"""RNA-RNA interactions on Rfam consensus positions (Interactions tab): interactions.json.

Production version of the RF00167 investigation in experiments/rf00167_interactions (map_interactions.py, analyse.py); the definitions
are the same and tests/test_interactions_rf00167.py checks that the family-wide numbers reproduce RF00167_interaction_conservation.csv.

Source: BGSU RNA 3D Hub FR3D annotations per PDB entry (base pairs and base stacking, CSV of unit-id, class, unit-id rows), fetched through
the cache; a seed directory (``<seed>/bgsu/{PDB}_{basepairs|stacking}.csv.gz`` + ``<seed>/PROVENANCE.json``) can provide them offline.
Unit id: PDB|Model|Chain|ResName|ResNum|Atom|AltId|InsCode|SymOp. Join: (pdb, author chain) -> instance; (author number, insertion code)
-> validated mapping residue -> label_seq_id -> Rfam state. Every mapped endpoint is checked against the instance's mmCIF (label_asym_id,
label_seq_id, residue name, model 1); a mismatch is counted and kept, never silently re-mapped.

Each interaction is written once, 5'->3' by label_seq_id (the FR3D class is flipped when a row is read from the 3' nucleotide: cWH <-> cHW,
s35 <-> s53), de-duplicated over the two row orders and alternate locations. FR3D classes are kept as annotated (no merging of classes).

Eligibility of an instance for a consensus position pair (i, j): the chain has FR3D annotations, both positions are match-state residues
with coordinates and a complete base (FR3D standard base heavy atoms; for modified residues the six pyrimidine-ring atoms), model 1.
Observed = at least one non-near FR3D interaction of the category between i and j in an eligible instance. Near interactions (FR3D 'n'
prefix) are counted separately and never as observed. BGSU's base-pair files contain no near base pairs; they are added from
near_basepairs.json (pinned local fr3d-python run, `rna_explorer near`) when present, tagged with their source and never used to decide
whether BGSU annotated a chain (annotation availability and gaps are BGSU-only).

Base pairs: paired = observed or near only (any class); `consensus` compares the position pair with the Rfam SS_cons: both = SS_cons pair
paired in at least one eligible chain, 2d = SS_cons pair paired in none (it gets a row even then), 3d = observed pair not in SS_cons. Insertion and flank endpoints are kept in the per-instance observations but are not
aggregated on consensus pairs. Missing, incomplete or unmapped nucleotides and unannotated chains are never counted as absence.
"""
from __future__ import annotations
import collections, csv, datetime as _dt, gzip, hashlib, io, json, logging, re
from pathlib import Path
from typing import Any, Optional

import requests

from . import __version__
from .cache import Cache
from .export import validate_json
from .near_basepairs import load_near_basepairs
from .pdbe import family_uses_chain_files, instance_cif

log = logging.getLogger(__name__)
BGSU_URL = "https://rna.bgsu.edu/rna3dhub/pdb/{pdb}/interactions/fr3d/{kind}/csv"
KINDS = {"basepairs": "basepair", "stacking": "stacking", "baseribose": "baseribose", "basephosphate": "basephosphate"}
# base-ribose and base-phosphate rows are directional (BGSU row order: base nucleotide, then the nucleotide whose ribose / phosphate is
# contacted); they are kept in that order (never flipped or re-ordered by sequence), and both directions between two nucleotides are distinct
DIRECTIONAL = {"baseribose", "basephosphate"}
# backbone atoms needed to evaluate a base-ribose / base-phosphate interaction on the contacted nucleotide (model 1)
SUGAR_ATOMS = frozenset("C1' C2' C3' C4' O4' O2'".split())
PHOSPHATE_ATOMS = (frozenset("P OP1 OP2 O5'".split()), frozenset("P O1P O2P O5'".split()))
# FR3D standard base heavy atoms (complete base = all present in model 1); modified residues: the six atoms of the pyrimidine ring
BASE_ATOMS = {"A": "N9 C8 N7 C5 C6 N6 N1 C2 N3 C4", "G": "N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4", "C": "N1 C2 O2 N3 C4 N4 C5 C6", "U": "N1 C2 O2 N3 C4 O4 C5 C6"}
RING = frozenset("N1 C2 N3 C4 C5 C6".split())
CANONICAL = {"AU", "UA", "GC", "CG", "GU", "UG"}
# per-position eligibility codes of an instance (string of length clen)
ELIGIBLE, INCOMPLETE, UNMODELLED, UNMAPPED = "E", "I", "M", "-"
# incomplete source annotation: a run of consecutive eligible positions without any BGSU annotation (any category, near included) at least
# this long is treated as "annotation unavailable" rather than as absence (e.g. RF00177 7PJS chain a, annotated only up to nt ~650)
GAP_MIN_POSITIONS, GAP_MIN_FRACTION = 15, 0.08


def annotation_gaps(eligibility: str, annotated: set[int]) -> list[list[int]]:
    """[start, end] Rfam ranges of long runs of eligible positions that carry no annotation at all (incomplete source annotation)"""
    el = [k + 1 for k, c in enumerate(eligibility) if c == ELIGIBLE]
    if not el:
        return []
    need = max(GAP_MIN_POSITIONS, GAP_MIN_FRACTION * len(el))
    gaps, run = [], []
    for p in el + [None]:
        if p is not None and p not in annotated:
            run.append(p); continue
        if len(run) >= need:
            gaps.append([run[0], run[-1]])
        run = []
    return gaps


def parse_unit(u: str) -> dict[str, Any]:
    f = u.split("|") + [""] * 9
    return {"unit": u, "pdb": f[0].upper(), "model": f[1], "chain": f[2], "resname": f[3], "num": int(f[4]), "alt": f[6], "ins": f[7], "sym": f[8]}


def flip(cls: str) -> str:
    """class of the same interaction written from the other nucleotide: cWH <-> cHW, s35 <-> s53 (cWW, s33, s55 unchanged)"""
    m = re.fullmatch(r"(n?)([ct])([WHS])([WHS])", cls)
    if m:
        return f"{m[1]}{m[2]}{m[4]}{m[3]}"
    m = re.fullmatch(r"(n?)s(\d)(\d)", cls)
    if m:
        return f"{m[1]}s{m[3]}{m[2]}"
    raise ValueError(f"unknown FR3D class {cls!r}")


def ss_cons_pairs(ss: str) -> set[tuple[int, int]]:
    st: dict[str, list[int]] = collections.defaultdict(list)
    pairs = set()
    for i, ch in enumerate(ss, 1):
        for o, c in ("()", "<>", "[]", "{}"):
            if ch == o:
                st[o].append(i)
            elif ch == c and st[o]:
                pairs.add((st[o].pop(), i))
    return pairs


# ----------------------------------------------------------------------------- sources
class BgsuSource:
    """FR3D CSVs per entry through the cache; optional seeding from a directory of previously downloaded gzip files."""

    def __init__(self, cache: Cache, seed_dir: Optional[Path] = None, offline: bool = False) -> None:
        self.cache, self.seed_dir, self.offline = cache, Path(seed_dir) if seed_dir else None, offline
        self.seed_prov: dict[str, Any] = {}
        if self.seed_dir and (self.seed_dir / "PROVENANCE.json").exists():
            self.seed_prov = json.load(open(self.seed_dir / "PROVENANCE.json"))

    def get(self, pdb: str, kind: str) -> tuple[Optional[str], dict[str, Any]]:
        url = BGSU_URL.format(pdb=pdb.upper(), kind=kind)
        path = self.cache._path(url, ".csv", "GET", None)
        seeded = None
        if not path.exists() and self.seed_dir:
            src = self.seed_dir / "bgsu" / f"{pdb.upper()}_{kind}.csv.gz"
            if src.is_file():
                path.write_bytes(gzip.decompress(src.read_bytes()))
                seeded = str(src)
        meta: dict[str, Any] = {"pdb_id": pdb.upper(), "kind": KINDS[kind], "url": url}
        if not path.exists() and self.offline:
            return None, {**meta, "status": "unavailable", "note": "offline and not cached"}
        try:
            text = self.cache.fetch_text(url, ".csv")
        except requests.RequestException as e:
            return None, {**meta, "status": "unavailable", "note": str(e)[:200]}
        p = self.seed_prov.get(url) or {}
        retrieved = p.get("retrieved") or _dt.datetime.fromtimestamp(path.stat().st_mtime, _dt.timezone.utc).date().isoformat()
        return text, {**meta, "status": "ok", "retrieved": retrieved, "sha256": hashlib.sha256(text.encode()).hexdigest(), "seeded_from": seeded or p.get("file")}


def chain_residues(cif_path: Path, label_asym_id: str) -> dict[int, dict[str, Any]]:
    """label_seq_id -> author chain, number, insertion code, name, base completeness (model 1) of one polymer chain."""
    import gemmi
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    out: dict[int, dict[str, Any]] = {}
    for ch in st[0]:
        for res in ch:
            if res.subchain != label_asym_id or res.label_seq is None:
                continue
            names = {a.name for a in res}
            need = set(BASE_ATOMS[res.name].split()) if res.name in BASE_ATOMS else RING
            out.setdefault(res.label_seq, {"chain": ch.name, "num": res.seqid.num, "ins": res.seqid.icode.strip(), "name": res.name, "complete": need <= names,
                                           "sugar": SUGAR_ATOMS <= names, "phosphate": any(ph <= names for ph in PHOSPHATE_ATOMS)})
    return out


SOURCE_BGSU, SOURCE_LOCAL = "bgsu", "fr3d_local"


# ----------------------------------------------------------------------------- build
def build_interactions(acc: str, cache_dir: Path, out_dir: Path, seed_dir: Optional[Path] = None, offline: bool = False, write: bool = True) -> dict[str, Any]:
    fam_dir = Path(out_dir) / "families" / acc
    family = json.load(open(fam_dir / "family.json"))
    instances = json.load(open(fam_dir / "instances.json"))
    clen = family["clen"]
    ref_el = fam_dir / "reference_structural_elements.json"
    region = json.load(open(ref_el))["track"] if ref_el.exists() else [None] * clen
    ss = ss_cons_pairs(family.get("ss_cons", ""))
    cache = Cache(cache_dir)
    chain_files = family_uses_chain_files(fam_dir)
    src = BgsuSource(cache, seed_dir, offline)
    near_file = load_near_basepairs(fam_dir)
    near_entries = (near_file or {}).get("entries", {})

    by_sn = {i["seq_name"]: i for i in instances}
    by_chain: dict[tuple[str, str], list[dict[str, Any]]] = collections.defaultdict(list)
    for i in instances:
        by_chain[(i["pdb_id"].upper(), i["auth_asym_id"])].append(i)
    maps = {sn: json.load(open(fam_dir / by_sn[sn]["mapping_file"])) for sn in by_sn}
    resmap = {sn: {(r["author_residue_number"], r["author_insertion_code"] or ""): r for r in m["residues"]} for sn, m in maps.items()}

    # coordinates: per instance label_seq_id -> residue facts (mmCIF used by the pipeline, model 1)
    coords: dict[str, dict[int, dict[str, Any]]] = {}
    for sn, inst in by_sn.items():
        try:
            coords[sn] = chain_residues(instance_cif(cache, inst, chain_files), inst["label_asym_id"])
        except Exception as e:   # noqa: BLE001 - recorded, the chain is then not comparable anywhere
            log.warning("coordinates unavailable for %s: %s", sn, e)
            coords[sn] = {}

    sources: list[dict[str, Any]] = []
    endpoint_status: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    verification: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    pair_kinds: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    rows_in_files: collections.Counter[str] = collections.Counter()
    obs: dict[tuple, dict[str, Any]] = {}
    inter_chain: dict[tuple, dict[str, Any]] = {}
    annotated: set[str] = set()
    entries_ok: dict[str, bool] = {}               # base-pair and stacking files (chain annotation availability, as validated)
    files_ok: dict[tuple[str, str], bool] = {}     # (pdb, category) -> file retrieved

    def resolve(p: dict[str, Any]) -> tuple[str, Optional[dict[str, Any]], Optional[dict[str, Any]]]:
        if p["model"] != "1":
            return "other_model", None, None
        if p["sym"] and p["sym"] not in ("1_555", "P_1"):
            return "symmetry_mate", None, None
        cands = by_chain.get((p["pdb"], p["chain"]), [])
        if not cands:
            return "other_chain_or_non_polymer", None, None
        if len(cands) > 1:
            return "ambiguous_chain", None, None
        inst = cands[0]
        r = resmap[inst["seq_name"]].get((p["num"], p["ins"]))
        if r is None:
            return "not_in_rna_polymer", inst, None
        return "mapped", inst, r

    local_endpoints: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    local_verification: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    local_kinds: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)

    def ingest(cat: str, raws: Any, source: str) -> int:
        """map and record BGSU-style rows (unit id, class, unit id) of one category; returns the number of rows read"""
        n_rows = 0
        for raw in raws:
            if len(raw) < 3:
                continue
            n_rows += 1
            if source == SOURCE_BGSU:
                rows_in_files[cat] += 1
            u1, cls, u2 = raw[0], raw[1], raw[2]
            ends = []
            for u in (u1, u2):
                p = parse_unit(u)
                st, inst, r = resolve(p)
                e = {"p": p, "inst": inst, "r": r, "status": st, "loc": None}
                if st == "mapped":
                    c = coords[inst["seq_name"]].get(r["polymer_seq_position"])
                    ok = c is not None and c["chain"] == p["chain"] and c["num"] == p["num"] and c["ins"] == p["ins"] and c["name"] == p["resname"]
                    e["check"] = "verified" if ok else "mismatch"
                    state = r["alignment_state"]
                    e["loc"] = r["rfam_consensus_pos"] if state == "match" else (f"{r['insertion_anchor']}^{r['insertion_offset']}" if state == "insert" else None)
                    st = {"match": "consensus", "insert": "insertion"}.get(state, "outside_alignment")
                (endpoint_status if source == SOURCE_BGSU else local_endpoints)[cat][st] += 1
                if "check" in e:
                    (verification if source == SOURCE_BGSU else local_verification)[cat][e["check"]] += 1
                ends.append(e)
            a, b = ends
            if a["r"] is None or b["r"] is None:
                (pair_kinds if source == SOURCE_BGSU else local_kinds)[cat]["unmapped_endpoint"] += 1
                continue
            if a["inst"]["seq_name"] != b["inst"]["seq_name"]:
                (pair_kinds if source == SOURCE_BGSU else local_kinds)[cat]["between_family_chains"] += 1
                k = (tuple(sorted((a["p"]["unit"], b["p"]["unit"]))) if cat not in DIRECTIONAL else (a["p"]["unit"], b["p"]["unit"])) + (cat, cls)
                if source == SOURCE_BGSU:
                    inter_chain.setdefault(k, {"category": cat, "fr3d_class": cls, "unit_1": a["p"]["unit"], "unit_2": b["p"]["unit"]})
                continue
            (pair_kinds if source == SOURCE_BGSU else local_kinds)[cat]["intra_chain"] += 1
            if cat not in DIRECTIONAL and a["r"]["polymer_seq_position"] > b["r"]["polymer_seq_position"]:
                a, b, cls = b, a, flip(cls)
            sn = a["inst"]["seq_name"]
            if cat not in DIRECTIONAL and source == SOURCE_BGSU:
                annotated.add(sn)
            k = (sn, cat, a["r"]["polymer_seq_position"], b["r"]["polymer_seq_position"], cls)
            rec = obs.get(k)
            if rec is None:
                obs[k] = {"seq_name": sn, "category": cat, "fr3d_class": cls, "a": a, "b": b, "alts": set(), "source": source}
                rec = obs[k]
            alt = a["p"]["alt"] or b["p"]["alt"]
            if alt:
                rec["alts"].add(alt)
        return n_rows

    near_status: dict[str, str] = {}
    for pdb in sorted({i["pdb_id"].upper() for i in instances}):
        ok_all = True
        for kind, cat in KINDS.items():
            text, meta = src.get(pdb, kind)
            files_ok[(pdb, cat)] = text is not None
            n_rows = 0
            if text is None:
                ok_all = ok_all and cat in DIRECTIONAL
            else:
                n_rows = ingest(cat, csv.reader(io.StringIO(text)), SOURCE_BGSU)
            sources.append({**meta, "n_rows": n_rows})
        entries_ok[pdb] = ok_all
        ne = near_entries.get(pdb)
        near_status[pdb] = "available" if ne and ne["status"] == "ok" else "unavailable"
        if ne and ne["status"] == "ok":
            ingest("basepair", [r for r in ne["rows"] if r[1].startswith("n")], SOURCE_LOCAL)

    # per-instance observations (all intra-chain interactions, including insertion and flank endpoints)
    observations: dict[str, list[list[Any]]] = {sn: [] for sn in sorted(by_sn)}
    for (sn, cat, l1, l2, cls), rec in sorted(obs.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][3], kv[0][4])):
        a, b = rec["a"], rec["b"]
        observations[sn].append([cat, cls, l1, l2, a["loc"], b["loc"], a["p"]["unit"], b["p"]["unit"],
                                 "verified" if a["check"] == b["check"] == "verified" else "mismatch", ";".join(sorted(rec["alts"])), rec["source"]])

    # per-instance eligibility and nucleotide identity on consensus positions
    inst_rows = []
    elig: dict[str, str] = {}
    backbone: dict[str, dict[str, str]] = {}
    for sn in sorted(by_sn):
        inst = by_sn[sn]
        codes, nts = [UNMAPPED] * clen, ["-"] * clen
        sugar, phos = [UNMAPPED] * clen, [UNMAPPED] * clen
        for r in maps[sn]["residues"]:
            if r["alignment_state"] != "match" or r["rfam_consensus_pos"] is None:
                continue
            k = r["rfam_consensus_pos"] - 1
            c = coords[sn].get(r["polymer_seq_position"])
            nts[k] = r["nucleotide_parent"] or "N"
            codes[k] = UNMODELLED if not (r["has_coordinates"] and c) else (ELIGIBLE if c["complete"] else INCOMPLETE)
            sugar[k] = UNMODELLED if not (r["has_coordinates"] and c) else (ELIGIBLE if c["sugar"] else INCOMPLETE)
            phos[k] = UNMODELLED if not (r["has_coordinates"] and c) else (ELIGIBLE if c["phosphate"] else INCOMPLETE)
        elig[sn] = "".join(codes)
        backbone[sn] = {"baseribose": "".join(sugar), "basephosphate": "".join(phos)}
        if not entries_ok[inst["pdb_id"].upper()]:
            ann, note = "unavailable", "BGSU FR3D annotation files could not be retrieved for this entry"
        elif sn not in annotated:
            ann, note = "unavailable", "no BGSU FR3D interaction records for this chain"
        else:
            ann, note = "available", None
        inst_rows.append({"seq_name": sn, "pdb_id": inst["pdb_id"].upper(), "auth_asym_id": inst["auth_asym_id"], "label_asym_id": inst["label_asym_id"],
                          "annotation": ann, "annotation_note": note, "eligibility": elig[sn], "sugar_eligibility": backbone[sn]["baseribose"], "phosphate_eligibility": backbone[sn]["basephosphate"],
                          "unavailable_categories": [c for c in KINDS.values() if ann != "available" or not files_ok[(inst["pdb_id"].upper(), c)]],
                          "nucleotides": "".join(nts), "n_interactions": len(observations[sn]),
                          "near_basepairs": near_status[inst["pdb_id"].upper()]})
    unavailable = {r["seq_name"]: set(r["unavailable_categories"]) for r in inst_rows}
    # incomplete BGSU annotation of a chain (long unannotated runs of eligible positions): not comparable there, never absent
    for r in inst_rows:
        annotated_pos = {q for o in observations[r["seq_name"]] if o[10] == SOURCE_BGSU for q in (o[4], o[5]) if isinstance(q, int)}
        r["annotation_gaps"] = annotation_gaps(r["eligibility"], annotated_pos) if r["annotation"] == "available" else []
    in_gap = {r["seq_name"]: {q for a, b in r["annotation_gaps"] for q in range(a, b + 1)} for r in inst_rows}
    ann_ok = {r["seq_name"] for r in inst_rows if r["annotation"] == "available"}

    # family-wide consensus pairs
    cons: dict[tuple[str, int, int], dict[str, list[str]]] = collections.defaultdict(lambda: collections.defaultdict(list))
    near: dict[tuple[str, int, int], dict[str, set[str]]] = collections.defaultdict(lambda: collections.defaultdict(set))
    for sn, rows in observations.items():
        for cat, cls, _l1, _l2, p1, p2, *_ in rows:
            if not (isinstance(p1, int) and isinstance(p2, int)) or p1 == p2:
                continue
            k = (cat, p1, p2)
            if cls.startswith("n"):
                near[k][sn].add(cls)
            else:
                cons[k][sn].append(cls)
    nt = {r["seq_name"]: r["nucleotides"] for r in inst_rows}

    def eligibility(sn: str, cat: str, i: int, j: int) -> str:
        """i, j: the two endpoints (for directional categories i = base nucleotide, j = nucleotide whose ribose / phosphate is contacted)"""
        if sn not in ann_ok or cat in unavailable[sn] or i in in_gap[sn] or j in in_gap[sn]:
            return "annotation_unavailable"
        a, b = elig[sn][i - 1], (backbone[sn][cat] if cat in DIRECTIONAL else elig[sn])[j - 1]
        if UNMAPPED in (a, b):
            return "unmapped"
        if UNMODELLED in (a, b):
            return "missing"
        if INCOMPLETE in (a, b):
            return "incomplete"
        return "eligible"

    # every SS_cons position pair gets a base-pair row, also when it is paired (true or near) in no chain ("2d")
    for i, j in ss:
        cons[("basepair", i, j)]
    pairs = []
    for (cat, i, j), per in sorted(cons.items()):
        st: collections.Counter[str] = collections.Counter()
        classes: collections.Counter[str] = collections.Counter()
        near_classes: collections.Counter[str] = collections.Counter()
        ent_p = set()
        nt_pairs: collections.Counter[str] = collections.Counter()
        ent_e, ent_o = set(), set()
        n_canonical = obs_inel = 0
        for sn in sorted(by_sn):
            e = eligibility(sn, cat, i, j)
            st[e] += 1
            cl = per.get(sn)
            if e != "eligible":
                obs_inel += bool(cl)
                continue
            ent_e.add(by_sn[sn]["pdb_id"].upper())
            if cl:
                st["observed"] += 1
                ent_o.add(by_sn[sn]["pdb_id"].upper())
                pair_nt = nt[sn][i - 1] + nt[sn][j - 1]
                nt_pairs[pair_nt] += 1
                for c in sorted(set(cl)):
                    classes[c] += 1
                n_canonical += "cWW" in cl and pair_nt in CANONICAL
            elif sn in near[(cat, i, j)]:
                st["near_only"] += 1
                for c in sorted(near[(cat, i, j)][sn]):
                    near_classes[c] += 1
                if cat == "basepair":
                    ent_p.add(by_sn[sn]["pdb_id"].upper())
            else:
                st["not_annotated"] += 1
        n, k = st["eligible"], st["observed"]
        ri, rj = region[i - 1], region[j - 1]
        in_ss = cat in ("basepair", "stacking") and (i, j) in ss
        structure_class = None if cat != "basepair" else "secondary" if in_ss else "tertiary" if (ri and rj and ri != rj) else "unclassified"
        kp = k + (st["near_only"] if cat == "basepair" else 0)
        consensus = None if cat != "basepair" else ("both" if kp else "2d") if in_ss else "3d"
        pairs.append({"id": f"{cat}:{i}-{j}", "category": cat, "pos_1": i, "pos_2": j, "directional": cat in DIRECTIONAL,
                      "consensus_nt": f"{family['consensus_sequence'][i - 1]}{family['consensus_sequence'][j - 1]}" if family.get("consensus_sequence") else None,
                      "region_1": ri, "region_2": rj, "in_ss_cons": in_ss, "structure_class": structure_class,
                      "n_observed": k, "n_eligible": n, "frequency": round(k / n, 4) if n else None,
                      "n_entries_observed": len(ent_o), "n_entries_eligible": len(ent_e),
                      "consensus": consensus, "n_paired": kp, "frequency_paired": round(kp / n, 4) if n else None,
                      "n_entries_paired": len(ent_o | ent_p) if cat == "basepair" else len(ent_o), "near_classes": dict(near_classes.most_common()),
                      "classes": dict(classes.most_common()), "nt_pairs": dict(nt_pairs.most_common()), "n_canonical": n_canonical,
                      "n_near_only": st["near_only"], "n_not_annotated": st["not_annotated"], "n_missing": st["missing"], "n_incomplete": st["incomplete"],
                      "n_unmapped": st["unmapped"], "n_annotation_unavailable": st["annotation_unavailable"], "n_observed_but_ineligible": obs_inel})

    ref = family.get("display_representative") or family.get("representative")
    ref_sn = ref["seq_name"] if ref else None
    out = {
        "schema_version": "1.4", "rfam_acc": acc, "clen": clen,
        "default_numbering_reference": {"seq_name": ref_sn, "rule": "family display representative (family.json display_representative)",
                                        "n_eligible_positions": elig.get(ref_sn, "").count(ELIGIBLE) if ref_sn else 0,
                                        "n_match_positions": sum(c != UNMAPPED for c in elig.get(ref_sn, "")) if ref_sn else 0},
        "position_regions": region,
        "definitions": {
            "observation": "[category, fr3d_class, label_seq_id_1, label_seq_id_2, rfam_1, rfam_2, unit_id_1, unit_id_2, mapping_check, altlocs]; base pairs and stacking "
                           "5'->3' by label_seq_id (class flipped when needed); base-ribose and base-phosphate in BGSU row order: endpoint 1 = base, endpoint 2 = "
                           "nucleotide whose ribose / phosphate is contacted (never re-ordered); "
                           "rfam_k = consensus position (match), 'anchor^offset' (insertion) or null (outside the aligned region)",
            "eligibility_codes": {ELIGIBLE: "match residue with coordinates and complete base", INCOMPLETE: "match residue with an incomplete base",
                                  UNMODELLED: "match residue without coordinates", UNMAPPED: "no match residue (deletion or outside the aligned region)"},
            "eligible": "chain annotated by FR3D and the category's file retrieved; base pairs and stacking: both positions match residues with coordinates and "
                        "complete bases; base-ribose: base nucleotide with a complete base and contacted nucleotide with a complete ribose (C1' C2' C3' C4' O4' "
                        "O2'); base-phosphate: complete base and a complete phosphate (P OP1 OP2 O5') on the contacted nucleotide (model 1)",
            "sugar_eligibility": "per consensus position, ribose completeness (codes as eligibility, I = incomplete ribose)",
            "phosphate_eligibility": "per consensus position, phosphate completeness (codes as eligibility, I = incomplete phosphate, e.g. a 5'-terminal residue without P)",
            "annotation_gaps": f"per chain, [start, end] Rfam ranges of >= {GAP_MIN_POSITIONS} (and >= {GAP_MIN_FRACTION:.0%} of the chain's eligible) consecutive eligible "
                               "positions without any BGSU annotation (incomplete source annotation); interactions touching them count as annotation unavailable",
            "pair_id": "category:pos_1-pos_2; base pairs and stacking pos_1 < pos_2; base-ribose / base-phosphate pos_1 = base, pos_2 = contacted backbone",
            "observed": "a non-near FR3D interaction of the category between the two positions in an eligible instance; near ('n') classes are counted as near_only",
            "structure_class": "base pairs only: secondary = position pair of the Rfam SS_cons; tertiary = other pair between two different validated structural "
                               "elements; unclassified = otherwise",
            "canonical": "n_canonical = observed instances annotated cWW with an AU, GC or GU nucleotide pair",
            "paired": "base pairs only: n_paired = observed + near only (an FR3D base pair of any class, true or near, between the two positions in an "
                      "eligible chain); other categories: n_paired = n_observed (near never counted)",
            "consensus": "base pairs only: both = Rfam SS_cons position pair paired in at least one eligible chain; 2d = SS_cons pair paired in no "
                         "eligible chain; 3d = observed base pair that is not an SS_cons pair; null for other categories",
            "near_classes": "near-only eligible chains per near FR3D class",
            "near_basepairs": "per chain: available = near base pairs computed by the pinned local fr3d-python run for its entry (near_basepairs.json); "
                              "unavailable = not computed (near base pairs cannot be observed there)",
            "observation_source": f"last observation field: {SOURCE_BGSU} (BGSU RNA 3D Hub file) or {SOURCE_LOCAL} (near base pair, local fr3d-python)",
        },
        "instances": inst_rows,
        "observations_files": {"path": "interactions_by_instance/{seq_name}.json", "content": "every intra-chain observation of the chain (definitions.observation)"},
        "pairs": pairs,
        "excluded": {"between_family_chains": sorted(inter_chain.values(), key=lambda r: (r["unit_1"], r["unit_2"], r["category"]))},
        "summary": {"rows_in_files": dict(rows_in_files),
                    "by_category": {c: {"rows_in_files": rows_in_files[c], "endpoint_status": dict(endpoint_status[c]), "endpoint_verification_against_mmcif": dict(verification[c]),
                                        "pair_kinds": dict(pair_kinds[c]), "n_intra_chain_interactions": sum(1 for v in observations.values() for o in v if o[0] == c and o[10] == SOURCE_BGSU),
                                        "n_consensus_pairs": sum(1 for q in pairs if q["category"] == c),
                                        "n_entries_with_file": sum(1 for (e, cc), ok in files_ok.items() if cc == c and ok)} for c in KINDS.values()},
                    "n_intra_chain_interactions": sum(1 for v in observations.values() for o in v if o[10] == SOURCE_BGSU), "n_consensus_pairs": len(pairs),
                    "n_instances": len(by_sn), "n_instances_annotated": len(ann_ok), "n_entries": len(entries_ok),
                    "consensus": dict(collections.Counter(q["consensus"] for q in pairs if q["consensus"])),
                    "near_basepairs_local": {"n_intra_chain_observations": sum(1 for v in observations.values() for o in v if o[10] == SOURCE_LOCAL),
                                             "endpoint_status": dict(local_endpoints["basepair"]), "endpoint_verification_against_mmcif": dict(local_verification["basepair"]),
                                             "pair_kinds": dict(local_kinds["basepair"])}},
        "sources": sources,
        "provenance": {"source": "BGSU RNA 3D Hub FR3D annotations (base pairs, base stacking, base-ribose, base-phosphate)", "url_template": BGSU_URL,
                       "coordinates": "PDBe ModelServer chain files" if chain_files else "PDBe updated mmCIF (model 1)",
                       "rfam_release": family.get("provenance", {}).get("rfam_release"), "mapping_generated": family.get("provenance", {}).get("generated"),
                       "pipeline_version": __version__, "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                       "method": "experiments/rf00167_interactions (RF00167 investigation, 2026-10-18)",
                       "near_basepairs": None if near_file is None else {"file": "near_basepairs.json", **{k: near_file["provenance"][k] for k in ("tool", "repository", "commit", "commit_date", "licence", "generated")},
                                                                         "n_entries_available": sum(v == "available" for v in near_status.values()), "n_entries": len(near_status)}},
    }
    # compact per-chain view for the frontend (the main file stays small for large RNAs): indices into `pairs` of every non-near observed
    # consensus interaction by class, near annotations on consensus positions, and the number of off-consensus interactions; the full
    # records (BGSU unit ids, label_seq_ids, mmCIF check, altlocs) are kept per chain in interactions_by_instance/<seq_name>.json
    pair_index = {(q["category"], q["pos_1"], q["pos_2"]): k for k, q in enumerate(pairs)}
    for r in inst_rows:
        observed: dict[str, list[int]] = collections.defaultdict(list)
        near_rows, off = [], 0
        for cat, cls, _l1, _l2, p1, p2, *_ in observations[r["seq_name"]]:
            if not (isinstance(p1, int) and isinstance(p2, int)) or p1 == p2:
                off += 1
            elif cls.startswith("n"):
                near_rows.append([cat, cls, p1, p2])
            else:
                k = pair_index[(cat, p1, p2)]
                if k not in observed[cls]:
                    observed[cls].append(k)
        r["observed"] = {c: sorted(v) for c, v in sorted(observed.items())}
        r["near"] = near_rows
        r["n_off_consensus"] = off
    validate_json(out, "interactions")
    per_instance = {sn: {"schema_version": out["schema_version"], "rfam_acc": acc, "seq_name": sn, "observation": out["definitions"]["observation"] + "; " + out["definitions"]["observation_source"],
                         "observations": rows}
                    for sn, rows in observations.items()}
    for v in per_instance.values():
        validate_json(v, "interactions_instance")
    if write:
        with open(fam_dir / "interactions.json", "w") as f:
            json.dump(out, f, separators=(",", ":"))
        d = fam_dir / "interactions_by_instance"
        d.mkdir(exist_ok=True)
        for old in d.glob("*.json"):
            if old.stem not in per_instance:
                old.unlink()
        for sn, v in per_instance.items():
            with open(d / f"{sn}.json", "w") as f:
                json.dump(v, f, separators=(",", ":"))
    return {**out, "observations": observations}
