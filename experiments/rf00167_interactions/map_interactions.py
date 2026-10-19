"""RF00167 RNA-RNA interaction investigation, stage 1: map BGSU/FR3D interaction records onto RNA chain instances and Rfam positions.

Inputs (read-only):
- BGSU FR3D per-entry CSVs (base pairs, stacking) cached by the research worktree's bp_fetch.py
  (~/Desktop/aggregated_rna_clustering/.../aggregated_secondary_structure_investigation/raw/bgsu, PROVENANCE.json there).
- This repository's validated per-chain mappings (data/families/RF00167/mappings) and instances.json.
- The PDBe updated mmCIF files already in the local cache (cache/www_ebi_ac_uk; git-ignored, not copied) for verification.

Unit ID (BGSU): PDB|Model|Chain|ResName|ResNum|Atom|AltId|InsCode|SymOp (trailing fields optional).
Join: (pdb, auth chain) -> instance; (author number, insertion code) -> mapped residue -> polymer position (label_seq_id) -> Rfam state.
Each endpoint is verified against _pdbx_poly_seq_scheme (label asym, seq_id, mon_id) and _atom_site (model 1).
Writes: endpoint_mapping.csv, RF00167_mapped_interactions.csv, base_completeness.csv, stage1_summary.json (this directory).
Usage: python3 experiments/rf00167_interactions/map_interactions.py   (from the repository root)
"""
import csv, glob, gzip, json, collections, re
from pathlib import Path
import gemmi

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
FAM = ROOT / "data/families/RF00167"
BGSU = Path.home() / "Desktop/aggregated_rna_clustering/experiments/rna_structural_clustering/results/aggregated_secondary_structure_investigation/raw"
CIF = ROOT / "cache/www_ebi_ac_uk"
# FR3D standard base heavy atoms (complete base = all present in model 1); modified residues: the six ring atoms of the pyrimidine ring
BASE = {"A": "N9 C8 N7 C5 C6 N6 N1 C2 N3 C4", "G": "N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4", "C": "N1 C2 O2 N3 C4 N4 C5 C6", "U": "N1 C2 O2 N3 C4 O4 C5 C6"}
RING = set("N1 C2 N3 C4 C5 C6".split())


def parse_unit(u):
    f = u.split("|") + [""] * 9
    return {"unit": u, "pdb": f[0].upper(), "model": f[1], "chain": f[2], "resname": f[3], "num": int(f[4]), "atom": f[5], "alt": f[6], "ins": f[7], "sym": f[8], "n_fields": len(u.split("|"))}


def flip(cls):
    """class of the same interaction written from the other nucleotide: cWH <-> cHW, s35 <-> s53 (s33, s55, cWW unchanged)"""
    m = re.fullmatch(r"(n?)([ct])([WHS])([WHS])", cls)
    if m:
        return f"{m[1]}{m[2]}{m[4]}{m[3]}"
    m = re.fullmatch(r"(n?)s(\d)(\d)", cls)
    if m:
        return f"{m[1]}s{m[3]}{m[2]}"
    raise ValueError(cls)


def main():
    instances = json.load(open(FAM / "instances.json"))
    by_chain = collections.defaultdict(list)
    for i in instances:
        by_chain[(i["pdb_id"].upper(), i["auth_asym_id"])].append(i)
    maps = {i["seq_name"]: json.load(open(FAM / "mappings" / f"{i['seq_name']}.json")) for i in instances}
    resmap = {sn: {(r["author_residue_number"], r["author_insertion_code"] or ""): r for r in m["residues"]} for sn, m in maps.items()}
    prov = json.load(open(BGSU / "PROVENANCE.json"))
    url_of = {v["file"].split("/")[-1]: (k, v["retrieved"], v["sha256"]) for k, v in prov.items() if v.get("file")}

    # mmCIF: poly_seq_scheme index and model-1 base completeness, per entry
    pss, present, complete = {}, set(), {}
    for pdb in sorted({i["pdb_id"].upper() for i in instances}):
        path = glob.glob(str(CIF / f"pdbe_entry-files_download_{pdb.lower()}_updated.cif-*"))[0]
        doc = gemmi.cif.read(path); b = doc.sole_block()
        t = b.find("_pdbx_poly_seq_scheme.", ["asym_id", "seq_id", "mon_id", "pdb_seq_num", "pdb_ins_code", "pdb_strand_id"])
        for r in t:
            ins = "" if r[4] in ("?", ".") else r[4]
            pss[(pdb, r[5], int(r[3]), ins)] = (r[0], int(r[1]), r[2])
        st = gemmi.make_structure_from_block(b); model = st[0]
        for ch in model:
            for res in ch:
                ins = res.seqid.icode.strip()
                key = (pdb, ch.name, res.seqid.num, ins)
                present.add(key)
                names = {a.name for a in res}
                ci = gemmi.find_tabulated_residue(res.name)
                parent = res.name if res.name in BASE else (ci.one_letter_code.upper() if ci and ci.one_letter_code.upper() in BASE else None)
                need = set(BASE[res.name].split()) if res.name in BASE else RING
                complete[key] = (parent, need <= names)

    status = collections.Counter(); verify = collections.Counter(); examples = collections.defaultdict(list)
    endpoint_rows, obs = [], {}
    rows_total = collections.Counter(); pair_kinds = collections.Counter(); reverse_check = collections.Counter()
    seen_rows = collections.defaultdict(set)

    def resolve(p):
        if p["model"] != "1":
            return "other_model", None, None
        if p["sym"] and p["sym"] not in ("1_555", "P_1"):
            return "symmetry_mate", None, None
        cands = by_chain.get((p["pdb"], p["chain"]), [])
        if not cands:
            # another chain of the entry: polymer (other RNA / protein / scaffold) or non-polymer
            return ("other_polymer_chain" if any(k[0] == p["pdb"] and k[1] == p["chain"] for k in pss) else "non_polymer_or_other_chain"), None, None
        if len(cands) > 1:
            return "ambiguous_chain", None, None
        inst = cands[0]
        r = resmap[inst["seq_name"]].get((p["num"], p["ins"]))
        if r is None:
            key = (p["pdb"], p["chain"], p["num"], p["ins"])
            return ("non_polymer_in_rna_chain" if key not in pss else "residue_not_in_mapping"), inst, None
        return "mapped", inst, r

    for f in sorted(glob.glob(str(BGSU / "bgsu" / "*.csv.gz"))):
        kind = "basepair" if "basepairs" in f else "stacking"
        src = url_of[Path(f).name]
        for raw in csv.reader(gzip.open(f, "rt")):
            if len(raw) < 3:
                continue
            rows_total[kind] += 1
            u1, cls, u2 = raw[0], raw[1], raw[2]
            seen_rows[kind].add((u1, cls, u2))
            res = []
            for u in (u1, u2):
                p = parse_unit(u)
                st, inst, r = resolve(p)
                v = None
                if st == "mapped":
                    key = (p["pdb"], p["chain"], p["num"], p["ins"])
                    s = pss.get(key)
                    ok = s is not None and s[0] == inst["label_asym_id"] and s[1] == r["polymer_seq_position"] and s[2] == p["resname"] and key in present
                    v = "verified" if ok else f"mismatch:{s}"
                    verify[v.split(':')[0]] += 1
                    state = r["alignment_state"]
                    pos = r["rfam_consensus_pos"] if state == "match" else (f"{r['insertion_anchor']}^{r['insertion_offset']}" if state == "insert" else None)
                    st2 = {"match": "consensus", "insert": "insertion"}.get(state, f"outside_alignment:{state}")
                    status[st2] += 1
                    if len(examples[st2]) < 3:
                        examples[st2].append({"unit": u, "seq_name": inst["seq_name"], "label_seq_id": r["polymer_seq_position"], "state": state, "rfam_pos": pos})
                    res.append({"p": p, "inst": inst, "r": r, "state": st2, "pos": pos, "verify": v})
                else:
                    status[st] += 1
                    if len(examples[st]) < 3:
                        examples[st].append({"unit": u, "class": cls, "partner": u2 if u == u1 else u1})
                    res.append({"p": p, "inst": inst, "r": None, "state": st, "pos": None, "verify": None})
            a, b = res
            if a["r"] is None or b["r"] is None:
                pair_kinds[("rna_rna" if a["state"] not in ("non_polymer_in_rna_chain", "non_polymer_or_other_chain") and b["state"] not in ("non_polymer_in_rna_chain", "non_polymer_or_other_chain") else "involves_non_polymer", "unmapped_endpoint")] += 1
                continue
            if a["inst"]["seq_name"] != b["inst"]["seq_name"]:
                pair_kinds[("rna_rna", "between_two_rf00167_chains")] += 1
                key = ("inter", kind, cls, a["p"]["unit"], b["p"]["unit"])
                obs.setdefault(key, (kind, cls, a, b, src, True))
                continue
            pair_kinds[("rna_rna", "intra_chain")] += 1
            # normalise: 5'->3' by label_seq_id; flip the class when the row is written from the 3' nucleotide
            if a["r"]["polymer_seq_position"] > b["r"]["polymer_seq_position"]:
                a, b, cls_n = b, a, flip(cls)
            else:
                cls_n = cls
            k = (kind, a["inst"]["seq_name"], a["r"]["polymer_seq_position"], b["r"]["polymer_seq_position"], cls_n)
            if k in obs:
                obs[k][4].add(a["p"]["alt"] or b["p"]["alt"] or "")
            else:
                obs[k] = [kind, cls_n, a, b, {a["p"]["alt"] or b["p"]["alt"] or ""}, src, cls]
    # reverse-order consistency: every row should have its mirror (u2, flip(cls), u1)
    for kind, rs in seen_rows.items():
        for (u1, c, u2) in rs:
            reverse_check[(kind, (u2, flip(c), u1) in rs)] += 1

    # mapped interactions table (intra-chain, normalised, de-duplicated over both orders and altlocs)
    out = []
    for k, v in obs.items():
        if k[0] == "inter":
            continue
        kind, cls_n, a, b, alts, src, cls_orig = v
        out.append({"seq_name": a["inst"]["seq_name"], "pdb_id": a["p"]["pdb"], "category": kind, "fr3d_class": cls_n, "near": cls_n.startswith("n"),
                    "unit_1": a["p"]["unit"], "unit_2": b["p"]["unit"], "auth_1": f"{a['p']['chain']}:{a['p']['resname']}{a['p']['num']}{a['p']['ins']}",
                    "auth_2": f"{b['p']['chain']}:{b['p']['resname']}{b['p']['num']}{b['p']['ins']}", "label_seq_1": a["r"]["polymer_seq_position"], "label_seq_2": b["r"]["polymer_seq_position"],
                    "nt_1": a["r"]["nucleotide_parent"], "nt_2": b["r"]["nucleotide_parent"], "state_1": a["state"], "state_2": b["state"], "rfam_pos_1": a["pos"], "rfam_pos_2": b["pos"],
                    "altlocs": ";".join(sorted(x for x in alts if x)), "mapping_check": "verified" if a["verify"] == b["verify"] == "verified" else f"{a['verify']}/{b['verify']}",
                    "source_url": src[0], "retrieved": src[1], "source_sha256": src[2]})
    out.sort(key=lambda r: (r["seq_name"], r["category"], r["label_seq_1"], r["label_seq_2"], r["fr3d_class"]))
    with open(OUT / "RF00167_mapped_interactions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    inter = [{"category": v[0], "fr3d_class": v[1], "unit_1": v[2]["p"]["unit"], "unit_2": v[3]["p"]["unit"]} for k, v in obs.items() if k[0] == "inter"]
    with open(OUT / "inter_chain_interactions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["category", "fr3d_class", "unit_1", "unit_2"]); w.writeheader(); w.writerows(inter)

    # base completeness of every mapped residue (eligibility input for stage 2)
    with open(OUT / "base_completeness.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["seq_name", "label_seq_id", "alignment_state", "rfam_pos", "has_coordinates", "base_complete"])
        for i in instances:
            for r in maps[i["seq_name"]]["residues"]:
                key = (i["pdb_id"].upper(), i["auth_asym_id"], r["author_residue_number"], r["author_insertion_code"] or "")
                pos = r["rfam_consensus_pos"] if r["alignment_state"] == "match" else (f"{r['insertion_anchor']}^{r['insertion_offset']}" if r["alignment_state"] == "insert" else "")
                w.writerow([i["seq_name"], r["polymer_seq_position"], r["alignment_state"], pos, r["has_coordinates"], bool(complete.get(key, (None, False))[1])])

    summary = {"rows_in_files": dict(rows_total), "endpoint_status": dict(status), "endpoint_verification_against_mmcif": dict(verify),
               "pair_kinds": {f"{a}|{b}": n for (a, b), n in pair_kinds.items()}, "reverse_order_rows_present": {f"{k}|{v}": n for (k, v), n in reverse_check.items()},
               "n_intra_chain_interactions_after_normalisation": len(out), "n_inter_chain_unique_rows": len(inter), "examples": examples,
               "n_instances": len(instances), "n_entries": len({i['pdb_id'] for i in instances})}
    json.dump(summary, open(OUT / "stage1_summary.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in summary.items() if k != "examples"}, indent=1))


if __name__ == "__main__":
    main()
