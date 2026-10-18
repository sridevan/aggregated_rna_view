"""Arpeggio pilot on deposited wwPDB archive mmCIF coordinates (asymmetric unit, model 1, no assembly or symmetry expansion).

Per PDB entry:
 1. target RNA chains = Rfam-mapped instances (family instances.json + mappings, or a POC numbering CSV for one chain);
 2. nonpolymer ligand instances read from the deposited mmCIF with gemmi (entity type non-polymer, water excluded),
    classified with rna_explorer/annotation_config.json (ligand / metal_ion / additive);
 3. pdbe-arpeggio run once per entry with every ligand instance selected (selection file, author ids);
 4. INTER contacts kept only when one side is a ligand instance and the other a residue of a target chain;
 5. atom-level records retained (chain, residue, atoms, distance, contact classes); residue-level and consensus-level aggregation;
 6. RNA residues mapped to (pdb_id, label_asym_id, label_seq_id) and Rfam consensus positions through the validated mapping;
 7. aggregation per RNA chain instance, with PDB-entry counts kept separately;
 8. ligand classes kept separate; 9. coordinate context recorded on every output.

Run with the Arpeggio venv python (needs gemmi):  <venv>/bin/python3 scripts/pilot.py --family RF00167 [--entries 1y26 3rkf]
                                               <venv>/bin/python3 scripts/pilot.py --entry 7k00 --chain A --mapping-csv ../../../poc_rf00177/results/ecoli_7k00_numbering_to_consensus.csv
"""
from __future__ import annotations
import argparse, csv, json, subprocess, sys, time
from collections import Counter, defaultdict
from pathlib import Path
import gemmi

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent
ROOT = PILOT.parents[2]
sys.path.insert(0, str(ROOT))
from rna_explorer.annotation_mapping import InstanceIndex, UnitId   # noqa: E402  (read-only use)
from rna_explorer.mapping import pp_category, pp_value              # noqa: E402

CONFIG = json.load(open(ROOT / "rna_explorer/annotation_config.json"))
WATER = {"HOH", "DOD", "H2O", "WAT"}
COORD_CONTEXT = {"coordinate_source": "wwPDB archive mmCIF, deposited coordinates (https://files.rcsb.org/download/{PDB}.cif)",
                 "model": 1, "assembly": "none: asymmetric unit as deposited", "symmetry_expansion": False,
                 "hydrogens": "added at run time by Arpeggio/OpenBabel; not deposited", "software": "pdbe-arpeggio 1.4.4, openbabel-wheel 3.1.1.23, gemmi 0.7.5",
                 "contact_definition": {"any": "Arpeggio INTER contact (includes 'proximal' atom pairs within the 5 A interacting cutoff)",
                                        "specific": "at least one non-proximal class (hbond, polar, weak_hbond, weak_polar, vdw, vdw_clash, ionic, metal_complex, aromatic, carbonyl, xbond, ...)"}}


def classify(ccd: str) -> str:
    for cls, ids in CONFIG["bound_molecule_classes"].items():
        if not cls.startswith("_") and ccd.upper() in {x.upper() for x in ids}:
            return cls
    return "ligand"


# ----------------------------------------------------------------------------- target chains and resolvers
class FamilyTargets:
    """Rfam-mapped chains of one entry from the production family build (read-only)."""
    def __init__(self, fam_dir: Path):
        self.idx = InstanceIndex(fam_dir)

    def chains(self, pdb: str):
        return [{"seq_name": r["seq_name"], "auth_asym_id": r["auth_asym_id"], "label_asym_id": r["label_asym_id"], "entity_id": r["entity_id"]} for r in self.idx.by_pdb.get(pdb, [])]

    def resolve(self, pdb: str, auth_chain: str, num: int, ins: str):
        ref = self.idx.resolve_unit(f"{pdb.upper()}|1|{auth_chain}|X|{num}" + (f"|||{ins}" if ins else ""))
        return {"seq_name": ref.seq_name, "label_asym_id": ref.label_asym_id, "label_seq_id": ref.label_seq_id, "residue_name": ref.residue_name,
                "alignment_state": ref.alignment_state, "rfam_consensus_pos": ref.rfam_consensus_pos, "pp_char": ref.pp_char, "pp_category": ref.pp_category,
                "has_coordinates": ref.has_coordinates, "mapping_status": ref.mapping_status}

    def observable(self, seq_name: str, pos: int) -> bool:
        return self.idx.observable(seq_name, pos)

    def clen(self) -> int:
        return self.idx.clen


class CsvTargets:
    """One chain mapped through a POC numbering CSV (author number -> consensus); label_seq_id from the deposited mmCIF."""
    def __init__(self, pdb: str, auth_chain: str, csv_path: Path, st: gemmi.Structure, clen: int):
        self.pdb, self.auth_chain, self._clen = pdb, auth_chain, clen
        self.rows = {int(r[next(k for k in r if k.startswith("ecoli_residue_number"))]): r for r in csv.DictReader(open(csv_path))}
        ch = st[0][auth_chain]
        self.label_asym = ch[0].subchain
        self.label_by_auth = {(res.seqid.num, res.seqid.icode.strip()): (res.label_seq, res.name) for res in ch}
        self.seq_name = f"{pdb}_?_{auth_chain}_{self.label_asym}"

    def chains(self, pdb: str):
        return [{"seq_name": self.seq_name, "auth_asym_id": self.auth_chain, "label_asym_id": self.label_asym, "entity_id": None}]

    def resolve(self, pdb: str, auth_chain: str, num: int, ins: str):
        if auth_chain != self.auth_chain:
            return {"seq_name": None, "label_asym_id": None, "label_seq_id": None, "residue_name": None, "alignment_state": None, "rfam_consensus_pos": None,
                    "pp_char": "", "pp_category": None, "has_coordinates": None, "mapping_status": "unresolved_identifier"}
        r = self.rows.get(num); lab = self.label_by_auth.get((num, ins))
        if r is None:
            return {"seq_name": self.seq_name, "label_asym_id": self.label_asym, "label_seq_id": lab[0] if lab else None, "residue_name": lab[1] if lab else None,
                    "alignment_state": None, "rfam_consensus_pos": None, "pp_char": "", "pp_category": None, "has_coordinates": lab is not None, "mapping_status": "outside_mapped_region"}
        st_ = r["alignment_state"]; cons = int(r["rfam_consensus_pos"]) if r["rfam_consensus_pos"] else None; pp = r["posterior_prob"] or ""
        status = "mapped" if st_ == "match" and lab else "mapped_no_coordinates" if st_ == "match" else st_
        return {"seq_name": self.seq_name, "label_asym_id": self.label_asym, "label_seq_id": lab[0] if lab else None, "residue_name": r["residue_name"],
                "alignment_state": st_, "rfam_consensus_pos": cons, "pp_char": pp, "pp_category": pp_category(pp) if pp else None, "has_coordinates": lab is not None, "mapping_status": status}

    def observable(self, seq_name: str, pos: int) -> bool:
        for num, r in self.rows.items():
            if r["alignment_state"] == "match" and r["rfam_consensus_pos"] and int(r["rfam_consensus_pos"]) == pos:
                return (num, "") in self.label_by_auth
        return False

    def clen(self) -> int:
        return self._clen


# ----------------------------------------------------------------------------- ligand inventory from the deposited mmCIF
def ligand_instances(st: gemmi.Structure) -> list[dict]:
    st.setup_entities()
    etype = {e.name: e.entity_type for e in st.entities}
    out = []
    for ch in st[0]:
        for res in ch:
            if res.name in WATER or res.het_flag != "H":
                continue
            ent = st.get_entity_of(ch.get_subchain(res.subchain)) if res.subchain else None
            if ent is not None and ent.entity_type != gemmi.EntityType.NonPolymer:
                continue   # branched or polymer het residues (modified nucleotides) are not ligand instances
            out.append({"ccd_id": res.name, "ligand_class": classify(res.name), "auth_asym_id": ch.name, "auth_seq_num": res.seqid.num, "ins_code": res.seqid.icode.strip(),
                        "label_asym_id": res.subchain, "entity_id": ent.name if ent else None, "n_atoms": len(res), "altlocs": sorted({a.altloc for a in res if a.altloc != "\0"})})
    return out


# ----------------------------------------------------------------------------- arpeggio
def run_arpeggio(arpeggio: Path, cif: Path, ligs: list[dict], out_dir: Path) -> tuple[Path, float]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{cif.stem}.json"
    sel = out_dir / "selection.txt"
    sel.write_text("".join(f"/{l['auth_asym_id']}/{l['auth_seq_num']}{l['ins_code']}/\n" for l in ligs))
    if out.exists():
        return out, 0.0
    t = time.time()
    r = subprocess.run([str(arpeggio), str(cif), "-sf", str(sel), "-o", str(out_dir)], capture_output=True, text=True)
    (out_dir / "arpeggio.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0:
        raise RuntimeError(f"arpeggio failed for {cif}: {r.stderr[-500:]}")
    return out, time.time() - t


# ----------------------------------------------------------------------------- per entry
def _exptl_method(cif: Path):
    try:
        v = gemmi.cif.read(str(cif)).sole_block().find_value("_exptl.method")
        return gemmi.cif.as_string(v) if v else None
    except Exception:
        return None


def process_entry(pdb: str, targets, cif: Path, arpeggio: Path, work: Path) -> dict:
    st = gemmi.read_structure(str(cif))
    ligs = ligand_instances(st)
    chains = targets.chains(pdb)
    target_auth = {c["auth_asym_id"]: c for c in chains}
    lig_by_key = {(l["auth_asym_id"], l["auth_seq_num"], l["ins_code"]): l for l in ligs}
    rec = {"pdb_id": pdb, "coordinate_context": {**COORD_CONTEXT, "file": cif.name, "n_models_in_file": len(st), "experimental_method": _exptl_method(cif),
                                                 "spacegroup": st.spacegroup_hm or None, "cell": [st.cell.a, st.cell.b, st.cell.c] if st.cell.a else None},
           "target_chains": chains, "ligand_instances": ligs, "arpeggio": {}, "atom_contacts": [], "residue_contacts": [], "per_chain": {}, "entry": {}}
    if not ligs:
        rec["arpeggio"] = {"status": "skipped_no_nonpolymer_ligands"}
        return rec
    out, secs = run_arpeggio(arpeggio, cif, ligs, work / pdb)
    data = json.load(open(out))
    rec["arpeggio"] = {"status": "ok", "seconds": round(secs, 2), "n_contacts_total": len(data), "selection_file": str((work / pdb / "selection.txt").relative_to(PILOT)),
                       "output": str(out.relative_to(PILOT))}
    # ---- atom-level filtering
    for c in data:
        if c["interacting_entities"] != "INTER":
            continue
        b, e = c["bgn"], c["end"]
        kb, ke = (b["auth_asym_id"], b["auth_seq_id"], b["pdbx_PDB_ins_code"].strip()), (e["auth_asym_id"], e["auth_seq_id"], e["pdbx_PDB_ins_code"].strip())
        if kb in lig_by_key and ke not in lig_by_key and e["auth_asym_id"] in target_auth and e["label_comp_id"] not in WATER:
            lig, rna, lk = b, e, kb
        elif ke in lig_by_key and kb not in lig_by_key and b["auth_asym_id"] in target_auth and b["label_comp_id"] not in WATER:
            lig, rna, lk = e, b, ke
        else:
            continue
        L = lig_by_key[lk]
        rec["atom_contacts"].append({"ligand_ccd": L["ccd_id"], "ligand_class": L["ligand_class"], "ligand_label_asym_id": L["label_asym_id"], "ligand_auth_asym_id": L["auth_asym_id"],
                                     "ligand_auth_seq_num": L["auth_seq_num"], "ligand_ins_code": L["ins_code"], "ligand_atom": lig["auth_atom_id"],
                                     "rna_auth_asym_id": rna["auth_asym_id"], "rna_auth_seq_num": rna["auth_seq_id"], "rna_ins_code": rna["pdbx_PDB_ins_code"].strip(),
                                     "rna_comp_id": rna["label_comp_id"], "rna_atom": rna["auth_atom_id"], "distance": c["distance"], "contact_classes": c["contact"], "type": c["type"]})
    # ---- residue level + mapping
    per_res: dict[tuple, dict] = {}
    for a in rec["atom_contacts"]:
        k = (a["ligand_label_asym_id"], a["rna_auth_asym_id"], a["rna_auth_seq_num"], a["rna_ins_code"])
        x = per_res.setdefault(k, {"ligand_ccd": a["ligand_ccd"], "ligand_class": a["ligand_class"], "ligand_label_asym_id": a["ligand_label_asym_id"], "ligand_auth_asym_id": a["ligand_auth_asym_id"],
                                   "ligand_auth_seq_num": a["ligand_auth_seq_num"], "rna_auth_asym_id": a["rna_auth_asym_id"], "rna_auth_seq_num": a["rna_auth_seq_num"], "rna_ins_code": a["rna_ins_code"],
                                   "rna_comp_id": a["rna_comp_id"], "n_atom_pairs": 0, "min_distance": None, "contact_classes": set(), "ligand_atoms": set(), "rna_atoms": set()})
        x["n_atom_pairs"] += 1; x["min_distance"] = a["distance"] if x["min_distance"] is None else min(x["min_distance"], a["distance"])
        x["contact_classes"].update(a["contact_classes"]); x["ligand_atoms"].add(a["ligand_atom"]); x["rna_atoms"].add(a["rna_atom"])
    for k, x in per_res.items():
        m = targets.resolve(pdb, x["rna_auth_asym_id"], x["rna_auth_seq_num"], x["rna_ins_code"])
        x.update(m); x["contact_classes"] = sorted(x["contact_classes"]); x["ligand_atoms"] = sorted(x["ligand_atoms"]); x["rna_atoms"] = sorted(x["rna_atoms"])
        x["specific"] = any(t != "proximal" for t in x["contact_classes"])
        rec["residue_contacts"].append(x)
    # ---- aggregation per RNA chain instance and per entry
    for c in chains:
        sn = c["seq_name"]
        mine = [x for x in rec["residue_contacts"] if x["seq_name"] == sn]
        pos: dict[int, dict] = {}
        for x in mine:
            if x["rfam_consensus_pos"] is None:
                continue
            p = pos.setdefault(x["rfam_consensus_pos"], {"ligand_instances": set(), "by_class": defaultdict(set), "by_ccd": defaultdict(set), "specific_ligand_instances": set(), "min_distance": None})
            lid = f"{x['ligand_ccd']}:{x['ligand_label_asym_id']}"
            p["ligand_instances"].add(lid); p["by_class"][x["ligand_class"]].add(lid); p["by_ccd"][x["ligand_ccd"]].add(lid)
            if x["specific"]:
                p["specific_ligand_instances"].add(lid)
            p["min_distance"] = x["min_distance"] if p["min_distance"] is None else min(p["min_distance"], x["min_distance"])
        rec["per_chain"][sn] = {"auth_asym_id": c["auth_asym_id"], "label_asym_id": c["label_asym_id"], "n_residue_contacts": len(mine),
                                "n_residue_contacts_on_consensus": sum(1 for x in mine if x["rfam_consensus_pos"] is not None),
                                "mapping_status_counts": dict(Counter(x["mapping_status"] for x in mine)),
                                "positions": {str(p): {"n_ligand_instances": len(v["ligand_instances"]), "ligand_instances": sorted(v["ligand_instances"]),
                                                       "by_class": {k: sorted(s) for k, s in v["by_class"].items()}, "by_ccd": {k: len(s) for k, s in v["by_ccd"].items()},
                                                       "n_specific_ligand_instances": len(v["specific_ligand_instances"]), "min_distance": v["min_distance"]} for p, v in sorted(pos.items())}}
    ent_pos: dict[int, dict] = defaultdict(lambda: {"chains": set(), "by_class": defaultdict(set), "by_ccd": defaultdict(set)})
    for x in rec["residue_contacts"]:
        if x["rfam_consensus_pos"] is not None:
            e = ent_pos[x["rfam_consensus_pos"]]; e["chains"].add(x["seq_name"]); e["by_class"][x["ligand_class"]].add(x["seq_name"]); e["by_ccd"][x["ligand_ccd"]].add(x["seq_name"])
    rec["entry"] = {"n_ligand_instances": len(ligs), "ligand_instances_by_class": dict(Counter(l["ligand_class"] for l in ligs)),
                    "n_target_chains": len(chains), "n_chains_with_any_contact": len({x["seq_name"] for x in rec["residue_contacts"]}),
                    "positions": {str(p): {"n_chains": len(v["chains"]), "by_class_n_chains": {k: len(s) for k, s in v["by_class"].items()}, "by_ccd_n_chains": {k: len(s) for k, s in v["by_ccd"].items()}}
                                  for p, v in sorted(ent_pos.items())}}
    return rec


# ----------------------------------------------------------------------------- family aggregation
def aggregate(records: list[dict], targets, out: Path) -> dict:
    clen = targets.clen()
    all_chains = [(r["pdb_id"], c["seq_name"]) for r in records for c in r["target_chains"]]
    evaluated = {c["seq_name"] for r in records if r["arpeggio"].get("status") == "ok" for c in r["target_chains"]}
    skipped = {c["seq_name"] for r in records if r["arpeggio"].get("status") != "ok" for c in r["target_chains"]}
    positions = []
    for pos in range(1, clen + 1):
        observable = {sn for _, sn in all_chains if targets.observable(sn, pos)}
        by_class_inst: dict[str, set] = defaultdict(set); by_class_entry: dict[str, set] = defaultdict(set); by_ccd_inst: dict[str, set] = defaultdict(set)
        spec_inst: dict[str, set] = defaultdict(set)
        for r in records:
            for sn, ch in r["per_chain"].items():
                p = ch["positions"].get(str(pos))
                if not p:
                    continue
                for cls, lids in p["by_class"].items():
                    by_class_inst[cls].add(sn); by_class_entry[cls].add(r["pdb_id"])
                    if p["n_specific_ligand_instances"]:
                        spec_inst[cls].add(sn)
                for ccd in p["by_ccd"]:
                    by_ccd_inst[ccd].add(sn)
        tracks = {}
        for cls in ("ligand", "metal_ion", "additive"):
            inst = by_class_inst.get(cls, set()); obs_eval = observable & evaluated
            tracks[cls] = {"n_instances": len(inst), "n_entries": len(by_class_entry.get(cls, set())), "n_instances_specific": len(spec_inst.get(cls, set())),
                           "denominators": {"counting_unit": "rna_chain_instance", "n_observable_instances": len(observable), "n_evaluated_instances": len(evaluated),
                                            "n_observable_and_evaluated": len(obs_eval), "n_positive_observable_and_evaluated": len(inst & obs_eval),
                                            "n_not_contacted": len(obs_eval - inst), "n_not_evaluated": len(observable - evaluated), "n_positive_not_observable": len(inst - observable)},
                           "frequency": round(len(inst & obs_eval) / len(obs_eval), 4) if obs_eval else None}
        positions.append({"rfam_consensus_pos": pos, "n_observable_instances": len(observable), "tracks": tracks, "by_ccd_n_instances": {k: len(v) for k, v in sorted(by_ccd_inst.items())}})
    agg = {"coordinate_context": COORD_CONTEXT, "n_entries": len(records), "n_chain_instances": len(all_chains), "n_instances_evaluated": len(evaluated), "n_instances_skipped": sorted(skipped),
           "entry_level_note": "entry counts are distinct PDB ids with >=1 contacting chain; chain counts are distinct Rfam-mapped chain instances; never summed",
           "positions": positions}
    json.dump(agg, open(out, "w"), indent=1)
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family"); ap.add_argument("--entries", nargs="*"); ap.add_argument("--entry"); ap.add_argument("--chain"); ap.add_argument("--mapping-csv", type=Path); ap.add_argument("--clen", type=int, default=1533)
    ap.add_argument("--arpeggio", type=Path, required=True); ap.add_argument("--out", type=Path, default=PILOT / "results"); ap.add_argument("--work", type=Path, default=PILOT / "work")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    records = []
    if a.family:
        targets = FamilyTargets(ROOT / "data/families" / a.family)
        pdbs = [p.lower() for p in a.entries] if a.entries else sorted(targets.idx.by_pdb)
        for pdb in pdbs:
            cif = PILOT / "cif" / f"{pdb.upper()}.cif"
            if not cif.exists():
                print("missing deposited mmCIF", cif); continue
            rec = process_entry(pdb, targets, cif, a.arpeggio, a.work)
            json.dump(rec, open(a.out / f"{pdb}.json", "w"), indent=1); records.append(rec)
            print(pdb, rec["arpeggio"].get("status"), rec["arpeggio"].get("seconds"), "ligands", len(rec["ligand_instances"]), "atom contacts", len(rec["atom_contacts"]), "residue contacts", len(rec["residue_contacts"]))
        aggregate(records, targets, a.out / f"{a.family}_chain_level_tracks.json")
    else:
        pdb = a.entry.lower(); cif = PILOT / "cif" / f"{pdb.upper()}.cif"
        st = gemmi.read_structure(str(cif))
        targets = CsvTargets(pdb, a.chain, a.mapping_csv, st, a.clen)
        rec = process_entry(pdb, targets, cif, a.arpeggio, a.work)
        json.dump(rec, open(a.out / f"{pdb}.json", "w"), indent=1)
        print(pdb, rec["arpeggio"], "ligands", len(rec["ligand_instances"]), "atom contacts", len(rec["atom_contacts"]), "residue contacts", len(rec["residue_contacts"]))
    # flat CSV of residue contacts
    rows = [{"pdb_id": r["pdb_id"], **{k: (";".join(v) if isinstance(v, list) else v) for k, v in x.items()}} for r in records or [rec] for x in r["residue_contacts"]]
    if rows:
        name = a.family or a.entry
        with open(a.out / f"{name}_residue_contacts.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
