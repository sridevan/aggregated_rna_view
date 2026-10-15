"""Step 6: validate residue identifiers against mmCIF _pdbx_poly_seq_scheme, coordinates via atom_site, continuity, insertion bookkeeping,
and agreement with Rfam's own hit coordinates (structures endpoint)."""
import csv, json, sys, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent)); from common import cif_category, atom_site_observed
DATA, INT, RES = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
chains = json.load(open(INT/"chains.json"))
rows_all = list(csv.DictReader(open(RES/"rf00177_residue_mapping.csv")))
cons_all = list(csv.DictReader(open(RES/"rf00177_consensus_columns.csv")))
summ = json.load(open(RES/"rf00177_summary.json"))
scores = {l.split()[1]: l.split() for l in open(INT/"rf00177_mapped.scores") if not l.startswith("#")}
rfam = {(m["pdb_id"].lower(), m["chain"]): m for m in json.load(open(DATA/"api/rf00177_rfam_structures.json"))["mapping"]}
report = {"chains": {}, "global": {}}
for ch in chains:
    key = f"{ch['pdb_id']}/{ch['auth_chain_id']}"; cif = DATA/f"cif/{ch['pdb_id']}_updated.cif"
    rows = sorted([r for r in rows_all if r["pdb_id"] == ch["pdb_id"] and r["auth_asym_id"] == ch["auth_chain_id"]], key=lambda r: int(r["polymer_seq_position"]))
    scheme = {int(r["seq_id"]): r for r in cif_category(cif, "_pdbx_poly_seq_scheme") if r["asym_id"] == ch["label_asym_id"] and r["entity_id"] == str(ch["entity_id"])}
    observed = atom_site_observed(cif, ch["label_asym_id"], ch["entity_id"])
    mism, obs_conflicts = [], []
    for r in rows:
        p = int(r["polymer_seq_position"]); s = scheme.get(p)
        if s is None: mism.append({"pos": p, "problem": "absent from _pdbx_poly_seq_scheme"}); continue
        chk = {"auth_asym_id": (s["pdb_strand_id"], r["auth_asym_id"]), "label_asym_id": (s["asym_id"], r["label_asym_id"]), "entity_id": (s["entity_id"], r["entity_id"]),
               "author_residue_number": (s["pdb_seq_num"], r["author_residue_number"]),
               "author_insertion_code": ("" if s["pdb_ins_code"] in (".", "?") else s["pdb_ins_code"], r["author_insertion_code"]), "residue_name(mon_id)": (s["mon_id"], r["residue_name"])}
        bad = {k: v for k, v in chk.items() if str(v[0]) != str(v[1])}
        if bad: mism.append({"pos": p, "mismatch": bad})
        cif_obs = s["auth_seq_num"] not in ("?", ".")   # scheme says modelled
        if cif_obs != (p in observed) or (p in observed) != bool(int(r["has_coordinates"])) or ((float(r["observed_ratio_api"] or 0) > 0) != (p in observed)):
            obs_conflicts.append({"pos": p, "scheme_modelled": cif_obs, "atom_site": p in observed, "csv": r["has_coordinates"], "api_ratio": r["observed_ratio_api"]})
    # continuity + insertion bookkeeping
    m = [(int(r["polymer_seq_position"]), int(r["rfam_consensus_pos"])) for r in rows if r["alignment_state"] == "match"]
    mono = all(b[1] > a[1] for a, b in zip(m, m[1:]))
    ins_ok, last_anchor, last_off, anchors_ok = True, None, 0, True; prev_cons = 0; anchor_after_deletion = 0
    for r in rows:
        st = r["alignment_state"]
        if st == "match": prev_cons = int(r["rfam_consensus_pos"]); last_anchor, last_off = None, 0
        elif st in ("insert", "local_end"):
            a, o = int(r["insertion_anchor"]), int(r["insertion_offset"])
            if not (prev_cons <= a): anchors_ok = False
            if a > prev_cons and (last_anchor != a): anchor_after_deletion += 1   # anchor is the model position preceding the insert column; may exceed the last *matched* residue when deletions intervene
            if (last_anchor == a and o != last_off + 1) or (last_anchor != a and o != 1): ins_ok = False
            last_anchor, last_off = a, o
            if r["rfam_consensus_pos"] != "": ins_ok = False
    # no residue of the mapped region left without a state; every consensus position accounted for exactly once
    cons = [c for c in cons_all if c["pdb_id"] == ch["pdb_id"] and c["auth_asym_id"] == ch["auth_chain_id"]]
    cpos = collections.Counter(int(c["rfam_consensus_pos"]) for c in cons)
    cons_complete = sorted(cpos) == list(range(1, summ["CLEN"] + 1)) and max(cpos.values()) == 1
    name = f"{ch['seq_name']}/{ch['map_start']}-{ch['map_end']}"; sc = scores[name]; sm = summ["per_seq"][name]; rf = rfam.get((ch["pdb_id"], ch["auth_chain_id"]), {})
    unobs_mapped = [int(r["polymer_seq_position"]) for r in rows if r["in_pdbe_mapped_region"] == "1" and r["has_coordinates"] == "0"]
    def runs(xs):
        out = []
        for x in xs:
            if out and x == out[-1][1] + 1: out[-1][1] = x
            else: out.append([x, x])
        return out
    report["chains"][key] = {
        "entity_id": ch["entity_id"], "label_asym_id": ch["label_asym_id"], "chain_length": len(ch["sequence"]), "scheme_rows": len(scheme),
        "mmcif_mismatches": mism, "n_mmcif_mismatches": len(mism), "observed_status_conflicts": obs_conflicts[:5], "n_observed_conflicts": len(obs_conflicts),
        "n_modelled_residues_atom_site": len(observed), "mapped_residues_without_coordinates": len(unobs_mapped), "unmodelled_runs_in_mapped_region": runs(unobs_mapped),
        "consensus_monotonic": mono, "insertion_offsets_consistent": ins_ok, "insertion_anchors_equal_preceding_consensus": anchors_ok,
        "consensus_columns_table_complete_and_unique": cons_complete, "insertion_sites_whose_anchor_follows_deleted_consensus": anchor_after_deletion,
        "states": dict(collections.Counter(r["alignment_state"] for r in rows)),
        "cmalign": {"cm_from": int(sc[3]), "cm_to": int(sc[4]), "trunc": sc[5], "bits": float(sc[6]), "avg_pp": float(sc[7])},
        "rfam_structures_endpoint": {"cm_start": rf.get("cm_start"), "cm_end": rf.get("cm_end"), "bits": rf.get("bit_score"), "pdb_start": rf.get("pdb_start"), "pdb_end": rf.get("pdb_end")},
        "cm_span_agrees_with_rfam": (int(sc[3]), int(sc[4])) == (rf.get("cm_start"), rf.get("cm_end")),
        "pdbe_region_agrees_with_rfam": (ch["map_start"], ch["map_end"]) == (rf.get("pdb_start"), rf.get("pdb_end")),
        "modified_in_mapped_region": {k: v["three_letter_code"] for k, v in ch["modified_positions_in_mapped_region"].items()},
        "author_numbering_offset_at_map_start": ch["map_auth_start"] - ch["map_start"],
        "author_numbering_continuous": all(int(rows[i+1]["author_residue_number"]) - int(rows[i]["author_residue_number"]) == 1 for i in range(len(rows)-1)),
        "species": {"api": ch["api_organism"], "mmcif": ch["cif_organism"], "mmcif_tax_id": ch["cif_tax_id"]},
        "label_vs_auth_chain_differ": ch["label_asym_id"] != ch["auth_chain_id"]}
cov = collections.Counter(int(c["rfam_consensus_pos"]) for c in cons_all if c["state"] == "match")
report["global"] = {"n_chains": len(chains), "residue_rows": len(rows_all), "total_mmcif_mismatches": sum(v["n_mmcif_mismatches"] for v in report["chains"].values()),
                    "total_observed_conflicts": sum(v["n_observed_conflicts"] for v in report["chains"].values()),
                    "consensus_positions_matched_in_all_8": sum(1 for v in cov.values() if v == 8), "consensus_positions_matched_in_none": [c for c in range(1, summ["CLEN"]+1) if c not in cov],
                    "consensus_positions_matched_only_by_bacteria_and_chloroplast": sum(1 for c in range(1, summ["CLEN"]+1) if cov.get(c, 0) == 3)}
json.dump(report, open(RES/"validation_report.json", "w"), indent=1)
for k, v in report["chains"].items():
    print(k, json.dumps({x: v[x] for x in ["n_mmcif_mismatches", "n_observed_conflicts", "n_modelled_residues_atom_site", "mapped_residues_without_coordinates", "unmodelled_runs_in_mapped_region", "consensus_monotonic", "insertion_offsets_consistent", "insertion_anchors_equal_preceding_consensus", "consensus_columns_table_complete_and_unique", "cm_span_agrees_with_rfam", "pdbe_region_agrees_with_rfam", "author_numbering_offset_at_map_start", "author_numbering_continuous", "label_vs_auth_chain_differ"]}))
    if v["mmcif_mismatches"]: print("   MISMATCH", v["mmcif_mismatches"][:3])
    if v["observed_status_conflicts"]: print("   OBS CONFLICT", v["observed_status_conflicts"][:3])
print("global", json.dumps(report["global"]))
