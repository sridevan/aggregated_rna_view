"""Validate residue mapping CSV: continuity, coverage, mmCIF identifier agreement, mapped-vs-full consistency."""
import csv, json, shlex, sys
from pathlib import Path
from collections import defaultdict
DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
chains = json.load(open(RES/"chains.json"))

def cif_loop(path, category):
    """Return list of dicts for an mmCIF loop_ category (simple tokenizer; handles quoted tokens, not multi-line ;-blocks)."""
    rows, keys, inloop = [], [], False
    with open(path) as f:
        lines = iter(f)
        for line in lines:
            if line.startswith("loop_"): inloop, keys = True, []; continue
            if inloop and line.startswith(category + "."): keys.append(line.split()[0].split(".",1)[1]); continue
            if inloop and keys:
                if line.startswith("#") or line.startswith("_") or line.startswith("loop_"):
                    if rows: return rows
                    inloop, keys = line.startswith("loop_"), []; continue
                toks = shlex.split(line.strip(), posix=True)
                if toks: rows.append(dict(zip(keys, toks)))
    return rows

report = {"chains": {}, "global": {}}
mapped = list(csv.DictReader(open(RES/"mapped_residue_mapping.csv")))
full = list(csv.DictReader(open(RES/"full_residue_mapping.csv")))
msum = json.load(open(RES/"mapped_summary.json"))
scores = {l.split()[1]: l.split() for l in open(RES/"rf00167_mapped.scores") if not l.startswith("#")}

for ch in chains:
    key = f"{ch['pdb_id']}/{ch['auth_chain_id']}"
    rows = [r for r in mapped if r["pdb_id"] == ch["pdb_id"] and r["auth_chain_id"] == ch["auth_chain_id"]]
    rows.sort(key=lambda r: int(r["polymer_seq_position"]))
    out = {}
    # --- 1. mmCIF agreement (pdbx_poly_seq_scheme) ---
    scheme = [r for r in cif_loop(DATA/f"{ch['pdb_id']}_updated.cif", "_pdbx_poly_seq_scheme")
              if r["asym_id"] == ch["label_asym_id"] and r["entity_id"] == str(ch["entity_id"])]
    by_seq = {int(r["seq_id"]): r for r in scheme}
    mismatches = []
    for r in rows:
        p = int(r["polymer_seq_position"]); s = by_seq.get(p)
        if s is None: mismatches.append((p, "missing in mmCIF scheme")); continue
        checks = {"auth_chain": (s["pdb_strand_id"], r["auth_chain_id"]),
                  "auth_num": (s["pdb_seq_num"], r["author_residue_number"]),
                  "ins_code": ("" if s["pdb_ins_code"] in (".", "?") else s["pdb_ins_code"], r["author_insertion_code"]),
                  "mon_id": (s["mon_id"], r["residue_name"])}
        bad = {k: v for k, v in checks.items() if str(v[0]) != str(v[1])}
        if bad: mismatches.append((p, bad))
        # unobserved in model coordinates per mmCIF (auth_seq_num '?')
    unobs_cif = sorted(int(s["seq_id"]) for s in scheme if s["auth_seq_num"] in ("?", "."))
    unobs_api = sorted(int(r["polymer_seq_position"]) for r in rows if float(r["observed_ratio"] or 0) == 0)
    partial_api = sorted((int(r["polymer_seq_position"]), r["observed_ratio"]) for r in rows if 0 < float(r["observed_ratio"] or 0) < 1)
    modres = sorted((int(r["polymer_seq_position"]), r["residue_name"]) for r in rows if r["residue_name"] not in ("A","C","G","U"))
    out["mmcif_scheme_rows"] = len(scheme); out["chain_length_api"] = len(ch["sequence"]); out["mmcif_mismatches"] = mismatches
    out["unobserved_residues_mmcif"] = unobs_cif; out["unobserved_residues_api"] = unobs_api; out["partially_observed_api"] = partial_api
    out["modified_residues"] = modres
    out["mapped_region_contains_unobserved"] = [p for p in unobs_cif if ch["map_start"] <= p <= ch["map_end"]]
    out["mapped_region_contains_modified"] = [m for m in modres if ch["map_start"] <= m[0] <= ch["map_end"]]
    # --- 2. continuity ---
    m = [(int(r["polymer_seq_position"]), int(r["rfam_consensus_pos"])) for r in rows if r["status"] == "match"]
    jumps = [(a, b) for a, b in zip(m, m[1:]) if b[1] - a[1] != b[0] - a[0]]
    mono = all(b[1] > a[1] for a, b in zip(m, m[1:]))
    out["consensus_monotonic"] = mono
    out["offset_changes"] = [{"between_polymer": (a[0], b[0]), "consensus": (a[1], b[1]), "deleted_consensus_positions": b[1]-a[1]-(b[0]-a[0])} for a, b in jumps]
    out["n_insert"] = sum(r["status"] == "insert" for r in rows); out["n_local_end"] = sum(r["status"] == "local_end" for r in rows)
    out["n_unexpected"] = sum(r["status"] not in ("match","insert","local_end","not_submitted_flank") for r in rows)
    # --- 3. coverage ---
    seqname = ch["seq_name"] + f"/{ch['map_start']}-{ch['map_end']}"
    sm = msum["per_seq"][seqname]
    n_mapped = ch["map_end"] - ch["map_start"] + 1
    out["mapped_region_len"] = n_mapped; out["residues_in_match_cols"] = sm["n_match"]
    out["fraction_mapped_residues_in_consensus"] = round(sm["n_match"] / n_mapped, 4)
    out["consensus_span"] = [sm["cons_min"], sm["cons_max"]]; out["consensus_positions_covered"] = sm["n_match"]
    out["fraction_of_CLEN_covered"] = round(sm["n_match"] / msum["CLEN"], 4)
    out["deleted_consensus_positions"] = sm["deleted_cons_positions"]
    sc = scores[seqname]
    out["cmalign"] = {"cm_from": int(sc[3]), "cm_to": int(sc[4]), "trunc": sc[5], "bit_score": float(sc[6]), "avg_pp": float(sc[7])}
    out["cmalign_span_equals_parsed_span"] = (int(sc[3]), int(sc[4])) == (sm["cons_min"], sm["cons_max"])
    low_pp = [(int(r["polymer_seq_position"]), r["posterior_prob"]) for r in rows if r["posterior_prob"] and r["posterior_prob"] not in ("*",) and r["posterior_prob"] < "9"]
    out["residues_with_pp_below_0.85"] = low_pp
    # --- 4. mapped vs full alignment consistency ---
    fr = {int(r["polymer_seq_position"]): r for r in full if r["pdb_id"] == ch["pdb_id"] and r["auth_chain_id"] == ch["auth_chain_id"]}
    diffs = []
    for r in rows:
        if r["status"] == "not_submitted_flank": continue
        p = int(r["polymer_seq_position"]); f = fr[p]
        if (r["rfam_consensus_pos"], r["status"]) != (f["rfam_consensus_pos"], f["status"]):
            diffs.append({"polymer_pos": p, "nt": r["nucleotide"], "mapped_aln": (r["status"], r["rfam_consensus_pos"]), "full_aln": (f["status"], f["rfam_consensus_pos"])})
    out["residues_assigned_differently_when_full_chain_aligned"] = len(diffs)
    out["examples_of_differences"] = diffs[:6]
    # --- 5. species provenance ---
    out["species"] = {"api": (ch["api_organism"], ch["api_tax_id"]), "mmcif": (ch.get("cif_organism"), ch.get("cif_tax_id")), "agree": str(ch["api_organism"]) == str(ch.get("cif_organism"))}
    out["label_vs_auth_chain_differ"] = ch["label_asym_id"] != ch["auth_chain_id"]
    out["author_numbering_offset"] = ch["map_auth_start"] - ch["map_start"]
    report["chains"][key] = out

# global: per consensus position, how many chains have a residue (the aggregated view)
cov = defaultdict(int)
for r in mapped:
    if r["status"] == "match": cov[int(r["rfam_consensus_pos"])] += 1
report["global"]["consensus_positions_with_any_residue"] = len(cov)
report["global"]["consensus_positions_with_all_5_chains"] = sum(1 for v in cov.values() if v == 5)
report["global"]["consensus_positions_never_covered"] = [i for i in range(1, 103) if i not in cov]
report["global"]["total_mmcif_mismatches"] = sum(len(v["mmcif_mismatches"]) for v in report["chains"].values())
json.dump(report, open(RES/"validation_report.json", "w"), indent=1)
print(json.dumps(report, indent=1))
