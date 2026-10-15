"""Parse a cmalign Stockholm file and emit a per-residue Rfam consensus mapping CSV.

RF line semantics (Infernal 1.1 user guide / cmalign man page):
  letter  -> consensus (match) column, numbered 1..CLEN
  '.'     -> insert column (residues here are insertions relative to consensus)
  '~'     -> local-end (EL) column (large local insertion; residues here are unaligned to the model)
Sequence-row gap characters: '-', '.', '_', '~' (all treated as gaps).
"""
import csv, json, re, sys
from pathlib import Path
from collections import defaultdict

sto_path, chains_json, out_prefix = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
GAPS = set("-._~")

def parse_stockholm(path):
    seqs, gr, gc = defaultdict(str), defaultdict(lambda: defaultdict(str)), defaultdict(str)
    order = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("# STOCKHOLM") or line == "//": continue
        if line.startswith("#=GC"):
            _, tag, txt = line.split(None, 2); gc[tag] += txt
        elif line.startswith("#=GR"):
            _, name, tag, txt = line.split(None, 3); gr[name][tag] += txt
        elif line.startswith("#"):
            continue
        else:
            name, txt = line.split(None, 1)
            if name not in seqs: order.append(name)
            seqs[name] += txt.strip()
    return order, seqs, gr, gc

order, seqs, gr, gc = parse_stockholm(sto_path)
RF, SS = gc["RF"], gc["SS_cons"]
L = len(RF)
assert all(len(s) == L for s in seqs.values()), "ragged alignment"
# column -> consensus position (1-based) or None; and type
col_type, col_cons = [], []
n = 0
for ch in RF:
    if ch == ".": col_type.append("insert"); col_cons.append(None)
    elif ch == "~": col_type.append("local_end"); col_cons.append(None)
    else: n += 1; col_type.append("match"); col_cons.append(n)
CLEN = n
# for insert / EL columns, which consensus position precedes them
prev_cons = []; last = 0
for c in col_cons:
    if c is not None: last = c
    prev_cons.append(last)

chains = {c["seq_name"]: c for c in json.load(open(chains_json))}
rows, cons_rows, summary = [], [], {}
for name in order:
    base, rng = name.rsplit("/", 1)
    s, e = map(int, rng.split("-"))
    ch = chains[base]
    res_by_num = {r["residue_number"]: r for r in ch["residues"]}
    aln, pp = seqs[name], gr[name].get("PP", "")
    pos = s - 1
    ins_counter = 0; per_res = {}
    cons_seen = {}
    for col, a in enumerate(aln):
        if a in GAPS:
            if col_type[col] == "match":
                cons_rows.append({"seq_name": name, "pdb_id": ch["pdb_id"], "auth_chain_id": ch["auth_chain_id"], "rfam_consensus_pos": col_cons[col],
                                  "aln_column": col+1, "status": "deletion", "polymer_seq_position": "", "nucleotide": "", "ss_cons": SS[col], "rf": RF[col]})
            if col_type[col] == "match": ins_counter = 0
            continue
        pos += 1
        r = res_by_num[pos]
        ctype = col_type[col]
        if ctype == "match":
            ins_counter = 0
            status, cpos, ins_info = "match", col_cons[col], ""
            cons_seen[cpos] = pos
            if a.islower(): status = "match_lowercase_unexpected"
        else:
            ins_counter += 1
            status, cpos = ctype, ""
            ins_info = f"{ctype}:after_cons_{prev_cons[col]}:#{ins_counter}"
        rows.append({"pdb_id": ch["pdb_id"], "entity_id": ch["entity_id"], "label_asym_id": ch["label_asym_id"], "auth_chain_id": ch["auth_chain_id"],
                     "polymer_seq_position": pos, "author_residue_number": r["author_residue_number"], "author_insertion_code": r["author_insertion_code"] or "",
                     "residue_name": r["residue_name"], "nucleotide": ch["sequence"][pos-1], "aligned_char": a,
                     "rfam_consensus_pos": cpos, "status": status, "insertion_info": ins_info,
                     "aln_column": col+1, "rf": RF[col], "ss_cons": SS[col], "posterior_prob": pp[col] if pp else "",
                     "observed_ratio": r.get("observed_ratio"), "in_pdbe_mapped_region": int(ch["map_start"] <= pos <= ch["map_end"]),
                     "organism_api": ch["api_organism"], "tax_id_api": ch["api_tax_id"], "organism_mmcif": ch.get("cif_organism"), "tax_id_mmcif": ch.get("cif_tax_id"),
                     "rfam_acc": "RF00167", "seq_name": name})
        if ctype == "match":
            cons_rows.append({"seq_name": name, "pdb_id": ch["pdb_id"], "auth_chain_id": ch["auth_chain_id"], "rfam_consensus_pos": cpos, "aln_column": col+1,
                              "status": "match", "polymer_seq_position": pos, "nucleotide": ch["sequence"][pos-1], "ss_cons": SS[col], "rf": RF[col]})
    assert pos == e, (name, pos, e)
    # residues of the chain not in the aligned range (flanks not submitted to cmalign)
    for rn in range(1, len(ch["sequence"]) + 1):
        if rn < s or rn > e:
            r = res_by_num[rn]
            rows.append({"pdb_id": ch["pdb_id"], "entity_id": ch["entity_id"], "label_asym_id": ch["label_asym_id"], "auth_chain_id": ch["auth_chain_id"],
                         "polymer_seq_position": rn, "author_residue_number": r["author_residue_number"], "author_insertion_code": r["author_insertion_code"] or "",
                         "residue_name": r["residue_name"], "nucleotide": ch["sequence"][rn-1], "aligned_char": "", "rfam_consensus_pos": "", "status": "not_submitted_flank",
                         "insertion_info": "", "aln_column": "", "rf": "", "ss_cons": "", "posterior_prob": "", "observed_ratio": r.get("observed_ratio"),
                         "in_pdbe_mapped_region": int(ch["map_start"] <= rn <= ch["map_end"]), "organism_api": ch["api_organism"], "tax_id_api": ch["api_tax_id"],
                         "organism_mmcif": ch.get("cif_organism"), "tax_id_mmcif": ch.get("cif_tax_id"), "rfam_acc": "RF00167", "seq_name": name})
    matched = [c for c in cons_seen]
    summary[name] = {"n_residues_aligned": e - s + 1, "n_match": len(cons_seen), "cons_min": min(matched), "cons_max": max(matched),
                     "n_insert": sum(1 for r in rows if r["seq_name"] == name and r["status"] == "insert"),
                     "n_local_end": sum(1 for r in rows if r["seq_name"] == name and r["status"] == "local_end"),
                     "n_deletions_within_span": sum(1 for c in range(min(matched), max(matched)+1) if c not in cons_seen),
                     "deleted_cons_positions": [c for c in range(min(matched), max(matched)+1) if c not in cons_seen]}

rows.sort(key=lambda r: (r["pdb_id"], r["auth_chain_id"], r["polymer_seq_position"]))
fields = list(rows[0].keys())
with open(out_prefix + "_residue_mapping.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
cons_rows.sort(key=lambda r: (r["pdb_id"], r["auth_chain_id"], r["rfam_consensus_pos"]))
with open(out_prefix + "_consensus_columns.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(cons_rows[0].keys())); w.writeheader(); w.writerows(cons_rows)
json.dump({"CLEN": CLEN, "alignment_length": L, "n_insert_columns": col_type.count("insert"), "n_local_end_columns": col_type.count("local_end"), "per_seq": summary},
          open(out_prefix + "_summary.json", "w"), indent=1)
print(f"{sto_path.name}: CLEN={CLEN} aln_len={L} insert_cols={col_type.count('insert')} EL_cols={col_type.count('local_end')}")
for k, v in summary.items(): print(" ", k, v)
