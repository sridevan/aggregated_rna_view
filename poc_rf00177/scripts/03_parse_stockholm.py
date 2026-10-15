"""Step 4/5: parse cmalign Stockholm -> per-residue mapping with explicit insertion anchor/offset, per-consensus-column table.

RF line: letter = consensus match column (numbered 1..CLEN); '.' = insert column; '~' = local-end (EL) column.
Residue states: match | insert | local_end ; insertions carry (anchor = preceding consensus position, 0 if before first; offset = 1..n).
Observed status comes from atom_site (coordinates), independently of alignment state.
"""
import csv, json, sys
from pathlib import Path
from collections import defaultdict
sys.path.insert(0, str(Path(__file__).parent)); from common import atom_site_observed
sto, chains_json, cif_dir, out_prefix = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
FAM = "RF00177"; GAPS = set("-._~")

def parse_stockholm(path):
    seqs, gr, gc, order = defaultdict(str), defaultdict(lambda: defaultdict(str)), defaultdict(str), []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("# STOCKHOLM") or line == "//": continue
        if line.startswith("#=GC"): _, tag, txt = line.split(None, 2); gc[tag] += txt
        elif line.startswith("#=GR"): _, name, tag, txt = line.split(None, 3); gr[name][tag] += txt
        elif line.startswith("#"): continue
        else:
            name, txt = line.split(None, 1)
            if name not in seqs: order.append(name)
            seqs[name] += txt.strip()
    return order, seqs, gr, gc

order, seqs, gr, gc = parse_stockholm(sto)
RF, SS = gc["RF"], gc["SS_cons"]; L = len(RF)
assert all(len(s) == L for s in seqs.values()) and len(SS) == L
col_type, col_cons, prev_cons = [], [], []; n = 0
for ch in RF:
    if ch == ".": col_type.append("insert"); col_cons.append(None)
    elif ch == "~": col_type.append("local_end"); col_cons.append(None)
    else: n += 1; col_type.append("match"); col_cons.append(n)
    prev_cons.append(n)   # consensus position preceding (or at) this column; 0 before first
CLEN = n
chains = {c["seq_name"]: c for c in json.load(open(chains_json))}
rows, cons_rows, summary = [], [], {}
for name in order:
    base, rng = name.rsplit("/", 1); s, e = map(int, rng.split("-")); ch = chains[base]
    res_by_num = {r["residue_number"]: r for r in ch["residues"]}
    observed = atom_site_observed(cif_dir/f"{ch['pdb_id']}_updated.cif", ch["label_asym_id"], ch["entity_id"])
    aln, pp = seqs[name], gr[name].get("PP", "")
    pos, ins_off, cons_seen, ins_after = s - 1, 0, {}, defaultdict(int)
    common = {"pdb_id": ch["pdb_id"], "entity_id": ch["entity_id"], "label_asym_id": ch["label_asym_id"], "auth_asym_id": ch["auth_chain_id"],
              "rfam_acc": FAM, "seq_name": name, "organism_api": ch["api_organism"], "organism_mmcif": ch.get("cif_organism"), "tax_id_mmcif": ch.get("cif_tax_id")}
    def resrow(p, **kw):
        r = res_by_num[p]
        d = dict(common, polymer_seq_position=p, author_residue_number=r["author_residue_number"], author_insertion_code=r["author_insertion_code"] or "",
                 residue_name=r["residue_name"], nucleotide_parent=ch["sequence"][p-1], is_modified=int(r["residue_name"] not in ("A","C","G","U")),
                 has_coordinates=int(p in observed), observed_ratio_api=r.get("observed_ratio"), in_pdbe_mapped_region=int(ch["map_start"] <= p <= ch["map_end"]))
        d.update(kw); return d
    for col, a in enumerate(aln):
        ctype = col_type[col]
        if a in GAPS:
            if ctype == "match":
                ins_off = 0
                cons_rows.append({"pdb_id": ch["pdb_id"], "auth_asym_id": ch["auth_chain_id"], "seq_name": name, "rfam_consensus_pos": col_cons[col], "aln_column": col+1,
                                  "state": "deletion", "polymer_seq_position": "", "nucleotide_parent": "", "residue_name": "", "has_coordinates": "", "posterior_prob": "", "ss_cons": SS[col], "rf": RF[col]})
            continue
        pos += 1
        if ctype == "match":
            ins_off = 0; cpos = col_cons[col]; cons_seen[cpos] = pos
            rows.append(resrow(pos, aligned_char=a, alignment_state="match", rfam_consensus_pos=cpos, insertion_anchor="", insertion_offset="",
                               aln_column=col+1, rf=RF[col], ss_cons=SS[col], posterior_prob=pp[col] if pp else ""))
            cons_rows.append({"pdb_id": ch["pdb_id"], "auth_asym_id": ch["auth_chain_id"], "seq_name": name, "rfam_consensus_pos": cpos, "aln_column": col+1, "state": "match",
                              "polymer_seq_position": pos, "nucleotide_parent": ch["sequence"][pos-1], "residue_name": res_by_num[pos]["residue_name"], "has_coordinates": int(pos in observed),
                              "posterior_prob": pp[col] if pp else "", "ss_cons": SS[col], "rf": RF[col]})
        else:
            ins_off += 1; anchor = prev_cons[col]; ins_after[anchor] += 1
            rows.append(resrow(pos, aligned_char=a, alignment_state=ctype, rfam_consensus_pos="", insertion_anchor=anchor, insertion_offset=ins_off,
                               aln_column=col+1, rf=RF[col], ss_cons=SS[col], posterior_prob=pp[col] if pp else ""))
    assert pos == e, (name, pos, e)
    for rn in range(1, len(ch["sequence"]) + 1):   # chain residues outside the aligned (PDBe-mapped) region
        if rn < s or rn > e:
            rows.append(resrow(rn, aligned_char="", alignment_state="not_submitted_flank", rfam_consensus_pos="", insertion_anchor="", insertion_offset="",
                               aln_column="", rf="", ss_cons="", posterior_prob=""))
    lo, hi = min(cons_seen), max(cons_seen)   # consensus positions outside the matched span are 'uncovered' (truncation), inside are 'deletion'
    for cr in cons_rows:
        if cr["seq_name"] == name and cr["state"] == "deletion" and not (lo <= int(cr["rfam_consensus_pos"]) <= hi): cr["state"] = "uncovered"
    my = [r for r in rows if r["seq_name"] == name and r["alignment_state"] != "not_submitted_flank"]
    ppv = [r["posterior_prob"] for r in my if r["posterior_prob"]]
    summary[name] = {"mapped_region_len": e - s + 1, "n_match": len(cons_seen), "n_insert": sum(r["alignment_state"] == "insert" for r in my),
                     "n_local_end": sum(r["alignment_state"] == "local_end" for r in my), "cons_span": [min(cons_seen), max(cons_seen)],
                     "n_deleted_in_span": sum(1 for c in range(min(cons_seen), max(cons_seen)+1) if c not in cons_seen),
                     "n_uncovered_outside_span": CLEN - (max(cons_seen) - min(cons_seen) + 1), "fraction_CLEN_matched": round(len(cons_seen)/CLEN, 4),
                     "n_insertion_sites": len(ins_after), "longest_insertion": max(ins_after.values()) if ins_after else 0,
                     "longest_insertion_anchor": max(ins_after, key=ins_after.get) if ins_after else None,
                     "pp_star_fraction": round(sum(p == "*" for p in ppv)/len(ppv), 4) if ppv else None,
                     "pp_below_0.7_count": sum(p not in ("*", "9", "8", "7") for p in ppv),
                     "n_mapped_with_coordinates": sum(r["has_coordinates"] for r in my), "n_match_with_coordinates": sum(r["has_coordinates"] for r in my if r["alignment_state"] == "match"),
                     "n_modified_in_mapped": sum(r["is_modified"] for r in my)}
rows.sort(key=lambda r: (r["pdb_id"], r["auth_asym_id"], r["polymer_seq_position"]))
with open(out_prefix + "_residue_mapping.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
cons_rows.sort(key=lambda r: (r["pdb_id"], r["auth_asym_id"], int(r["rfam_consensus_pos"])))
with open(out_prefix + "_consensus_columns.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(cons_rows[0].keys())); w.writeheader(); w.writerows(cons_rows)
json.dump({"CLEN": CLEN, "alignment_length": L, "n_insert_columns": col_type.count("insert"), "n_local_end_columns": col_type.count("local_end"), "per_seq": summary}, open(out_prefix + "_summary.json", "w"), indent=1)
print(f"CLEN={CLEN} aln_len={L} insert_cols={col_type.count('insert')} EL_cols={col_type.count('local_end')}")
for k, v in summary.items(): print(" ", k, json.dumps(v))
