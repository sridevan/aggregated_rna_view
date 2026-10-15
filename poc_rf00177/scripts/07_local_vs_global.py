"""Stability test: compare consensus assignment of every residue between local (default) and global (-g) cmalign modes."""
import csv, json, sys, collections
loc = {(r["pdb_id"], int(r["polymer_seq_position"])): r for r in csv.DictReader(open(sys.argv[1]))}
glo = {(r["pdb_id"], int(r["polymer_seq_position"])): r for r in csv.DictReader(open(sys.argv[2]))}
out = {}
for pdb in sorted({k[0] for k in loc}):
    keys = sorted(k for k in loc if k[0] == pdb and loc[k]["alignment_state"] != "not_submitted_flank")
    same_match = diff_match = loc_match_only = glo_match_only = both_nonmatch = 0; examples = []
    for k in keys:
        a, b = loc[k], glo[k]
        if a["alignment_state"] == "match" and b["alignment_state"] == "match":
            if a["rfam_consensus_pos"] == b["rfam_consensus_pos"]: same_match += 1
            else:
                diff_match += 1
                if len(examples) < 5: examples.append((k[1], a["rfam_consensus_pos"], b["rfam_consensus_pos"]))
        elif a["alignment_state"] == "match": loc_match_only += 1
        elif b["alignment_state"] == "match": glo_match_only += 1
        else: both_nonmatch += 1
    out[pdb] = {"n": len(keys), "match_same_consensus": same_match, "match_different_consensus": diff_match, "match_in_local_only": loc_match_only,
                "match_in_global_only": glo_match_only, "nonmatch_in_both": both_nonmatch, "examples_different": examples,
                "local": {"insert": sum(loc[k]["alignment_state"] == "insert" for k in keys), "local_end": sum(loc[k]["alignment_state"] == "local_end" for k in keys)},
                "global": {"insert": sum(glo[k]["alignment_state"] == "insert" for k in keys), "local_end": sum(glo[k]["alignment_state"] == "local_end" for k in keys)}}
json.dump(out, open(sys.argv[3], "w"), indent=1)
for k, v in out.items(): print(k, json.dumps(v))
