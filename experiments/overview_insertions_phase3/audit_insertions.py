"""Phase 3 audit: how the RF00167 per-chain mappings represent insertions, deletions, flanks and missing coordinates (read-only)."""
import json, glob, collections, sys
from datetime import datetime, timezone
from pathlib import Path

ACC = sys.argv[1] if len(sys.argv) > 1 else "RF00167"
root = Path(__file__).resolve().parents[2]
files = sorted(glob.glob(str(root / f"data/families/{ACC}/mappings/*.json")))
states, columns = collections.Counter(), collections.Counter()
runs, flanks, checks = [], collections.Counter(), collections.Counter()
for f in files:
    d = json.load(open(f)); sn = Path(f).stem
    R = sorted(d["residues"], key=lambda r: r["polymer_seq_position"])
    states.update(r["alignment_state"] for r in R); columns.update(c["state"] for c in d["consensus_columns"])
    clen = d["provenance"]["cm_clen"]
    for r in R:
        if r["alignment_state"] == "not_submitted_flank":
            flanks["n"] += 1; flanks["with_coordinates"] += r["has_coordinates"]
            checks["flank_has_no_anchor"] += r["insertion_anchor"] is None and r["rfam_consensus_pos"] is None
        if r["alignment_state"] in ("insert", "local_end"):
            checks["insert_has_no_consensus_pos"] += r["rfam_consensus_pos"] is None
    cur = None
    for r in R:
        if r["alignment_state"] not in ("insert", "local_end"):
            cur = None; continue
        last = cur["residues"][-1] if cur else None
        if cur and cur["state"] == r["alignment_state"] and cur["anchor"] == r["insertion_anchor"] and r["polymer_seq_position"] == last["label_seq_id"] + 1 and r["insertion_offset"] == last["offset"] + 1:
            cur["residues"].append({"label_seq_id": r["polymer_seq_position"], "author": f'{r["author_residue_number"]}{r["author_insertion_code"]}', "nt": r["nucleotide_parent"], "offset": r["insertion_offset"], "has_coordinates": r["has_coordinates"], "pp": r["posterior_prob_char"]})
        else:
            a = r["insertion_anchor"]
            cur = {"seq_name": sn, "state": r["alignment_state"], "anchor": a, "left": a if a >= 1 else None, "right": a + 1 if a < clen else None,
                   "residues": [{"label_seq_id": r["polymer_seq_position"], "author": f'{r["author_residue_number"]}{r["author_insertion_code"]}', "nt": r["nucleotide_parent"], "offset": r["insertion_offset"], "has_coordinates": r["has_coordinates"], "pp": r["posterior_prob_char"]}]}
            runs.append(cur)
    # the following match after each run lies at a consensus position > anchor (anchor = preceding consensus column)
    for run in [x for x in runs if x["seq_name"] == sn]:
        nxt = next((r for r in R if r["alignment_state"] == "match" and r["polymer_seq_position"] > run["residues"][-1]["label_seq_id"]), None)
        checks["next_match_after_anchor"] += bool(nxt and nxt["rfam_consensus_pos"] > run["anchor"])
for run in runs: run["sequence"] = "".join(x["nt"] for x in run["residues"])
out = {"rfam_acc": ACC, "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"), "n_mapping_files": len(files),
       "residue_alignment_states": dict(states), "consensus_column_states": dict(columns), "flanks": dict(flanks), "checks": dict(checks),
       "n_insertion_runs": len(runs), "n_instances_with_insertions": len({r["seq_name"] for r in runs}),
       "distinct_events": dict(collections.Counter(f'{r["state"]} anchor {r["anchor"]} ({r["left"]}|{r["right"]}) {r["sequence"]}' for r in runs)),
       "terminal_runs": [r for r in runs if r["left"] is None or r["right"] is None], "runs": runs}
json.dump(out, open(Path(__file__).with_name("audit_insertions.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "runs"}, indent=1))
