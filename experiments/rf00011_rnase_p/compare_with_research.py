"""Phase 1 check: our selection build of RF00005 (data/families/RF00011) against the research build of the same 50 chains
(~/Desktop/aggregated_rna_clustering/.../data/rf00011_50, read-only). Compares every residue (state, consensus position, insertion
anchor/offset, coordinates) and every consensus column per chain. Writes compare_with_research.json next to this script."""
import json
from pathlib import Path

ours = Path(__file__).resolve().parents[2] / "data/families/RF00011"
theirs = Path.home() / "Desktop/aggregated_rna_clustering/experiments/rna_structural_clustering/data/rf00011_50/families/RF00011"
K = ("polymer_seq_position", "alignment_state", "rfam_consensus_pos", "insertion_anchor", "insertion_offset", "has_coordinates", "author_residue_number", "residue_name")
oi, ti = {i["seq_name"]: i for i in json.load(open(ours / "instances.json"))}, {i["seq_name"]: i for i in json.load(open(theirs / "instances.json"))}
out = {"n_ours": len(oi), "n_research": len(ti), "only_ours": sorted(set(oi) - set(ti)), "only_research": sorted(set(ti) - set(oi)), "chains": {}}
for sn in sorted(set(oi) & set(ti)):
    a, b = json.load(open(ours / "mappings" / f"{sn}.json")), json.load(open(theirs / "mappings" / f"{sn}.json"))
    ra, rb = [tuple(r.get(k) for k in K) for r in a["residues"]], [tuple(r.get(k) for k in K) for r in b["residues"]]
    ca, cb = [(c["rfam_consensus_pos"], c["state"], c["polymer_seq_position"]) for c in a["consensus_columns"]], [(c["rfam_consensus_pos"], c["state"], c["polymer_seq_position"]) for c in b["consensus_columns"]]
    out["chains"][sn] = {"residues_identical": ra == rb, "n_residue_differences": sum(x != y for x, y in zip(ra, rb)) + abs(len(ra) - len(rb)),
                         "columns_identical": ca == cb, "validation_ok": oi[sn]["validation_ok"], "n_match": oi[sn]["n_match"], "n_insert": oi[sn]["n_insert"],
                         "n_match_with_coordinates": oi[sn]["n_match_with_coordinates"], "organism": oi[sn].get("organism"), "pdbe_2d_available": oi[sn].get("pdbe_2d_available")}
c = out["chains"].values()
out["summary"] = {"chains_compared": len(out["chains"]), "residues_identical": sum(x["residues_identical"] for x in c), "columns_identical": sum(x["columns_identical"] for x in c),
                  "validation_ok": sum(bool(x["validation_ok"]) for x in c), "pdbe_2d_available": sum(x["pdbe_2d_available"] is True for x in c)}
json.dump(out, open(Path(__file__).with_name("compare_with_research.json"), "w"), indent=1)
print(json.dumps({k: out[k] for k in ("n_ours", "n_research", "only_ours", "only_research", "summary")}, indent=1))
for sn, x in out["chains"].items():
    if not (x["residues_identical"] and x["columns_identical"]): print("DIFF", sn, x)
