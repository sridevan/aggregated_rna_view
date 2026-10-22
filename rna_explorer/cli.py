"""Command line: python -m rna_explorer build RF00167 [--out data] [--cache cache] [--cpu 2] [--max-entries N]"""
from __future__ import annotations
import argparse, json, logging, sys
from pathlib import Path
from .config import DEFAULT_CACHE_DIR, DEFAULT_OUT_DIR
from .pipeline import build_family


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="rna_explorer")
    sub = p.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="precompute mappings, validation, statistics and frontend JSON for one Rfam family")
    b.add_argument("acc")
    b.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR)
    b.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    b.add_argument("--cpu", type=int, default=2)
    b.add_argument("--max-entries", type=int, default=None)
    b.add_argument("--selection", type=Path, default=None, help="CSV of curated chains (pdb_id, label_asym_id): map only these, with per-chain ModelServer coordinates")
    b.add_argument("-v", "--verbose", action="store_true")
    o = sub.add_parser("orient", help="compute reference-orientation transforms for a built family (orientation.json)")
    o.add_argument("acc"); o.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR); o.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    n = sub.add_parser("annotate", help="build annotation files (reference elements, Motif Atlas loops, bound molecules, per-instance records, tracks) for a built family")
    n.add_argument("acc"); n.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR); n.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    n.add_argument("--raw-dir", type=Path, default=None, help="directory of previously downloaded raw responses used to seed the cache (e.g. annotation_assessment/cache)")
    n.add_argument("--offline", action="store_true", help="never contact the network; fail sources that are not cached")
    n.add_argument("-v", "--verbose", action="store_true")
    sa = sub.add_parser("sasa", help="per-nucleotide SASA of each isolated RNA chain (FreeSASA, Shrake-Rupley, probe 1.4 A) aggregated on consensus positions (sasa.json, sasa_tracks.json)")
    sa.add_argument("acc"); sa.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR); sa.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    sa.add_argument("--instances", nargs="*", default=None, help="restrict to these seq_names (pilot runs)"); sa.add_argument("--pilot-dir", type=Path, default=None, help="write outputs here instead of the family directory")
    sa.add_argument("-v", "--verbose", action="store_true")
    pk = sub.add_parser("pockets", help="fpocketR predicted pockets for every chain instance (isolated x86-64 env), mapped to consensus positions (pockets.json, pocket_tracks.json)")
    pk.add_argument("acc"); pk.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR); pk.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    pk.add_argument("--experiment-dir", type=Path, required=True, help="directory for inputs, work files and raw fpocket outputs (kept out of data/)")
    pk.add_argument("--instances", nargs="*", default=None); pk.add_argument("--no-reuse", action="store_true", help="recompute even when a run output already exists"); pk.add_argument("--no-write", action="store_true")
    pk.add_argument("-v", "--verbose", action="store_true")
    md = sub.add_parser("modifications", help="nonstandard nucleotide components of the deposited RNA polymers mapped to consensus positions (modifications.json, modification_tracks.json)")
    md.add_argument("acc"); md.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR); md.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    md.add_argument("--pdbe", choices=["cache", "off"], default="cache", help="PDBe modified_AA_or_NA cross-check (cached; 'off' skips it)")
    md.add_argument("--nakb", choices=["cache", "off"], default="cache", help="NAKB/BGSU component enrichment (cached; 'off' skips it)")
    md.add_argument("--seed-dir", type=Path, default=None, help="investigation raw directory (pdbe/batch_all_entries.json, nakb/...) used to seed the cache instead of new requests")
    md.add_argument("--offline", action="store_true", help="never contact the network; unavailable optional sources are recorded as such")
    md.add_argument("-v", "--verbose", action="store_true")
    sd = sub.add_parser("seqdiff", help="reported sequence differences of the family RNA chains (mmCIF _struct_ref_seq_dif) mapped to consensus positions (sequence_differences.json, sequence_difference_tracks.json)")
    sd.add_argument("acc"); sd.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR); sd.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    sd.add_argument("--offline", action="store_true", help="never contact the network; entries whose mmCIF is not cached are recorded as source unavailable")
    sd.add_argument("-v", "--verbose", action="store_true")
    sc = sub.add_parser("clustering", help="structural groups for the Structures tab from precomputed research matrices (structure_clustering.json); no comparison is recomputed")
    sc.add_argument("acc"); sc.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR)
    sc.add_argument("--research-root", type=Path, default=Path.home() / "Desktop/aggregated_rna_clustering", help="research worktree root (read-only)")
    sc.add_argument("-v", "--verbose", action="store_true")
    it = sub.add_parser("interactions", help="BGSU/FR3D RNA-RNA base pairs and stacking mapped to consensus positions with eligibility (interactions.json)")
    it.add_argument("acc"); it.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR); it.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    it.add_argument("--seed-dir", type=Path, default=None, help="directory with bgsu/{PDB}_{basepairs|stacking}.csv.gz and PROVENANCE.json used to seed the cache (e.g. experiments/rf00167_interactions/raw)")
    it.add_argument("--offline", action="store_true", help="never contact the network; entries whose files are not cached are recorded as annotation unavailable")
    it.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO if getattr(a, "verbose", False) else logging.WARNING, format="%(levelname)s %(message)s")
    if a.cmd == "orient":
        from .cache import Cache
        from .orientation import build_orientation
        from .pdbe import family_uses_chain_files, instance_cif
        cache = Cache(a.cache)
        fam_dir = a.out / "families" / a.acc
        chain_files = family_uses_chain_files(fam_dir)
        res = build_orientation(fam_dir, lambda inst: instance_cif(cache, inst, chain_files), fam_dir / "orientation.json")
        q = {}
        for v in res["transforms"].values(): q[v["fit_quality"]] = q.get(v["fit_quality"], 0) + 1
        print(json.dumps({"acc": a.acc, "reference": res["reference"]["seq_name"], "instances": len(res["transforms"]), "fit_quality": q}, indent=1))
        return 0
    if a.cmd == "pockets":
        from .pockets import build_pockets
        res = build_pockets(a.acc, a.cache, a.out, a.experiment_dir, a.instances, write=not a.no_write, reuse_existing_runs=not a.no_reuse)
        print(json.dumps({k: res[k] for k in ("acc", "n_instances", "status_counts", "n_pockets_total", "total_seconds")}, indent=1))
        return 0
    if a.cmd == "clustering":
        from .structure_clustering import build_structure_clustering
        res = build_structure_clustering(a.acc, a.out, a.research_root)
        print(json.dumps({"acc": a.acc, **{k: res["dataset"][k] for k in ("n_mapped_instances", "n_clustered_instances", "n_clustered_entries", "n_excluded_instances")},
                          **{m: {"k": v["clustering"]["cut"]["k"], "cut_height": v["clustering"]["cut"]["height"], "groups": [g["n_instances"] for g in v["clustering"]["groups"]],
                                 "representatives": [g["representative"] for g in v["clustering"]["groups"]]} for m, v in res["metrics"].items()}}, indent=1))
        return 0
    if a.cmd == "interactions":
        from .interactions import build_interactions
        res = build_interactions(a.acc, a.cache, a.out, a.seed_dir, a.offline)
        print(json.dumps({"acc": a.acc, **res["summary"], "default_numbering_reference": res["default_numbering_reference"]}, indent=1))
        return 0
    if a.cmd == "seqdiff":
        from .sequence_differences import build_sequence_differences
        res = build_sequence_differences(a.acc, a.cache, a.out, a.offline)
        print(json.dumps({"acc": res["acc"], "counts": res["counts"]}, indent=1))
        return 0
    if a.cmd == "modifications":
        from .modifications import build_modifications
        res = build_modifications(a.acc, a.cache, a.out, a.pdbe, a.nakb, a.seed_dir, a.offline)
        print(json.dumps({"acc": res["acc"], "pdbe": res["pdbe"], "nakb": res["nakb"]} | res["counts"], indent=1))
        return 0
    if a.cmd == "sasa":
        from .sasa import build_sasa
        res = build_sasa(a.acc, a.cache, a.out, a.instances, write=a.pilot_dir is None)
        if a.pilot_dir is not None:
            a.pilot_dir.mkdir(parents=True, exist_ok=True)
            json.dump(res["sasa"], open(a.pilot_dir / "sasa.json", "w"), indent=1); json.dump(res["tracks"], open(a.pilot_dir / "sasa_tracks.json", "w"), indent=1)
        print(json.dumps({k: res[k] for k in ("acc", "n_measured", "n_failed", "status_counts", "n_positions_with_value", "max_median_A2")}, indent=1))
        return 0
    if a.cmd == "annotate":
        from .annotations import build_annotations
        summary = build_annotations(a.acc, a.cache, a.out, a.raw_dir, a.offline)
        print(json.dumps({k: summary[k] for k in ("rfam_acc", "n_instances", "n_entries", "n_unresolved_identifiers", "file_sizes_bytes")} | {"loops": summary["loops"]["n_loops"], "ligands": summary["ligands"]["n_bound_molecule_records"]}, indent=1))
        return 0
    if a.cmd == "build":
        print(json.dumps(build_family(a.acc, a.cache, a.out, a.cpu, a.max_entries, selection=a.selection), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
