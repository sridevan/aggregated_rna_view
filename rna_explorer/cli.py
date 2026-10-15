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
    b.add_argument("-v", "--verbose", action="store_true")
    o = sub.add_parser("orient", help="compute reference-orientation transforms for a built family (orientation.json)")
    o.add_argument("acc"); o.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR); o.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    n = sub.add_parser("annotate", help="build annotation files (reference elements, Motif Atlas loops, bound molecules, per-instance records, tracks) for a built family")
    n.add_argument("acc"); n.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR); n.add_argument("--cache", type=Path, default=DEFAULT_CACHE_DIR)
    n.add_argument("--raw-dir", type=Path, default=None, help="directory of previously downloaded raw responses used to seed the cache (e.g. annotation_assessment/cache)")
    n.add_argument("--offline", action="store_true", help="never contact the network; fail sources that are not cached")
    n.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO if getattr(a, "verbose", False) else logging.WARNING, format="%(levelname)s %(message)s")
    if a.cmd == "orient":
        from .cache import Cache
        from .orientation import build_orientation
        from .pdbe import get_updated_cif
        cache = Cache(a.cache)
        fam_dir = a.out / "families" / a.acc
        res = build_orientation(fam_dir, lambda pdb: get_updated_cif(cache, pdb), fam_dir / "orientation.json")
        q = {}
        for v in res["transforms"].values(): q[v["fit_quality"]] = q.get(v["fit_quality"], 0) + 1
        print(json.dumps({"acc": a.acc, "reference": res["reference"]["seq_name"], "instances": len(res["transforms"]), "fit_quality": q}, indent=1))
        return 0
    if a.cmd == "annotate":
        from .annotations import build_annotations
        summary = build_annotations(a.acc, a.cache, a.out, a.raw_dir, a.offline)
        print(json.dumps({k: summary[k] for k in ("rfam_acc", "n_instances", "n_entries", "n_unresolved_identifiers", "file_sizes_bytes")} | {"loops": summary["loops"]["n_loops"], "ligands": summary["ligands"]["n_bound_molecule_records"]}, indent=1))
        return 0
    if a.cmd == "build":
        print(json.dumps(build_family(a.acc, a.cache, a.out, a.cpu, a.max_entries), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
