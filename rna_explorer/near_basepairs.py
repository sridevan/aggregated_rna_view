"""Near base pairs of the family chains from a pinned local fr3d-python run (near_basepairs.json), read by `interactions`.

Why: BGSU's per-entry base-pair CSVs (the production source of true base pairs) contain no near ('n'-prefixed) classes, because fr3d-python
writes near base pairs only when the category 'near' is requested (experiments/sscons_vs_3d/NEAR_CHECK.md). They are therefore computed
locally with BGSU's own annotation code:

    NA_pairwise_interactions.py -c basepair,near --chain <family author chains> <PDB>.cif

on the PDBe updated mmCIF of each entry. The mmCIF is streamed into a temporary directory, annotated and deleted (never stored). Only near
base pairs whose two residues both belong to a family chain of the entry are kept, as BGSU-style rows (unit id, class, unit id); true base
pairs still come from BGSU (only their count is recorded, for comparison). The fr3d-python checkout must be at a fixed commit, recorded with
the run date, interpreter and package versions in the provenance. Entries already annotated with the same commit and chain set are reused
unless `force`.
"""
from __future__ import annotations
import concurrent.futures as cf, datetime as _dt, hashlib, json, logging, shutil, subprocess, tempfile
from pathlib import Path
from typing import Any, Optional

import requests

from . import __version__
from .export import validate_json

log = logging.getLogger(__name__)
CIF_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}_updated.cif"
REPOSITORY = "https://github.com/BGSU-RNA/fr3d-python"
SCRIPT = "fr3d/classifiers/NA_pairwise_interactions.py"
LICENCE = "fr3d-python: license.txt is MIT (Copyright 2026 BGSU RNA Group); pyproject.toml declares Apache-2.0 (both permissive; discrepancy upstream)"
PACKAGES = ("numpy", "scipy", "biopython", "mmcif-pdbx")


def tool_info(fr3d_dir: Path, python: str) -> dict[str, Any]:
    """commit (refusing a dirty checkout), interpreter and package versions of the fr3d-python installation"""
    git = lambda *a: subprocess.run(["git", "-C", str(fr3d_dir), *a], capture_output=True, text=True, check=True).stdout.strip()
    commit = git("rev-parse", "HEAD")
    if git("status", "--porcelain", "--untracked-files=no"):
        raise RuntimeError(f"fr3d-python checkout {fr3d_dir} has local modifications; use a clean checkout of a fixed commit")
    code = "import importlib.metadata as m, sys, json; print(json.dumps({'python': sys.version.split()[0], **{p: m.version(p) for p in %r}}))" % (PACKAGES,)
    versions = json.loads(subprocess.run([python, "-c", code], capture_output=True, text=True, check=True).stdout)
    return {"tool": f"fr3d-python {SCRIPT}", "repository": REPOSITORY, "commit": commit, "commit_date": git("log", "-1", "--format=%cs"),
            "categories": "basepair,near", "environment": versions, "licence": LICENCE}


def is_near(cls: str) -> bool:
    return cls.startswith("n")


def run_entry(pdb: str, chains: list[str], fr3d_dir: Path, python: str, tmp_root: Optional[Path]) -> dict[str, Any]:
    """annotate one entry (family chains only); the mmCIF is deleted with its temporary directory"""
    d = Path(tempfile.mkdtemp(dir=tmp_root))
    try:
        cif, h = d / f"{pdb.upper()}.cif", hashlib.sha256()
        url = CIF_URL.format(pdb=pdb.lower())
        with requests.get(url, stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(cif, "wb") as f:
                for block in r.iter_content(1 << 20):
                    f.write(block); h.update(block)
        cmd = [python, str(Path(fr3d_dir) / SCRIPT), "-i", str(d), "-o", str(d), "-c", "basepair,near", "--chain", ",".join(chains), cif.name]
        p = subprocess.run(cmd, cwd=fr3d_dir, env={"PYTHONPATH": str(fr3d_dir), "PATH": "/usr/bin:/bin"}, capture_output=True, text=True, timeout=7200)
        outs = sorted(d.glob("*basepair*.txt"))
        if not outs:
            raise RuntimeError(f"fr3d-python wrote no base-pair file (exit {p.returncode}): {p.stderr[-300:]}")
        want = set(chains)
        near, n_true = set(), 0
        for f in outs:
            for line in open(f):
                x = line.rstrip("\n").split("\t")
                if len(x) < 3 or not all(u.split("|")[2:3] and u.split("|")[2] in want for u in (x[0], x[2])):
                    continue
                if is_near(x[1]):
                    near.add((x[0], x[1], x[2]))
                else:
                    n_true += 1
        return {"status": "ok", "url": url, "retrieved": _dt.date.today().isoformat(), "cif_sha256": h.hexdigest(), "n_near_rows": len(near),
                "n_true_rows_local": n_true, "rows": [list(r) for r in sorted(near)]}
    finally:
        shutil.rmtree(d, ignore_errors=True)


def build_near_basepairs(acc: str, out_dir: Path, fr3d_dir: Path, python: str, workers: int = 4, tmp_dir: Optional[Path] = None,
                         force: bool = False, write: bool = True) -> dict[str, Any]:
    fam_dir = Path(out_dir) / "families" / acc
    instances = json.load(open(fam_dir / "instances.json"))
    chains: dict[str, set[str]] = {}
    for i in instances:
        chains.setdefault(i["pdb_id"].upper(), set()).add(i["auth_asym_id"])
    info = tool_info(Path(fr3d_dir), python)
    path = fam_dir / "near_basepairs.json"
    old = json.load(open(path)) if path.exists() and not force else {}
    reuse = old.get("provenance", {}).get("commit") == info["commit"]
    entries: dict[str, dict[str, Any]] = {}
    todo = []
    for pdb, ch in sorted(chains.items()):
        prev = old.get("entries", {}).get(pdb) if reuse else None
        if prev and prev["status"] == "ok" and prev["chains"] == sorted(ch):
            entries[pdb] = prev
        else:
            todo.append(pdb)
    log.info("%s: %d entries, %d to annotate (%d reused)", acc, len(chains), len(todo), len(entries))

    def job(pdb: str) -> tuple[str, dict[str, Any]]:
        try:
            return pdb, run_entry(pdb, sorted(chains[pdb]), Path(fr3d_dir), python, tmp_dir)
        except Exception as e:   # noqa: BLE001 - recorded as unavailable; downstream these chains have no near base pairs
            return pdb, {"status": "failed", "note": str(e)[:300], "rows": []}

    with cf.ThreadPoolExecutor(max(1, workers)) as ex:
        for pdb, res in ex.map(job, todo):
            entries[pdb] = {"chains": sorted(chains[pdb]), **res}
            log.info("%s %s %s", acc, pdb, res["status"] if res["status"] != "ok" else f"{res['n_near_rows']} near")
    out = {"schema_version": "1.0", "rfam_acc": acc,
           "definitions": {"rows": "[unit_id_1, fr3d_class, unit_id_2] as written by fr3d-python (BGSU unit ids); near (n-prefixed) base-pair classes "
                                   "only, both residues in family chains of the entry; all models and symmetry operators as annotated (filtered "
                                   "downstream like the BGSU files)",
                           "n_true_rows_local": "true (non-near) base-pair rows of the same run between family chains; not used downstream (BGSU is "
                                                "the source of true base pairs), kept for comparison"},
           "provenance": {**info, "coordinates": "PDBe updated mmCIF (streamed to a temporary directory, deleted after annotation)",
                          "url_template": CIF_URL, "pipeline_version": __version__,
                          "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                          "method": "experiments/sscons_vs_3d/NEAR_CHECK.md (validated: local true pairs identical to BGSU in 298 of 321 chains)"},
           "summary": {"n_entries": len(entries), "n_ok": sum(e["status"] == "ok" for e in entries.values()),
                       "n_failed": sum(e["status"] != "ok" for e in entries.values()), "n_near_rows": sum(len(e["rows"]) for e in entries.values())},
           "entries": dict(sorted(entries.items()))}
    validate_json(out, "near_basepairs")
    if write:
        with open(path, "w") as f:
            json.dump(out, f, separators=(",", ":"))
    return out


def load_near_basepairs(fam_dir: Path) -> Optional[dict[str, Any]]:
    p = Path(fam_dir) / "near_basepairs.json"
    return json.load(open(p)) if p.exists() else None
