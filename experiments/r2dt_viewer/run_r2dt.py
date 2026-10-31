"""Exploration: R2DT (2d-3d-viewer branch) 2D layouts for every chain of a family, validated against our mappings.

Runs ``r2dt.py pdb_2d_3d <PDB> --chain <auth> --format cif`` in the ``rnacentral/r2dt:develop`` image with the
``2d-3d-viewer`` branch code mounted over it (the image predates ``pdb_2d_3d``), keeps only the viewer JSON
(``api.json``, ``fr3d.json``, ``lbn.json``; structure files stay in ``work/``) and accepts a layout only when every
drawn position equals our mapping residue at the same ``label_seq_id`` (base, author number). Modified nucleotides,
which R2DT draws as ``N``, are allowed; dropped positions are recorded.

    python3 experiments/r2dt_viewer/run_r2dt.py RF00167 --r2dt-src <R2DT checkout on 2d-3d-viewer>

Writes ``data/families/<ACC>/r2dt/<seq_name>/{api,fr3d,lbn}.json`` for accepted chains and
``data/families/<ACC>/r2dt/index.json`` (all chains, accepted or not, with reasons and provenance).
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
IMAGE = "rnacentral/r2dt:develop"


def run_chain(pdb: str, chain: str, src: Path, work: Path, mode: str) -> tuple[Path, str]:
    name = f"{pdb.upper()}_{chain}"
    out = work / name
    log = work / f"{name}.log"
    if not (out / "viewer" / "api.json").exists():
        shutil.rmtree(out, ignore_errors=True)
        cmd = ["docker", "run", "--rm",
               "-v", f"{src}/r2dt.py:/rna/r2dt/r2dt.py", "-v", f"{src}/utils:/rna/r2dt/utils",
               "-v", f"{src}/data/viewer:/rna/r2dt/data/viewer", "-v", f"{work}:/rna/r2dt/output",
               "-w", "/rna/r2dt", IMAGE, "./r2dt.py", "pdb_2d_3d", pdb.upper(), f"output/{name}",
               "--chain", chain, "--format", "cif", "--mode", mode]
        with log.open("w") as fh:
            subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, check=False)
    return out / "viewer", log.read_text(errors="ignore") if log.exists() else ""


def validate(api: dict, residues: list[dict]) -> dict:
    ours = {r["polymer_seq_position"]: r for r in residues}
    labels = api["label_seq_ids"][1:-1]
    auth = api["auth_seq_ids"][1:-1]
    seq = api["sequence"]
    unobserved = set(api.get("unobserved_label_seq_ids", []))
    problems, modified_as_n = [], []
    for k, lab in enumerate(labels):
        r = ours.get(lab)
        if r is None:
            problems.append(f"{lab}: not in mapping"); continue
        if seq[k] != r["nucleotide_parent"]:
            if seq[k] in "NX" and r["is_modified"]:   # R2DT draws modified residues as N (X when unresolved)
                modified_as_n.append(lab)
            else:
                problems.append(f"{lab}: base {seq[k]} != {r['nucleotide_parent']}")
        # R2DT falls back to label numbering for the author number of unresolved residues (display only)
        if lab not in unobserved and auth[k] != r["author_residue_number"]:
            problems.append(f"{lab}: auth {auth[k]} != {r['author_residue_number']}")
        if (lab in unobserved) == r["has_coordinates"]:
            problems.append(f"{lab}: observed flag differs")
    drawn = set(labels)
    not_drawn = sorted(set(ours) - drawn)
    not_drawn_matched = [p for p in not_drawn if ours[p]["alignment_state"] == "match"]
    return {"ok": not problems, "n_drawn": len(labels), "n_chain": len(ours), "problems": problems[:20],
            "n_problems": len(problems), "modified_drawn_as_N": modified_as_n, "not_drawn": not_drawn,
            "not_drawn_consensus_matched": not_drawn_matched}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("family")
    ap.add_argument("--r2dt-src", type=Path, required=True)
    ap.add_argument("--mode", default="auto", choices=["auto", "templated", "templatefree"])
    ap.add_argument("--limit", type=int)
    a = ap.parse_args()
    fam = ROOT / "data" / "families" / a.family
    work = HERE / "work" / a.family
    work.mkdir(parents=True, exist_ok=True)
    dest = fam / "r2dt"
    dest.mkdir(exist_ok=True)
    instances = json.loads((fam / "instances.json").read_text())[: a.limit]
    commit = subprocess.run(["git", "-C", str(a.r2dt_src), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    digest = subprocess.run(["docker", "image", "inspect", IMAGE, "--format", "{{index .RepoDigests 0}}"],
                            capture_output=True, text=True).stdout.strip()
    rows = []
    for i, inst in enumerate(instances, 1):
        seq_name = inst["seq_name"]
        viewer, log = run_chain(inst["pdb_id"], inst["auth_asym_id"], a.r2dt_src.resolve(), work, a.mode)
        m = re.search(r"Matched template: (\S+)", log)
        row = {"seq_name": seq_name, "pdb_id": inst["pdb_id"], "auth_asym_id": inst["auth_asym_id"],
               "template": m.group(1) if m else None,
               "layout": "templated" if m else ("templatefree" if "templatefree" in log else None)}
        if not (viewer / "api.json").exists():
            row.update(status="failed", reason=(log.strip().splitlines() or ["no log"])[-1][:300])
        else:
            api = json.loads((viewer / "api.json").read_text())
            residues = json.loads((fam / inst["mapping_file"]).read_text())["residues"]
            v = validate(api, residues)
            row.update(status="accepted" if v["ok"] else "rejected", validation=v)
            if "Could not map to full sequence" in log:
                row["r2dt_warning"] = "Could not map to full sequence"
            out = dest / seq_name
            if v["ok"]:
                out.mkdir(exist_ok=True)
                for f in ("api.json", "fr3d.json", "lbn.json"):
                    if (viewer / f).exists():
                        shutil.copy(viewer / f, out / f)
            else:
                shutil.rmtree(out, ignore_errors=True)
        rows.append(row)
        print(f"[{i}/{len(instances)}] {seq_name} {row['status']} {row.get('template')}", flush=True)
    index = {
        "family": a.family,
        "provenance": {"tool": "R2DT pdb_2d_3d", "r2dt_repo": "https://github.com/r2dt-bio/R2DT",
                       "r2dt_branch": "2d-3d-viewer", "r2dt_commit": commit, "docker_image": digest,
                       "mode": a.mode, "structure_format": "cif (downloaded by R2DT)",
                       "generated": dt.date.today().isoformat(),
                       "acceptance": "every drawn position equals the mapping residue at the same label_seq_id "
                                     "(base, author number when observed, observed flag); modified nucleotides drawn as N or X allowed"},
        "counts": {s: sum(r["status"] == s for r in rows) for s in ("accepted", "rejected", "failed")},
        "chains": rows,
    }
    (dest / "index.json").write_text(json.dumps(index, indent=1) + "\n")
    print(index["counts"])


if __name__ == "__main__":
    main()
