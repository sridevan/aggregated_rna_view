"""Is the PDBe 2D topology already the R2DT template layout? For every chain with both an accepted R2DT layout and a PDBe topology file,
compare the drawn positions of the same label_seq_id after removing the best translation and uniform scale (the two canvases differ).

    python3 experiments/r2dt_viewer/compare_pdbe.py RF00167 ... > experiments/r2dt_viewer/pdbe_vs_r2dt.json

PDBe files (https://www.ebi.ac.uk/pdbe/static/entry/<pdb>_<entity>_<AUTH>.json) are cached in work/pdbe/.
"""
from __future__ import annotations

import json
import math
import statistics
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = Path(__file__).resolve().parent / "work" / "pdbe"


def pdbe_topology(pdb: str, entity: int, auth: str) -> dict | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{pdb}_{entity}_{auth.upper()}.json"
    if not f.exists():
        try:
            with urllib.request.urlopen(f"https://www.ebi.ac.uk/pdbe/static/entry/{f.name}", timeout=60) as r:
                f.write_bytes(r.read())
        except urllib.error.HTTPError:
            f.write_text("null")
    return json.loads(f.read_text())


def points(api: dict) -> dict[int, tuple[float, float]]:
    """label_seq_id -> drawn point (end point of each nucleotide's path, as the viewer draws it)."""
    out = {}
    for k, lab in enumerate(api["label_seq_ids"][1:-1]):
        p = api["svg_paths"][k + 2].replace("M", " ").replace("L", " ").replace(",", " ").split()
        if lab is not None and len(p) >= 4:
            out[int(lab)] = (float(p[2]), float(p[3]))
    return out


def similarity_residual(a: dict, b: dict) -> tuple[int, float] | None:
    """median distance after the best translation + uniform scale of b onto a, in units of a's median nucleotide spacing"""
    common = sorted(set(a) & set(b))
    if len(common) < 10:
        return None
    ax, ay = statistics.fmean(a[c][0] for c in common), statistics.fmean(a[c][1] for c in common)
    bx, by = statistics.fmean(b[c][0] for c in common), statistics.fmean(b[c][1] for c in common)
    num = sum((a[c][0] - ax) * (b[c][0] - bx) + (a[c][1] - ay) * (b[c][1] - by) for c in common)
    den = sum((b[c][0] - bx) ** 2 + (b[c][1] - by) ** 2 for c in common) or 1
    s = num / den
    d = [math.dist((a[c][0] - ax, a[c][1] - ay), (s * (b[c][0] - bx), s * (b[c][1] - by))) for c in common]
    sp = statistics.median(math.dist(a[c], a[c + 1]) for c in common if c + 1 in a) or 1
    return len(common), statistics.median(d) / sp


def family(acc: str) -> dict:
    fam = ROOT / "data" / "families" / acc
    idx = json.loads((fam / "r2dt" / "index.json").read_text())
    inst = {i["seq_name"]: i for i in json.loads((fam / "instances.json").read_text())}
    rows = []
    for r in idx["chains"]:
        if r["status"] != "accepted":
            continue
        i = inst[r["seq_name"]]
        p = pdbe_topology(i["pdb_id"], i["entity_id"], i["auth_asym_id"])
        if not p:
            rows.append({"seq_name": r["seq_name"], "pdbe": False, "r2dt_layout": r["layout"], "r2dt_template": r["template"]}); continue
        res = similarity_residual(points(json.loads((fam / "r2dt" / r["seq_name"] / "api.json").read_text())), points(p))
        rows.append({"seq_name": r["seq_name"], "pdbe": True, "r2dt_layout": r["layout"], "r2dt_template": r["template"],
                     "shared_positions": res[0] if res else 0, "median_residual_in_nt_spacings": round(res[1], 3) if res else None})
    both = [x for x in rows if x["pdbe"] and x["median_residual_in_nt_spacings"] is not None]
    return {"family": acc, "accepted_r2dt": len(rows), "with_pdbe_topology": len(both),
            "same_layout_as_pdbe (<0.1 spacing)": sum(x["median_residual_in_nt_spacings"] < 0.1 for x in both),
            "different_layout (>0.5 spacing)": sum(x["median_residual_in_nt_spacings"] > 0.5 for x in both),
            "r2dt_only (no PDBe topology)": [x["seq_name"] for x in rows if not x["pdbe"]], "chains": rows}


if __name__ == "__main__":
    print(json.dumps([family(a) for a in sys.argv[1:]], indent=1))
