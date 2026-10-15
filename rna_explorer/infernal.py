"""cmalign runner and score-file parser (local mode by default; global mode is not used, see poc_rf00177 report)."""
from __future__ import annotations
import subprocess
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Optional, Sequence


@dataclass
class ScoreRow:
    name: str
    length: int
    cm_from: int
    cm_to: int
    trunc: str
    bit_score: float
    avg_pp: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def cmalign_version(cmalign: Path | str) -> str:
    out = subprocess.run([str(cmalign), "-h"], capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        if "INFERNAL" in line:
            return line.strip("# ").strip()
    return "unknown"


def parse_scores(path: Path | str) -> dict[str, ScoreRow]:
    rows: dict[str, ScoreRow] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            rows[t[1]] = ScoreRow(t[1], int(t[2]), int(t[3]), int(t[4]), t[5], float(t[6]), float(t[7]))
    return rows


def run_cmalign(cmalign: Path | str, cm_path: Path | str, fasta: Path | str, out_sto: Path | str, scores_path: Path | str,
                cpu: int = 2, extra: Sequence[str] = ()) -> dict[str, ScoreRow]:
    cmd = [str(cmalign), "--cpu", str(cpu), "-o", str(out_sto), "--sfile", str(scores_path), *extra, str(cm_path), str(fasta)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"cmalign failed ({r.returncode}): {r.stderr[-2000:]}")
    return parse_scores(scores_path)


def cmalign_options(cpu: int, extra: Sequence[str] = ()) -> str:
    return " ".join(["--cpu", str(cpu), *extra]) or "(defaults)"
