"""Stockholm parsing and RF-column typing (validated semantics: letter = match column, '.' = insert column, '~' = local-end column)."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

GAPS = frozenset("-._~")


@dataclass
class Alignment:
    order: list[str] = field(default_factory=list)
    seqs: dict[str, str] = field(default_factory=dict)
    gr: dict[str, dict[str, str]] = field(default_factory=lambda: defaultdict(dict))
    gc: dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.gc.get("RF", "")) if self.gc.get("RF") else (len(next(iter(self.seqs.values()))) if self.seqs else 0)


def parse_stockholm(path: Path | str) -> Alignment:
    """Parse a (possibly multi-block) Stockholm file, concatenating sequence, #=GR and #=GC lines."""
    aln = Alignment()
    seqs: dict[str, str] = defaultdict(str)
    gr: dict[str, dict[str, str]] = defaultdict(lambda: defaultdict(str))
    gc: dict[str, str] = defaultdict(str)
    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line or line.startswith("# STOCKHOLM") or line == "//":
                continue
            if line.startswith("#=GC"):
                _, tag, txt = line.split(None, 2)
                gc[tag] += txt.strip()
            elif line.startswith("#=GR"):
                _, name, tag, txt = line.split(None, 3)
                gr[name][tag] += txt.strip()
            elif line.startswith("#"):
                continue
            else:
                name, txt = line.split(None, 1)
                if name not in seqs:
                    aln.order.append(name)
                seqs[name] += txt.strip()
    aln.seqs = dict(seqs)
    aln.gr = {k: dict(v) for k, v in gr.items()}
    aln.gc = dict(gc)
    L = len(aln.gc.get("RF", ""))
    if L and any(len(s) != L for s in aln.seqs.values()):
        raise ValueError("ragged Stockholm alignment")
    return aln


@dataclass
class ColumnMap:
    rf: str
    ss_cons: str
    col_type: list[str]            # 'match' | 'insert' | 'local_end' per column
    col_cons: list[Optional[int]]  # consensus position (1..CLEN) for match columns
    prev_cons: list[int]           # consensus position preceding (or at) the column; 0 before the first match column
    clen: int


def type_columns(rf: str, ss_cons: str = "") -> ColumnMap:
    col_type: list[str] = []
    col_cons: list[Optional[int]] = []
    prev: list[int] = []
    n = 0
    for ch in rf:
        if ch == ".":
            col_type.append("insert"); col_cons.append(None)
        elif ch == "~":
            col_type.append("local_end"); col_cons.append(None)
        else:
            n += 1
            col_type.append("match"); col_cons.append(n)
        prev.append(n)
    if ss_cons and len(ss_cons) != len(rf):
        raise ValueError("SS_cons length differs from RF")
    return ColumnMap(rf, ss_cons or "." * len(rf), col_type, col_cons, prev, n)
