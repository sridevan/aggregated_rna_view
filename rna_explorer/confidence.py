"""Alignment confidence statistics per consensus position.
Sequence level: over deduplicated unique mapped sequences. Structure level: over chain instances, with coordinate availability.
Raw counts and sums are stored; no composite score."""
from __future__ import annotations
from typing import Any
from .mapping import AlignedSequence, ChainMapping, PP_BINS, PP_THRESHOLDS, pp_category, pp_value
from .stockholm import ColumnMap


def _empty(clen: int) -> dict[str, list[Any]]:
    return {"n_match": [0] * clen, "n_deletion": [0] * clen, "n_uncovered": [0] * clen, "pp_n": [0] * clen, "pp_sum": [0.0] * clen,
            "pp_min": [None] * clen, "pp_high": [0] * clen, "pp_medium": [0] * clen, "pp_low": [0] * clen}


def _add_pp(acc: dict[str, list[Any]], i: int, ch: str) -> None:
    v = pp_value(ch)
    if v is None:
        return
    acc["pp_n"][i] += 1
    acc["pp_sum"][i] += v
    acc["pp_min"][i] = v if acc["pp_min"][i] is None else min(acc["pp_min"][i], v)
    acc[f"pp_{pp_category(ch)}"][i] += 1


def _finish(acc: dict[str, list[Any]], n: int) -> dict[str, Any]:
    clen = len(acc["n_match"])
    acc["coverage"] = [round(acc["n_match"][i] / n, 4) if n else None for i in range(clen)]
    acc["pp_mean"] = [round(acc["pp_sum"][i] / acc["pp_n"][i], 4) if acc["pp_n"][i] else None for i in range(clen)]
    return acc


def consensus_statistics(cm: ColumnMap, uniques: dict[str, AlignedSequence], chains: list[ChainMapping], provenance: dict[str, Any]) -> dict[str, Any]:
    clen = cm.clen
    seq = _empty(clen)
    for a in uniques.values():
        for c in a.consensus:
            i = c.consensus_pos - 1
            seq[f"n_{c.state}"][i] += 1
            if c.state == "match":
                _add_pp(seq, i, c.pp_char)
    st = _empty(clen)
    st["n_match_with_coordinates"] = [0] * clen
    st["n_instances_with_region"] = [0] * clen
    for m in chains:
        for c in m.cons_rows:
            i = c.rfam_consensus_pos - 1
            st[f"n_{c.state}"][i] += 1
            if c.state == "match":
                _add_pp(st, i, c.posterior_prob)
                if c.has_coordinates:
                    st["n_match_with_coordinates"][i] += 1
    seq_out = _finish(seq, len(uniques))
    st_out = _finish(st, len(chains))
    st_out["structural_coverage"] = [round(st["n_match_with_coordinates"][i] / len(chains), 4) if chains else None for i in range(clen)]
    seq_out["n_unique_sequences"] = len(uniques)
    st_out["n_instances"] = len(chains)
    return {"rfam_acc": provenance.get("rfam_acc"), "clen": clen, "rf": list(c for c, t in zip(cm.rf, cm.col_type) if t == "match"),
            "ss_cons": list(c for c, t in zip(cm.ss_cons, cm.col_type) if t == "match"),
            "sequence_level": seq_out, "structure_level": st_out,
            "pp_bins": PP_BINS, "pp_display_thresholds": PP_THRESHOLDS,
            "notes": ["sequence_level counts each unique mapped sequence once; structure_level counts every chain instance",
                      "pp values are Infernal posterior-probability bin midpoints; categories are provisional display bands (high >= 0.95, medium >= 0.80)",
                      "structural_coverage = fraction of chain instances with a match residue that has coordinates at the position"],
            "provenance": provenance}
