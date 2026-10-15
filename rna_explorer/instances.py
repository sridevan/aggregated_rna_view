"""Instance table rows (one per mapped chain) and representative-structure selection."""
from __future__ import annotations
from typing import Any, Optional
from .infernal import ScoreRow
from .mapping import ChainMapping


def instance_row(m: ChainMapping, search_doc: Optional[dict[str, Any]], rfam_hit: Optional[dict[str, Any]], clan_hits: list[tuple[str, float]],
                 score: Optional[ScoreRow], acc: str, validation_ok: Optional[bool]) -> dict[str, Any]:
    info, s = m.info, m.summary
    best = max(clan_hits, key=lambda x: x[1]) if clan_hits else None
    d = search_doc or {}
    first = lambda k: (d.get(k) or [None])[0] if isinstance(d.get(k), list) else d.get(k)
    return {"pdb_id": info.pdb_id, "entity_id": info.entity_id, "label_asym_id": info.label_asym_id, "auth_asym_id": info.auth_asym_id, "seq_name": info.seq_name,
            "sequence_key": info.sequence_key, "mapping_file": f"mappings/{info.seq_name}.json",
            "organism": first("organism_scientific_name") or info.cif_organism or info.api_organism, "tax_id": first("tax_id") or info.cif_tax_id or info.api_tax_id,
            "organism_source": "pdbe_search" if first("organism_scientific_name") else ("mmcif" if info.cif_organism else ("pdbe_molecules" if info.api_organism else None)),
            "molecule_name": info.molecule_name, "experimental_method": first("experimental_method"), "resolution": first("resolution"), "release_year": first("release_year"), "title": d.get("title"),
            "chain_length": len(info.sequence), "map_start": info.map_start, "map_end": info.map_end, "map_auth_start": info.map_auth_start, "map_auth_end": info.map_auth_end,
            "mapped_length": info.map_end - info.map_start + 1, "n_modified_in_mapped": s["n_modified_in_mapped"],
            "rfam_bit_score": rfam_hit.get("bit_score") if rfam_hit else None, "rfam_evalue": rfam_hit.get("evalue_score") if rfam_hit else None,
            "rfam_cm_start": rfam_hit.get("cm_start") if rfam_hit else None, "rfam_cm_end": rfam_hit.get("cm_end") if rfam_hit else None,
            "clan_families_hit": [{"acc": f, "bits": b} for f, b in sorted(clan_hits, key=lambda x: -x[1])], "clan_best_family": best[0] if best else None,
            "is_clan_best": (best[0] == acc) if best else None,
            "cmalign": score.to_dict() if score else None,
            "n_match": s["n_match"], "n_insert": s["n_insert"], "n_local_end": s["n_local_end"], "n_deleted_in_span": s["n_deleted_in_span"], "n_uncovered": s["n_uncovered_outside_span"],
            "consensus_span": s["cons_span"], "fraction_clen_matched": s["fraction_CLEN_matched"], "n_mapped_with_coordinates": s["n_mapped_with_coordinates"],
            "n_match_with_coordinates": s["n_match_with_coordinates"], "pp_star_fraction": s["pp_star_fraction"], "pp_category_counts": s["pp_category_counts"],
            "validation_ok": validation_ok}


def choose_representative(rows: list[dict[str, Any]]) -> Optional[dict[str, Any]]:
    """Documented rule: clan-best (or unknown), most consensus-matched residues with coordinates, then best (lowest) resolution, then earliest pdb_id."""
    cands = [r for r in rows if r.get("is_clan_best") in (True, None)]
    if not cands:
        return None
    return min(cands, key=lambda r: (-r["n_match_with_coordinates"], r["resolution"] if isinstance(r.get("resolution"), (int, float)) else 99.0, r["pdb_id"]))


def choose_display_representative(rows: list[dict[str, Any]], min_coord_fraction: float = 0.95, min_clen_fraction: Optional[float] = None) -> Optional[dict[str, Any]]:
    """Frontend demonstration structure. Documented rule, evaluated in this order:
    eligible = clan-best (or unknown) AND X-ray diffraction AND PDBe 2D topology available AND an identifiable natural species
    (species-level name resolved from NCBI taxonomy, not 'synthetic construct') AND >= 95 % of consensus-matched residues have coordinates
    AND model coverage >= the family median of fraction_clen_matched (so a short fragment never wins on resolution alone);
    ranking = best (lowest) resolution, then most matched residues with coordinates, then pdb_id, then auth chain.
    Falls back to None when nothing is eligible (the frontend then uses `representative`)."""
    if not rows:
        return None
    med = sorted(r["fraction_clen_matched"] for r in rows)[len(rows) // 2] if min_clen_fraction is None else min_clen_fraction
    def eligible(r: dict[str, Any]) -> bool:
        return (r.get("is_clan_best") in (True, None) and (r.get("experimental_method") or "").lower().startswith("x-ray") and r.get("pdbe_2d_available") is True
                and bool(r.get("organism_species")) and r.get("organism_species") != "synthetic construct"
                and r["n_match"] > 0 and r["n_match_with_coordinates"] / r["n_match"] >= min_coord_fraction and r["fraction_clen_matched"] >= med)
    cands = [r for r in rows if eligible(r)]
    if not cands:
        return None
    best = min(cands, key=lambda r: (r["resolution"] if isinstance(r.get("resolution"), (int, float)) else 99.0, -r["n_match_with_coordinates"], r["pdb_id"], r["auth_asym_id"]))
    return {**best, "selection_rule": "x-ray; PDBe 2D topology available; natural species; >=95% matched residues with coordinates; model coverage >= family median; best resolution",
            "eligible_count": len(cands), "min_clen_fraction_used": med}
