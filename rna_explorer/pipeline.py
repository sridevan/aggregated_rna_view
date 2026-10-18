"""End-to-end precomputation for one Rfam family: fetch -> instances -> unique sequences -> cmalign -> mapping -> validation -> statistics -> export."""
from __future__ import annotations
import datetime as _dt
import csv
import hashlib
import json
import logging
from collections import Counter
from pathlib import Path
from typing import Any, Optional
from . import __version__
from .cache import Cache
from .config import DEFAULT_CACHE_DIR, DEFAULT_OUT_DIR, find_cmalign
from .confidence import consensus_statistics
from .export import write_family
from .infernal import ScoreRow, cmalign_options, cmalign_version, run_cmalign
from .instances import choose_display_representative, choose_representative, instance_row
from .mapping import AlignedSequence, ChainMapping, build_chain_mapping, map_aligned_sequence
from .mmcif import atom_site_observed, entity_source, poly_seq_scheme
from .ncbi import taxon
from .pdbe import get_chain_cif, pdbe_2d_available, scheme_from_residue_listing
from .pdbe import get_entry_nucleic_mappings, get_molecules, get_nucleic_mappings, get_residue_listing, get_updated_cif, search_rfam_entities
from .rfam import get_clan_members, get_cm, get_family, get_release, get_structures, structures_by_chain
from .sequences import ChainInstance, chain_instance_from_pdbe, unique_sequences, write_fasta
from .stockholm import parse_stockholm, type_columns
from .validation import validate_chain

log = logging.getLogger("rna_explorer")


def align_uniques(cmalign: Path, cm_path: Path, fasta: Path, work: Path, cpu: int) -> tuple[dict[str, AlignedSequence], Any, dict[str, Any]]:
    sto, sfile = work / "aligned.sto", work / "aligned.scores"
    scores = run_cmalign(cmalign, cm_path, fasta, sto, sfile, cpu=cpu)
    aln = parse_stockholm(sto)
    cm = type_columns(aln.gc["RF"], aln.gc.get("SS_cons", ""))
    aligned = {name[4:]: map_aligned_sequence(name, aln.seqs[name], aln.gr.get(name, {}).get("PP", ""), cm) for name in aln.order}
    return aligned, cm, {k[4:]: v.to_dict() for k, v in scores.items()}


OVERRIDES_PATH = Path(__file__).with_name("family_overrides.json")


def apply_display_override(acc: str, rows: list[dict[str, Any]], automatic: Optional[dict[str, Any]], path: Path = OVERRIDES_PATH) -> Optional[dict[str, Any]]:
    """The display representative from `family_overrides.json` when the family has one (recorded choice with its reason), else the
    automatic rule's pick. The automatic pick is kept in the result for comparison."""
    ov = json.load(open(path)).get(acc, {}).get("display_representative") if path.exists() else None
    if not ov:
        return automatic
    row = next((r for r in rows if r["seq_name"] == ov["seq_name"]), None)
    if row is None:
        raise ValueError(f"{acc}: display_representative override {ov['seq_name']} is not among the built instances")
    return {**row, "selection_rule": f"override ({ov.get('chosen_by', 'user')}, {ov.get('date', '')}): {ov['reason']}",
            "automatic_choice": automatic["seq_name"] if automatic else None}


def read_selection(path: Path) -> list[tuple[str, str, Optional[str]]]:
    """Curated chain selection: (pdb_id, label_asym_id, auth_asym_id or None) per row. Columns `pdb_id` and `label_asym_id` (or the 16S
    curation's `ssu_label_asym_id`); `auth_asym_id` / `ssu_auth_chain` optional; other columns are metadata and ignored."""
    rows = list(csv.DictReader(open(path)))
    out: list[tuple[str, str, Optional[str]]] = []
    for r in rows:
        asym = r.get("label_asym_id") or r.get("ssu_label_asym_id")
        if not r.get("pdb_id") or not asym:
            raise ValueError(f"{path.name}: every row needs pdb_id and label_asym_id (or ssu_label_asym_id)")
        out.append((r["pdb_id"].strip().lower(), asym.strip(), (r.get("auth_asym_id") or r.get("ssu_auth_chain") or "").strip() or None))
    if len(set(out)) != len(out):
        raise ValueError(f"{path.name}: duplicate chains")
    return out


def build_family(acc: str, cache_dir: Path = DEFAULT_CACHE_DIR, out_dir: Path = DEFAULT_OUT_DIR, cpu: int = 2,
                 max_entries: Optional[int] = None, cmalign: Optional[Path] = None, selection: Optional[Path] = None) -> dict[str, Any]:
    """Map every PDBe-mapped chain of the family (default), or only the chains of a curated `selection` CSV. A selection build reads each
    chain's coordinates from the PDBe ModelServer (one chain, first model) instead of the whole-entry updated mmCIF, and rebuilds the
    validator's residue scheme from the PDBe residue listing (the chain file has no _pdbx_poly_seq_scheme); everything else is identical."""
    cache = Cache(cache_dir)
    cmalign = cmalign or find_cmalign()
    fam_dir = out_dir / "families" / acc
    work = fam_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    # --- Rfam ---
    fam = get_family(cache, acc)
    cm_info = get_cm(cache, acc)
    release = get_release(cache)
    rfam_hits = structures_by_chain(get_structures(cache, acc))
    clan = get_clan_members(cache, acc)
    clan_hits: dict[tuple[str, str], list[tuple[str, float]]] = {}
    for member in (clan or [acc]):
        for (pdb, chain), h in structures_by_chain(get_structures(cache, member)).items():
            clan_hits.setdefault((pdb, chain), []).append((member, float(h["bit_score"])))
    # --- PDBe ---
    selected = read_selection(selection) if selection else None
    mappings = {} if selected else get_nucleic_mappings(cache, acc)
    pdb_ids = sorted({p for p, _, _ in selected}) if selected else (sorted(mappings)[:max_entries] if max_entries else sorted(mappings))
    molecules = get_molecules(cache, pdb_ids)
    search = search_rfam_entities(cache, acc)
    instances: list[ChainInstance] = []
    cif_paths: dict[str, Path] = {}
    skipped: list[dict[str, Any]] = []
    listings: dict[str, Any] = {}
    for pdb, asym, _auth in (selected or []):
        listing = listings.setdefault(pdb, get_residue_listing(cache, pdb))
        segs = [g for g in get_entry_nucleic_mappings(cache, pdb).get("Rfam", {}).get(acc, {}).get("mappings", []) if g["struct_asym_id"] == asym]
        if len(segs) != 1:
            skipped.append({"pdb_id": pdb, "label_asym_id": asym, "reason": f"{len(segs)} {acc} segments for this chain in PDBe"}); continue
        seg = segs[0]
        ent = {int(e["entity_id"]): e for e in molecules.get(pdb, [])}.get(int(seg["entity_id"]))
        if ent is None:
            skipped.append({"pdb_id": pdb, "label_asym_id": asym, "reason": "entity not in molecules"}); continue
        try:
            cif = cif_paths[f"{pdb}:{asym}"] = get_chain_cif(cache, pdb, asym)
            instances.append(chain_instance_from_pdbe(pdb, seg, ent, listing, atom_site_observed(cif, asym, seg["entity_id"]), None))
        except ValueError as e:
            skipped.append({"pdb_id": pdb, "label_asym_id": asym, "reason": str(e)})
    for pdb in ([] if selected else pdb_ids):
        cif_paths[pdb] = get_updated_cif(cache, pdb)
        listing = get_residue_listing(cache, pdb)
        ents = {int(e["entity_id"]): e for e in molecules.get(pdb, [])}
        for seg in mappings[pdb]:
            ent = ents.get(int(seg["entity_id"]))
            if ent is None:
                skipped.append({"pdb_id": pdb, "segment": seg, "reason": "entity not in molecules"}); continue
            try:
                observed = atom_site_observed(cif_paths[pdb], seg["struct_asym_id"], seg["entity_id"])
                instances.append(chain_instance_from_pdbe(pdb, seg, ent, listing, observed, entity_source(cif_paths[pdb], seg["entity_id"])))
            except ValueError as e:
                skipped.append({"pdb_id": pdb, "segment": seg, "reason": str(e)})
    log.info("%s: %d entries, %d instances, %d skipped", acc, len(pdb_ids), len(instances), len(skipped))
    # --- alignment of unique mapped sequences ---
    uniques = unique_sequences(instances)
    fasta = work / "unique_sequences.fasta"
    write_fasta(uniques, fasta)
    aligned, cm, scores = align_uniques(cmalign, Path(cm_info.path), fasta, work, cpu)
    if cm.clen != cm_info.clen:
        raise RuntimeError(f"alignment CLEN {cm.clen} != CM CLEN {cm_info.clen}")
    generated = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    provenance = {"rfam_acc": acc, "rfam_release": release, "cm_name": cm_info.name, "cm_checksum": cm_info.checksum, "cm_clen": cm_info.clen, "cm_build_date": cm_info.build_date,
                  "infernal_version": cmalign_version(cmalign), "cmalign_options": cmalign_options(cpu), "cmalign_mode": "local (default)", "pipeline_version": __version__,
                  "generated": generated, "pdbe_mapping_source": f"PDBe API v2 nucleic_mappings/{acc}", "rfam_hit_source": f"rfam.org/family/{acc}/structures"}
    if selected:
        provenance |= {"pdbe_mapping_source": "PDBe API v2 nucleic_mappings/{pdb_id} (per entry)",
                       "selection": {"file": str(selection), "sha256": hashlib.sha256(Path(selection).read_bytes()).hexdigest(), "n_chains": len(selected)},
                       "coordinate_source": "PDBe ModelServer per chain (first model); validation residue scheme from the PDBe residue_listing API"}
    # --- per-chain mapping, validation ---
    chains: list[ChainMapping] = []
    validation: list[dict[str, Any]] = []
    for inst in instances:
        m = build_chain_mapping(inst, aligned[inst.sequence_key], cm, acc)
        scheme = (scheme_from_residue_listing(listings[inst.pdb_id], inst.entity_id, inst.label_asym_id, inst.auth_asym_id) if selected
                  else poly_seq_scheme(cif_paths[inst.pdb_id], inst.label_asym_id, inst.entity_id))
        sc = scores.get(inst.sequence_key)
        v = validate_chain(m, scheme, cm.clen, rfam_hits.get((inst.pdb_id, inst.auth_asym_id)), ScoreRow(**sc) if sc else None)
        chains.append(m)
        validation.append(v)
    stats = consensus_statistics(cm, aligned, chains, provenance)
    rows = [instance_row(m, search.get((m.info.pdb_id, m.info.entity_id)), rfam_hits.get((m.info.pdb_id, m.info.auth_asym_id)),
                         clan_hits.get((m.info.pdb_id, m.info.auth_asym_id), []), ScoreRow(**scores[m.info.sequence_key]) if m.info.sequence_key in scores else None,
                         acc, v["ok"]) for m, v in zip(chains, validation)]
    # species-level organism names (NCBI taxonomy; strain/subspecies ids resolved to their species) and PDBe 2D topology availability
    tax_cache: dict[str, Optional[dict[str, Any]]] = {}
    for r in rows:
        t = r.get("tax_id")
        if t is not None and str(t) not in tax_cache:
            tax_cache[str(t)] = taxon(cache, t)
        info = tax_cache.get(str(t)) if t is not None else None
        r["organism_species"] = info["species_name"] if info else None
        r["species_tax_id"] = info["species_tax_id"] if info else None
        r["organism_rank"] = info["rank"] if info else None
        r["pdbe_2d_available"] = pdbe_2d_available(cache, r["pdb_id"], r["entity_id"], r["auth_asym_id"])
    methods = Counter(r["experimental_method"] or "unknown" for r in rows)
    rep = choose_representative(rows)
    display_rep = apply_display_override(acc, rows, choose_display_representative(rows))
    sl, st = stats["sequence_level"], stats["structure_level"]
    family = {"rfam_acc": acc, "rfam_id": fam.get("id") or cm_info.name, "description": fam.get("description") or cm_info.description,
              "comment": fam.get("comment"), "rfam_url": f"https://rfam.org/family/{acc}", "clen": cm_info.clen,
              "consensus_sequence": "".join(stats["rf"]), "ss_cons": "".join(stats["ss_cons"]), "clan_members": clan,
              "curation": fam.get("curation"), "counts": {"pdb_entries": len({r["pdb_id"] for r in rows}), "instances": len(rows), "unique_mapped_sequences": len(uniques),
                                                        "instances_clan_best": sum(1 for r in rows if r["is_clan_best"] in (True, None)), "instances_validation_ok": sum(1 for r in rows if r["validation_ok"]),
                                                        "skipped_segments": len(skipped)},
              "experimental_methods": dict(methods), "representative": rep, "display_representative": display_rep,
              "coverage_summary": {"sequence_coverage_mean": _mean(sl["coverage"]), "structural_coverage_mean": _mean(st["structural_coverage"]),
                                   "positions_with_sequence_coverage": sum(1 for x in sl["n_match"] if x), "positions_with_structural_coverage": sum(1 for x in st["n_match_with_coordinates"] if x),
                                   "positions_clen": cm_info.clen},
              "skipped": skipped, "provenance": provenance}
    alignments = {k: {"bit_score": v["bit_score"], "avg_pp": v["avg_pp"], "trunc": v["trunc"], "cm_from": v["cm_from"], "cm_to": v["cm_to"]} for k, v in scores.items()}
    write_family(fam_dir, family, rows, chains, stats, validation, provenance, alignments)
    summary = {"acc": acc, "entries": len(pdb_ids), "instances": len(rows), "unique_sequences": len(uniques), "validation_ok": family["counts"]["instances_validation_ok"],
               "skipped": len(skipped), "out_dir": str(fam_dir)}
    json.dump(summary, open(fam_dir / "build_summary.json", "w"), indent=1)
    return summary


def _mean(xs: list[Any]) -> Optional[float]:
    v = [x for x in xs if isinstance(x, (int, float))]
    return round(sum(v) / len(v), 4) if v else None
