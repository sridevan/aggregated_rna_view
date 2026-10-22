"""Per-nucleotide solvent-accessible surface area (SASA) of isolated experimental RNA chains, aggregated on Rfam consensus positions.

Method
  * Coordinates: the PDBe updated mmCIF already used by the pipeline (cache), first model only. The chain is selected by label_asym_id
    (gemmi subchain), so only the polymer atoms of that RNA chain enter the calculation: no protein, other chains, ligands, ions or water.
  * Atoms: heavy atoms only (H/D removed). Alternate locations: for each atom name within a residue the location with the highest
    occupancy is kept (ties: first altloc label in alphabetical order); every other altloc copy is dropped, never double counted.
  * Radii: element-based van der Waals radii (Bondi 1964; P 1.80 as in NACCESS): C 1.70, N 1.55, O 1.52, P 1.80, S 1.80, Se 1.90,
    F 1.47, Cl 1.75, Br 1.85, I 1.98. Atoms of any other element make the nucleotide 'unsupported' (never a silent zero radius).
  * Algorithm: FreeSASA Shrake-Rupley, probe radius 1.4 Å, 1000 test points per atom (freesasa.calcCoord: coordinates + radii only,
    so FreeSASA's protein-oriented residue classifiers are not involved).
  * Nucleotide SASA = sum of its atoms' SASA; chain total = sum over nucleotides (identity holds by construction and is tested).
  * Completeness: expected heavy atoms per residue come from the wwPDB Chemical Component Dictionary (PDBe pdbechem_v2 CIF, RCSB ligand CIF only as fallback; cached),
    excluding atoms with pdbx_leaving_atom_flag = Y (OP3). Status: complete | incomplete (missing heavy atoms listed; the 5'-terminal
    phosphate is the usual case) | unsupported (CCD unavailable or atom element without radius) | no_coordinates (residue in the
    sequence but not in the model). Only 'complete' nucleotides enter the primary aggregated track.
  * Aggregation: one measurement per RNA chain instance per consensus match position; insert / local-end / flank residues are kept
    per instance and never aggregated; statistics are median (primary), mean, min, max, sd (n>1), quartiles (numpy linear).
"""
from __future__ import annotations
import datetime as _dt
import json
import logging
from collections import Counter
from pathlib import Path
from typing import Any, Optional
import numpy as np
from . import __version__
from .annotation_mapping import InstanceIndex
from .cache import Cache
from .export import validate_json
from .pdbe import family_uses_chain_files, instance_cif

log = logging.getLogger("rna_explorer.sasa")
SCHEMA_VERSION = "1.0"
PROBE_RADIUS = 1.4
N_POINTS = 1000
RADII: dict[str, float] = {"C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80, "SE": 1.90, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}
#: wwPDB Chemical Component Dictionary definition: PDBe first, RCSB only when PDBe cannot provide it (same atoms, leaving flags and parents for
#: all 41 components used so far, checked 2026-10-22)
CCD_URLS = ("https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{ccd}.cif", "https://files.rcsb.org/ligands/view/{ccd}.cif")


def method_record() -> dict[str, Any]:
    import freesasa, gemmi   # noqa: E401
    try:
        from importlib.metadata import version
        fs_version = version("freesasa")
    except Exception:   # pragma: no cover
        fs_version = "unknown"
    return {"software": {"freesasa": fs_version, "gemmi": gemmi.__version__, "numpy": np.__version__, "rna_explorer": __version__},
            "algorithm": "Shrake-Rupley (FreeSASA calcCoord)", "probe_radius_A": PROBE_RADIUS, "n_test_points": N_POINTS,
            "radii_A": RADII, "radii_source": "Bondi (1964) van der Waals radii; P 1.80 (NACCESS)", "hydrogens": "excluded",
            "altloc_policy": "per atom name within a residue keep the highest-occupancy alternate location (ties: first altloc label); other copies dropped",
            "chain_isolation": "polymer atoms of the selected label_asym_id only (gemmi subchain); no other chains, ligands, ions or water",
            "model": "first model of the deposited/updated mmCIF", "completeness_reference": "wwPDB CCD heavy atoms minus pdbx_leaving_atom_flag=Y atoms (PDBe pdbechem_v2 CIF; RCSB ligand CIF only when PDBe is unavailable)",
            "coordinate_source": "PDBe updated mmCIF (same files as the mapping pipeline)"}


# ----------------------------------------------------------------------------- chemical component dictionary
class CcdExpectedAtoms:
    def __init__(self, cache: Cache) -> None:
        self.cache = cache
        self._memo: dict[str, Optional[dict[str, Any]]] = {}

    def get(self, ccd: str) -> Optional[dict[str, Any]]:
        """heavy atoms expected in a polymer residue of this component (leaving atoms removed) or None when the CCD is unavailable."""
        if ccd in self._memo:
            return self._memo[ccd]
        import gemmi
        res = None
        for url in (u.format(ccd=ccd) for u in CCD_URLS):   # PDBe first; RCSB only when PDBe cannot provide a readable definition
            try:
                block = gemmi.cif.read_string(self.cache.fetch_text(url, ".cif")).sole_block()
                tab = block.find("_chem_comp_atom.", ["atom_id", "type_symbol", "pdbx_leaving_atom_flag"])
                atoms = [(gemmi.cif.as_string(r[0]), gemmi.cif.as_string(r[1]).upper(), gemmi.cif.as_string(r[2])) for r in tab]
                if not atoms:
                    raise ValueError("no _chem_comp_atom rows")
                heavy = [a for a in atoms if a[1] not in ("H", "D")]
                res = {"ccd_id": ccd, "expected_atoms": sorted(a for a, _, lv in heavy if lv != "Y"), "leaving_atoms": sorted(a for a, _, lv in heavy if lv == "Y"),
                       "elements": {a: e for a, e, _ in heavy}, "type": gemmi.cif.as_string(block.find_value("_chem_comp.type") or "?"),
                       "parent": gemmi.cif.as_string(block.find_value("_chem_comp.mon_nstd_parent_comp_id") or "?"), "source_url": url}
                break
            except Exception as e:   # network failure, unknown component, parse error: try the next source
                log.warning("CCD %s unavailable from %s: %s", ccd, url, e)
        if res is None:
            log.warning("CCD %s unavailable from every source", ccd)
        self._memo[ccd] = res
        return res


# ----------------------------------------------------------------------------- structure
def load_chain(cif_path: Path, label_asym_id: str) -> dict[str, Any]:
    """Heavy atoms of one RNA chain (first model) with the altloc policy applied. Returns residues in sequence order."""
    import gemmi
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    model = st[0]
    sub = model.get_subchain(label_asym_id)
    if len(sub) == 0:
        raise ValueError(f"no subchain {label_asym_id} in {cif_path.name}")
    ent = st.get_entity_of(sub)
    residues = []
    for res in sub:
        by_name: dict[str, list[Any]] = {}
        for atom in res:
            if atom.element.name.upper() in ("H", "D"):
                continue
            by_name.setdefault(atom.name, []).append(atom)
        atoms = []
        n_alt_dropped = 0
        for name, copies in by_name.items():
            if len(copies) > 1:
                copies = sorted(copies, key=lambda a: (-a.occ, a.altloc))
                n_alt_dropped += len(copies) - 1
            a = copies[0]
            atoms.append({"name": name, "element": a.element.name.upper(), "x": a.pos.x, "y": a.pos.y, "z": a.pos.z, "occupancy": a.occ, "altloc": a.altloc if a.altloc != "\0" else "", "b_iso": a.b_iso})
        residues.append({"label_seq_id": res.label_seq, "author_residue_number": res.seqid.num, "insertion_code": res.seqid.icode.strip(), "residue_name": res.name,
                         "atoms": atoms, "n_altloc_copies_dropped": n_alt_dropped})
    return {"pdb_id": st.name.lower(), "label_asym_id": label_asym_id, "n_models_in_file": len(st), "entity_type": ent.entity_type.name if ent else None,
            "polymer_type": ent.polymer_type.name if ent else None, "residues": residues}


def compute_atom_sasa(residues: list[dict[str, Any]]) -> tuple[list[float], list[str]]:
    """Per-atom SASA for all supported atoms of the chain (one FreeSASA run). Returns (areas aligned with atom order, unsupported elements)."""
    import freesasa
    coords: list[float] = []
    radii: list[float] = []
    unsupported: list[str] = []
    for res in residues:
        for a in res["atoms"]:
            r = RADII.get(a["element"])
            if r is None:
                unsupported.append(a["element"])
                a["radius"] = None
                continue
            a["radius"] = r
            coords.extend([a["x"], a["y"], a["z"]]); radii.append(r)
    params = freesasa.Parameters({"algorithm": freesasa.ShrakeRupley, "probe-radius": PROBE_RADIUS, "n-points": N_POINTS})
    result = freesasa.calcCoord(coords, radii, params) if radii else None
    areas = [result.atomArea(i) for i in range(result.nAtoms())] if result else []
    # write back
    k = 0
    for res in residues:
        for a in res["atoms"]:
            if a["radius"] is None:
                a["sasa"] = None
            else:
                a["sasa"] = areas[k]; k += 1
    return areas, unsupported


# ----------------------------------------------------------------------------- per instance
def measure_instance(idx: InstanceIndex, cache: Cache, inst: dict[str, Any], ccd: CcdExpectedAtoms, chain_files: bool = False) -> dict[str, Any]:
    cif = instance_cif(cache, inst, chain_files)   # SASA is computed for the isolated chain either way
    chain = load_chain(cif, inst["label_asym_id"])
    warnings: list[str] = []
    if chain["polymer_type"] != "Rna":
        warnings.append(f"polymer type {chain['polymer_type']} (expected Rna)")
    areas, unsupported_elements = compute_atom_sasa(chain["residues"])
    if unsupported_elements:
        warnings.append(f"atoms with unsupported elements excluded from the calculation: {dict(Counter(unsupported_elements))}")
    by_seq = {r["label_seq_id"]: r for r in chain["residues"]}
    mapping = idx.mapping(inst["seq_name"])
    nts: list[dict[str, Any]] = []
    status_counts: Counter = Counter()
    for mres in mapping["residues"]:
        i = int(mres["polymer_seq_position"])
        base = {"label_seq_id": i, "author_residue_number": mres["author_residue_number"], "insertion_code": mres["author_insertion_code"] or "",
                "residue_name": mres["residue_name"], "nucleotide_parent": mres["nucleotide_parent"], "is_modified": bool(mres["is_modified"]),
                "alignment_state": mres["alignment_state"], "rfam_consensus_pos": mres["rfam_consensus_pos"], "pp_char": mres["posterior_prob_char"] or ""}
        r = by_seq.get(i)
        if r is None or not r["atoms"]:
            nts.append({**base, "status": "no_coordinates", "sasa_A2": None, "n_observed_heavy_atoms": 0, "n_expected_heavy_atoms": None, "missing_atoms": None, "extra_atoms": None})
            status_counts["no_coordinates"] += 1
            continue
        if r["residue_name"] != mres["residue_name"]:
            warnings.append(f"residue name mismatch at label_seq_id {i}: structure {r['residue_name']} vs mapping {mres['residue_name']}")
        exp = ccd.get(r["residue_name"])
        observed = sorted(a["name"] for a in r["atoms"])
        unsupported_here = [a["element"] for a in r["atoms"] if a["radius"] is None]
        if exp is None or unsupported_here:
            status = "unsupported"
            missing = None if exp is None else sorted(set(exp["expected_atoms"]) - set(observed))
            extra = None if exp is None else sorted(set(observed) - set(exp["expected_atoms"]) - set(exp["leaving_atoms"]))
        else:
            missing = sorted(set(exp["expected_atoms"]) - set(observed))
            extra = sorted(set(observed) - set(exp["expected_atoms"]) - set(exp["leaving_atoms"]))
            status = "complete" if not missing else "incomplete"
        sasa = None if status == "unsupported" and unsupported_here else float(sum(a["sasa"] for a in r["atoms"] if a["sasa"] is not None))
        nts.append({**base, "status": status, "sasa_A2": round(sasa, 3) if sasa is not None else None, "n_observed_heavy_atoms": len(observed),
                    "n_expected_heavy_atoms": len(exp["expected_atoms"]) if exp else None, "missing_atoms": missing, "extra_atoms": extra,
                    "n_altloc_copies_dropped": r["n_altloc_copies_dropped"],
                    "incomplete_reason": ("5prime_phosphate_absent" if missing and set(missing) <= {"P", "OP1", "OP2"} else "missing_heavy_atoms") if status == "incomplete" else None,
                    "unsupported_reason": (f"elements without radius: {sorted(set(unsupported_here))}" if unsupported_here else "chemical component definition unavailable") if status == "unsupported" else None,
                    "ccd_type": exp["type"] if exp else None})
        status_counts[status] += 1
    # structure residues not in the mapping (should not happen: mapping covers the whole polymer)
    extra_struct = [i for i in by_seq if i not in {n["label_seq_id"] for n in nts}]
    if extra_struct:
        warnings.append(f"{len(extra_struct)} structure residues absent from the mapping: {extra_struct[:5]}")
    chain_total = float(sum(a for a in areas))
    return {"seq_name": inst["seq_name"], "pdb_id": inst["pdb_id"], "entity_id": inst["entity_id"], "label_asym_id": inst["label_asym_id"], "auth_asym_id": inst["auth_asym_id"],
            "source_file": cif.name, "n_models_in_file": chain["n_models_in_file"], "polymer_type": chain["polymer_type"],
            "n_atoms_used": len(areas), "chain_total_sasa_A2": round(chain_total, 3), "sum_nucleotide_sasa_A2": round(float(sum(n["sasa_A2"] or 0 for n in nts)), 3),
            "status_counts": dict(status_counts), "warnings": warnings, "nucleotides": nts}


# ----------------------------------------------------------------------------- aggregation
def _stats(values: list[float]) -> dict[str, Any]:
    v = np.asarray(values, dtype=float)
    return {"median": round(float(np.median(v)), 2), "mean": round(float(v.mean()), 2), "min": round(float(v.min()), 2), "max": round(float(v.max()), 2),
            "sd": round(float(v.std(ddof=1)), 2) if len(v) > 1 else None, "q1": round(float(np.percentile(v, 25)), 2), "q3": round(float(np.percentile(v, 75)), 2), "n": int(len(v))}


def aggregate(measurements: list[dict[str, Any]], clen: int) -> list[dict[str, Any]]:
    positions = []
    for pos in range(1, clen + 1):
        contrib: list[tuple[str, float, str]] = []
        n_mapped_coords = 0; excl: Counter = Counter()
        for m in measurements:
            seen = False
            for n in m["nucleotides"]:
                if n["rfam_consensus_pos"] != pos or n["alignment_state"] != "match":
                    continue
                if seen:   # one measurement per instance per position (defensive; a match position occurs once per chain)
                    continue
                seen = True
                if n["status"] == "no_coordinates":
                    continue
                n_mapped_coords += 1
                if n["status"] == "complete" and n["sasa_A2"] is not None:
                    contrib.append((m["seq_name"], n["sasa_A2"], n["residue_name"]))
                else:
                    excl[n["status"]] += 1
        rec: dict[str, Any] = {"rfam_consensus_pos": pos, "n_contributing_instances": len(contrib), "n_mapped_instances_with_coordinates": n_mapped_coords,
                               "n_excluded_incomplete": excl.get("incomplete", 0), "n_excluded_unsupported": excl.get("unsupported", 0),
                               "nucleotide_identity": dict(Counter(c[2] for c in contrib)), "contributing_instances": sorted(c[0] for c in contrib)}
        rec.update({"median_sasa_A2": None, "mean_sasa_A2": None, "min_sasa_A2": None, "max_sasa_A2": None, "sd_sasa_A2": None, "q1_sasa_A2": None, "q3_sasa_A2": None})
        if contrib:
            s = _stats([c[1] for c in contrib])
            rec.update({"median_sasa_A2": s["median"], "mean_sasa_A2": s["mean"], "min_sasa_A2": s["min"], "max_sasa_A2": s["max"], "sd_sasa_A2": s["sd"], "q1_sasa_A2": s["q1"], "q3_sasa_A2": s["q3"]})
        positions.append(rec)
    return positions


# ----------------------------------------------------------------------------- build
def build_sasa(acc: str, cache_dir: Path, out_dir: Path, instances: Optional[list[str]] = None, write: bool = True) -> dict[str, Any]:
    fam_dir = out_dir / "families" / acc
    idx = InstanceIndex(fam_dir)
    chain_files = family_uses_chain_files(fam_dir)
    cache = Cache(cache_dir)
    ccd = CcdExpectedAtoms(cache)
    rows = [r for r in idx.instances if not instances or r["seq_name"] in set(instances)]
    generated = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    measurements = []
    failures = []
    for r in rows:
        try:
            measurements.append(measure_instance(idx, cache, r, ccd, chain_files))
            log.info("%s total %.1f A2 %s", r["seq_name"], measurements[-1]["chain_total_sasa_A2"], measurements[-1]["status_counts"])
        except Exception as e:
            log.error("%s failed: %s", r["seq_name"], e)
            failures.append({"seq_name": r["seq_name"], "error": str(e)})
    mp = idx.mapping(idx.instances[0]["seq_name"])["provenance"]
    prov = {"rfam_acc": acc, "rfam_release": mp.get("rfam_release"), "cm_checksum": mp.get("cm_checksum"), "cm_clen": idx.clen, "mapping_pipeline_version": mp.get("pipeline_version"),
            "mapping_generated": mp.get("generated"), "generated": generated, "subset": sorted(instances) if instances else None}
    sasa = {"schema_version": SCHEMA_VERSION, "rfam_acc": acc, "clen": idx.clen, "measurement": {"name": "solvent-accessible surface area of the isolated RNA chain", "unit": "A^2", "level": "nucleotide (sum of its heavy atoms)"},
            "method": method_record(), "ccd_components": {k: v for k, v in ccd._memo.items()},
            "counts": {"n_instances_requested": len(rows), "n_instances_measured": len(measurements), "n_failed": len(failures),
                       "nucleotide_status_counts": dict(sum((Counter(m["status_counts"]) for m in measurements), Counter())),
                       "n_nucleotides_with_altloc_copies_dropped": sum(1 for m in measurements for n in m["nucleotides"] if n.get("n_altloc_copies_dropped"))},
            "failures": failures, "instances": measurements, "provenance": prov}
    positions = aggregate(measurements, idx.clen)
    medians = [p["median_sasa_A2"] for p in positions if p["median_sasa_A2"] is not None]
    tracks = {"schema_version": SCHEMA_VERSION, "rfam_acc": acc, "clen": idx.clen, "track_id": "sasa_median", "category": "Solvent accessibility",
              "measurement": {"name": "Median SASA", "unit": "A^2", "statistic": "median over contributing RNA chain instances", "context": "isolated RNA chain (RNA-only accessibility)"},
              "aggregation_rules": {"contribution": "one nucleotide measurement per RNA chain instance per Rfam consensus match position", "eligibility": "status 'complete' (all CCD heavy atoms observed; leaving atoms excluded)",
                                    "excluded": "incomplete and unsupported nucleotides are counted but not aggregated; insert/local-end/flank residues never aggregated; residues without coordinates never counted as zero",
                                    "normalisation": "none (absolute A^2)"},
              "colour_domain_suggestion_A2": {"min": 0.0, "max": float(np.ceil(max(medians) / 50.0) * 50) if medians else None, "rule": "0 to the maximum median rounded up to the next 50 A^2"},
              "summary": {"n_positions_with_value": len(medians), "median_of_medians_A2": round(float(np.median(medians)), 2) if medians else None, "max_median_A2": round(max(medians), 2) if medians else None,
                          "n_instances_measured": len(measurements)},
              "positions": positions, "method": sasa["method"], "provenance": prov}
    if write:
        validate_json(sasa, "sasa"); validate_json(tracks, "sasa_tracks")
        fam_dir.mkdir(parents=True, exist_ok=True)
        with open(fam_dir / "sasa.json", "w") as f: json.dump(sasa, f, indent=1)
        with open(fam_dir / "sasa_tracks.json", "w") as f: json.dump(tracks, f, indent=1)
    return {"acc": acc, "n_measured": len(measurements), "n_failed": len(failures), "status_counts": sasa["counts"]["nucleotide_status_counts"], "n_positions_with_value": len(medians),
            "max_median_A2": tracks["summary"]["max_median_A2"], "sasa": sasa, "tracks": tracks}
