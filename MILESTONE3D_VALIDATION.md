# Milestone 3D validation: RNA solvent accessibility (SASA) annotations

Date: 2026-10-16. New module `rna_explorer/sasa.py`, CLI `python3 -m rna_explorer sasa RF00167` (optional deps `pip install freesasa gemmi numpy`, listed as `[sasa]` extra), schemas `sasa.schema.json` and `sasa_tracks.schema.json`, outputs `data/families/RF00167/sasa.json` (3.8 MB, per-instance nucleotides) and `sasa_tracks.json` (149 KB). Pilot: `experiments/sasa_pilot/PILOT_REPORT.md`.

## Method (documented in every output under `method`)
Isolated RNA chain (polymer atoms of the instance's label_asym_id, first model of the PDBe updated mmCIF; no proteins, other chains, ligands, ions or water); heavy atoms; altloc = highest occupancy per atom name; element radii C 1.70 / N 1.55 / O 1.52 / P 1.80 Å (Bondi, P as NACCESS; S, Se, halogens defined); FreeSASA 2.2.1 Shrake–Rupley via `calcCoord`, probe 1.4 Å, 1000 points; nucleotide SASA = sum of its atoms. Completeness from the wwPDB CCD (leaving atoms excluded): complete | incomplete (missing atoms listed, e.g. 5′ phosphate) | unsupported (no CCD or unknown element; never a zero radius) | no_coordinates (never zero SASA). Independent numpy Shrake–Rupley check on 1Y26: 0.06 % total, ≤ 1.1 % per nucleotide.

## Aggregation (`sasa_tracks.json`)
One measurement per RNA chain instance per Rfam consensus match position; only `complete` nucleotides contribute; insert/local-end/flank residues stay per instance; statistics median (primary), mean, min, max, sd, Q1, Q3 (numpy linear percentiles), contributing-instance count, mapped-with-coordinates count, excluded counts by reason, nucleotide identity distribution, contributing instance list. No normalisation by instance count; absolute Å².

RF00167: 84/84 instances measured, 6,133 complete / 46 incomplete / 422 no-coordinate nucleotides, 0 unsupported; all 102 positions have a median (median of medians 171.7 Å², max 307.1 Å²); consensus 53: 80 instances, median 134.8 Å² (U 78, C 2). Colour domain 0–350 Å² (max median rounded up to 50, stored in the file).

## Frontend
New collapsible category **Solvent accessibility** with one track **Median SASA (Å²)**: fixed-size marks (same geometry as loops/ligands), continuous light→dark sequential scale over 0–350 Å² (`ScaleDef` kind `quantity`; a legend gradient labelled in Å² is generated), grey marks for positions without a qualifying measurement, measured zero would be the lightest colour. Tooltip: "Rfam residue N / Median SASA: x Å² / Based on k RNA instances"; pinned details add mean/min/Q1/Q3/max/sd, counts and nucleotide identities. Category ⓘ uses the requested text. Download menu gains "Solvent accessibility, per-nucleotide (JSON)" and "Solvent accessibility, aggregated (JSON)" (`<ACC>_sasa.json`, `<ACC>_sasa_tracks.json`, HEAD-checked; a family without SASA shows neither the category nor working downloads). Loop and ligand tracks untouched.

## Tests
- `tests/test_sasa.py` (synthetic gemmi structure): chain isolation (protein chain and water excluded), hydrogens removed, altloc policy, isolated-atom area = 4π(r+1.4)² within 1 %, unsupported element flagged (not zero), completeness/incomplete/no_coordinates statuses, sum identity, aggregation (one per instance, exclusion of incomplete, insert not aggregated, median/quartiles/sd, reproducibility, duplicate-position guard), CCD parsing.
- `tests/test_sasa_rf00167.py`: schema validity, 84 measured, status counts, finite/non-negative, sum identity per chain, no zero for missing residues, consensus 53 values, denominator identities, ordered quantiles.
- Playwright `sasa-track.spec.ts`: category and single track, mark geometry equal to other tracks, colour monotonic with median, Å² legend with domain, missing vs zero colours, axis alignment, tooltip text and counts, pinned details, category popover (no track icon), both downloads valid for the family, loop/ligand tracks unchanged. Screenshot `docs/screenshots/m3d_annotations_sasa.png`.

## Limitations
RNA-only accessibility of the deposited chain (no biological assembly, no protein or ligand context); 5′-terminal residues usually lack the phosphate and are therefore excluded from the primary track; modified nucleotides contribute with their own chemistry (identity recorded); Shrake–Rupley with 1000 points has ~1 % per-nucleotide numerical noise; no relative SASA.
