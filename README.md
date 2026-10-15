# RNA Structure Explorer

Prototype for exploring experimentally determined RNA structures aggregated by Rfam family, using Rfam covariance-model match-state positions as the family coordinate system.

Status: **Milestone 1 (data preparation), Milestone 2 (three-tab frontend with synchronised 2D/3D/consensus selection) and Milestone 3B (annotation backend: reference elements, Motif Atlas loops, bound molecules, consensus tracks) implemented**; see `MILESTONE1_REVIEW.md`, `RNA_VISUALIZATION_REVIEW.md`, `MILESTONE2_VALIDATION.md`, `MILESTONE2_1_VALIDATION.md`, `MILESTONE2_2_VALIDATION.md`, `annotation_assessment/` (Milestone 3A) and `MILESTONE3B_VALIDATION.md`.

## Layout

| Path | Content |
|---|---|
| `rna_explorer/` | Python package: retrieval (Rfam, PDBe, mmCIF), sequences, Infernal, mapping, validation, confidence statistics, export, pipeline, CLI |
| `rna_explorer/schemas/` | JSON Schemas for `family.json`, `instances.json`, `consensus_stats.json`, `mappings/*.json`, and the annotation files (`reference_structural_elements`, `loops`, `ligands`, `tracks`, `annotation_instance`) |
| `rna_explorer/annotation_config.json` | reviewable bound-molecule classification (ligand / metal ion / additive / water), default-track classes and artefact rules |
| `annotation_assessment/` | Milestone 3A data assessment (`ANNOTATION_ASSESSMENT.md`, `PROPOSED_ANNOTATION_SCHEMA.md`, raw source cache used to seed annotation builds) |
| `tests/` | pytest suite incl. offline regression tests against the RF00167 and RF00177 proofs of concept |
| `data/families/<ACC>/` | pipeline output for the frontend (static JSON + CSV twins + validation + `work/` alignment) |
| `cache/` | cached downloads (safe to delete; rebuilt on demand) |
| `tools/infernal-1.1.5/bin/` | cmalign binary built from source (Infernal 1.1.5) |
| `poc_rf00167/`, `poc_rf00177/` (+ `bgsu_validation/`) | original proofs of concept and validation reports (unchanged) |
| `DESIGN.md` | Milestone 1 design: architecture, schemas, reuse map |
| `frontend/` | React/TypeScript/Vite app (see `frontend/README.md`, `frontend/ARCHITECTURE.md`) |
| `experiments/rna2d_poc/` | isolated PDBe RNA Viewer + PDBe Mol* experiment (1y26) |
| `docs/screenshots/` | Milestone 2 screenshots |

## Run

    pip install requests jsonschema pytest        # python >= 3.10
    python3 -m rna_explorer build RF00167 --cpu 4 -v
    python3 -m rna_explorer orient RF00167      # orientation.json: Kabsch superposition onto the display representative
    python3 -m rna_explorer annotate RF00167 --raw-dir annotation_assessment/cache   # annotations/, reference_structural_elements.json, loops.json, ligands.json, tracks.json
    python3 -m pytest -q

`cmalign` is found via `RNA_EXPLORER_CMALIGN`, then `tools/infernal-1.1.5/bin/cmalign`, then `PATH`.

## Pipeline (per family)

1. Rfam: family metadata, CM (CLEN, checksum, build date), release, `structures` hits (bit score, E-value, model span), clan members for clan competition.
2. PDBe: `nucleic_mappings/<ACC>` segments, molecules (batch), residue listings, search-API entity metadata (organism, method, resolution), updated mmCIF.
3. Chain instances: PDBe-mapped region extracted from the parent-base sequence; coordinates from `atom_site`; species from mmCIF source categories.
4. Deduplicate identical mapped sequences; one `cmalign` run (local mode) on the unique sequences.
5. Per-chain mapping: match residues get a consensus position; insert and local-end residues get an anchor (preceding consensus position, 0 before the first) and offset, never a consensus position; deleted vs uncovered consensus positions; residues outside the mapped region are `not_submitted_flank`.
6. Validation against `_pdbx_poly_seq_scheme`, coordinate cross-check, continuity, insertion bookkeeping, span agreement with Rfam.
7. Confidence statistics per consensus position: sequence level (unique sequences) and structure level (chain instances, coordinates), raw counts, PP bin midpoints, provisional display bands (>= 0.95 high, >= 0.80 medium). No composite score.
8. Export: `family.json`, `instances.json`, `consensus_stats.json`, `mappings/<pdb>_<entity>_<auth_chain>.json`, `validation.json`, CSV twins.

Annotation step (`annotate`, Milestone 3B): BGSU reference element labels (family file + per-chain transfer), BGSU RNA 3D Motif Atlas loops and motif groups (loop download CSV + motifs page), PDBe `binding_sites` bound molecules; all identifiers resolved through the mapping files to `(pdb_id, label_asym_id, label_seq_id)` and consensus positions; outputs `annotations/<seq_name>.json`, `reference_structural_elements.json`, `loops.json`, `ligands.json`, `tracks.json`, `annotation_summary.json`. `--raw-dir` seeds the cache from previously downloaded responses; `--offline` forbids network access.

Representative structure rule: clan-best instance with the most consensus-matched residues that have coordinates, then best resolution, then lowest PDB id.
