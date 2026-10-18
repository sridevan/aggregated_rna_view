# Aggregated RNA 3D View

PDBe-KB-style aggregation of experimentally determined RNA structures per Rfam family, using Rfam CM match-state positions as the family coordinate system. Development family: RF00167 (purine riboswitch).

Full status, decisions, commands and open questions: **`PROJECT_HANDOVER.md`** (read it first). Milestone reports: `MILESTONE1_REVIEW.md`, `MILESTONE2*_VALIDATION.md`, `annotation_assessment/ANNOTATION_ASSESSMENT.md`, `MILESTONE3B_VALIDATION.md`, `MILESTONE3C_VALIDATION.md`, `MILESTONE3D_VALIDATION.md`, `MILESTONE3E_VALIDATION.md`, `MILESTONE3F_VALIDATION.md` (nonstandard nucleotides; completes the Annotations milestone), `MILESTONE4_STRUCTURES_VALIDATION.md` (RF00167 structural clustering, Structures tab; plan in `experiments/structures_clustering_integration/REPORT.md`). How chains are mapped to Rfam positions, with limitations: `docs/RFAM_PDB_MAPPING.md`.

## Layout
- `rna_explorer/` Python backend (`build`, `orient`, `annotate` commands; JSON schemas in `rna_explorer/schemas/`; ligand classification in `annotation_config.json`).
- `data/families/<ACC>/` generated JSON consumed by the frontend (symlinked from `frontend/public/families`).
- `frontend/` React 19 + TypeScript + Vite + Tailwind; pdbe-molstar (3D) and vendored PDBe RNA Viewer (2D).
- `tests/` pytest; `frontend/tests/` Playwright; vitest beside the sources.
- `annotation_assessment/` Milestone 3A assessment and the raw source cache that seeds `annotate`.

## Commands
    python3 -m rna_explorer build RF00167 --cpu 4 -v
    python3 -m rna_explorer orient RF00167
    python3 -m rna_explorer annotate RF00167 --raw-dir annotation_assessment/cache
    python3 -m rna_explorer sasa RF00167            # needs freesasa, gemmi, numpy
    python3 -m rna_explorer pockets RF00167 --experiment-dir experiments/fpocketr_full_RF00167   # needs tools/fpocketr-env (x86-64 under Rosetta)
    python3 -m rna_explorer modifications RF00167 --seed-dir experiments/modifications_RF00167/raw   # cached mmCIF; PDBe/NAKB optional
    python3 -m rna_explorer clustering RF00167          # structure_clustering.json from the research matrices (read-only ~/Desktop/aggregated_rna_clustering)
    python3 -m pytest -q
    cd frontend && npm run dev   # npm test, npm run test:e2e

## Rules
- Consensus positions come only from the validated mappings in `data/families/<ACC>/mappings/`; never re-align or force insertions, local ends or unmapped residues onto consensus positions.
- Every aggregated number carries an explicit numerator and denominator; distinguish "not annotated by the source" and "source unavailable" from biological absence; never divide by the total instance count indiscriminately.
- Keep provenance (source URL, retrieval date, Rfam/CM/Motif Atlas releases, pipeline and config versions) on every output; retain flagged or excluded observations in the data.
- Internal identifiers are `label_asym_id` + `label_seq_id`; author numbering is for display and BGSU unit ids only.
- Stay within the requested milestone scope; backend changes need tests and schema validation, frontend changes need Playwright/vitest coverage.
- Git repository on `main` (no remote yet). Commit only when asked. `.gitignore` excludes scratch `work/` folders, downloaded structure files under `experiments/` (never store full mmCIFs) and regenerable RF00177 outputs over 10 MB plus its per-chain `annotations/`.
