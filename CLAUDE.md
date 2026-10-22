# Aggregated RNA 3D View

PDBe-KB-style aggregation of experimentally determined RNA structures per Rfam family, using Rfam CM match-state positions as the family coordinate system. Development family: RF00167 (purine riboswitch).

Full status, decisions, commands and open questions: **`PROJECT_HANDOVER.md`** (read it first). Milestone reports: `MILESTONE1_REVIEW.md`, `MILESTONE2*_VALIDATION.md`, `annotation_assessment/ANNOTATION_ASSESSMENT.md`, `MILESTONE3B_VALIDATION.md`, `MILESTONE3C_VALIDATION.md`, `MILESTONE3D_VALIDATION.md`, `MILESTONE3E_VALIDATION.md`, `MILESTONE3F_VALIDATION.md` (nonstandard nucleotides; completes the Annotations milestone), `MILESTONE4_STRUCTURES_VALIDATION.md` (RF00167 structural clustering, Structures tab; plan in `experiments/structures_clustering_integration/REPORT.md`). How chains are mapped to Rfam positions, with limitations: `docs/RFAM_PDB_MAPPING.md`. UX refinements after the UX review (annotation units and denominators, missing vs absent, cross-tab action, the ligand-denominator correction in `tracks.json` `by_ccd_entries`, compact comparison header): `UX_REFINEMENTS_VALIDATION.md`. Structures tab Structural Dissimilarity mode (TPPL seriation precomputed in `structure_clustering.json`, vendored reference `rna_explorer/ordering_similarity.py`): `STRUCTURAL_DISSIMILARITY_VALIDATION.md`. Interactions tab (BGSU FR3D base pairs and stacking on Rfam positions, eligibility-based conservation, status in the displayed structure, two-chain comparison with base-anchored local superposition; basis `experiments/rf00167_interactions/`): `INTERACTIONS_TAB_VALIDATION.md`. Mol* viewer molecular context (Complete-structure layers, Show ligands with a 4 Å spatial association): `MOLSTAR_CONTEXT_VALIDATION.md`. Structures tab and Interactions comparison moved from 3Dmol.js to the same Mol* engine (shared `components/molstar/MolstarSuperposition.tsx`): `MOLSTAR_UNIFICATION_VALIDATION.md`. All tooltip text: `docs/TOOLTIPS.md`. Ligands tab (RNA-contacting ligand gallery with PDBe CCD depictions, ligand-contact residue and fpocketR pocket tracks moved from Annotations into a shared `TrackBrowser`): `LIGANDS_TAB_VALIDATION.md`; data basis `experiments/ligand_annotation_investigation/REPORT.md`. BGSU merged loops supersede their two parts in the loop tracks (inferred from nucleotide sets; all records kept): `MERGED_LOOPS_VALIDATION.md`. Sixth family RF00011 (RNase P class B, 50 curated cryo-EM chains, default 9OWQ A by override): `experiments/rf00011_rnase_p/REPORT.md`.

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
    python3 -m rna_explorer interactions RF00167 --seed-dir experiments/rf00167_interactions/raw   # interactions.json + interactions_by_instance/ (BGSU FR3D; seeded, cached)
    python3 -m rna_explorer interactions RF00177          # 16S; per-chain files git-ignored (> 10 MB); same for RF00059, RF00162, RF00005
    python3 -m pytest -q
    cd frontend && npm run dev   # npm test, npm run test:e2e

## Rules
- Consensus positions come only from the validated mappings in `data/families/<ACC>/mappings/`; never re-align or force insertions, local ends or unmapped residues onto consensus positions.
- Every aggregated number carries an explicit numerator and denominator; distinguish "not annotated by the source" and "source unavailable" from biological absence; never divide by the total instance count indiscriminately.
- Keep provenance (source URL, retrieval date, Rfam/CM/Motif Atlas releases, pipeline and config versions) on every output; retain flagged or excluded observations in the data.
- Internal identifiers are `label_asym_id` + `label_seq_id`; author numbering is for display and BGSU unit ids only.
- Stay within the requested milestone scope; backend changes need tests and schema validation, frontend changes need Playwright/vitest coverage.
- Git repository on `main`, remote `origin` = git@github.com:sridevan/aggregated_rna_view.git. Commit and push only when asked. The clustering worktree (`~/Desktop/aggregated_rna_clustering`, branch `research/rna-clustering`) shares this repository: never run history-rewriting commands that touch its branch. `.gitignore` excludes scratch `work/` folders, downloaded structure files under `experiments/` (never store full mmCIFs) and regenerable RF00177 outputs over 10 MB plus its per-chain `annotations/`.
