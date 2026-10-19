# RF00011 (bacterial RNase P RNA class B, 50 curated chains): family page build

Date: 2026-10-19.

The selection is the research 50-chain subset (`data/selections/RF00011_50_chains.csv`), copied from the read-only research worktree
with `label_asym_id`/`auth_asym_id` added (`data/selections/README.md`).
- **Research basis:** `RF00011_RNASEP_B_CLUSTERING_RESULTS.md`, the 50-chain section, in
  `~/Desktop/aggregated_rna_clustering/experiments/rna_structural_clustering`.
- **User decisions (2026-10-19):**
  - default structure 9OWQ A;
  - discrepancy k = 2 (stable);
  - TM k = 2 (provisional);
  - no fpocketR pockets;
  - no curated state labels in the app for now.

## Data

| Measure | Value |
|---|---|
| Chains / entries | 50 / 50, all *Geobacillus stearothermophilus* P RNA, cryo-EM 2.69–3.9 Å |
| Studies | Degenhardt et al. 2026 (Nature) conformers; Lee et al. 2026 (Nat Commun); unpublished holoenzyme ± tRNA entries from the Lee group |
| Unique mapped sequences | 4 (wild type; tetraloop mutant 9OY2/9OY3/9OY4) |
| Chains passing validation | 50 / 50; 0 skipped |
| Model length | 365; matched positions 362–364, all with coordinates |
| PDBe 2D topology | 24 of 50 chains |

This is one molecule in many states, not family diversity. The group sizes reflect the subset design (stratified selection), not
populations.

## Phase 1: mapping

**Build:** `python3 -m rna_explorer build RF00011 --cpu 4 -v --selection data/selections/RF00011_50_chains.csv`. This is the
per-chain ModelServer coordinate mode.

**Agreement with the research:** `compare_with_research.py` (results in `compare_with_research.json`) shows all 50 chains identical to
the research build `data/rf00011_50`, for every residue and every consensus column.

**Default structure:** `rna_explorer/family_overrides.json` sets it to **9OWQ A** (user, 2026-10-19).
- All chains are cryo-EM, so the automatic display rule (X-ray first) finds no chain; the coverage representative is 9OWY A.
- 9OWQ A (Lee 2026, RNA alone + pre-tRNA, closed, 2.78 Å) is the best-resolution chain with a PDBe 2D topology and a published study.
- BGSU's RF00011 reference annotation file is also on 9OWQ A.

## Phase 2: orientation, SASA, nonstandard nucleotides, annotations, interactions

**Orientation** (whole-chain superposition on 9OWQ A): 11 good, 14 moderate, 25 poor.
- **Open chains:** all 21 are poor (RMSD 7.7–11 Å). The S-domain moves (P12 by 24–42 Å), so a whole-chain fit onto the closed reference
  cannot be good. This is the conformational change, not a mapping problem.
- **Closed chains:** 4 Degenhardt conformers are also poor (4.5–6.6 Å): 9ZDV (RNA alone), and 9ZIL, 9ZIU and 9ZIY (with RnpA). These
  are probably the closed conformer subgroup the research's TM k = 4 cut isolates (peripheral variation, not interpreted).

**SASA:** 50 of 50 chains measured; 20,848 nucleotides complete, 2 incomplete. 364 positions have a value.

**Nonstandard nucleotides:** none (in-vitro transcribed RNA).

**Annotations:**
- **Reference elements** (BGSU file on 9OWQ A): 27 elements on 323 consensus positions.
- **Loops:** 1,303 (HL 500, IL 554, J3 149, J4 50, J6 50). The J6 loops are in the junction track (the J5+ fix from RF00005).
- **Bound molecules:** metal ions only (Mg²⁺ in 30 entries, Ca²⁺ in 20), so the ligand-binding track is empty. "Show ligands" in 3D is
  disabled for the default chain: no small molecule lies within 4 Å.

**Fix:** BGSU's RF00011 reference file lists 39 nucleotides with an empty element label. `annotate` failed on them with an IndexError.
- `annotation_sources.parse_reference_annotations` now skips such rows: they are no annotation.
- Test: `tests/test_annotations.py::test_parsers`.

**Interactions** (BGSU FR3D): all 50 chains annotated, 46,047 intra-chain interactions, 1,239 consensus pairs.

| Category | Consensus pairs |
|---|---|
| Base pairs | 225 |
| Stacking | 556 |
| Base–ribose | 247 |
| Base–phosphate | 211 |

The default numbering reference is 9OWQ A. With 365 positions the circular diagram is off (limit 300), so the table and 3D are shown.

## Structural groups

**Import:** `python3 -m rna_explorer clustering RF00011` reads the research run `results/rf00011_50/T362` (50 chains, 362 shared positions,
1,225 pairs). The config entry is in `rna_explorer/structure_clustering_config.json`.

| Metric | k | Groups | Validation |
|---|---|---|---|
| Geometric discrepancy | 2 | 29 (closed) / 21 (open); cut 0.387; representatives 9OWN A, 9ZI7 A | stable |
| TM-score (1 − TM) | 2 | 42 / 8; representatives 9OV3 A, 9ZF9 A | provisional |

**Discrepancy:**
- The groups match the curated closed/open labels exactly (29/29 closed, 21/21 open).
- The tetraloop mutants 9OY2/3/4 fall in the open group.
- Evidence:
  - the largest merge-height gap (0.311 → 0.463);
  - the best silhouette (0.554);
  - position-bootstrap median and 5th percentile ARI 1.00;
  - structural-element jackknife minimum ARI 1.00.

**TM-score:** k = 2 does not reproduce the closed/open split.

**No floor line:** the subset has no duplicate pairs (exact duplicate entries were collapsed during selection). The nearest empirical
reference is the 9OVB/9OVC pair (Mg²⁺ 5 vs 10 mM, discrepancy 0.069).

## App

- **Registration:** `AVAILABLE_ACCS`, `FAMILY_NAMES` ("RNase P RNA (class B)") and the home-page blurb in `frontend/src/lib/route.ts`.
- **Thumbnail:** `public/thumbnails/RF00011.webp` (9OWQ A, 417 nt).
  - The thumbnail script's 3Dmol helpers were removed with the Mol* unification. They are restored, for thumbnails only, in
    `frontend/scripts/thumbnail3dmol.ts`, so all thumbnails keep one style.
- **Family Overview viewer:** the "(none nearby)" and "(ions/additives only)" notes beside Show ligands moved into its tooltip. With
  "(none nearby)" the toolbar wrapped onto two lines for this family and for tRNA 1EHZ.
- **Git-ignored:**
  - `data/families/RF00011/loops.json` (10.2 MB);
  - `data/families/RF00011/sasa.json` (11 MB);
  - the per-chain `annotations/`.

  As for RF00177, these are regenerable and deployed but not committed.
- **Tests:** `frontend/tests/family-rf00011.spec.ts` (3); `home.spec.ts` lists six families.

## Reproduce

    python3 -m rna_explorer build RF00011 --cpu 4 -v --selection data/selections/RF00011_50_chains.csv
    python3 experiments/rf00011_rnase_p/compare_with_research.py
    python3 -m rna_explorer orient RF00011
    python3 -m rna_explorer annotate RF00011
    python3 -m rna_explorer sasa RF00011
    python3 -m rna_explorer modifications RF00011
    python3 -m rna_explorer interactions RF00011
    python3 -m rna_explorer clustering RF00011
    cd frontend && npm run thumbnails -- RF00011     # with npm run dev running
