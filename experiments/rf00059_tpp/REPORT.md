# RF00059 (TPP riboswitch) — family page build

## Phase 1: Rfam mapping (2026-10-17)

Command: `python3 -m rna_explorer build RF00059 --cpu 4 -v`. The model is Rfam 15.1 "TPP" (TPP riboswitch aptamer, THI element), CLEN 105. The research worktree reports the same build: 38 chains, 28 entries, 5 unique sequences.

| Measure | Value |
|---|---|
| PDB entries | 28 |
| RNA chains | 29 *E. coli*, 8 *A. thaliana*, 1 recorded as *O. iheyensis*; 38 in total |
| Unique mapped sequences | 5 |
| Chains passing validation | 38 / 38 |
| Skipped | 0 |
| Log warnings | 0 |

All chains have PDBe 2D topology. 36 are X-ray and 2 cryo-EM.

| Group | Chains | Mapping pattern |
|---|---|---|
| *E. coli* thiM | 27 | CM 6–93, 77 match, 3-nt insert GCG after consensus 24 |
| 3K0J E/F | 2 | Same as thiM, but an 8-nt engineered insert (UUGCACUC) at 24 |
| *A. thaliana* thiC | 8 | CM 10–89, 68 match, no inserts |
| 9C6J, 9C6K (cryo-EM, TPP riboswitch embedded in a group IIC intron) | 2 | CM 1–105, 94 match |

**Ligand-free structure:** 8F4O A/B (apo, 3.1 Å). Coordinates exist for 71 and 64 of 77 positions.

**Coordinate coverage:** many *E. coli* chains lack some base coordinates; for example, 2HOP has 70 of 77. This is the base-geometry limitation the research describes.

**Flags:**
- 9C6J A has coordinates at only 13 of 94 mapped positions, from a scaffold-focused map. The research excludes it from clustering at T = 64.
- 9C6J's species is recorded as *Oceanobacillus iheyensis*, while 9C6K's is *E. coli*. This is the source's entity organism for a chimeric intron–riboswitch construct; it is kept as the source states and flagged, not corrected.

**Reference:** 2GDI X, which is both the display and the coverage representative (X-ray 2.05 Å, *E. coli*, 77 of 105 matched).

## Phases 2–3: orientation, SASA, modified nucleotides, annotations, pockets (2026-10-17)

All commands exited cleanly.

**Orientation** (reference 2GDI X):

| Fit quality | Chains | Which |
|---|---|---|
| good | 25 | *E. coli* thiM |
| moderate | 10 | the 8 *A. thaliana* thiC homologs (2.0–2.1 Å), 2HOP (2.42 Å, pyrithiamine), 4NYG (2.36 Å) |
| poor | 2 | 8F4O A (10.7 Å) and 8F4O B (7.7 Å), the apo chains |
| insufficient | 1 | 9C6J A (13 positions with coordinates) |

The apo chains fit poorly because P5 and P3 are displaced (Lee et al. 2023). Their 3D placement against the reference frame is offset accordingly, and the Structures tab shows a "poor fit" warning chip for them.

**SASA:**
- 38 of 38 chains measured; 3,172 nucleotides complete, 132 incomplete, 612 without coordinates.
- 90 of 105 positions have a value; the uncovered positions are outside every construct.

**Modified nucleotides:** 24 records, all terminal or flank artefacts of in-vitro transcription.

| CCD | Records | Location |
|---|---|---|
| GTP | 8 | 5′-terminal, consensus 6 |
| CCC | 8 | 3′-terminal cytidine 2′,3′-cyclic phosphate, consensus 93 |
| A23 | 8 | outside the mapped region |

**Annotations:**
- The BGSU named-element file is unavailable (source has none), recorded as such.
- 284 loops.
- 75 bound-molecule records; additives PG4 and PGE are classified as additives by the generic config, and no config change was needed.
- Ligands: TPP (9 entries) plus 17 analogues and fragments: VIB, GMI, TPS, BFT, 218, PYI, D2X, 2QB, 2QC, SVN, HPA, KWU, KX3, KXC, KWL.
- Metals: MG, MN, CA, K, NA, BA and IRI.
- Entries without a ligand: 8F4O (apo) and 9C6J (intron scaffold).
- **TPP-contact positions** (at least half of the TPP entries): 15, 33, 35, 36, 61, 62, 76, 77, 81, 82.

**Pockets** (38 chains, 264 s): 37 chains with pockets (66 in total); none in 8F4O B.
- The top-ranked pocket sits on the TPP site: positions 15, 16, 32, 33, 35 and 36 in 92 % of chains, and 77 in 86 %.
- The P5 side of the site (81, 82) is rarely covered.

## Structural-group candidates

Computed from research run T64: 37 chains, 42 positions, 9C6J excluded. Script: `experiments/rf00162_structural_groups/candidates.py RF00059 T64`. Output: `structural_groups/candidates.json`.

**Discrepancy:**
- k = 2 separates the apo pair 8F4O A/B from 35 bound chains. Its merge-height gap is 0.828 (0.257 → 1.085), the largest by far; its silhouette is 0.843, the best.
- Research bootstrap ARI at k = 2: median 1.0, 5th percentile 1.0. The T64 versus T68 run sensitivity is 1.0.
- k = 3 adds the 8 *A. thaliana* thiC homologs, but on a gap of only 0.027, with a bootstrap 5th percentile of 0.18.

**TM-score:**
- k = 2 gives the same apo split, with the largest gap (0.121) and the best silhouette (0.540). The T64 versus T68 run sensitivity is 1.0.
- There is no bootstrap for TM.

**Floor:** 0.127 for discrepancy and 0.208 for 1 − TM, each from 15 duplicate pairs. It is stored in the per-run summary and is now read by the importer.

## Structural groups (user decision, 2026-10-17)

The user chose k = 2 for both metrics: geometric discrepancy "stable", TM-score "provisional". The groups are apo 8F4O A/B versus 35 bound chains. The rationales are in `rna_explorer/structure_clustering_config.json`; the output is `data/families/RF00059/structure_clustering.json`, which passes schema validation.

| Measure | Discrepancy | TM-score |
|---|---|---|
| Chains clustered | 37 | 37 |
| Excluded | 9C6J (13 of 64 positions) | 9C6J |
| Shared positions | 42 | 42 |
| Gap at the cut | 0.828 | 0.121 |
| Cophenetic correlation | 0.998 | 0.906 |
| Duplicate floor (15 pairs) | 0.127 | 0.208 |

**Importer fix:** the per-run summary layout stores the floor, cophenetic correlation and per-k partitions under `runs.<primary run>.metrics`. The importer now reads them; earlier it treated the floor as absent. RF00162 keeps no floor line through an explicit `omit_floor` setting with the research's reason; its 2-pair value, 0.060, stays in the evidence. RF00167 is unchanged.

## Front end and tests

**Front end:** RF00059 is added to the available families (`src/lib/route.ts`, "TPP riboswitch") and the 404 page. No other change was needed.

**Browser review of every tab:** no page errors.
- The "+3" insertion marker shows at 24|25.
- The duplicate floor line is drawn.
- The apo group's two chains carry the "poor fit" warning, because their displaced arms make the superposition poor.
- 9C6J is explained as not clustered.
- The analogue and fragment ligand rows show; additives do not.

**Tests:** `tests/test_structure_clustering.py::test_rf00059_apo_vs_bound_groups` (backend, 60 passing) and `tests/family-rf00059.spec.ts` with 3 browser tests.
