# Clustering research inventory (read-only audit)

Worktree: `/Users/sria/Desktop/aggregated_rna_clustering` (branch `research/rna-clustering`, HEAD 9e94873 = same initial commit as main).
Everything under audit is **untracked**: `experiments/rna_structural_clustering/`, `INTERIM_REPORT_TASKS1-4.md`, `data/families/RF00059/`, `data/families/RF00162/`.
The main checkout (`/Users/sria/Desktop/aggregated_rna_view/data/families/`) contains only RF00167. Below, `E/` = `experiments/rna_structural_clustering/`.
Nothing was modified or re-run; pytest was **not** run. Values marked (V) were verified by reading result files; (R) = taken from reports only.

## Conventions

- Instance id: the production `seq_name` = `<pdb>_<entity>_<auth_asym>_<label_asym>` (e.g. `4fe5_1_B_A` = auth chain B, label A; `1y26_1_X_A`; 16S `8bf7_52_2_ZA`). This was checked against `data/families/RF00167/instances.json`. The docstring in `E/lib/rsc/benchmark.py:2` says the same. Report prose ("4FE5/B") uses author chains.
- Residues are addressed internally by label_seq_id, through `mapping["residues"][*].polymer_seq_position`. Only `alignment_state == "match"` with `has_coordinates` counts (`benchmark.py:36`).
- Coordinates: PDBe updated mmCIF in `cache/www_ebi_ac_uk/...` (gitignored), first model and first altloc (`E/lib/rsc/coords.py`). For 16S, per-chain model-server files in `cache/model_server/<pdb>_<label_asym>.cif` take precedence (`benchmark.py:17-27`).
- Mapping provenance is copied into every manifest: Rfam 15.1, Infernal 1.1.5, pipeline 0.1.0. Software: rsc 0.1.0, numpy 2.4.2, python 3.11.9, scipy 1.15.3, sklearn 1.9.0, US-align 20261012.

## A. Per-family datasets (V from `results/<fam>/<run>/manifest.json` and `instances.json`)

| Family | Mapped chains / entries | Primary run | Clustered chains / entries / unique seqs | Excluded (reason) | Shared positions n (span) |
|---|---|---|---|---|---|
| RF00167 | 84 / 66 (`data/families/RF00167`) | `T64` (≥64 observed consensus positions; aptamer core 19–87) | 76 / **61** / 27 | 8: 2XO1, 3FO4, 9NBC (63 obs), 6P2H (55), 9NBH ×4 (26; scaffold only). Reason string `"n_observed X < 64"` | 57 (20–83) |
| RF00059 TPP | 38 / 28 (worktree `data/families/RF00059`) | `T64` | 37 / 27 / 5 | 1: 9C6J (13 obs, scaffold-focused map) | 42 (10–89) |
| RF00162 SAM-I | 29 / 28 (worktree `data/families/RF00162`) | `T91` | 29 / 28 / 23 | 0 | 86 (1–108) |
| RF00177 16S | research subset of 22 / 20 (`E/data/rf00177_subset`, built by `scripts/rf00177_subset_build.py`; the production POC has ~1817 chains) | `T1481` | 22 / 20 / 10 (R) | 0 | 1475 (5–1524) |

- **Entry count:** CLUSTERING_RESULTS.md and the Phase-2 prompt say RF00167 T64 has "60 entries". The manifest has 61 distinct PDB IDs: 66 entries minus 2XO1, 3FO4, 6P2H, 9NBC and 9NBH. This is a small inconsistency in the report.
- **Sensitivity runs (all with full matrices):**
  - RF00167: T63 (79 chains, 56 positions), T65 (74, 58; drops 5E54 A and 5SWD A), T67 (69, 59).
  - RF00059: T68 (36, 48; drops 8F4O B), T70 (28, 54; drops the thiC chains), T74 (24, 59; no apo).
  - RF00162: `Tte_only` (26, 90; 3V7E C/D and 4KQY excluded by name).
  - RF00177: T1500 (21, 1485), T1506 (16, 1505).
- **How positions are excluded** (`benchmark.py:93-99` `select`):
  - Eligibility requires `len(observed) >= T`. Chain length is not used.
  - The shared set S is the intersection of every eligible chain's **frame-complete** positions. A position is frame-complete when all parent heavy base atoms are present, plus the glycosidic N9/N1, C4′ and C3′ (`benchmark.py:43`).
  - S can optionally be restricted to a core. RF00167 uses 19–87 by default in `rf00167_benchmark.py`, because the 9NDD scaffold occupies CM 1–18 and 88–102. The generic `family_benchmark.py` runs used no core (`core_restriction: null`).
  - Positions lost to base geometry:
    - RF00167 T64: 22 (5C7W Z base) and 68 (9NBC/9NDD backbone-only U).
    - RF00059 T64: 37, 50, 56, 65, 67, 68, 76, 81 (U46/C55 bases unmodelled in most thiM entries).
    - RF00162: 9 and 14 (3V7E).
  - Per-instance coverage files: `E/results/coverage_scan/<fam>_instance_coverage.csv`, `<fam>_position_observation.csv`, `<fam>_threshold_scan.json` (keys `threshold_scan[*].{threshold_min_observed, n_eligible, excluded, n_shared_frame_complete, shared_frame_complete_positions, positions_lost_to_base_geometry,...}`).
- **US-align standard mode instance sets** (whole chains, independent of T): RF00167 79 (union at T63), RF00059 37, RF00162 29, RF00177 22.

## B. Metrics and matrix files

All matrices are CSV, written by `write_matrix`. The header is `seq_name,<id1>,...`; each row starts with the id followed by values to 6 decimals. Row and column order equal `manifest.instances` (tested). The order is the sorted order of `instances.json`, not the dendrogram order.

Per run, the directory `E/results/<fam>/<run>/matrices/` contains (V, symmetry checked):

| File | Meaning | Symmetric | Diagonal |
|---|---|---|---|
| `disc_base.csv` | BGSU geometric discrepancy D, base centres, w = 1 | yes | 0 |
| `disc_gly.csv` | D with glycosidic N9/N1 centres (R3DMCS variant; same base frames) | yes | 0 |
| `center_rmsd_base.csv` | RMSD of superposed base centres (location only). D_pos = this/√n (verified to 6e-7) | yes | 0 |
| `rmsd_c4.csv`, `rmsd_c3.csv` | Kabsch RMSD over S (Å) | yes | 0 |
| `ufix_tm_user.csv` | US-align `-I` fixed correspondence, TM normalised by user L = n (`-u n`) — **symmetric "TM-fixed"** | yes | 1 |
| `ufix_tm_L1row.csv` | same run, TM normalised by the length of the **row** structure (directional, whole chain length) | **no** | 1 |
| `ufix_rmsd.csv` | US-align RMSD over the n fixed pairs under the TM-optimal superposition | yes | 0 |
| `ufix_lali.csv` | aligned length (= n for every off-diagonal pair; diagonal 0) | yes | 0 |

- **`<run>/pairs_long.csv`** columns: `seq_name_1, seq_name_2, n_positions, disc_base, disc_gly, center_rmsd_base, rmsd_c4, rmsd_c3, max_angle_deg, max_angle_position, ufix_lali, ufix_rmsd, ufix_tm_norm_L1, ufix_tm_norm_L2, ufix_tm_norm_n, ufix_command` (upper triangle only).
- **`<run>/per_position.npz`** contains:
  - `positions` (n,)
  - `pair_a`, `pair_b` (n_pairs,), indices into `names`
  - `names` (N,)
  - `angles_deg` (n_pairs, n): per-position base rotation angle under the base-centre superposition
  - `c4_dev` (n_pairs, n): per-position C4′ deviation **after the C4′ Kabsch** superposition (not base-centre displacement)
- **`<run>/manifest.json`** contains: family, threshold, (core_restriction in family runs), n_eligible, n_excluded, `excluded[{seq_name,n_observed,reason}]`, shared_positions, instances, metric definitions, `usalign_fixed_pairs_with_lali_not_equal_n` (0 for every run), software, mapping_provenance, coordinates, generated.

US-align standard mode, `E/results/<fam>/usalign_standard/`:

| File | Meaning |
|---|---|
| `tm_avg.csv` | `-a T`, TM normalised by the **average length**; symmetric; diagonal 1 |
| `tm_norm_by_row.csv` | TM normalised by the row structure's length; asymmetric (max \|M − Mᵀ\| = 0.40 in RF00167) |
| `rmsd.csv`, `lali.csv` | RMSD and aligned length (both vary per pair) |
| `pairs_long.csv` | `L1, L2, lali, rmsd, seq_id, tm_norm_L1, tm_norm_L2, tm_norm_avg, command` |
| `manifest.json` | mode, version, instances |
| `raw/*.txt` | full US-align outputs (gitignored, present locally) |

- **Defect:** `rf00177/usalign_standard/pairs_long.csv` contains only the header. The `--reuse-standard` path overwrote it with an empty table (`family_benchmark.py:41-57`). The matrices and the 231 raw outputs still exist.

Derived distances:

- **1 − TM:** never saved as a matrix file. It is computed in memory by the clustering scripts:
  - `one_minus_tm_fixed = 1 − ufix_tm_user`
  - `one_minus_tm_std = 1 − tm_avg` (restricted to the run's instances)
  
  It appears only in `clustering_summary.json`, the `control_pairs_*.csv` files and the figures.
- **D_pos and D_ori matrices** (CSV, same format):
  - RF00167: `E/results/rf00167/orientation_ablation/{D_T64,D_pos_T64,D_ori_T64}.csv`. Verified: D_T64 equals disc_base exactly, and D² = D_pos² + D_ori² to 1e-6.
  - RF00059: `E/results/rf00059/clustering/{D_pos_T64,D_ori_T64}.csv`.
  - RF00162: `E/results/rf00162/clustering/{D_pos_T91,D_ori_T91}.csv`.
  - None for RF00177.
- **RF00177 region metrics** (CSV, 4 decimals) in `E/results/rf00177/rad_analysis/`: `body_rmsd_c4`, `body_disc`, `head_disp_bodyframe_c4`, `head_disp_bodyframe_base`, `head_disc_bodyframe`, `head_rot_angle_deg`, `head_only_rmsd_c4`, plus `per_structure_head_rotation.csv`.

## C. Geometric discrepancy as implemented

Code: `E/lib/rsc/base_frames.py`, `E/lib/rsc/discrepancy.py`, `E/lib/rsc/benchmark.py`.

**Base frame** (`base_frames.py:65-83`, `base_frame`):
- Heavy base atoms of the **parent** base (from the mapping's `nucleotide_parent`), matched by name against the FR3D standard base coordinates. These are in `E/data/reference/fr3d_standard_bases.json`, extracted from fr3d-python `NAbasecoordinates` at commit 288f98c.
- Kabsch fit (`best_transformation`, lines 44-62) of observed atoms onto the standard ones. Row-vector convention: (x − mean)·R lies in the standard frame.
- Base centre = mean of the atoms if complete; otherwise `mean_r − U·mean_s`. Glycosidic centre = N9 (A, G) or N1 (C, U).
- Needs at least 3 atoms (`MIN_ATOMS`). The benchmark additionally requires completeness (`benchmark.py:43`).

**Superposition** (`best_transformation`):
- SVD Kabsch on the centred sets. A reflection is rejected when det ≈ −1.
- Returns U, mean1, mean2, rmsd and sse. The convention is (set1 − mean1)·U ≈ set2 − mean2.
- The superposition is unweighted and done on the per-nucleotide **centres only**; base frames are not used for it.

**Discrepancy** (`discrepancy.py:23-32`, `matrix_discrepancy`):
- U = Kabsch(centres1 → centres2).
- L = sse of the superposed centres.
- θ_i = `angle_of_rotation(U @ R2_i @ R1_iᵀ)` (lines 19-20; arccos((tr − 1)/2), clipped).
- O = Σθ_i² in radians.
- **D = √(L + w·O)/n with w = 1**: 1 rad is weighted like 1 Å. The docstring gives units as "Å per nucleotide", which is effectively a mixed Å/rad quantity.
- Asserts n > 2; the FR3D n = 2 formula is not implemented.
- Returns `discrepancy, n, location_sse, center_rmsd, orientation_sse, angles_rad, rotation`. **The rotation U is returned but discarded** by `fixed_metrics` (`benchmark.py:144-152`). The translation (the means) is not returned at all.

**Decomposition:** D² = D_pos² + D_ori², with D_pos = √L/n = center_RMSD/√n and D_ori = √O/n. This is implemented in `scripts/rf00167_orientation_ablation.py:3-5,35` and `family_clustering.py:144-147`.

**Correspondence:**
- Pairs are the same consensus positions in both chains, in the same sorted order of S (`features`, `benchmark.py:138-141`). The same S is used for every pair in a run, so n is constant within a matrix.
- Insertions, local ends, unmapped residues and positions outside S are never used.
- Missing residues or atoms cannot occur inside S by construction: S requires frame-completeness in every eligible chain. Chains lacking a position are excluded by T, or the position is dropped from S.
- `check_mapping_consistency` (`benchmark.py:65-78`) asserts that the mmCIF comp_id at the mapped label_seq_id equals the mapping's residue_name, and that the frame parent equals the mapping parent. It runs before every run (`assert not problems`).

**Comparison frame:** none is persisted. D is frame-invariant: each pair is superposed on its own centroids.

**Modified residues:** handled through the parent's standard atom names. IU and CCC are kept; 50N and 50L have no parent and therefore no frame.

**n-dependence:** D ∝ 1/√n, so values cannot be compared across families or runs (`discrepancy.py:10-11`; METHODS_COMPARISON.md §1.5).

**Validation** (GEOMETRIC_DISCREPANCY_VALIDATION.md; `tests/test_discrepancy.py`):

| Check | Result |
|---|---|
| fr3d-python test, 1S72 case 1 (expected 0.9402) | 0.94017 |
| fr3d-python test, 1S72 case 2 (expected 1.1414) | 1.14143 |
| FR3D Matlab 1GID base centres (316 nucleotides) | max Δ 4e-14 Å |
| FR3D Matlab 1GID rotation matrices | max Δ 6.4e-8 |
| Self-discrepancy | 0 |
| Symmetry D(a,b) = D(b,a) | to 1e-9 |
| 1/√2 scaling with duplicated lists | exact |

- **Not yet done:** an independent R3DMCS server comparison (the comparison page returned an error) and a BGSU `map_across_chains` vs Rfam correspondence check for RF00167.
- **Documented deviations:** stricter completeness; glycosidic variant also computed; unit weights only.

## D. US-align

**Binary:** `E/tools/bin/USalign` (local build artefact, gitignored). Built by `E/tools/build_usalign.sh` from pylelab/USalign commit `746ad8a2a4dea9f317075cc482df0ecb5440bdea`, version string "20261012", with `c++ -std=c++11 -O3 -ffast-math`.

**Inputs:** single-chain PDB files written by `write_pdb` (`benchmark.py:102-116`) to `results/<fam>/pdb/<seq_name>.pdb` (gitignored). Each file has chain A, residue number = label_seq_id and one-letter parent residue names. Only residues with C3′ are written.

**Standard mode** (`family_benchmark.py:46`):
- Command: `USalign a.pdb b.pdb -mol RNA -a T`.
- Default C3′ representative atom; sequence-independent structural alignment; whole chains.
- No `-mm` or `-ter` options (single chain per file). No `-TMscore` option (default TM-align/RNA-align search).
- Normalisations kept: L1, L2 and the average (`-a T`).
- Unequal lengths: the aligned subset differs per pair; length and construct differences enter the score. 9NDD is the first standard-mode split in RF00167.

**Fixed correspondence ("TM-fixed")** (`family_benchmark.py:94-96`):
- Command: `-mol RNA -I <pair>.fasta -u n`.
- The FASTA, built by `fixed_alignment_fasta` (`benchmark.py:119-135`), aligns only the shared consensus positions; every other residue is gapped.
- US-align keeps this alignment and optimises the superposition for TM-score over those pairs. The output literally prints "User-specified initial alignment", but the aligned length equals n for every pair in every run (manifest counter = 0; tested).
- Normalisations kept:
  - user n (`ufix_tm_user`; symmetric, d0 from n, e.g. 2.01 Å at n = 57)
  - L1 and L2 (directional, whole chain lengths)
- RMSD is over the n pairs under the TM-optimal (not Kabsch) superposition.

**Coverage vs D:** the TM-fixed matrices cover exactly the same instances and positions as D. Standard mode covers a slightly larger set for RF00167 (79 vs 76).

**Not saved:** `-m` (rotation matrix) and `-do` were **not** used in the benchmark runs. They were used only in the smoke tests `E/results/usalign_smoke/4fe5B_vs_{5e54A,2eesA}/{standard,fixed}_full.txt`, which contain matrices for those two pairs. `-outfmt 2` was used only in the smoke tests; benchmark outputs are parsed from the default text output by regex (`parse_usalign`, `benchmark.py:155-175`).

**Recommendation in the reports:** keep TM-score as a sequence- or fold-level descriptor and an optional mapping-independent cross-check, not as the state metric. Rationale: d0 saturation ranks class above state (RF00167 AUC 0.76 fixed / 0.85 standard).

## E. Clustering

**Method:** `scipy.cluster.hierarchy.linkage(squareform(M), "average")` (UPGMA) on each dissimilarity.

- Scripts: `scripts/rf00167_clustering.py` (RF00167, T64/63/65/67) and `scripts/family_clustering.py` with `data/curated/rf00059_config.json` / `rf00162_config.json`.
- RF00177 is clustered inside `scripts/rf00177_rad_analysis.py`.
- No preset k. Cut sweeps use `fcluster(Z, k, "maxclust")`: k = 2–10 for RF00167, 2–8 for the family script, 2–? for RF00177 (not checked).
- There is one distance cut, the **noise-floor cut** at height = 95th percentile of the duplicate floor. The floor is defined as NCS copies plus same sequence and same ligand, both bound.

**Metrics clustered** (`dissimilarities()` in both scripts): disc_base, disc_gly, rmsd_c4, ufix_rmsd, **one_minus_tm_fixed**, **one_minus_tm_std**, rmsd_std. The family script also clusters D_pos and D_ori. **TM-based clustering therefore exists**, as summary statistics and figures, alongside D.

**Outputs per metric** in `E/results/<fam>/clustering/clustering_summary.json`:
- Location: `thresholds.<T>.metrics.<metric>` for RF00167, `runs.<tag>.metrics.<metric>` for the family script.
- `cophenetic_correlation`
- `noise_floor{n,median,p95,max}`, `noise_floor_cut{height,n_clusters,singletons}`
- `by_pair_category`, `all_pairs`
- `partitions.<k>{silhouette, ARI_state, ARI_sequence_group, (ARI_class), clusters{id:{n,states,groups,members}}}`
- `outliers_top5`
- `positive_vs_bound_separability` (family script)

Also present: metric agreement (Spearman plus Mantel p), partition ARI between metrics, threshold or run sensitivity ARI, stability, localisation, and control pairs.

**Defect — member lists are truncated** (12 in rf00167, `rf00167_clustering.py:143`; 10 in the family script, `family_clustering.py:111`). Full group assignments are therefore **not saved** for RF00167, RF00059 or RF00162. They are deterministic and recomputable from the matrices. The RF00177 summary stores full member lists (`rad_analysis_summary.json:clustering.<metric>.<k>`).

**Not saved anywhere:** linkage matrices Z, leaf order (`leaves_list` is used only for heatmap figures), and medoids or representatives (none computed).

**Figures (PNG):**
- `clustering/figures/dendrogram_<T|run>_<metric>.png` and `heatmap_..._<metric>.png` for disc_base, rmsd_c4, one_minus_tm_fixed and one_minus_tm_std
- `metric_scatter_*.png`, `per_position_*.png`
- RF00177: `rad_analysis/figures/dendrogram_{disc_base,rmsd_c4,one_minus_tm_fixed,one_minus_tm_std,head_disp_bodyframe_c4,head_rot_angle_deg}.png`
- RF00167 orientation ablation: `orientation_ablation/figures/dendrogram_{D,D_pos,D_ori}.png`

**Stability** (D and rmsd_c4 only; **none for TM**):
- Element jackknife (leave one structural element out; ARI at k = 2/3/4).
- Position bootstrap: 100 resamples for RF00167, 50 for the family script.
- RF00167 also has k = 4 co-clustering matrices: `clustering/coclustering_bootstrap_k4_{disc_base,rmsd_c4}.csv`.

**Labels** (literature, kept separate from the clusters), in `E/data/curated/<fam>_instance_labels.csv`:
- Columns: seq_name, pdb_id, label_asym_id, auth_asym_id, state, state_note, ligand, (ligand_class), riboswitch_class, sequence_group, sequence_key, mutations / organism, method, resolution, release_year, n_observed.
- RF00167 states: bound 66, apo 3 (5E54 A/B, 5SWD B), intermediate_bound 1 (5SWD A), ligand_free_holo_like 2 (3IVN), transitional_bound 4 (6VWT/6VWV).
- RF00059: bound 35, apo 2 (8F4O A/B), scaffold_focused_map 1.
- RF00162: bound 28, ligand_free_holo_like 1 (3IQP).
- RF00177: `rf00177_primary_set.csv` with RAD angles.

**Key partitions (disc_base, average linkage), V from the summaries:**

| Family | k = 2 | k = 3 / 4 | Silhouette | Floor cut | Cophenetic |
|---|---|---|---|---|---|
| RF00167 T64 | {5E54 A, 5E54 B, 5SWD A, 5SWD B} vs 72 | k3: +{3IVN A/B, 9NDD ×4}; k4: +3DS7 B singleton | 0.59 / 0.46 / 0.44 | height 0.124 → 18 clusters, 7 singletons | 0.985 |
| RF00059 T64 | {8F4O A, B} vs 35 | k3: +thiC (8); k4: +2HOP singleton | 0.84 / 0.58 / 0.54 | 0.127 → 7 clusters, 4 singletons | 0.998 (R) |
| RF00162 T91 | {3V7E C, D} vs 27 | k3: +4KQY; k4: +{5FK4, 5FJC, 2YGH, 4B5R} | 0.57 / 0.56 / 0.35 | 0.060 → 19 clusters, 13 singletons (floor has n = 2; the report says it is unusable) | — |
| RF00177 T1481 | swivelled 7 (RAD 13.1–20°) vs 15 | next: 8BF7 | — | — | — |

TM-fixed partitions differ from D:
- RF00167 TM-fixed k = 2: Guanine-II + xanthine-II (17) vs 59. Floor cut 0.131 → 14 clusters, 5 singletons.
- RF00167 TM-std k = 2: 9NDD ×4 vs rest.

## F. Superposition outputs

- **No rotation or translation matrices are saved for any pair or onto any reference.** U is computed per pair in `matrix_discrepancy` and in the Kabsch calls, then discarded. US-align ran without `-m`, except the two smoke-test pairs.
- The RF00177 body-core superposition onto 4V9D-BA is computed in `rf00177_rad_analysis.py:60-92`. Only the angles, signs and displacements are kept.
- Residue correspondence is **implicit**: the shared consensus positions are in `manifest.shared_positions`, and each chain's label_seq_id comes from its validated mapping. Per-pair `-I` FASTA files exist in `results/<fam>/<run>/usalign_fixed/*.fasta` (gitignored, present locally).
- No per-instance table of label_seq_id per shared position is written (unverified whether any helper does).
- The production `orientation.json` (main repo) is a separate, existing Kabsch-to-reference output and is not used by this research.

## G. Validated vs provisional; key conclusions

**Validated:**
- The discrepancy reimplementation against FR3D and fr3d-python.
- Output integrity tests: symmetry, identity order, lali = n, and recomputation of random pairs.
- The user approved Dataset B on 2026-10-16. Session record, verbatim intent: "Minimum 64 … all 76 eligible … intersection of 57 consensus positions … average-linkage … sensitivity (63, 65 and 67) … Do not create a separate maximum-coverage dataset".

**Provisional** (the reports' own words):
- CLUSTERING_RESULTS.md: "Provisional conclusion for RF00167 … Cluster boundaries below the apo/bound split remain provisional."
- INTERIM: "Cross-family picture so far (provisional)".
- MULTIFAMILY marks findings as (C) computational, (L) literature or (I) interpretation; its §5 recommendation is "not implemented".

**Key conclusions:**
- **Position vs orientation:** "The apo/bound separation and the clustering are carried by the position term … orientation term contributes ~5 % of D²". The advantage comes from the base centre, not the orientation term. An erratum was appended to CLUSTERING_RESULTS.md. Orientation helps only for analog-induced base flips (9LKC/9LKF, TPP 3D2G vs 3D2V) and is noisy at loop bases: flip-prone positions 52, 57, 69, 70, 77; position 69 is syn/anti ambiguous.
- **Floors** (duplicate 95th percentile of D):

  | Family | Floor | Notes |
  |---|---|---|
  | RF00167 | 0.124 | n = 27 pairs; resolution-stratified: ≤2 Å ≤ 0.07, 2.5–3 Å up to 0.24 |
  | RF00059 | 0.127 | n = 15 pairs |
  | RF00162 | 0.060 | n = 2 pairs; use near-wt bound max 0.071 or k-turn p95 0.167 instead |
  | RF00177 | ≈ 0.011–0.016 | isomorphous pairs |

  Quote: "thresholds must be derived per family from its own duplicates".
- **n-dependence:** "discrepancy values scale as 1/√n … rules out cross-family thresholds".
- **RF00167:** D separates every apo/bound pair from every bound/bound pair (min 0.531 vs max 0.479, AUC 1.00). TM AUC is 0.76 (fixed) and 0.85 (standard). The secondary splits are unstable (bootstrap k3 median ARI 0.70, p05 0.15).
- **TPP:** apo 8F4O vs bound D 1.01–1.13 vs ≤ 0.33; all metrics find it.
- **SAM-I:** correctly no global split (AUC 0.63); 3IQP vs 3IQR D 0.102 localised to A46.
- **16S:** all whole-chain metrics give ρ ≈ 0.82 with |Δ head rotation|. The body-frame head angle gives ρ 0.89–0.91, and r = 0.989 per structure vs RAD. Recommendation: region-frame views for large RNAs.
- **Recommendations for the app** (CLUSTERING_RECOMMENDATIONS.md §7–8; MULTIFAMILY §5):
  - Fixed-correspondence base-centre discrepancy as the default, with L and O stored separately.
  - Per-position profile alongside every pair.
  - Floor line on dendrograms; cut sweep, never a fixed k; label only stable splits.
  - Colour by sequence group; literature state as an overlay.
  - US-align standard mode as an optional cross-check.
  - "Discrete 'state' groups should not be shown without literature labels."

## H. Gaps for a production Structures tab

1. Nothing is in the production backend. Everything is research scripts. RF00059 and RF00162 mappings exist only in the worktree; RF00177 is only a 22-chain research subset (`E/data/rf00177_subset`), not a production family. Every output is uncommitted.
2. No saved linkage matrix, leaf order, full partition membership (RF00167/59/162 truncated), medoids or representatives, or explicit cut heights per family beyond the floor-cut statistic.
3. No saved superposition transforms (needed for Mol* superposition), and no common reference frame. US-align runs lacked `-m`.
4. No 1 − TM matrix files. TM clustering exists only as summaries and figures, with no TM stability analysis. Standard-mode TM covers a different instance set (RF00167 79).
5. Per-position base-centre displacement under the D superposition is not stored. `c4_dev` uses the C4′ Kabsch superposition; angles use the base-centre superposition.
6. JSON outputs are not schema-validated or in the `data/families/<ACC>/` frontend layout. The manifests carry provenance but not CM checksum per run consistently (RF00167 manifests do). No numerator/denominator fields beyond n.
7. Open items:
   - R3DMCS server cross-check
   - BGSU correspondence check for RF00167
   - Reduced-core T55 run for 6P2H/9NBC (proposed, not run; only `coverage_scan/discrepancy_smoke_T55.json`)
   - Optional exclusion of P1 positions 20–21 (proposed)
   - Resolution-stratified floor (proposed)
   - 5C7W treatment (kept)
8. Data defect: `rf00177/usalign_standard/pairs_long.csv` contains only the header.
9. The RF00167 report says 60 entries where the data give 61.
10. Scaling: US-align standard mode on 16S took about 1 h for 231 pairs (R), so the full RF00177 (~1800 chains) is impractical without subsetting.

## I. Tests (`E/tests/`; not run in this audit)

The reports state 40 pass (RF00167 + discrepancy); a count of the parametrised cases gives 43 including the family tests.

- **`test_discrepancy.py`** (5):
  - fr3d-python expected values (2 cases)
  - self = 0 and symmetry
  - 1/√2 scaling
  - 1GID frames vs FR3D Matlab (centre and rotation < 1e-3; skipped if atoms are missing)
- **`test_rf00167_benchmark_outputs.py`** (35):
  - Shared positions equal the coverage scan and lie within 19–87 (T63/64/65/67).
  - Matrix properties for 7 matrices × 4 T: ids equal the manifest order, shape, no NaN, symmetric, diagonal 0/1, TM in (0,1], lali = n.
  - US-align standard tm_avg symmetric with diagonal 1.
  - Mapping consistency, plus 25 random pairs recomputed (disc_base, rmsd_c4 < 1e-5).
  - Excluded instances recorded and complete.
- **`test_family_benchmark_outputs.py`** (3, parametrised over RF00059 T64, RF00162 T91, RF00177 T1481):
  - Same matrix checks.
  - Selection reproduces the manifest.
  - Mapping consistency.
  - 10 random disc_base pairs recomputed.
- **Not tested:** clustering outputs, the orientation ablation, the RAD analysis, and US-align standard pairs_long.

## Session record

`~/.claude/projects/-Users-sria-Desktop-aggregated-rna-clustering/` contains only `a9cad1d0-….jsonl`, `subagents/` and `tool-results/`; there is no `memory/`. The user messages contain:
- the task prompts
- the Dataset B approval (2026-10-16T03:09Z)
- the Phase 2 prompt, which asked TM-score to be "investigate[d for] an appropriate transformation into structural dissimilarity"

No user decision on cut heights, a TM normalisation choice, or medoids was found.

Related memory: `~/.claude/projects/-Users-sria-Desktop-aggregated-rna-view/memory/clustering-investigation.md` agrees with the reports.
