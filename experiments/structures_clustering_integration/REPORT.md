# Structures tab: structural clustering audit and integration plan

Date: 2026-10-16. This was an audit and planning task. Nothing was implemented, recomputed or modified in either worktree, in the demo project, or in any Claude Code session record. The only repository writes are this file and its two appendices:

- `APPENDIX_A_antibody_demo_audit.md`: the demo source audit, with line numbers.
- `APPENDIX_B_clustering_inventory.md`: the research-output inventory, with paths, formats and fields.

**Labels used throughout:**

| Label | Meaning |
|---|---|
| **[V]** | verified in code or data during this audit |
| **[H]** | historical intent, from session records or commit messages |
| **[S]** | scientific conclusion of the clustering reports |
| **[R]** | recommendation needing implementation |

The two screenshots supplied by the user are the visual target:

- **Kir2.2, PDB-CPX-119152:** 11 assemblies, full heatmap and dendrogram.
- **Human haemoglobin, PDB-CPX-154652:** 341 assemblies, collapsed group view.

---

## 1. Executive summary

- **The research data are good enough for the default metric.** Complete, tested pairwise matrices exist for all four benchmark families [V]:
  - geometric discrepancy D, using base centres plus FR3D base frames;
  - fixed-correspondence US-align TM-score on the same instances and positions;
  - several RMSD variants.

  The discrepancy implementation is validated against FR3D reference values [V/S].
- **The clustering outputs the tab needs are missing** [V]. No family has saved linkage matrices, leaf order, full group memberships, medoids, 1 − TM files or superposition transforms. TM-based clustering exists only as summary statistics and figures, with no stability analysis.
- **Phase 1 is a small backend command.** It reads the approved research matrices, runs average-linkage clustering, cuts it, computes medoids and writes schema-validated JSON. It needs no new similarity computation, but the tests and provenance rules of this project require it to sit in `rna_explorer/`.
- **The demo is a strong design and algorithm reference, but it cannot be dropped in** [V]. It is a different stack: React 18 JSX, plain CSS, 3Dmol.js for 3D rather than Mol*, and no tests. Its clustering runs in the browser. Its UI logic ports well:
  - average linkage with a k-cut and medoids;
  - the SVG dendrogram and collapsed group overview;
  - the groups panel and its rules;
  - figure exports.

  The 3D view must be rebuilt on the app's existing `PdbeMolstarViewer`.
- **The 50-instance threshold already exists in the demo** as `DETAIL_MAX = 50` [V]. Under it, RF00167 (76 clustered chains) uses the collapsed view. RF00059 (37), RF00162 (29) and the RF00177 subset (22) use the full heatmap.
- **3D superposition can work immediately for RF00167.** The app's `orientation.json` already superposes every instance onto the family frame by C4′ Kabsch [V]. It covers all 76 clustered chains: 48 good, 26 moderate, and 2 poor (the apo chains 5E54 B and 5SWD B). Group overlays can therefore be built without new transforms. Base-centre transforms consistent with D are a later refinement.
- **The cut level must be decided before Phase 1** (§23). The research approved the dataset and the method, but only the first split is called stable [S].

  | Family | Stable first split |
  |---|---|
  | RF00167 | apo 5E54/5SWD versus 72 others |
  | RF00059 | apo 8F4O versus 35 |

  A cut at the duplicate "floor" gives 18 groups for RF00167, 7 of them singletons.
- **Only RF00167 can be integrated without first building family data in this checkout.** RF00059 and RF00162 have production mappings only in the research worktree. RF00177 is a research subset of 22 chains, not a production family.

## 2. Verified Git worktrees and project locations [V]

| Path | Branch / HEAD | Role |
|---|---|---|
| `/Users/sria/Desktop/aggregated_rna_view` | `main` @ 9e94873 | Main application. Recent work is uncommitted. |
| `/Users/sria/Desktop/aggregated_rna_clustering` | `research/rna-clustering` @ 9e94873 | Research worktree of the same repository. Untracked: `INTERIM_REPORT_TASKS1-4.md`, `data/families/RF00059/`, `data/families/RF00162/`, `experiments/rna_structural_clustering/` (272 MB). |
| `/Users/sria/Desktop/antibody_antigen_interface` | own repository, HEAD `acd238a` "Structural groups: F-ATP synthase page" | Demo. The page is in `app-react/src/groups/` and `app-react/src/states/`. |

All research work is untracked; nothing on the research branch is committed beyond the shared initial commit.

Research root, referred to as E below: `/Users/sria/Desktop/aggregated_rna_clustering/experiments/rna_structural_clustering`. It contains:

- **Reports:** `README.md`, `CLUSTERING_RESULTS.md` (with an erratum), `CLUSTERING_RECOMMENDATIONS.md`, `GEOMETRIC_DISCREPANCY_VALIDATION.md`, `METHODS_COMPARISON.md`, `COVERAGE_ASSESSMENT.md`, `ORIENTATION_CONTRIBUTION.md`, `RF00167_OUTLIER_ANALYSIS.md`, `RF00059_TPP_RESULTS.md`, `RF00162_SAMI_RESULTS.md`, `RF00177_16S_RAD_RESULTS.md`, `MULTIFAMILY_CLUSTERING_BENCHMARK.md`, and the root-level `INTERIM_REPORT_TASKS1-4.md`.
- **Code:** `lib/rsc`, `scripts/`, `tests/`.

## 3. Review of `antibody_antigen_interface` [V unless marked]

Full detail is in Appendix A. Key facts:

- **Stack** (`app-react/package.json`): React ^18.3.1, Vite ^5.4, JavaScript/JSX and one hand-written `src/styles.css`. There is no TypeScript, no Tailwind and no tests. The only workflow deploys to GitHub Pages.
- **3D:** 3Dmol.js 2.5.2 loaded from a CDN (`index.html:8`) and used by `src/states/SuperpositionViewer.jsx`. The Mol* component in the repository is not used by this page.
- **Data:** one `public/<complex>-similarity/instance_similarity.json` per complex, holding labels, the 1 − TM matrix, the RMSD matrix, assemblies with metadata, ligands, mutations and modified residues. Pre-superposed backbone CIFs sit in `assemblies/<asm>.cif`.
  - The generator is `scripts/build_instance_similarity.py`. It uses distance = 1 − (TM₁ + TM₂)/2 (line 518), refuses incomplete matrices, and applies the US-align matrices to coordinates at build time.
- **Clustering in the browser** (`src/groups/clustering.js`):
  - average linkage (UPGMA) with a deterministic tie-break;
  - the number of groups is fixed per page in `src/main.jsx` (`nGroups` 3, 3, 6 or 2);
  - the cut height is the midpoint between the lowest merge undone and the highest kept;
  - the medoid is the member with the lowest mean distance to the others;
  - groups are numbered by size.
- **Layout** (`src/groups/StructuralGroupsApp.jsx`): four panels, plus Associated features.
  - The threshold is `DETAIL_MAX = 50` (line 50, used at line 234). It is duplicated as `MAX_LABELLED` in `DissimilarityHeatmap.jsx:76` and `HEATMAP_EXPORT_MAX` in `exportFigures.js:23`.
  - At or below 50: an HTML-table heatmap (viridis, grey diagonal, zoom) with `Dendrogram.jsx` in the left gutter (dashed cut line, clickable group bars).
  - Above 50: `GroupOverview.jsx`, with one wedge per group, a proportional bar and "cut at X".
- **Groups panel** (`GroupPanels.jsx`):
  - colours `['#E8685A','#B0327A','#8C4A1E','#3FA7D6']`, with grey for single assemblies and any group past the fourth;
  - `SMALL_GROUP = 5`;
  - three groups shown, the rest folded when there are more than four;
  - the largest group selected by default;
  - within-group and across-group mean 1 − TM and their ratio, from `BlockSummary.jsx` `summariseSelection`.
- **3D:** group members are overlaid up to `OVERLAY_MAX = 8` (`src/groups/StructuralGroupsApp.jsx:39`). Larger groups show the medoid alone. "Compare representatives" shows two medoids. There are no manual toggles (removed in `ad026fb`), and no viewer cleanup on unmount.
- **Associated features** (`src/groups/features.js`):
  - counted per PDB entry, comparing the group against all other groups pooled;
  - chips shown when the difference is at least 15 percentage points (`ENRICH_PP`);
  - ligands shown in both directions, mutations and modified residues only when enriched;
  - residue-1 substitutions excluded.

## 4. Associated Claude Code session [H]

The records are at `~/.claude/projects/-Users-sria-Desktop-antibody-antigen-interface/`: three jsonl transcripts, the largest 22 MB, plus `memory/`. They were read with read-only text search; nothing was resumed or modified. Details and citations are in Appendix A §8. Points relevant to this plan:

**Declined.** Fisher exact tests with Benjamini–Hochberg grading of feature chips. The user wanted plain per-PDB-entry counting without interpretation.

**Removed** (do not reproduce):

- the "Group N only" within-group view; a stale comment remains in `main.jsx`;
- manual 3D visibility toggles;
- the "Bound macromolecules" row.

The overlay cap was also reduced from 10 to 8 after a colour-vision check.

**Raised but not fixed** (worth fixing in the RNA version):

- percentages hide counts;
- the fixed 15-point rule lets any 1-of-3 feature through;
- two-group pages show mirror-image chips;
- absence is unexplained: shared versus unannotated;
- viridis min/max stretching compresses most of the range;
- the analysis notebook is private;
- whole-assembly TM-score cannot resolve substates.

**Not read.** Two Claude Docs written in that session (scientific questions; UX decisions).

The research worktree's session (`~/.claude/projects/-Users-sria-Desktop-aggregated-rna-clustering/`) has no `memory/` folder. Its user messages contain the approval of "Dataset B" for RF00167 (T ≥ 64, 76 chains, 57 positions, average linkage, sensitivity at T 63/65/67). They contain no decision on cut heights, TM normalisation or medoids.

## 5. Reusable frontend components and dependencies [R]

| Demo source (`app-react/src/…`) | Decision | Changes for the RNA app |
|---|---|---|
| `groups/clustering.js` (UPGMA, k-cut, medoid, numbering) | **Port the algorithm to the Python backend**; keep a small TS reader | Provenance and the backend-test rule require precomputed, schema-validated groups. Use the same tie-break and the same numbering by size so results match the demo. |
| `groups/Dendrogram.jsx` | **Adapt** (port to TSX) | Read the stored linkage and leaf order; labels are `seq_name` shown as "4FE5 B"; colours from the shared group palette; Tailwind. |
| `groups/GroupOverview.jsx` (collapsed view) | **Adapt** (TSX) | Same, with "cut at" showing the metric and its unit. |
| `states/DissimilarityHeatmap.jsx` | **Adapt** | Keep the table core, viridis, zoom and labels. Drop drag-to-drill and the demo's own measure switch, which the metric selector replaces. Add explicit min/max from the data. |
| `groups/GroupPanels.jsx` (`GroupsCard`, selected group card) | **Adapt** | "assemblies" becomes "RNA instances" plus "PDB entries". Show metric-specific within-group and across-group statistics with numerator and denominator. |
| `states/BlockSummary.jsx` `summariseSelection` | **Port into the backend** | Group statistics become precomputed fields. |
| `groups/ExportMenu.jsx`, `groups/exportFigures.js` | **Adapt** | Reuse the SVG/PNG/CSV export pattern; files use descriptive names with the Rfam accession. |
| `groups/RepresentativeComparison.jsx` | **Adapt** | The "Compare with" control drives the Mol* view, not 3Dmol. |
| `states/SuperpositionViewer.jsx` (3Dmol) | **Reimplement** on `frontend/src/components/PdbeMolstarViewer.tsx` | The app already uses Mol* with `orientation.json` transforms (`lib/orientation.ts` `toMat4`). It needs multi-structure loading, per-structure colours, a legend, the 8-structure cap and `plugin.dispose()` on unmount. |
| `groups/features.js` | **Reimplement the content; keep the layout and per-entry counting** | Mutations, UniProt data and the residue-1 rule do not apply. Use RNA ligand classes from `ligands.json`, `modifications.json`, species and method. Show counts as well as percentages, and separate "not annotated" from "source unavailable". |
| `scripts/build_instance_similarity.py` | **Do not reuse** | Its assumptions are protein-assembly specific: whole-assembly US-align, ModelServer coordinates, assembly IDs. The RNA matrices already exist. |

**Protein- and assembly-specific assumptions that must not carry over:**

- assemblies as the instance unit (RNA uses chain instances, `label_asym_id`);
- 1 − mean TM of whole assemblies as the only metric;
- UniProt mutations and residue-1 substitutions;
- PISA interfaces;
- ligands truncated to the top 6;
- features pooled per entry across assemblies without saying so.

**No new dependencies are needed.** The heatmap and dendrogram are plain SVG or HTML, and Mol* is already installed. 3Dmol.js must not be added.

## 6. RNA clustering output inventory [V]

Full detail is in Appendix B. The primary runs are approved for RF00167 and used in the reports for the other families. The `seq_name` values match the main app; for RF00167 all 76 are present in `instances.json`.

| Family | Mapped chains / entries | Run | Clustered chains / entries | Excluded (reason) | Shared positions n |
|---|---|---|---|---|---|
| RF00167 | 84 / 66 | `results/rf00167/T64` | 76 / 61 | 8 below T = 64 observed positions: 2XO1, 3FO4, 9NBC (63 each), 6P2H (55), 9NBH A–D (26 each; scaffold construct) | 57 (aptamer core) |
| RF00059 | 38 / 28 | `results/rf00059/T64` | 37 / 27 | 9C6J (13 observed) | 42 |
| RF00162 | 29 / 28 | `results/rf00162/T91` | 29 / 28 | none | 86 |
| RF00177 | 22 / 20 (research subset `E/data/rf00177_subset`) | `results/rf00177/T1481` | 22 / 20 | none within the subset | 1,475 |

- **RF00167 entry count:** the research reports say 60 entries, but the manifest has 61 distinct PDB IDs [V]. Use the computed value in production.
- **Each run directory holds:**
  - `manifest.json`: instances in matrix order, `shared_positions`, and `excluded[{seq_name, n_observed, reason}]`.
  - `matrices/*.csv`: a header row of IDs and one row per ID, to 6 decimals. The files are `disc_base`, `disc_gly`, `center_rmsd_base`, `rmsd_c4`, `rmsd_c3`, `ufix_tm_user`, `ufix_tm_L1row`, `ufix_rmsd` and `ufix_lali`.
  - `pairs_long.csv` and `per_position.npz`.
- **Other inputs:**
  - Standard-mode US-align: `results/<fam>/usalign_standard/`, which for RF00167 covers 79 chains, not 76.
  - Clustering summaries: `results/<fam>/clustering/clustering_summary.json` (cophenetic correlation, floor, k-sweep silhouette and ARI, truncated member lists, stability).
  - Literature labels: `E/data/curated/<fam>_instance_labels.csv`.
- **Defect:** `rf00177/usalign_standard/pairs_long.csv` contains only a header. The matrices are intact.
- **Validated versus exploratory:**
  - Validated: the discrepancy implementation, the output-integrity tests (43 tests in `E/tests`), and the RF00167 Dataset B choice.
  - Exploratory or provisional: every split below the first, the TM clustering, the orientation-only analyses, and RF00177.

## 7. Geometric discrepancy methodology [V/S]

- **Implementation:** `E/lib/rsc/base_frames.py:65-83` and `E/lib/rsc/discrepancy.py:23-32` (`matrix_discrepancy`, `best_transformation` with SVD Kabsch and reflections rejected).
- **Per nucleotide:** a base centre and a base rotation, fitted to the FR3D standard bases (fr3d-python commit 288f98c).
- **Correspondence:** the same sorted list of shared Rfam consensus match positions for every pair. Only chains with all n positions observed (complete base, C4′ and C3′) are clustered, so there is no missing-residue handling inside a pair. Insertions and construct extensions are excluded by construction.
- **Superposition:** Kabsch on the n base centres.
- **Formula:**
  - L = Σ‖cᵢ − (U c′ᵢ + t)‖²;
  - θᵢ = rotation angle of U·R′ᵢ·Rᵢᵀ, and O = Σθᵢ²;
  - D = √(L + O) / n, with weight 1, so 1 radian counts like 1 Å.
- **Properties:** symmetric with a zero diagonal (tested). It is a mixed-unit dissimilarity, not an RMSD.
- **Decomposition** (`ORIENTATION_CONTRIBUTION.md`): D² = D_pos² + D_ori², where D_pos = base-centre RMSD/√n.
  - D_pos alone reproduces the RF00167 partitions (ρ = 0.999, AUC 1.00); orientation is about 5 % of D² for apo/bound pairs [S].
  - The earlier claim that orientation explains D's advantage over C4′ RMSD was retracted. The advantage comes from using base centres [S].
  - Orientation adds local signal for ligand-analog base flips, for example 9LKC/9LKF versus 9LJN, and TPP U35. It is contaminated by syn/anti modelling choices at flip-prone positions 52, 57, 69, 70 and 77, so it should be reported separately and not relied on [S].
- **Comparability:** D scales with 1/√n, so values are not comparable across families (n = 57, 42, 86, 1,475). Floors, legends and cut heights must be family-specific.
- **Validation:** the implementation matches the fr3d-python test values and the FR3D Matlab 1GID data. The independent R3DMCS server check is pending (`GEOMETRIC_DISCREPANCY_VALIDATION.md`).

## 8. TM-score methodology and normalisation [V]

- **Binary:** `E/tools/bin/USalign`, version 20261012, commit 746ad8a.
- **Inputs:** single-chain PDB files with residue number = `label_seq_id`, using C3′ atoms.
- **Fixed mode** (`ufix_*`, the recommended alternative): `-mol RNA -I pair.fasta -u n`. The correspondence is the same n shared positions as D, and the aligned length equals n for every pair (tested). `ufix_tm_user` is normalised by n and is therefore symmetric. Instances and positions are identical to D.
- **Standard mode** (`usalign_standard`): `-mol RNA -a T` on whole chains, giving `tm_avg` (symmetric, average normalisation) and `tm_norm_by_row` (directional). RF00167 covers 79 chains, a different set from D.
- **No transforms:** `-m` was not used, so no rotation matrices exist except for two smoke-test pairs.
- **No 1 − TM files:** the transformation to a dissimilarity was only done in memory for clustering. The demo uses 1 − (TM₁ + TM₂)/2.

[R] Use **TM-score (fixed correspondence, normalised by the n shared positions)** as the alternative metric, with dissimilarity 1 − TM for clustering and dendrograms.

- It is symmetric, it needs no new normalisation, and it uses exactly the same instances and correspondence as D. Groups under the two metrics are therefore directly comparable.
- Show TM-score itself, where higher means more similar, in tooltips and the heatmap legend. State the dendrogram axis as "1 − TM-score". Never call it RMSD or discrepancy.
- Standard-mode TM can be offered later as a separate "whole-chain" option with its own instance set.

## 9. Metric-specific clustering readiness [V]

| Output needed per metric | D (`disc_base`) | 1 − TM-fixed |
|---|---|---|
| Complete pairwise matrix | ✓ all four families | ✓ TM matrix all four; 1 − TM is a trivial transform |
| Linkage matrix / leaf order | ✗ not saved | ✗ |
| Cut / group assignments | partial: summaries with member lists truncated to 10–12; full lists only for RF00177 | partial: summaries only |
| Medoids / representatives | ✗ | ✗ |
| Group statistics | ✗ (only global silhouette) | ✗ |
| Stability (bootstrap, leave-one-out) | ✓ (RF00167, RF00059, RF00162) | ✗ |

**Preprocessing still required** [R]:

1. Read the matrices.
2. Run scipy average linkage, the same as the research and the demo.
3. Store the linkage, leaf order and cophenetic correlation.
4. Cut at the agreed level (§23).
5. Compute groups numbered by size, medoids and statistics.
6. Do the same for 1 − TM.

TM stability analysis is optional for the first release, but the tab must say that TM groups are less validated. Group IDs under the two metrics are unrelated. On switching metric, the UI selects the new largest group.

## 10. Structural-group assignments and representatives [R]

- **Representative = group medoid** under the active metric: the member with the lowest mean distance to the other members. It is deterministic, with ties broken by sorted `seq_name`. This matches the demo and is a real clustered instance.
  - The representative can change when the metric changes; the UI must re-derive it.
  - A singleton's representative is itself.
- **The existing display representative stays as it is.** It is a family-level choice for Family Overview (4FE5 B, chosen by an explicit rule), not a group medoid.
- **Optional quality tie-break**, applied only among members within 5 % of the best mean distance: higher structural coverage, then better resolution. Document it in the JSON `method` block if adopted.
- **Group statistics:** n instances, n PDB entries, within-group mean and max dissimilarity, the mean across groups (this group versus all others) and their ratio. All carry metric units: D is dimensionless (Å/rad mix); 1 − TM is dimensionless.

## 11. 50-instance visualisation strategy [R]

- Apply `DETAIL_MAX = 50` to the number of instances in the **selected metric's clustering**. It is a visualisation limit only; clustering always uses all instances.
- **At or below 50** (RF00059 37, RF00162 29, RF00177 22): full heatmap, instance dendrogram, group bars and labels.
- **Above 50** (RF00167 76): collapsed group overview with wedges, sizes, representatives and the cut line, with no full heatmap by default.
- **Drill-down:** when the selected group has between 2 and 50 members, show "View this group's pairwise matrix". It is the submatrix in the group's own leaf order, reusing the heatmap component with a within-group dendrogram taken from the stored linkage.
  - This reintroduces in clean form what the demo removed ("Group N only"), so its revisit needs your sign-off.
  - Groups above 50 keep the summary. Members are never truncated to the first 50.
- **Mode can differ between metrics.** D and TM-fixed share the instance set, so for these four families they do not. A future whole-chain TM option could differ, for example 79 chains for RF00167.

## 12. Existing Structures table functionality [V]

Source: `frontend/src/components/InstanceTable.tsx` (87 lines) and `StructuresTab.tsx`.

- **Columns:** PDB ID, Chain (author), Title (truncated with a tooltip), Organism (species in italics), Method, Res. (Å), Region (`map_start–map_end`, nt), Seq. cov. %, Struct. cov. %, Avg PP, Rfam bits. Several have ⓘ tooltips.
- **Behaviour:**
  - free-text search, click-to-sort on every column, 10 rows per page with Previous/Next;
  - a selected row tied to the app-wide instance selection (`state/selection.tsx`), with automatic paging to the selection and a notice when the selection is hidden by the search;
  - keyboard selection.
- **Rest of the tab:** one `PdbeMolstarViewer` for the selected instance with the family-frame orientation, and three `PendingPanel` placeholders.
- **Tests:** `structures-selection.spec.ts`, `orientation.spec.ts`, `refinements.spec.ts`, `annotations.spec.ts` and `screenshots.spec.ts` use `tab-structures`, `structures-current` and `instance-search`. These test IDs must keep working.
- **Not present:** structure links and column filters other than the search.

## 13. Compact and expanded table design [R]

- **Compact mode** (default, about 60/40 with the Structural groups panel):
  - columns: RNA instance (`4FE5 B` with a medoid tag), Structural group (badge), Title (truncated), Res. (Å) and Method;
  - search, sort and pagination kept;
  - order "matching the dendrogram" as an extra sort option and the default, as in the demo;
  - non-clustered instances listed last with a grey "not clustered" badge and a tooltip giving the reason, for example "63 of 64 required positions observed". They stay visible in the family.
- **Expanded mode** ("Expand table"):
  - full width, with the groups panel hidden but its selection kept;
  - all 11 existing columns plus the group badge, with the same sort, search and page;
  - "Compact view" restores the layout;
  - the panels below (pairwise, Mol*, features) never hide.
- **Badges are interactive:** clicking one selects that group. Rows of the selected group get the group-colour left border, as in the demo.
- **Row click** keeps the existing app-wide instance selection. Initially it only highlights the row and shows that instance in the Mol* panel's "Displayed structures" list (§18).
- **Responsive:** below the `lg` breakpoint, stack the table above the groups panel. Below about 768 px, use the expanded table automatically and keep horizontal scrolling inside the table, as today.

## 14. Structural groups panel design [R]

- **Group list:** colour dot, "Group k", "n RNA instances", and "largest" / "small" tags (small = fewer than 5).
  - Show three groups and fold the rest ("N more groups · M single instances"). A selected folded group is pinned into view.
  - Palette: the demo's four colours plus grey for singletons, checked against the app's teal, orange and magenta annotation colours.
- **Selected group card:**
  - n instances and n PDB entries, the representative (medoid tag);
  - within-group mean and max, the across-group mean and the ratio, labelled per metric: "Mean D (dimensionless)" or "Mean 1 − TM-score";
  - an ⓘ explaining lower = more similar (D) and higher TM = more similar.
- **"Clustering details":** method, linkage, cut rule and value, cophenetic correlation, instances clustered and excluded with reasons, and provenance.
- **Default:** the largest group. On a metric change, select the new largest group.

## 15. Heatmap and dendrogram integration strategy [R]

- **Panel header:** "Pairwise structural dissimilarity" with an ⓘ, a metric selector (Geometric discrepancy default; TM-score) and a Download menu.
- **Colour:** viridis with explicit family-specific min/max from the data, and a legend reading "More similar ← dissimilarity → More different".
  - For TM, show 1 − TM in cells and TM in tooltips. Alternatively show TM on a reversed scale with its own legend.
  - Either way the legend states the metric and direction.
- **Interactions:** hovering a cell shows "A × B: D = 0.124 (57 positions)". Clicking a dendrogram group bar or a wedge selects the group. Selected-group outline. Zoom −/fit/+. Heatmap cells are not clickable in the first release, as in the demo.
- **Cut line:** the dashed cut at the stored height, labelled with its value.
- **Floor marker (recommended):** the family-specific duplicate floor, labelled "duplicate floor (95th pct)", as recommended by `CLUSTERING_RECOMMENDATIONS.md`.
- **Downloads:** dendrogram SVG/PNG, pairwise CSV in long format, the instance-group table CSV, and the heatmap SVG/PNG when there are at most 50 instances. Figures are drawn from data. Remove the demo's notebook link or replace it with the downloadable JSON.

## 16. Mol* superposition readiness [V → R]

- **No transforms are saved in the research outputs** [V], and the TM matrices carry none.
- **The app already has transforms onto one frame** [V]. `data/families/RF00167/orientation.json` holds a Kabsch C4′ fit of every instance onto 4FE5 B on shared consensus positions. All 76 clustered chains are included: 48 good, 26 moderate, and 2 poor (5E54 B 4.63 Å, 5SWD B 4.59 Å, which are apo outliers). `lib/orientation.ts` turns them into Mol* matrices.
- **Recommended phase 4a:** show the medoid; overlay all members for groups of 2–8, each in its own colour, with a "Displayed structures" legend; for larger groups show the medoid plus a "Show more members" picker (up to 8 total). Compare representatives overlays two medoids. Fullscreen and Centre view.
  - **Caveat for the info text:** C4′ superposition onto 4FE5 B is not the base-centre fit that defines D, so the overlay is for visual comparison only.
- **Phase 4b (optional precomputation):** for each group, a base-centre Kabsch fit of every member onto its medoid over the n shared positions, stored as `group_superpositions`. This needs coordinates only, not a new similarity run. It makes the overlay consistent with D and allows inspecting a group's spread around its medoid.
- **Viewer work in `PdbeMolstarViewer.tsx`:** it currently handles one structure per plugin.
  - It needs a multi-structure mode: load N structures, transform each, colour per structure (cartoon or ribbon of the mapped region), dispose cleanly, keep the 60 s watchdog.
  - Cap at 8 structures.
  - Lazily load coordinates from PDBe as today; there is no backend file service.
- **Compatibility:** the existing family frame keeps 3D views in Family Overview, the Structures tab and group overlays consistent.

## 17. Associated-features integration [R]

**Unit:** chain instance, the clustering unit. Compare inside the selected group with all other clustered instances. Show both counts and percentages on each chip, for example "HPA 9/30 (30%) vs 3/46 (7%)". Show no significance grading, as the user requested in the demo session.

| Feature | Source (main app) | Unit and caveat |
|---|---|---|
| Ligands (default-track classes) | `ligands.json` `bound_molecules` | Entry/entity-level PDBe annotation. Count as "chain in an entry whose entity has ligand X annotated", labelled "entry-level annotation". Never present as an observed chain contact. Chain-specific evidence is a separate future layer (`experiments/chain_specific_ligands/`). |
| Nonstandard nucleotides | `modifications.json` | Chain-specific; on-axis and off-axis components can both be listed. |
| Species | `instances.json` (`organism_species`) | Chain-level. |
| Experimental method, resolution bins | `instances.json` | Entry-level, repeated per chain. |
| Construct length, sequence and structural coverage | `instances.json` / mappings | Chain-level; show as distributions (median, IQR) rather than chips. |
| Predicted pockets at the ligand site | `pockets.json` | Chain-specific; optional, labelled as predicted. |

**Display rule:**

- show chips where the absolute difference is at least 15 percentage points **and** the group has at least 3 positive instances;
- keep "show all" so nothing is hidden;
- one-direction display on two-group pages to avoid mirror chips;
- say "none annotated" versus "source unavailable" explicitly.

**Interpretation:** no causal or state language. The panel description reads "Features of RNA instances in this group compared with the other clustered instances".

## 18. Coordinated interaction model [R]

A single `StructuresState` (React context, mirroring `state/selection.tsx`) holds:

```
{ metric: 'discrepancy' | 'tm', groupId, tableMode: 'compact' | 'expanded', sort, query, page, compareGroupId?, shownMembers[] }
```

The existing app-wide `selection.instance` stays separate.

| Event | Effect |
|---|---|
| Initial load | metric discrepancy; compact table; largest group selected; pairwise view full or collapsed by n; Mol* shows the medoid, plus members if 2–8; features for the selected group |
| Select group (groups panel, badge, dendrogram bar or wedge) | Update every panel: card, row highlight, heatmap outline, representative, Mol*, features. Do not change the sort, search or page, except jumping to the group's first row when the order matches the dendrogram. |
| Expand / compact table | Only `tableMode` changes. |
| Switch metric | Load that metric's clustering; select its largest group; clear `compareGroupId` and `shownMembers`; keep `selection.instance`. |
| Click a row | Set `selection.instance` (as today). If the instance is in the selected group and not shown, add it to Mol* (up to the cap). No multi-select and no residue-level linking. |

## 19. Proposed production JSON schemas [R]

These follow the conventions of the existing data files: `data/families/<ACC>/`, `schema_version`, `provenance`, and a JSON schema per file in `rna_explorer/schemas/`.

**`structure_clustering.json`** holds family-level metadata, the instance set and per-metric clustering. The matrices are stored once.

```json
{
  "schema_version": "1.0", "rfam_acc": "RF00167",
  "dataset": { "n_mapped_instances": 84, "n_clustered_instances": 76, "n_clustered_entries": 61,
    "inclusion_rule": "observed (complete base + C4' + C3') at >= 64 of the aptamer-core consensus positions",
    "shared_positions": [20, 21, 22],
    "excluded": [{ "seq_name": "6p2h_1_A_A", "n_observed": 55, "reason": "below threshold 64" }],
    "instances": ["1y26_1_X_A", "..."] },
  "default_metric": "discrepancy",
  "metrics": {
    "discrepancy": { "label": "Geometric discrepancy", "kind": "dissimilarity", "unit": "dimensionless (A and rad, weight 1)",
      "definition": "...", "matrix_file": "structure_matrices/discrepancy.json",
      "clustering": { "method": "agglomerative", "linkage": "average", "distance": "D",
        "linkage_matrix": [[0, 1, 0.012, 2]], "leaf_order": [5, 2],
        "cophenetic_correlation": 0.0, "cut": { "rule": "k", "k": 0, "height": 0.0 }, "floor": { "value": 0.124, "rule": "95th pct of duplicate pairs" },
        "groups": [{ "group_id": 1, "members": ["..."], "n_instances": 0, "n_entries": 0, "medoid": "4fe5_1_B_A",
          "within": { "mean": 0.0, "max": 0.0, "n_pairs": 0 }, "across": { "mean": 0.0, "n_pairs": 0 }, "small": false }],
        "assignments": { "1y26_1_X_A": 1 }, "stability": { "bootstrap_ari_median": 0.0, "source": "..." } } },
    "tm": { "label": "TM-score (fixed correspondence)", "kind": "similarity", "distance": "1 - TM", "normalisation": "n shared positions",
      "usalign": { "version": "20261012", "options": "-mol RNA -I <pair.fasta> -u n" }, "matrix_file": "structure_matrices/tm.json",
      "clustering": { "...": "same structure as above, independent group ids" } } },
  "provenance": { "research_run": "experiments/rna_structural_clustering/results/rf00167/T64", "research_commit": null, "matrices_sha256": {},
    "mapping_pipeline_version": "0.1.0", "rfam_release": "15.1", "cm_checksum": "1292453307", "generated": "..." }
}
```

**`structure_matrices/<metric>.json`** holds the values in `dataset.instances` order, upper triangle only, rounded to 4 decimals:

```json
{ "metric": "discrepancy", "instances": [], "upper_triangle": [] }
```

For RF00167 that is 2,850 values, about 25 KB, so the files are small enough to load eagerly.

**Other outputs:**

- The instance metadata stays in `instances.json`. Group assignments are read from `structure_clustering.json`, so nothing is duplicated.
- **`structure_superpositions.json`** (phase 4b): for each metric and group, the medoid as reference and, for each member, rotation, translation, n matched positions and RMSD. Until then the frontend uses `orientation.json`.
- **Associated features** are computed in the frontend from existing files, which avoids duplicating them. Alternatively a precomputed `structure_group_features.json` can list per feature: `{feature_type, id, observation_unit, in_group: {k, n}, outside: {k, n}, source}`. The precomputed form is preferred if the counting rules need backend tests.

## 20. Benchmark-family readiness [V]

| | RF00167 | RF00059 | RF00162 | RF00177 |
|---|---|---|---|---|
| Family data in this checkout | ✓ full, with annotations | ✗ mappings in the research worktree only | ✗ same | ✗ research subset only (22 E. coli 70S chains), not a production family |
| D matrix | ✓ 76 | ✓ 37 | ✓ 29 | ✓ 22 |
| TM-fixed matrix | ✓ 76 | ✓ 37 | ✓ 29 | ✓ 22 |
| Metric-specific clustering (stored) | ✗ | ✗ | ✗ | ✗ |
| Group assignments (full) | ✗ (truncated) | ✗ (truncated) | ✗ (truncated) | ✓ in the summary |
| Representatives | ✗ | ✗ | ✗ | ✗ |
| Dendrogram data | ✗ | ✗ | ✗ | ✗ |
| View mode | collapsed (76 > 50) | full (37) | full (29) | full (22) |
| Drill-down | depends on cut; with k = 2 the large group has 72 members, so summary only | n/a | n/a | n/a |
| Superposition | ✓ via `orientation.json` (C4′) | needs `orient` after a build | needs `orient` | needs a family build |
| Associated features | ✓ ligands, modifications, metadata | needs `annotate` and `modifications` | same | same |

**Integration order:**

1. RF00167 first.
2. RF00059 to validate the full heatmap mode. It needs `build`, `orient`, `annotate` and `modifications` in this checkout, using the same production commands that built it in the research worktree.
3. RF00162 next.
4. RF00177 later. It needs a production-scale strategy and a decision on the 16S subset.

## 21. Incremental implementation phases [R]

1. **Clustering data.**
   - New command `python3 -m rna_explorer clustering <ACC> --research-run <dir>`.
   - It copies nothing from the research code. It reads the research matrices and manifests, verifies the instance IDs against `instances.json` and the mappings, and checks matrix symmetry and diagonals.
   - It runs average linkage for D and 1 − TM-fixed, then cuts, computes medoids and group statistics, and writes `structure_clustering.json` plus matrices, with JSON schemas.
   - Tests: synthetic linkage and medoid cases; a regression test on the RF00167 first split (5E54/5SWD ×2 versus 72) and the 8 exclusions.
2. **Table and groups.** Compact and expanded modes, group badges, groups panel, selected-group card, "not clustered" handling, `StructuresState`. Playwright tests for selection, expansion, metric switching and existing test IDs.
3. **Pairwise visualisation.** Port the dendrogram, collapsed overview and heatmap to TSX, with legends, the cut and floor lines, drill-down for groups of 50 or fewer, downloads, and vitest and Playwright coverage. Validate full mode on RF00059 after its family build.
4. **Mol*.** Multi-structure mode in `PdbeMolstarViewer`, group overlay up to 8, medoid display, Compare representatives, legend, fullscreen and centre. Optionally (4b) per-group base-centre superpositions.
5. **Associated features.** Per-instance counting with explicit k/n, the display rule, and source and unit labels.

## 22. Estimated development effort [R]

Estimates assume the current pace and the existing test conventions.

| Phase | Effort | Main uncertainty |
|---|---|---|
| 1 Clustering data (RF00167, both metrics) | 1–1.5 days | the cut decision (§23) |
| 1′ Build RF00059 (and RF00162) family data in this checkout | 0.5 day each | network and cache; mapping review |
| 2 Table + groups panel | 1.5–2 days | integration with the existing selection and tests |
| 3 Heatmap / dendrogram / overview / drill-down / downloads | 2–3 days | porting and layout at both thresholds |
| 4a Mol* multi-structure overlay | 2–3 days | `PdbeMolstarViewer` multi-structure refactor, performance |
| 4b Base-centre group superpositions | 1 day | — |
| 5 Associated features | 1–1.5 days | ligand counting unit wording |
| **Total** | **about 9–13 working days** | |

## 23. Implementation risks and unresolved decisions

**Decisions needed before Phase 1:**

1. **The cut rule per family and metric.** Options:
   - **(a) Fixed k from the stable split.** RF00167 k = 2 gives 4 versus 72. It is scientifically safest, but the largest group is not drillable.
   - **(b) The duplicate-floor height.** RF00167 gives 18 groups, 7 of them singletons. It is detailed but includes provisional splits.
   - **(c) A k chosen by silhouette from the stored sweep.**

   The demo uses a hand-set k per page. **Recommendation:** store the full linkage and use a per-family configured cut. Start with the floor cut for display, draw the floor line, and mark groups below the stable split as "provisional" in the details. The user must confirm.
2. **The TM variant.** The recommendation is fixed correspondence, 1 − TM_n, with the same 76 chains. Standard mode has 79 chains and a different correspondence.
3. **Whether to reintroduce per-group drill-down,** which the demo removed.
4. **Whether to show the floor line and "provisional" markers.**
5. **Wording for excluded instances,** "not clustered: incomplete core", and whether they get a row badge or a footnote.

**Risks:**

- **Research code and data are uncommitted.** The 272 MB research directory exists only in an untracked worktree. Production provenance should record file hashes, and the research outputs should be committed or archived before relying on them.
- **The RF00167 entry count is inconsistent:** 61 in the manifest versus 60 in the reports.
- **D values are not comparable across families.** Legends and cut heights must be per family.
- **Two apo chains have poor C4′ orientation fits.** Overlays in their group will look worse than D implies until phase 4b.
- **`PdbeMolstarViewer` multi-structure support is new code.** It risks memory pressure and a slow first load beyond 8 structures; keep the cap.
- **TM groups lack stability analysis.** Label them as less validated.
- **RF00177 is not a production family,** and its `usalign_standard/pairs_long.csv` holds only a header.
- **Discrepancy validation is incomplete:** the independent R3DMCS check is pending.
