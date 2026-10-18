# Milestone 4: RF00167 structural clustering in the Structures tab

Date: 2026-10-16. Scope: RF00167 only.

The plan is in `experiments/structures_clustering_integration/REPORT.md`. This milestone integrated precomputed research matrices and adapted the reference demo's components. No structural comparison was recomputed. The research worktree and `antibody_antigen_interface` were only read.

## Backend

Command:

    python3 -m rna_explorer clustering RF00167 [--research-root ~/Desktop/aggregated_rna_clustering]

**Files:**

- Module: `rna_explorer/structure_clustering.py`.
- Cut configuration with the rationale for each cut: `rna_explorer/structure_clustering_config.json`.
- Schema: `rna_explorer/schemas/structure_clustering.schema.json`.
- Output: `data/families/RF00167/structure_clustering.json`, 110 KB, holding both metrics, upper-triangle matrices and provenance including the sha256 of each source matrix.

**Inputs:**

- `results/rf00167/T64/manifest.json`, giving 76 instances, 57 shared positions, threshold T = 64 and 8 exclusions.
- `matrices/disc_base.csv`, the geometric discrepancy D.
- `matrices/ufix_tm_user.csv`, the US-align fixed-correspondence TM-score normalised by n = 57. It is symmetric with diagonal 1, and this was checked.
- `clustering/clustering_summary.json`, the stability evidence, which is copied into the output.

The standard-mode TM files, which cover 79 chains, are not used.

**Checks before writing:**

- The matrix order equals the manifest order.
- Values are finite and symmetric, and the diagonal is correct.
- Every clustered instance is a mapped instance with an observed match residue at each of the 57 positions in this checkout's mappings.
- Clustered plus excluded instances cover exactly the 84 mapped instances.

**Algorithm.** Average linkage (UPGMA), ported from the demo's `clustering.js` with the same tie-breaking and drawing order. A unit test confirms it gives scipy's heights and partitions.

- The cut undoes the k − 1 highest merges. The stored height is the midpoint.
- Groups are numbered by size and named "Structural group N".
- The representative is the medoid, with ties going to the earliest instance in matrix order.
- Stored per group: within-group mean and maximum, the across-group mean and the ratio.
- Stored per metric: silhouette for k = 2–10, the cophenetic correlation, and whether the cut is the largest gap of the tree.

## Clustering results

| | Geometric discrepancy (default) | TM-score (1 − TM for clustering) |
|---|---|---|
| Cophenetic correlation | 0.9845 (= research) | 0.9562 (= research) |
| Fixed cut | k = 2 at 0.495 (merges 0.397 / 0.594); largest gap of the tree | k = 4 at 0.327 (merges 0.283 / 0.372); largest gap of the tree |
| Groups | 72 + 4 | 51 + 17 + 4 + 4 |
| Representatives | 6UC8 B; 5E54 B | 6UC8 B; 9V4V A; 5E54 B; 9NDD A |
| Status | **stable** | **provisional** |
| Duplicate floor | 0.124 | 0.131 (1 − TM) |

**Discrepancy evidence.** Its silhouette is the best for k = 2–10, at 0.590. The partition holds at every alternative coverage threshold (63, 65, 67) and under every structural-element jackknife, with ARI 1.0, and in 100 position bootstraps, with median and 5th percentile ARI 1.0. The 4-member group is 5E54 A/B and 5SWD A/B.

**TM evidence.** k = 4 holds at all three alternative coverage thresholds (ARI 1.0); k = 2 does not hold at threshold 67. The research ran no bootstrap or jackknife for TM, and it describes TM as a sequence/fold-level descriptor. The TM groups are therefore labelled provisional in the interface and the JSON.

**Data counts.** The 76 clustered instances come from 61 PDB entries; the family has 84 instances from 66 entries. The 8 instances not clustered have fewer than 64 observed core positions: 2XO1 A, 3FO4 A and 9NBC A (63 each), 6P2H A (55), and 9NBH A–D (26 each).

## Frontend

| Part | Source in `antibody_antigen_interface/app-react/src/` | Here |
|---|---|---|
| Group colours, structure palette, overlay cap 8, detail limit 50, small group under 5 | `groups/clustering.js`, `StructuralGroupsApp.jsx` | `src/lib/structureGroups.ts` |
| Collapsed group dendrogram | `groups/GroupOverview.jsx` (+ CSS) | `src/components/structures/GroupOverview.tsx`. Added the floor line, representative label and selected-instance mark. |
| Instance dendrogram | `groups/Dendrogram.jsx` | `structures/Dendrogram.tsx`. Added the floor line and instance mark. |
| Heatmap | `states/DissimilarityHeatmap.jsx` (viridis, grey diagonal) | `structures/PairwiseHeatmap.tsx`. One fixed scale per metric from 0 to the family maximum; drag-to-drill omitted. |
| Groups panel | `groups/GroupPanels.jsx` (`GroupsCard`) | `structures/GroupsCard.tsx`, with metric-specific labels and clustering details including the status and rationale. |
| 3D view | `states/SuperpositionViewer.jsx` (3Dmol.js 2.5.2 from the CDN) | `structures/AlignmentViewer.tsx`, with legend, visibility, add member and compare representatives in `StructuresTab.tsx`. |
| Associated features | `groups/features.js`, `AssociatedFeatures` | `structures/AssociatedFeatures.tsx` |

**Layout:**

- **A.** The existing `InstanceTable` with a compact column set (PDB ID, Chain, Structural group, Title, Method, Res.), beside the Structural groups panel. The new optional props are group badges, the Not clustered badge, a dendrogram-order sort and the "selected group only" filter.
  - "Expand table" shows every existing column at full width and hides the panel. Search, sort, page, metric and group are kept, because the table keeps its place in the component tree.
  - Below the `lg` breakpoint the panels stack.
- **B.** "Pairwise structural comparison" and the "3D alignment view".
  - The metric selector reads Geometric discrepancy | TM-score, with discrepancy as the default.
  - With 76 clustered instances, which is above 50, the view is the collapsed group dendrogram with the cut and floor lines.
  - Drill-down is offered only for a selected group of 2–50 members; its subtree is taken from the stored linkage. For the 72-member group the panel says why the matrix is not drawn.
  - Downloads: the JSON, long-format pairwise CSV for the active metric, group and not-clustered CSV, and an SVG of the current figure.
- **C.** Associated features.

**State.** `state/structures.tsx` sits above the tabs. It holds the metric and a selected group per metric, so group IDs are never carried across metrics and a new metric opens on its largest group. It also holds the table mode, the representative comparison and the displayed and hidden structures. The selected instance remains the app-wide selection.

**3D alignment:**

- Each chain is fetched from the PDBe ModelServer for its `label_asym_id` and model 1.
- It is placed in the family frame with the existing `orientation.json` transform. That is the Kabsch fit on C4′ onto 4FE5 B: 48 good, 26 moderate and 2 poor among the 76 clustered chains.
- 3Dmol.js 2.5.2 could not parse ModelServer mmCIF, even unmodified, so atoms are parsed from `atom_site` and written as PDB records. The first alternate location is kept, and polymer residues are written as ATOM.
- By default the view shows every member for groups of 2–8, otherwise the representative.
- Up to 8 structures could be shown (5 per structural group since the 2026-10-17 change below). Each has a visibility checkbox and can be removed. "Compare representatives" overlays two medoids.
- 5E54 B and 5SWD B are flagged as poor fits. The ⓘ text explains that this overlay is not the base-centre fit behind D.

**Associated features:**

- The unit is the clustered chain instance, compared with the other clustered instances; not-clustered instances are excluded.
- Each chip shows k/n for both sides.
- Ligands are labelled entity-level PDBe annotations, from the `entity_chains` of the default-track bound molecules.
- Nonstandard nucleotides are chain-specific, from `modifications.json`, at any position.
- Species and method are compared as categories; resolution, structural coverage and chain length as medians.
- The display rule is a difference of at least 15 points, with "Show all features". No tests of significance.

**Replaced.** The Structures tab's single-instance Mol\* view and the three placeholder panels are gone. The Family Overview Mol\* view is unchanged, and tests that checked Mol\* now do so there.

## Tests

- **pytest, `tests/test_structure_clustering.py`** (4 tests):
  - UPGMA against scipy (heights, partitions, cophenetic correlation) and determinism;
  - ties, cut midpoint, medoid tie-breaking, group statistics;
  - matrix validation errors;
  - RF00167 regression: counts 84/66/76/61/8, the exact excluded IDs, group sizes and representatives for both metrics, rebuilding the tree from the stored matrix, medoids, contiguous groups in leaf order, cophenetic correlations equal to the research values, cut height and labels.
- **vitest, `src/lib/structureGroups.test.ts`:**
  - tree, matrix and display helpers; feature counting;
  - atom_site parsing and PDB writing — this caught a bug that stripped the prime from C4′;
  - RF00167 clusterings and the 1 − TM distances.
- **Playwright, `tests/structures-clustering.spec.ts`** (7 tests):
  - totals and exclusions; default metric, group and medoid;
  - the collapsed view and the 72-member summary; drill-down of group 2 (16 cells, dendrogram order) and coordinated selection;
  - 3D content, poor-fit flags and visibility;
  - TM switching: its own groups and representative, the largest-group fallback, a 17×17 drill-down showing TM-score values with a 1 − TM legend, and the discrepancy selection restored;
  - compact/expanded persistence, badges and the group filter;
  - compare representatives, the 8-structure cap, a not-clustered selection, feature counts;
  - a missing clustering file, and a narrow screen.
- **Updated:** `structures-selection`, `orientation`, `refinements`, `annotations` and `screenshots`, where Mol\* checks moved to the Family Overview and column expectations now cover both table modes.

Screenshots: `docs/screenshots/m4_structures_default.png`, `m4_structures_expanded.png`, `m4_structures_group2_drilldown.png`, `m4_structures_3d_alignment.png` and `m4_structures_tm.png`.

## Caveats and limitations

- Groups are descriptive. The discrepancy split is the only stable one, and the TM groups are provisional.
- D is dimensionless and depends on n, so values are comparable only within RF00167.
- The 3D superposition is C4′ onto the family reference, not the base-centre fit of D. It is the same for both metrics, and two chains have poor fits.
- 3Dmol.js loads from jsdelivr at run time, as in the demo, so the 3D panel needs internet access, as Mol\* and PDBe already do.
- The research outputs remain uncommitted in the research worktree. The provenance records their paths and sha256 hashes. Committing or archiving them is a separate decision.
- Ligand features are entity-level and may not reflect contacts with every chain. The nonstandard-nucleotide feature includes terminal components outside the mapped region, for example GTP.
- RF00059, RF00162 and RF00177 are not integrated.

## Final test run (2026-10-16)

| Suite | Result |
|---|---|
| pytest | 54 passed |
| vitest | 23 passed |
| Playwright | 51 passed |

Fixed during the run:

- With group badges, the expanded table was wider than 1400 px. The truncated Title column is now 11rem when groups are shown, so every column fits again.
- `refinements.spec.ts` now finds tooltips by their labels instead of their position, because new ⓘ icons shifted the indices.

## Panel layout refinement (2026-10-17, frontend only)

**Pairwise structural comparison panel:**

- One line above the figure, for example "76 RNA instances · 2 structural groups · Fixed cut: 0.495". TM-score adds a "provisional" tag.
- One sentence below the figure:
  - collapsed view: "Groups with more than 50 instances are shown in collapsed form.";
  - group drill-down: "Pairwise D (or 1 − TM-score) within Structural group N, in dendrogram order.";
  - full matrix: "… between all clustered RNA instances".
- The drill-down button sits beside that sentence.
- The metric definition, linkage and cut rationale, the duplicate floor, the 50-instance rule and the interpretation caveat moved into the panel's ⓘ popover, one paragraph each.

**Heights:**

- Both panels share one grid row of equal height. The pairwise figure keeps its natural size and is never stretched.
- The 3D viewer flexes into the remaining height, with a minimum of 240 px on desktop and 320 px when stacked. Its ResizeObserver refits the scene.
- The representative-comparison controls fit on one row; the button reads "Compare".
- Nothing about clustering, data, viewer behaviour, the table or features changed.

**Tests:** `structures-clustering.spec.ts` has a new layout test and updated text expectations. It checks:

- the caption, the note, and the popover content;
- equal outer heights for both metrics, in collapsed and drill-down views, and during a representative comparison;
- the viewer growing and its canvas following the box, fullscreen and centre;
- the legend never clipped, and narrow-screen stacking without overflow.

Results: vitest 23 passed, Playwright 52 passed. Screenshots: `docs/screenshots/m4_structures_panels_aligned.png` and `m4_structures_panels_aligned_drilldown.png`.

## Table and viewer refinements (2026-10-17, frontend only)

**Table:**

- "Region" is replaced by "Length (nt)": the mapped length, with numeric sorting. Hovering a value shows its range, as "RNA region: 3–69".
- Methods are shown as X-ray, EM or NMR in both modes, from `shortMethod` in `lib/format.ts`. The full value is kept in the data and shown on hover.

**3D panel:**

- The viewer minimum is 340 px on desktop, and the canvas is about 588 × 346 px in the default view.
- "Displayed structures" became a one-row footer of chips: visibility, colour, label, ◆ medoid, fit warning, remove, "+ member…" and "+ selected".
- The compare box is tighter.
- To balance the heights, the collapsed overview's minimum total row height is 320 px; this is the demo's own stretch for pages with few groups.
- Panel heights are equal in every view tested. Measured heights:

| View | Height |
|---|---|
| D default | 485 / 485 px |
| D group 2 | 511 / 511 px |
| TM default | 488 / 488 px |
| TM 17×17 drill-down | 951 / 951 px |

**Tests:** a length and method test, and a canvas size of at least 300 × 300 in the layout test. Playwright: 53 passed before this change; the affected specs, 18, pass after it.

## Structural groups panel formatting (2026-10-17, display only)

**Headings by metric:**

| Metric | Headings |
|---|---|
| Geometric discrepancy | "Within-group geometric discrepancy", "Across-group geometric discrepancy" |
| TM-score | "Within-group 1 − TM-score", "Across-group 1 − TM-score" |

**Values:**

- Rows read "Mean" and "Maximum", shown to two decimals.
- The ratio is shown to one decimal. It is computed from the stored, unrounded means, so the shown figures may not divide exactly to it; the ⓘ text says so.
- The stored values and the JSON are unchanged.

**New Playwright test.** For every group under both metrics it recomputes the within-group and across-group mean, the maximum and the ratio from the stored matrix, using 1 − TM for TM-score. It checks them against the JSON to 1e-5 and against the displayed two-decimal text.

The Structures spec passes (10 tests). Screenshots: `docs/screenshots/m4_groups_panel_discrepancy.png` and `m4_groups_panel_tm.png`.

## Tooltip simplification (2026-10-17, text only)

Tooltips are now one or two sentences in the selected metric's terms, with no bare "D":

- **Coverage:** counts, threshold and shared positions taken from the data.
- **Representative:** the medoid definition.
- **Within-group and across-group:** "Mean pairwise geometric discrepancy / 1 − TM-score …", with no pair counts and no "Not a statistical test".
- **Structural groups:** meaning plus the interpretation caveat.
- **Pairwise comparison:** what the metric measures and what the dotted floor line means.
- **3D alignment:** visual superposition only; poor fits may look offset.

The axis, legend and cell tooltips also say "geometric discrepancy" (or "discrepancy" in the narrow dendrogram axis).

**Moved into "Clustering details":** the metric definition, linkage and cut rationale, the floor definition, pair counts and the ratio note, the 50-instance rule, the list of not-clustered instances with reasons, and the interpretation caveat.

**Tests:** a new Playwright test checks exact tooltip texts for both metrics, the absence of "D", pair counts and "statistical test", the content of Clustering details, and that every popover stays inside the viewport.

| Suite | Result |
|---|---|
| vitest | 24 passed |
| Playwright | 55 passed |

**3D alignment tooltip (2026-10-17).** It now reads: "RNA structures are superposed using C4′ atoms at shared Rfam positions. This visual alignment differs from that used to calculate geometric discrepancy. Poor-quality alignments are flagged."

- The full description moved to Clustering details: Kabsch onto the family reference frame of the Family Overview, the same for both metrics, and the poor-fit threshold of C4′ RMSD above 4 Å.
- The tooltip test checks the exact text, a width of at most 340 px, and that the popover is placed over its own panel.
- Screenshot: `docs/screenshots/m4_tooltip_3d_alignment.png`.

## Representative comparison path and group colours in 3D (2026-10-17)

**Comparison path.** In comparison mode the stored-tree route between the two groups is drawn in heavy ink, with other branches dimmed. The route uses `pathBetween`, ported from the demo, and runs up to the lowest common ancestor and back down.

- Both compared groups get a marker in their group colour, and the compared group has a dashed outline. The selected group and the selected-instance mark are kept.
- The instance-level dendrogram marks the compared representatives' leaves and draws the route only when both are in the drawn tree. A group drill-down contains only its own representative.
- Leaving comparison removes the route and restores the selected group's default display. The selected group itself never changes.

**Colours.** A group's representative always takes the group colour from the shared palette, matching the dendrogram, the groups panel and the table badges. Other members take alternating lighter and darker shades of their own group's colour (`groupShade` / `structureColours`). Each swatch in the legend carries its colour.

**Tests:**

- The route equals the path computed from the stored linkage, for discrepancy 1–2 and TM-score 1–4, 3–2 and 2–4.
- Colours match the group chips for a representative alone, a 4-member group, comparisons, and both metrics; the table badge colour matches too.
- Exiting comparison restores the previous state; visibility, removal and centre still work.
- The drill-down marks only its own group's representative.
- Vitest covers `pathBetween` and the colour rule.

| Suite | Result |
|---|---|
| vitest | 26 passed |
| Playwright | 58 passed |

Screenshots in `docs/screenshots/`: `m4_compare_g1_rep_alone.png`, `m4_compare_g2_rep_alone.png`, `m4_compare_g1_vs_g2.png`, `m4_compare_tm_g2_vs_g3.png`.

## Associated features refinement (2026-10-17, frontend only)

**Counting audit, before the change:**

- Already correct: each instance was counted once per feature; group membership came from the active metric's assignments; not-clustered instances were excluded; medians skipped missing values.
- **Corrected — mixed regions.** Nonstandard nucleotides inside and outside the Rfam-mapped region were pooled. For example, the terminal GTP was counted as 3/17 for TM group 2. They are now two separate tables, mapped region first.
- **Corrected — ligand availability.** Ligand "absence" was not tied to source availability. The denominator is now the instances whose PDBe entry-entity was evaluated, from `ligands.json` `entry_entity_sources`. All 66 were evaluated, so the counts are unchanged; a 0 is described as "not annotated by PDBe".
- **Corrected — unit mixing.** The header line mixed PDB-entry counts with instance counts; it is gone.
- **Corrected — medians.** They did not report their contributing instances; each median now carries n in its tooltip and in `data-n`.

**Presentation:**

- Counts-first tables: Ligand associations, Nonstandard nucleotides (Rfam-mapped region), Nonstandard nucleotides (outside the mapped region), Species, Experimental method (X-ray / EM / NMR), and Continuous properties as medians.
- Each cell is k/n; the percentage appears only on hover.
- Removed: arrows, green and orange colouring, the 15-point filter, and "Show all features".
- Long lists show 6 rows with a "Show all N" control, ordered by total count.

**Columns:**

- Discrepancy: Structural group 1 (n = 72) and Structural group 2 (n = 4), always in group order. The selected column gets a light tint and nothing else changes.
- TM-score: four provisional groups, shown as "Structural group N (selected)" against "Other clustered RNA instances", with a "provisional groups" tag.

One descriptive note sits at the top. The tooltips use the requested wording, including entity-level ligand provenance and that instances from the same entry are not independent.

**Tests:**

- vitest: `categoricalRows` (deduplication; unknown status excluded from n) and `medianOf` (missing values skipped, n reported).
- New Playwright test: recomputes the ligand, mapped and unmapped nonstandard, species and median values from `instances.json`, `ligands.json` and `modifications.json` for both groups under discrepancy and the first two groups under TM-score. It also checks there are no arrows, enrichment colours or 15-point controls, the column order and selection, the TM layout, and the tooltips.

| Suite | Result |
|---|---|
| vitest | 27 passed |
| Playwright | 60 passed |

Screenshots: `docs/screenshots/m4_features_D_group1.png`, `m4_features_D_group2.png`, `m4_features_TM_group1.png`, `m4_features_TM_group3.png`.

## Associated features: collapsed by default, five rows per category (2026-10-17)

**Collapsing:**

- The panel is collapsed by default. Collapsed, it shows only the heading, the selected-group badge and a chevron.
- The whole header is one toggle button with `aria-expanded` and `aria-controls`, following the Annotations categories, and works from the keyboard.
- The open state lives in the shared Structures state, so it is kept across group, metric and tab switches. The badge updates without opening the panel.

**Expanded view:**

- The ⓘ icon and the "provisional groups" tag sit beside the descriptive note.
- Each category shows at most 5 rows, with its own "Show all N" / "Show fewer" control (`aria-expanded`). Categories of 5 or fewer show every row with no control.
- Order: total supporting instances across all clustered instances, descending, then the identifier. Selection doesn't change it; the metric changes it only through that metric's dataset.

**Tests:**

- New Playwright test: default collapsed state with exactly heading + badge + chevron and height under 60 px; badge updates without opening; keyboard toggle; state kept across groups, metrics and tabs; 5-row defaults (ligands 28, species 6); no control for 4 and 2 rows; independent expansion; ordering by totals; unchanged order and counts on group change; TM ordering; collapsing again.
- The recomputation test now opens the panel and every category before comparing counts.
- The Length sort test now waits for the re-render; it had one timing failure.

Results: vitest 27, Playwright 61 passed. Screenshots: `m4_features_collapsed.png`, `m4_features_D_group1.png`, `m4_features_D_group2.png`, `m4_features_TM_group1.png`, `m4_features_TM_group3.png`.

## Associated features follow the viewing mode (2026-10-17)

**Modes:**

- **Single-group mode** (default): one column for the selected group, with n given as RNA instances. Rows are the features observed in that group.
- **Comparison mode:** derived from the existing representative-comparison state; Associated features keeps no state of its own. It shows exactly the two compared groups, selected first as in the comparison controls, and the union of their observed features.
- The header shows one badge, or "Group A vs Group B". Selecting another group or switching metric ends comparison, as before, and returns to single-group mode.
- The panel's collapse state and per-category expansion are kept.

**Section order (user request):** Ligand associations | Species, then Nonstandard (mapped) | Nonstandard (outside the mapped region), then Experimental method | Continuous properties. The two nonstandard tables now share a row.

Counting rules are unchanged.

**Tests:**

- The recomputation test checks single-group mode for both groups and both metrics: exactly one column, and every row with support. It checks comparison mode for D 1–2 and TM 2–4 and 3–1: exactly the two columns in order, recomputed counts and medians, and the return to single mode.
- The row-limit test checks that group 2 shows only its own ligand row.
- The Length sort test now waits on `aria-sort`, after two timing failures.

| Suite | Result |
|---|---|
| vitest | 27 passed |
| Playwright | 61 passed |

Screenshots: `m4_features_D_group1.png`, `m4_features_D_group2.png`, `m4_features_D_compare_g1_g2.png`, `m4_features_TM_compare_g2_g4.png`, `m4_features_collapsed.png`.

## Family Overview residue focus (2026-10-17)

`frontend/src/components/PdbeMolstarViewer.tsx` now uses Mol*'s built-in structure focus, the same mechanism as the PDBe entry pages. pdbe-molstar already registers that behaviour (`selectInteraction`), so nothing new was added.

**On selecting a resolved residue:**

- `managers.structure.focus.setFromLoci` shows the nucleotide plus its 5 Å surroundings in ball-and-stick, with non-covalent contacts.
- The residue keeps the app's selection colour through `visual.select`.
- The camera closes in on it with `focusLoci` (extraRadius 7 Å, minRadius 12 Å, 250 ms transition), keeping the view direction.

**Representation and state:**

- The Rfam-region representation is now cartoon only, with bases drawn by the cartoon, instead of cartoon plus all-atom sticks. The selection therefore stands out; complete-structure mode is unchanged.
- Each selection replaces the focus, so representations never accumulate. Clearing the selection, or picking a position without coordinates, a deletion or an uncovered position, clears the focus and does not move the camera.
- Switching instance clears the focus before loading and re-focuses that instance's own residue for the same consensus position.
- Reset orientation / Fit restores the whole-structure view and keeps the selection.

**Tests:** `tests/residue-focus.spec.ts` covers:

- the focus target equals the mapped residue, with tracks, 2D and Mol* in agreement (76 → C56);
- the camera centring and moving much closer than the whole view;
- neighbouring, far and repeated selections (a single target);
- reset, and clearing;
- instance switching (4FE5 to 1Y26);
- a deleted position (1Y26 64), no coordinates (2XNW 19) and uncovered (2B57 19), none with a focus or camera move.

The Length sort test now retries its page-back step.

| Suite | Result |
|---|---|
| Playwright | 63 passed, plus the 15-test Structures spec passing twice under parallel load |
| vitest | 27 passed |

Screenshots: `docs/screenshots/m5_overview_default.png`, `m5_overview_residue_focus.png`.

## Family Overview track reorganisation, Phase 1 (2026-10-17)

`frontend/src/components/ConsensusTracks.tsx`:

- **Primary tracks:** three are always visible, in order: Sequence coverage, Structural coverage, then the selected RNA chain's state per position. Scales, denominators, tooltips and state semantics are unchanged. The insertion triangles already on the chain track are untouched.
- **Reference and alignment details:** this section sits directly below them and is collapsed by default, showing only the heading and a chevron. It is a button with `aria-expanded` / `aria-controls`.
- **Expanded content:**
  - a one-line consensus-model summary: Rfam accession and id, 102 match positions, clan members, and an ⓘ with the RF / WUSS / mean PP explanations;
  - the Alignment confidence (mean PP) track;
  - RF and SS_cons as one character per consensus column, taken unchanged from `family.json`.
- **Axis:** one shared Rfam position axis sits below the tracks and the optional section.
- **Alignment:** each part is its own SVG with the same viewBox width and scaling, so all columns align in both states. Every strip shows the selected column and accepts clicks.
- **Removed:** the standalone "Consensus model" card in `FamilyOverview.tsx`. Its content is now only in the reference section.
- **Unchanged:** selection, 2D and Mol* synchronisation, and instance switching.

**Tests:**

- New `tests/overview-tracks.spec.ts`:
  - exactly three primary tracks by default, the reference section collapsed and showing only its heading, the order tracks → heading → axis, tooltips kept;
  - expanded content: the confidence track, RF and SS_cons equal to the family strings, the consensus model shown once;
  - column centres equal within 0.5 px across all strips for positions 1/50/76/102, the axis tick aligned in both states, the selection backdrop in every strip;
  - a click in the reference strip selects; collapsing keeps the selection;
  - an instance change updates the chain track;
  - keyboard toggle;
  - a 600 px viewport with no overflow and aligned strips.
- Updated `overview-layout.spec.ts` (card removed) and `refinements.spec.ts` (confidence tooltips are now in the reference section).

| Suite | Result |
|---|---|
| vitest | 27 passed |
| Playwright | 68 passed |

Screenshots: `docs/screenshots/m6_tracks_default.png`, `m6_tracks_reference_expanded.png`, `m6_overview_selected_position.png`, `m6_tracks_narrow.png`.

**Known:** at narrow widths the SVG tracks scale down as before, while the HTML section heading keeps its size.

## Family Overview Phase 2: per-instance RNA loops (2026-10-17)

**Phase 1 refinements:**

- The visible consensus-model line was removed; its content is now in the ⓘ next to "Reference and alignment details".
- Rows were renamed "Consensus sequence (RF)" and "Secondary structure (SS_cons)".

**Data.** The audit is in `experiments/overview_loops_phase2/` and found no issues; membership is taken from `loops.json`.

**Loop tracks:**

- "RNA loop annotations" sits below the selected-chain track and above the reference section. It is collapsed by default and loads `loops.json` on first opening; its open state is kept across instance changes.
- Three tracks: Hairpin loops (HL), Internal loops (IL), Junction loops (JL = J3 + J4).
- Blocks:
  - One block per original strand of each loop, spanning that strand's consensus match positions in the selected chain (`lib/loops.ts`).
  - Adjacent loops stay separate blocks; a discontinuous loop shows one block per segment.
  - Overlapping loops of one track go to separate lanes, so J3 and J4 sit in two JL lanes.
  - Category colours, rounded ends, small insets.
- Hover outlines every block of that loop and fades the others. The tooltip reads "ID · type · n nucleotides · k segments".
- Loops with no consensus position are listed under "Not on the consensus grid" and are selectable.
- Instance-specific and binary: no family frequencies. The 102-column grid is unchanged and no insertion indicators were added.

**Whole-loop selection:**

- A click selects the loop's original membership as a group (selection kind `group`): all strands, plus members with insert, flank or no consensus state.
- The viewers receive exactly the members in this chain that have coordinates. Members without coordinates stay in the membership; members in another chain are reported as "k in another chain (not shown)".
- 2D: all members are coloured with one call to the viewer's non-zooming selection path, so the user's zoom and pan are preserved. Single nucleotides now also keep the zoom.
- Mol*: the structure focus of all resolved members, and one camera `focusLoci` over their combined extent.
- Status line, for example "Internal loop · IL_4FE5_001 · 2 segments · 6 nucleotides", with the unit ids on hover. The residue list in the selection panel is folded for loops.
- A new selection replaces the old one; Clear removes both kinds.
- An instance change clears a loop selection, because a loop belongs to its instance. A consensus position is still kept.
- Collapsing the section keeps the selection.
- Group positions now count only members of the current instance, so another chain's member is never translated onto this chain.
- Ordinary column, 2D and 3D clicks still select single nucleotides.

**Tests:**

- vitest `src/lib/loops.test.ts`: instance filtering, blocks per strand, adjacent IDs kept distinct, discontinuous loops, lanes for overlaps, inserts, flanks, other-chain and missing-coordinate members, whole-loop membership excluding intervening residues.
- Playwright `tests/overview-loops.spec.ts` (full run: vitest 32, Playwright 72 passed), all values recomputed from `loops.json`:
  - every 4FE5 B loop's blocks, J3/J4 lanes, category colours, grid edges within 1 px, hover;
  - IL, HL, J3 and J4 selection with exact membership in the viewer ids, 2D, Mol* selection, focus and camera; intervening stem excluded; status and title; 2D zoom unchanged after a wheel zoom; folded member list;
  - replacement, collapse, clear, single nucleotides;
  - instance switch clearing the loop; inserts (IL_9IWF_001, residues 59/60); off-grid and inter-chain loops (9NBH).
- `overview-tracks.spec.ts` updated for the Phase 1 refinements.

Screenshots (`docs/screenshots/m7_*`): sections collapsed, loop tracks expanded, selected hairpin / internal / junction loop, J3/J4 overlap lanes, 3D focus on the complete loop, 2D loop highlight with zoom kept.

**Limitations:**

- RF00167 has no loop nucleotides without coordinates, so that path is covered by unit tests only.
- Inter-chain members cannot be shown in the single-chain 2D view and 3D selection.
- Mol*'s focus view also shows the 5 Å surroundings, so a large junction selection looks busy (resolved in the refinements below).

## Family Overview loop selection and viewer refinements (2026-10-17, frontend only)

No loop data, mappings or other tabs changed. The Structures tab's 3Dmol viewer is untouched.

**Simplified Mol\* display (Family Overview):**

- The focus representation shows only the selected nucleotides as heavy-atom ball-and-stick, in a uniform dark selection colour.
- The 5 Å surroundings and the non-covalent contact lines are switched off through the focus behaviour's own parameters.
- Water is hidden by visibility in both view modes. The component and coordinates stay in the structure.
- The RNA cartoon remains the structural context.
- This replaces the earlier "residue plus 5 Å surroundings and contacts" focus for single nucleotides too.

**Structured loop tooltip:** a hover card replaces the one-line SVG title. It shows the category with its code, the original BGSU loop ID, the nucleotide and segment counts, and members in another chain where present. A "Motif" line appears only when `loops.json` carries an explicit BGSU motif-group assignment. That line gives the group ID with the BGSU annotation text, plus the Motif Atlas release and the assignment basis. 456 of 496 RF00167 loops have such an assignment; 40 have none and show no motif line. The card is clamped inside the viewport.

**Toggle-off:**

- Clicking any block of the active loop again clears it.
- Clicking the active consensus column, residue or 2D/3D nucleotide again clears it as well.
- Each such deselection, and Clear selection, increments a reset counter in the selection state. Both viewers then restore their default views: the fitted 2D topology with zoom and pan reset, and the Reset orientation / Fit camera in Mol\*.
- Selecting a different loop or nucleotide replaces the selection without an intermediate reset.
- An instance change clears an incompatible loop. The newly loaded structure gets its default view as before.

**Selected-loop outline:** a solid red outline is drawn inside every block of the selected loop ID, so block widths, positions and gaps never change. Hover uses a thin dashed dark outline. Category fills are unchanged.

**Defect found and fixed:** in Phase 2, changing a selection reset the user's 2D zoom. The 2D viewer's own clear-selection call is a "click" deselect, which also resets zoom. The Phase 2 test did not catch it, because the viewer disables wheel zoom, so its "zoomed" state was never zoomed. Selection changes now clear colours with a non-click deselect, and the tests zoom through the viewer's zoom-in button.

**Tests:**

- New `frontend/tests/overview-loop-refinements.spec.ts` with five tests covering checklist items 1 to 15.
- New `frontend/src/state/selection.test.ts` for the toggle and reset reducer, and new tooltip tests in `lib/loops.test.ts`.
- Updated existing tests:
  - overview-loops: outline, tooltip and zoom button;
  - residue-focus: Clear restores the default camera;
  - overview-layout: radius tolerance, since Mol\* re-derives it from scene bounds;
  - two timing checks now poll, in single-selection and structures-clustering.
- Full run: vitest 41, Playwright 78 passed.

**Screenshots** (`docs/screenshots/`):
- `m8_3d_simplified_internal_loop.png`
- `m8_3d_simplified_junction_loop.png`
- `m8_loop_tooltip.png`
- `m8_red_outline_internal_loop.png`
- `m8_overlap_selected_vs_hover.png`
- `m8_selected_hairpin.png`
- `m8_selected_internal_loop.png`
- `m8_selected_junction_loop.png`
- `m8_deselected_default_views.png`

`m7_2d_loop_zoom_kept.png` was regenerated with a real zoom.

## Family Overview consistent RNA colour (2026-10-17, display only)

The Family Overview Mol\* viewer now draws every instance's RNA in one fixed slate grey-blue, `#64748B`. The colour is passed to the viewer as a prop from the structure explorer. Nothing derives it from the PDB ID, chain, instance or structural group.

- **Rfam region mode:** the region cartoon is built with a uniform slate colour theme.
- **Complete structure mode:** the polymer component is recoloured to slate. Ligands and ions keep their element colours.
- **Selection:** nucleotides and whole loops now use violet `#7C3AED` instead of `#6250D6`, for better contrast with slate. The Mol\* overpaint and the focus target use the same violet. Both are identical for every instance.
- **Unchanged:** the red loop-block outline, the HL, IL and JL colours, the track selection shading, the 2D viewer colours, the Structures tab's 3Dmol viewer, and all coordinates, camera, focus and representation types.

**Tests:** `frontend/tests/overview-structure-colour.spec.ts`. It checks the representation colour parameters and also classifies the rendered pixels by hue.
- **Default colour:** 4FE5 B, 1Y27 X and 2EET A are slate in both view modes, after instance switches and after Reset orientation.
- **Selection:** a selected nucleotide and the J3_1Y27_001 loop render violet over slate.
- **Clear:** clearing leaves no violet.
- **Other test updates:** four older checks in `annotations.spec.ts` now match the selection hook partially, since it also records the colour. One refinement test now waits for the complete-structure state.
- **Full run:** vitest 41, Playwright 81 passed.

**Screenshots:**
- `m9_colour_1y27_X_default.png`
- `m9_colour_2eet_A_default.png`
- `m9_colour_2eet_A_complete.png`
- `m9_colour_2eet_A_selected_junction.png`
- `m9_colour_2eet_A_selected_junction_overview.png`

The m8 refinement screenshots were regenerated with the new colours.

## Family Overview Phase 3: insertion markers (2026-10-17, frontend only)

The audit is in `experiments/overview_insertions_phase3/README.md`, with `audit_insertions.py` and `audit_insertions.json`. The existing mappings already carry every needed field, so the mapping pipeline and the data are unchanged.

RF00167 has 17 genuine insertion runs, all the same event: `CA` between consensus 77 and 78, in the 9IWF, 9IWG, 9LJN, 9LKC, 9LKE, 9LKF, 9LMC and 9V4U–9V51 chains. It has no terminal insert states and no local-end residues. The 635 terminal "not_submitted_flank" residues are not insertions.

**Display:**
- An amber `+N` label sits on the column boundary in the selected-chain track, in the gap above the row, with a tick through the row. There is no extra column.
- The track viewBox, coverage-bar geometry and axis are identical between a chain with a marker (9IWF A) and one without (4FE5 B). Loop-block geometry is unchanged by marker interactions.
- The hit target is larger than the label. Markers are keyboard-focusable buttons with an aria-label.
- The old insertion triangles were replaced.
- An "i" in the track names the coordinate system with the requested sentence. Legend & help explains all five per-position states. The existing state colours and patterns are unchanged.

**Interaction:**
- Hovering shows "Insertion · +2 nt", "Between Rfam positions 77 and 78" and "Sequence: CA". Terminal runs get one-sided wording.
- Clicking opens details: the sequence, a table of label_seq_id, author number with insertion code, and coordinates, plus provenance. A click never selects a consensus position.
- The details offer "Show in 2D and 3D". It highlights the inserted residues as an instance-bound group with no consensus position, so the viewers address them by label_seq_id.
- Clicking the open marker again closes the details and clears that selection, restoring the default views. Escape, the close button and an outside click also close the details.
- Switching instance clears the details, any marker tooltip, and the insertion selection.
- Loop membership is unchanged: IL_9IWF_001 still includes residues 59–60.

**Tests:**
- `src/lib/insertions.test.ts`: 8 tests. They cover real mappings (9IWF A, 9LJN A, flank-only chains, and all 84 chains against their summaries) and synthetic terminal, local-end, split and aggregated cases.
- `tests/overview-insertions.spec.ts`: 5 tests covering markers, position, grid invariance, tooltip, details, view and toggle, column click beside the marker, instance switching, loops, responsive widths, accessibility and legend.

**Screenshots:**
- `m10_insertion_marker_9iwf_A.png`
- `m10_marker_closeup.png`
- `m10_insertion_tooltip.png`
- `m10_insertion_details.png`
- `m10_no_insertions_4fe5_B.png`
- `m10_alignment_with_loops_9iwf_A.png`

**Limitation:** RF00167 exercises only one insertion event. Terminal, local-end, multi-run and collision handling are covered by unit tests only.

## Family Overview tooltip wording and presentation (2026-10-17, frontend only)

No data, calculations, mappings or selection behaviour changed. Values come from the same `consensus_stats.json` and mapping fields as before.

**Per-position tooltips:** the native SVG titles became one shared hover card with a bold heading, a main value and an optional muted line. It is at most 340 px wide, rendered in a portal so no container clips it, and kept inside the viewport.

| Track | Heading | Main value | Muted line |
|---|---|---|---|
| Sequence coverage | `Sequence coverage · Position 41` | `29/29 unique sequences (100%)` | Deleted and not-covered counts, only when nonzero |
| Structural coverage | `Structural coverage · Position 40` | `80/84 RNA instances (95.2%)` | `84 instances contain the nucleotide in sequence`, only when it differs from the observed count |
| Selected chain | `Selected chain · Position N` | `Coordinates available`, `Coordinates unavailable`, `Deleted in chain` or `Outside aligned region` | A short explanation for unavailable and deleted positions |
| Alignment confidence | `Alignment confidence · Position N` | `Mean PP 0.987 (high)` | Number of aligned unique sequences |
| RF / SS_cons | `<row> · Position N` | The character | None |

Percentages keep one decimal, without a trailing `.0`.

**Information tooltips:**
- Each now has a heading and one or two sentences. This uses a new optional heading on InfoTip, which other tabs don't pass, so their popovers are unchanged.
- The selected-chain text is "Shows nucleotide availability at Rfam match positions. +N markers indicate insertions between positions."
- The loop and reference texts are as specified.
- The removed detail moved to Legend & help. That covers loop segments and lanes, the red selection outline, click again to clear, and members without a consensus position. It also covers RF, SS_cons/WUSS, mean PP, and the model with its clan members.
- The loop-block and insertion tooltips are unchanged.

**Fix:** dashed "no coordinates" cells were hoverable only on their stroke, because they have no fill. They now respond over the whole cell.

**Tests:** new `tests/overview-tooltips.spec.ts` with 3 tests; values come from the JSON. `overview-tracks`, `overview-insertions` and `refinements` were updated for the new wording. Full run: vitest 49, Playwright 91 passed (the refinements spec passed on its rerun after the update).

**Screenshots:**
- `m11_tooltip_sequence_coverage.png`
- `m11_tooltip_structural_coverage.png`
- `m11_tooltip_selected_chain_unavailable.png`
- `m11_tooltip_selected_chain_info.png`
- `m11_tooltip_loop_annotations_info.png`
- `m11_tooltip_reference_info.png`

## Family Overview simplified controls and Viewer details (2026-10-17, frontend only)

**Toolbar:** the toolbar now holds only the instance selector and Clear selection, side by side on desktop and wrapping on narrow screens.
- The "Consensus position" input and "Select" button are removed. The selection state and its `toggle` and `select` API are unchanged.
- Positions are selected from the tracks, the 2D viewer, the 3D viewer, the loop blocks and the insertion markers.
- The idle sentence ("No nucleotide selected: …") is gone. The status line appears only for an active selection and keeps the existing summaries: nucleotide, loop with ID, type, segments and nucleotide count, and insertion.
- Clear selection is disabled when nothing is selected.

**Viewer details:** the paragraph under the viewers became a "Viewer details" section, collapsed by default, using the same style as Legend & help.
- It holds the specified 2D and 3D sentences.
- A muted technical line covers the viewer names, the author chain and `label_asym_id`, addressing by `label_seq_id`, the Rfam-region default, coordinates only, no water or contact lines, click-again and Clear restoring views, and per-instance nucleotides without a consensus position.

**Tests:**
- New `tests/overview-controls.spec.ts` with 4 tests: removed controls, no idle text, Clear enabled and disabled, selection from tracks, 2D, 3D, loops and insertion markers, Viewer details, and layout at 1400, 700 and 420 px.
- The 13 specs that used the input now click the consensus column through a shared helper in `tests/helpers.ts`, which is the real user path.
- Three idle-text assertions now expect no status line.
- Full run: vitest 49, Playwright 96 passed; type-check clean.

**Screenshots:**
- `m12_default_no_selection.png`
- `m12_active_loop_selection.png`
- `m12_viewer_details_expanded.png`

## Family Overview: Mol* footer removed, Viewer details simplified (2026-10-17, frontend only)

**Footer:** the alignment overlay at the bottom of the Mol\* panel is removed ("Reference frame: superposed on the family reference, RMSD … · reference structure"). It was an overlay, so no layout space is left behind. The camera fit and its overlay margins are unchanged, as are the transforms, the camera and selection. Mol\*'s orientation axes stay on, and a test checks this.

**Viewer details** stays collapsed by default and now has two compact sections:
- **Molecular viewers:** the specified text, with no chain or residue identifiers, water or contact-line wording.
- **3D alignment:** read from `orientation.json` for the selected instance. It shows the reference structure, marked "(this structure)" for 4FE5 B, and whether the structure is superposed or shown as deposited. It also shows the RMSD, the matched atom pairs with their atom types, and the alignment quality.

Chain and residue-numbering conventions remain in `frontend/ARCHITECTURE.md`.

**Tests:**
- `overview-controls`: no footer, axes on, collapsed and expanded states, exact viewer text, no internal identifiers, alignment for 4FE5 B and 1Y27 X from the JSON, and selection with reset.
- `orientation` and `structures-selection`: they read the alignment from Viewer details instead of the footer, covering 1Y26 (1.90 Å, good), 4FE5 (this structure) and 9NBH (as deposited, insufficient).
- Full run: vitest 49, Playwright 97 passed; type-check clean.

**Screenshots:**
- `m13_viewers_without_footer.png`
- `m13_viewer_details_4fe5_B.png`
- `m13_viewer_details_1y27_X.png`

## Structures tab: the selected instance follows into the 3D alignment view (2026-10-17, frontend only)

Before this change, selecting an instance in the Structures tab updated only the shared selection and the highlights; the 3D view kept showing the selected group. Now a table row, heatmap leaf or dendrogram leaf also brings the instance into 3D:
- **Already displayed:** the structures are unchanged and the emphasis moves to it.
- **In another group:** that group is selected, which also leaves a comparison. A small group shows all its members; a large group shows its representative plus the instance.
- **Same group, not yet displayed:** it is added to the overlay, up to 5 structures of one group; when that group is full, it is shown with the representative.
- **Not clustered:** it is shown on its own, with the caption "… is not clustered; shown on its own" and the existing reason note.

**Emphasis:** the instance is drawn with a thicker cartoon, the other structures fade to 35 % opacity, and its "Displayed" chip is outlined in blue. This restyles the models in place without reloading. Group colours are unchanged, so in closely superposed same-group overlays the emphasis is subtle. The Family Overview instance dropdown still doesn't change the Structures 3D view.

**Tests:** `tests/structures-follow.spec.ts`, 2 tests on the discrepancy groups from the JSON (group 1 with 72 members, group 2 with 4, and an excluded instance). The existing Structures specs still pass (22 tests).

**Screenshot:** `m14_structures_selected_instance_in_3d.png`

**Overlay limit (2026-10-17):** the limit is now 5 structures per structural group, representative included, down from 8 overall. The 8 came from the PDBe-KB demo and had no scientific basis. Same-group structures are told apart by shades of the group colour, which stay distinguishable for about five. Groups with up to 5 members are shown in full by default; the RF00167 defaults are unchanged. When a group is full, the add control reads "Max 5 per group", and a newly selected instance is shown with the representative. Tests: a unit test for the per-group count and the default, an extended `structures-follow` test (fill to 5, then a sixth), and the updated cap test in `structures-clustering`.

## Structures table: default order by the active structural group (2026-10-17, frontend only)

**Default order:** when clustering data is present, the table opens sorted by the "Structural group" column.
- The active group comes first, with its representative (medoid) on the top row of page 1.
- The rest of that group follows in dendrogram order, then the other groups in dendrogram order, then non-clustered instances.
- Selecting another group, from a chip, badge or the dendrogram, moves its representative to the top and returns to page 1. The metric switch uses that metric's groups.
- It is an ordinary sort: the column header shows it as ascending, and other headers re-sort as before. Sorting the group column ascending again restores the default and page 1.
- Without clustering data, the table keeps the PDB ID order.

**Selection:**
- The representative is not auto-selected; the shared selected instance, 4FE5 B by default, stays highlighted wherever it falls.
- Opening the tab always shows page 1.
- A row selected in the same step as a group change keeps its page visible, so a clicked row never disappears. Clearing a search keeps the selected row's page.

**Fix:** the group column was sorted as text, so 10 came before 9. It is now numeric, which also corrects the old dendrogram-order sort.

**Tests:**
- New `tests/structures-table-order.spec.ts` with 3 tests: entering the tab, group 2, TM group 2, user sort and return to the default, a row selected deep in a group, and no clustering data.
- New unit test `src/components/InstanceTable.test.ts`.
- `refinements` and `structures-selection` are updated where they assumed the old first-visit rule of opening on the selected row's page in PDB ID order.

**Screenshots:**
- `m15_table_default_group_order.png`
- `m15_table_group2_order.png`

## Family Overview 3D controls: Clear selection replaces Focus selection (2026-10-17, frontend only)

The Mol\* control bar now reads "Rfam region | Complete structure · Reset orientation / Fit · Clear selection". The new button calls the same shared clear as the toolbar's Clear selection: it removes the nucleotide, loop or insertion selection from the tracks, 2D and 3D, and both viewers return to their default views. The 2D view is refitted with zoom reset, and Mol\* goes back to the Reset orientation / Fit camera. It is disabled when nothing is selected. "Focus selection" was removed, since a selection is already focused when it is made.

**Tests:** a new test in `overview-controls` covers the disabled state, a J3 loop with 2D zoom, both default views after the click, and a single nucleotide. `refinements` was updated, since it no longer clicks Focus selection.

## Family Overview: collapsible, zoomable "Sequence and annotations" (2026-10-17, frontend only; prepares RF00177)

**Shared scale** (`ConsensusTracks.tsx`):
- Every strip is an SVG in track units drawn at a fixed 1.2 px per unit, so text and row heights never stretch.
- Zoom changes only `cell`, the units per consensus position. Every strip uses the same `x(pos) = LEFT + (pos − 1) · cell`, and the number of positions comes from the loaded family.
- Insertion markers, loop blocks, the selection shading and click-to-position all use that one scale. Clicks convert pointer coordinates through the SVG transform, so scrolling and zoom never shift the selected position.

**Layout:**
- One horizontal scroll container holds all strips and the axis.
- Each strip's labels sit in a sticky left column, and the loop and reference headings are sticky too, so names stay in view while the tracks scroll together.
- The page itself never scrolls sideways.

**Controls:** a "Sequence and annotations" heading with a collapse chevron, and −, + and Fit all buttons.

| Behaviour | Detail |
|---|---|
| Default | Expanded and at Fit all. RF00167 looks as before, at about 10 px per position. |
| Zoom levels | 1, 1.5, 2, 3, 4, 6, 8, 11, 16, 24, 32 px per position, never below Fit all |
| Zoom anchor | Zooming keeps the selection, or the view centre, in place |
| Collapse | Keeps the selection and never reloads the 2D or 3D viewer |
| Selection elsewhere | A selection from 2D, 3D or a loop that lies off-screen is scrolled into view, without re-zooming |

**Level of detail:**

| Element | Rule |
|---|---|
| RF and SS_cons letters | Shown from 7 px per position; below that, a hover target per position and a "zoom in to read" label |
| +N insertion markers | Shown from 4 px per position |
| Gaps between columns | Shrink to zero at overview scale |
| Axis labels | Step from 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, at least 60 px apart |
| Loop blocks | Remain one block per original segment |

**Long sequences:** checked with a synthetic 1,533-position fixture that exists only in `tests/long-sequence-fixture.spec.ts`. It serves RF00167 files with every per-position array tiled, under the test accession RF99999; it is not RF00177 data.
- First render took 0.75 s, and each zoom step 0.1–0.3 s.
- At Fit all, every position is present with no letters and labels every 100.
- At letter scale, the strips stay aligned in mid-model, and a click selects the column under the pointer.
- No new dependency and no virtualisation were needed.

**Tests:**
- New `tests/tracks-zoom.spec.ts` with 4 tests covering default fit and numbering, zoom steps, and alignment at minimum, default and maximum zoom. It also covers collapse with the selection and the same Mol\* instance, zoomed clicks, scrolling a 3D selection into view, zoomed discontinuous loops, insertion anchoring at maximum zoom, and a narrow window.
- New `tests/long-sequence-fixture.spec.ts`.
- `tests/helpers.ts` now reads each strip's geometry (`data-cell`, `data-left`). The specs that hard-coded 9 units per position use it, and geometry comparisons are expressed in columns.
- Full run: vitest 55, backend 60, type-check clean, Playwright 120 of 124 passed. The 4 failures were outdated test expectations: three download-label checks left from the Download menu change, and one test reading the chain name from the data strip. All 4 were updated and pass on rerun (8 tests in those specs).

**Screenshots:**
- `m17_rf00167_default.png`
- `m17_rf00167_zoomed.png`
- `m17_rf00167_collapsed.png`
- `m17_long_fixture_fit_all.png`
- `m17_long_fixture_zoomed.png`

## Position numbering (Family Overview, display only) — 2026-10-17

**What it does.** A "Position numbering" selector in the *Sequence and annotations* heading: **Rfam consensus** (default) or the author numbering of any RNA chain of the family. It is independent of the active instance (Instance selector) and changes labels only.

**Unchanged.** Rfam match positions remain the only internal coordinate:
- Selections, 2D/3D mappings, loops, insertions and aggregation all still run on Rfam positions.
- No backend, schema or data file was changed (checked: no file under `data/` or `rna_explorer/` is newer than the work).
- Downloads are untouched.

**Source.** The reference chain's existing validated mapping file (`mappings/<seq_name>.json`). The lookup is consensus column → match residue → `residue_name` + `author_residue_number` + `author_insertion_code` (`frontend/src/lib/numbering.ts`).
- No offsets and no re-alignment.
- A position with no match residue in the reference (deletion, outside the aligned region) gets no axis label and is never filled in.
- Rfam insertion residues never get a column.
- The active chain's mapping is reused when it is the reference; otherwise only the reference mapping file is fetched.

**Axis.** Labels are the reference's author numbers at the Rfam columns:
- They are round numbers at the zoom-dependent step.
- They are thinned so labels never overlap.
- Stretches with irregular numbering or insertion codes are filled with real mapped residues.
- The axis caption names the convention, e.g. "4V9D BA · author numbering".
- Zoom, scroll position and column geometry are unchanged by a switch.

**Tooltips.** The selected-chain tooltip shows these rows:
- Rfam position.
- Author residue: the active chain's own identifier, or Unavailable / Not present / Unmapped.
- `<ref> reference residue`, when the reference is another chain.
- Status: Coordinates available / Coordinates unavailable / Deletion / Outside aligned region.

The coverage, confidence, RF and SS_cons tooltips say "Rfam position N" and add the reference residue. +N insertion tooltips show the run's deposited author identifiers. The selection summary reads `Rfam position 50 → 2XNW A G46 · 1Y26 X reference G46 · highlighted and centred in 2D/3D`.

**Tests.**
- **`src/lib/numbering.test.ts` (6).** Synthetic fixture covering an insertion code, an author-number gap, a modified residue, a missing author number, a deletion, uncovered positions and an Rfam insertion residue. It also tests tick selection for continuous, irregular/insertion-code and 1,533-position numbering.
- **`tests/position-numbering.spec.ts` (3, real data; all expected values read from the mapping files).**
  - RF00167, first part: default numbering; reference 1Y26 X with active 4FE5 B; axis labels equal 1Y26 author numbers; zoom, scroll and column x unchanged; no 2D/3D reload, with only the reference mapping fetched.
  - RF00167, tooltip and summary checks:
    - Rows at 50, 64 (reference deletion) and 5 (reference outside the alignment).
    - Instance switched to 2XNW A with the reference kept; 19 is mapped without coordinates; 64 is a deletion; 5 is outside the aligned region.
    - The selection summary is correct and 2D still selects the residue.
    - When the reference is the active chain there is a single author row.
    - Switching back to Rfam leaves the selection unchanged.
  - RF00059: 9C6J A insertion codes (304, 306, 306C, 306F…) on the axis and in tooltips.
  - RF00177: 4V9D BA axis sparse at Fit all, follows zoom, no overlapping labels.
- **Updated wording** in `overview-tooltips`, `overview-tracks`, `overview-controls`, `annotations`, `single-selection` and `tracks-zoom`.
- **Results:** vitest 70, backend 70, Playwright 136 (134 in the full run; the 2 failures were old wording, fixed and re-run).

**Annotations tab (shared choice).** The choice lives in `state/numbering.tsx` (`NumberingProvider` above the tabs; one `NumberingSelect` control), so it survives tab switches and the reference mapping is loaded once.
- The Annotations tab shows the same selector beside its zoom controls.
- Its axis shows the reference's author numbers at the Rfam columns (same tick rule, shared zoom window).
- Mark tooltips and the pinned details add `<ref> reference residue: …`.
- Marks, aggregation and downloads stay on Rfam positions.
- Test: `position-numbering.spec.ts`, fourth test: shared choice in both directions; ticks equal the 1Y26 author numbers and sit on their Rfam columns; tooltip and pinned details; zoom; back to Rfam.
