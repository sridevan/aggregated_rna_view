# Structures tab: Structural Dissimilarity mode (2026-10-18)

The Structures tab has two modes, chosen with a segmented control at the top of the tab:

- **Structural Dissimilarity** (first button, default): the full geometric discrepancy matrix, with rows and columns in tree-penalized path length (TPPL) seriation order.
- **Clustering**: unchanged.

Without a stored seriation the tab shows Clustering only.

No structural comparison, clustering result or superposition was recomputed or changed.

## Seriation (backend, precomputed)

**Reference implementation:**
- `ordering_similarity.py` from https://github.com/sridevan/ordering_analysis (commit `37f63c9`), vendored byte for byte as `rna_explorer/ordering_similarity.py` (sha256 `1215f580…`, checked by a test).
- The algorithm is TPPL (Aliyev & Zirbel 2023, EJOR 305(2):617–629), adapted from BGSU RNA-3D-correspondence.
- It is called as in the reference `run.py`: `treePenalizedPathLength(dist, 10, 39873)`, which runs `orientPath` itself. The local `antibody_antigen_interface/scripts/build_instance_similarity.py` calls the same function with the same seed.

**Wrapper (`rna_explorer/seriation.py`).** It adds no logic to the algorithm. It only:
- sorts the instance identifiers before the call, so the order never depends on how the caller listed them;
- refuses an incomplete, negative, asymmetric or non-zero-diagonal matrix (see "Missing values" below);
- restores Python's global `random` state after the reference seeds it;
- requires scipy. Without it the reference silently drops the tree penalty and returns a different order, so scipy is now in a new `clustering` extra in `pyproject.toml`.

**Where it runs:** in `python3 -m rna_explorer clustering <ACC>`, on the discrepancy matrix (`metrics.discrepancy.seriation` in `structure_clustering.json`; schema extended and the field required). It stores:
- `order` (seq_names; the same order applies to rows and columns);
- `repetitions`;
- `seed`;
- `missing_pairs`;
- the reference record.

**Why precomputed:** the result depends on Python's seeded `random` module. Reimplementing it in the browser would not reproduce the reference order, and it would rerun on every load. The browser only reads the order.

**Regeneration:**
- The unmodified pipeline first reproduced all five families' `structure_clustering.json` exactly, apart from timestamps.
- After the change, the only difference is the new `seriation` record.

**Checks:**
- For every family, the stored order equals a direct call of the cloned reference on the same matrix.
- It also equals the order recomputed from the stored, rounded values.
- It is unchanged when the input is shuffled, and identical on repeat runs.

| Family | Instances | Labels on the Y axis |
|---|---|---|
| RF00167 | 76 | no |
| RF00162 | 29 | yes |
| RF00059 | 37 | yes |
| RF00005 | 60 | no |
| RF00177 | 60 | no |

**Missing values:**
- The clustering pipeline already rejects non-finite matrices (`check_matrix`), and none of the five families has a missing pair (`missing_pairs` = 0).
- The reference would impute a missing value as the largest observed value. A missing comparison is not a maximal difference, so the wrapper refuses such input instead of ordering on an invented value.
- In the browser, a missing or non-finite stored value stays `null`. It is drawn hatched, labelled "not measured" in the legend and tooltip, and never coloured.

## Frontend

**Mode state** (`state/structures.tsx`, above the tabs, so it is kept while the family page is open and across tab visits):
- `view` defaults to Structural Dissimilarity.
- `dissimDisplayed` holds this mode's own 3D content: `null` means the representative alone, otherwise one instance or one pair, never more than two. It is separate from the Clustering overlay (`displayed` and `hidden`), so neither mode leaves structures in the other.
- The 3D viewer caches coordinates per chain, so switching modes or structures downloads nothing twice.

**Default selection:** the first instance of the stored TPPL order (`seriation.order[0]`; 9LKF A in RF00167). It is table row 1, diagonal cell 1 and the only structure in 3D. The representative is not forced to the front. The family's designated representative only keeps its "Representative" label tag.

**Matrix:**
- `dissimilarityMatrix()` in `lib/structureGroups.ts` builds it once per loaded dataset, in the stored order, from `values_upper` and `dataset.instances`.
- It checks that the order lists every clustered instance exactly once.
- Switching modes neither refetches the file nor rebuilds the matrix (both checked in the browser tests).

**Layout in Structural Dissimilarity mode:**
- Hidden: the structural-groups panel, the dendrogram and heatmap panel (with its metric selector), Associated features (per group) and the table's Structural group column.
- The instance table uses the full width (the expanded layout).

**Table (`InstanceTable`, `seriation` prop):**
- **Default order ("Seriation (default)"):** exactly the stored TPPL order (`seriation.order`), applied to the whole filtered list before paging. The heatmap uses the same order.
- **Instances without a seriation position:** in RF00167 these are the 8 not-clustered instances, which have no pairwise values. They get no artificial position: they follow the ordered instances, by PDB ID and chain, with a "not in matrix" tag. They can still be selected (shown alone in 3D, no diagonal), and their metadata stays in the table. Nothing is listed twice or left out.
- **Column sorting:** a column header click sorts the table by that column, as before. The heatmap order, its values and the selection do not change.
- **Order by:** a compact select beside the search, "Seriation (default)" or the active column sort with its arrow, driven by the same sort state as the header arrows. Choosing Seriation restores the stored order without recomputing anything and keeps the selected instance or pair and the 3D view.
- **Sort note:** while a column sort is active, "Table sorted by Res. (Å) · Heatmap ordered by structural similarity".
- **Search and paging:** unchanged. A search keeps the sort mode, and paging never changes the selection or the 3D view.

**Selection, by seq_name only (never row or matrix positions):**
- **Table row or diagonal click:** that RNA alone in 3D. The row is selected and its diagonal checked, and any pair is cleared. A diagonal hover still shows the instance's metadata.
- **Below-diagonal click:** exactly that pair, replacing any previous one. The cell is outlined, both diagonals are checked and both table rows marked; the "Selected:" line and the toolbar read "A ↔ B".
- **Upper-triangle click:** nothing.
- **Selected instance hidden by a search:** the existing "hidden by the current search" notice; nothing else is selected.
- **Reset:** the first seriation instance alone, selected and checked, with the pair statistics cleared. The table sort, search and page are kept: Reset does not change the page even when that row is on another page.

**3D panel: one structure or exactly two, in both modes.** It is the same panel, viewer and superposition (the `orientation.json` transforms onto the family frame) in both modes.

- **One model drives it.** `StructuresTab` builds a `ViewerModel` holding:
  - the displayed structures (one or two);
  - a caption;
  - the stored pair value, if any;
  - which structures count as representatives;
  - what Reset does.
- **Title:** "3D structure" for one structure, "3D superposition" for two.
- **Labels:** below the viewer, read-only and colour-coded, in the viewer's colours. One structure reads "Displayed: ■ 6UC8 B · Representative"; a pair reads "Comparing: ■ A ↔ ■ B". Poor-fit flags stay on the labels.
- **Removed in both modes:** the member dropdown, "+ selected", the visibility checkboxes, the remove buttons and the per-structure emphasis.
- **Kept:** Fullscreen, Centre view, rotation and zoom.
- **Downloads:** the viewer caches coordinates per chain and cancels a load when the selection changes, so a quick series of clicks leaves only the last selection.

**State (`state/structures.tsx`):**
- **Clustering mode:** `displayed` is null (the selected group's representative) or one instance. A pair exists only through `compareGroup` (two representatives). The `hidden` list is gone.
- **Structural Dissimilarity mode:** `dissimDisplayed` is null (the family representative), one instance or one pair.
- The selected structural group (`groupByMetric`) and the app-wide selected instance stay separate from what is displayed.
- The old overlay helpers (`OVERLAY_MAX`, `groupCount`, `defaultDisplayed`) were removed.

**Clustering mode:**
- **Default:** the largest group selected and its representative alone, captioned "Structural group 1 · Showing representative".
- **Group selection** (groups panel, dendrogram or badge): that group's representative alone, ending any comparison.
- **Table row, heatmap cell or dendrogram leaf:** that RNA alone, captioned "Structural group 1 · Showing member 4FE5 B". Its group becomes the selected group. A non-clustered instance is shown alone with its reason.
- **Compare representatives** (unchanged control above the viewer): exactly the two representatives, never the member that was shown. The caption reads "Comparing representatives of structural groups 1 and 2", with the stored value. "Back to group" or another group returns to one structure.
- **Reset:** the selected group's representative alone, with the group kept.

**Structural Dissimilarity mode:**
- **Default:** the first seriation instance alone, with the caption "Select a pair from the heatmap to compare structures."
- **Table row or diagonal cell:** that RNA alone.
- **Below-diagonal cell:** exactly that pair, captioned "Comparing A and B", with the stored value.
- **Reset:** the first seriation instance.

**Heatmap (`components/structures/DissimilarityHeatmap.tsx`):**
- Same viridis scale as the clustering heatmap, from 0 to the family's maximum discrepancy, with a legend and a grey diagonal.
- Hover and click are handled once on the SVG and the cells are memoised.
- **Hover:** off-diagonal cells show the pair, the stored value and the shared positions. Diagonal cells show the instance's PDB ID and chain, title, method, resolution and organism; absent fields are omitted.
- **Diagonal click:** pins the same details. It closes on Escape, ×, an outside click or a scroll, and does not change the 3D view.
- **Below-diagonal click:** shows exactly that pair, replacing any previous one.
- **Upper-triangle click:** nothing.
- **Displayed structures:** checked on the diagonal, as a blue cell with a white check: the representative or chosen instance alone, or both instances of a pair. A pair's cell is also outlined, with its row and column in a light outline and the mirrored cell dashed.
- **Y-axis labels:** only for 50 or fewer instances.

## Tests

- **Backend: `tests/test_seriation.py`:**
  - the vendored file's checksum;
  - the reference's own test dataset (50 uniform points, seed 2276393) gives the reference order at 100 and at 10 repetitions;
  - determinism, a true permutation, independence of input order, and the row/column correspondence;
  - the global random state left untouched;
  - incomplete or invalid input refused;
  - each family's stored order valid and reproducible from the stored values.
- **Vitest: `src/lib/dissimilarityMatrix.test.ts`:** the order applied to rows and columns, each pair on the right instances, symmetry and diagonal, null for missing values (never the maximum), and rejection of an order that does not list every instance once.
- **Playwright: `tests/structures-dissimilarity.spec.ts` (6 tests):**
  - **Modes:** the default and the button order; the hidden panels and controls; Clustering's multi-structure overlay, members and remove kept intact across switches; the representative alone on returning; no refetch or rebuild; the mode kept across tabs.
  - **Default and sorting:**
    - row 1, diagonal 1 and the 3D content are the first seriation instance;
    - the full table order across all pages, with the not-in-matrix instances after;
    - row 5 gives diagonal 5;
    - PDB ID and resolution sorts change only the table, with the Order by value, header arrow and sort note consistent;
    - selecting after a sort checks the right diagonal;
    - search and paging keep the sort mode and the selection;
    - Seriation restores the exact order and keeps a selected pair;
    - no matrix rebuild.
  - **Matrix:** values and order against the JSON.
  - **Labels:** shown at 29 instances, hidden at 60 and 76.
  - **Selection sync:**
    - hover above, below and on the diagonal;
    - a diagonal click giving that RNA alone with its row selected;
    - upper-triangle click ignored;
    - two successive pairs, each exactly two structures in the one viewer, with title, labels, colours, the stored value, the cell and diagonal highlights, and both table rows marked;
    - a table row giving one structure, and a not-clustered row giving one structure with no diagonal;
    - Reset back to the first seriation instance, with the table sort and page untouched.
  - **Narrow window:** 760 px.
- **Playwright: `tests/structures-follow.spec.ts` (rewritten, 3 tests):**
  - Clustering default, another group, table members (replacing, never adding), a group's own representative and a non-clustered instance, each checked for viewer content, title, labels and colours, and caption;
  - Compare representatives (not the shown member), leaving a comparison by group change or Back, and Reset from a member and from a comparison;
  - rapid consecutive selections, and mode switches never showing more than two structures.
- **Other specs:**
  - existing Structures specs now select Clustering explicitly, and the missing-data specs check that no view switch is offered;
  - tests of the old overlay (members, visibility, remove, emphasis) now check the single-structure or pair equivalents;
  - the comparison-header height limit is now 80 px, because the header has the comparison caption as well as the metric line.

Screenshots: `docs/screenshots/m20_dissimilarity_default.png`, `m20_structures_dissimilarity_mode.png`, `m20_dissimilarity_hover_pair.png`, `m20_dissimilarity_diagonal_hover.png`, `m20_dissimilarity_labels_rf00162.png`, `m21_clustering_member.png`, `m21_clustering_compare.png`, `m21_dissimilarity_pair.png`, `m21_dissimilarity_single.png`, `m22_table_sorted_by_resolution.png`.

## Limitations

- **Metric:** only the geometric discrepancy matrix is seriated, as specified. TM-score has no Structural Dissimilarity view.
- **Groups:** structure colours in the 3D view and its labels still come from the structural group (as in Clustering mode). The matrix, its order and the table in this mode do not use the groups.
- **Instances shown:** the matrix contains only the clustered instances. Instances excluded from clustering have no pairwise values.
- **Large ribosomes:** for RF00177, superposing two 16S chains loads two large structures, exactly as the existing comparison does.
