# Audit: "Structural groups" page in antibody_antigen_interface (read-only)

Repo: `/Users/sria/Desktop/antibody_antigen_interface` (git, branch main, HEAD `acd238a`).
App: `/Users/sria/Desktop/antibody_antigen_interface/app-react`. All paths below are relative to
`app-react/` unless absolute. Everything here was verified by reading code unless marked
**[not verified]**.

## 0. Corrections to the brief (important)

1. **The 3D view is NOT Mol\***. The structural-groups page uses **3Dmol.js 2.5.2** loaded from CDN in
   `index.html:8` (`https://cdn.jsdelivr.net/npm/3dmol@2.5.2/build/3Dmol-min.js`), via
   `src/states/SuperpositionViewer.jsx`. A Mol*/MVS viewer exists (`src/components/MolstarViewer.jsx`)
   but is not used by the groups page. Memory note `viewer-feedback-next-steps.md`: "The viewer is
   3Dmol.js ... NOT Mol* (Mol* is the kept production backup)."
2. **Clustering is computed in the browser**, not at build time (`src/groups/clustering.js`). The
   build script does emit an upstream `cluster_id`/`n_clusters` (auto_gap) but the groups page ignores it.
3. **Not TypeScript, not Tailwind.** Plain JSX + one hand-written stylesheet `src/styles.css` (1024 lines,
   no dark mode: 0 `prefers-color-scheme`).
4. **No charting library** for these panels: heatmap, dendrogram and overview are hand-written SVG / HTML
   table / canvas. (`recharts` is a dependency but used elsewhere, not by `src/groups`.)
5. **No tests at all** in the repo (no vitest, playwright, jest, pytest files; only a Pages deploy workflow).
6. There is **no per-structure visibility toggle** any more (removed in `ad026fb`).

## 1. Framework, tooling, organisation

`app-react/package.json`:
- deps: `react ^18.3.1`, `react-dom ^18.3.1`, `recharts ^2.12.7`
- devDeps: `vite ^5.4.0`, `@vitejs/plugin-react ^4.3.1`
- scripts: `sync-data` (`scripts/sync-data.mjs`, run as `predev`/`prebuild`), `dev`, `build`, `preview`.
- `vite.config.js`: base `/complex-interface-demo/` for build (GitHub Pages), `/` for dev, port 5173.
- 3D: 3Dmol.js via CDN global `window.$3Dmol` (index.html:8). No npm Mol*/pdbe-molstar.
- CSS: `src/styles.css`, BEM-ish prefixes `cs-` (similarity page), `sg-` (structural groups), `bs-` (block summary chips).

Routing: hash routes in `src/main.jsx` (no router lib). Groups pages, `main.jsx:48-93`:
`#kir22-groups` (11 asm, nGroups 3), `#atpsynthase-groups` (16, 3), `#human-hb-groups` (341, 6),
`#atcase-groups` (58, 2). Config object: `{ basePath, complexId, title, organism, similarityHash,
grouping: { method: 'Agglomerative hierarchical clustering', linkage: 'Average', nGroups } }`.
nGroups is hand-picked per complex (comments in main.jsx explain each choice).

Component organisation (`src/groups/`, plus shared pieces from `src/states/` and `src/components/`):

| File | Lines | Role |
|---|---|---|
| `src/groups/StructuralGroupsApp.jsx` | 557 | Page shell, state, table, layout, downloads, notebook hand-off |
| `src/groups/clustering.js` | 159 | UPGMA, tree cut, medoid, group colours, `pathBetween` |
| `src/groups/GroupPanels.jsx` | 292 | `GroupsCard` (group list + selected-group card + clustering details), `AssociatedFeatures` |
| `src/groups/features.js` | 66 | Per-entry feature counting (`groupFeatures`) |
| `src/groups/Dendrogram.jsx` | 149 | Leaf dendrogram in heatmap gutter (<=50) |
| `src/groups/GroupOverview.jsx` | 210 | Collapsed wedge dendrogram (>50) |
| `src/groups/RepresentativeComparison.jsx` | 108 | Compare-representatives control/panel/colours |
| `src/groups/ExportMenu.jsx` | 53 | Download dropdown |
| `src/groups/exportFigures.js` | 326 | SVG/PNG/CSV figure builders (drawn from data) |
| `src/states/DissimilarityHeatmap.jsx` | 764 | Shared heatmap (table <=50 / canvas >50), viridis, zoom |
| `src/states/BlockSummary.jsx` | 320 | `summariseSelection` (within/between stats), `ENRICH_PP`, notes |
| `src/states/SuperpositionViewer.jsx` | 293 | 3Dmol viewer, PCA orientation, fullscreen/centre, legend |
| `src/components/Hint.jsx`, `SortIcon.jsx`, `Pager.jsx` | | info-icon popups, sort glyph, pager |

State management: plain React `useState`/`useMemo` in `StructuralGroupsApp` (`selectedId`, `sort`,
`page`, `groupOnly`, `compareWith`, `comparing`, `nbOpen`, `panelSize`, `extrasH`), props drilled down.
No context, no store. One selection for the whole page, always a group id (`StructuralGroupsApp.jsx:142-146`).
Card heights synchronised via `ResizeObserver` (`Measured`, `onSize`, `extrasRef`).

## 2. Data model

Loading: one fetch, `StructuralGroupsApp.jsx:176-184`:
`${BASE}${basePath}/instance_similarity.json`; requires `heatmap.metrics.tmscore`.
CIFs: `${BASE}${basePath}/assemblies/<assembly_id>.cif` (`SuperpositionViewer.jsx:8`).
Staged into `public/<name>-similarity/` by `scripts/sync-data.mjs:73+` from `data/processed/<CX>/`.
(`api_data/` is the separate aggregated-interface API contract; not used by this page.)

`instance_similarity.json` (example `public/atpsynthase-similarity/instance_similarity.json`):
```
complex_id: 'PDB-CPX-106364'
reference_assembly: '6re1_1'                     # superposition frame
method: { score_type:'combined', linkage:'average', cut_mode:'auto_gap', n_clusters:4,
          source:'pdbe_complex_clustering/clustering.py (notebook defaults)',
          coordinate_source:'PDBe model server assembly coordinates', atoms_shipped:'backbone' }
subset: null
complex: { name, complex_portal_id, oligomeric_state, symmetry, representative:'6rdb_1', representative_in_set }
data_notes: [str...], data_notes_scope: 'mixed'
assemblies: [ { assembly_id:'6rdb_1', pdb_id:'6rdb', structure_title, resolution:2.8, exp_method:'em',
                cluster_id:2, packing_group:1, superposes_with_reference:true, mean_dissimilarity:0.6712,
                ligands:[{comp:'MG', copies:5, name:'MAGNESIUM ION'}],          # top 6 per assembly
                mutations:[{label:'K266R', type:'Conflict'}],
                modified:[...comp codes],
                transform:{ rotation:[9], translation:[3], source:'6rdb_1-6re1_1_transform.txt',
                            derivation:'direct'|'inverted'|'identity' } } ]
heatmap: { labels:[ids, matrix order],
           ordering:'tree-penalised path-length seriation (matrix only, pre-clustering)',
           dendrogram_order:[...], default_metric:'tmscore',
           coverage:{ pairs, metrics:{ shape:{pairs,shown}, tmscore:{...}, rmsd:{pairs,shown,per_pair} } },
           rmsd:{ label:'Backbone RMSD after superposition (US-align)', cell_label:'RMSD', unit:'Å', matrix:[[...]] },
           metrics:{ shape:{label,cell_label,order,matrix}, tmscore:{label:'Structural (1 − TM-score)',
                     cell_label:'1 − TM-score', order:[seriation], matrix:[[n×n]]} } }
```
What the groups page reads: `heatmap.labels`, `heatmap.metrics.tmscore.matrix`, `heatmap.rmsd`,
`assemblies[*]` (ids, pdb_id, title, resolution, exp_method, ligands, mutations, modified). It does
NOT read `cluster_id`, `dendrogram_order`, `metrics.tmscore.order`, or `transform`.
**Leaf order, linkage tree, group assignments, cut height, medoids and within/across stats are all
computed client-side** (see below). There is no separate linkage/dendrogram file.

Generator: `/Users/sria/Desktop/antibody_antigen_interface/scripts/build_instance_similarity.py` (1521 lines).
- Inputs: per-pair US-align sidecars in gitignored `PDB-CPX-<id>/<asmA>/` dirs:
  `<A>-<B>_scores.csv` (`length_structure_1/2, aligned_length, rmsd, tm_score_1, tm_score_2`),
  `<A>-<B>_transform.txt` (TM-align-format rotation), `<A>-<B>_cc.csv` (chain correspondence).
  Produced upstream (not by this repo) **[how US-align was invoked, e.g. -mm flags, not verified]**.
  Shape (Zernike+spectral) scores come from `~/Desktop/pdbe_complex_clustering` (clustering.py,
  merged TSV / DuckDB).
- `tmscore_matrix()` (:492-527): `d = 1 - (tm_score_1 + tm_score_2)/2`, i.e. **average of the two
  length normalisations**, clamped >= 0; refuses incomplete matrices (returns None).
- `rmsd_matrix()` (:530): US-align RMSD, 2 dp.
- `seriation_order()` (:578): notebook's tree-penalised path-length seriation (used by similarity page only).
- `cluster_matrix()` (:395-410): auto_gap cut via notebook `find_optimal_num_clusters`, clamped to
  [2, n-1]; emitted as provenance only (spec says auto_gap gave k=4 on ATCase vs true 2, k=18 on human Hb).
- Superposition: `transform_to_reference()` (:618) reads the direct file or inverts the reverse one
  (R' = R^T, t' = -R^T t); `write_superposed()` (:699) writes **pre-transformed, backbone-only
  (N/CA/C/O/OXT, no HETATM), first-model-only** CIFs. Coordinates fetched from PDBe ModelServer with
  RCSB fallback (`fetch_assembly`, :640). Reference = medoid of the shape matrix (`--reference auto`)
  or given.
- Annotations: ligands from the assembly CIF HETATM (`ligands_of`, :759, top 6, waters excluded);
  mutations/modified from PDBe `mutated_AA_or_NA` / `modified_AA_or_NA` POST endpoints
  (`entry_annotations`, :1101), label regex `[A-Z][0-9]+[A-Z]`, author numbering with fallback.

Client-side clustering, `src/groups/clustering.js`:
- `averageLinkage(matrix)` (:15-58): naive O(n^3) UPGMA over a flat Float64Array of 2n-1 nodes;
  strict `<` tie-break (first pair met) for determinism; children ordered larger-first, then lowest
  index (:42-45), which fixes the leaf order (largest group at top).
- `cutTree(root, k)` (:63-79): undo the k-1 highest merges; `cutHeight` = midpoint between the lowest
  merge undone and the highest kept (this is the "cut at 0.094" figure).
- `medoid(members, matrix)` (:83-93): member with minimum mean distance to other members (ties: first).
- `structuralGroups(labels, matrix, {nGroups})` (:117-139): returns `{root, order, groups, groupOf,
  cutHeight, maxHeight}`; groups sorted by size desc, ids 1..k, `from/to` = leaf-order span,
  `representative` = medoid.
- `pathBetween(root, a, b)` (:144-159): node set on the tree route between two nodes (compare mode).

Within/across stats: `summariseSelection` in `src/states/BlockSummary.jsx:94-166`:
`within` = mean over ordered member pairs excl. diagonal; `withinMax`; `between` = mean over
member x non-member; `rmsdWithin`/`rmsdMax` from the RMSD matrix; `entries` = distinct pdb_id.
Ratio in `GroupPanels.jsx:94-98` = between/within **using the 4-dp rounded figures** (commit d4d22c2).

## 3. Structural groups (representation, selection, colours)

- Group object: `{ id, name:'Group N', color, node, members, from, to, representative }`.
- Default selection = largest (id 1): `StructuralGroupsApp.jsx:202-203` (`selectGroup(1)` on clustering
  change); card text "Largest group selected by default." (`GroupPanels.jsx:148`).
- Colours: `clustering.js:112-113`: `GROUP_COLORS = ['#E8685A','#B0327A','#8C4A1E','#3FA7D6']`
  (coral, magenta, rust, light cyan) given to groups with >= 2 members in size order;
  `UNPLACED = '#6b7480'` grey for every singleton and for multi-member groups past the 4th.
  Rationale comments: must not resemble viridis (OKLab dE checks), validated with the dataviz palette
  validator; never used for heatmap cells.
- "small" rule: `SMALL_GROUP = 5` (`GroupPanels.jsx:30`); tag "single" if 1 member, "small" if < 5
  (:133-137); amber note in card with number of pairs (:164-173); `AssociatedFeatures` gets
  `small` and adds "read the counts rather than the percentages".
- "largest" tag on group 1 (`LARGEST_NOTE`).
- Folding: `MAX_LISTED = 3` (:25); folds only if groups > 4 (:69); selected hidden group is appended to
  the folded list; fold never opens automatically; label (:80-88) e.g. "3 more groups · all single
  assemblies" / "2 single, 1 small".
- Selected card: name + dot, "N assemblies · M PDB entries", Representative (mono id + "Medoid"
  badge), Within-group (mean 4dp, max 3dp), Across-group (mean 4dp, "X.x× the within-group mean"),
  `<details>` "Clustering details" (Method, Linkage, Groups, Input + provenance note).
- Table (`StructuralGroupsApp.jsx:261-281, 361-438`): `PAGE_SIZE = 10`; default order = leaf order
  ("Order: matching the dendrogram"); 3-state column sort (asc, desc, off) with leaf order as tiebreak;
  "match the dendrogram again" link; on group select, jumps to page containing `g.from` (:200);
  rows of selected group get `sg-row-sel`, rows shown in 3D get a colour swatch; medoid tag `cs-rep`;
  group cell is a clickable chip; assembly id links to PDBe entry page. "selected group only"
  checkbox exists **only when n > 50** (`overview`, :369-379).

## 4. Heatmap + dendrogram

- `DETAIL_MAX = 50` (`StructuralGroupsApp.jsx:50`); `overview = nAll > DETAIL_MAX` (:234); <= 50 renders
  `DissimilarityHeatmap` with `leftPanel` = `Dendrogram` (:490-506); > 50 renders `GroupOverview`
  (:507-515). Separate constants with the same value: `MAX_LABELLED = 50` in
  `DissimilarityHeatmap.jsx:76` (labels + DOM table vs canvas) and `HEATMAP_EXPORT_MAX = 50` in
  `exportFigures.js:23`.
- Heatmap tech: <= 50 an HTML `<table>` with per-cell `<td>` (aria labels, roving tabindex);
  > 50 a canvas with a 1-px-per-cell offscreen tile (:375-486). On the groups page only the table path
  is reached. Cell size solved from card width (:161-191). Props used by groups page:
  `dragSelect={false}`, `pickable={() => false}`, `onPick` no-op, `bands` (white outline per group,
  CSS `.cs-hm-band`), `block` (dark outline of selected group, `.cs-hm-mark`), `leftGutter=124`,
  `colorOf` (row labels/rings for structures in 3D).
- Colour scale: viridis, 10 matplotlib stops linearly interpolated (:35-50), scaled 0..max of the full
  matrix; diagonal neutral grey `#d7dae0`; legend ramp with numeric end ticks and "More similar ←
  dissimilarity → More different" + info hint (:669-698). Known limitation (spec): min/max stretch;
  percentile clipping proposed, not implemented.
- Hover: tooltip with pair ids, 1 − TM-score (3 dp) and RMSD; diagonal shows title/method/resolution.
  Cells not clickable on the groups page.
- Zoom: Zoom − / fit / + (factor 1.6, max 40 px cell), centre-preserving scroll (:244-262). Drag-to-drill
  disabled on groups page.
- `Dendrogram.jsx`: SVG; root left, leaves right aligned to rows; linear x in merge height; 8 px leaf
  stub; grey branches, selected group's subtree in ink; dashed cut line with `<title>` "Cut at …";
  clickable/keyboardable group bars (16 px) in group colour with rotated name if tall; axis shows
  only 0 and max; compare mode marks the two leaves and draws the tree route in ink and dims the rest.
- `GroupOverview.jsx` (>50): tree drawn down to group nodes only; each multi-member group a
  `<polygon>` wedge whose point is at the group's top merge height; singletons a line; row height
  clamped `MIN_ROW 34` to `MAX_ROW 132`, stretched to `MIN_ROWS_TOTAL 300`; colour strip, name +
  count, proportional bar (n/total); whole row clickable; dashed cut line with text
  "cut at {fmt(cutHeight)}" (:156-163); axis label "{quantity} at merge"; compare mode marks rows and route.
- Downloads (`StructuralGroupsApp.jsx:310-317`, `exportFigures.js`): Dendrogram SVG/PNG (always full
  labelled horizontal tree, even for 341), Assembly instances CSV
  (`<CX>_structural_group_assemblies.csv`: assembly_id, pdb_id, structural_group,
  group_representative, structure_title, experimental_method, resolution_angstrom), Pairwise
  dissimilarity CSV (`<CX>_pairwise_tm_dissimilarity.csv`, long format: assembly_a, assembly_b,
  tm_dissimilarity, tm_score, rmsd_angstrom, group_a, group_b), Heatmap SVG/PNG only if n <= 50.
  Figures are generated as SVG strings from data (not DOM capture); PNG via Image→canvas at 3x, max
  side 12000. Failures surface via `onError` notice.
- Notebook: "Explore clustering in notebook" button toggles a **placeholder** box (notebook repo is
  private) with "Download the inputs (JSON)" hand-off: `{complex_id, input:{measure,label,labels,matrix},
  clustering:{method, linkage, n_groups, cut_height}, result:{leaf_order, groups:[{name,
  representative, members}]}}` (:287-306, 467-474).

## 5. 3D alignment view

- Component: `src/states/SuperpositionViewer.jsx`, 3Dmol.js. Props `{basePath, entries:[{assembly_id,
  color}], height}`.
- Coordinates: fetched local static CIFs, **already superposed at build time** onto one reference
  assembly (TM-align/US-align transforms; backbone only). No runtime superposition, no matrix applied.
  Text cached per assembly_id in a ref Map.
- Orientation: PCA of CA coordinates of the first loaded structure (Jacobi eigen, :27-98), applied by
  rewriting mmCIF text (`rotateCif`, :105-119), reused for the whole page; sign pinned deterministically.
- Style: `cartoon: { color, style:'oval', ribbon:true, arrows:false, thickness:0.4, width:0.4 }`
  (one colour per model). `fit` = `zoomTo()` + `zoom(0.82)`; ResizeObserver refits.
- What is shown (`StructuralGroupsApp.jsx:248-259`): `OVERLAY_MAX = 8` (:39). Group size 2..8: all
  members, representative first, colours from `SERIES` (8 colours, first six Okabe-Ito, :42-43);
  group > 8 or singleton: representative only in the group colour. No manual toggles.
- Compare representatives (`RepresentativeComparison.jsx`): select another group (default first other
  group), button switches viewer to the two medoids in group colours (grey groups drawn `#2E3440`);
  panel shows TM-score, 1 − TM-score (from matrix) and RMSD (US-align pairwise); "← Back to group view";
  any group selection exits compare. Caveat (memory): overlay uses each structure's build-time fit to the
  common reference, not a pairwise refit, while the RMSD shown is the pairwise US-align value.
- Buttons: Fullscreen (Fullscreen API on wrapper, resize on change), Centre view (resize + fit).
- Legend: "Displayed structures" overlay list with swatches.
- Cleanup: effect sets `cancelled` flag only; **no viewer dispose on unmount** (viewerRef never
  cleared; 3Dmol canvas left to GC with the DOM node). Fullscreen listener and ResizeObserver are removed.

## 6. Associated features

- Code: `src/groups/features.js` (`featuresOf`, `groupFeatures`), render `GroupPanels.jsx:227-292`.
- Sources: per-assembly `ligands` (CIF HETATM comps, top 6 by copies; note ligands_of truncates to 6),
  `mutations` (PDBe mutated_AA_or_NA), `modified` (PDBe modified_AA_or_NA).
- Counting: **per PDB entry**; an entry's features = union over its assemblies; an entry with
  assemblies on both sides counts on both sides (features.js:3-6, 32-38). nIn/nOut = entries inside /
  outside (outside = all other groups pooled); assembly counts also reported in the note when different.
- Shares: `inShare = a/nIn`, `outShare = c/nOut`, `delta = in − out`.
- Threshold: `|delta| >= ENRICH_PP = 0.15` (`BlockSummary.jsx:65`), display only, sorted by |delta|.
  Ligands show both directions; mutations and modified residues only `delta > 0` (features.js:60-62).
- ▲/▼: CSS pseudo-elements, `styles.css:744-747` (`.bs-up::before '▲ '`, `.bs-down::before '▼ '`),
  class chosen by `delta > 0`.
- Exclusion: position-1 substitutions `/^[A-Z]1[A-Z]$/` (initiator Met) removed; count of affected
  assemblies shown "· N with a position-1 substitution, excluded".
- Empty states: "nothing differs by 15% points or more"; modified "none deposited for this complex";
  `noRest` → "This group is the whole set".
- Evidence display: hover title with ligand full name and "a/nIn entries in this group / c/nOut entries
  outside"; help popups LIGAND_NOTE / MUTATION_NOTE / MODIFIED_NOTE. No significance tests (by user
  decision).

## 7. Tests

None. No vitest/jest/playwright configs or spec files anywhere in the repo; `.github/workflows/deploy.yml`
only builds and deploys. Memory: "features.js, no tests"; verification was done by driving the page
manually / ad hoc Playwright scripts **[not committed]**.

## 8. Session history and design decisions (relevant only)

Sources: `~/.claude/projects/-Users-sria-Desktop-antibody-antigen-interface/ca58ea99-...jsonl`
(main groups session, 2026-09-29 to 10-01), `memory/structural-groups-prototype.md`,
`memory/no-statistical-interpretation-layer.md`, `memory/similarity-measured-limits.md`,
`memory/ui-text-conventions.md`, local gitignored `spec_revised.md` ("Assembly-similarity view" §§).

Git (`git log --oneline`): `065e0f2` prototype pages (Kir2.2, ATCase, human Hb) → `d4d22c2` card figures
agree with ratio → `1553995` small flag, per-entry features, drop group-only view, pairs CSV →
`ad026fb` 3D view fixed by the group (8 max, no toggles) → `930e169` grey = unplaced, manual fold →
`acd238a` F-ATP synthase page.

Design decisions:
- Terminology: always "Structural group", never "conformational state"/"open/closed"; heading
  "Associated features" not "Enrichment" (memory). GROUP_NOTE sentence reused everywhere.
- TM-score only on the groups page; no measure switch (`MEASURE='tmscore'`, :58).
- Clustering deliberately temporary: UPGMA, k fixed per route; user: "Do not spend time trying to find
  the mathematically optimal clustering strategy." Re-cutting/other methods → notebook.
- One selection model: group only; no drag, no manual subsets (user, 2026-09-29).
- Groups card holds group metadata and separation only; annotations only in Associated features
  (user removed "Associated" row and "Bound macromolecules" row, d4d22c2).
- >50 rule (user 2026-09-29/30): no N×N matrix on the page above 50; wedge overview instead; downloads exempt.
- k for human Hb = 6 because average linkage peels outliers first (k=3 → 337/3/1; k=6 → 174/161/3/1/1/1).

Revisions / abandoned ideas:
- 3D: originally "show all if <= 10, else medoid, with manual toggles within group" (user 19:09 on
  09-29) → reduced to 8 colours (validator) and toggles removed (ad026fb; user: "anything more should
  be notebook exploration").
- "Group N only" within-group matrix/tree view for small groups of large complexes: built, then
  removed 2026-09-30 (user: "users could use notebook for this"). Stale comment remains in
  `main.jsx` atcase route ("open in full in the group-only view").
- Fisher exact + Benjamini-Hochberg graded chips (supported/suggestive/weak): proposed with numbers
  (e.g. Kir2.2 PIO 5/6 vs 0/5 p=0.015, q=0.21), **declined** by user ("lets not help to interpret it
  but use pdb-entry counting").
- Similarity page: earlier "five superposed structures" cap (spec) — different page.
- Bound macromolecules row: added, then removed (no dataset carries the field).

Known issues / unfinished:
- Associated-features weaknesses the assistant listed on 2026-09-30 07:49 and the user did not act on
  (in session): two sides not comparable (1/3 vs 40/295); fixed 15-pp rule favours small groups (a
  1-of-3 always clears); "outside" is pooled, two-group pages show mirror chips; absence is
  unexplained (shared vs unannotated); framing suggests cause. Proposed fixes (counts on the chip,
  minimum count >= 2, collapse mirror chips, "all entries share ligands") **not implemented**.
- Notebook link is a placeholder (private repo).
- Study-level redundancy cannot be normalised (no study id).
- Viridis min/max stretch compresses (human Hb middle 50% in 13% of ramp; Pol II outlier) —
  percentile clipping not implemented.
- `ligands_of(top=6)` truncates ligand lists, so rare ligands in large assemblies can be missed
  (inference from code; not discussed in session).
- Compare overlay uses fits to the common reference, not pairwise refits.
- Whole-assembly TM-score cannot resolve substates (ATP synthase rotor dilution; human Hb T vs T-high).
- O(n^3) browser UPGMA; comment says fine at n<=50 but it also runs for 341 (ms, after a Float64Array fix).
- Two Claude Docs (scientific questions; UX decisions for UX specialist and team lead) were written
  in the session **[content not read; links not verified]**.

## 9. Reuse recommendations for the RNA app (React 19 + TS + Tailwind 4 + pdbe-molstar 3.12)

The RNA app has per-family `instances.json` (84 chains for RF00167), consensus mappings, Kabsch C4′
orientation onto a display representative (`orientation.json`), and a Structures tab placeholder;
no pairwise matrix exists yet (PROJECT_HANDOVER.md §6).

| Component | Recommendation | Changes needed |
|---|---|---|
| `clustering.js` (UPGMA, cutTree, medoid, pathBetween) | **Reuse directly** (port to TS, add vitest) | Type the node; consider moving clustering to the Python backend so groups, cut height and medoids are in a schema-validated JSON with provenance (RNA project rules: provenance on every output, backend tests). Keep the browser version only for display-time recomputation, or drop it. Tie-break determinism must be kept. |
| Group colours / UNPLACED rule | **Reuse** | Re-check against the RNA app's existing palettes (track colours, 2D viewer) and dark mode; define as CSS tokens. |
| `GroupsCard` | **Adapt** | Unit is chain/instance (label_asym_id) not "assembly"; "PDB entries" count stays relevant; numbers need explicit numerator/denominator; within/across from the RNA distance (and its definition, e.g. intersection of observed consensus positions per pair). Rewrite in Tailwind. |
| Instances table | **Adapt / reimplement** in TS+Tailwind | Columns: PDB, chain (auth for display, label internal), organism/taxid (memory: species per chain), ligand, resolution, method; link to the RNA app's own instance selection rather than PDBe complex tab. Integrate with existing `state/selection.tsx` (`group` selection kind already planned). |
| `DissimilarityHeatmap` | **Adapt** (extract the table/canvas core + viridis + zoom; drop drag-drill, pickable, measure toolbar) | 84 chains for RF00167 is > 50, so on RF00167 the overview path would apply; consider whether the RNA tab wants the full matrix (canvas path exists and is solid) — decision needed. Remove complex-specific notes. |
| `Dendrogram.jsx` | **Reuse with light adaptation** | Pure SVG, prop-driven; port to TS; label "assemblies" → instances. |
| `GroupOverview.jsx` | **Reuse with light adaptation** | Same; wording. |
| `exportFigures.js` + `ExportMenu` | **Reuse** | Rename columns (`assembly_id`→instance id with label_asym_id + author chain), add provenance header rows (Rfam release, pipeline/config versions) to CSV/SVG. |
| `SuperpositionViewer` (3Dmol) | **Reimplement** on pdbe-molstar | RNA app already has `PdbeMolstarViewer.tsx` with `TransformStructureConformation` and camera fit; load multiple chains (PDBe/BinaryCIF URLs or local) and apply per-instance Kabsch transforms from `orientation.json` (or pairwise ones) at runtime rather than shipping pre-transformed backbone CIFs. Per-structure colour via overpaint/uniform colour theme; cap 8 with the validated SERIES palette; explicit `plugin.dispose()` on unmount (antibody viewer has none). Representation: cartoon/tube for RNA backbone. PCA orientation not needed (orientation already solved). |
| `RepresentativeComparison` | **Adapt** | Logic is generic (pair lookup, colours); only the viewer hook changes. Make clear which superposition the overlay uses (common reference vs pairwise). |
| `AssociatedFeatures` / `features.js` | **Adapt substantially** | Replace protein features: ligands → RNA app's classified ligands (`ligands.json`, annotation_config classes, exclude additives by class not by display threshold only); mutations (UniProt/PDBe mutated_AA, initiator-Met rule) → irrelevant for RNA (drop, or replace with sequence variants at consensus positions); modified residues → `modifications.json`; add apo/holo, organism, construct (scaffold) and method. Must distinguish "not annotated"/"source unavailable" from absence and show numerator/denominator on the chip, not only on hover (also addresses the open weaknesses above). Counting per PDB entry is reusable. |
| `summariseSelection` | **Reuse the stats part** | Drop ligand/mutation enrichment parts; keep within/between/max. |
| Build pipeline (`build_instance_similarity.py`) | **Reimplement** | Protein/assembly-specific: PDBe complexes, assemblies, PISA-free but assembly ModelServer fetch, packing groups, US-align sidecars from an upstream corpus. For RNA: compute pairwise US-align (`-mol RNA`) or Kabsch/RMSD on consensus-matched C4′ from the validated mappings, with explicit per-pair overlap denominators; record TM normalisation choice (average of both, as here, or normalised by consensus length). |

Protein/antibody/assembly assumptions to remove: "assembly" vocabulary and `<pdb>_<asm>` ids,
PDBe complex links, `exp_method` vocabulary (fine), UniProt-based mutation labels and position-1 rule,
complex `data_notes`, oligomeric state, packing groups, backbone N/CA/C/O atom filter (RNA would use
P/C4′ or all-atom), CA-based PCA orientation.
