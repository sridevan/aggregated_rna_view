# Milestone 2 validation report: RNA Structure Explorer frontend

Date: 2026-10-14. Backend data: `data/families/RF00167/` (unchanged since Milestone 1). Frontend: `frontend/` (React 19, TypeScript, Vite 8, Tailwind 4, pdbe-molstar 3.12.0, vendored PDBe RNA Viewer 0.3.0).

## What was built

- Three tabs: Family Overview (counts, methods, model length, coverage, consensus RF/SS_cons, representative with rule, Mol* of the representative), Structures (searchable/sortable instance table with organism, method, resolution, mapped region, sequence and structural coverage, average PP, Rfam bits; Mol* of the selected instance; three explicit "unavailable" panels for similarity, clustering/dendrogram and groups), Annotations (consensus tracks, instance selector, selection panel, PDBe RNA 2D viewer and PDBe Mol* side by side on desktop, stacked below `lg`).
- Shared selection state in React (`SelectionContext`): consensus positions are family coordinates, preserved when the instance changes; residues without a consensus coordinate (insert, local end, flank) are selectable per instance (`instanceResidues`) and cleared on instance change. Both viewers only ever receive residues with coordinates.
- Identifier handling: 2D viewer is fed the author chain (its PDBe data are keyed by author chain) and selected by label_seq_id; Mol* is queried by label_asym_id + label_seq_id; author numbers and insertion codes are displayed only.
- Feedback-loop avoidance: viewers are driven from React state with `doNotReport` (2D) and `visual.select` (3D, emits no click events); incoming clicks are translated and dispatched to the reducer, which re-applies to both viewers.
- Cleanup: listeners removed, 2D container emptied, Mol* plugin cleared after any in-flight load settles; one 2D instance mounted at a time (Annotations tab only).

## Automated checks (`npm run test:e2e`, Playwright, Chromium headless, 7 passed in 60 s; `npm test`, vitest, 3 passed)

| # | Requirement | Test | Result |
|---|---|---|---|
| 1 | PDB 1y26 loads in both viewers | `1y26 loads in both viewers…`: 2D status ready with 71 nucleotide glyphs, Mol* status ready, both reporting pdb 1y26 | pass |
| 2 | RF00167 consensus 50 selects polymer 34 | same test: selection panel shows `50 → G34 (auth 46, PP *)`; unit test also checks 9ndd/A (polymer 64) | pass |
| 3 | Both viewers highlight the correct nucleotide | 2D element `rnaview_1y26_34` carries the selection fill; Mol* received `{label_asym_id A, label_seq_id 34}` | pass |
| 4 | Clicking a nucleotide in 2D updates 3D | DOM click on `rnaview_1y26_20` emits `PDB.RNA.viewer.click` -> consensus 36 -> Mol* selection [20] | pass |
| 5 | Clicking a residue in 3D updates consensus | `PDB.molstar.click` event with pdbe-molstar's payload (`seq_id 10, auth_asym_id X`) -> consensus 26 -> 2D selection [10]; an event for chain B is ignored | pass |
| 6 | Switching instances preserves the selection | positions {50, 36} kept after switching 1y26 -> 4tzx; re-translated (4tzx: 50 -> G34 auth 46, Mol* [20, 34]) | pass |
| 7 | Deleted / uncovered / unobserved handled | 9nbh/A: 39 reported deleted, 50 reported "no coordinates", only 1 -> U16 highlighted (Mol* [16]); 6p2h/A: 1 reported uncovered, 39 and 50 highlighted | pass |
| + | Family statistics equal backend JSON; table rows equal instances.json (84); search filters; 3 pending panels | `family statistics…`, `structures table…` | pass |

Real-browser checks (Chrome, this session): consensus 50 and 36 highlighted in 2D (orange G) and 3D (violet, focused); a real mouse click on the Mol* canvas in the isolated experiment produced a genuine `PDB.molstar.click` (seq_id 34) that round-tripped to consensus 50 (see `RNA_VISUALIZATION_REVIEW.md`).

Screenshots: `docs/screenshots/m2_family_overview.jpg`, `docs/screenshots/m2_structures.jpg`, `docs/screenshots/m2_annotations_1y26_sync.jpg`.

## Defects found and fixed during validation

- 2D wrapper drove a viewer instance whose `render()` had not completed when instances were switched quickly (library accessed `this.options` -> React tree crashed). Fixed by tracking the rendered instance and guarding calls.
- PDBe RNA Viewer throws when the FR3D base-pair file is missing (e.g. 9nbh, 6p2h). The wrapper now fetches the two PDBe static files itself and shows an explicit "No PDBe 2D topology available" state; tests accept `ready|unavailable` for such chains.
- Mol* clicks are dispatched on the render target; the listener now sits on the wrapper element so bubbling and synthetic events both work.
- Creating a new pdbe-molstar plugin per structure on one element produced "Invalid data cell" rejections; the wrapper now keeps one plugin per mount and switches structures with `visual.update`, serialising loads.

## Known limitations

- Playwright's headless Chromium still logs two "Invalid data cell" page errors per Mol* load that do not occur in desktop Chrome and do not affect the tests; they look WebGL/headless related and were not chased further.
- 2D topology is available only for chains PDBe has processed (missing for the 2025 entries 9nbh, 9ndd and for 6p2h in this family); the UI states this explicitly.
- Insert/local-end residues are selectable per instance but not from the consensus tracks (they have no consensus coordinate); only via 2D/3D clicks.
- Node 20.14 is below Vite 8's stated minimum (20.19); dev and build work but print a warning.
- Structural similarity, clustering, Motif Atlas and ligand tracks, R2DT family-level 2D: not implemented (pending panels / later milestones).
