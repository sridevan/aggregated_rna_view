# Milestone 3E validation: predicted RNA binding pockets (fpocketR) in Annotations

Date: 2026-10-16. Frontend integration of the completed fpocketR analysis (`experiments/fpocketr_full_RF00167/REPORT.md`, `data/families/RF00167/pockets.json`, `pocket_tracks.json`). No rerun of fpocketR, no change to pocket calculations, SASA, loop or ligand annotations, mappings or the clustering worktree.

## Implementation
- `api/data.ts` `loadPocketTracks` (optional, memoised; a family without the file shows no category), `state/useAnnotations.ts` loads it with the other annotation datasets, `types/backend.ts` `PocketTracksDataset`.
- `AnnotationsTab.tsx`: new collapsible category **Predicted binding pockets** (placed between Ligand binding sites and Solvent accessibility) with exactly one track **Pocket occurrence**, built from `pocket_tracks.json` → `tracks.all_score_ge_0` (fpocketR's published default filter, recommended by the family-wide report; the top-ranked and Score ≥ 0.4 variants remain in the downloadable file). Marks use the shared fixed geometry and a new purple frequency `ScaleDef` (legend "Predicted pocket occurrence", five bins); value = `n_supporting_instances / n_evaluable_instances`. Grey marks when no instance is evaluable at the position; empty when evaluated and in no pocket. Tooltip: "Rfam residue N / In a predicted pocket in k of n evaluable RNA instances". Pinned details: supporting entries, distinct pockets, fpocket Score range, unevaluable count and reasons, filter. Category ⓘ: fpocketR / fpocket 4.0.3 with published RNA parameters on ligand-free chains; geometric prediction, not experimental binding; zero = no pocket passed the size criteria.
- Download menu: "Predicted pockets, per-structure (JSON)" (`<ACC>_pockets.json`) and "Predicted pockets, aggregated (JSON)" (`<ACC>_pocket_tracks.json`), HEAD-checked like the others; seven items in total.

## Validation (RF00167)
- `pocket-track.spec.ts`: single track, category order, mark count = positions with support, geometry identical to other tracks, colour monotonic with frequency, data-value equals the file's frequency, legend, tooltip text with 61 of 80 at consensus 56, pinned details, popover content and non-collapsing click, no track-level icon, no grey marks (every position evaluable), downloads with valid JSON for the family, loop/ligand/SASA tracks unchanged. Updated `annotation-tracks.spec.ts` counts (4 categories, 4 header icons, 7 download items) and `screenshots.spec.ts`. Screenshot: `docs/screenshots/m3e_annotations_pockets.png`.
- Backend: 41 pytest (incl. `test_pockets.py`, `test_pockets_rf00167.py`).

## Limitations
Only the default filter is visualised; 19 of 84 instances contributed zeros (evaluated, no pocket), which the popover explains; the track reflects ligand-stripped deposited chains and fpocketR's size criteria (purine pocket at the 42-alpha-sphere edge), so a zero is not evidence of a missing cavity.

## Popover refinement (2026-10-16, later)
- Category ⓘ texts shortened. Solvent accessibility: "Per-nucleotide SASA calculated for isolated experimental RNA chains and aggregated by Rfam position. Colours indicate median SASA (Å²), reflecting nucleotide exposure." Predicted binding pockets: "RNA pockets predicted using fpocketR and mapped to Rfam positions. Colours indicate pocket occurrence frequency. No prediction does not necessarily mean a binding pocket is absent." Versions, parameters and caveats stay in the milestone documents and the downloadable JSON. Loop and ligand texts unchanged.
- `InfoTip.tsx` (shared by all category headers and the Family Overview cards) now renders its popover through a React portal into `document.body` with fixed positioning: opens below the icon, flips above when the space below is insufficient, is clamped horizontally with 14 px viewport padding, max width 320 px, max height = available space with internal scrolling, and is repositioned on window resize and on any scroll (capture listener). Hover/focus show it transiently; click pins it (second click, Escape or an outside click closes it); clicks never reach the category toggle.
- Tests: `info-popover.spec.ts` (texts; below placement at the top; flip above when the icon sits near the bottom; adjacency to the icon on the reported side; horizontal clamp at the right edge; resize and scroll repositioning; collapsed category; pin/unpin/outside/Escape; 390×500 viewport with internal scrolling and keyboard focus). Screenshots: `docs/screenshots/m3e_popover_pockets.png`, `m3e_popover_sasa.png`.

## Terminology refinement (frontend only)

User-facing "evaluable/unevaluable" wording was removed from the Annotations tab. Backend fields, schemas, datasets, calculations, colours, geometry and category order are unchanged; internal identifiers such as the `unevaluable` mark state remain.

| Place | Before | After |
|---|---|---|
| Loop tooltip | `Type: Hairpin loop (HL)` / `Observed in 66 RNA instances` | `Hairpin loop (HL)` / `Observed in 66 of 80 RNA instances (82.5%)` |
| Ligand tooltip | `Ligand: HPA — Hypoxanthine` / `Observed in N PDB entries` | `HPA — Hypoxanthine` / `Observed in N of M PDB entries (x%)` (M = PDB entries containing the bound ligand, the existing denominator) |
| Pocket tooltip | `In a predicted pocket in 61 of 80 evaluable RNA instances` | `Pocket occurrence: 61 of 80 RNA instances (76.3%)` |
| SASA tooltip | unchanged | `Median SASA: … Å²` / `Based on N RNA instances` |
| Pinned details | "observable, source-evaluated instances", "instances not evaluable here", "entry-entities observable and evaluated" | plain counts, "Not counted: k RNA instances (no coordinates: k)", plus "Frequencies are calculated from instances where annotation presence or absence can be determined." and the source-specific rule (BGSU overlap; successful fpocketR runs without a pocket count as absent, failed runs and missing residues are not counted; ligand counting per PDB entry) |
| Legend | "no evaluable structure / …", "evaluated, not annotated" | "presence or absence cannot be determined at the position (SASA: no qualifying measurement)", "annotation absent at the position" |

Numerators and denominators in tooltips are the same fields that set the colours (`n_positive_observable_and_evaluated`/`n_observable_and_evaluated` for loops; `n_supporting_instances`/`n_evaluable_instances` for pockets; supporting entries / `n_entries_in_default_track` for ligands). Percentages use at most one decimal (`pct`, vitest `src/components/AnnotationsTab.test.ts`). Playwright `tests/terminology.spec.ts` checks that no banned wording is visible in the tab, tooltips, pinned details, popovers or Download menu, and that `data-value` still equals the backend frequency.
