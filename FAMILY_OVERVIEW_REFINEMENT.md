# Family Overview and 3D viewer refinements: validation report

Date: 2026-10-15. Frontend only; the Rfam consensus mapping, Infernal alignment, orientation (Kabsch) dataset and annotation data model are unchanged.

## Layout (Family Overview)

| Row | Content |
|---|---|
| 1 | "About this family" (full width): Rfam long description, link to the Rfam entry |
| 2 | Six summary cards in one row on desktop (2 / 3 / 6 columns at narrow / medium / large widths), values and counting conventions unchanged |
| 3 | Molecular view card (chip MOLECULAR VIEW, heading "Representative RNA — 4FE5 chain B", or "Selected RNA — …" when another instance is selected): full-width Mol* viewer, 560 px tall, Rfam-mapped region by default, white background, Structure Tools hidden, RNA-only / complete toggle, Reset orientation / Fit and Focus selection |
| 4 | 40 % column: Representative structure (chip STRUCTURE; prominent PDB id, author chain and italic species; compact label/value list for method, resolution, Rfam model coverage, matched residues with coordinates, PDBe 2D topology; selection rule in an info tooltip) above Experimental methods (chip STATISTICS; bars styled like the species chart, counts and percentages of chain instances, conventions in a tooltip). 60 % column: Species distribution (chip DIVERSITY; toggle top-right; "5 distinct identified species · 60 of 66 PDB entries have species-level assignments"; provenance and counting conventions in a tooltip; labelled horizontal bars 14 px high with counts at the end and percentages on hover; legend directly beneath). Cards are top-aligned and sized by content (the equal-height requirement was withdrawn) |
| 5 | Consensus model (chip REFERENCE), full-width, collapsible, closed by default |

All sections use the shared `SectionCard` (eyebrow chip, 18 px semibold heading, description, right-aligned actions, collapsible with `aria-expanded`). Chips are neutral grey, uppercase 10 px. On narrow screens the stat cards wrap, the viewer keeps its height and scales in width, the metadata and diversity cards stack, species labels stack above their bars, and there is no horizontal overflow (tested at 390 px).

## Camera framing

Previous behaviour: Mol*'s reset places the camera at a distance derived from the bounding *sphere* radius (`Camera.targetDistance`), so an elongated RNA occupied only ~36 % × 54 % of the canvas. Two further problems were found: `Camera.Reset` discards the position in a supplied snapshot and re-derives it from the radius, and the loader's own reset request is resolved on a later frame, overriding anything applied immediately after the load.

New behaviour (`PdbeMolstarViewer.resetView`):
1. After the structure, the precomputed orientation transform and the visible representation (region or complete) are ready, every *visible* atom is projected onto the fixed view axes of the reference frame (right, up, towards-camera); hidden components never enter the calculation.
2. The camera target is the centre of the projected bounding box; the distance is the exact perspective requirement: the maximum over atoms of `max(|y|/tan(fov/2), |x|/(tan(fov/2)·aspect)) + depth-towards-camera`, with a 6 % margin, using the camera's fov and the real canvas aspect ratio, and a 40 px band reserved top and bottom for the overlay controls.
3. The snapshot is applied with `Camera.SetSnapshot` after two animation frames (so pending loader resets have resolved) and verified once; automatic camera resets stay disabled (`manualReset`).
4. The view direction and up vector are constants of the reference frame, so instances in the common orientation differ only in position/distance; highlighting never touches the camera; Reset orientation / Fit re-runs the same fit; Focus selection remains a separate action.

Measured on the rendered canvas (non-white pixels, overlays excluded, 1400 px viewport): Family Overview 53 % × 78 %, Structures 71 % × 77 %, Annotations 75 % × 66 %; applied camera distance equals the computed one (e.g. 69.0 Å); no atom reaches the overlay bands.

## Tests

- vitest: 12 passed (translation, orientation convention, species aggregation, SectionCard).
- Playwright: 20 passed, including the new `overview-layout.spec.ts`: desktop row order, six cards on one row, 40/60 top-aligned columns, chip and heading sizes, representative metadata and tooltip, no horizontal overflow; mobile stacking and label visibility at 390 px; camera occupancy 60–95 % of the limiting dimension on 4FE5, 1Y26 and 1Y27 with identical view direction, selection leaving the camera untouched, reset after a user rotation restoring the exact view. Existing suites (synchronised selection, orientation, structures selection, table) unchanged and passing.
- Python tests untouched (no backend change).

Screenshots: `docs/screenshots/m2_1_family_overview.png`, `docs/screenshots/m2_1_structures.png`, `docs/screenshots/m2_1_annotations.png`, `docs/screenshots/m2_2_orientation_*.png`.
