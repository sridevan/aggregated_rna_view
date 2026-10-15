# Milestone 2.2 validation report: representative structure, flash-free viewer, common 3D orientation

Date: 2026-10-15. Residue mapping, Infernal alignment and residue-correspondence code untouched (10 → 13 pytest with the new orientation tests; 3 → 6 vitest).

## 1. Single representative

The Family Overview shows one card, "Representative structure": **4FE5 chain B** (*Bacillus subtilis* guanine riboswitch, X-ray 1.32 Å, 67/67 matched residues with coordinates, PDBe 2D topology available), chosen by the existing display-representative algorithm (`rna_explorer/instances.py::choose_display_representative`, rule printed on the card). The coverage-based `representative` field stays in `family.json` (used as a fallback when no instance satisfies the display rule) but is no longer rendered. The rule is family-agnostic.

## 2. Mol* background transition

Cause: `visual.update(options, fullLoad = true)` calls `plugin.clear()`, which re-initialises the canvas with Mol*'s default (dark) renderer settings; the white `bgColor` was only re-applied after the new structure loaded, so a black frame was visible between the two.
Fix (`frontend/src/components/PdbeMolstarViewer.tsx`):
- One plugin per mounted viewer; structures are switched with `visual.update({ moleculeId }, false)` (loads the new entry without clearing) followed by `hierarchy.remove(<previous trajectories>)`, so the canvas and its settings are never re-created.
- `bgColor` white at `render()`, re-asserted with `canvas.setBgColor` after render and after every switch (setting camera props re-derives renderer defaults).
- A white overlay covers the canvas until the new structure, orientation and view are applied (`data-status="loading"`), so no unconfigured canvas is visible.
- Structure Tools panel hidden (`hideControls`, no reserved space); synchronised highlighting unchanged.
Verified by the browser test "switching structures keeps a white canvas…" over three consecutive switches (4FE5 → 1Y26 → 1Y27 → 3RKF): overlay present and white while loading, renderer background `#ffffff` and exactly one structure in the hierarchy after each switch.

## 3. Orientation dataset (`data/families/RF00167/orientation.json`, schema 1.0)

Produced by `python3 -m rna_explorer orient RF00167` (`rna_explorer/orientation.py`), regenerable independently of the residue mapping.

- Reference: 4FE5 chain B (label asym A). Targets: all 84 instances.
- Correspondence: consensus positions that are match states in both chains, have coordinates in both, and PP `*` (≥ 0.95) in both (relaxed to ≥ `9` only when fewer than 10 pairs remain; never needed except for the 9NBH chains, which still have 0 pairs). Atom per pair: C4′, falling back to P then C1′ only if C4′ is absent in either residue (all 80 fits used C4′ exclusively).
- Coordinates: first model and first alternate location of the PDBe updated mmCIF (`mmcif.atom_site_coords`).
- Fit: Kabsch (SVD) rotation + translation, no scaling; reflections rejected (det R = +1 in all 80 transforms).
- Stored per target: rotation (row-major 3×3), translation, RMSD, number of pairs, shared consensus counts, reference/target coverage (pairs ÷ matched residues with coordinates), PP filter used, atoms used, excluded counts (no coordinates in target/reference, low PP, missing atom), paired consensus positions, determinant, fit quality (good ≤ 2 Å and ≥ 20 pairs; moderate ≤ 4 Å; poor; insufficient < 6 pairs → no transform), and a `check` point (one C4′ before/after) used by the frontend convention test.
- **Convention:** `x_ref = R · x_target + t` (column vectors). Mol*: `Mat4.ofRows([[R00,R01,R02,t0],[R10,R11,R12,t1],[R20,R21,R22,t2],[0,0,0,1]])` applied with `TransformStructureConformation` (`transpose: false`), inserted below the structure node exactly as Mol*'s own superposition UI does. The target structure is moved into the reference frame; the reference transform is the identity. Verified numerically: the vitest `orientation.test.ts` rebuilds the Mat4 for all 80 applicable transforms and reproduces the backend `check` point to 1e-6 Å with `Vec3.transformMat4`.

Results (RF00167, 84 instances): 51 good, 27 moderate, 2 poor (5E54 chain B 4.63 Å, 5SWD chain B 4.59 Å: genuinely different conformations, reported not forced), 4 insufficient (9NBH chains A–D: the aptamer core has no coordinates, 0 shared positions with coordinates). RMSD median 1.82 Å.

| Target | pairs | RMSD (Å) | quality | note |
|---|---|---|---|---|
| 4FE5/B | 67 | 0.000 | good | identity (`is_identity: true`, rotation = I, translation = 0) |
| 1Y26/X | 67 | 1.898 | good | adenine riboswitch; per-pair deviation median 1.46 Å, max 4.16 Å at consensus 19 (chain end) |
| 1Y27/X | 67 | 1.444 | good | |
| 3RKF/A | 67 | 1.002 | good | |
| 6P2H/A | 53 | 2.768 | moderate | partial model coverage (consensus 26–80): reference coverage 0.79 |
| 9NDD/A | 61 | 2.315 | moderate | 4 pairs excluded for low PP |
| 9NBH/A | 0 | – | insufficient | core unmodelled; shown as deposited |

Reproducibility: recomputing the 1Y26 fit from raw coordinates outside the pipeline gives the identical rotation and RMSD 1.8984. Chain identity: coordinates are read by label_asym_id + entity (1Y26 label A = author X; 71/71 residues have C4′).

## 4. Mol* integration

- The transform is inserted once per loaded structure below the structure node (updated in place if present); the region component is created below the transform node so both display modes are transformed; pdbe-molstar's label-id selection works through Mol*'s hierarchy (`structure.transform`) unchanged.
- Camera: Mol*'s automatic camera reset is disabled (`camera.manualReset`), and after every load/mode change the camera is set from an explicit snapshot (centre of the visible region, fixed view direction (0,0,−1) and up (0,1,0), distance from the bounding radius). Note: `Camera.getFocus(target, radius, up, dir)` only hemisphere-matches the current direction (`Vec3.matchDirection`), so it cannot give a fixed orientation; the explicit snapshot is required.
- Highlighting never moves the camera; "Focus selection" and "Reset orientation / Fit" are separate buttons. A badge reports the superposition quality or "as deposited".
- Limitation found: pdbe-molstar re-shows its canvas buttons after `visual.update` (cosmetic), and logs one "reading 'transform'" rejection when the previous structure is removed (harmless, after load completion).

## 5. Tests

| Suite | Result |
|---|---|
| pytest | 13 passed: Kabsch identity, known rotation recovery and convention, reflection rejection, RMSD reproducibility, < 3 points error, pair selection with missing atoms / fallback atom / low PP relaxation / identifiers preserved, quality bands, display-representative rule |
| vitest | 6 passed: translation library; orientation convention (identity reference, 80 check points, insufficient fits not applied) |
| Playwright | 14 passed: Milestone 2/2.1 suites (updated for the single card) + `orientation.spec.ts`: single card, 3 consecutive switches with white overlay and white renderer, one structure loaded; transform applied for 1Y26 (RMSD 1.90, 67 pairs, badge), camera direction unchanged by selection and identical between 1Y26 and 1Y27, consensus 50 → label_seq_id 34 (1Y26) / 33 (1Y27) after the transform, reference identity badge, 9NBH shown as deposited, user rotation then Reset orientation restores the fixed direction; screenshots |

Screenshots (common orientation): `docs/screenshots/m2_2_orientation_{4fe5,1y26,1y27,3rkf,2xnw}.png`, `docs/screenshots/m2_2_family_overview.png`; the three-tab screenshots `m2_1_*.png` were regenerated.

## 6. Remaining limitations
- Orientation is relative to the display representative; a change of representative (or Rfam release) requires regenerating `orientation.json`.
- Partial-coverage instances are superposed on their shared positions only; their screen position depends on the centre of their own visible region.
- Insufficient fits (9NBH) are shown as deposited with a badge; no alternative reference is attempted.
- The extra Mol* canvas buttons after a structure switch and the pdbe-molstar internal rejection noted above are not fixed (library behaviour).
