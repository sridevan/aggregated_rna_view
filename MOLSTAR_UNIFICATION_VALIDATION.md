# One 3D engine: Mol* in the Structures tab and the Interactions comparison

Date: 2026-10-19. Scope: the 3D viewers of the Structures tab (`AlignmentViewer`) and of the Interactions tab's Compare two
(`InteractionPairViewer`). Until now both used 3Dmol.js 2.5.2 (CDN), while the Family Overview used PDBe Mol*. Both now run on the
same PDBe Mol* (pdbe-molstar 3.x / molstar 5.8) engine as the Family Overview. Data, transforms, colours, fits and every control are
unchanged.

## 1. Design

- **Shared scene: `frontend/src/components/molstar/MolstarSuperposition.tsx`.**
  - One `PDBeMolstarPlugin` per mount, rendered without a molecule. It has a white background and hides the Mol* controls, its canvas
    buttons and its highlight box; the viewer keeps its own Fullscreen, Centre view and residue label.
  - **Per chain:**
    - one BinaryCIF fetch from the PDBe ModelServer (`label_asym_id`, model 1), cached per chain for the page session;
    - parsed in Mol*, then placed by a `TransformStructureConformation` node whose matrix is the `orientation.json` transform
      (x′ = R x + t, the family frame), optionally preceded by an extra rigid transform;
    - drawn as a uniform-colour cartoon.
  - **Updates without refetching:** changing a transform, the emphasis or the focused residues updates state nodes only. No chain is
    refetched or re-parsed.
  - **Camera:** the same convention as the Family Overview, a fixed view direction in the family frame fitted to the displayed chains
    or to the focused residues.
  - Mol*-free helpers (`residueLabel`, ModelServer URL, fetch cache with retry) live in `frontend/src/lib/molstarScene.ts`, so vitest
    can test them without loading Mol*. The scene is lazy-loaded.
- **Structures tab:** `AlignmentViewer` is a thin wrapper. It keeps the one-or-two-structure view, group colours, emphasis, the
  representative comparison and the `window.__rnaExplorer.alignment` / `alignmentProbe` test hooks.
- **Interactions Compare two:** `InteractionPairViewer` keeps the base-anchored local fit, now applied as B's extra matrix.
  - The fit is the Kabsch fit on the base heavy atoms of the interaction's first nucleotide, read in the family frame.
  - Changing the selected interaction only updates that matrix.
  - The selected nucleotides are drawn as ball-and-stick in each structure's colour over faded cartoons.
- **Removed:** the 3Dmol.js CDN script, the CIF text parser and the PDB writer that 3Dmol needed. 3Dmol could not read ModelServer
  CIF directly.

## 2. Problems found and fixed during the change

- **StrictMode double mount:** the second mount called `createRoot()` on the container of the torn-down plugin. Each plugin now
  renders into its own child element, skips rendering if already unmounted, and is disposed on unmount.
- **Canvas size:** Mol* only follows window resizes, so the canvas did not grow with the panel (drill-down). A `ResizeObserver` on
  the container now calls `handleResize`.
- **Status:** for one frame, `data-structures` named the new chain while `data-status` still said `ready`. The status is now derived
  from the entries the scene was last built for.
- **ModelServer:** it occasionally returns HTTP 502. Fetches retry twice on 5xx and network errors, with an 0.8 s and 1.6 s backoff;
  4xx is not retried.
- **Hover probe (test hook):** it now projects the residue's O3′ atom, Mol*'s nucleic-acid cartoon trace atom, so the probe point lies
  on the tube.
  - The label under the pointer can belong to a residue in front of the probed one; the tests check the label's format and chain, not
    the exact residue.

## 3. Results

Load times on a local dev server (Playwright, Chromium):

| Case | Ready after |
|---|---|
| RF00167, one structure (first mount) | 0.24–0.55 s |
| RF00167, add the second representative | 0.09–0.12 s |
| RF00177 16S, one structure | 1.5–2.0 s |
| RF00177 16S, add the second representative | 0.6–1.5 s |
| Interactions Compare two (5E54 B vs 5SWE X) | 2.8 s |

- **Local fit:** on 5E54 B vs 5SWE X, base pair 26–56, the fit gives RMSD 0.014 Å on 8 base atoms. The C1′–C1′ distances are the same
  before and after the fit, and each pair is fetched once.
- **Hover:** labels match the format `5E54.B A52 (author)`. Multi-letter author chains are shown in full (`4V9D.BA …`).
- **Screenshots** (`docs/screenshots/`):
  - `m25_structures_rf00167_one.png`
  - `m25_structures_rf00167_compare.png`
  - `m25_structures_hover.png`
  - `m25_structures_rf00177_compare.png`
  - `m25_interactions_compare_chains.png`
  - `m25_interactions_compare_local_fit.png`

## 4. Tests

- **vitest:** 99 passed. `AlignmentViewer.test.ts` keeps the `residueLabel` test. The 3Dmol-only parser and writer tests in
  `structureGroups.test.ts` were removed with that code.
- **Playwright:** full suite 184 passed (21.9 min) after the fixes in §2.
  - Re-run after the Mol* UI hiding: `structures-hover`, `structures-clustering` and `interactions` 38/38.
  - `structures-hover.spec.ts` now also checks that Mol*'s highlight box and canvas buttons stay hidden.
  - `structures-clustering.spec.ts`: the canvas must match the viewer box within 5 px, because of the 1 px border plus subpixel
    rounding (previously 3 px).
- **Home page:** the credit now reads "3D: PDBe Mol*". `home.spec.ts` still checks that no viewer script or coordinates are loaded there.

## 5. Follow-up after review: toolbars, Show ligands, ligand style

- **Toolbars:** in every Mol* viewer the controls moved out of the canvas into a toolbar strip above it. In the Family Overview viewer
  "Reset orientation / Fit" became "Reset view", and the shown ligand names moved from the label to the tooltip.
- **Family Overview viewer, layout side effects:** the "Loading …" note moved from the toolbar onto the canvas, so the toolbar height
  (and therefore the canvas size the camera is fitted to) never changes while loading. The viewer now follows its canvas area with a
  `ResizeObserver`. The Structures hover test scrolls the viewer into view before probing, because the toolbar moved it 35 px lower.
- **Ligand style (all viewers):** sticks with equal atom and bond radii (`sizeAspectRatio` 1). Carbons are in the RNA chain's or
  structure's colour; the other atoms are CPK.
- **Show ligands in the superposition viewer** (Structures tab, Interactions Compare two):
  - Off by default; nothing is fetched until it is ticked.
  - For each displayed chain, the ModelServer `residueSurroundings` query (label_asym_id, 4 Å, model 1, BinaryCIF) is fetched once and
    placed with the chain's own matrix, so the ligands follow the orientation and the local base fit.
  - In it Mol* selects small-molecule residues with any atom within 4 Å of the chain. Ions, waters, saccharides, lipids and the family's
    `additive` / `metal_ion` CCDs are left out, exactly as in the Family Overview.
  - Checked results: 6UC8 B ANG; 5E54 B none (apo); 9C6K X TPP; 8F4O A none (apo); 5SWE X ADE; 4V9D AA and 7OTC a none. 16S ligands load
    in about 1.8 s.
  - Test: `tests/structures-ligands.spec.ts`.
- **Tests after the follow-up:** full Playwright suite 186 passed and 1 failed. The failure was `refinements.spec.ts`, which read the
  background of the viewer's first child, now the toolbar; it now reads the canvas host and passes (3/3). vitest 99 passed.

## 6. Limitations

- The local fit and the whole-chain superposition are unchanged: visual comparisons in the family frame, not the base-centre fit
  behind the geometric discrepancy.
- The Complete-structure layers exist only in the Family Overview viewer; the superposition viewer shows chains and, on request, their
  nearby ligands. The two viewers share the engine and conventions, not one component.
- The static home-page thumbnails were rendered offline with 3Dmol.js. `thumbnails/manifest.json` keeps that as provenance.

## Structures tab: thin backbone trace and a contrasting pair of colours (2026-10-22)

**Scope** (user decision, revising the first version of the same day):
- Only the **Structures tab** 3D viewer draws polymers as a thin round backbone trace, like the PDBe-KB superposition viewer.
- The Family Overview viewer and the Interactions comparison are restored to the Mol* cartoon used before.
- `MolstarSuperposition` takes `polymerStyle`: 'cartoon' by default, 'trace' from `AlignmentViewer`, exposed as `data-polymer-style`.

**Trace settings** (`frontend/src/lib/molstarStyle.ts`, Mol* cartoon):
- visuals: `polymer-trace` and `polymer-gap` only;
- `aspectRatio` 1 and `nucleicProfile` elliptical;
- `sizeFactor` 0.35 (0.5 emphasis).

**How the trace was chosen:**
- Candidates were rendered on 4FE5 B + 1Y27 X (RF00167) and 16S 4V9D Y + 4U27 A: current cartoon, trace-only cartoon, uniform putty and backbone cylinders.
- Trace-only cartoon and uniform putty look the same; backbone cylinders look angular.
- At size 0.2 the 16S trace becomes hairlines.
- Images: `docs/screenshots/m20_rna_style_candidates.png`, `m20_rna_trace_thickness.png`.

### Colours of a superposed pair

**The rule** (`viewerColours`, `lib/structureGroups.ts`):
- Superposition in the Structures tab is always a pair, so the pair takes a fixed contrasting pair by display order: **violet `#a569ff`** and **dark amber `#bf8400`**.
- This applies to two group representatives, two members of one group, and a pair picked in the dissimilarity heatmap.
- A single structure keeps its group colour (representative) or group shade (member).
- The "Comparing" swatches under the viewer use the same colours.

**Why these colours.** The pair must not read as a structural group or as a dissimilarity value. It was chosen by an sRGB search for colours away from every colour already used in the tab:
- the group colours `#E8685A`, `#B0327A`, `#8C4A1E`, `#3FA7D6`;
- unplaced `#6b7480`;
- the heatmap selection blue `#1d4ed8`;
- the viridis ramp of the dissimilarity heatmap.

**Distances** (OKLab ×100; colour-blindness simulated with Machado 2009):

| | To the reserved colours | Under deutan / protan / tritan simulation |
|---|---|---|
| Violet | ≥ 20 | ≥ 10 |
| Amber | ≥ 12.6 | 5.0 to the coral group colour under deutan |
| The pair itself | 33 | ≥ 12.6 |

- The amber/coral case is covered by the swatch and label beside each structure.
- Both colours have ≥ 3:1 contrast on white, and all dataviz palette-validator checks pass.

**Rejected pairs:**
- vermilion/blue: ΔE 7.5 to the coral group, 5.9 to viridis;
- black/amber: black fails the lightness and chroma checks.

**Images:** `docs/screenshots/m20_structures_pair_compare.png` (group 2 vs group 1 with group chips and dendrogram) and `m20_structures_pair_heatmap.png` (a pair picked in the heatmap).

**Tests:**
- vitest: `viewerColours`.
- Playwright (86 passed: Structures, Interactions, Family Overview and family specs):
  - pairs take `#a569ff,#bf8400` and the trace style;
  - a single structure keeps its group colour;
  - the Interactions comparison stays `cartoon`.
