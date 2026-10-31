# R2DT 2D layouts (2d-3d-viewer branch): coverage on six families

Exploration on branch `explore/r2dt-2d-3d-viewer`, 2026-10-30. The R2DT developer suggested the new 2D/3D viewer
(https://github.com/r2dt-bio/R2DT/tree/2d-3d-viewer, commit `191a13d9`). Question: should our 2D panel use it, and for which families?

## What the R2DT 2D/3D viewer is

- `r2dt.py pdb_2d_3d <PDB> --chain X` downloads the structure, runs FR3D, draws the chain on a matched template (Rfam, GtRNAdb, RNase P,
  CRW/RiboVision) or template-free, and writes a `viewer/` folder: `api.json` + `fr3d.json` (the **same input contract** as the PDBe static
  topology files our `PdbRnaViewer2D` already consumes, label_seq_id based), `lbn.json` (layered dot-bracket), the structure file, the
  PDBe RNA Viewer 0.3.0 plugin and a 4.3k-line glue script `r2dt-2d-3d-viewer.js` that drives its own pdbe-molstar 3.12 pane.
- The plugin is the one we vendor. Ours is a newer upstream build (DOM look-ups scoped to `svg.rnaTopoSvg`); R2DT's carries one
  local patch: a filled circle for non-canonical cWW pairs (e.g. U-U) instead of the canonical line.
- The glue is one global instance per page (`R2DTViewer.create`, "v1 limitation: only one per page"), patches the plugin's DOM, installs a
  fetch shim and owns its Mol* pane. It does not fit our shared selection and shared Mol* scene, so it was **not** adopted; only the
  layouts are used.
- Rfam templates exist for ~3,800 families, but R2DT blacklists rRNA, tRNA, RNase P and tmRNA families from the Rfam search and draws them
  from dedicated libraries (GtRNAdb per isotype, RNase P, CRW/RiboVision). `pdb_2d_3d` has no `--force_template` (the `draw` command has).

## Method

- `run_r2dt.py <ACC> --r2dt-src <checkout>`: per chain in `instances.json`, `pdb_2d_3d <PDB> --chain <auth> --format cif --mode auto` in
  `rnacentral/r2dt:develop` (sha256:098154a4…, arm64) with the branch's `r2dt.py`, `utils/`, `data/viewer/` mounted over it (the image
  predates `pdb_2d_3d`). Raw outputs in `work/` (git-ignored, 4.3 GB of mmCIFs; deletable).
- **Acceptance gate** (the project rule that positions come only from validated mappings): every drawn position must equal our mapping
  residue at the same `label_seq_id`: parent base, author number (when observed), observed/unobserved flag. Modified residues drawn as `N`
  (`X` when unresolved) are allowed. Accepted chains are copied to `data/families/<ACC>/r2dt/<seq_name>/{api,fr3d,lbn}.json`;
  `r2dt/index.json` lists **every** chain with status, template, validation details and provenance (R2DT commit, image digest, mode, date).
- `analyse.py` → `coverage.json`: counts, templates, rejection causes, and layout consistency (distance between the drawn points of the same
  Rfam match position in every pair of chains on one template, after removing a translation).
- `compare_pdbe.py` → `pdbe_vs_r2dt.json`: for accepted chains, the PDBe static topology (`/pdbe/static/entry/<pdb>_<entity>_<AUTH>.json`,
  fetched live, cached in `work/pdbe/`) vs the R2DT layout, same label_seq_id after the best translation + uniform scale, in units of the
  median nucleotide spacing. < 0.5 spacing = same layout.

## Results (321 chains)

| family | chains | accepted | template | same layout as PDBe | different from PDBe | accepted, **no PDBe topology** | rejected |
|---|---|---|---|---|---|---|---|
| RF00167 purine | 84 | 70 | RF00167 68, template-free 2 (3DS7 A/B) | 58 | 0 | **12** | 14 |
| RF00162 SAM-I | 29 | 26 | RF00162 26 | 26 | 0 | 0 | 3 |
| RF00059 TPP | 38 | 26 | RF00059 26 | 0 | **26** | 0 | 12 |
| RF00005 tRNA | 60 | 39 | 25 GtRNAdb isotype templates (E_Phe 9, A_iMet 3, …) + RF00005 2 | 25 | 5 | **9** | 21 |
| RF00011 RNase P B | 50 | 50 | RNAseP_b_G_stearothermophilus_3D 50 | 24 | 0 | **26** | 0 |
| RF00177 16S | 60 | 29 | EC_SSU_3D 29 | 21 | 0 | **8** | 31 |
| **total** | **321** | **240** | | **154** | **31** | **55** | **81** |

- **PDBe already uses R2DT templates.** Where both exist, 154 of 185 chains have the same layout (residual 0.02–0.3 spacing; the small
  non-zero values are template revisions). The exceptions are all 26 TPP chains (PDBe draws a compact layout, not the RF00059 template)
  and 5 tRNAs.
- **The real gain is coverage: 55 accepted chains PDBe has no 2D for**, half of them RNase P (26 of 50; R2DT draws all 50 on one
  conventional *G. stearothermophilus* diagram).
- **Consistency within a template is excellent:** median pairwise offset ≈ 0 px for RF00162, RF00011, RF00177 and E_Phe tRNAs, 0.4 px TPP,
  5.5 px purine (spacing 30 px). One outlier: 9NBC A (160-nt scaffold construct, distorted RF00167 drawing, 185 px).
- **tRNA has no family-wide layout**: 25 isotype templates; chains align only within an isotype.
- TPP is a judgement call: the R2DT/Rfam template is the conventional orientation but sparse at panel size (screenshots `2hoj_pdbe.png`
  vs `2hoj_r2dt.png`); PDBe's compact drawing reads better at 560 px.
- Accepted with gaps: 7 RF00177 chains miss their last 1–8 3′ positions in 2D; tRNA 3′ CCA is not drawn where the template ends earlier.
  Those positions cannot be highlighted in 2D (3D unaffected). Modified residues drawn as `N`: 5 RF00167, 9 RF00005, 13 RF00177 chains.

## R2DT failure modes (all caught by the gate; candidates for feedback to R2DT)

1. **Warned but still written, positions shifted.** "Could not map to full sequence: Resolved count (n) != dot-bracket length (m)" and
   the viewer files are produced anyway (RF00167 10, RF00059 10: 9NBH A–D, 5SWD A, 9V4V/9V4Z/9V51 A, 9LMC B, 6P2H A).
2. **Alternate conformations counted as extra nucleotides, no warning.** 5FK3 A residue 9 altloc A/B → one inserted nucleotide; 9NDD A–D
   residues 141–144 → 4 inserted (164 drawn for 160). RF00167 4, RF00162 3, RF00005 3 (1FIR, 2TRA, 9SFQ; altloc residues checked in the mmCIF and matching the
   insertion positions).
3. **Unresolved residues dropped on the library-template path, no warning.** The full-sequence/mask expansion does not reach
   RiboVision/GtRNAdb drawings (9SRO A: 31 unresolved residues, 1,501 drawn as 1–1501, shifted after residue 80). All 31 RF00177 rejections
   and most RF00005 rejections are chains with internal unresolved residues.
4. **5′ residue outside the tRNA template dropped and the rest renumbered** (1QTQ, 5E6M, 7XJZ, 7XK0: label 1 is missing, every position
   shifts by one).
5. Cosmetic: modified residues drawn as `N`/`X`; unresolved residues get label numbers as author numbers (5E54 A 66–71 instead of 78–83);
   tRNA 3′ CCA not drawn.
6. Request: expose `--force_template` in `pdb_2d_3d` so a family's chains can all be drawn on its own template (3DS7 A/B match nothing).

## Built on the branch

- `experiments/r2dt_viewer/run_r2dt.py`, `analyse.py`, `compare_pdbe.py`; outputs `coverage.json`, `pdbe_vs_r2dt.json`, `screenshots/`.
- `data/families/<ACC>/r2dt/` for all six families (240 accepted chains, ~14 MB, RF00177 7.7 MB) + `index.json`.
- Frontend: `PdbRnaViewer2D` takes an optional `layoutUrl` (R2DT `api.json`/`fr3d.json` instead of the PDBe files; `data-layout`
  attribute); `StructureExplorer` shows a **2D layout: PDBe | R2DT** toggle when the family has `r2dt/index.json` (PDBe default; the note
  names the template, or says the chain has no accepted R2DT layout and PDBe is shown); `loadR2dtIndex` in `api/data.ts`.
- `frontend/tests/r2dt-layout.spec.ts` (5 tests, pass): toggle on 4FE5 B with selection mirrored by label_seq_id and PDBe restored;
  toggle persists across instances; TPP 2HOJ A, RNase P 9OV3 A (PDBe unavailable, R2DT ready), tRNA 1EHZ A. The plugin's
  `createBPDropdown` console error during re-render also occurs on `main` (instance switching).

## Recommendation

1. **Use R2DT as a fallback, not a replacement**: show the R2DT layout automatically when PDBe has no topology and the chain passed the
   gate (+55 chains, 26 of them RNase P); keep PDBe elsewhere, since it is the same template for 154/185 chains. Drop the toggle, or keep it
   only where layouts differ (TPP, 5 tRNAs).
2. **Never use R2DT output without the gate**: three of its failure modes are silent and shift positions.
3. Before a pipeline command (`rna_explorer r2dt`, schema, tests, provenance) wait for `pdb_2d_3d` in a tagged image, or pin commit +
   digest as here.
4. Send failure modes 1–4 and the `--force_template` request to the R2DT developer (step 5); fixes to 2 and 3 would recover most of the
   81 rejected chains.

## Reproduce

    git clone --branch 2d-3d-viewer https://github.com/r2dt-bio/R2DT.git <r2dt>   # commit 191a13d9
    docker pull rnacentral/r2dt:develop
    python3 experiments/r2dt_viewer/run_r2dt.py RF00167 --r2dt-src <r2dt>        # cached per chain in work/
    python3 experiments/r2dt_viewer/analyse.py RF00167 RF00162 RF00059 RF00005 RF00011 RF00177 > experiments/r2dt_viewer/coverage.json
    python3 experiments/r2dt_viewer/compare_pdbe.py RF00167 RF00162 RF00059 RF00005 RF00011 RF00177 > experiments/r2dt_viewer/pdbe_vs_r2dt.json
    cd frontend && npx playwright test tests/r2dt-layout.spec.ts
