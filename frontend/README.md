# Aggregated RNA 3D View (RNA Structure Explorer): frontend

React 19 + TypeScript + Vite 8 + Tailwind CSS 4, PDBe Mol* (`pdbe-molstar` 3.12, npm) for 3D and the PDBe RNA Viewer 0.3.0 (vendored build in `public/vendor/pdb-rna-viewer`, Apache-2.0, see `VENDORED.txt`) for PDB-chain 2D.
Reads only the static JSON produced by the Python pipeline (`data/families/<ACC>/`, exposed as `/families/<ACC>/...` through the `public/families` symlink). No Python, Infernal or API server at runtime; coordinates and 2D topology are fetched from PDBe.

## Setup and run

    cd frontend
    npm install                      # Node 20.14 works (Vite prints a version warning); if the build complains about a missing rolldown binding:
                                     #   npm install --no-save @rolldown/binding-darwin-arm64@$(node -p "require('rolldown/package.json').version")
    npm run dev                      # http://127.0.0.1:5173
    npm run build                    # dist/ incl. a copy of data/families (npm run copy-data)
    npm run preview                  # serve dist/ on :4173
    npm test                         # vitest: consensus translation unit tests (use data/families/RF00167)
    npx playwright install chromium  # once
    npm run test:e2e                 # Playwright: 7 browser tests (starts the dev server; needs network access to PDBe)

## Layout

| Path | Content |
|---|---|
| `src/App.tsx` | providers, banner (app title; `family.description` as heading; `family.rfam_acc` badge linking to `family.rfam_url`) + three tabs |
| `src/state/selection.tsx` | authoritative selection: consensus positions (family coordinates, preserved across instances) + instance-only residues |
| `src/state/family.tsx`, `src/state/useMapping.ts` | start-up JSON, lazy per-chain mapping |
| `src/lib/consensus.ts` | consensus -> residues (`translateConsensus`), residue -> consensus (`residueToConsensus`), per-position states |
| `src/components/ConsensusTracks.tsx`, `TrackLegend.tsx` | axis, sequence coverage, structural coverage, PP confidence, instance state track; legend |
| `src/components/InstanceTable.tsx`, `InfoTip.tsx` | paginated/searchable/sortable table (opens on the page of the selected instance, flags a selection hidden by the search), accessible metric tooltips |
| `src/components/PdbeMolstarViewer.tsx` | 3D: one plugin per mount, flash-free structure switching, precomputed orientation transform (orientation.json) + fixed camera direction, Rfam-region / complete modes, white background, hidden controls, label-id selection (no auto-focus), Focus / Reset orientation actions |
| `src/lib/orientation.ts` | R,t -> Mol* Mat4 (convention test in `orientation.test.ts`) |
| `src/components/PdbRnaViewer2D.tsx` | 2D: author chain for data, label_seq_id selection, click -> label_seq_id; one instance at a time |
| `src/components/AnnotationsTab.tsx` | wires selection, translation and both viewers |
| `tests/*.spec.ts` (17) | browser tests (see MILESTONE2_VALIDATION.md, MILESTONE2_1_VALIDATION.md); `screenshots.spec.ts` regenerates docs/screenshots |

Window hook `window.__rnaExplorer` (selection, translation, last Mol* selection) exists for tests and debugging only.
