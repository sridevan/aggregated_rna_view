# Isolated experiment: PDBe RNA Viewer (2D) + PDBe Mol* (3D) driven by Rfam consensus positions (PDB 1y26, RF00167)

Not part of the main frontend. Run from the repository root:

    python3 -m http.server 8765 --bind 127.0.0.1
    open http://127.0.0.1:8765/experiments/rna2d_poc/index.html

- `vendor/pdb-rna-viewer/` git clone of PDBeurope/pdb-rna-viewer 0.3.0 (Apache-2.0), built with
  `npx tsc --incremental --skipLibCheck && NODE_OPTIONS=--openssl-legacy-provider npx webpack --mode production && npx gulp`
  (the stock `npm run build` fails on Node 20: OpenSSL "digital envelope" error and a @types/trusted-types redeclaration).
- `node_modules/pdbe-molstar/` pdbe-molstar 3.12.0 (Apache-2.0) bundle `build/pdbe-molstar-plugin.js`.
- `index.html` loads our mapping `data/families/RF00167/mappings/1y26_1_X_A.json`, translates consensus -> label_seq_id, and wires
  `PDB.RNA.viewer.click` (2D -> consensus -> 3D) and `PDB.molstar.click` (3D -> consensus -> 2D).
- `data/` the PDBe static inputs the viewer fetches (`1y26_1_X.json` topology, `1y26_X_basepair.json` FR3D pairs), build logs, screenshot.
