# RNA visualisation components: review and recommendation (Milestone 2 input)

Date: 2026-10-14. Sources inspected directly: PDBeurope/pdb-rna-viewer (git clone, source and built bundle), pdbe-molstar 3.12.0 (npm package type definitions and source), R2DT documentation (docs.r2dt.bio usage page), RiboVision2 (Crossref/PubMed record, About page, GitHub repository listing). Isolated proof of concept in `experiments/rna2d_poc/` (PDB 1y26, chain X, RF00167).

## 1. Tool comparison

| | PDBe RNA Viewer (`pdb-rna-viewer`) | R2DT | RiboVision2 |
|---|---|---|---|
| What it is | Lit/D3 web component drawing the PDBe-precomputed 2D layout of one PDB RNA chain, with FR3D base pairs | Pipeline/web service that draws a 2D diagram for any RNA sequence from a template library (Rfam, CRW, RiboVision, GtRNAdb, RNase P) | Web server for rRNA: phylogenetic/structural data on 1D/2D/3D views of ribosomes |
| Scope | Per PDB chain (PDB-specific numbering) | Per sequence, family-level templates | rRNA only (SSU/LSU), ribosome complexes |
| Inputs | `pdbId`, `entityId`, `chainId`; optional pre-fetched `apiData`, `FR3DData` | FASTA (+ optional dot-bracket); CLI `r2dt.py draw`, `r2dt.py rfam draw RF00162 ...`, EBI Web Services API | PDB entries, user alignments (uploaded to server), CSV |
| Data source | `https://www.ebi.ac.uk/pdbe/static/entry/{pdb}_{entity}_{CHAIN}.json` (svg paths, sequence, `label_seq_ids`, `auth_seq_ids`, `pdb_ins_codes`, `unobserved_label_seq_ids`) and `{pdb}_{CHAIN}_basepair.json` (FR3D: `seq_id1/2` = label_seq_id, `3d_id1/2` = auth, `bp` class, `crossing`) | Templates + Infernal/R-scape/RNArtist layouts; outputs `svg/`, `fasta/` (dot-bracket), `json/` (RNA 2D JSON schema), `tsv/metadata.tsv` (template used) | Server-side; repository bundles `pdbe-rna-viewer` and `pdb-topology-viewer` folders |
| Residue numbering | `label_seq_id` for all public methods and emitted events; author numbers only displayed | Ordinal positions of the input sequence (`svg2json.py` "ordinal numbering"); no documented link to Rfam CM consensus positions | not documented in accessible sources |
| Programmatic highlight | `selectResidue(label_seq_id)`, `selectResidueRange(a,b)`, `highlightResidue`, `highlightResidueRange`, `clearSelection`, `clearHighlight` (all with `doNotReport` to suppress echo events) | none (static SVG/JSON) | n/a |
| Events | emits `PDB.RNA.viewer.click` / `mouseover` / `mouseout` / `deselect` with `eventData = {pdbId, entityId, label_seq_ids}`; with `subscribeEvents: true` listens to `protvista-click/mouseover/mouseout` (`detail.start/end` = label_seq_id range) and `PDB.molstar.click/mouseover/mouseout` (`eventData.residueNumber`, filtered by `auth_asym_id === chainId`) | n/a | n/a |
| Secondary structure annotations | FR3D base pairs (cWW, tHS, ... with `crossing` flag for pseudoknots), rendered as dots/lines; "Base Pairings List" and nested-only filter in the UI | predicted/template base pairs in dot-bracket | n/a |
| Install | not on npm (`pdb-rna-viewer` 404); build from git. Stock build fails on Node 20 (OpenSSL + TS types); works with `--skipLibCheck` and `NODE_OPTIONS=--openssl-legacy-provider` | Docker/Singularity/bare metal; or EBI API | hosted; MIT per About page |
| License | Apache-2.0 (EMBL-EBI) | Apache-2.0 | MIT (web server); bundled viewer licences not stated |
| Published description | PDBe component (RiboVision2 repository embeds it) | Sweeney et al. Nat Commun 2021; NAR 2025 | McCann, Meade, Banerjee, Penev, Williams, Petrov, "RiboVision2: A Web Server for Advanced Visualization of Ribosomal RNAs", J Mol Biol 436(17):168556, 2024, PMID 39237196 |

Journal article note: the ScienceDirect page (PII S0022283624001517) returned HTTP 403, so only the Crossref record, the PubMed abstract and the server's About page were used. The abstract states that RiboVision2 "visualize[s] phylogenetic, structural, and evolutionary properties of ribosomal RNAs simultaneously at the levels of primary, secondary, and three-dimensional structure"; the detailed methods (layout generation, alignment mapping, 3D viewer) are not in any accessible source and are therefore not described here. The repository file listing shows it embeds the PDBe RNA viewer and a Python backend; the 3D viewer is not identified.

## 2. PDBe RNA Viewer: findings relevant to integration

- **Numbering matches our mapping exactly.** The topology JSON for 1y26 chain X lists `label_seq_ids` 1-71 and `auth_seq_ids` 13-83, identical to `mappings/1y26_1_X_A.json` (polymer positions 1-71, author 13-83). The viewer's whole API is label_seq_id based, so `consensus -> polymer_seq_position` from our mapping is the only translation needed.
- **PDBe Mol* (`pdbe-molstar` 3.12) uses the same identifiers.** `visual.select({data: [{label_asym_id, label_seq_id, color, focus}]})` is the current (non-deprecated) form; `struct_asym_id`/`start_residue_number` are deprecated aliases. Its click/hover events carry `seq_id` (= label_seq_id), `auth_seq_id`, `auth_asym_id`, `entity_id`, `ins_code`, and `residueNumber` (copy of `seq_id`).
- **Compatibility caveat between the two PDBe components:** the RNA viewer filters Mol* events by `auth_asym_id === chainId` while selecting by label_seq_id. This is consistent for 1y26 (auth X, label A) because the filter uses the author chain, but any integration must pass the author chain as `chainId` and the label asym to Mol*. Our mapping files carry both.
- **Unobserved residues:** the topology JSON lists `unobserved_label_seq_ids` and the viewer draws them in a distinct colour and suppresses their events; our `has_coordinates` flag gives the same information and must be the gate for 3D highlighting.
- **Insert / local-end residues:** the 2D layout is per PDB chain, so every residue has a 2D position regardless of CM state. Residues that have no consensus position (insert, local end) can still be shown and clicked in 2D; the reverse mapping returns `state != match` and the consensus selection should not change (the POC ignores them).
- **Modified residues:** the viewer draws the one-letter sequence from PDBe; no special handling needed.
- **Global coupling:** events are dispatched on `document` and the viewer keeps module-level static state (`UiActionsService.selected/highlighted`, DOM ids like `${pdbId}-rnaTopologyTooltip`). Two instances of the same PDB id on one page would collide; one viewer per page is safe. React integration must create the plugin in an effect, call `render(container, options)` once, and remove listeners on unmount.
- **Build:** no npm release; vendor the built bundle (plugin JS + CSS, ~95 kB + 5 kB) or build from the fork in CI with the two Node 20 workarounds. Apache-2.0 permits redistribution with attribution.

## 3. Proof of concept (PDB 1y26, RF00167), `experiments/rna2d_poc/index.html`

All five steps requested were exercised in Chrome against the local static server (screenshot: `experiments/rna2d_poc/data/poc_1y26_sync_screenshot.jpg`):

1. 2D display: PDBe RNA viewer rendered 71 nucleotides with FR3D pairs for 1y26 chain X.
2. 3D display: PDBe Mol* loaded 1y26.
3. Consensus selection: consensus 50 -> mapping -> polymer 34 (author G46), state `match`, has coordinates.
4. Highlight in both viewers: `rna.selectResidue(34)` coloured G34 in 2D; `mol.visual.select({label_asym_id:'A', label_seq_id:34, focus:true})` coloured and focused the residue in 3D (screenshot).
5. Bidirectional sync with real events: a DOM click on the 2D nucleotide for label_seq_id 20 emitted `PDB.RNA.viewer.click` with `label_seq_ids [20]`, reverse-mapped to consensus 36 and selected in 3D; a real mouse click on the Mol* canvas emitted `PDB.molstar.click` with `seq_id 34`, `auth_asym_id X`, reverse-mapped to consensus 50 and selected in 2D. A synthetic Mol* event for residue 10 behaved the same. Deletion (consensus 64) and uncovered (consensus 1) positions were correctly reported as not highlightable.

Observed quirk: after programmatic `clearSelection(undefined, true)` + `selectResidue`, previously selected 2D text elements kept a `fill` attribute; the viewer's visible selection state was correct in screenshots, but a React wrapper should track its own selected set and call `clearSelection` per residue rather than rely on the viewer's internal map.

## 4. R2DT as a family-level layout generator

R2DT draws one sequence at a time against a template; for an Rfam family it uses the Rfam template (`r2dt.py rfam draw RF00167 ...`, choosing the layout engine with the fewest overlaps). Its JSON gives per-nucleotide coordinates with ordinal numbering of the input sequence; it does not emit Rfam CM consensus positions. A family-level 2D would therefore be produced by (a) drawing the Rfam consensus (RF) sequence itself as the input, which yields one coordinate per consensus position by construction, or (b) drawing each instance and re-mapping ordinals through our mapping. Option (a) is the only one that gives a true family coordinate system; insertions have no slot in it, which matches our numbering policy. This is a precomputation step (Docker/EBI API), not a frontend component, and is deferred.

## 5. Recommendation

- **Structures tab:** keep Mol* for 3D. Switching the Milestone 2 viewer from raw Mol* (`createPluginUI`) to `pdbe-molstar` is recommended: it gives the label-id selection API, the standard PDBe look, and the click events the RNA viewer already understands, with the same Apache-2.0 licence. This replaces a component that does not exist yet (the MolstarViewer wrapper was not written), so no existing component is replaced.
- **Annotations tab:** add the PDBe RNA viewer as an optional "2D (PDB chain)" panel beside the consensus tracks and Mol*, driven by the same `SelectionContext` (consensus positions -> `translateConsensus` -> label_seq_id). It is PDB-specific, so it changes when the instance changes while the consensus selection persists, exactly like the 3D panel. It also delivers FR3D base-pair annotations for free, which Milestone 3 can read from the same `{pdb}_{chain}_basepair.json`.
- **Custom 2D viewer:** not necessary for PDB-chain 2D. A custom (or R2DT-derived) family-level 2D of the consensus is still needed later for the R2DT-style view planned in the project brief, because no existing component draws a layout in Rfam consensus coordinates.

## 6. Proposed implementation plan (after approval)

1. Vendor `pdb-rna-viewer` 0.3.0 build output into `frontend/public/vendor/` (JS + CSS, with LICENSE) or add the git repo as a build-time dependency with the two Node 20 fixes; add `pdbe-molstar` from npm.
2. `PdbeMolstarViewer.tsx`: wraps `PDBeMolstarPlugin.render` in an effect; props `{pdbId, labelAsymId, selection: {label_seq_id[]}}`; uses `visual.select` with `label_asym_id` + `label_seq_id`; subscribes to `PDB.molstar.click` on its container and reports `label_seq_id` upward (not to `document`).
3. `PdbRnaViewer2D.tsx`: wraps `PdbRnaViewerPlugin.render` with `subscribeEvents: false` (so our state, not the viewer, is the source of truth); applies the current translation with `selectResidueRange`/`selectResidue(..., doNotReport=true)`; listens to `PDB.RNA.viewer.click` and reports `label_seq_ids` upward; one instance mounted at a time.
4. `AnnotationsTab`: both panels consume `translateConsensus(mapping, selection.positions).highlightable`; clicks from either viewer go through `residueToConsensus` (new helper in `lib/consensus.ts`, inverse of the mapping) and only `match` residues update `SelectionContext`.
5. Tests: unit tests for the inverse mapping; a Playwright/vitest-browser smoke test is optional since the POC page already demonstrates the integration manually.
6. Later (Milestone 3+): FR3D pair track from the PDBe basepair JSON; family-level 2D via R2DT drawing of the consensus sequence.
