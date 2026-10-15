# Aggregated RNA 3D View: project handover

Last updated: 2026-10-15 (end of Milestone 3B). Owner: asridevan86@gmail.com. Working directory `/Users/sria/Desktop/aggregated_rna_view` (macOS arm64, Python 3.11.9 via pyenv, Node 20.14 via nvm). **Not yet a git repository** (see §7).

Goal: a PDBe-KB-style aggregated view of experimentally determined RNA structures per Rfam family, using Rfam covariance-model (CM) match-state positions as the family coordinate system. RF00167 (purine riboswitch) is the development family; RF00177 (16S rRNA) was used in the proofs of concept.

## 1. Architecture

### Backend (`rna_explorer/`, Python, `python3 -m rna_explorer <cmd>`)
| Module | Role |
|---|---|
| `cache.py` | file cache for all HTTP downloads under `cache/<host>/` (hash-named); builds are offline once warm |
| `rfam.py`, `pdbe.py`, `mmcif.py`, `ncbi.py` | retrieval: Rfam family/CM/structures/clan; PDBe nucleic_mappings, molecules (POST, form content-type), residue listing, search API (`rfam_accession:` for organism), updated mmCIF; NCBI taxonomy (strain→species) |
| `sequences.py`, `infernal.py`, `stockholm.py` | chain instances from PDBe-mapped regions, unique-sequence dedup, `cmalign` (local mode) on uniques, Stockholm parsing (RF/SS_cons/PP) |
| `mapping.py`, `validation.py`, `confidence.py` | residue ↔ consensus mapping (states `match`/`insert`/`local_end`/`not_submitted_flank`; consensus columns `match`/`deletion`/`uncovered`; PP bins `'*'`=0.975, `'0'`=0.025, `'1'..'9'`=0.1..0.9; display bands high ≥ 0.95, medium ≥ 0.80), validation vs `_pdbx_poly_seq_scheme`, per-position statistics |
| `instances.py` | instance rows, `choose_representative` (coverage) and `choose_display_representative` (rule below) |
| `orientation.py` | Kabsch superposition on C4′ onto the display representative (`x_ref = R·x_target + t`, det = +1) |
| `export.py`, `schemas/*.schema.json` | JSON writers validated with jsonschema; CSV twins |
| `pipeline.py`, `cli.py` | `build`, `orient`, `annotate` commands |
| `annotation_sources.py`, `annotation_mapping.py`, `annotations.py`, `annotation_config.json` | Milestone 3B annotation backend (§3) |
| `tools/infernal-1.1.5/bin/cmalign` | Infernal built from source; found via `RNA_EXPLORER_CMALIGN`, then this path, then `PATH` |

Commands:

    pip install requests jsonschema pytest
    python3 -m rna_explorer build RF00167 --cpu 4 -v          # fetch → instances → cmalign → mapping → validation → stats → export
    python3 -m rna_explorer orient RF00167                     # orientation.json
    python3 -m rna_explorer annotate RF00167 --raw-dir annotation_assessment/cache   # annotation files; --offline forbids network
    python3 -m pytest -q                                       # 29 tests

### Data (`data/families/RF00167/`, served to the frontend via symlink `frontend/public/families → ../../data/families`)
`family.json` (description, comment, `rfam_url`, clen 102, consensus sequence, SS_cons, clan members, counts, `representative` = 9ndd/A, `display_representative` = 4fe5 chain B), `instances.json` (84 rows incl. organism, species, tax id, `pdbe_2d_available`), `consensus_stats.json`, `validation.json`, `orientation.json`, `mappings/<seq_name>.json` (seq_name = `<pdb>_<entity>_<auth_chain>_<label_chain>`), `residue_mapping.csv`, `consensus_columns.csv`, `work/` (alignment, ignored by git), and the 3B annotation files `reference_structural_elements.json`, `loops.json`, `ligands.json`, `tracks.json`, `annotation_summary.json`, `annotations/<seq_name>.json`.

### Frontend (`frontend/`, React 19 + TypeScript + Vite 8 + Tailwind 4)
- 3D: `pdbe-molstar` 3.12 (`pdbe-molstar/lib/viewer`, needs `sass-embedded`; rolldown binding installed with `--no-save @rolldown/binding-darwin-arm64`). `components/PdbeMolstarViewer.tsx`: one plugin per mount, structure switch via `visual.update(…, false)` + hierarchy removal, white background, orientation transform via `TransformStructureConformation`, region/complete modes, custom camera fit (fixed view axes, exact perspective distance, margin 1.06, 40 px bands, `Camera.SetSnapshot` after two frames), "Reset orientation / Fit" and "Focus selection", 60 s watchdog + Retry.
- 2D: vendored PDBe RNA Viewer 0.3.0 in `public/vendor/pdb-rna-viewer` (Apache-2.0), `components/PdbRnaViewer2D.tsx` (prefetches PDBe topology/FR3D, explicit unavailable state, SVG fit).
- State: `state/family.tsx`, `state/selection.tsx` (consensus positions preserved across instances; `instanceResidues` for insert/local-end residues); `lib/consensus.ts` translate/inverse; `lib/species.ts`; `lib/orientation.ts`.
- Tabs: Family Overview (About → six stat cards → molecular view of the display representative, 560 px → 40/60 top-aligned row: representative rule + methods | species chart → collapsible consensus reference), Structures (`InstanceTable.tsx`, 10/page, auto-page to selection), Annotations (`ConsensusTracks.tsx`, `TrackLegend.tsx`, synchronised 2D/3D).
- Tests: `npm test` (vitest, 12), `npm run test:e2e` (Playwright, 20 specs in `frontend/tests/`); screenshots in `docs/screenshots/`. `frontend/ARCHITECTURE.md` predates pdbe-molstar and still names `MolstarViewer.tsx`; treat this file and `frontend/README.md` as the current description.
- Known harmless console noise: pdbe-molstar "Invalid data cell" (headless), one "reading 'transform'" rejection on structure switch.

### Other directories
`poc_rf00167/`, `poc_rf00177/` (+ `bgsu_validation/`): original proofs of concept and reports. `annotation_assessment/`: Milestone 3A assessment, CSV samples, `cache/{bgsu,pdbe}` raw responses, `scripts/assess.py`. `experiments/rna2d_poc/`: isolated 2D/3D experiment. `docs/screenshots/`. Design/validation documents: `DESIGN.md`, `MILESTONE1_REVIEW.md`, `RNA_VISUALIZATION_REVIEW.md`, `MILESTONE2_VALIDATION.md`, `MILESTONE2_1_VALIDATION.md`, `MILESTONE2_2_VALIDATION.md`, `FAMILY_OVERVIEW_REFINEMENT.md`, `annotation_assessment/ANNOTATION_ASSESSMENT.md`, `annotation_assessment/PROPOSED_ANNOTATION_SCHEMA.md`, `MILESTONE3B_VALIDATION.md`.

## 2. Completed milestones

| Milestone | Delivered | Key decisions and validation |
|---|---|---|
| POCs | RF00167 (5 chains) and RF00177 (8 chains) mapping POCs; BGSU comparison | vs BGSU `align_chains`: precision 1.0, recall 0.991. BGSU unit id = `PDB|model|chain|comp|authnum|||inscode`; BGSU SeqtoUnitMapping seqpos == label_seq_id |
| 1 (backend) | `rna_explorer` package, RF00167 build: 66 entries, 84 instances, 29 unique sequences, 84/84 validated; clan RF00167/RF01510/RF03165 | consensus position only for match states; insertions carry anchor + offset; no composite confidence score; species recorded per chain |
| 2 (frontend) | three tabs, shared selection in consensus coordinates, Mol* + PDBe RNA Viewer synchronised | identifiers are `label_asym_id` + `label_seq_id`; only residues with coordinates are sent to viewers |
| 2.1 (refinement) | Family Overview layout, stat cards, species chart (5 species, 60/66 entries, B. subtilis 34), instance table paging, display representative | display representative rule: clan-best, X-ray, PDBe 2D available, species-level natural organism, ≥ 95 % matched residues with coordinates, model coverage ≥ family median (0.657), best resolution → **4FE5 chain B** (1.32 Å, 53 eligible). "Use as selected instance" button removed by user request; third-row cards size naturally and align at the top |
| 2.2 (orientation) | `orientation.json` (schema 1.0), consistent 3D orientation across instances, camera fit | Kabsch on C4′ atoms: 51 good / 27 moderate / 2 poor (5e54/B, 5swd/B ≈ 4.6 Å) / 4 insufficient (9nbh); 1y26 1.90 Å, 1y27 1.44 Å, 3rkf 1.00 Å |
| 3A (annotation assessment) | source investigation, coverage tables, mapping issue list, proposed schema | 0 unresolved identifiers across 12,430 references; ambiguities listed in §3 |
| 3B (annotation backend) | sources, mapping, builder, config, 5 schemas, CLI `annotate`, 16 tests, RF00167 outputs | counts reproduce 3A (80 / 496 / 456 / 29 / 37); see `MILESTONE3B_VALIDATION.md` |

## 3. RNA annotations (Milestone 3A/3B)

Sources and identifiers
- **BGSU reference structural elements**: family file `https://rna.bgsu.edu/data/nucleotide_annotations/RF00167_annotations.txt` (one curated reference chain 4FEO|1|B, 67 nt, labels Helix P1, J1/2, Helix P2, L2, J2/3, Helix P3, L3, J3/1) and per-chain transfer `correspondence/nucleotide_annotation?chain=PDB%7C1%7CCHAIN` (author chain id, model 1; 9NBH returns empty). Consensus spans: P1 19–25/81–87, J1/2 26–28, P2 29–35/43–49, L2 36–42, J2/3 50–57, P3 58–63/73–78, L3 65–71, J3/1 79–80.
- **BGSU RNA 3D Motif Atlas**: loop membership `rna3dhub/loops/download/{PDB}` (CSV) and loop→motif-group assignment from the HTML page `rna3dhub/pdb/{PDB}/motifs` (no JSON endpoint). "Valid loops" rows carry the group either as a direct **member** (annotation cell) or via a **homologous**/**geometric** mapping through a release loop (5th "Mapped info" cell). "Problematic loops" carry a reason ("poor quality numbers", "complementary") and no group. Release id from `rna3dhub/motifs/release/{hl,il,j3,j4}/current` (4.15 on 2026-10-15). Do not run FR3D locally.
- **PDBe bound molecules**: `api/v2/pdb/entry/binding_sites/{pdb}/{entity}`: entity-level; `startIndex`/`endIndex` = label_seq_id; `bestChainId` = label_asym_id (metadata only); residues pooled per CCD id with `boundMoleculeId` per range. No Arpeggio atom-level data retrieved.
- Canonical residue id `(pdb_id, label_asym_id, label_seq_id)`; mapping statuses `mapped`, `mapped_no_coordinates`, `insert`, `local_end`, `outside_mapped_region`, `unresolved_identifier`.

Aggregation rules (implemented in `rna_explorer/annotations.py`)
1. `reference_structural_elements.json`: one categorical family-level track derived from the reference chain through our mapping; BGSU per-chain labels kept per instance for highlighting and checked for agreement (5,283/5,283), never counted as observations.
2. `loops.json`: 496 loops (HL 234, IL 82, J3 114, J4 66); one loop id = one observation even across chains (4 multi-chain IL in 9NBH/9NDD); strands preserved; nucleotides outside consensus retained with null position (68 inserts, 168 scaffold/flank). Motif groups: 456 assigned (78 member, 198 homologous, 180 geometric), 29 distinct; 40 unassigned (33 problematic + 7 valid).
3. J3 and J4 are separate tracks; `junction_loops_union` is a set union. The FR3D hairpin HL_56046.1 on the J2/3 strand (50–57) remains a hairpin and is not reclassified.
4. `ligands.json`: 166 records, 37 CCD ids (30 ligand, 5 metal ion, 2 additive), counting unit (pdb_id, entity_id); chain attribution `resolved` only for single-chain entities (139) else `unresolved` (27, 12 multi-chain entities); contacts deduplicated per label_seq_id with merged bound-molecule ids. Classification and artefact rules in `rna_explorer/annotation_config.json`: default track = class `ligand`; `ccd_is_polymer_residue` and `contacts_only_polymer_termini` (window 2) flagged and excluded from the default track (GTP in 9IWF/9IWG/9LKE/9LMC). 63 entries in the default track; 3IVN, 5E54, 9NBH have none. Pocket: consensus 25, 26, 55, 80, 81, 56, 51 contacted in 54–62 of 64 observable evaluated entry-entities.
5. `tracks.json`: per position `observable`, reference element, HL/IL/J3/J4/union/ligand tracks each with counts, `by_motif_group`, `loop_ids`, PP confidence and explicit denominators: observable, source-evaluated, positive, not-annotated-by-source, source-unavailable; `frequency = positive ∧ observable ∧ evaluated / observable ∧ evaluated` (null when 0). Never divide by 84.
6. `annotations/<seq_name>.json`: every polymer residue with reference element, loops, bound molecules (`attribution: chain | entity_unresolved`), plus per-source status `evaluated | unavailable`.

Schemas: `rna_explorer/schemas/{reference_structural_elements,loops,ligands,tracks,annotation_instance}.schema.json`, version 1.0. Documentation: `annotation_assessment/PROPOSED_ANNOTATION_SCHEMA.md` (implemented schema), `MILESTONE3B_VALIDATION.md`.

Decisions awaiting scientific review: (a) GTP terminal-contact exclusion (config flag); (b) whether homologous/geometric group mappings should count like direct membership; (c) presentation of entity-level contacts on multi-chain entities; (d) `loops.json` size (3.9 MB) for the viewer.

## 4. Outstanding scientific question: chain-specific RNA–ligand contacts
Current limitation: PDBe `binding_sites` pools residues per ligand over all copies in the entity, so for the 12 multi-copy entities (3DS7, 3IVN, 3RKF ×4, 5E54, 5SWD, 6VWT, 6VWV, 9IWF, 9IWG, 9LMC, 9NBH ×4, 9NDD ×4) the contacting chain is unknown; `bestChainId` is kept as metadata only and instance counts are lower bounds. Options to investigate:
- PDBe bound-molecule / ligand-site endpoints (`api/pdb/entry/binding_sites/{pdb}` v1 per chain, `api/v2/pdb/entry/bound_molecules/{pdb}`, `pdbe/graph-api/pdb/bound_ligand_interactions/{pdb}/{bm_id}` which carries Arpeggio interactions per bound molecule with chain ids).
- RCSB Data API (`core/nonpolymer_entity_instance`, `rcsb_target_neighbors` with `struct_asym_id` / `label_seq_id`).
- Running Arpeggio locally on the cached updated mmCIF files (atom-level, interaction types; heaviest option).
Whatever source is chosen, keep entry-entity counts as the primary statistic and add chain-level contacts as a separate, provenance-labelled layer; never convert `entity_unresolved` records into chain contacts without new evidence.

## 5. Next milestone: 3C, interactive annotation tracks
Frontend only; consume the 3B files without changing the backend. Suggested scope:
- Load `tracks.json` and `reference_structural_elements.json` at start; `loops.json`/`ligands.json`/`annotations/<seq_name>.json` lazily (per selected instance).
- Add tracks to `ConsensusTracks.tsx`: reference elements (categorical lane, labelled as transferred from 4FEO), hairpin, internal, J3, J4 (separate lanes; optional union lane), ligand-binding sites (default classes; toggle for metals/additives and flagged molecules). Show numerator/denominator and group breakdown on hover; drill-down to loop ids and CCD ids.
- Synchronise with 2D/3D: clicking a track segment selects consensus positions (existing `SelectionContext`); per-instance highlighting uses `annotations/<seq_name>.json` residues (reference element, loop membership, contacts with attribution caveat).
- Distinguish visually: not observable / not annotated by the source / source unavailable; entity-level ligand contacts vs chain-resolved.
- Extend TS types (`src/types/backend.ts`) from the new schemas; add vitest + Playwright coverage and screenshots to `docs/screenshots/`.

## 6. Future clustering investigation (not started)
Structures tab still shows `PendingPanel` placeholders for similarity matrix, clustering/dendrogram and groups. Planned investigation:
- BGSU geometric discrepancy (from `rna3dhub` correspondence / equivalence classes) as one similarity source.
- US-align TM-score / RMSD on the consensus-matched C4′ (or all-atom) sets; the Kabsch machinery in `orientation.py` and the per-instance C4′ extraction can be reused.
- Similarity restricted to the intersection of observed Rfam consensus positions per pair (explicit denominators, as in the tracks).
- Multi-family benchmarks (RF00167, RF00177 and additional families); RAD-based 16S rRNA validation on RF00177: RAD is the Ribosome Analysis Database (Whitford group, Northeastern University, https://radtool.rc.northeastern.edu/), which lists nearly every PDB ribosome structure with small-subunit rotation/tilt angles, LSU listings and alignment metadata computed with RADtool (open-source VMD plugin; downloads at `/downloads.html`, full database at `/data.html`, webtools "Find Similar Structures" / "Compare Two" at `/webtools.html`). Use it as an independent reference for 16S rRNA (RF00177) similarity/clustering results, building on `poc_rf00177/bgsu_validation/`.
- Do this in a separate git worktree (§7) so the main app stays stable.

## 7. Pending technical tasks
- `git init` the project; review `.gitignore` (currently `cache/`, `data/families/*/work/`, `__pycache__/`, `.pytest_cache/`; frontend has its own for `node_modules`, `dist`, `test-results`, `playwright-report`, `.vite`). Decide whether `data/families/` (15 MB, incl. 3B files ~9 MB), `tools/infernal-1.1.5/` (20 MB binary build), `annotation_assessment/cache/` and `docs/screenshots/` are committed; `.DS_Store` should be ignored.
- First commit, then a separate worktree (e.g. `git worktree add ../aggregated_rna_view-clustering clustering`) for the clustering research.
- `frontend/ARCHITECTURE.md` refresh (pdbe-molstar, PdbRnaViewer2D, SelectionContext details).

## 8. Important constraints
- Scientific provenance on every output: source URL, retrieval date, Rfam release and CM checksum, Motif Atlas release, pipeline/config versions.
- Explicit numerators and denominators; "not annotated by the source" is never "confirmed absence"; missing source data is never a negative.
- Annotation uncertainty is preserved (entity-level ligands, transferred reference labels, unassigned/problematic loops, flagged artefacts); nothing is silently discarded.
- The validated Rfam mappings in `data/families/RF00167/mappings/` are the single coordinate authority; annotations are resolved through them, never re-aligned, and insertions/local ends/unmapped residues are never forced onto consensus positions.
- Identifiers: `label_asym_id` + `label_seq_id` internally; author numbering only for display and for BGSU unit ids.
- Frontend and backend changes are kept to the milestone's scope; representative/display rules and layout decisions above were confirmed by the user.

## Memory
Claude Code auto-memory for this project lives in `~/.claude/projects/-Users-sria-Desktop-aggregated-rna-view/memory/` (`rfam-mapping-poc.md`, `species-per-chain.md`, `annotation-sources-rf00167.md`).
