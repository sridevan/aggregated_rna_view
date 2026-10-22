# Aggregated RNA 3D View: project handover

Last updated: 2026-10-16 (Milestone 3C UI reorganisation done; chain-specific ligand investigation and Arpeggio pilot done). Owner: asridevan86@gmail.com. Working directory `/Users/sria/Desktop/aggregated_rna_view` (macOS arm64, Python 3.11.9 via pyenv, Node 20.14 via nvm). Git repository on `main` since 2026-10-18, pushed to https://github.com/sridevan/aggregated_rna_view (see §7). **Sections 2–6 describe the state on 2026-10-16; see §0 for what changed since.**

Goal: a PDBe-KB-style aggregated view of experimentally determined RNA structures per Rfam family, using Rfam covariance-model (CM) match-state positions as the family coordinate system. RF00167 (purine riboswitch) is the development family; RF00177 (16S rRNA) was used in the proofs of concept.

## 0. Status update (2026-10-18)
Done since the sections below were written (details in the linked reports; facts per family in the auto-memory files):
- Annotations milestone completed: `MILESTONE3D_VALIDATION.md`, `MILESTONE3E_VALIDATION.md`, `MILESTONE3F_VALIDATION.md` (SASA, fpocketR pockets, nonstandard nucleotides).
- Milestone 4, structural clustering and the Structures tab: `MILESTONE4_STRUCTURES_VALIDATION.md` (also records the later UI sections: Annotations zoom, Position numbering). §6 below is superseded.
- Further families, all live: RF00162 SAM-I, RF00059 TPP, RF00177 16S (curated 60 chains via `build --selection data/selections/RF00177_60_chains.csv`; ModelServer per-chain coordinates; no pockets). Reports in `experiments/rf00162_phase1/`, `experiments/rf00059_tpp/`, `experiments/rf00177_16s/`.
- Deployed on Netlify: https://aggregated-rna-view.netlify.app/ (home page with representative thumbnails from `npm run thumbnails`) and /{rfam_acc} family pages; `cd frontend && npm run deploy`.
- Mapping methodology and limitations: `docs/RFAM_PDB_MAPPING.md`.

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
| `sasa.py`, `pockets.py`, `modifications.py` | 3D SASA, 3E fpocketR pockets, 3F nonstandard nucleotides (deposited mmCIF polymer components; PDBe cross-check, CCD/NAKB parents) |
| `tools/infernal-1.1.5/bin/cmalign` | Infernal built from source; found via `RNA_EXPLORER_CMALIGN`, then this path, then `PATH` |

Commands:

    pip install requests jsonschema pytest
    python3 -m rna_explorer build RF00167 --cpu 4 -v          # fetch → instances → cmalign → mapping → validation → stats → export
    python3 -m rna_explorer orient RF00167                     # orientation.json
    python3 -m rna_explorer annotate RF00167 --raw-dir annotation_assessment/cache   # annotation files; --offline forbids network
    python3 -m rna_explorer modifications RF00167 --seed-dir experiments/modifications_RF00167/raw   # modifications.json, modification_tracks.json
    python3 -m pytest -q                                       # 29 tests

### Data (`data/families/RF00167/`, served to the frontend via symlink `frontend/public/families → ../../data/families`)
`family.json` (description, comment, `rfam_url`, clen 102, consensus sequence, SS_cons, clan members, counts, `representative` = 9ndd/A, `display_representative` = 4fe5 chain B), `instances.json` (84 rows incl. organism, species, tax id, `pdbe_2d_available`), `consensus_stats.json`, `validation.json`, `orientation.json`, `mappings/<seq_name>.json` (seq_name = `<pdb>_<entity>_<auth_chain>_<label_chain>`), `residue_mapping.csv`, `consensus_columns.csv`, `work/` (alignment, ignored by git), and the 3B annotation files `reference_structural_elements.json`, `loops.json`, `ligands.json`, `tracks.json`, `annotation_summary.json`, `annotations/<seq_name>.json`.

### Frontend (`frontend/`, React 19 + TypeScript + Vite 8 + Tailwind 4)
- 3D: `pdbe-molstar` 3.12 (`pdbe-molstar/lib/viewer`, needs `sass-embedded`; rolldown binding installed with `--no-save @rolldown/binding-darwin-arm64`). `components/PdbeMolstarViewer.tsx`: one plugin per mount, structure switch via `visual.update(…, false)` + hierarchy removal, white background, orientation transform via `TransformStructureConformation`, region/complete modes, custom camera fit (fixed view axes, exact perspective distance, margin 1.06, 40 px bands, `Camera.SetSnapshot` after two frames), "Reset orientation / Fit" and "Focus selection", 60 s watchdog + Retry.
- 2D: vendored PDBe RNA Viewer 0.3.0 in `public/vendor/pdb-rna-viewer` (Apache-2.0), `components/PdbRnaViewer2D.tsx` (prefetches PDBe topology/FR3D, explicit unavailable state, SVG fit).
- State: `state/family.tsx`, `state/selection.tsx` (one active selection: `position` | `residue` | `group` (future BGSU loops, exact membership + segments); consensus position preserved across instances; selected instance held separately); `lib/consensus.ts` translate/inverse; `lib/species.ts`; `lib/orientation.ts`.
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
| 3D (SASA, 2026-10-16) | `rna_explorer/sasa.py`: isolated-chain per-nucleotide SASA (FreeSASA Shrake–Rupley, probe 1.4 Å, element radii, CCD-based completeness, altloc policy) aggregated as median per consensus position; `sasa.json`/`sasa_tracks.json`; Annotations category "Solvent accessibility" with one quantitative track and two downloads | 84/84 instances, 6,133 complete nucleotides, independent numpy check 0.06 %; see `MILESTONE3D_VALIDATION.md`, `experiments/sasa_pilot/PILOT_REPORT.md` |
| 3E (predicted pockets UI, 2026-10-16) | Annotations category "Predicted binding pockets" with one "Pocket occurrence" frequency track from `pocket_tracks.json` (filter all Score ≥ 0), two downloads | `MILESTONE3E_VALIDATION.md`, screenshot `docs/screenshots/m3e_annotations_pockets.png` |
| 3F (nonstandard nucleotides, 2026-10-16; final Annotations category) | `rna_explorer/modifications.py` + CLI `modifications`; `modifications.json`/`modification_tracks.json`; Annotations category "Nonstandard nucleotides" with one fixed-colour "Occurrence" track, two downloads; residue tooltip flips above near the viewport bottom | 13 observations (7 CCD, 9 instances, 8 entries), 6 on the axis at 22/32/35/84/87, 7 terminal off-axis (incl. 23G 3IVN); occurrence = nonstandard / determined identity (coordinates not required); see `experiments/modifications_RF00167/REPORT.md`, `MILESTONE3F_VALIDATION.md` |
| 4 (Structures tab, RF00167 clustering, 2026-10-16) | `rna_explorer/structure_clustering.py` + CLI `clustering` (UPGMA ported from the antibody_antigen_interface demo, fixed cuts in `structure_clustering_config.json`, medoids) → `structure_clustering.json`; Structures tab rebuilt after the PDBe-KB Complexes prototype: compact/expanded table + Structural groups panel, collapsed dendrogram (76 > 50) with drill-down for groups ≤ 50, metric selector (geometric discrepancy default, TM-score), 3Dmol.js alignment view on orientation.json transforms, associated features | D: 72 + 4 (stable, cut 0.495); TM: 51/17/4/4 (provisional, cut 0.327); 76 of 84 instances clustered; see `MILESTONE4_STRUCTURES_VALIDATION.md` |
| 3C (UI reorganisation, 2026-10-16) | Family Overview = shared instance selector + compact coverage/instance tracks + 2D + 3D (`StructureExplorer.tsx`); single-nucleotide selection with Mol* highlight + camera focus and 2D mirror; Annotations = independent PDBe-KB-style fixed-size marks coloured by frequency ("RNA loop annotations" HL/IL/JL, Ligand binding sites per CCD), scrollable categories, concise tooltips, category-header ⓘ popovers, Download dropdown (loops/ligands/tracks JSON); `ScaleDef` per track ready for pocket/RSA tracks | 26 Playwright + 12 vitest pass; see `MILESTONE3C_VALIDATION.md`, screenshots `docs/screenshots/m3c_*.png` |

## 3. RNA annotations (Milestone 3A/3B)

Sources and identifiers
- **BGSU reference structural elements**: family file `https://rna.bgsu.edu/data/nucleotide_annotations/RF00167_annotations.txt` (one curated reference chain 4FEO|1|B, 67 nt, labels Helix P1, J1/2, Helix P2, L2, J2/3, Helix P3, L3, J3/1) and per-chain transfer `correspondence/nucleotide_annotation?chain=PDB%7C1%7CCHAIN` (author chain id, model 1; 9NBH returns empty). Consensus spans: P1 19–25/81–87, J1/2 26–28, P2 29–35/43–49, L2 36–42, J2/3 50–57, P3 58–63/73–78, L3 65–71, J3/1 79–80.
- **BGSU RNA 3D Motif Atlas**: loop membership `rna3dhub/loops/download/{PDB}` (CSV) and loop→motif-group assignment from the HTML page `rna3dhub/pdb/{PDB}/motifs` (no JSON endpoint). "Valid loops" rows carry the group either as a direct **member** (annotation cell) or via a **homologous**/**geometric** mapping through a release loop (5th "Mapped info" cell). "Problematic loops" carry a reason ("poor quality numbers", "complementary") and no group. Release id from `rna3dhub/motifs/release/{hl,il,j3,j4}/current` (4.15 on 2026-10-15). Do not run FR3D locally.
- **PDBe bound molecules**: `api/v2/pdb/entry/binding_sites/{pdb}/{entity}`: entity-level; `startIndex`/`endIndex` = label_seq_id; `bestChainId` = label_asym_id (metadata only); residues pooled per CCD id with `boundMoleculeId` per range. No Arpeggio atom-level data retrieved.
- Canonical residue id `(pdb_id, label_asym_id, label_seq_id)`; mapping statuses `mapped`, `mapped_no_coordinates`, `insert`, `local_end`, `outside_mapped_region`, `unresolved_identifier`.

Aggregation rules (implemented in `rna_explorer/annotations.py`)
1. `reference_structural_elements.json`: one categorical family-level track derived from the reference chain through our mapping; BGSU per-chain labels kept per instance for highlighting and checked for agreement (5,283/5,283), never counted as observations.
2. `loops.json`: 496 loops (HL 234, IL 82, J3 114, J4 66); one loop id = one observation even across chains (4 multi-chain IL in 9NBH/9NDD); strands preserved; nucleotides outside consensus retained with null position (68 inserts, 168 scaffold/flank). Motif groups: 456 assigned (78 member, 198 homologous, 180 geometric), 29 distinct; 40 unassigned (33 problematic + 7 valid).
3. J3 and J4 are separate tracks; `junction_loops_union` is a set union. **Merged loops (2026-10-21, user decision):** BGSU also emits a merged loop where two loops share one canonical closing pair and are linked across it; the download does not flag it, so a loop whose nucleotide set is exactly the union of two other non-problematic loops sharing exactly two nucleotides is marked `decomposition.role = merged` and its parts `superseded_by`. The HL/IL/J tracks and the Family Overview loop rows count and draw only loops not superseded; all records stay in `loops.json` (`superseded_loop_ids` per track position). So the FR3D hairpin HL_56046.1 on the J2/3 strand (closed by the isolated pair G46·C53) is now a superseded part of the merged J3, and an internal loop BGSU merged into a junction is counted as part of that junction. Previous position (HL and J3 both kept): `experiments/loop_overlap_pos53/SCIENTIFIC_REVIEW.md`; validation: `MERGED_LOOPS_VALIDATION.md`.
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
**Investigated 2026-10-15** (`experiments/chain_specific_ligands/`, read `RECOMMENDATIONS.md`): PDBe graph-api `bound_molecule_interactions`/`bound_ligand_interactions` give chain-specific Arpeggio interactions but only for ligand copies in PDBe's *preferred assembly* (3RKF: chain B only; `bestChainId` can be wrong); RCSB `rcsb_ligand_neighbors` gives every copy's residues within 5 Å with label ids; `ligand_monomers.struct_asym_id` joins the two. Local Arpeggio works (≈1 s/entry) but is not required. **Arpeggio pilot on deposited archive mmCIF** (`experiments/chain_specific_ligands/arpeggio_pilot/`): 66 RF00167 entries in 68 s, chain-level tracks, symmetry-omission check (409 symmetry-image contacts, all ions/additives/terminal GTP); 7K00 70S 16S chain works (paromomycin at h44) but took 6.7 h with all 329 ligands selected → ligand pre-filter required for ribosomes. Proposed: RCSB proximity as the chain-level evidence layer, PDBe classes where available, v2 entity-level counts kept as primary until validated on all 66 entries.
Whatever source is chosen, keep entry-entity counts as the primary statistic and add chain-level contacts as a separate, provenance-labelled layer; never convert `entity_unresolved` records into chain contacts without new evidence.

## 5. Milestone 3C (done 2026-10-16) and what remains
Implemented as described in `MILESTONE3C_VALIDATION.md`. Done later: solvent accessibility (3D). Not yet done (candidates for 3E): predicted binding pockets (feasibility `experiments/fpocketr_feasibility/REPORT.md`; approved pilot done `experiments/fpocketr_pilot/REPORT.md`: isolated env `tools/fpocketr-env` (fpocket 4.0.3, fpocketR 1.3.6, x86-64 under Rosetta), holo chains 4FE5/1Y26/1Y27 give exactly the purine pocket (consensus 25,26,51,54,55,56,79,80,81; Known by fpocketR's own criteria), apo 5E54 shifted, 3IVN/9NDD none; full RF00167 run done: `rna_explorer/pockets.py`, `data/families/RF00167/pockets.json` + `pocket_tracks.json`, report `experiments/fpocketr_full_RF00167/REPORT.md` (65/84 chains with a pocket, 58/75 holo chains recover the purine site as pocket 1, 17 holo chains yield none; recommend one track, filter all Score ≥ 0); frontend track done, see 3E), chain-specific ligand evidence layer.
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
- Git: done 2026-10-18 (`main`, pushed to GitHub `sridevan/aggregated_rna_view`). Before the first push the history was cleaned of raw downloads (POC mmCIFs/API dumps, RF00167 `work/`); these stay on disk and are ignored. `.gitignore` excludes caches, scratch `work/` folders, downloaded structure files under `experiments/` (full mmCIFs are never stored), the regenerable RF00177 outputs over 10 MB (`loops.json`, `sasa.json`, `residue_mapping.csv`) and RF00177 `annotations/`.
- `frontend/ARCHITECTURE.md` refresh (pdbe-molstar, PdbRnaViewer2D, SelectionContext details).

## 8. Important constraints
- SASA (3D): RNA-only, deposited chain, first model; only `complete` nucleotides enter the median track; never normalise by instance count or treat missing residues as zero.

- Scientific provenance on every output: source URL, retrieval date, Rfam release and CM checksum, Motif Atlas release, pipeline/config versions.
- Explicit numerators and denominators; "not annotated by the source" is never "confirmed absence"; missing source data is never a negative.
- Annotation uncertainty is preserved (entity-level ligands, transferred reference labels, unassigned/problematic loops, flagged artefacts); nothing is silently discarded.
- The validated Rfam mappings in `data/families/RF00167/mappings/` are the single coordinate authority; annotations are resolved through them, never re-aligned, and insertions/local ends/unmapped residues are never forced onto consensus positions.
- Identifiers: `label_asym_id` + `label_seq_id` internally; author numbering only for display and for BGSU unit ids.
- Frontend and backend changes are kept to the milestone's scope; representative/display rules and layout decisions above were confirmed by the user.

## Memory
Claude Code auto-memory for this project lives in `~/.claude/projects/-Users-sria-Desktop-aggregated-rna-view/memory/` (`rfam-mapping-poc.md`, `species-per-chain.md`, `annotation-sources-rf00167.md`).
