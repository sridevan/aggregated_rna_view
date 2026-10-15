# Milestone 3B validation: RNA annotation backend (RF00167)

Date: 2026-10-15. Scope: backend only (`rna_explorer/annotation_sources.py`, `annotation_mapping.py`, `annotations.py`, `annotation_config.json`, five JSON schemas, CLI `annotate`). No frontend change. Build command:

    python3 -m rna_explorer annotate RF00167 --raw-dir annotation_assessment/cache

The build resolves every external identifier through the existing validated mapping files (`mappings/<seq_name>.json`); Infernal is not rerun. Raw responses from Milestone 3A seed the cache, so the build is reproducible offline (`--offline`). Motif Atlas release pages were fetched once (release 4.15).

## 1. Data coverage

| Source | Unit | Evaluated | Unavailable | Notes |
|---|---|---|---|---|
| BGSU reference family file | family | 1 file, 4FEO chain B, 67 nt | – | 67 consensus positions 19–87 labelled; 0 reference nucleotides off-consensus; 0 conflicts |
| BGSU per-chain element transfer | instance | 80 / 84 | 4 (9NBH A–D: empty response) | 5,283 labels, 5,283 agree with the family track, 0 disagree |
| Motif Atlas loops + motifs page | PDB entry | 66 / 66 (84 instances) | 0 | 496 loops; 0 loops on the motifs page missing from the download |
| PDBe binding_sites | (entry, entity) | 66 / 66 | 0 | 166 bound-molecule records, 37 CCD ids |

Instances/entries with an evaluated source but no annotation at a position are reported as "not annotated by the source"; the 9NBH chains are "source unavailable" for reference elements and are never counted as negatives.

## 2. Annotation counts (vs Milestone 3A)

| Quantity | 3A | 3B | Comment |
|---|---|---|---|
| Instances with reference-element labels | 80 | 80 | identical |
| Motif Atlas loops | 496 | 496 | HL 234, IL 82, J3 114, J4 66 |
| Loops with motif group | 456 | 456 | **new**: 78 `member`, 198 `homologous`, 180 `geometric` (mapped through a release loop) |
| Loops without group | 40 | 40 | 33 `problematic` (31 "poor quality numbers", 2 "complementary") + 7 valid but unassigned |
| Distinct motif groups | 29 | 29 | HL 10, IL 7, J3 9, J4 3 |
| Multi-chain loops | 4 | 4 | IL in 9NBH/9NDD scaffolds, one observation each |
| Loop nucleotide memberships | 5,655 | 5,655 | mapped 5,419, insert 68, outside region 168, unresolved 0 |
| CCD ids | 37 | 37 | ligand 30, metal ion 5, additive 2 |
| Entries with a default-track ligand | 63 | 63 | GTP (4 entries) now flagged and excluded, see §6 |
| Binding-residue contacts | 1,492 rows | 1,370 contacts | **deviation explained**: 3A rows were per bound-molecule range and per chain; 3B deduplicates per (entry, entity, CCD, label_seq_id) and merges bound-molecule ids (1,358 mapped + 12 outside region) |
| IL instances / entries | 50 / 46 | 58 / 48 | 3A counted only chains with consensus-mapped IL nucleotides; 3B lists every chain participating in an IL, including scaffold/insert-only loops (20 loops have no consensus position) |

## 3. Mapping validation
- 12,430 external nucleotide references (5,283 + 5,655 + 1,492 ranges expanded), 0 unresolved identifiers (`annotation_summary.json: unresolved_identifiers`).
- Insertions (68 loop nucleotides in 9IWF/9IWG/9LMC/9LJN/9LKC), local ends (0) and nucleotides outside the mapped region (168 loop nucleotides on the 9NBH/9NDD scaffolds and 6P2H flanks; 12 contacts) are retained in `loops.json`, `ligands.json` and the per-instance files with `rfam_consensus_pos: null` and are excluded from consensus aggregation.
- The BGSU per-chain transfer agrees with our own reference transfer at all 5,283 labelled nucleotides, confirming that the 3A mapping and BGSU's Rfam-based correspondence coincide for this family.
- PDBe `startIndex`/`endIndex` are label_seq_id; `bestChainId` is label_asym_id (3A validation reused).

## 4. Aggregation rules (as implemented)
1. Reference elements: single categorical track derived from 4FEO|1|B through the mapping; instance labels kept for highlighting only.
2. Loops: one loop id = one observation; per position, counts of distinct loops, instances (chains with a consensus-mapped nucleotide) and entries; `by_motif_group` incl. `unassigned`; PP confidence of contributing nucleotides.
3. J3 and J4 separate; `junction_loops_union` is a set union (e.g. consensus 25: J3 75 instances, J4 66, union 79).
4. The FR3D hairpin HL_56046.1 on the J2/3 strand (50–57) stays a hairpin loop; J3 covers the same positions. Both are reported.
5. Ligands: counting unit (pdb_id, entity_id); chain attribution resolved only for single-chain entities (139 records) and unresolved for multi-chain entities (27 records); default track = class `ligand` minus flagged artefacts; metals/additives retained with counts in `by_class_n_entry_entities`.
6. Denominators per position: observable, source-evaluated, positive, not-annotated-by-source, source-unavailable; frequency = positive / (observable ∧ evaluated), null when 0.

Example, consensus 25 (ligand pocket): 62 of 64 observable and evaluated entry-entities have a ligand contact (frequency 0.9688); consensus 37 (L2): 80 hairpin-loop instances of 80 observable (HL_59138.3 52, HL_35704.1 21, HL_55605.1 2, unassigned 5).

## 5. Tests
`python3 -m pytest -q`: 29 passed (13 pre-existing + 16 new).
- `tests/test_annotations.py` (synthetic family, fake fetcher): identifier parsing and resolution incl. insertion/no-coordinate/unknown chain; reference transfer and instance availability; discontinuous strands and insertions in loops; multi-chain loop = one observation; J3/J4 overlap union; motif-group counts and homologous assignment; entity-level attribution, contact deduplication and class exclusion; entry-level ligand denominators; loop denominators and missing sources (frequency null vs 0); schema validation and traceability.
- `tests/test_annotations_rf00167.py`: 3A counts, schema validation of all outputs, every track count recomputed from `loops.json`/`ligands.json`, denominator identities at all 102 positions.

## 6. Known limitations and decisions for scientific review
1. **GTP at the 5′ terminus** (9IWF, 9IWG, 9LKE, 9LMC): contacts only polymer residue 1–2 (consensus 19–20). Flagged `contacts_only_polymer_termini` and excluded from the default track by config (`exclude_from_default_track: true`). Reviewers may prefer to keep it visible.
2. **Motif-group assignment types**: 378 of 456 assignments are `homologous`/`geometric` mappings onto a release loop, not direct group membership. Tracks currently count all three types alike; a per-type filter is available in `loops.json`.
3. **Problematic loops** (33) carry no motif group by Atlas policy; they are counted as `unassigned`, not as absent loops.
4. **Reference track** is alignment-transferred (4FEO); observed base-pairing per structure is not part of this milestone.
5. **Entity-level ligands**: for 12 multi-chain entities, per-instance records show contacts on every chain of the entity with `attribution: entity_unresolved`; frontend should present them as entity-level, not chain-level.
6. **Classification lists** are generic (`annotation_config.json`); NCO (cobalt hexammine) is classed as a metal ion, ACT/BR as additives. No water records in this endpoint.
7. **Motif Atlas release** is read from the release page at build time (4.15); motif-group ids are version-suffixed and may change between releases.
8. BGSU's per-chain endpoint takes the author chain id and model 1 only.

## 7. Output file sizes (RF00167)

| File | Size |
|---|---|
| `reference_structural_elements.json` | 28 KB |
| `loops.json` | 3.9 MB (496 loops with full nucleotide membership) |
| `ligands.json` | 0.7 MB |
| `tracks.json` | 0.6 MB (102 positions × 7 tracks with loop id lists) |
| `annotations/*.json` | 4.0 MB total, 84 files (≈ 47 KB each) |
| `annotation_summary.json` | 9 KB |

## 8. Readiness for frontend integration
- Tracks are position-indexed (1..CLEN) and carry explicit numerators/denominators, group breakdowns and loop ids for drill-down; per-instance files give residue-level labels for synchronized 2D/3D highlighting through the existing consensus translation.
- `loops.json` is large because every nucleotide carries its full ResidueRef; the frontend should load it lazily or a trimmed variant can be exported in 3C.
- All files validate against `rna_explorer/schemas/*.schema.json` (schema_version 1.0).
