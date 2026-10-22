# Reported sequence differences of RNA chains from mmCIF: method

This document describes how we extract, **per polymer chain**, the sequence differences that a PDB entry records between each polymer's deposited sequence and its reference sequence. It is written so the same process can be rebuilt independently, for example to extend the PDBe API `mutated_AA_or_NA` to nucleic acids.

Reference implementation: `rna_explorer/sequence_differences.py`.
- Part A, `extract_chain_differences`: independent of Rfam and of this project.
- Part B, `place_on_family` and the build function: our Rfam layer.

Tests: `tests/test_sequence_differences.py`. Investigation: `experiments/mutations_investigation/REPORT.md`.

## Why this is needed

The PDBe API `POST /api/v2/pdb/entry/mutated_AA_or_NA` currently reports protein chains only. Its description says: "at present it does not provide information about mutated nucleotides in RNA or DNA chains".
- For the 289 entries of our six RNA families it returned 640 records, all on protein chains.
- Entries with known RNA changes, such as 3FO4 (C74U), 4LX6 and 9IWG, return HTTP 404.

The information is deposited. It is in the mmCIF category `_struct_ref_seq_dif`, which the API already uses for proteins: for 3V7E and 3K0J the API records equal the UniProt-referenced protein rows of that category exactly.

## Part A: per-chain extraction from one entry (polymer-type agnostic)

### Input

The entry's mmCIF; we use the PDBe updated mmCIF, `https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}_updated.cif`. Categories read:

| Category | Items | Use |
|---|---|---|
| `_struct_ref_seq_dif` | `align_id`, `pdbx_pdb_strand_id`, `seq_num`, `pdbx_auth_seq_num`, `pdbx_pdb_ins_code`, `mon_id`, `db_mon_id`, `pdbx_seq_db_seq_num`, `details`, `pdbx_ordinal` | the difference rows |
| `_struct_ref_seq` | `align_id`, `ref_id`, `pdbx_strand_id` | alignment → reference; fallback strand |
| `_struct_ref` | `id`, `entity_id`, `db_name`, `db_code`, `pdbx_db_accession` | reference → entity; reference database |
| `_pdbx_poly_seq_scheme` | `asym_id`, `entity_id`, `seq_id`, `mon_id`, `pdb_strand_id`, `pdb_seq_num`, `pdb_ins_code` | strand → label chain; identity and numbering check |
| `_entity_poly` | `entity_id`, `type` | polymer type (filtering only) |

Treat `?`, `.` and empty as missing. Columns other than `_struct_ref_seq_dif.align_id` may be absent from a file and must be read as optional.

### Algorithm, for every `_struct_ref_seq_dif` row

1. **Entity:** `align_id` → `_struct_ref_seq.ref_id` → `_struct_ref.entity_id`.
   - Keep `db_name`, `db_code` and `pdbx_db_accession`: the differences are relative to this reference.
   - If unresolved: `chain_status = entity_unresolved`.
2. **Author strand:** `pdbx_pdb_strand_id`, else `_struct_ref_seq.pdbx_strand_id` of the alignment.
3. **Label chain (`struct_asym_id`):** look up (`entity_id`, strand) in `_pdbx_poly_seq_scheme` (`entity_id`, `pdb_strand_id`) → `asym_id`.
   - If not found: `chain_status = chain_unresolved`.
   - **Never assign the row to another chain**, even when the entity has only one other copy.
4. **Residue:**
   - `label_seq_id = seq_num`;
   - `deposited = mon_id`, `reference = db_mon_id`;
   - author identifiers: `pdbx_auth_seq_num` and `pdbx_pdb_ins_code`;
   - `reference_seq_num = pdbx_seq_db_seq_num`.
5. **Change kind:**
   - `substitution`: deposited and reference both present.
   - `added_residue`: deposited only, e.g. an added 5′ G for transcription.
   - `absent_residue`: reference only, no `seq_num`; a reference residue missing from the construct.
6. **Classification:** `details`, verbatim. Do not change case or wording. Observed values on RNA are "Engineered mutation", "Conflict", "variant", "Insertion", "deletion" and "SEE REMARK 999"; on proteins also "Expression tag", "Linker", "Initiating methionine" and others.
7. **Identity and numbering check** against `_pdbx_poly_seq_scheme` at (`asym_id`, `seq_id = label_seq_id`):
   - `agrees`: `mon_id` equals the deposited residue, and `pdb_seq_num` / `pdb_ins_code` equal the author identifiers when the row gives them.
   - `identity_conflict`: different residue.
   - `author_number_conflict`: same residue, different author numbering.
   - `not_in_scheme`: no such position.
   - `no_deposited_residue`: an absent residue.
   - Keep rows that fail the check, with the status.

**Output per row:**
- entry and entity, polymer type;
- `auth_asym_id` and `label_asym_id`, `label_seq_id`, author number and insertion code;
- reference and deposited residues, change kind;
- classification;
- reference database, code and accession;
- chain status and check status;
- `pdbx_ordinal`.

**Mapping to the API's fields:**

| API field | Source |
|---|---|
| `entity_id` | resolved entity |
| `chain_id` | author strand |
| `struct_asym_id` | label chain |
| `residue_number` | `label_seq_id` |
| `author_residue_number`, `author_insertion_code` | author identifiers |
| `chem_comp_id` | deposited residue |
| `mutation_details.from` | reference residue |
| `mutation_details.to` | deposited residue |
| `mutation_details.type` | `details` |

`absent_residue` rows have no `residue_number` and would need a representation that the current API schema does not have.

### Copies

Each copy of an entity has its own rows (one `_struct_ref_seq` alignment per strand). In our data every multi-copy entity lists identical rows for each copy: 6VWT, 6VWV, 9IWF, 9IWG and 8F4O, chains A and B. Check this per entity and report copies that disagree, but do not copy rows between chains.

### Interpretation (important for any consumer)

- **Relative to the depositor's reference.** Differences are measured against the reference sequence in `_struct_ref`, often a genomic GenBank region. That is not a canonical wild type and not a family consensus.
- **Labels are free text.** In our data:
  - "Conflict" (5C7U) and "variant" (6VWT/6VWV) rows are a redesigned P1 helix, i.e. construct design.
  - Added transcription nucleotides are labelled "Engineered mutation" (2B57) or "Insertion" (8F4O).
  - Present the label as deposited, not as a verified experimental mutation.
- **Absence is not evidence.** An entity whose `_struct_ref` points to its own PDB entry (112 of 122 RNA entities in three of our families) records no differences, whatever was engineered. A chain without rows has no recorded difference; that doesn't mean it is unchanged.

## Part B: our Rfam layer (optional for PDBe)

**Filtering.** Keep rows whose (entry, entity) is a family RNA entity and whose resolved label chain is a family chain. Rows of other polymers are counted by polymer type in `counts.rows_on_other_entities_by_polymer_type` and not shown.

**Placement.** Look up `label_seq_id` in that chain's validated residue → Rfam mapping. A row gets a consensus position only when all of these hold:
- the Part A check is `agrees`;
- the mapping's `residue_name` equals the deposited residue;
- its author number and insertion code equal the row's;
- the residue is a match state.

Otherwise the position stays null, with one of these statuses:
- `no_consensus_position:<insert|not_submitted_flank|local_end|…>`;
- `not_in_mapping`;
- `no_deposited_residue`;
- `mapping_identity_conflict`;
- `mapping_author_conflict`.

These rows are listed, not dropped.

**Aggregation**, per consensus position and classification:
- observation count;
- distinct RNA chain instances;
- distinct PDB entries;
- `reference>deposited` change counts and change kinds;
- the observation ids.

No frequency is computed, because absence is not a negative.

**Display.** In the Annotations tab, one presence track per classification found on the axis, with the classification text as the label. Pinned details list every observation; each opens its own chain in the Family Overview and highlights the position there.

## Validation (2026-10-22)

**Three families with cached mmCIF:**
- RF00167: 85 rows; 77 on consensus positions, 8 flank rows.
- RF00162: 2 rows, a `deletion` with no residue and a `Conflict` at an insert.
- RF00059: 8 `Insertion` rows, 2 on consensus positions.
- All 94 rows with a deposited residue agree with both the scheme and our mapping (identity and author numbering).
- No row is unresolved, and copies are consistent.
- Results for all six families (192 RNA rows in 289 entries; none unresolved; no copy disagreement; no identity or numbering conflict): `SEQUENCE_DIFFERENCES_VALIDATION.md`.

**Spot checks:**
- 3FO4 A C74U → consensus 80 (Engineered mutation).
- 4LX6 A U74C → consensus 80.
- 9IWG A/B: 7 rows each.
- 2B57 A G16 added → consensus 20.
