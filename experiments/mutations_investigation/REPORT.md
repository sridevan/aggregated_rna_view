# RNA mutation annotations: investigation (2026-10-22)

**Outcome: stopped before implementation.**
- The requested source, PDBe `mutated_AA_or_NA`, contains **no RNA records** for any of the six families.
- The underlying mmCIF data do contain RNA sequence differences, and they map perfectly to our chains and Rfam positions.
- However, their classifications are free text, often describe construct design rather than mutations, and are recorded for only a few entries. Absence is uninformative.
- A "Mutations" track built on them needs a decision on scope and wording first (§6).

No application code or data were changed.

**Files (all in this folder):**
- `fetch.py`: queries the API; raw responses in `raw/`.
- `analyse.py` → `records_classified.json`: each API record against the family chains.
- `mmcif_seqdif.py` → `mmcif_seqdif_family_rna.json`: `_struct_ref_seq_dif` rows of family RNA entities in the cached PDBe updated mmCIF, mapped through the validated family mappings.

## 1. PDBe `mutated_AA_or_NA` (POST, all entries of the six families, 2026-10-22)

| Family | Entries | Entries with records | Records | On a family RNA chain |
|---|---|---|---|---|
| RF00167 | 66 | 0 | 0 | 0 |
| RF00162 | 28 | 1 | 2 | 0 |
| RF00059 | 28 | 1 | 8 | 0 |
| RF00005 | 59 | 23 | 420 | 0 |
| RF00011 | 50 | 16 | 32 | 0 |
| RF00177 | 58 | 39 | 178 | 0 |

- **All 640 records are on other chains (proteins).** None has a nucleotide `chem_comp_id`.
- **Classifications seen** (`mutation_details.type`):

  | Classification | Records |
  |---|---|
  | Expression tag | 375 |
  | Conflict | 185 |
  | Linker | 30 |
  | Modified residue | 13 |
  | Engineered mutation | 9 |
  | See remark 999 | 7 |
  | Variant | 7 |
  | Initiating methionine | 5 |
  | Insertion | 5 |
  | Acetylation | 4 |

- **Direct check.** Entries whose mmCIF has RNA sequence-difference rows (3FO4, 4LX6, 9IWG) return HTTP 404 "Requested endpoint does not contain any data".
- **Conclusion.** The endpoint does not cover nucleic-acid chains in these data, so the track cannot be built from it.
- **Confirmed by PDBe's endpoint description** (provided by the user, 2026-10-22): "This call provides a list of modelled instances of mutated amino acids in proteins in a given set of PDB entries. (Note that at present it does not provide information about mutated nucleotides in RNA or DNA chains, but it would do so in near future.)" Nucleic-acid chains are excluded by design for now.

## 2. The mmCIF source: `_struct_ref_seq_dif` (RF00167, RF00162, RF00059; cached PDBe updated mmCIF)

**What the rows are.** Each row describes a residue that differs between the deposited sequence (`mon_id`) and the entity's reference sequence (`_struct_ref`, `db_mon_id`). For the family RNA entities:
- 122 of 122 have a `_struct_ref`: 112 refer to the PDB entry itself, 10 to GenBank.
- Differences can only be recorded against an external reference. The 112 self-referenced entities therefore carry no rows, whatever their construct.

**Volume:** 95 rows on family RNA entities in 11 entries (RF00167 85, RF00059 8, RF00162 2).

| Entry | Family | Reference | Chains | Rows per chain | Classification (`details`, verbatim) | What they are |
|---|---|---|---|---|---|---|
| 3FO4 | RF00167 | PDB:3FO4 | A | 1 | Engineered mutation | C74U (consensus 80): the guanine→adenine specificity switch |
| 2B57 | RF00167 | GB:BACYACA | A | 3 | Engineered mutation | C74U (cons. 80); two added 5′ G/A (cons. 20–21; reference residue "–") |
| 4LX5 | RF00167 | GB:AE016796 | A | 6 | Engineered mutation ×4, SEE REMARK 999 ×2 | substitutions at cons. 32, 46, 51, 55; added 5′ G and 3′ C (flanks) |
| 4LX6 | RF00167 | GB:AE016796 | A | 7 | Engineered mutation ×5, SEE REMARK 999 ×2 | as 4LX5 plus C74U (cons. 80) |
| 9IWF | RF00167 | GB:CP011058.1 | A, B | 4 + 4 | Engineered mutation | P1 stem ends (cons. 20–21, 85–86) |
| 9IWG | RF00167 | GB:CP011058.1 | A, B | 7 + 7 | Engineered mutation | P1 ends plus cons. 40, 44, 67 |
| 5C7U | RF00167 | GB:X83878.1 | A | 6 | Conflict | P1 stem ends (cons. 19–21, 85–87): construct design reported as conflict |
| 6VWT, 6VWV | RF00167 | GB:CP009262.1 | A, B | 10 + 10 each | variant | P1 stem redesign (cons. 19–22, 84–88, one flank) |
| 8F4O | RF00059 | GB:CP026935.2 | A, B | 4 + 4 | Insertion | added residues (reference "–"): 5′ G, 3′ GCA; one at cons. 93 |
| 4KQY | RF00162 | GB:CP003695 | A | 2 | Conflict ×1, deletion ×1 | C→A at an insert (no consensus position); one reference U absent from the construct (no residue) |

**Mapping (validated family mappings):**
- All 94 rows with a deposited residue match our mapping in both the nucleotide (`mon_id` = mapping `residue_name`) and the author number with insertion code.
- 79 map to a consensus position.
- 15 have none: 14 in non-submitted flanks, 1 at an insertion (4KQY).
- 1 deletion has no residue at all.
- Multi-chain entities have a row for every chain (6VWT, 6VWV, 9IWF, 9IWG, 8F4O: A and B), so chain-specific handling works.
- No ambiguous mapping was found.

## 3. Classification issues (why not to implement as specified)

1. **The source is free text.** `details` takes depositor or annotator wording: "Engineered mutation", "Conflict", "variant", "Insertion", "deletion", "SEE REMARK 999" (inconsistent case). It is not a controlled vocabulary.
2. **The labels don't separate experimental mutations from construct design.**
   - 5C7U's "Conflict" and 6VWT/6VWV's "variant" rows are a redesigned P1 helix.
   - 9IWF's "Engineered mutation" rows are the same kind of P1 change.
   - Only 3FO4, 2B57, 4LX6 (C74U) and the 4LX5/4LX6 core substitutions are mutations in the functional sense.
3. **Insertions and deletions are mixed with substitutions.**
   - In 2B57, added 5′ nucleotides are "Engineered mutation" with no reference nucleotide.
   - 8F4O has "Insertion" rows, and 4KQY has a "deletion" with no deposited residue.
   - A from→to display needs to show these as "added" and "absent", not as substitutions.
4. **Differences are relative to the reference sequence chosen by the depositor,** often a genomic GenBank region. They are not relative to the Rfam consensus or a defined wild type, and they can differ between entries of the same RNA.
5. **Coverage is sparse and absence is uninformative.**
   - Only 11 of 122 family entries have any row.
   - Most engineered RNA constructs in these families reference themselves (PDB) and record nothing, e.g. many C74U guanine-riboswitch structures.
   - An empty track position must never read as "no mutation".

## 4. API fields (for reference)

**Fields:** `entity_id`, `chain_id` (author), `struct_asym_id` (label), `author_residue_number`, `author_insertion_code`, `residue_number` (= label_seq_id), `chem_comp_id`, `mutation_details.{from,to,type}`.

**mmCIF equivalents in `_struct_ref_seq_dif`:**
- `pdbx_pdb_strand_id` (author chain);
- `seq_num` (label_seq_id);
- `pdbx_auth_seq_num`, `pdbx_pdb_ins_code`;
- `mon_id` (deposited), `db_mon_id` (reference);
- `details` (classification).

The label chain comes from `_pdbx_poly_seq_scheme` (entity + strand → asym), and the entity from `_struct_ref_seq.ref_id` → `_struct_ref.entity_id`.

## 5. Not checked

- tRNA, 16S and RNase P mmCIFs are not cached and were not downloaded (16S files are large). Their API results already show no RNA records.
- Whether PDBe plans nucleic-acid coverage for `mutated_AA_or_NA`.

## 6. Options

| Option | Implementation |
|---|---|
| A. Do not add a Mutations track yet | None |
| B. Add "Reported sequence differences" from mmCIF | New backend step like `modifications`: PDBe updated mmCIF `_struct_ref_seq_dif` per family chain, mapped through the validated mappings; one subtrack per verbatim `details` value. |
| C. As B, curated | The same rows with a reviewed per-entry label (functional mutation vs construct design), recorded with its rationale. |

**Option B** needs these rules:
- Substitutions, added residues and absent residues are shown separately.
- Markers are counted per RNA instance, with chain-specific highlighting.
- Off-axis rows (flanks, insertions, deletions) are listed, not dropped.
- The copy says differences are relative to the depositor's reference sequence and that absence is not a negative.
- Tracks are named after the source wording, not "mutations".

**Option C** is more meaningful for users, but it is manual curation that has to be repeated for every family.

**Recommendation:** B, with the "Reported sequence differences" naming and caveats, if the track is wanted now; A otherwise. Not as specified, because the specified API has no RNA data.
