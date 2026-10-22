# Reported sequence differences: validation (2026-10-22)

**Scope.** A new "Reported sequence differences" category in the Annotations tab. It shows the sequence differences that each PDB entry records for the family RNA chains, read per chain from the PDBe updated mmCIF (`_struct_ref_seq_dif`) and mapped to Rfam positions.

**Related documents:**
- Method, written so it can be rebuilt at PDBe: `docs/SEQUENCE_DIFFERENCES_METHOD.md`.
- Investigation, including why the PDBe API cannot be used yet: `experiments/mutations_investigation/REPORT.md`.

**Why not the API.** The PDBe API `mutated_AA_or_NA` covers proteins only at present; its own description says so. It returned no RNA record for any family.

## 1. Entries tested and classifications observed

All 289 entries were read (streamed for RF00005, RF00011 and RF00177, never stored). Classifications are kept exactly as deposited.

| Family | Entries read | RNA rows | On Rfam positions | Off-axis | Unresolved | Copies disagreeing | Classifications (rows) | Rows on other polymers (not shown) |
|---|---|---|---|---|---|---|---|---|
| RF00167 | 66/66 | 85 | 77 | 8 | 0 | 0 | variant 40, Engineered mutation 35, Conflict 6, SEE REMARK 999 4 | – |
| RF00162 | 28/28 | 2 | 0 | 2 | 0 | 0 | deletion 1, Conflict 1 | polypeptide 2 |
| RF00059 | 28/28 | 8 | 2 | 6 | 0 | 0 | Insertion 8 | polypeptide 8 |
| RF00005 | 59/59 | 24 | 17 | 7 | 0 | 0 | modified residue 14, Conflict 4, Insertion 3, Engineered mutation 2, Expression tag 1 | polypeptide 412, polyribonucleotide 16 |
| RF00011 | 50/50 | 73 | 32 | 41 | 0 | 0 | Conflict 46, Engineered mutation 24, Expression tag 3 | polypeptide 32, polyribonucleotide 40 |
| RF00177 | 58/58 | 0 | 0 | 0 | 0 | 0 | – | polypeptide 99, polyribonucleotide 36 |

**Examples, checked against the source rows** (reference → deposited, author numbering, Rfam position):

| Entry | Family | Chain | Change | Classification | Rfam position | Note |
|---|---|---|---|---|---|---|
| 3FO4 | RF00167 | A | C74U | Engineered mutation | 80 | The guanine→adenine specificity switch. |
| 2B57 | RF00167 | A | C74U | Engineered mutation | 80 | |
| 2B57 | RF00167 | A | G16 added | Engineered mutation | 20 | |
| 4LX6 | RF00167 | A | U74C | Engineered mutation | 80 | Reverse direction to 3FO4 and 2B57. |
| 4LX6 | RF00167 | A | U28G, G42C, U47C, U51C | Engineered mutation | 32, 46, 51, 55 | |
| 6VWT, 6VWV | RF00167 | A and B | 10 rows each | variant | 19–22, 84–88 | P1 stem redesign; each copy separately. |
| 5C7U | RF00167 | B | A15G, C16G, U17A, G79U, G80C, U81C | Conflict | 19–21, 85–87 | Construct P1, labelled as conflict. |
| 9IWF / 9IWG | RF00167 | A and B | 4 and 7 rows per chain | Engineered mutation | | |
| 1EHZ | RF00005 | A | G10→2MG, U16→H2U, … | modified residue | | tRNA modifications recorded as differences; also in the Nonstandard nucleotides track, legitimately. |
| 8B0X | RF00005 | Z | C34U, U36C | Conflict | 33, 35 | Anticodon loop. |
| 3L0U | RF00005 | A | C3G, G70C | Engineered mutation | 3, 68 | |
| 9OY2, 9OY3, 9ZDV, … | RF00011 | A | C126U | Conflict | 113 | |
| 4KQY | RF00162 | A | C63A | Conflict | – | At an insertion. |
| 4KQY | RF00162 | A | U absent | deletion | – | No deposited residue. |
| 8F4O | RF00059 | A and B | 4 added residues each | Insertion | 1 at 93, 3 in flanks | |

## 2. API fields and limitations

The mmCIF items and their API equivalents are listed in the method document.

**Limitations of the data, not of the processing:**
- **Labels are free text.** Case and wording vary; for example "variant" vs "Variant", and "modified residue" on tRNA.
- **Differences are relative to the depositor's reference sequence** (here GenBank, or the entry itself), not to the Rfam consensus.
- **Labels mix construct design with mutations:** P1 redesign appears as "Conflict" or "variant", and added transcription nucleotides as "Engineered mutation", "Insertion" or "Expression tag".
- **Most entries record nothing.** For example, 112 of 122 RNA entities in RF00167, RF00162 and RF00059 reference their own PDB entry. So an empty position is never shown as "unchanged".

## 3. Mapping

- **Resolution:** every RNA row resolved to an entity and to the label chain of its own strand. None was unresolved, and none was moved to another chain.
- **Identity and numbering:** every row with a deposited residue agrees with `_pdbx_poly_seq_scheme` at its `label_seq_id` and with our validated mapping, on residue identity and on author number with insertion code. There are no conflicts.
- **Rows without an Rfam position:** 64 rows across the families. They are listed in the data and in the off-axis note under the tracks:
  - non-submitted flanks: 61;
  - insertions: 2 (4KQY C63A; one in 1EHZ);
  - absent residue: 1 (4KQY).
- **Copies:** every multi-copy entity reports the same rows for each copy.

## 4. Implementation

**Backend.** `rna_explorer/sequence_differences.py`, run with `python3 -m rna_explorer seqdiff <ACC>`.
- Part A, `extract_chain_differences`: the per-chain extraction from one mmCIF, independent of Rfam.
- Part B: placement through the chain's validated mapping (`place_on_family`), aggregation, and the copy-consistency check.
- Outputs: `sequence_differences.json` and `sequence_difference_tracks.json`, with new schemas.
- Unchanged: `ligands.json`, `tracks.json`, `modifications.json` and every other output.

**Frontend.** "Reported sequence differences" category between RNA loop annotations and Nonstandard nucleotides (`seqDiffTracks`).
- Tracks: one per classification found on the axis, named verbatim and ordered by RNA instances. Each has a fixed colour per label and a presence legend; there is no frequency.
- Tooltips: classification, changes as reference→deposited, RNA instances and entries, and up to three observations ("3FO4 A C74U").
- Pinned details: every observation as a button that opens that RNA chain in the Family Overview and highlights the position there. The Overview note says "Recorded for this RNA chain".
- Off-axis rows are listed under the categories; families with none on the axis say so.
- Two new downloads.

## 5. Validation

**pytest:** `tests/test_sequence_differences.py`, 4 tests.
- A synthetic mmCIF covers: substitution on two copies, added residue, absent residue, identity conflict, unresolved strand and a protein row.
- Placement statuses.
- The copy-consistency check.
- RF00167 regression: counts, 3FO4 C74U at 80, 6VWT both chains, consensus 80 = C>U ×2, U>C ×1.

**vitest (115 passed):**
- `seqDiffTracks`: tracks only for classifications on the axis, fixed colours, no mark where nothing is reported.
- Tooltip lines and per-chain observations.
- Labels for substitution, added and absent residues, insertion codes and modified residues.

**Playwright:**
- `annotations-seqdiff.spec.ts` (4 tests):
  - category order and tracks equal the data file;
  - legend and info text;
  - off-axis note;
  - the C74U tooltip and pinned observations, with "4LX6 A" opening chain 4LX6 A at C74;
  - two copies as two instances (6VWT B opens chain B);
  - downloads;
  - numbering changes labels only;
  - SAM-I shows its "none at Rfam positions" note.
- Full suite: 208 passed and 1 failed (a SAM-I Structures-tab 3D test); it passed when re-run, together with the Ligands-tab spec (14 passed).
- pytest total: 111 passed.
- Updated for the new category: `annotation-tracks`, `annotations-ui`, `annotations-download` and `modification-track` (category, legend and download counts).

## 6. Limitations and possible improvements

- **Switch to PDBe's own data when available.** Once `mutated_AA_or_NA` covers nucleic acids, the source can be swapped and checked against these 192 rows.
- **Optional curation.** A per-entry note could separate functional mutations from construct design (option C in the investigation).
- **Absent residues have no position on the axis.** They are listed only.
- **The 16S family records no differences on its chains.** The category is hidden there and a note says so.
