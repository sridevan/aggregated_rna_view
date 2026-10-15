# BGSU RNA 3D Hub correspondence comparison for RF00177 mappings

Date: 2026-10-14. Inputs: BGSU `correspondence/align_chains` (pairwise and 10-chain, with `column`), BGSU `rna3dhub/rest/SeqtoUnitMapping`, our RF00177 mappings (`poc_rf00177/results/`, Infernal 1.1.5 cmalign local mode, Rfam 15.1 RF00177 CM, CKSUM 2093381478). Original RF00167 and RF00177 results were not modified; the documented pair (5J7L, 1J5E) was aligned with the same scripts into `bgsu_validation/results/repro/`.

## Step 1: BGSU response format

- One nucleotide per line, tab-separated unit ids, `NULL` where a chain has no observed nucleotide in that alignment column. Unit id = `PDB|model|chain|comp|number[|atom|alt|insertion_code|symmetry]`; `comp` is the mmCIF component id (3-letter for modified residues), `number` is the author residue number, and insertion codes appear in field 8 (e.g. `1J5E|1|A|C|190|||C`).
- With `column` as the first pseudo-chain, the first field is the alignment column number. Columns run to 2908 for 16S, are shared by all pairwise and multi-chain requests (13,539 residues seen in several files, 0 column inconsistencies), i.e. BGSU serves slices of one fixed family-wide multiple alignment.
- BGSU columns are not Rfam consensus positions numerically, but the 1533 columns that hold any of our match-state residues map one-to-one and monotonically onto the 1533 RF00177 consensus positions, with no conflicts and no column mixing match with insert/local-end residues (`bgsu_column_to_consensus.csv`). The other columns are insert columns. This is consistent with BGSU aligning every chain with cmalign to the same RF00177 model.
- Only observed nucleotides have unit ids; unobserved residues are `NULL`, so BGSU correspondences exist only between residues with coordinates.
- `SeqtoUnitMapping` lists every experimental sequence position with its unit id or `NULL`; its sequence position equals mmCIF `label_seq_id` for all ten chains tested (12,194 observed units, 0 differences).

## Identifier normalisation

Every BGSU unit id in all files resolved to exactly one `_pdbx_poly_seq_scheme` row through (pdb_strand_id, pdb_seq_num, pdb_ins_code), with the BGSU component equal to `mon_id` (modified residues included): 0 unresolved, 0 ambiguous, 0 component mismatches across 17,171 pairwise correspondences and the 10-chain alignment. Differences attributable to identifier mapping: none.

## Comparison dataset (all available through the API)

| Pair | Purpose |
|---|---|
| 5J7L/AA vs 1J5E/A | documented BGSU example, reproduced with our pipeline (E. coli vs Thermus) |
| 7K00/A vs 4Y4O/1a, 5J7L/AA, 1J5E/A | bacterial; 5J7L is a second E. coli structure with identical sequence |
| 7K00/A vs 6ERI/BA | bacterial vs chloroplast |
| 7K00/A vs 7P2E/A, 5MRC/aa, 6YW5/aa, 9GYT/2, 8APN/B3 | bacterial vs mitochondrial (human 12S, yeast 15S, Neurospora, Brassica, Polytomella fragment) |
| 5MRC vs 6YW5, 7P2E vs 9GYT, 6ERI vs 9GYT, 9GYT vs 6YW5 | organellar vs organellar, to probe correspondences inside our unnumbered regions |

Our correspondence for a pair = two residues that both occupy match states at the same consensus position. Observed subset = both residues have coordinates (the only domain BGSU can report).

## Quantitative results (`pair_metrics.csv`, `comparison_summary.json`)

| Pair | BGSU pairs | ours (all) | ours (observed) | matching | precision | recall | ours-only (observed) | BGSU-only | BGSU-only category |
|---|---|---|---|---|---|---|---|---|---|
| 5J7L vs 1J5E | 1484 | 1488 | 1484 | 1484 | 1.000 | 1.000 | 0 | 0 | |
| 7K00 vs 1J5E | 1479 | 1496 | 1479 | 1479 | 1.000 | 1.000 | 0 | 0 | |
| 7K00 vs 4Y4O | 1473 | 1496 | 1473 | 1473 | 1.000 | 1.000 | 0 | 0 | |
| 7K00 vs 5J7L | 1519 | 1508 | 1500 | 1500 | 1.000 | 0.988 | 0 | 19 | insert/insert, identical E. coli sequences |
| 7K00 vs 6ERI | 1469 | 1476 | 1469 | 1469 | 1.000 | 1.000 | 0 | 0 | |
| 7K00 vs 7P2E | 694 | 694 | 694 | 694 | 1.000 | 1.000 | 0 | 0 | |
| 7K00 vs 5MRC | 1283 | 1331 | 1279 | 1279 | 1.000 | 0.997 | 0 | 4 | insert/insert at anchor 208 |
| 7K00 vs 6YW5 | 1295 | 1339 | 1295 | 1295 | 1.000 | 1.000 | 0 | 0 | |
| 7K00 vs 9GYT | 1391 | 1398 | 1391 | 1391 | 1.000 | 1.000 | 0 | 0 | |
| 7K00 vs 8APN | 367 | 367 | 367 | 367 | 1.000 | 1.000 | 0 | 0 | |
| 5MRC vs 6YW5 | 1262 | 1279 | 1220 | 1220 | 1.000 | 0.967 | 0 | 42 | 40 local-end/local-end, 2 insert/insert |
| 7P2E vs 9GYT | 712 | 694 | 694 | 694 | 1.000 | 0.975 | 0 | 18 | local-end/local-end at anchor 1444 |
| 6ERI vs 9GYT | 1388 | 1388 | 1388 | 1388 | 1.000 | 1.000 | 0 | 0 | |
| 9GYT vs 6YW5 | 1355 | 1310 | 1286 | 1286 | 1.000 | 0.949 | 0 | 69 | 68 local-end/local-end, 1 insert/insert |
| **Total (14 pairs)** | **17171** | **17264** | **17019** | **17019** | **1.000** | **0.991** | **0** | **152** | 26 insert/insert, 126 local-end/local-end |

- Agreement among experimentally observed residues: 17,019 of 17,019 of our observed match-state pairs are in BGSU (100 %). Our 245 further pairs involve at least one residue without coordinates; BGSU cannot report those.
- High-confidence pairs (both residues PP "*"): 14,761, agreement 100 %. Low-confidence pairs: 2,258, agreement 100 %. Pairs inside our flagged 50-position windows (low identity, low PP, sparse coverage, long insertion, pairing inconsistency): 2,785, agreement 100 %. Region-by-region (`region_windows.csv`, 31 windows per pair): no window with agreement below 1.0.
- Not a single shift (match-match pair with different consensus positions) occurred in any pair. Differences attributable to Infernal settings or Rfam versions: none detectable; the exact equality of pair counts for 7P2E (694) and 8APN (367), the one-to-one column-to-consensus mapping, and the absence of shifts imply BGSU uses the same RF00177 model length (1533) and a local alignment mode equivalent to ours. BGSU does not publish its Infernal/Rfam versions in the response; this remains an assumption.

## Step 4: unnumbered regions (`unnumbered_regions_bgsu.csv`, `unnumbered_bgsu_pairs_detail.csv`)

Our five states were kept distinct: ordinary insert residues (anchor+offset), local-end residues (anchor+offset, separate state), deleted consensus positions (inside the aligned span), uncovered consensus positions (outside the span), and residues without coordinates (independent flag).

- Bacterial reference vs organellar: BGSU gives **no** partner for any of our local-end residues (7P2E 249, 5MRC 248, 6YW5 421, 9GYT 353, 8APN 10) and none for organellar insert residues, except 4 residues of 5MRC at anchor 208 paired with 7K00 insert residues at the same anchor and offset (0 of 4 nucleotides identical).
- Identical sequences: 7K00 vs 5J7L (both E. coli K-12 rrnB) get 19 insert/insert pairs from BGSU, same anchor and same offset, 19/19 identical nucleotides. These are correct by sequence identity, but they arise from column coincidence, not from the model.
- Organellar vs organellar: BGSU pairs 40 (5MRC-6YW5), 18 (7P2E-9GYT) and 68 (9GYT-6YW5) residues that are local-end in both chains, all at the same anchor (1444, 420, 82) and in every case either the same offset counted from the left or the same offset counted from the right (exactly half each for 7P2E-9GYT and 9GYT-6YW5), with nucleotide identity 45 %, 28 % and 21 %. This is the signature of Infernal's split insert placement in a multiple alignment (first half left-justified, second half right-justified) rather than a model-based alignment: the run lengths differ (e.g. 18 vs 71 at anchor 1444, 37 vs 59, 66 vs 113), the CM has no states for these residues, and BGSU's own documentation says cross-domain alignments "may not be as reliable". No different reference model or alternative alignment strategy is involved; it is the same cmalign output read column-wise.
- Deleted and uncovered consensus positions have no residue and therefore no BGSU unit on either side; BGSU represents them exactly as we do (NULL / absent).
- Residues without coordinates never appear in BGSU; our mapping keeps them with a consensus position and `has_coordinates = 0`, which is the correct superset for a sequence-level view and must be filtered for a coordinate-level view.

## Step 6: numbering scheme review

The scheme (Rfam consensus position for match states; anchor+offset for ordinary insertions; explicit unnumbered local-end residues with anchor+offset; preserved polymer and author numbering; alignment provenance) reproduces BGSU's correspondences completely wherever BGSU and the CM agree there is a correspondence, and it refrains from asserting the 152 column-coincidence pairs that BGSU reports inside insert and local-end regions. No additional family-wide numbering scheme is supported by these results. Two additions would help motif aggregation without changing numbering:

1. A per-residue `bgsu_column` (or equivalently the full alignment column of our own multiple alignment) so that two residues in the same insert column can be flagged as "co-columnar" without being called equivalent. This is the only information BGSU provides beyond our scheme, and it is weak.
2. A per-pair confidence attribute when two chains are compared: both match + PP, both match + flagged window, co-columnar insert, co-columnar local end.

Entity-level alignment for duplicated chains (5J7L BA, 4Y4O 2a) remains sufficient: BGSU's within-entry alignments are identical for both copies.

## Answers

1. Yes. For conventional bacterial 16S (5J7L/1J5E documented example, 7K00 vs 4Y4O/1J5E/5J7L) every one of our observed correspondences is in BGSU and every BGSU correspondence is ours, except the 19 insert-residue pairs between the two identical E. coli sequences.
2. Bacterial vs organellar pairs do not disagree at all (precision and recall 1.0 for 6ERI, 7P2E, 6YW5, 9GYT, 8APN; 4 insert pairs in 5MRC). Organellar vs organellar pairs differ only inside shared local-end regions (anchors 82, 420, 1444), where BGSU reports 126 column-coincidence pairs.
3. No. BGSU's extra correspondences in our unnumbered regions are same-anchor, left- or right-justified offset coincidences with near-random nucleotide identity, not model-supported alignments.
4. Yes. Local-end residues should remain unnumbered in the viewer; they can be shown as "unaligned region after consensus N (k residues)" and, if desired, co-columnar residues in other structures can be listed as weak candidates, clearly labelled.
5. Yes for consensus-numbered residues: identifier mapping is exact, correspondence agrees with an independent implementation at 100 % on observed residues, and confidence flags are available per residue and per window. Aggregation must filter on coordinates for structure-derived annotations (ligand contacts, motifs) and must treat insert/local-end residues as per-structure, non-aggregatable unless explicitly co-columnar.

## Caveats

- BGSU is treated as an independent implementation of the same method, not as ground truth; the identical model length and alignment behaviour mean the two implementations share Infernal's systematic errors.
- BGSU Infernal/Rfam versions are not exposed; agreement shows they are compatible for these chains but could diverge after an Rfam seed update.
- Only ten chains, one family. Families with more variable seeds (group II introns, RNase P) or where BGSU uses a different family (cross-domain SSU) were not tested.
