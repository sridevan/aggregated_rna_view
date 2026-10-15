# RF00177 (bacterial SSU rRNA) consensus residue mapping: validation report

Date: 2026-10-14. Follows the RF00167 proof of concept (`poc_rf00167/results/VALIDATION_REPORT.md`); scripts reuse its design (same Stockholm parsing semantics, same mmCIF checks) with a generalised mmCIF reader and an explicit insertion representation.

## Provenance

| Item | Value |
|---|---|
| Rfam release | 15.1 (FTP README; CM header DATE 2024-07-10) |
| Covariance model | RF00177 SSU_rRNA_bacteria, CLEN 1533, W 1587, STATES 4758, NODES 1197, CKSUM 2093381478, built with Infernal 1.1.4 (cmbuild), seed 99 sequences (Gutell/Nawrocki, PMID 11869452) |
| Infernal | 1.1.5 (Sep 2023), built from source; `cmalign` default local mode, `--cpu 4`; global mode (`-g`) run separately for the stability test |
| PDBe inputs | `api/v2/nucleic_mappings/RF00177`, `pdb/entry/molecules` (POST batch), `pdb/entry/residue_listing`, updated mmCIF files |
| Rfam inputs | `family/RF00177/cm`, `family/RF00177/structures` (JSON), `family/RF00177/alignment/stockholm` (seed), `database_files/pdb_full_region.txt.gz`, structures endpoints of the other SSU-clan families (RF01959, RF01960, RF02542, RF02545) |
| PDBe search API | `search/pdb/select?q=rfam_accession:RF00177` for organism and taxonomy per entity (the molecules endpoint returns an empty `source` for most ribosome entries) |

## Step 1: RF00177 mappings in PDBe

| Metric | Value |
|---|---|
| PDB entries | 1473 |
| Mapped RNA instances (chains) | 1817 |
| Unique (entry, entity) pairs | 1494 |
| Unique entity sequences (md5) | 233 |
| Unique mapped-region sequences | 210 |
| Mapped length min / median / max | 60 / 1519 / 1808 nt |
| Length bins | 0-199: 10; 300-499: 10; 900-1099: 107 (mammalian mito 12S); 1300-1499: 18; 1500-1599: 1650; 1600-1899: 22 |
| Distinct organisms | 58 (T. thermophilus HB8 620, E. coli 458 + K-12 195, H. sapiens 90, M. pneumoniae 51, ...) |
| Entities with >1 mapped chain | 316 (e.g. two 70S copies per asymmetric unit); the alignment is identical for all copies of an entity |
| Instances with modified nucleotides in the entity | 769 |
| Author insertion codes at segment boundaries | 0 (but present inside 4y4o, see below) |
| Rfam structures endpoint vs PDBe list | identical: same 1473 entries, same 1817 chains, same segment boundaries |
| Rfam bit score distribution (structures endpoint) | <100: 11; 100-300: 123; 300-600: 11; 600-1200: 4; >=1200: 1668 |

Taxonomy question raised during the run: membership in a "bacterial" family is decided by score against the model, not by taxonomy. Organellar rRNAs (mitochondrial 12S/15S/16S/18S-named, chloroplast 16S) are all RF00177 by score. There are no archaeal chains in the list. There are four eukaryotic cytosolic 18S chains (9b0q, 9b0s; human) that RF00177 also hits (456 bits) while RF01960 (eukaryotic SSU) scores 1669; these are the only instances where RF00177 is not the best family within the SSU clan (CL00111: RF00177, RF01959, RF01960, RF02542, RF02545). The Rfam FTP table `pdb_full_region.txt` is not clan-competed at all (almost every bacterial chain also hits RF01959/RF01960/RF02542); the family `structures` endpoint is clan-competed except for these four. Column `rf00177_is_clan_best` in `rf00177_instances.csv` records this. 4v4b (yeast, molecule name "18S Ribosomal RNA") scores 1479 bits to RF00177 and 633 to RF01960: its deposited sequence is bacterial 16S despite the name.

## Step 2: selected structures and rationale

| PDB/chain | Entity | label_asym | Organism (mmCIF) | Chain len | Mapped region (polymer) | Author numbering | Rfam bits | Why |
|---|---|---|---|---|---|---|---|---|
| 7k00/A | 1 | A | Escherichia coli | 1542 | 1-1542 | 1-1542 | 1578 | Reference bacterial 16S at 2.0 A; E. coli numbering is the community frame; 11 modelled modified nucleotides |
| 4y4o/1a | 32 | FA | Thermus thermophilus HB8 | 1521 | 1-1519 | 1-1542 with insertion codes | 1508 | Second classic bacterium (X-ray 2.3 A); author numbering follows E. coli with gaps and insertion codes 189A-189I, 1030A-D |
| 6eri/BA | 34 | HA | Spinacia oleracea (chloroplast) | 1481 | 2-1481 | 3-1482 | 1517 | Chloroplast 16S; RF00177 mapping confirmed (model 1-1523) |
| 7p2e/A | 1 | A | Homo sapiens (mito 12S) | 955 | 4-949 | 651-1596 (mtDNA coordinates) | 251 | Mitochondrial SSU explicitly mapped to RF00177; reduced rRNA, many deletions |
| 5mrc/aa | 75 | WB | Saccharomyces cerevisiae (mito 15S) | 1649 | 8-1634 | 8-1634 | 260 | Divergent, expanded mito rRNA; large insertions, lowest full-length confidence |
| 6yw5/aa | 35 | KA | Neurospora crassa (mito) | 1864 | 46-1853 | 46-1853 | 597 | Longest mapped region in the family (1808 nt) |
| 9gyt/2 | 87 | IC | Brassica oleracea (mito 18S) | 1766 | 7-1766 | 9-1848 | 1232 | Plant mito rRNA: long insertions in a confidently aligned background; 2.11 A; 6 modified nucleotides |
| 8apn/B3 | 62 | JB | Polytomella magna (fragmented mtSSU-3) | 379 | 1-379 | 13-391 | 126 | Fragmented rRNA piece covering model 488-980 only: partial coverage, both ends truncated |

Mitochondrial assignments were verified per entry (`nucleic_mappings/{pdb}`): human, pig, bovine, yeast, Neurospora, Brassica, Arabidopsis, Polytomella and Chlamydomonas mito SSU rRNAs are RF00177; Trypanosoma mito SSU (6hiv) is RF02545 and was not selected. Not aligned but noted: 6ydp/6ydw (pig) deposit one 16,770-nt entity that carries 12S, 16S and nine tRNA hits on a single chain.

## Step 3: sequences and normalisation

PDBe's `sequence` field already gives parent bases for modified residues (`pdb_sequence` shows e.g. `(PSU)`); 34 modified residues lie inside the mapped regions (7k00 11, 4y4o 12, 9gyt 6, 7p2e 5) and all were aligned as their parent base to a consensus match column. No non-ACGU letters occurred. The 3-letter code is retained per residue (`residue_name`, `is_modified`). Full chain sequences are kept in `intermediate/rf00177_full_sequences.fasta`; cmalign input is the PDBe-mapped region only (RF00167 lesson).

## Step 4: alignment summary (local mode)

| Chain | cmalign bits | avg PP | trunc | model span | match | insert | local-end | deleted in span | uncovered | model coverage | PP "*" fraction | residues PP<0.7 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7k00/A | 1571.5 | 0.990 | no | 1-1533 | 1516 | 26 | 0 | 17 | 0 | 98.9 % | 0.958 | 13 |
| 4y4o/1a | 1576.6 | 0.986 | no | 1-1533 | 1512 | 7 | 0 | 21 | 0 | 98.6 % | 0.932 | 14 |
| 6eri/BA | 1535.0 | 0.986 | 3' | 1-1523 | 1478 | 2 | 0 | 45 | 10 | 96.4 % | 0.958 | 26 |
| 7p2e/A | 264.8 | 0.887 | 3' | 9-1522 | 694 | 3 | 249 | 820 | 19 | 45.3 % | 0.630 | 143 |
| 5mrc/aa | 647.6 | 0.858 | no | 1-1533 | 1346 | 33 | 248 | 187 | 0 | 87.8 % | 0.658 | 344 |
| 6yw5/aa | 650.4 | 0.911 | no | 1-1533 | 1355 | 32 | 421 | 178 | 0 | 88.4 % | 0.768 | 214 |
| 9gyt/2 | 1233.3 | 0.972 | 3' | 1-1528 | 1399 | 8 | 353 | 129 | 5 | 91.3 % | 0.896 | 49 |
| 8apn/B3 | 126.8 | 0.922 | 5'&3' | 488-980 | 367 | 2 | 10 | 126 | 1040 | 23.9 % | 0.728 | 31 |

Runtime: 2.1 s wall for all eight sequences (`--cpu 4`), peak memory 210 MB. The alignment has 2633 columns: 1533 match, 84 insert (`.`), 1016 local-end (`~`). Model spans agree with Rfam's structures endpoint for 7 of 8 chains; 9gyt ends at 1528 in cmalign versus 1533 in Rfam's cmsearch hit (different algorithm; the last five consensus positions are "uncovered" here). cmalign and cmsearch bit scores differ (e.g. 5mrc 648 vs 260) because cmalign scores the forced alignment of the full region in a different mode; per-chain confidence should be taken from one consistently chosen source.

Every residue of every mapped region received exactly one state: match (consensus position), insert (anchor + offset), or local-end (anchor + offset), and every consensus position appears exactly once per chain in `rf00177_consensus_columns.csv` as match, deletion (inside the aligned span) or uncovered (outside the span, i.e. truncation or fragment).

## Step 5: insertion numbering

Representation used in `rf00177_residue_mapping.csv`: `alignment_state` in {match, insert, local_end, not_submitted_flank}; `rfam_consensus_pos` only for match; `insertion_anchor` = consensus position of the nearest match column to the left in the RF line (0 before the first); `insertion_offset` = 1..n within the run. Author insertion codes are a separate column and never interact with this.

- Synthetic test (`intermediate/anchor0_test.*`): 7k00 region with 5 nt prepended and 6 nt appended. Result: prefix residues anchor 0, offsets 1-5; suffix residues anchor 1533, offsets 1-6; all 1516 core match residues keep exactly the same consensus positions as in the real run.
- Real data: 141 insertion sites across the eight chains; 17 sites with >= 20 residues; longest: 225 nt after consensus 1130 in 9gyt (polymer 1102-1326, author 1121-1406, all with coordinates), 122 nt after 1022 in 6yw5 (no coordinates), 120 nt after 82 in 6yw5, 113 nt after 420 in 6yw5, 102 nt after 464 in 5mrc. Full list: `rf00177_longest_insertions.csv`.
- Anchor after deleted consensus positions: when consensus positions are deleted immediately before an insertion, the anchor is the model position preceding the insert column, which can be higher than the last *matched* consensus position of that chain (15 such sites in 7p2e, 4 in 5mrc, 6 in 6yw5, 3 each in 9gyt and 8apn). This is correct in model coordinates and the validator checks it, but consumers must not assume "anchor = previous matched residue".
- Stability around long insertions: the consensus positions flanking the five longest insertions are identical between local and global alignment modes for 9gyt and 5mrc at those sites, while global mode (`-g`) reassigns 38-55 match residues in 5mrc, 6yw5 and 7p2e elsewhere and collapses 7p2e completely (30 bits, 5' misassigned by 440 positions). Bacterial and chloroplast chains are mode-independent (0, 0 and 1 residues differ). Local mode with explicit local-end states is therefore the correct production setting; a residue in a local-end region must not be given a consensus position. Details: `local_vs_global_stability.json`.
- Even E. coli has insertions relative to the seed consensus: 26 residues at four sites (author 72-75 after consensus 71; 94-96 after 89; 201-211 after 208; 461-468 after 458), and 17 consensus positions deleted (123, 178-184, 195-201, 451, 1125). Spinach and Thermus delete the same 177-201 block. Consensus positions 177-210 are a seed-definition artefact region (variable helix) and should be treated as low reliability for cross-species aggregation.

## Step 6: residue identifier validation

All 11,157 residue rows (eight full chains, flanks included) were compared with `_pdbx_poly_seq_scheme`: entity_id, label_asym_id, auth_asym_id (pdb_strand_id), label_seq_id, pdb_seq_num, pdb_ins_code and mon_id. Result: 0 mismatches. Coordinate presence was determined from `atom_site` (any atom for the label_seq_id) and cross-checked with `auth_seq_num == '?'` in the scheme and PDBe `observed_ratio`: 0 conflicts.

Residues with and without coordinates (never treated as deletions):

| Chain | modelled / chain length | mapped residues without coordinates | runs |
|---|---|---|---|
| 7k00/A | 1519 / 1542 | 23 | 1; 205-213; 841-845; 1535-1542 |
| 4y4o/1a | 1500 / 1521 | 19 | 1-3; 78-84; 1511-1519 |
| 6eri/BA | 1481 / 1481 | 0 | |
| 7p2e/A | 955 / 955 | 0 | |
| 5mrc/aa | 1501 / 1649 | 141 | 15 runs, largest 421-477 and 564-599 (both inside local-end insertions) |
| 6yw5/aa | 1435 / 1864 | 384 | 126-236, 327-358, 563-667, 1195-1328 (the large insertions are mostly unmodelled) |
| 9gyt/2 | 1766 / 1766 | 0 | |
| 8apn/B3 | 379 / 379 | 0 | |

Identifier edge cases found: (1) 4y4o uses E. coli-style author numbering with 34 jumps and 21 residues carrying insertion codes 189A-189I, 1030A-1030D; 17 of those are consensus match residues and 4 are Infernal insertions, which demonstrates that author insertion codes and alignment insertions are independent. (2) 7p2e author numbers are mitochondrial-genome coordinates (651-1596, offset 647). (3) 9gyt author numbering has 13 jumps. (4) 6eri, 4y4o, 5mrc, 6yw5, 9gyt and 8apn have label_asym_id different from auth_asym_id (two-letter label IDs such as HA, WB, KA, IC, JB). (5) Entities with several chains: 4y4o has two 16S copies (1a, 2a), one aligned.

## Step 7: bacterial vs organellar comparison

| Chain | identity to E. coli at shared consensus positions | base pairs evaluable (SS_cons) | pair consistency (WC/GU) | flagged 50-position windows |
|---|---|---|---|---|
| 7k00/A | 1.000 | 455 | 0.952 | 0 |
| 4y4o/1a | 0.763 | 457 | 0.950 | 0 |
| 6eri/BA | 0.761 | 446 | 0.937 | 2 |
| 9gyt/2 | 0.737 | 423 | 0.924 | 6 |
| 6yw5/aa | 0.574 | 403 | 0.916 | 13 |
| 7p2e/A | 0.548 | 204 | 0.917 | 18 |
| 5mrc/aa | 0.518 | 395 | 0.922 | 17 |
| 8apn/B3 | 0.569 | 78 | 0.897 | 5 |

Correspondence checks independent of alignment score (`reference_correspondence.json`, `rf00177_window_confidence.csv`, `ecoli_7k00_numbering_to_consensus.csv`):

- Universally conserved motifs located by sequence in E. coli and read back at the same consensus positions in every chain: 530 loop (E. coli 515-533 = consensus 505-523) is 84-100 % identical in all eight; 790 loop (785-797 = 775-787) 92-100 %; decoding site h44 (1494-1507 = 1485-1498) 100 % in bacteria, chloroplast and Brassica, 93 % in human 12S, 57-71 % in yeast and Neurospora (known divergence); 3'-terminal anti-Shine-Dalgarno (1530-1542 = 1521-1533) is present only in E. coli, Thermus and partially spinach, and absent (deleted/uncovered) in all mitochondrial rRNAs, as expected biologically. 8apn does not reach these positions.
- Modified nucleotides land on the same consensus positions in E. coli and Thermus (PSU516 -> 506, m7G527 -> 517, m2G966/m5C967 -> 956/957, m4Cm1402 -> 1393, m5C1407 -> 1398, m3U1498 -> 1489, m6,6A1518/1519 -> 1509/1510): an independent confirmation that the consensus frame transfers correctly between bacteria.
- Human 12S (7p2e): 820 consensus positions deleted, in eight runs >= 20 positions (consensus 62-347, 403-487, 574-747, 819-860, 986-1035, 1106-1175, 1233-1286, 1427-1456; E. coli 62-350, 406-497, 584-757, ...). These are the known missing helices of mammalian mito-rRNA, not alignment failures; 249 residues sit in local-end regions (all with coordinates) and have no consensus position.
- Low-confidence regions (window flags: sparse coverage, identity to E. coli < 0.5, mean PP < 0.8, pair consistency < 0.7, >= 20 inserted residues): none for E. coli and Thermus; 2 windows for spinach (51-100 low PP, 801-850); 6 for Brassica (dominated by its three large insertions); 13-18 for yeast, Neurospora and human 12S. Consensus 151-200 in Neurospora and 601-650 in the Polytomella fragment are flagged for pairing inconsistency as well as low identity, i.e. the alignment there is doubtful.
- Coverage of consensus positions: 254 positions have a residue in all eight chains; 60 positions only in the three bacterial/chloroplast chains; none uncovered by all.

## Step 8: insertion profile

Files: `rf00177_insertion_profile.csv` (long), `rf00177_insertion_profile_wide.csv` (anchor x instance, plus seed statistics), `rf00177_insertion_sites_summary.csv`, `rf00177_longest_insertions.csv`, `rf00177_insertion_profile.png`.

Family level (Rfam seed, 99 sequences): 112 anchors carry an insertion in at least one seed sequence, none in >= 50 % of them, longest seed insertion 13 nt. Structural observations: 141 anchors in eight chains, 17 shared by >= 2 chains, 29 anchors not seen in the seed at all, 17 sites with >= 20 residues. The organellar insertions (57-225 nt) are an order of magnitude longer than anything in the seed because the seed is bacterial; the family statistics therefore do not predict the structural observations and the two must be kept separate in any aggregated view.

## Modified nucleotides and numbering problems

- Modified nucleotides: handled through PDBe's parent-base sequence; all 34 in mapped regions aligned to match columns. Untested: modified residues whose parent is unknown (one-letter N/X) and residues PDBe lists as multiple-residue positions with microheterogeneity.
- Author numbering is non-contiguous (4y4o, 9gyt), offset (7p2e by 647, 8apn by 12, 6eri by 1), and carries insertion codes (4y4o). Polymer position (label_seq_id) remains the only safe key.
- The PDBe molecules endpoint returns an empty `source` for most ribosome entries; the search API and the mmCIF `_entity_src_nat` agree for all eight chains and should be the organism source.

## Comparison with RF00167

| Aspect | RF00167 (purine riboswitch) | RF00177 (16S) |
|---|---|---|
| Insert columns in mapped-region alignment | 0 | 84 columns, 113 insert residues |
| Local-end columns | 0 | 1016 columns, 1281 residues, in five organellar chains |
| Deletions | 2-4 per chain | 17-820 per chain |
| Truncation flags | none | 3' (6eri, 7p2e, 9gyt), 5'&3' (8apn) |
| Author insertion codes | none | 21 residues in 4y4o |
| Residues in mapped region without coordinates | 72/98 in 9nbh | up to 384 (6yw5); large insertions are often unmodelled |
| Full-chain vs region alignment | diverged | not repeated; local vs global compared instead |
| mmCIF identifier mismatches | 0 / 434 | 0 / 11,157 |
| cmalign cost | 0.01 s/seq | ~0.3 s/seq, 210 MB |

## Schema assessment

The RF00167 schema needs these additions, all already implemented in the RF00177 CSVs:

1. `alignment_state` with four values (match, insert, local_end, not_submitted_flank) instead of a free-text insertion string.
2. `insertion_anchor` (0..CLEN) and `insertion_offset` (1..n) as integer columns.
3. `has_coordinates` from atom_site, separate from `observed_ratio`, and never conflated with deletions.
4. Per-consensus-column table with three non-residue states: match, deletion (inside aligned span), uncovered (outside span: truncation or fragment).
5. Per-chain alignment provenance: cmalign mode, bit score, truncation flag, model span, avg PP, plus Rfam's own bit score/E-value from the structures endpoint, plus clan-best flag.
6. `is_modified` and `residue_name` alongside `nucleotide_parent`.
7. Entity-level alignment reused for all chains of the entity (316 multi-chain entities in this family).
8. Rfam release, CM checksum and Infernal version stored with every mapping.

## Answers to the expected-conclusion questions

1. Can the strategy handle bacterial 16S rRNA? Yes. E. coli and Thermus align with 98.6-98.9 % model coverage, mode-independent assignments, perfect identifier agreement, and conserved motifs and modified nucleotides land on the same consensus positions.
2. Are insertions correctly identified and represented? Yes. Insert (`.`) and local-end (`~`) columns are parsed from the RF line; anchor/offset bookkeeping was verified by the validator on all 141 sites and by the synthetic anchor-0/trailing test; the longest insertions (up to 225 nt) are reported with their PDB residues and coordinate status.
3. Does numbering remain reliable across divergent sequences? Partly. Conserved cores (530 loop, 790 loop, h44, modified-nucleotide sites) are consistent across all chains including mitochondrial 12S. Divergent organellar chains (bits 250-650, PP "*" fraction 0.63-0.77) have 13-18 flagged windows each where assignments are mode-dependent and should be shown with confidence, not as facts. Local-end residues have no consensus position by design.
4. Can we distinguish insertions, deletions, missing coordinates and unmapped residues? Yes, each is a distinct state in the two tables, and the validator confirms completeness and uniqueness per chain. Missing coordinates co-occur with large insertions (6yw5, 5mrc) but are recorded independently.
5. Are schema changes required before generalising? Yes, the eight items above. Additionally the pipeline must apply clan competition (exclude 9b0q/9b0s-type cases) and carry per-residue PP and per-window flags so that the front end can grey out low-confidence regions.

## Unexpected behaviour (documented, not worked around)

- cmalign and cmsearch bit scores differ substantially for divergent chains (5mrc 648 vs 260).
- 9gyt model span ends at 1528 (cmalign) vs 1533 (Rfam cmsearch).
- The Rfam `pdb_full_region.txt` FTP table lacks all 165 instances released in 2026; the `structures` endpoint is current.
- The PDBe molecules endpoint `source` field is empty for most ribosome entries.
- 4v4b is labelled "18S Ribosomal RNA" but its sequence is bacterial 16S by score.
- 6ydp/6ydw deposit an entire mitochondrial genome-length RNA entity (16,770 nt) carrying multiple Rfam hits on one chain.
