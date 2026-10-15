# RF00167 (purine riboswitch) consensus residue mapping: proof-of-concept validation report

Date: 2026-10-14. Infernal 1.1.5 (built from source). Covariance model: RF00167 from rfam.org (Rfam 15 build, CLEN 102, DATE 2024-07-04).
Source of chain list: PDBe API `nucleic_mappings/RF00167` (66 entries, 84 mapped segments, all entity 1, mapped lengths 55-98 nt).

## Chains selected

| PDB | Entity | label_asym | auth chain | Chain length | PDBe mapped region (polymer) | author numbering | Mapped nt | Why chosen |
|---|---|---|---|---|---|---|---|---|
| 6p2h | 1 | A | A | 69 | 6-60 | 22-76 | 55 | shortest mapped segment; 5'-terminal GTP modified residue |
| 1y26 | 1 | A | X | 71 | 3-69 | 15-81 | 67 | canonical adenine riboswitch; auth chain differs from label asym |
| 3ds7 | 1 | B | B | 67 | 1-67 | 215-281 | 67 | author numbering offset by +214 in chain B only |
| 4fe5 | 1 | A | B | 67 | 1-67 | 15-81 | 67 | auth chain B but label asym A |
| 9nbh | 1 | A | A | 160 | 16-113 | 16-113 | 98 | longest segment, covers the full model; cryo-EM scaffold construct |

## Alignment results (mapped regions only, `cmalign` default local mode, Stockholm output)

| Chain | CM from-to | trunc | bit score | avg PP | residues in consensus columns | insert residues | local-end residues | deleted consensus positions |
|---|---|---|---|---|---|---|---|---|
| 6p2h/A | 26-80 | no | 36.3 | 1.000 | 55/55 | 0 | 0 | none |
| 1y26/X | 19-87 | no | 63.4 | 0.999 | 67/67 | 0 | 0 | 64, 72 |
| 3ds7/B | 19-87 | no | 42.6 | 0.999 | 67/67 | 0 | 0 | 64, 72 |
| 4fe5/B | 19-87 | no | 64.9 | 0.999 | 67/67 | 0 | 0 | 64, 72 |
| 9nbh/A | 1-102 | no | 40.9 | 0.986 | 98/98 | 0 | 0 | 39, 51, 64, 72 |

Every residue submitted was placed in a consensus (match) column. No `.` insert columns and no `~` local-end columns were produced for the mapped regions, so those code paths were exercised only by the secondary full-chain alignment (see edge case 2), where the parser handled 2 insert columns and 102 local-end columns correctly.

## Validation checks

1. **mmCIF identifier agreement.** Every residue row in the CSV (434 rows over 5 chains, including flanks) was compared with `_pdbx_poly_seq_scheme` in the PDBe updated mmCIF: label asym id, author strand id, polymer `seq_id`, `pdb_seq_num`, `pdb_ins_code`, `mon_id`. Result: 0 mismatches. The PDBe residue-listing API and the mmCIF agree exactly, including the 3ds7 chain B offset numbering and the 6p2h GTP residue. No insertion codes occur in any of the 84 RF00167 segments, so insertion-code handling is implemented but untested.
2. **Continuity.** Consensus positions are strictly monotonic along each chain. Offset changes occur only where a consensus position is deleted (e.g. polymer 47 to 48 in 1y26 jumps consensus 63 to 65). Deleted consensus positions 64 and 72 are shared by all five chains, i.e. these two model positions will have no structural observations from any of these entries. The per-column table `mapped_consensus_columns.csv` records deletions explicitly.
3. **Coverage.** 100 % of each mapped region lands in consensus columns. Model coverage per chain ranges from 54 % (6p2h, CM 26-80) to 96 % (9nbh). Across the 5 chains, all 102 consensus positions have at least one residue and 51 positions have a residue from all five chains.
4. **Agreement with the PDBe/Rfam segment boundaries.** The cmalign model span for each mapped region (26-80, 19-87, 1-102) is a clean hit with trunc = "no", consistent with the segments having been derived from a cmsearch hit with the same model. The parsed span equals the `--sfile` span in all cases.
5. **Posterior probabilities.** All residues have PP >= 0.85 except 9nbh polymer position 65 (PP bin 8). In default (truncation-aware) mode the first seven 9nbh residues get PP "9" rather than "*"; with `--notrunc` the residue placement is identical but PP becomes "*". PP is carried per residue in the CSV.

## Edge cases identified

1. **Mapped region with no coordinates (9nbh).** Polymer positions 31-102 of every 9nbh chain (72 residues, essentially the entire aptamer core) are absent from `atom_site`; only 26 of the 98 mapped residues are modelled. The PDBe mapping is SEQRES-based and does not know this. Sibling entry 9nbc has 97 unmodelled residues; 9ndd has none. Any aggregated structural view must join the mapping with observed status (PDBe `observed_ratio` or mmCIF `auth_seq_num == '?'`), otherwise consensus positions 16-91 will appear "covered" by 9nbh with nothing to display.
2. **Alignment is input-dependent (6p2h).** Aligning the full 69-nt chain instead of the mapped 55-nt region changes the assignment of polymer residues 59-60 (author 75-76): they become insertions, and consensus positions 79-80 are instead taken by residues 61-62. Both alignments are legitimate local optima; the full-chain version also reports a 5' truncation and a slightly higher bit score (36.5 vs 36.3). For 9nbh the full-chain alignment is catastrophic (bit score -13.2, 102 residues in local-end columns, only 58 matched). Conclusion: always align the Rfam/PDBe-defined segment (or a cmsearch hit), never the whole chain, and record which region was aligned.
3. **Author vs label identifiers.** 1y26 (auth X / label A) and 4fe5 (auth B / label A) have different author and label chain IDs; 3ds7 chain B author numbering starts at 215. Author numbering is not a stable key even within an entry. Polymer sequence position (label_seq_id) must be the primary key; author number + insertion code are display attributes.
4. **Species provenance is inconsistent.** For 3ds7 and 4fe5 the PDBe API reports "synthetic construct" and "Bacillus subtilis" respectively, while the wwPDB mmCIF `_pdbx_entity_src_syn` has organism and taxonomy ID as "?" (details: "mutagenesis study", "T7 RNA transcript"). 3ds7 is in fact the B. subtilis xpt-pbuX guanine riboswitch but is labelled synthetic. Three of five chains are "synthetic construct" (taxid 32630). Species for in vitro-transcribed riboswitch constructs is a curation question, not a data field that can be trusted directly; the CSV carries both API and mmCIF values.
5. **Modified residues.** 6p2h position 1 is GTP (PDBe `pdb_sequence` "(GTP)GG...", one-letter parent G). It lies outside the mapped region here, but modified residues inside a mapped region require the parent one-letter code for cmalign and the 3-letter code in the output. Unknown/non-standard residues (N/X) were not encountered and are untested.
6. **Entity vs chain.** Alignment is a property of the entity sequence; 3ds7 (2 chains), 9nbh (4 chains) and 3rkf (4 chains) share one sequence per entity. Align once per entity and fan out to chains, joining per-chain observed status. This also guarantees identical consensus assignment across copies.
7. **Insert and local-end columns.** Correctly parsed (`.` = insert, `~` = local end, lowercase residues in those columns, `-`/`.`/`~` as gaps) but only seen in the full-chain alignment. When aligning well-defined segments with this model they appear to be rare; they will be common for families with variable-length loops and must be kept in the output as "insert after consensus N, index k".
8. **Partially observed terminal residues.** First/last residues frequently have `observed_ratio` 0.85-0.96 (missing terminal atoms). Harmless for mapping, relevant for per-atom features.
9. **Model versioning.** Consensus numbering is only meaningful relative to the CM build (Rfam 15, CLEN 102). A seed change in a future Rfam release renumbers everything; mappings must be stored with the Rfam release and CM checksum (CKSUM 1292453307) and regenerated per release. It was not verified that the PDBe mapping and this CM come from the same Rfam release, but the span agreement suggests they do.

## Unresolved scientific assumptions

- The Rfam consensus column is treated as a stable, biologically meaningful coordinate for aggregating structures. Deleted consensus positions (64, 72 in all 5 chains) and insert residues show it is an approximation of homology, not a one-to-one residue correspondence.
- cmalign's single optimal-accuracy alignment is taken as the truth. Alternative near-optimal alignments (6p2h 3' end) are not represented; PP values partially capture this.
- The PDBe nucleic_mappings segment boundaries are taken as authoritative. They are SEQRES-based cmsearch hits and ignore whether residues are modelled.
- Species of synthetic constructs: unresolved. Options are to display sample source with the taxid caveat, or to derive the "natural origin" from the Rfam full-alignment hit, which was not attempted here.
- Modified-nucleotide handling inside mapped regions is assumed to work via the parent one-letter code; not tested.

## Recommendation

The approach is reliable enough to scale, with conditions. Identifier agreement with mmCIF was perfect, continuity and coverage behaved exactly as the alignment predicts, and cmalign cost is negligible (about 0.01 s per sequence; 84 segments for this family). Conditions for the full pipeline:

1. Align the PDBe/Rfam mapped segment per entity, never the full chain; record the aligned range and the cmalign bit score, truncation flag and PP per residue.
2. Key everything on (pdb_id, entity_id, label_asym_id, polymer position); carry author chain, author number and insertion code as attributes.
3. Join with observed status from the residue listing or mmCIF before any aggregation or display.
4. Keep insert and local-end residues in the output with their anchor consensus position; keep deletions in a per-consensus-column table.
5. Store Rfam release and CM checksum with every mapping; regenerate on Rfam release change.
6. Treat species as provisional; store both API and mmCIF values and flag disagreement.
7. Before scaling to all families, test a family with long variable insertions (e.g. a group II intron or an rRNA family) and one with insertion-coded residues, since neither occurs in RF00167.
