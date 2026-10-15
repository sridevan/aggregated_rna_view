# RF00167 Rfam consensus mapping proof of concept

Self-contained; does not touch any production code. Rerun with Infernal 1.1.x `cmalign` on PATH:

    python3 -I scripts/01_build_inputs.py data results          # FASTA + chain metadata (incl. species) from PDBe API JSON + mmCIF
    cmalign -o results/rf00167_mapped.sto --sfile results/rf00167_mapped.scores data/RF00167.cm results/mapped_regions.fa
    cmalign -o results/rf00167_full.sto   --sfile results/rf00167_full.scores   data/RF00167.cm results/full_sequences.fa
    python3 -I scripts/02_parse_stockholm.py results/rf00167_mapped.sto results/chains.json results/mapped
    python3 -I scripts/02_parse_stockholm.py results/rf00167_full.sto   results/chains.json results/full
    python3 -I scripts/03_validate.py data results

Deliverables (results/):
- rf00167_mapped.sto            Stockholm alignment of the PDBe-mapped regions (primary)
- rf00167_full.sto              Stockholm alignment of full chains (secondary; exercises '.' insert and '~' local-end columns)
- mapped_residue_mapping.csv    per-residue mapping (primary deliverable)
- mapped_consensus_columns.csv  per (chain, consensus position) incl. deletions
- full_*.csv                    same, for the full-chain alignment
- validation_report.json        machine-readable validation output
- VALIDATION_REPORT.md          human-readable report, edge cases, recommendation

data/: raw PDBe API JSON, PDBe updated mmCIF files, RF00167.cm (Rfam 15).
