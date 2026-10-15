# RF00177 (bacterial SSU rRNA) Rfam consensus mapping proof of concept

Extends `poc_rf00167`. Raw downloads in `data/` (api/, cif/, cm/), intermediates in `intermediate/`, deliverables in `results/`.
Needs Infernal 1.1.x `cmalign` on PATH and matplotlib.

    python3 -I scripts/00_fetch_molecules_batch.py data        # molecules JSON for all 1473 entries (POST batch)
    python3 -I scripts/01_instances.py data results            # rf00177_instances.csv + summary (organism via PDBe search API)
    python3 -I scripts/01b_clan_competition.py data results    # Rfam bits/E-value/model span + SSU clan competition
    python3 -I scripts/02_build_inputs.py data intermediate    # mapped-region FASTA, chains.json (species from mmCIF), modified-base normalisation
    cmalign --cpu 4 -o intermediate/rf00177_mapped.sto --sfile intermediate/rf00177_mapped.scores data/cm/RF00177.cm intermediate/rf00177_sequences.fasta
    python3 -I scripts/03_parse_stockholm.py intermediate/rf00177_mapped.sto intermediate/chains.json data/cif results/rf00177
    python3 -I scripts/04_validate.py data intermediate results
    python3 -I scripts/06_reference_check.py results
    python3 -I scripts/05_insertion_profile.py data results intermediate
    # optional: cmalign -g ... intermediate/rf00177_mapped_global.sto ; scripts/07_local_vs_global.py

Deliverables (results/): rf00177_instances.csv, rf00177_sequences.fasta, rf00177_mapped.sto, rf00177_residue_mapping.csv,
rf00177_consensus_columns.csv, rf00177_insertion_profile.csv (+ _wide, _sites_summary, longest_insertions), rf00177_insertion_profile.png,
rf00177_window_confidence.csv, ecoli_7k00_numbering_to_consensus.csv, reference_correspondence.json, local_vs_global_stability.json,
validation_report.json, VALIDATION_REPORT.md.
