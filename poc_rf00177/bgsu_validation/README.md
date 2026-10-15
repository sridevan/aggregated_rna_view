# BGSU RNA 3D Hub correspondence validation (RF00177)

data/      BGSU API downloads (pairwise/multi align_chains with column, SeqtoUnitMapping), selection for the documented pair
scripts/   01_normalise_bgsu.py (unit id -> mmCIF polymer position), 02_compare.py (metrics, discrepancies, regions, unnumbered regions), 03_columns_and_unnumbered.py (BGSU column -> consensus; nature of pairs inside insert/local-end regions)
results/   bgsu_correspondences_normalised.csv, bgsu_alignments_normalised.csv, bgsu_normalisation_stats.json, repro/ (our mapping of 5J7L/1J5E), pair_metrics.csv, comparison_summary.json, discrepancies.csv, region_windows.csv, unnumbered_regions_bgsu.csv, unnumbered_bgsu_pairs_detail.csv, bgsu_column_to_consensus.csv, BGSU_COMPARISON_REPORT.md

Rerun: python3 -I ../scripts/02_build_inputs.py ../data intermediate data/selection_repro.json ; cmalign ... ; python3 -I ../scripts/03_parse_stockholm.py ... results/repro/rf00177 ;
       python3 -I scripts/01_normalise_bgsu.py data ../data/cif results ; python3 -I scripts/02_compare.py results ../results ; python3 -I scripts/03_columns_and_unnumbered.py results ../results
