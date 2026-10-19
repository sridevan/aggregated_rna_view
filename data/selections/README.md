# Curated chain selections

A selection limits `python3 -m rna_explorer build <ACC> --selection <csv>` to the listed chains. The required columns are `pdb_id` and `label_asym_id` (or the 16S file's `ssu_label_asym_id`); `auth_asym_id` / `ssu_auth_chain` are optional. Any other columns are curated metadata, carried along but not used by `build`.

| File | Source (read-only research worktree `~/Desktop/aggregated_rna_clustering`, branch research/rna-clustering) | Copied |
|---|---|---|
| `RF00177_60_chains.csv` | `experiments/rna_structural_clustering/data/curated/rf00177_expanded_set.csv`: the approved 60-chain *E. coli* 16S set, 22 benchmark chains plus 38 additions (`results/rf00177/dataset_expansion_investigation/investigation_report.md`, `RF00177_EXPANDED_60_RESULTS.md`) | 2026-10-17, unchanged (md5 identical) |
| `RF00177_60_states.csv` | `.../data/curated/rf00177_expanded_states.csv`: state categories from deposition titles and deposited polymer entities | 2026-10-17, unchanged |
| `RF00005_60_chains.csv` | `experiments/rna_structural_clustering/data/curated/trna_demo_set.csv`: the 60-chain tRNA demo set (the proposed set with 9PS7 NC replaced by 8Z1F A and 1F7U A by 2AKE A; `results/trna_dataset_investigation/investigation_report.md`, `RF00005_DEMO_CLUSTERING_RESULTS.md`) | 2026-10-18, unchanged (sha256 identical) |
| `RF00011_50_chains.csv` | `experiments/rna_structural_clustering/data/curated/rf00011_selection_50.csv`: the 50-chain RNase P class B subset (24 distinct-condition structures + 26 stratified Degenhardt conformers; `scripts/rnasep_b_select50.py`, `RF00011_RNASEP_B_CLUSTERING_RESULTS.md`) | 2026-10-19; source sha256 fda32b12…; columns `label_asym_id` and `auth_asym_id` added from the research build's `data/rf00011_50/.../instances.json` (all A/A), every other column unchanged |

**RAD values:** `rad_*` columns are RAD (Ribosome Analysis Database) values for the 22 benchmark chains only; they are empty for the additions. Head-rotation values for the additions were computed by the research and are not in these files.

## Deployed copy

The site copy (`npm run build` → `frontend/dist/families`) leaves out pipeline-only files that the app never fetches and the download menu does not offer: `work/`, the per-chain `annotations/` files and the CSV tables. They stay in `data/families/`.
