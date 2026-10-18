# Curated chain selections

A selection limits `python3 -m rna_explorer build <ACC> --selection <csv>` to the listed chains. The required columns are `pdb_id` and `ssu_label_asym_id` (label asym id); `ssu_auth_chain` is optional. Any other columns are curated metadata, carried along but not used by `build`.

| File | Source (read-only research worktree `~/Desktop/aggregated_rna_clustering`, branch research/rna-clustering) | Copied |
|---|---|---|
| `RF00177_60_chains.csv` | `experiments/rna_structural_clustering/data/curated/rf00177_expanded_set.csv`: the approved 60-chain *E. coli* 16S set, 22 benchmark chains plus 38 additions (`results/rf00177/dataset_expansion_investigation/investigation_report.md`, `RF00177_EXPANDED_60_RESULTS.md`) | 2026-10-17, unchanged (md5 identical) |
| `RF00177_60_states.csv` | `.../data/curated/rf00177_expanded_states.csv`: state categories from deposition titles and deposited polymer entities | 2026-10-17, unchanged |

**RAD values:** `rad_*` columns are RAD (Ribosome Analysis Database) values for the 22 benchmark chains only; they are empty for the additions. Head-rotation values for the additions were computed by the research and are not in these files.

## Deployed copy

The site copy (`npm run build` → `frontend/dist/families`) leaves out pipeline-only files that the app never fetches and the download menu does not offer: `work/`, the per-chain `annotations/` files and the CSV tables. They stay in `data/families/`.
