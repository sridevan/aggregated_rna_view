Chain-specific RNA–ligand interaction investigation (2026-10-15). Read `RECOMMENDATIONS.md` first.
- `DATA_SOURCE_ASSESSMENT.md` endpoints tested, response formats, availability
- `IDENTIFIER_VALIDATION.md` mapping of every source to (pdb_id, label_asym_id, label_seq_id) and consensus positions
- `ARPEGGIO_FEASIBILITY.md` local Arpeggio proof of concept
- `CONTACT_COMPARISON.csv` 499 residue-level contact rows for 1Y26, 3RKF, 9NDD, 9NBH from five sources (`cache/comparison_summary.txt` = per-ligand summary)
- `scripts/fetch_sources.sh` (retrieval), `scripts/compare_contacts.py` (comparison; read-only use of `rna_explorer.annotation_mapping`), `scripts/run_arpeggio_poc.sh`
- `cache/` 145 raw responses and Arpeggio outputs (4.6 MB)
Nothing in `rna_explorer/`, `frontend/` or `data/` was modified.

`arpeggio_pilot/` (2026-10-15/16): Arpeggio on deposited wwPDB mmCIF (asymmetric unit, no assembly/symmetry): all 66 RF00167 entries (68 s), 7K00 70S 16S chain (6.7 h, see runtime caveat), chain-level tracks, symmetry omission check. Read `arpeggio_pilot/ARPEGGIO_PILOT.md`.
