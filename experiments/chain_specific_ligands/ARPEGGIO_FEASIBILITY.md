# Local Arpeggio feasibility

Setup (2026-10-15, macOS arm64, Python 3.11.9): `python3 -m venv`; `pip install pdbe-arpeggio openbabel-wheel gemmi` → pdbe-arpeggio 1.4.4, openbabel-wheel 3.1.1.23 (prebuilt wheel, no conda/system Open Babel needed), gemmi 0.7.5, BioPython. Install time ≈ 1 min. Script: `scripts/run_arpeggio_poc.sh`.

| Aspect | Finding |
|---|---|
| Input | PDBe updated mmCIF from the project cache; selection `-s /CHAIN/RESNUM/` (author ids) or `RESNAME:<CCD>` |
| Runtime | 1Y26 (1 ligand): 0.46 s wall, 67 MB; 3RKF (4 DX4 copies in one run): 1.4 s, 107 MB. 66 entries would take well under a minute |
| Output | `<pdb>.json`: list of atom-pair contacts `{bgn, end, contact[], distance, interacting_entities, type}` with `auth_asym_id`, `auth_seq_id`, `auth_atom_id`, `label_comp_id`, `pdbx_PDB_ins_code`; `interacting_entities` ∈ INTER, INTRA_NON_SELECTION, SELECTION_WATER, NON_SELECTION_WATER, WATER_WATER. 1Y26: 769 contacts, 260 INTER |
| Identifiers | author only; mapped through `InstanceIndex.resolve_unit` (author chain + number + insertion code) → 48/48 rows resolved to label ids and consensus |
| Chain and ligand selection | explicit; all copies of a ligand can be processed in one run (3RKF A/B/C/D) → covers the copies PDBe skips |
| Classifications | with the PDBe updated mmCIF (no hydrogens; `type_symbol` present) the default run yields mostly `proximal`, some `vdw`, `vdw_clash`, `polar`, `hbond`; **no atom-plane / plane-plane** records, and `-mh` changed nothing. PDBe's deployment reports aromatic, CARBONPI/CATIONPI, FF/OF stacking etc. for the same ligand, so PDBe uses additional preprocessing (protonation, ring perception) that the bare CLI does not reproduce here. Residue sets nevertheless agree: 1Y26 ADE local = PDBe's 7 residues + 2 proximal-only (consensus 54, 79) |
| Modified nucleotides | Arpeggio reads components by `label_comp_id`; standard atom typing covers A/C/G/U; modified residues depend on Open Babel perception (not tested; RF00167 has few) |
| Assemblies and symmetry | runs on the file as given (asymmetric unit, model 1). Symmetry-mate contacts need an assembly file (PDBe `*_assembly.cif` or ModelServer) as input |
| Waters | included (SELECTION_WATER records); filterable |

Verdict: local Arpeggio is **not necessary** for chain-specific residue identification (PDBe + RCSB together cover every ligand instance), but it is cheap and the only way to obtain atom-level, typed interactions for ligand copies outside PDBe's preferred assembly or for symmetry-aware analysis, provided the hydrogen/ring preprocessing gap is resolved (compare against PDBe's `bound_ligand_interactions` for processed copies as a calibration set).
