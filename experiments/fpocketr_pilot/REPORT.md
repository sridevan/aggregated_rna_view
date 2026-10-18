# fpocketR pilot on RF00167 (Apple M3 Pro, isolated x86-64 environment)

Date: 2026-10-16. Follows `experiments/fpocketr_feasibility/REPORT.md`. Scope: environment, 6 RNA chain instances, mapping validation, known-site comparison, filter sensitivity. No production pipeline, no frontend change, no modification of SASA, loop/ligand annotations, mappings or the clustering worktree. Nothing installed outside `tools/micromamba` and `tools/fpocketr-env`.

Files: `environment/ENVIRONMENT.md`, `environment/fpocketr-env.yml`, `environment/package_list.txt` (reproducible environment); `scripts/prepare_inputs.py`, `scripts/analyze_pilot.py`; `inputs/` (PDB inputs + `*_residue_map.csv` + `input_validation.json`); `runs/<name>_clean_out/` (complete fpocket/fpocketR output per run incl. `pockets/pocketN_atm.pdb`, `pocketN_vert.pqr`, `*_clean_info.txt`, `*_pocket_characteristics.csv`, figures) with logs; `runs_repeat/` (reproducibility), `runs_timing/`; `pockets_pilot.json` (mapped pockets, comparison, sensitivity).

## 1. Environment (verified)
Project-local x86-64 micromamba 2.9.0 running under the pre-existing Rosetta 2; environment `tools/fpocketr-env` with **fpocket 4.0.3** (conda-forge `h8a118fb_2`, Mach-O x86_64), **fpocketR 1.3.6** (PyPI), Python 3.11.17 (x86_64), pymol-open-source 2.5.0, prody 2.3.0, rdkit 2023.09.6, scipy 1.12.0, trimesh 4.2.0, biopython 1.83, numpy 1.26.4, pandas 2.1.4, seaborn 0.13.2, setuptools 79.0.1. All versions equal the Bioconda recipe pins; no newer fpocket was substituted. Two deviations from the README route, both documented in `environment/ENVIRONMENT.md`: pinned Python dependencies were installed as conda-forge binaries because pip could not build prody 2.3.0 from source under Rosetta (pip then installed fpocketR with `--no-deps`; `pip check` clean), and setuptools < 80 was added because prody 2.3.0 imports `pkg_resources`. `python -m fpocketR --help` confirms the parameter names and defaults: `-m` 3.0 (min α-sphere radius), `-M` 5.70 (max), `-i` 42 (min α-spheres per pocket), `-D` 1.65 (clustering distance), `-A` 3, `-p` 0.0, `--qualityfilter` 0.0. Practical finding: fpocketR 1.3.6 fails after fpocket when `--out` is an absolute path (`analysis.strip('/')` in `pocket.py` removes the leading slash → "analysis directory does not exist"); a relative `--out` works.

## 2. Inputs
Isolated chains written by `scripts/prepare_inputs.py` (project Python + gemmi, read-only) from the cached PDBe updated mmCIF: all polymer residues of the instance's `label_asym_id` (construct extensions kept), heavy atoms, altloc = highest occupancy per atom name, model 1, modified residues written as ATOM records with their CCD names, author chain id and author numbering; a `*_residue_map.csv` records `pdb_chain, pdb_resnum, pdb_icode, resname ↔ label_asym_id, label_seq_id` per residue, and the written PDB was re-read and compared (`roundtrip_ok` true for all 10 files). Holo variants append only the ligand copy nearest the chain for fpocketR's own Known/Novel classification (`--ligand`); fpocketR removes it before pocket detection.

| Instance | Method / Å | Residues | Heavy atoms | Altloc copies dropped | Modified | Ligand (holo variant) |
|---|---|---|---|---|---|---|
| 4FE5 B (representative) | X-ray 1.32 | 67 (auth 15–81) | 1,422 | 0 | – | HPA |
| 1Y26 X | X-ray 2.1 | 71 (13–83) | 1,499 | 0 | – | ADE |
| 1Y27 X | X-ray 2.4 | 68 (14–81) | 1,456 | 0 | GDP (5′), CCC | GUN |
| 5E54 A (apo) | X-ray 2.3 | 65 (13–77), 6 nt unmodelled | 1,383 | 0 | – | none |
| 3IVN A (apo, U65C mutant) | X-ray 2.8 | 69 (5–73) | 1,466 | 0 | – | none |
| 9NDD A (cryo-EM scaffold construct) | EM 2.67 | 160 (1–160) | 3,394 | 87 | – | OXG (copy A/201, 2.2 Å from chain) |

No numbering problems, no insertion codes, no unsupported residues; GDP and CCC were accepted by fpocketR (kept as polymer).

## 3. Runs and reproducibility
Command per input: `python -m fpocketR -pdb <file> -c <chain> -m 3.0 -M 5.7 -i 42 -D 1.65 -A 3 -p 0.0 -dpi 10 -y -o ../runs [-l <CCD>]`. All 10 runs exit 0 in 4.7–6.3 s wall each (fpocket alone 0.48 s, 21 MB RSS; the fpocketR Python/PyMOL stage dominates: 7.7 s, 292 MB RSS at `-dpi 10`). Repeat runs of 4FE5 B and 1Y26 X: pocket characteristics CSV identical; `*_clean_info.txt` identical for 4FE5, for 1Y26 only the `Volume` line differed (146.2 vs 151.2 Å³) — fpocket estimates pocket volume by Monte-Carlo sampling, so volume carries ~±3 % run-to-run noise while scores, α-sphere counts, SASA and lining atoms are deterministic. **[inferred from fpocket's documented Monte-Carlo volume]**

## 4. Predictions

| Run | Pockets (Score ≥ 0) | Rank 1 Score / Druggability | Volume Å³ | α-spheres | Pocket SASA Å² | Shape | Lining nucleotides (author) | Consensus positions |
|---|---|---|---|---|---|---|---|---|
| 4FE5 B | 1 | 0.55 / 0.001 | 140 | 43 | 34 | balanced (NPR 0.59/0.76) | 21, 22, 47, 50, 51, 52, 73, 74, 75 | 25, 26, 51, 54, 55, 56, 79, 80, 81 |
| 4FE5 B + HPA | 1 | 0.55 (Type **Known**: pocket overlap 1.0, ligand overlap 1.0, centre 0.47 Å) | 140 | 43 | 34 | balanced | same | same |
| 1Y26 X | 1 | 0.53 / 0.001 | 150 | 48 | 39 | balanced | 21, 22, 47, 50, 51, 52, 73, 74, 75 | 25, 26, 51, 54, 55, 56, 79, 80, 81 |
| 1Y26 X + ADE | 1 | 0.53 (**Known**: 1.0 / 1.0 / 0.55 Å) | 150 | 48 | 39 | balanced | same | same |
| 1Y27 X | 1 | 0.53 / 0.001 | 150 | 46 | 40 | balanced | 21, 22, 47, 50, 51, 52, 73, 74, 75 | 25, 26, 51, 54, 55, 56, 79, 80, 81 |
| 1Y27 X + GUN | 1 | 0.53 (**Known**: 1.0 / 1.0 / 0.54 Å) | 150 | 46 | 40 | balanced | same | same |
| 5E54 A (apo) | 1 | 0.50 / 0.13 | 190 | 50 | 62 | rod-like | 23, 24, 25, 46, 47, 49, 52, 53, 73 | 27, 28, 29, 50, 51, 53, 56, 57, 79 |
| 3IVN A (apo) | 0 | – | | | | | | |
| 9NDD A (cryo-EM) | 0 | – | | | | | | |
| 9NDD A + OXG | 0 | – | | | | | | |

Every holo-type chain yields exactly one pocket, and it is the purine-binding pocket: fpocketR's own criteria classify it as Known with complete alpha-sphere and ligand overlap and the pocket centre within 0.5 Å of the ligand. All fpocket "Druggability" scores are near zero, as expected for the polar, purine-sized pocket (80.8 % polar atoms, mean α-sphere radius 3.4 Å); fpocket's Score (0.5–0.55) is the ranking quantity.

## 5. Mapping validation
`PocketNT` in the characteristics CSV equals the residue set parsed from `pockets/pocketN_atm.pdb` for every pocket (`nt_match` true); the lining-atom file lists 26 atoms for 4FE5 pocket 1 (e.g. 7 atoms of A21, 5 of U22, 1 of U47). Author residues map through `*_residue_map.csv` and independently through `InstanceIndex.resolve_unit` to identical `label_seq_id`s (`map_consistent` true) and to consensus positions; e.g. 4FE5 B: auth 21 → label 7 → consensus 25; 47 → 33 → 51; 75 → 61 → 81. No pocket nucleotide fell on insert, flank or scaffold residues in this pilot (`n_insert_or_flank_nt` 0). Identifier correspondence is therefore reversible and lossless for these chains (no insertion codes, single-character chain ids).

## 6. Comparison with experimental ligand contacts
Reference: chain-level Arpeggio contacts on the deposited coordinates (`experiments/chain_specific_ligands/arpeggio_pilot/results`), specific interactions (non-proximal) and any contact (≤ 5 Å), for the same chain and ligand.

| Chain | Pocket positions (consensus) | Specific contacts | Any contacts | Recall vs specific | Precision vs any | Jaccard vs specific |
|---|---|---|---|---|---|---|
| 4FE5 B / HPA | 25 26 51 54 55 56 79 80 81 | 25 26 51 55 56 80 81 | 25 26 51 54 55 56 79 80 81 | 1.00 | 1.00 | 0.78 |
| 1Y26 X / ADE | same | same | same | 1.00 | 1.00 | 0.78 |
| 1Y27 X / GUN | same | same | same | 1.00 | 1.00 | 0.78 |

The predicted pocket is exactly the set of nucleotides within 5 Å of the bound purine; the two "extra" positions (54, 79: C50/U/A73 in riboswitch numbering) are the stacking neighbours that Arpeggio also lists as proximal contacts. Apo structures: 5E54 A gives one lower-scoring, rod-like pocket displaced onto the J2/3 strand and J1/2 (consensus 27–29, 50–57, 79) sharing only positions 51 and 56 with the ligand site — consistent with the published observation that apo pockets "do not always perfectly align" with holo sites; 3IVN A (apo U65C mutant, 2.8 Å) gives no pocket. 9NDD A (cryo-EM, 2.67 Å, 160-nt scaffold, 87 altloc copies dropped) gives no pocket even with the 8-oxoguanine present in the input — the binding site is occupied in the model, but fpocketR removed the ligand and still found no cavity meeting the 42-α-sphere/polar criteria; whether this reflects map resolution, model geometry (B-factors/altlocs) or the construct cannot be decided from one case and parameters were not tuned.

## 7. Filter sensitivity
Per instance, nucleotides were deduplicated across pockets (set union). Because every pilot chain produced at most one pocket with Score 0.50–0.55, the three alternatives are identical here:

| Filter | Pockets retained (4FE5, 1Y26, 1Y27, 5E54, 3IVN, 9NDD) | Consensus positions with any occurrence | Occurrence at the family pocket (25, 26, 55, 80, 81 / 51 / 56) |
|---|---|---|---|
| all Score ≥ 0 | 1, 1, 1, 1, 0, 0 | 15 | 3/6 instances (3/5 at 51), 4/6 at 56 |
| top-ranked only | identical | 15 | identical |
| Score ≥ 0.4 | identical | 15 | identical |

Denominators: 6 instances observable at every position except 51 (5: unmodelled in 5E54 A). The pilot therefore cannot discriminate the filters; the family-wide run must record all pockets so that the three variants can be compared where multiple pockets occur (expected in scaffold constructs). Construct sensitivity is visible already: the 160-nt scaffold chain contributes zeros across the aptamer while the 67–71-nt aptamer chains contribute the pocket.

## 8. Limitations and caveats
- Only 6 chains; three holo chains are closely related constructs with identical numbering, so the perfect agreement partly reflects redundancy.
- Pocket volume is stochastic (±3 %); use Score, α-sphere count and lining atoms for comparisons, report volume as approximate.
- fpocket "Druggability" is protein-trained and uninformative here; fpocketR ranking relies on Score plus geometry.
- The holo-with-ligand variant only adds fpocketR's Known/Novel classification; detection itself is always ligand-free, so a tightly bound ligand leaves a cavity that a true apo structure may not have (5E54 illustrates the apo shift).
- Absolute `--out` paths break fpocketR 1.3.6; multi-character chain ids and insertion codes were not exercised; the x86-64-under-Rosetta environment is validated only by internal consistency, not yet against a Linux x86-64 run.
- Cryo-EM/scaffold case (9NDD) yielded nothing: expect "no pocket" to be a frequent, structure-quality-dependent outcome that must be reported as an evaluated zero, distinct from failed runs.

## 9. Recommendation
Proceed to all 84 RF00167 instances with the same environment, inputs and parameters (`-dpi 10`, relative `--out`, one run per chain, holo `--ligand` variant for single-chain entities with a default-track ligand), expected runtime ≈ 10 min. Record every pocket (section D.7 of the feasibility report) and compute the three filter variants before choosing the primary definition; keep "evaluated, no pocket" (3IVN, 9NDD) separate from failures in the denominators; flag pockets that include incomplete nucleotides (from `sasa.json` status) or terminal/scaffold residues; and run one Linux x86-64 cross-check of 4FE5 B when such a host is available.
