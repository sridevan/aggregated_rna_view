# SASA pilot report (RF00167, 6 representative RNA chains)

Date: 2026-10-16. Outputs: `sasa.json`, `sasa_tracks.json` (pilot subset), `independent_check_1y26_1_X_A.json`, scripts `rna_explorer/sasa.py` (production module, run with `python3 -m rna_explorer sasa RF00167 --instances … --pilot-dir experiments/sasa_pilot`) and `independent_check.py`.

## Method
FreeSASA 2.2.1 (`freesasa.calcCoord`: coordinates and radii supplied explicitly, so no protein residue classifier is involved), Shrake–Rupley, probe 1.4 Å, 1000 test points per atom. Heavy atoms only; element van der Waals radii C 1.70, N 1.55, O 1.52, P 1.80 (Bondi; P as NACCESS); S/Se/halogens defined for modified components. Chain = polymer atoms of the instance's `label_asym_id` in the PDBe updated mmCIF (gemmi 0.6.7 subchain), first model, so proteins, other RNA/DNA chains, ligands, ions and water are absent by construction. Alternate locations: per atom name the highest-occupancy location is kept (ties: first altloc label). Expected heavy atoms per residue: wwPDB CCD (RCSB ligand CIF, cached) minus `pdbx_leaving_atom_flag = Y` atoms (OP3); standard A/C/G/U expect 22/20/23/20 atoms.

## Test structures and results

| Instance | Why | Atoms | Chain total (Å²) | Nucleotides complete / incomplete / no coords | Per-nt range (Å²) | Modified |
|---|---|---|---|---|---|---|
| 4FE5 B (representative) | X-ray 1.32 Å | 1,422 | 11,051.7 | 66 / 1 / 0 | 74.5–308.8 | – |
| 1Y26 X | X-ray 2.1 Å, independent structure | 1,499 | 11,933.2 | 70 / 1 / 0 | 76.7–316.8 | – |
| 1Y27 X | modified nucleotides | 1,456 | 11,592.3 | 68 / 0 / 0 | 76.1–434.2 | GDP (5′ flank, 28/28 atoms, 434 Å²), CCC (match, 23/23) |
| 5C7W C | two modified match residues, 3.22 Å | 1,428 | 11,174.8 | 67 / 0 / 0 | 77.9–329.9 | 50N (23/23), 50L (23/23; leaving atom O8 handled) |
| 5E54 A | residues without coordinates | 1,383 | 11,469.3 | 65 / 0 / 6 | 76.0–357.6 | – |
| 9NDD A | cryo-EM, 160-nt scaffold construct, altlocs | 3,394 | 26,849.5 | 158 / 2 / 0 | 39.9–333.3 | – (87 altloc copies dropped) |

- Sum of per-nucleotide SASA equals the chain total within 0.004 Å² in every chain (identity by construction; asserted in tests).
- All values finite and non-negative; totals and ranges are physically plausible (a free nucleotide has roughly 300–450 Å²; buried helix residues 40–80 Å²; 5′-terminal GDP with free triphosphate-like group is the maximum).
- Incomplete nucleotides: the 5′-terminal residue lacking P/OP1/OP2 in 4FE5, 1Y26 and 9NDD (`incomplete_reason = 5prime_phosphate_absent`), and 9NDD U80 missing its whole base (7 atoms, `missing_heavy_atoms`). They are kept in `sasa.json` with the missing-atom list and excluded from the aggregated track.
- Modified nucleotides: all CCDs resolved (GDP, CCC, 50N, 50L here; GTP, IU, 23G in the full set), typed by element, complete → included with their own identity recorded in `nucleotide_identity`. No `unsupported` nucleotide occurred.
- Residue names in the structure matched the mapping in every case (no warnings).

## Independent numerical check
`independent_check.py` re-implements Shrake–Rupley in plain numpy (golden-spiral points, 2000 per atom, brute-force neighbour test) on the same 1,499 atoms and radii of 1Y26 X: chain total 11,926.1 vs FreeSASA 11,933.2 Å² (relative difference 0.06 %); per-nucleotide median relative difference 0.4 %, maximum 1.1 % (2.0 Å² on a 300 Å² nucleotide). The two implementations share no code, so this validates the FreeSASA call, radii and atom bookkeeping independently of the library.

## Suitability
No methodological problem found; the pipeline was run on all 84 instances (`data/families/RF00167/sasa.json`, `sasa_tracks.json`): 6,133 complete, 46 incomplete (38 missing 5′ phosphate, 8 missing heavy atoms), 422 without coordinates, 0 unsupported, 0 failures; every one of the 102 consensus positions has a median value (80 contributing instances at consensus 53, median 134.8 Å²). Colour domain for the track: 0–350 Å² (maximum median 307 Å² rounded up to 50).
