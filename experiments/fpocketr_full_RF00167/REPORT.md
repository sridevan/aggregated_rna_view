# fpocketR across all RF00167 RNA chain instances

Date: 2026-10-16. Analysis only: no frontend change, no change to SASA, loop/ligand annotations, Rfam mappings or the clustering worktree. Module `rna_explorer/pockets.py` (CLI `python3 -m rna_explorer pockets RF00167 --experiment-dir experiments/fpocketr_full_RF00167`) reuses the pilot's input, parsing and mapping logic; outputs `data/families/RF00167/pockets.json` (525 KB) and `pocket_tracks.json` (160 KB), schemas `pockets.schema.json` / `pocket_tracks.schema.json`. Raw fpocket/fpocketR outputs, inputs and residue maps stay in `experiments/fpocketr_full_RF00167/{inputs,work/runs}`; `repeat/` holds the reproducibility rerun; `analysis.json` holds every number below (`scripts/analyze_full.py`).

## A. Execution summary
- Environment: `tools/fpocketr-env` (fpocket 4.0.3 conda-forge h8a118fb_2, fpocketR 1.3.6, Python 3.11.17 x86_64, PyMOL 2.5.0) under Rosetta 2, unchanged since the pilot; relative `--out` workaround kept.
- Parameters: `-m 3.0 -M 5.7 -i 42 -D 1.65 -A 3 -p 0.0 -dpi 10 -y`, `--chain <auth chain>`; ligand-free single-chain inputs only (the pilot's holo inputs were used solely for fpocketR's Known/Novel check and are not part of this run).
- Inputs: 84/84 instances prepared and round-trip validated; 738 alternate-location copies dropped (87 per 9NBH/9NDD chain, 42 in 5SWD A); modified residues kept (GTP ×4, CCC ×2, GDP, 50N, 50L, IU); no insertion codes, no multi-character chain ids, no warnings.
- Runs: **84 attempted, 65 with ≥ 1 pocket, 19 with zero pockets, 0 failures**; 426 s total, 4.6–7.9 s per run (median 5.0 s; fpocket itself < 0.5 s, the rest is PyMOL/figures).
- Mapping checks: for all 70 pockets `PocketNT` = lining-atom residues, residue map = resolver (`label_seq_id`), 0 unresolved nucleotides, one membership per nucleotide; insert/flank nucleotides kept with `rfam_consensus_pos: null` (6 pockets contain such nucleotides, 4 lie entirely on the 9NBH scaffold).

## B. Pocket statistics (70 pockets, 65 instances)
- Pockets per instance: 0 ×19, 1 ×60, 2 ×5 (3RKF D, 4FEN B, 4FEO B, 5C7W C, 9NBC A). fpocketR therefore behaves as a one-pocket detector on this family.
- Score: min 0.11, Q1 0.455, median 0.53, Q3 0.56, max 0.62; 59 pockets ≥ 0.4, 11 below.
- Volume (Monte-Carlo): median 190 Å³ (130–700); nucleotides per pocket median 9 (4–16); shapes: 47 balanced, 15 rod-like, 8 disc-like; fpocket "Druggability" ≈ 0 throughout.
- Most frequent consensus positions by distinct-pocket count: 56 (62), 26/51/79 (59), 55 (57), 25/80/81 (56), 54 (54), then 27/28 (23), 57 (21), 29/53 (13), 50 (9).
- Alpha spheres of the pocket overlapping the ligand site in holo chains: min 42, Q1 46, median 51 — the purine pocket sits at the lower edge of fpocketR's `-i 42` size criterion (18 of 70 pockets have ≤ 45 alpha spheres).

## C. Aggregated pocket occurrence (filter "all Score ≥ 0", primary candidate)
Denominator per position = instances whose match residue has coordinates (run succeeded for all): 79–80 at the pocket positions, 75 at consensus 51 (unmodelled in several chains). The independent reconstruction from pocket records reproduces every numerator/denominator.

| Consensus | supporting / evaluable | frequency | entries | distinct pockets |
|---|---|---|---|---|
| 25 | 56 / 79 | 0.71 | 49 | 56 |
| 26 | 58 / 80 | 0.73 | 51 | 59 |
| 51 | 59 / 75 | 0.79 | 51 | 59 |
| 54 | 54 / 80 | 0.68 | 48 | 54 |
| 55 | 57 / 80 | 0.71 | 50 | 57 |
| 56 | 61 / 80 | 0.76 | 52 | 62 |
| 79 | 58 / 80 | 0.73 | 49 | 59 |
| 80 | 56 / 80 | 0.70 | 49 | 56 |
| 81 | 56 / 79 | 0.71 | 49 | 56 |

Variable region (0.1–0.5): 27, 28, 29 (J1/2), 50, 53, 57, 58 (J2/3 strand), fed by a second, lower-scoring pocket in 4FEN/4FEO/3RKF D and by shifted single pockets in the apo chains (5E54 A/B, 5SWD A/B), 6P2H and the 2025 entries (9IWF, 9LJN, 9LKC, 9LKE, 9LMC, 9V4U, 9V51) whose pocket extends into J1/2 and J2/3. Sporadic positions (< 0.1, one or two instances): 21, 23, 30–32, 35–38, 46–49, 52, 59, 70, 71, 78, 83, 85 (second pockets of 4FEN/4FEO at the P2 stem, 5C7W pocket 2 on the P1 terminus, 9NBC pocket 2 in L2/L3). 36 of 102 positions have any occurrence; 66 are evaluated zeros.

## D. Filter sensitivity

| Filter | qualifying pockets | instances with ≥ 1 | positions non-zero | known-site frequency (25…81) | positions ≥ 0.5 |
|---|---|---|---|---|---|
| all Score ≥ 0 | 70 | 65 | 36 | 0.68–0.79 | the 9 pocket positions only |
| top-ranked only | 65 | 65 | 20 | identical to "all" (0.68–0.79) | identical |
| Score ≥ 0.4 | 59 | 56 | 26 | 0.60–0.69 | identical set |

The known site is identical under "all" and "top1" because every instance's top pocket is the ligand pocket wherever one exists; the second pockets (5 instances, Scores 0.11–0.43) add only sporadic positions (P1 terminus, P2 stem, L2/L3, J1/2) with ≤ 3 instances each, so low-ranked pockets introduce diffuse but low-frequency annotations rather than a competing site. Score ≥ 0.4 removes 11 pockets, 9 of which are the ligand pocket itself in lower-quality chains (3DS7 B 0.39, 5C7W C 0.33, 6UC7 B 0.37, 9LKC 0.37, 9V51 0.25, 5E54 B 0.39, 6P2H 0.36, 9NBH A/C) and lowers the site frequency by ~0.08 without changing its shape: it trades sensitivity for nothing specific here. Denominators are identical across filters (zero-pocket runs count as zeros).

## E. Experimental ligand-site comparison
Reference per chain: chain-level Arpeggio contacts of the ligand-class molecule on the deposited coordinates (`experiments/chain_specific_ligands/arpeggio_pilot/results`), "any" = ≤ 5 Å, "specific" = non-proximal classes; the pilot set 25, 26, 51, 54, 55, 56, 79, 80, 81 is used only as a reference frame.
- 75 chains have a chain-specific ligand contact set (holo); 9 have none (apo: 3IVN A/B, 5E54 A/B, 5SWD B, 9NBH A–D).
- Holo recovery: **58 of 75** chains have a pocket overlapping the contact set, and in all 58 it is the top-ranked pocket; residue-level recall of the specific contacts median 1.0 (min 0.25), precision against the ≤ 5 Å set median 1.0 (min 0.44), Jaccard median 1.0 (30 chains exact), best-pocket Score median 0.54, volume median 175 Å³. Lower Jaccard values (0.29–0.69) come from pockets that extend into J1/2 (27–29) and J2/3 (57, 58) in the 2025 guanine-II/Vibrio variants (9IWF, 9IWG, 9Lxx, 9V4x), in 3RKF A and 3LA5, and from 5SWD A whose single pocket is shifted away from the site.
- **17 holo chains yield no pocket at all**, including high-resolution xpt structures (4FEJ 1.5 Å, 2XNZ 1.59 Å, 3FO6 1.9 Å, 2EET 1.95 Å, 6UC9 1.94 Å, 3GOG 2.1 Å), 4TZX, the 3.0 Å "transitional" 6VWT/6VWV chains, 5C7U, 9V4V and all four 9NDD cryo-EM chains (fpocket's info file is empty: no pocket passed fpocket's own criteria). These are evaluated zeros in the track. Given the 42–51 alpha-sphere counts of recovered pockets, the most plausible cause is that the ligand-free purine cavity falls just below the 42-alpha-sphere minimum in these models **[inferred, not tested: parameters were not varied]**; mutants (A21G/U75C, C74U, A24U) and lower resolution appear among them, but 4FEJ at 1.5 Å shows that resolution alone is not the explanation.
- Apo chains: 3IVN A/B (2.8 Å, U65C) no pocket; 5E54 A/B and 5SWD B one pocket each (Scores 0.39–0.50) displaced to J1/2 + J2/3 (27–29, 50, 51, 53, 56, 57, 79) sharing 3–4 positions with the reference frame; 9NBH A–D one pocket each entirely on the scaffold (consensus null). No ligand-induced conformational change is inferred from this; apo pockets are simply geometrically different and lower scoring, consistent with the published apo behaviour.

## F. Construct and structural variability
- Method: X-ray 60/75 chains with pockets (mean 0.85 pockets), cryo-EM 5/9 (0.67; 9NBC has the site pocket, 9NBH scaffold-only pockets, 9NDD none).
- Resolution bands: ≤ 2.0 Å 22/27 chains with pockets; 2.0–2.5 Å 23/25; 2.5–3.0 Å 17/25; > 3.0 Å 3/7. Detection drops at > 3 Å but is not monotonic (see 4FEJ).
- Constructs: aptamer-only chains (≤ 72 nt, n = 75) give the same known-site frequencies as the full set (0.71–0.80), so the 160-nt scaffold chains (9NBC/9NBH/9NDD) do not distort the profile; they contribute 4 scaffold-only pockets (insert/flank nucleotides, excluded from consensus aggregation) and one L2/L3 pocket (9NBC pocket 2).
- Modified nucleotides: 6/7 chains with modified residues have a pocket; GDP/GTP 5′ residues did not generate pockets.
- Flags: 4 pockets include an incomplete nucleotide (5′-terminal, from `sasa.json`), 4 touch a chain terminus (5C7W pocket 2 on P1, 9NBH scaffold pockets), 6 include non-consensus nucleotides; none of these affect the known-site frequencies.

## G. Scientific interpretation
Across homologous experimental purine riboswitch structures, fpocketR predicts essentially one pocket per chain, and where it predicts one it is the purine-binding pocket in 58 of 65 cases (the remaining 7 are the apo/shifted pockets and one scaffold case). The aggregated occurrence profile is therefore sharp: nine consensus positions at 0.68–0.79, a shoulder of J1/2–J2/3 positions at 0.1–0.3, and evaluated zeros elsewhere. What the data do **not** show: (i) the ~25 % of chains without a pocket at the site are not evidence of an absent cavity but of a geometric criterion (≥ 42 alpha spheres of radius 3.0–5.7 Å) that the purine pocket only just meets; (ii) pocket occurrence is a geometric statement on ligand-stripped coordinates, not experimental binding; (iii) apo-state behaviour is based on 9 chains from 4 entries; (iv) parameters were validated by the authors on drug-like ligands (QED ≥ 0.3), so a family whose natural ligands are purines is a demanding test rather than a representative one.

## H. Recommendation for the frontend (not implemented)
- One aggregated track "Pocket occurrence" in a category "Predicted binding pockets" is appropriate: the profile is stable, reconstructible and interpretable, with explicit per-position denominators.
- Default filter: **all pockets with Score ≥ 0** (fpocketR's published default). It is identical to top-1 at every consistently predicted position, keeps the low-scoring ligand pockets that Score ≥ 0.4 would drop, and its only cost is a handful of ≤ 3-instance sporadic positions that the colour scale renders faintly; expose the other two variants in the downloadable `pocket_tracks.json` only.
- Legend: the existing frequency scale (fraction of evaluable instances). Tooltip: "Rfam residue 56 / In a predicted pocket in 61 of 80 evaluable RNA instances" (+ distinct pockets when > instances). Pinned details: entries, pocket ids, score range.
- Category ⓘ caveats: fpocketR 1.3.6 / fpocket 4.0.3 with the published RNA parameters on ligand-free single chains; geometric prediction, not experimental binding; zero means no pocket passed the size criteria, not absence of a cavity; frequency denominators exclude unmodelled residues; 19 chains contributed only zeros.
- Downloads: `pockets.json` (per-structure pockets with lining atoms/nucleotides, scores, volumes, flags) and `pocket_tracks.json` (three filters).
