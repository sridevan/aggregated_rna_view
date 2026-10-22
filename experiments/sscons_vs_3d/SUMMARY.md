# Rfam SS_cons vs 3D base pairs: summary (2026-10-22)

> **Superseded in part by `FOLLOWUP_REPORT.md` (2026-10-22):** 16S is E. coli only; only RNase P 324–342 is sheared (325–341 / 326–340 are unpaired); TPP 64–83 is a G·G cWH alternative, not a register shift; SAM-I 50–64 vs 50–66 is an alignment-placement problem; low agreement of 3QSY C and 6AHU M is mostly FR3D annotation sensitivity.

The full report is `REPORT.md` in this folder; data and script are alongside it (`analyse.py`, the CSVs, `summary.json`). This is an investigation only: no app changes.

## What was compared

For each Rfam consensus base pair, every chain where both positions have a modelled residue with a complete base and FR3D annotations. Pairs inside a stretch with no FR3D annotations (≥ 15 positions) are left out. Support = chains forming a Watson–Crick (cWW) pair / those eligible chains.

| Family | Pairs with support ≥ 0.9 | Pairs with support < 0.5 | Median per-chain agreement |
|---|---|---|---|
| Purine riboswitch (RF00167) | 18/22 | 0 | 0.95 |
| SAM-I (RF00162) | 24/33 | 4 | 0.93 |
| TPP (RF00059) | 19/28 | 3 | 0.90 |
| tRNA (RF00005) | 17/21 | 0 | 0.95 |
| RNase P B (RF00011) | 72/97 | 8 | 0.89 |
| 16S (RF00177) | 375/462 | 32 | 0.92 |

## Consensus pairs not formed as Watson–Crick in 3D

- **Wrong pair type.** Some consensus pairs form a different, non-Watson–Crick pair:
  - SAM-I 46–69 is a sheared (tSH) pair in 26 of 29 chains;
  - RNase P has a run of tandem sheared G·A pairs (324–342, 325–341, 326–340), plus 129–162 and 7–361;
  - 16S: 139–171 and 712–723.
- **Shifted by one position.** The partner is one position off:
  - SAM-I 50–64: 50 pairs with 66 in about 96% of chains;
  - TPP 20–28 / 21–27 (shifted in about 30% of chains) and 64–83 (64 pairs with 82);
  - 16S 372–385: 372 pairs with 386 in all 60 chains;
  - 16S 648–737.
- **16S pairs where a position pairs elsewhere**, e.g. 1228 with 1328, 1520 with 807, 1149 with 1172; a few others have no partner at all. These may be places where Rfam's consensus differs from the E. coli / T. thermophilus structures, or where the alignment placed columns against the 3D pairing. Candidates, not proven Rfam errors.
- **A 16S hairpin (75–86 to 78–83)** forms in only 23–49% of chains.

## 3D pairs missing from the consensus

Mostly the expected tertiary structure, which a nested secondary structure cannot express:
- the purine L2–L3 kissing loop (41–66 in 80 of 80 chains);
- the SAM-I pseudoknot (25–28 with 77–80);
- the tRNA tertiary pairs (18–54, the Levitt pair 15–46, 8–14).

## Chains with low agreement

Failures sit in specific helices, not across whole chains, so they point to real structural states rather than mapping errors:
- tRNA bound in complexes: 3QSY 0.43, 6AHU 0.48, 1J2B with its remodelled D arm;
- the Guanine-II riboswitch 9V50 and the scaffold construct 9NBC;
- 16S chains lacking the often-unformed hairpin.

## Caveats

- **Circularity.** Chains are placed on Rfam positions with the Rfam model, which is built on the consensus structure. Agreement is an upper bound; disagreement is the more informative signal.
- **Near pairs are missing.** The per-chain files contain no FR3D "near" base pairs, so that status could never be observed.
- **Two families are dominated by a few species.** RNase P is a single RNA; 16S is mostly two species.

## Recommendation (not built)

Add a consensus-structure layer to the Interactions tab:
- each consensus pair coloured by support, shown as "formed as Watson–Crick in x / y chains" with a status breakdown;
- a toggle for the recurrent 3D-only pairs;
- a per-chain "consensus pairs formed: x / y" line.

No new data sources are needed.

**Open decisions:**
- Non-canonical Watson–Crick pairs (G·A, U·U): count them as support (as here) or show them separately?
- 16S: summarise its 462 pairs by helix rather than show each pair?
