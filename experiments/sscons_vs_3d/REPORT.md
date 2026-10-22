# Rfam SS_cons base pairs vs base pairs observed in 3D (investigation, 2026-10-22)

> **Superseded in part by `FOLLOWUP_REPORT.md` (2026-10-22):** 16S is E. coli only; only RNase P 324–342 is sheared (325–341 / 326–340 are unpaired); TPP 64–83 is a G·G cWH alternative, not a register shift; SAM-I 50–64 vs 50–66 is an alignment-placement problem; low agreement of 3QSY C and 6AHU M is mostly FR3D annotation sensitivity.

**Question.** How well does each family's Rfam consensus secondary structure (SS_cons) agree with the base pairs that FR3D annotates in the experimentally determined chains? And which recurrent 3D pairs does SS_cons not contain?

**Status.** Investigation only. No app, backend or data changes.

**Script.** `analyse.py`. Run it with `python3 experiments/sscons_vs_3d/analyse.py [ACC...]`.

**Outputs, all in this folder:**
- `ss_cons_pair_support.csv`: one row per SS_cons pair.
- `helix_support.csv`: one row per helix.
- `recurrent_3d_pairs_not_in_ss_cons.csv`: recurrent 3D pairs absent from SS_cons.
- `chain_agreement.csv`: one row per chain.
- `summary.json`: per-family summary.

## 1. Method

### Inputs

Existing outputs only, read unchanged:
- `family.json` `ss_cons`;
- `interactions.json`:
  - per-chain eligibility, nucleotides and annotation gaps;
  - the aggregated `pairs`, with frequencies over eligible chains;
- `interactions_by_instance/<chain>.json`: each FR3D observation of a chain on Rfam positions.

These come from the validated mappings (match states only) and BGSU FR3D. Releases and provenance are in each family's `interactions.json`.

### SS_cons pairs

- **Parsing.** WUSS brackets `()<>[]{}`, with pseudoknot letters also handled.
- **None of the six families has pseudoknot letters.** So every Rfam consensus structure here is nested-only.

### Status of each SS_cons pair (i, j) in each chain

**Which chains count.** A chain counts when both positions are eligible: a match residue with coordinates, a complete base, and FR3D available.

**Gap exclusion.** Pairs where either position lies in one of the chain's annotation gaps are not evaluable. These are counted in `n_in_annotation_gap` and left out. An annotation gap is a run of 15 or more eligible positions without any BGSU annotation; for example, 7PJS 652–1531.

**Statuses:**

| Status | Meaning |
|---|---|
| `cww_canonical` | FR3D cWW, AU/GC/GU |
| `cww_noncanonical` | FR3D cWW with another combination (G·A, U·U, …) |
| `other_class` | another Leontis–Westhof family between i and j (e.g. tSH) |
| `near_only` | only an FR3D near pair |
| `shifted` | no i–j pair, but i or j forms a cWW with a neighbour of the other partner (register shift candidate) |
| `unpaired` | none of these |

**Support.** Support = (canonical + non-canonical cWW) / eligible chains. The canonical share is reported separately.

**Per-chain agreement.** Agreement is the same ratio over the chain's eligible SS_cons pairs.

### Recurrent 3D pairs not in SS_cons

**Selection.** From `interactions.json` `pairs`, base pairs with `in_ss_cons = false`, frequency ≥ 0.5 and ≥ 3 eligible chains.

**Kinds:**
- `pseudoknot (crosses SS_cons)`;
- `helix_extension`: a cWW stacked on an SS_cons pair;
- `shift`: shares a position with an SS_cons pair of support < 0.5;
- `tertiary (nested)`: everything else.

### Caveats

1. **Circularity.** Chains are placed on Rfam positions by cmalign with the family CM, which is built on SS_cons.
   - The alignment is therefore biased towards SS_cons, so agreement here is an upper bound.
   - Disagreement is more informative than agreement.
   - A wrong SS_cons pair can also pull an alignment column. Section 4 lists the candidates.
2. **Near pairs are not available.** The per-chain `near` lists contain base–ribose, base–phosphate and stacking near-interactions, never near base pairs. So `near_only` is always 0, which does not mean there were none.
3. **Missing observations.** "Unpaired" means no FR3D base pair between the two consensus positions. It can also be caused by:
   - a local alignment shift that the CM forced;
   - modified bases that FR3D does not pair;
   - low resolution: many RF00011 and RF00177 chains are cryo-EM.
4. **Family composition.**
   - RF00177 is 60 selected 16S chains, mostly from E. coli and T. thermophilus.
   - RF00011 is 50 chains of one RNA (Gst P RNA).
   - Per-pair support therefore mostly reflects a few species, not the whole family.
5. **Rfam SS_cons is a consensus.** For some families it was built without 3D data, and for 16S it is a large alignment-derived model. Disagreement is not by itself an error in Rfam.

## 2. Summary per family

| Family | SS_cons pairs | Chains | Pairs with support ≥ 0.9 | Pairs with support < 0.5 | cWW share of all pair observations | Median chain agreement | Chains with agreement < 0.8 | 3D-only recurrent pairs (pk / tertiary / extension / shift) |
|---|---|---|---|---|---|---|---|---|
| RF00167 purine | 22 | 84 | 18 | 0 | 0.947 | 0.95 | 2 | 7 / 9 / 1 / 0 |
| RF00162 SAM-I | 33 | 29 | 24 | 4 | 0.900 | 0.926 | 2 | 8 / 3 / 0 / 1 |
| RF00059 TPP | 28 | 38 | 19 | 3 | 0.877 | 0.895 | 3 | 8 / 3 / 1 / 2 |
| RF00005 tRNA | 21 | 60 | 17 | 0 | 0.916 | 0.952 | 7 | 5 / 1 / 1 / 0 |
| RF00011 RNase P B | 97 | 50 | 72 | 8 | 0.865 | 0.887 | 8 | 20 / 22 / 3 / 1 |
| RF00177 16S | 462 | 60 | 375 | 32 | 0.898 | 0.921 | 7 | 81 / 99 / 10 / 21 |

**Overall agreement.**
- For the riboswitches and tRNA, SS_cons agrees with the structures for about 90–95% of pair observations.
- The pairs that fail are concentrated, not scattered:
  - helix-terminal pairs that form a non-Watson–Crick pair instead;
  - a few register shifts;
  - one 16S hairpin that is often unformed.
- Isolated (lone) SS_cons pairs are not the problem: they have median support 0.97–1.0, except TPP (0.73, one pair).

## 3. SS_cons pairs not formed as cWW in 3D (support < 0.5)

### RF00167 purine and RF00005 tRNA

None.

### RF00162 SAM-I

| Pair | Support | What the structures show |
|---|---|---|
| 46–69 (AU in all 29) | 0.00 | tSH in 26 of 29 chains: a kink-turn-type sheared pair, not Watson–Crick. |
| 50–64 (AA in 25 of 28) | 0.11 | Not paired. Position 50 pairs cWW with 66 in about 96% of chains. A likely register difference at this pair of SS_cons. |
| 53–61 | 0.33 | Only 3 eligible chains. |
| 86–92 | 0.00 | Only 1 eligible chain. |

### RF00059 TPP

| Pair | Support | What the structures show |
|---|---|---|
| 20–28 and 21–27 | 0.29, 0.27 | About 30% of chains are shifted: 27 pairs with 22 (frequency 0.55), and the two chains disagree. The register in this short stem varies. |
| 64–83 (GU / GC) | 0.00 | 64 pairs with 82 (cWH, 0.94): a one-position shift in the 3D pairing. |

### RF00011 RNase P B (all 50 chains are Gst P RNA)

| Pair | Support | What the structures show |
|---|---|---|
| 7–361 GU | 0.00 | tSH / cWS. |
| 129–162 GA | 0.00 | tSH (0.68). |
| 290–305 AU | 0.38 | cWH (0.98 when not cWW). |
| 324–342, 325–341, 326–340 (G·A / A·G tandem) | 0.00 | Sheared G·A tandem (tSH), as expected for a G·A motif. Rfam writes these as Watson–Crick brackets. |
| 110–115 UG, 141–150 GU | 0.38, 0.30 | Mostly unpaired. |

### RF00177 16S (32 pairs)

**Unformed hairpin.**
- The stem 75–86 / 76–85 / 77–84 / 78–83 forms in only 23–49% of chains.
- It is the same block in every species; check whether this is a cryo-EM resolution effect.

**Consistent register shifts:**
- 372–385: 372 pairs with 386 in all 60 chains.
- 648–737: shifted in 57 of 59 chains.
- 649–736: shifted in 30 of 59 chains.

**Formed as a non-cWW pair:**

| Pair | Observed class |
|---|---|
| 139–171 GG | cWH in 55 of 60 |
| 680–687 GU | tSH-type, 46 of 59 |
| 712–723 GG | 53 of 59 |
| 986–1035 AC | 26 of 59 |

**Never paired, with the position pairing elsewhere:**

| Pair | What position i or j pairs with instead |
|---|---|
| 13–21 | 13 is cWH with 905 (56 of 60) and tHW with 20 |
| 37–537 | 37 is tWW with 394 (58 of 60) |
| 572–748 | 572 is tWH with 749 |
| 934–1228 | 1228 is cWW with 1328 (58 of 59) |
| 1043–1196 | 1043 is cWH with 1047 |
| 1149–1169 | 1149 is cWW with 1172 (57 of 59) |
| 1497–1520 | 1520 is cWW with 807 (59 of 59) |

**No partner observed at all:**
- 596–622;
- 950–965, 951–964;
- 1015–1026;
- 1054–1183;
- 1150–1168.

**Interpretation.** These are candidates for:
- an SS_cons that differs from the E. coli / T. thermophilus 3D structure;
- or local alignment columns placed by the CM against the 3D pairing.

They are not proven Rfam errors; see caveat 1. Distinguishing the two would need a structure-aware alignment check for the few chains concerned.

## 4. Recurrent 3D pairs absent from SS_cons (frequency ≥ 0.5)

Every family has some. Most are the expected tertiary structure, which a nested secondary structure cannot contain.

**RF00167 purine.**
- The L2–L3 kissing-loop pseudoknot:
  - 41–66 cWW (80 of 80);
  - 42–65 cWW (75 of 80);
  - 37–71 and 38–70 tWH.
- The J1/2–J2/3 base triples around the ligand:
  - 26–56 cWW (0.91);
  - 27–50 tWS;
  - 54–81 cWS;
  - 50–57.
- A P1 extension 14–92 … 17–89, present only in the 4–8 scaffold chains.

**RF00162 SAM-I.**
- The pseudoknot between L2 and J3/4: 25–80, 26–79, 27–78 and 28–77 cWW.
- 12–44, 12–71, 30–73 and 23–74.
- 50–66, the shift partner of 50–64.

**RF00005 tRNA.** The classic tertiary pairs:
- G18–Ψ55-type 18–54 cWW;
- the Levitt pair 15–46 tWW;
- 8–14 tWH, 21–45, 8–20, 52–56 tWH;
- 25–43 extending the anticodon stem.

**RF00059 TPP.**
- 32–36 tWH and 23–26 cWW.
- 64–82 cWH, the shift partner of 64–83.
- Several tSS / cSS contacts crossing SS_cons.
- 15–51 G·A cWW, extending a helix.

**RF00011 and RF00177.**
- RF00011: 20 pseudoknot, 22 tertiary, 3 extension and 1 shift pairs.
- RF00177: 81, 99, 10 and 21 respectively.
- Examples are among the 16S pairs listed above (1228–1328, 1149–1172, 807–1520), which are classic long-range 16S pairs (e.g. the central pseudoknot region).

## 5. Chains with low agreement (< 0.8)

The failures are at helix ends and in specific helices, not across whole chains.

| Chain | Agreement | Where it fails | Likely reason |
|---|---|---|---|
| RF00167 9NBC A | 0.68 | P1, P2, P3 ends | Scaffold construct (RNA scaffold attached to the aptamer). |
| RF00167 9V50 A | 0.65 | P1, P2, P3 ends | Guanine-II riboswitch, a divergent class within the family. |
| RF00162 3GX6 / 3GX7 A | 0.79 | P2 13–42 (0 of 4 in 3GX6), 50–64 | U34 variant structures. |
| RF00059 2HOP A, 4NYA B, 4NYG A | 0.65–0.79 | 17–31, 61–86, 67–80 | E. coli thiM with analogue ligands; cause not checked. |
| RF00005 3QSY C | 0.43 | Acceptor, D and anticodon stems | Initiator tRNA in a translation factor complex. |
| RF00005 6AHU M | 0.48 | Acceptor, D and anticodon stems | Mature tRNA in human RNase P, cryo-EM. |
| RF00005 6ZM6 | 0.71 | — | Mitochondrial tRNA in the mitoribosome. |
| RF00005 1J2B | 0.76 | D stem 0 of 4 | tRNA in the archaeosine transglycosylase complex, known to remodel the D arm (λ-form). |
| RF00177 4V9C, 4V9O, 6LKQ, 7SSD | 0.71–0.79 | The 72–89 hairpin, 9–25, 32–542 | Mostly the often-unformed 75–86 hairpin. |

Low agreement therefore mostly marks real structural states, such as remodelled tRNA in complexes or a divergent riboswitch class, or construct/resolution effects. It does not mark mapping errors: none of these chains fails across all helices.

## 6. Conclusions

1. **The SS_cons ↔ 3D comparison works and is informative** with the data we already have. No new downloads are needed.
2. **For the riboswitches and tRNA, SS_cons is well supported.** The disagreements are biologically meaningful:
   - sheared or kink-turn pairs written as Watson–Crick (SAM-I 46–69, RNase P G·A tandem);
   - register shifts (TPP 20–28 / 21–27 and 64–83; SAM-I 50–64 vs 50–66; 16S 372–385).
3. **3D-only pairs are dominated by the expected tertiary structure:** kissing loops, pseudoknots and tRNA tertiary pairs. This is exactly what a consensus secondary structure cannot express.
4. **Per-chain agreement flags remodelled or non-standard states** (tRNA in complexes, Guanine-II, scaffolds), which is useful for users.

## 7. Recommendation for the app (not implemented)

**Add a consensus-structure layer to the Interactions tab.**
- Each SS_cons pair is drawn and coloured by support, with the numerator and denominator.
- The status breakdown (canonical, non-canonical cWW, other class, shifted, unpaired) appears in the tooltip.
- A toggle shows the recurrent 3D-only pairs with their kind.
- Per chain, "SS_cons pairs formed: x / y eligible" is shown in the displayed structure.

**Data and labelling.**
- This needs no new sources: the analysis runs on `interactions.json` and the per-chain files.
- Generating a small `ss_cons_support.json` per family would carry the provenance and schema.
- The circularity caveat should be stated in the info text.

**Before building it, decide:**
1. Treat non-canonical cWW (G·A, U·U) as support (as here), or show it separately.
2. The status threshold wording: say "not formed as Watson–Crick in x / y chains", never "Rfam error".
3. Whether 16S (462 pairs) needs a helix-level summary instead of per-pair display.
