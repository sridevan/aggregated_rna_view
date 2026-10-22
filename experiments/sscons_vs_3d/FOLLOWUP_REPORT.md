# Follow-up: Rfam SS_cons versus observed 3D base pairs (validation, 2026-10-22)

> **Update (2026-10-22): see `NEAR_CHECK.md`.** Near base pairs are not unavailable: BGSU's base-pair export omits them, but fr3d-python produces them with `-c basepair,near`. They cover about half of the SS_cons observations without a true pair, and most of the "WC-like" ones in section 3.

**What this is.** A validation of `REPORT.md` / `SUMMARY.md` before any consensus-structure layer is built. This is investigation only: the application, backend and data are unchanged.

**Scripts.** All outputs are written to `followup/`.

| Script | What it produces |
|---|---|
| `followup.py` | Three support measures, weighting, alternative partners, completeness, outliers and classification attributes, from existing data only |
| `geometry.py` | Watson–Crick-edge geometry of every SS_cons pair observation without an FR3D base pair, in five families. PDBe updated mmCIF is streamed into memory and never stored. |
| `geometry_16s_hairpin.py` | The same geometry measure for the 16S 75–86 hairpin, in 7 chains |

**Inputs.** Unchanged:
- `family.json`, `interactions.json` and `interactions_by_instance/`;
- `instances.json`;
- the validated mappings;
- the cmalign output `work/aligned.sto`;
- cached BGSU FR3D CSVs (Rfam 15.1).

**Positions** are Rfam consensus positions. Residues are given as name + author number of the chain.

## Corrections to the first report

| Statement in REPORT.md / SUMMARY.md | Corrected |
|---|---|
| 16S is "mostly E. coli and T. thermophilus" | **All 59 chains with an organism are E. coli** (1 chain has none). The 15 distinct sequences fall into one 95%-identity cluster. RF00011 is likewise one RNA (Gst P RNA, one cluster). |
| RNase P 324–342 / 325–341 / 326–340 form a "sheared G·A tandem" | Only **324–342** is mostly sheared (tSH in 37 of 50 chains). 325–341 and 326–340 have no FR3D pair in 40 and 43 of 50 chains, and their bases are > 5 Å apart in every chain measured. |
| TPP 64–83 is "a one-position shift" | No Watson–Crick register shift. G at 64 forms a **cWH** pair with G at 82, and the nucleotide at 83 is unannotated or incomplete in most chains (section 2). |
| Low agreement of 3QSY C and 6AHU M reflects remodelled tRNA in complexes | Mostly **FR3D annotation sensitivity**: the stems are geometrically paired but not annotated (section 5). The exception is the C1·A72 mismatch of 3QSY. |
| SAM-I 50–64 vs 50–66 is a "register difference" | An **alignment-placement (correspondence) problem**: posterior probability is 0.05–0.15 at 64 (section 2). |

## 1. Revised support: three measures

All three measures use the same evaluable chains: both positions eligible, and not in an annotation gap.

- **any:** any non-near FR3D base pair between i and j.
- **cWW:** a cWW pair.
- **canonical:** a cWW pair of AU, GC or GU.

| Family | Evaluable SS_cons pairs | Support ≥ 0.9 (any / cWW / canonical) | Support < 0.5 (any / cWW / canonical) | Pooled observations (any / cWW / canonical) |
|---|---|---|---|---|
| RF00167 purine | 22 | 18 / 18 / 17 | 0 / 0 / 0 | 0.953 / 0.948 / 0.921 |
| RF00162 SAM-I | 33 | 24 / 24 / 24 | 3 / 4 / 4 | 0.932 / 0.900 / 0.900 |
| RF00059 TPP | 23 (5 never evaluable) | 20 / 19 / 19 | 3 / 3 / 3 | 0.887 / 0.876 / 0.874 |
| RF00005 tRNA | 21 | 17 / 17 / 16 | 0 / 0 / 0 | 0.920 / 0.916 / 0.908 |
| RF00011 RNase P B | 97 | 73 / 72 / 72 | 4 / 8 / 8 | 0.898 / 0.865 / 0.865 |
| RF00177 16S | 455 (7 never evaluable) | 376 / 375 / 366 | 29 / 32 / 45 | 0.906 / 0.898 / 0.870 |

**How the measures differ.** They differ only for a few pairs. In 16S, 13 consensus pairs are consistently non-canonical cWW in E. coli (e.g. G·A, U·U). These count under cWW but not under canonical.

**SS_cons pairs whose partners interact, but not through Watson–Crick edges.** The expected partners are present in these cases; the geometry differs.

| Family | SS_cons pair | Nucleotides | Evaluable chains | Any FR3D pair | cWW | Observed classes |
|---|---|---|---|---|---|---|
| RF00162 | 46–69 | A·U (all) | 29 | 0.90 | 0.00 | tSH 26, none 3 |
| RF00011 | 324–342 | G·A | 50 | 0.74 | 0.00 | tSH 37, none 13 |
| RF00011 | 129–162 | G·A | 50 | 0.68 | 0.00 | tSH 33, tSW 1, none 16 |
| RF00011 | 7–361 | G·U | 50 | 0.54 | 0.00 | cWS 27, none 23 |
| RF00011 | 290–305 | A·U | 50 | 0.98 | 0.38 | cWH 30, cWW 19 |
| RF00177 | 139–171 | G·G | 60 | 0.92 | 0.00 | cWH 55 |
| RF00177 | 712–723 | G·G | 59 | 0.90 | 0.00 | tWW 53 |
| RF00177 | 680–687 | G·U | 59 | 0.78 | 0.00 | tSH 46 |
| RF00011 | 325–341 / 326–340 | G·A / A·G | 50 | 0.20 / 0.14 | 0.00 | mostly none; bases > 5 Å apart (40 of 40 and 42 of 43 measured) |

**Conclusion for SAM-I 46–69 and the RNase P sheared pairs.** The low cWW support reflects a different interaction geometry, not missing partners:
- For SAM-I 46–69, the A and U are paired, as a trans Sugar-edge / Hoogsteen pair, in 90% of chains.
- In RNase P, the same holds only for 324–342 and 129–162. 325–341 and 326–340 are genuinely not paired in these structures.

## 2. Alternative pairing partners

For each case, all chains were checked, using:
- residues and author numbers from the mapping;
- cmalign posterior probabilities (PP);
- insertions inside the span;
- the base-pair network of i, j and the alternative partner;
- the local cmalign rows (`followup/alt_partner_alignment.txt`).

**PP codes** (Infernal): `*` ≥ 0.95, `9` 0.85–0.95, `7` 0.65–0.75, `4` 0.35–0.45, `1` 0.05–0.15. High PP means the CM is confident; it does not show that the alignment is correct.

| Case | Representative chains | Residues (cons i / j / alt) | PP (i j alt) | SS_cons pair | Alternative pair | Verdict |
|---|---|---|---|---|---|---|
| SAM-I 50–64 vs 50–66 | 2GIS A, 3GX5 A and 23 other chains with the same construct sequence (T. tengcongensis, recorded as Caldanaerobacter subterraneus) | A49 / A52 / U54 | 4 1 7 | not formed (0 of 25); A·A is not complementary | A49–U54 cWW in 25 of 25 | **Correspondence problem.** Consensus 51–63 is almost entirely deleted in this sequence (`AACCA----`, `a---AAUGG`), and PP at 64 is 0.05–0.15. Three chains with confident alignments (PP 8–\*) form 50–64 cWW: 3V7E C/D, a T. tengcongensis construct with A49–U86 at 50–64, and 4KQY A, B. subtilis G50–C73. The observed A49–U54 pair is real; its consensus assignment is uncertain. |
| TPP 64–83 vs 64–82 | 2GDI Y, 4NYA A (E. coli); 2CKY A, 3D2G A (A. thaliana) | E. coli G60 / U79 / G78; A. thaliana G48 / C67 / G66 | \*\*\* | 0 of 23 evaluable, including canonical G48–C67 in A. thaliana; bases > 5 Å apart in 23 of 23 | G·G **cWH** in 33 of 38 | **Observed structural alternative; not a correspondence problem.** Confident alignment, no insertions, two species agree. The nucleotide at 83 has no FR3D annotation in 22 chains and an incomplete base in 12. No Watson–Crick pair is formed in this register or the next. |
| 16S 372–385 vs 372–386 | 4YBB AA, 7K00 A, 4V9D AA (all 60 E. coli) | U375 / G388 / A389 | \*\*\* | 0 of 60 | U375–A389 cWW in 60 of 60; A389 also tWH with consensus 369 in 51 (a triple) | **Bulge-position ambiguity.** SS_cons pairs 372–385 (U·G) and places a bulged A (conserved in RF) at 386. E. coli forms U375–A389 and bulges G388. The two alignments are equally consistent with the sequence, and with one species the structures cannot show which is consensus-wide. The structural alternative is credible and reproducible. |
| 16S 648–737 vs 648–738 | the same E. coli chains | C658 / A747 / G748 | \*\*\* | 0 of 59; C·A is not complementary | C658–G748 cWW in 57 of 59; U659 pairs A747 (649–737) in 30 and A746 (649–736, SS_cons) in 23 | **Bulge-position ambiguity.** SS_cons puts the bulge at 738 (RF `u`), which E. coli fills with a paired G748. In E. coli the unpaired nucleotide is A746 or A747, and this varies between structures. The structural alternative is credible. |
| TPP 20–28 / 21–27 (vs 21–28 / 22–27) | 2GDI Y, 4NYG A (E. coli); 2CKY A (A. thaliana) | E. coli C24 U25 U26 / A34 A35 | \*\*9 | E. coli: C24·A35 is a mismatch, 0 formed; A. thaliana (U16–A23, G17–C22) forms the consensus pairs | E. coli U25–A35 (10 of 29) and U26–A34 (11 of 29) | **Correspondence ambiguity plus a partly unformed stem.** The CM keeps a C·A mismatch rather than a bulge, although the E. coli sequence fits the shifted register. In 25 chains the stem is not paired at all (bases apart). |

**Conclusion.**
- Two of the four requested cases are *bulge-placement ambiguities* in a single species: 16S 372 and 16S 648.
- One is a *low-confidence alignment placement*: SAM-I 50.
- One is a *genuine non-Watson–Crick alternative*: TPP 64.
- None is evidence that Rfam's consensus structure is wrong for the family. Resolving the 16S cases needs 3D structures from other species, such as T. thermophilus 30S (not in our 60-chain set).

## 3. Annotation completeness

**Near base pairs are unavailable at the source.**
- Not one of the 286 cached BGSU base-pair files contains an `n`-prefixed class.
- All 286 stacking files and all base-ribose files do contain near classes.
- So BGSU publishes no near base pairs. "Near only" cannot be assessed, and the absence of near annotations carries no biological meaning.

**The 15-position gap threshold is not hiding shorter gaps.** Runs of eligible consensus positions without any BGSU annotation, of any category:

| Family | Run length (number of runs) |
|---|---|
| RF00167 | 1 (77) |
| RF00162 | 1 (69) |
| RF00059 | 1 (68), 2 (2) |
| RF00005 | 1 (68), 2 (11), 3 (1) |
| RF00011 | 1 (177), 2 (4) |
| RF00177 | 1 (912), 2 (61), 3 (8), 4 (1), 5 (1), 879 (1: 7PJS, the known gap) |

- No runs of 6–14 positions exist.
- Only 2 SS_cons observations (16S) touch a run of 3–5.
- Single unannotated residues are usually flipped-out or disordered nucleotides, e.g. TPP E. coli U79 at 83.

**The real completeness issue is FR3D sensitivity inside annotated regions.**
- **Method.** For every SS_cons observation with no FR3D base pair (815 in the five non-16S families), I measured:
  - the Watson–Crick-edge N–N distance (purine N1 / pyrimidine N3);
  - the C1′–C1′ distance.
- **Calibration.** 98% of SS_cons pairs that FR3D *does* call cWW (1,121 of 1,146 in the same chains) meet "WC-like" (N–N ≤ 3.6 Å and C1′–C1′ 8.5–12.5 Å).
- **Result:**

  | Geometry of the 815 observations without an FR3D pair | Count |
  |---|---|
  | WC-like (paired, not annotated) | 176 (22%) |
  | Distorted (N–N ≤ 5 Å) | 248 (30%) |
  | Apart | 391 (48%) |

- **WC-like share of complementary pairs by resolution:**

  | Resolution | WC-like share |
  |---|---|
  | ≤ 2.5 Å | 10 of 50 |
  | 2.5–3.0 Å | 79 of 271 |
  | 3.0–3.5 Å | 70 of 235 |
  | > 3.5 Å | 10 of 45 |

  This is not a simple resolution effect.

- **Effect if WC-like observations were counted as formed:**

  | Family | Pooled cWW support | Pairs ≥ 0.9 |
  |---|---|---|
  | RF00167 | 0.948 → 0.963 | 18 → 19 |
  | RF00162 | 0.900 → 0.917 | 24 → 28 |
  | RF00059 | 0.876 → 0.884 | 19 → 19 |
  | RF00005 | 0.916 → 0.942 | 17 → 19 |
  | RF00011 | 0.865 → 0.885 | 72 → 74 |

  - The largest single change: RNase P 110–115 (U·G closing a 4-nt loop) goes from 0.38 to 0.96; 29 of its 31 "unpaired" chains are WC-like.

**Modified nucleotides.**
- They are eligible when the base ring is complete. SS_cons pairs with a modified end are cWW in 138 of 147 observations.
- **BGSU base-pair gap for pseudouridine.** In 9 of 62 entries containing Ψ, BGSU lists Ψ stacking but no Ψ base pair at all:
  - 1B23, 1EHZ, 3CW5, 8CTH, 8UU5, 9B4Z, 9SFQ, 9SRE, 9TI5.
  - In 1EHZ the 2′-O-methylated C32 and G34 and wybutosine 37 (OMC, OMG, YYG) are likewise absent from the base-pair file.
  - Consequence: tRNA 30–38 (A31–Ψ39) is "unpaired" in 1B23, 1EHZ, 1TTT, 8CTH and 8UU5. The frequency of the 3D-only G18–Ψ55 pair (18–54) is probably depressed.
  - Cause unknown; it is a source annotation gap, not absence.

**Missing residues and incomplete bases** are excluded from the denominators. Here are the SS_cons pair × chain combinations excluded:

| Family | Evaluable | Incomplete base | Unmodelled | Not mapped (deletion or local end) |
|---|---|---|---|---|
| RF00167 | 1,587 | 2 | 106 | 153 |
| RF00162 | 822 | 1 | 0 | 134 |
| RF00059 | 793 | 17 | 32 | 222 |
| RF00005 | 1,256 | 0 | 2 | 2 |
| RF00011 | 4,850 | 0 | 0 | 0 |
| RF00177 | 27,262 | 0 | 36 | 422 (and 7PJS gap) |

## 4. Dataset redundancy

Sequence clusters are greedy, at ≥ 95% identity over shared consensus positions.

| Family | Chains | Entries | Unique sequences | Sequence clusters (95%) | Species | Largest cluster (chains) |
|---|---|---|---|---|---|---|
| RF00167 | 84 | 66 | 29 | 15 | 6 (B. subtilis 38) | 23 |
| RF00162 | 29 | 28 | 23 | 6 | 3 (T. tengcongensis 27) | 13 |
| RF00059 | 38 | 28 | 5 | 4 | 3 (E. coli 29) | 18 |
| RF00005 | 60 | 59 | 43 | 38 | 14 | 5 |
| RF00011 | 50 | 50 | 4 | **1** | 1 | 50 |
| RF00177 | 60 | 58 | 15 | **1** | 1 (E. coli) | 60 |

**Instance-weighted vs cluster-weighted cWW support** (each sequence cluster has weight 1):
- The band counts barely change: pairs ≥ 0.9 move by at most 1–2 per family, and pairs < 0.5 are unchanged.
- Individual pairs can move a lot:

  | Pair | Instance-weighted | Cluster-weighted | Note |
  |---|---|---|---|
  | RF00167 64–72 | 0.889 | 0.641 | 18 chains, 3 clusters |
  | RF00162 83–95 | 0.897 | 0.653 | |
  | RF00162 50–64 | 0.107 | 0.233 | |
  | RF00059 16–49 | 0.730 | 0.611 | |
  | RF00005 30–38 | 0.850 | 0.900 | |

**For RF00011 and RF00177**, every support value describes **how reproducibly one RNA is modelled across structures and functional states**. It says nothing about evolutionary conservation. For the riboswitches, a single heavily studied sequence still dominates:
- B. subtilis xpt: 23 chains in one cluster;
- T. tengcongensis SAM-I: 27 of 29 chains;
- E. coli thiM: 29 of 38 chains.

## 5. Outliers

Unannotated SS_cons pairs were measured as in section 3.

| Chain | Context (from the entry) | Agreement | What the structure shows | Attribution |
|---|---|---|---|---|
| 3QSY C | E. coli initiator tRNA with archaeal aIF2; X-ray 3.2 Å | 0.43 | 9 cWW annotated. Of 12 unannotated SS_cons pairs, 9 are WC-like or nearly so (N–N 2.5–3.7 Å; e.g. G10–C25 2.5 Å, G4–C69 3.1 Å). C1·A72 (1–70): 24.8 Å apart. G2–C71: 7.9 Å. C28–G42: 15.5 Å. | **Mostly annotation sensitivity.** C1·A72 is a sequence feature: the C1·A72 mismatch of bacterial initiator tRNA is well documented. 2–69 and 27–41 are genuinely open in this model. |
| 6AHU M | Human mature tRNA in human RNase P; cryo-EM 3.66 Å | 0.48 | All 11 unannotated SS_cons pairs are complementary, N–N 3.0–5.1 Å, C1′–C1′ 10–12 Å. | **Annotation sensitivity**; the stems are present. No evidence for a structural difference. |
| 1J2B C | P. horikoshii tRNA(Val) with archaeosine tRNA-guanine transglycosylase; title "λ-form tRNA"; X-ray 3.3 Å | 0.76 | D stem 10–24 to 13–21: bases 13–20 Å apart. D-stem nucleotides pair with other positions (24↔45, 23↔insert 45^1, 22↔46). | **Genuine structural difference.** Consistent with the λ-form named in the entry, described by Ishitani et al. (Cell 2003). |
| 9V50 A | S. aureus Guanine-II riboswitch; X-ray 3.25 Å | 0.65 | Rfam bit score 43.7; only 67% of consensus positions matched; 18–88 unmapped. 59–77 is A·C (mismatch). The other unannotated pairs are complementary, N–N 2.9–3.9 Å. | **Sequence divergence** (a distinct riboswitch class per the entry title), **plus annotation sensitivity.** |
| 9NBC A | Synthetic 8-oxoguanine aptamer on an RNA scaffold; cryo-EM 2.87 Å | 0.68 | 18–88, 19–87 and 64–72 unmodelled. 35–43 is U·U. Other unannotated pairs are complementary but distorted (N–N 3.5–5.4 Å). | **Construct and sequence**, plus annotation sensitivity / distortion. |
| 16S hairpin 75–86 … 78–83 | E. coli G79 A80 A81 G82 / C87 U88 U89 C90; fully complementary in 57 of 60 | — | All 4 cWW in 7 chains, none in 18, unmodelled in 5. In 4YBB (X-ray 2.1 Å), 76–85, 77–84 and 78–83 are WC-like (2.8–3.2 Å) but unannotated. In 7K00 and 8AYE (cryo-EM 1.96–1.98 Å), 78–83 is open (7.0–7.6 Å) and 75–86 / 77–84 are distorted (4.2–4.9 Å). | **Both annotation sensitivity and genuine local variation** in the hairpin tip. Not explained by resolution. No functional interpretation is claimed. |

## 6. Classification of recurrent 3D-only pairs

**Which pairs were assessed.** 307 pairs not in SS_cons, with frequency ≥ 0.5 over ≥ 3 eligible chains. Each was given independent attributes (`followup/extra_pair_attributes.csv`).

**Problems with the first report's kinds:**
1. **They are not mutually exclusive.** 196 of 307 pairs have more than one of these attributes:
   - crosses SS_cons;
   - shares a nucleotide with an SS_cons pair;
   - stacks on an SS_cons pair;
   - both positions unpaired in SS_cons.

   The first report's fixed precedence hid this.
2. **"Pseudoknot" for any crossing pair is too broad.** Most crossing pairs are single tertiary contacts, often non-Watson–Crick, for example:
   - TPP 15–35 tSS;
   - tRNA Levitt pair 15–46 tWW and 8–14 tWH;
   - SAM-I 12–44 cSS.

   Only stacked runs of crossing cWW pairs are pseudoknot helices:
   - RF00167 41–66 / 42–65 (L2–L3 kissing loop);
   - RF00162 25–80 to 28–77.
3. **"Shift" (sharing a nucleotide with a poorly supported SS_cons pair) conflates base triples with replacements.** Per-chain co-occurrence separates them (counts are both formed / only the 3D pair formed):

   | Kind | 3D-only pair | SS_cons pair sharing a nucleotide | Both formed | Only the 3D pair formed |
   |---|---|---|---|---|
   | Base triple | RF00167 54–81 | 25–81 | 70 | 1 |
   | Base triple | RF00162 12–44 | 44–71 | 28 | 0 |
   | Base triple | RF00059 13–60 | 13–53 | 34 | 0 |
   | Replacement | SAM-I 50–66 | 50–64 | 0 | 25 |
   | Replacement | TPP 64–82 | 64–83 | 0 | 21 |
   | Replacement | TPP 22–27 | 21–27 | 0 | 11 |

   In 16S, 21 pairs are replacements and 63 are triples or bifurcations.
4. **"Helix extension" only looks one pair deep.** RF00167 14–92, 15–91 and 16–90 continue the extended P1 through 17–89. They are found in scaffold chains only (4–8 eligible). They were labelled "tertiary".
5. **"Tertiary (nested)" includes local pairs that SS_cons lacks**, such as cWW across a hairpin loop: TPP 23–26, RF00167 50–57.

**Proposed scheme.**
- Attributes stay separate in the data:
  - geometry class;
  - crosses SS_cons;
  - per-chain co-occurrence with any SS_cons pair sharing a nucleotide;
  - stacked continuation of an SS_cons helix;
  - crossing-helix membership (≥ 2 stacked crossing cWW pairs).
- One display label is derived from them, by stated precedence:
  1. **Replaces an SS_cons pair:** the 3D pair forms mostly in chains where the SS_cons pair does not.
  2. **Base triple / bifurcation with an SS_cons pair.**
  3. **Helix extension** (chained).
  4. **Crossing helix (pseudoknot).**
  5. **Crossing tertiary contact.**
  6. **Nested pair outside SS_cons helices.**
- The Leontis–Westhof class is always shown.

## 7. Recommendations for the consensus-structure layer (not implemented)

1. **Status per SS_cons pair and chain.** Replace "unpaired" with explicit, mutually exclusive states:
   - **cWW formed** (canonical or non-canonical);
   - **other FR3D pair** (with its class);
   - **alternative partner** (name the partner and the consensus position);
   - **no FR3D pair, WC-like geometry** ("paired in the model, not annotated by FR3D");
   - **no pair** (bases distorted or apart);
   - **not evaluable**: missing, incomplete, unmapped, annotation gap, or **base-pair annotation unavailable for this modified nucleotide**.

   The geometry state needs a small backend step: two distances per pair from the mmCIF the pipeline already reads for eligibility. It needs tests and a schema.
2. **Show the three measures** (any pair / cWW / canonical) with explicit numerators and denominators. Use cWW as the colour and the others in the tooltip.
3. **Carry redundancy with every support value:** n chains, n sequence clusters, n species.
   - Offer the cluster-weighted value as a secondary figure.
   - For single-cluster families (RF00011, RF00177), state "one RNA: reproducibility across structures, not conservation".
4. **Flag low alignment confidence.** Mark pairs where either end has PP < 0.5 in a chain (e.g. SAM-I 50–64), so that a correspondence problem is not shown as a structural observation.
5. **Show alternative partners only where co-occurrence shows replacement.** Word them as observations:
   - "In 25 of 28 chains, position 50 pairs with 66 instead; alignment confidence at 64 is low";
   - "E. coli forms U375–A389 and bulges G388; consensus bulge at 386".

   Never call a disagreement an Rfam error.
6. **Use the attribute-based classification of section 6** for 3D-only pairs. Restrict "pseudoknot" to crossing helices.
7. **Before release:**
   - check the Ψ annotation gap with BGSU (9 entries);
   - decide whether to include a structure from another species (e.g. T. thermophilus 30S) to test the 16S bulge cases.

## Files (`followup/`)

| File | Contents |
|---|---|
| `support_three_measures.csv` | Per-pair support under the three measures, plus sequence-, cluster- and species-weighted cWW support |
| `focus_pair_classes.csv` | Class distributions of the geometry examples |
| `alt_partner_chains.csv`, `alt_partner_alignment.txt` | The alternative-partner cases |
| `unpaired_causes.csv`, `unpaired_geometry.csv`, `unannotated_runs.csv` | Completeness diagnostics |
| `redundancy.csv` | Redundancy per family |
| `outlier_helices.csv`, `hairpin_16s.csv`, `hairpin_16s_geometry.csv` | Outliers |
| `extra_pair_attributes.csv` | Classification attributes of 3D-only pairs |
| `followup_summary.json` | Summary numbers |

**Species-weighted support** is in the CSV but is not meaningful for 16S, where one chain has no organism and gets half the weight.
