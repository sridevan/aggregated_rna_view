# Scientific review of the Interactions results: TPP riboswitch (RF00059) and bacterial 16S rRNA (RF00177)

Date: 2026-10-19.

**Data.** `data/families/<ACC>/interactions.json` and `interactions_by_instance/`, the production datasets built by
`rna_explorer interactions`. All inputs were used read-only.

**Reproduction.**
```
python3 experiments/interactions_review/review.py
```
This writes every table cited below, in about 1 s. The shared helpers are in `common.py`. Status calls use the production definitions:
- **present:** a non-near FR3D annotation;
- **absent:** both endpoints eligible but not annotated;
- **not comparable:** otherwise.

## Summary

| | TPP riboswitch (RF00059) | 16S rRNA (RF00177) |
|---|---|---|
| Chains / entries | 38 / 28 (E. coli thiM ×27, A. thaliana thiC ×8, scaffold / intron-embedded cryo-EM ×2, apo ×2) | 60 / 58 (E. coli; 70S, 30S and translocation complexes) |
| Mapping | every endpoint verified against the mmCIF | every endpoint verified against the mmCIF (688,932) |
| Consensus helices | all Rfam helices present at 97–100%, except an optional hairpin absent from every structure and one pair that differs in 3D (§1.2) | 441 of 462 SS_cons pairs observed; median chain helix coverage 95.5% |
| Literature check | apo 8F4O vs holo reproduces the Lee et al. 2023 changes (L5–P3 contact lost, P5 uncoiled, P3 register shift, J2/4 disordered) | central pseudoknot, 530-loop pseudoknot and 570/866 pairs in ≥ 97% of chains; the decoding-centre switch is recovered qualitatively |
| Main source of variation | sequence (E. coli vs A. thaliana); ligand state (apo vs holo) | functional state (decoding centre); candidate head-rotation-associated contacts |
| Data problems | sparse annotation in 2HOP (3.3 Å) and 4NYG (3.05 Å) | truncated BGSU annotation of 7PJS chain a (stops near nt 650); **fixed** in the pipeline (§2.1) |

**Overall.** The interaction data are scientifically sound for both families. Conserved architecture is recovered at high frequency, and
the known conformational changes appear where the literature places them. Most variation follows sequence or a documented functional
state, not noise. One data defect was found and has been fixed (§2.1). Three interpretation caveats should be visible to users: annotation coverage,
crystal-form effects, and resolution confounding.

---

## 1. TPP riboswitch (RF00059)

### 1.1 Data integrity

**Coverage** (`tpp_chain_coverage.csv`). The median consensus helix coverage is 91%.
- **Lowest:** 2HOP (pyrithiamine, 3.3 Å, 13/20), 4NYA B (17/22), 4NYG (thiamine, 3.05 Å, 15/19).
- **Not evaluable:** 9C6J (intron-embedded, 13 of 94 positions modelled) has no evaluable helix pair.
- **Truncation:** no chain has a truncated annotation (`tpp_truncation.csv` is empty).

**Mapping.** All 23,808 endpoints are verified against the mmCIF (`interactions.json` summary).

### 1.2 Conserved architecture (`tpp_consensus_pairs.csv`)

Consensus stems are named S1–S7 by order along the Rfam consensus structure, with 2GDI (E. coli thiM, Serganov et al. 2006) numbering.

**Present at 97–100% of eligible chains**
- **P1:** S1, G11–U88 … U14–A85.
- **P2:** S2, C15–G51 … G18–C48.
- **Lower P3:** S4, G21–C38, C22–G37, C23–G36.
- **P4:** S6.
- **P5:** S7.

**Exceptions, each with an explanation**
- **S5 (Rfam 38–48, an optional hairpin in the consensus):** absent from every structure. No chain has residues there; all structured TPP
  riboswitches lack this extension. This reflects consensus variability, not missing annotation.
- **Rfam pair G60–U79 (S6):** never annotated (0 of 23 evaluable chains). Instead, G60–G78 (cWH) is present in 33 of 35 chains. The
  structures have U79 kinked out to form the pyrophosphate pocket, with G60 pairing G78 (Lee et al. 2023). Here Rfam's consensus pairing
  differs from the 3D structure.
- **Upper P3 (C24–A35, U25–A34):** only 10/35 and 9/34. This is construct-dependent; see §1.4.

**Conserved tertiary contacts.** These non-SS_cons base pairs are present in ≥ 80% of eligible chains:
- U39–A43 tWH (J3/2 turn);
- G19–G42 tSS;
- G60–G78 cWH;
- G17–A56 tSS;
- C50–A84 tSS and A56–G83 tHS (the J2/4 A-minor tetrad region);
- G16–A84 cSS.

These are the junction contacts that hold the two sensor arms around the three-way junction.

### 1.3 Apo vs holo: test of Lee et al. 2023 (`tpp_apo_vs_holo.csv`)

The comparison uses the 16 holo chains with exactly the 8F4O sequence (2HO\*, 4NYB/C/D/G, 7TD\*, 7TZT/U) against apo 8F4O A and B, in 8F4O
numbering.

| Reported change (Lee et al. 2023) | Data |
|---|---|
| L5–P3 tertiary interaction lost | A69 base–ribose contacts to C23 (13/16 holo) and C38 (9/16), C22–A70 cSS (11/16), G21–U71 cSS (9/16) and the A69/A70 stack (16/16) are all **absent in apo A**. Apo B: not comparable (L5 not modelled). |
| P5 uncoils: C63•G76 and G60•G78 lost; C59-region pair kept | C63–G76 (13/14 holo) and G60–G78 (14/14) are **absent in both apo chains**; U59–A80 (15/16) is **present in both**. |
| P3 register shift: U25•A35 → U25•A34; C24 pairs A35; U26•G33 wobble | Apo A has **U25–A34, C24–A35 and U26–G33 cWW**; none occurs in holo chains of this sequence. |
| P2 and J2/4 tetrad partly kept; J2/4 (A53•A84) transient | G16–C50, G17–C49 and A56–G83 are present in both apo chains. A53–A84 (8/15 holo) is not comparable in apo A (A53 unmodelled) and absent in apo B. |
| P1 partially disordered | Apo A has positions 10–11 unmodelled ("M"), so they are not comparable rather than absent. |

**Every change reported in the paper that our data can test is recovered from the FR3D annotations.** Where the paper reports disorder,
our pipeline correctly calls the affected interactions "not comparable", never "absent".

### 1.4 Crystal-form effects in identical sequences

Two E. coli thiM constructs share the P3–L5 sequence and differ only at the P1 ends:
- **51ac69:** 2GDI, 4NYA, 7TZR/S;
- **acac3e:** 2HO\*, 4NYB/C/D/G, 7TD\*, 7TZT/U, 8F4O.

Even so, the upper P3 pairs C24–A35 and U25–A35 and the A69/C24 stack are annotated only in 51ac69 chains, and never in acac3e holo
chains. In acac3e, the L5–P3 contact is made through A69 base–ribose contacts instead. This is a crystal-form or construct effect in how
upper P3 and L5 pack. It is not ligand-related and not a mapping error. Comparisons across the two constructs will show these as
"only in" differences.

### 1.5 Variation among identical-sequence holo chains (`tpp_holo_variable.csv`)

Seven base pairs vary among the 16 same-sequence holo chains: G16–A84, G21–U71, C22–A70, C49–A56, A53–A84, G65–C74, G66–C73. All lie in
the L5–P3 contact, J2/4, or P5, the regions that close around the ligand.
- **Low-coverage chains:** the chains lacking them include the two sparsely annotated ones, 2HOP (pyrithiamine, 3.3 Å) and 4NYG
  (thiamine without pyrophosphate, 3.05 Å). Only those two lack the P5 pairs G65–C74 and G66–C73.
- **Well-resolved chains:** absences also occur in some well-resolved chains with complete coverage (e.g. 2HOJ, TPP, 2.5 Å, lacks
  G21–U71 and A53–A84).
- **Conclusion:** ligand chemistry (pyrophosphate-less or analogue ligands) and resolution are confounded here. These pairs should be read
  as dynamic contacts, not as ligand signatures.

### 1.6 E. coli vs A. thaliana (`tpp_species_differences.csv`)

There are 40 interactions whose frequency differs by ≥ 50 points between the species (≥ 4 eligible chains each). Of these:
- 34 have an endpoint whose nucleotide differs between the species;
- 38 lie within one nucleotide of a sequence difference;
- the remaining two are G21–U71 (13/26 vs 8/8) and a base–phosphate contact.

The species differences are therefore subclass and sequence effects (P3/L3, J3/2–P2, J2/4, P1 ends), consistent with the RF00167
finding that interaction classes follow sequence.

---

## 2. 16S rRNA (RF00177)

### 2.1 Data defect: truncated BGSU annotation of 7PJS chain a

**The problem**
- **Source annotation:** BGSU's FR3D files for 7PJS have only 526 base-pair rows (263 pairs) for chain a (16S), against about 1,150–1,330
  rows for other 16S chains. All rows lie between nucleotides 1 and about 651, although the chain is fully modelled (1,514 eligible
  positions, cryo-EM 2.35 Å).
- **What the pipeline does with it:** it treats the unannotated nucleotides as eligible. So 7PJS is counted **absent** for 2,343
  interactions (`16s_truncation_impact.csv`).
- **Effect on frequencies:** they fall by a median of 1.5 points (at most 3.6 points).
- **Effect on 7PJS itself:** its helix coverage is 38%, and any comparison involving 7PJS reports false "only in" differences for the
  3′ two thirds.

**Scope.** No other chain in any of the five families shows this pattern (`*_truncation.csv`).

**Fixed (pipeline schema 1.3, 2026-10-19).** `rna_explorer interactions` now records, per chain, `annotation_gaps`: runs of consecutive
eligible positions with no BGSU annotation of any category that are at least 15 positions long and at least 8% of the chain's eligible
positions.
- **Gap found:** RF00177 7PJS chain a, Rfam 652–1531. No other chain in the five families has one.
- **Effect on the family:** interactions touching a gap count as "annotation unavailable" for that chain. 2,344 RF00177 pairs moved 7PJS
  from eligible to unavailable (e.g. U13–A915 went from 56/60 to 56/59). The other families' pairs are byte-identical.
- **In the tab:** these interactions are "Not comparable — BGSU annotation incomplete for this chain (no annotation at Rfam 652–1531)".
- **Coverage:** 7PJS now shows 169 of 177 evaluable helix pairs (95.5%) instead of 38%.
- **Tests:** pytest (gap rule, and only 7PJS affected) and vitest (status and reason).

### 2.2 Coverage

- **Median:** consensus helix coverage is 95.5%.
- **Excluding 7PJS, the lowest** are 4V9C (3.3 Å, 74%), 8UZ3 (3.2 Å, 75%) and 4V9O PG (2.9 Å, 77%).
- **Effect on group 2:** sparsely annotated chains are slightly over-represented in the head-rotated group (median 93.7% vs 96.1%).

### 2.3 Conserved architecture (`16s_known_tertiary.csv`, `16s_conserved_long_range.csv`)

**Classic long-range interactions.** These are recovered in ≥ 97% of eligible chains (59 chains, 7PJS excluded):

| Interaction | E. coli numbering | Chains |
|---|---|---|
| Central pseudoknot | U17–A918, C18–G917, A19–U916 | 58/59, 59/59, 57/59 |
| 530-loop pseudoknot | G505–C526, G506–C525, C507–G524 (Woese & Gutell 1989; Powers & Noller 1991) | 57/59, 58/59, 59/59 |
| 570/866 pairs | G570–C866, U571–A865 | 58/59, 59/59 |

These are long-range or pseudoknotted, so the tab correctly classes them as tertiary (they are not in Rfam's nested SS_cons).

**Further conserved contacts.** 35 long-range (|i − j| ≥ 50) non-SS_cons base pairs are present in ≥ 90% of chains, mostly A-minor and
sheared contacts. Examples: A26–G558 tHS, G21–A914 tSH, G299–G566 cWH, G922–A1398 tSS. Two need comment:
- **C1407–G1494:** a canonical h44 pair, but outside Rfam's SS_cons, because the consensus leaves that position unpaired. C1407 is m5C,
  correctly mapped.
- **C817–G1529:** annotated in every chain and present in the raw BGSU file (e.g. 7K00). I found no literature description, so it is
  reported as observed but unverified.

### 2.4 Decoding centre (`16s_decoding_centre.csv`)

A1492 and A1493 of h44 switch between an intrahelical state and an extruded state that monitors the codon–anticodon helix together with
G530 (Ogle et al. 2001). Aminoglycosides bound to h44 favour the extruded state. The annotations separate the states qualitatively:
- **Extruded-state contact A1492 → C518 ribose:** only in chains with an A-site tRNA or an h44-binding aminoglycoside, namely 7B5K, 8BF7,
  8SYL (amikacin), 9B4Z (arbekacin) and 9B51 (G418).
- **Intrahelical contacts A1408–A1492 cWW and A1493–G1494 cSH:** only in chains with an empty A site (8B0X, and the tetracycline complexes
  9PIH and 9SX2, plus A1493–G1494 in 9PII and 9PIJ).
- **The A1492/A1493 stack:** present in 38 of 57 evaluable chains in both states, so it does not discriminate.

**Limits.** Several A-site-tRNA chains (8R6C, 9TMI, 9O2Y, 9QQQ) show neither state's diagnostic contact, so the interaction data are an
indicator, not a classifier, of decoding state. The two 30S mRNA-delivery complexes (9GUP, 9GUV) are not comparable at h44. Our 16S
dataset contains intra-16S contacts only; codon–anticodon contacts are inter-chain and excluded.

### 2.5 Structural groups: head rotation (`16s_group_differences.csv`)

Group 2 (10 chains: EF-G translocation intermediates, chimeric states, an open-head 30S) vs group 1 (49): 10 base pairs differ by ≥ 40
points. All are noncanonical minor-groove or sheared contacts that are less frequent in group 2:
- **Head domain:** G1089–A1170, A1093–C1109, U1307–U1330, A996–C1045.
- **Body, 530-loop and h18 regions:** G404–A498, C503–A510, A510–G542.
- **Central pseudoknot surroundings:** A16–A1080 and G15–U920.

The differences persist when both groups are restricted to chains with ≥ 90% helix coverage (e.g. G1089–A1170 38/41 vs 2/7; C503–A510
34/41 vs 1/7). They are candidate conformational differences associated with head swivel and translocation. However, group 2 has only 7
well-annotated chains, a lower median resolution (3.15 vs 2.8 Å), and these minor-groove contacts are resolution-sensitive. Targeted
geometric checks are needed before biological claims.

---

## 3. Recommendations

1. ~~Fix the 7PJS truncation (§2.1)~~: **done** (annotation gaps, schema 1.3).
2. **Done (§13.11 of the validation file): annotation context wherever "only in" or "absent" is interpreted.** The coverage line in Compare two already does this; adding
   resolution next to coverage would help, since sparse annotation (2HOP, 4NYG, 4V9C, 8UZ3) explains part of the variation.
3. **Done (§13.11): warn about crystal-form or construct effects** when two chains share a sequence but come from different constructs (TPP §1.4). This
   could be a note when the two compared chains have different `sequence_key`s, or the same sequence from different crystal forms.
4. **Near base pairs (fr3d-python)** would help most in the 16S head and in TPP L5/J2/4, where strained minor-groove contacts drive the
   variation.

## 4. Sources

- Lee H-K, …, Wang Y-X, Stagno JR (2023). *Structure* 31:848–859. Crystal structure of *Escherichia coli* thiamine pyrophosphate-sensing
  riboswitch in the apo state. doi:10.1016/j.str.2023.05.003; open access PMC10335363. PDB 8F4O.
- Serganov A, Polonskaia A, Phan AT, Breaker RR, Patel DJ (2006). *Nature* 441:1167–1171. doi:10.1038/nature04740. PDB 2GDI.
- Ogle JM, Brodersen DE, Clemons WM Jr, Tarry MJ, Carter AP, Ramakrishnan V (2001). *Science* 292:897–902. Recognition of cognate transfer
  RNA by the 30S ribosomal subunit. doi:10.1126/science.1060612.
- Powers T, Noller HF (1991). *EMBO J* 10:2203–2214. A functional pseudoknot in 16S ribosomal RNA. PMC452909. The 505–507/524–526 pairing
  was proposed by Woese & Gutell (1989).
