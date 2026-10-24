# Ribosome–drug inventory (2026-10-23)

**Question.** Is there enough ribosome–drug structural data, with usable PDBe contacts, to find new patterns (cross-resistance, shared hotspots, species selectivity)? Script: `inventory.py`. Output: `out/inventory.json`, one row per ligand copy.

## Sources
- **Entries:** PDBe search (`rfam_accession`) for four Rfam families.
  - **Bacterial:** RF00177 (SSU) and RF02541 (LSU). Both also contain the human mitochondrial 12S and 16S rRNAs.
  - **Eukaryotic:** RF01960 (SSU) and RF02543 (LSU).
  - **Result:** 2,948 entries with these rRNAs.
  - **Why the search index:** the PDBe API `nucleic_mappings` endpoint did not answer within 5 minutes for these families.
- **Drug classes:** assigned by hand to 248 chemical components (CCD IDs; `CLASSES` in the script). Ions, buffers, modified nucleotides, cofactors, polyamines and SAM/SAH are left out.
- **Contacts:** per ligand copy, from PDBe API `bound_molecules` and `bound_ligand_interactions`.
  - Contacts of type proximal, vdw, vdw_clash and clash are not counted.
  - "ok" means at least one contact with an rRNA chain of the four families.

## Counts
- **Overall:** 612 entries contain a classed drug, giving 1,330 ligand copies.
  - 1,179 copies are usable: 568 entries with contacts to rRNA.
  - 133 have no PDBe contacts (HTTP 404).
  - 18 contact only proteins or other RNAs.

| Class | Copies | Usable (entries) | SSU / LSU | Bacterial / eukaryotic / mitochondrial | ≤ 2.5 Å |
|---|---|---|---|---|---|
| Aminoglycoside (2-DOS) | 670 | 599 (146) | 314 / 285 | 268 / 322 / 9 | 132 |
| Tetracycline | 168 | 124 (49) | 44 / 80 | 53 / 55 / 16 | 41 |
| Eukaryotic inhibitors | 108 | 104 (91) | 8 / 96 | 5 / 99 / 0 | 16 |
| Macrolide/ketolide | 80 | 80 (78) | 0 / 80 | 79 / 1 / 0 | 18 |
| Other PTC/tunnel drugs | 51 | 49 (47) | 10 / 39 | 30 / 19 / 0 | 9 |
| Other SSU drugs | 51 | 37 (17) | 10 / 27 | 35 / 2 / 0 | 26 |
| Streptogramin | 31 | 31 (31) | 0 / 31 | 22 / 0 / 9 | 9 |
| Phenicol | 28 | 25 (25) | 0 / 25 | 25 / 0 / 0 | 4 |
| Spectinomycin | 26 | 20 (18) | 18 / 2 | 20 / 0 / 0 | 12 |
| Oxazolidinone | 20 | 20 (20) | 0 / 20 | 19 / 0 / 1 | 8 |
| Streptomycin | 19 | 19 (19) | 19 / 0 | 14 / 0 / 5 | 0 |
| Hygromycin B | 19 | 19 (17) | 19 / 0 | 7 / 12 / 0 | 8 |
| Lincosamide | 19 | 16 (16) | 0 / 16 | 16 / 0 / 0 | 6 |
| Pleuromutilin | 18 | 16 (16) | 0 / 16 | 15 / 0 / 1 | 2 |
| Negamycin, kasugamycin, orthosomycin | 22 | 20 | | bacterial | |

- **Species (usable copies):** Thermus thermophilus 326, human 179, S. cerevisiae 163, E. coli 152, Candida albicans 104, C. auris 55, Acinetobacter baumannii 34, Deinococcus radiodurans 25, Leishmania 16.
- **Missing PDBe contacts:** 126 of the 133 copies with no contacts are cryo-EM structures released from 2021 on (for example 7M4U–7M4Y eravacycline, 7P7Q/7P7R spectinomycin). This is a PDBe coverage gap worth raising internally.
- **Copies on non-rRNA chains** (LSU counts for aminoglycosides and tetracyclines): these are partly secondary sites, for example H69 and other 23S sites. They need site assignment on Rfam positions before any per-site claim.

## First look: contacts with the nucleotides modified by resistance enzymes
E. coli and Thermus copies only, since both use E. coli numbering (498 copies). Each cell is copies with a non-vdW contact **to the modified atom** / **to the nucleotide**, out of n.

| Class | n | A2058 N6 (Erm) | A2503 C8 (Cfr) | A1408 N1 (NpmA) | G1405 N7 (ArmA, RmtB) |
|---|---|---|---|---|---|
| Aminoglycoside (2-DOS) | 261 | – | – | 76 / 105 | 18 / 103 |
| Macrolide/ketolide | 47 | 10 / 46 | 31 / 41 | – | – |
| Lincosamide | 11 | 7 / 10 | 9 / 11 | – | – |
| Phenicol | 11 | 0 / 1 | 6 / 9 | – | – |
| Pleuromutilin | 9 | – | 9 / 9 | – | – |
| Streptogramin (mostly group A) | 21 | – | 0 / 21 | – | – |
| Oxazolidinone | 9 | – | 0 / 2 | – | – |

**What this already shows**
- **Partial match with the resistance phenotypes:**
  - Erm (MLS_B group): macrolides and lincosamides contact A2058.
  - Cfr (PhLOPSA group): phenicols, lincosamides and pleuromutilins contact A2503 C8.
- **Contacts alone don't predict resistance:**
  - Most macrolides contact A2503 C8 (31 of 41), yet Cfr does not confer resistance to 14- and 15-membered macrolides.
  - Cfr does confer resistance to oxazolidinones and group A streptogramins, yet these contact C8 in 0 copies.
  - Erm dimethylates N6, but only 10 of 46 macrolide copies contact N6 itself.
- **Implication:** a useful analysis has to model the added methyl group (steric clash or displacement within a few ångström of the modified atom) rather than count contacts. That is a testable, structure-wide question with known resistance spectra as the answer key.

## Data are sufficient for
1. **Cross-resistance from structure.** Model the methyl groups added by Erm, Cfr, NpmA, ArmA and RmtB on every copy and compare the predicted clashes with the known resistance spectra: MLS_B, PhLOPSA, and the aminoglycoside methyltransferases' spectra.
2. **Species selectivity.**
   - Aminoglycosides have about 270 bacterial and about 320 eukaryotic copies (yeast, Candida, human cytosolic), plus some mitochondrial ones.
   - Tetracyclines have about 50 bacterial, 55 eukaryotic and 16 mitochondrial copies.
   - These sets need the cross-family correspondence (bacterial vs eukaryotic SSU, for example via US-align).
3. **Shared PTC/tunnel hotspots.** About 230 usable copies across 8 classes, with fingerprints computed within sites defined on Rfam positions.

## Still needed
- **Rfam-position mapping:** a mapping for the drug-bound chains, because non-E. coli species have their own numbering.
- **Copy deduplication:** Thermus crystals hold two ribosomes per asymmetric unit.
- **Ligand quality filters:** fit to density, for example Q-score for cryo-EM.
- **Missing contacts:** contacts for the 133 copies without PDBe contacts.

## Follow-up (2026-10-24): resistance-methylated ribosomes already in the PDB
**Question.** Before modelling methyl groups, what do the experimental structures of methylated ribosomes show? Script: `methylated.py`. Output: `out/methylated_pairs.json`.

**Finding the entries.** PDBe API `modified_AA_or_NA` for the 1,130 rRNA entries that contain a methylated A or G, filtered to the resistance positions. Native modifications at other positions (m⁶₂A1518/1519, m²A2503, m⁷G527 and others) were ignored.

| Enzyme | Modification (CCD) | Entries | With a drug |
|---|---|---|---|
| Erm | m⁶₂A2058 (MA6) | 6XHV, 7RQ9, 8FC4, 8FC5, 8FC6, 8UD7, 8VTX | erythromycin, azithromycin, telithromycin, macrolone MCX-128, iboxamycin, cresomycin |
| Cfr | m²,⁸A2503 (8AH) | 7LVK, 8G29, 8G2A, 7S1K, 8G2B, 8G2C, 8UD8 | radezolid (E. coli cryo-EM), iboxamycin, tylosin, cresomycin |
| NpmA | m¹A1408 (1MA) | 9B4Z, 9B51 | arbekacin, G418 (E. coli cryo-EM) |
| ArmA/RmtB | m⁷G1405 | 7CPJ only | none; a stalled E. coli ribosome with no resistance enzyme, so treated as a modelling artefact |

- **Controls:** nearly every drug complex has an unmethylated twin with the same drug from the same paper or lab. 9B51 (G418) has none.
- **Systems:** Thermus crystals hold two ribosomes, and both copies are reported.
- **Coordinates:** PDBe ModelServer `residueSurroundings`, 10 Å around each drug copy.
- **Secondary sites:** aminoglycoside copies more than 15 Å from A1408 are secondary sites. They are counted but not compared.

**Per drug copy** (copy 1 / copy 2 for Thermus; distances in Å):

| Enzyme | Drug | Ideal methyl → drug in the unmethylated complex | Real methyl → drug in the methylated complex | Shift of the modified base | Drug shift, RMSD / max | Drug B ÷ pocket B, methylated vs unmethylated |
|---|---|---|---|---|---|---|
| Erm | erythromycin | 3.4 / 3.2 | 2.9 / **1.4** | 0.8 / 0.9 | 0.7 / 2.7 | 1.61, 1.37 vs 1.21, 1.12 |
| Erm | azithromycin | 3.5 / 3.0 | 2.9 / 2.9 | 0.6 / 0.3 | 0.8 / 2.4 | 1.75, 1.59 vs 1.34, 1.18 |
| Erm | telithromycin | 3.4 / 3.5 | 2.9 / 2.9 | 0.6 / 0.6 | 0.6 / 2.8 | 1.38, 1.24 vs 0.98, 0.96 |
| Erm | macrolone MCX-128 | 3.7 / 3.3 | 3.0 / 3.0 | 0.7 / 0.6 | 0.5 / 1.8 | 1.66, 1.37 vs 1.13, 1.10 |
| Erm | iboxamycin | **1.6 / 1.5** | 2.9 / 2.9 | **1.7 / 2.1** | 0.7, 0.5 / 1.2 | 1.25, 1.04 vs 1.04, 0.90 |
| Erm | cresomycin | **1.5 / 1.6** | 2.8 / 2.8 | **1.6 / 2.3** | 0.3 / 0.6 | 0.98, 0.98 vs 0.97, 1.00 |
| Cfr | iboxamycin | **1.5 / 1.7** | 2.9 / 2.9 | 1.0 / 1.0 | 0.4 / 0.8 | 1.04, 1.00 vs 1.04, 0.90 |
| Cfr | tylosin | **2.0 / 2.2** | 3.2 / 3.2 | 0.3 / 0.5 | 0.8 / 3.0 | 1.19, 1.13 vs 0.93, 0.96 |
| Cfr | cresomycin | **1.5 / 1.7** | 2.9 / 2.9 | 0.4 / 0.9 (U2585 1.4–1.9) | 0.4–0.5 / 1.1 | 1.26, 1.11 vs 0.97, 1.00 |
| Cfr | radezolid (EM) | **2.3** | 3.1 | 0.7 (A2062 3.9, U2585 2.0) | 0.2 / 0.4 | 0.98 vs 1.00 |
| NpmA | arbekacin (EM) | **1.3** (N3 of ring I) | 3.4 | 1.1 | 0.7 / 1.9 (N3) | 1.09 vs 1.00 |
| NpmA | G418 (EM) | no control | 3.1 | – | – | – |

**What this shows**
1. **Placing a methyl is easy and accurate.** Ideal in-plane methyls built on the parent atom lie 0.01–0.3 Å from the deposited methyl carbons in every structure.
2. **A clash test on the unmethylated complex gets resistance backwards.**
   - Macrolides lie 3.0–3.7 Å from where Erm's methyls go, so no steric clash, yet Erm confers strong resistance to them.
   - Iboxamycin and cresomycin, which were designed to escape Erm and Cfr, clash hardest (about 1.5 Å).
   - The ribosome makes room by moving the methylated base: A2058 by 1.6–2.3 Å for the lincosamide-type drugs, and A2503 by 0.3–1.0 Å.
   - So the methyl-clash model proposed above would mispredict. It is dropped.
3. **Every drug still binds the methylated ribosome in these structures.** Crystal soaks and cryo-EM samples use high drug concentrations, so binding at saturation says little about affinity.
   - **Only indirect sign of weaker binding:** the drug's B-factor relative to its pocket rises for the drug classes the enzyme defeats. It rises 0.25–0.5 for the macrolides on Erm, and 0.1–0.3 for tylosin and cresomycin on Cfr and for iboxamycin on Erm.
   - **Where it stays flat:** cresomycin on Erm.
   - **Cryo-EM B-factors** aren't comparable in the same way (9B51 ratio 0.34).
4. **Model issue:** in 8FC4 copy 2, erythromycin O8 is 1.38 Å from methyl C9, with full occupancy and no alternate locations. That is an unresolved steric overlap in the deposited model; copy 1 is at 2.9 Å.
5. **Coverage is narrow.** There are 12 drug complexes from 6 papers, and most are next-generation drugs designed to escape resistance. No structure exists for:
   - the classic Cfr targets: chloramphenicol, linezolid, pleuromutilins, streptogramin A;
   - lincosamides on Cfr;
   - ArmA/RmtB with an aminoglycoside.

**Reading.** The PDB already holds the methylated structures, and comparing them with their unmethylated twins is more informative than modelling methyls. These structures show the ribosome adapting: the methylated base moves and the drug shifts at its edges. They don't show the drug being excluded. Explaining the resistance spectra therefore needs affinity data, not geometry alone.

## Bridge pilot (2026-10-24): one SSU axis across kingdoms
**Question.** The novelty review named one blocking risk: Rfam keeps bacterial and eukaryotic SSU rRNA in separate families (RF00177, RF01960), so there is no shared position axis until a bridge between them is built and shown to work. Script: `bridge_pilot.py`. Output: `out/bridge_pilot.json`.

**Set.** One SSU rRNA chain from each of 8 aminoglycoside complexes, chosen for resolution and species spread:

| Chain | Species | Group | Drug |
|---|---|---|---|
| 7K00 A | E. coli (1.98 Å) | bacterial | paromomycin |
| 2VQE A | T. thermophilus (2.5 Å) | bacterial | paromomycin |
| 9VXM AA | human 12S (2.9 Å) | mitochondrial | neomycin |
| 9VXN S2 | human 18S (2.51 Å) | cytosolic | neomycin |
| 5NDV 2 | S. cerevisiae (3.3 Å) | cytosolic | paromomycin |
| 9G6J A | C. albicans (2.15 Å) | cytosolic | G418 |
| 9FXO S1 | L. major (2.25 Å) | cytosolic | paromomycin |
| 9SUL AA | T. brucei (2.1 Å) | cytosolic | G418 |

**Two independent routes to E. coli numbering.** Both end at an RF00177 consensus position, which 7K00 A turns into an E. coli number.
- **Covariance model:** each chain aligned with `cmalign` to its own family, and every cytosolic chain also to RF00177 (cross-family).
- **Structure:** US-align (`-mol RNA`, C3′) onto 7K00 A; a residue takes the position of the 7K00 residue it pairs with, if within 4 Å after superposition.

| Chain | Own-family fit (bits / avg pp) | Cross-family fit to RF00177 | US-align TM | Routes agree (same / jointly mapped) |
|---|---|---|---|---|
| 7K00 A | 1571.5 / 0.99 | – | 1.00 | 1489 / 1489 |
| 2VQE A | 1566.4 / 0.985 | – | 0.96 | 1368 / 1390 |
| 9VXM AA | 260.3 / 0.887 | – | 0.95 | 614 / 640 |
| 9VXN S2 | 1661.4 / 0.969 | 429.0 / 0.911 | 0.82 | 928 / 985 |
| 5NDV 2 | 1806.6 / 0.968 | 533.6 / 0.895 | 0.78 | 838 / 877 |
| 9G6J A | 1815.9 / 0.978 | 535.1 / 0.891 | 0.82 | 1076 / 1108 |
| 9FXO S1 | 1223.5 / 0.907 | 440.5 / 0.895 | 0.75 | 941 / 992 |
| 9SUL AA | 1199.6 / 0.925 | 420.5 / 0.900 | 0.74 | 932 / 992 |

**The bridge passes its test.**
1. **Known determinants come out right in all 8 chains, with no hand-fixing.** Position 1408 is A in the bacteria and in the human mitochondrial 12S (its own A1492), and G in every cytosolic chain: the aminoglycoside selectivity switch. Positions 1490 and 1491 map to the mitochondrial deafness nucleotides A1555 and C1556. Positions 1492–1496 recover A1557–C1561 (mito), A1824–C1828 (human 18S) and A1755–C1759 (yeast).
2. **The two routes agree on 94–98% of jointly mapped residues**, so each checks the other.
3. **Each chain's single decoding-site aminoglycoside copy is recovered by both routes**, agreeing on 7–9 of the contacted positions. The routes differ mostly at 1406 and 1407.

**Limits found**
- **The cross-family fit is weak** for cytosolic chains: 420–535 bits against about 1,570 for the bacteria. It is usable but should be the check, with the structural route primary.
- **Mitochondrial 12S fits RF00177 at only 260 bits**, and US-align pairs 807 of 949 residues. Mitochondria need per-structure checks, not pooling.
- **Only 55–65% of eukaryotic residues map** to bacterial positions, as expected from the expansion segments. Drug sites lie in the conserved core, so this does not block the analysis.
- **Position 1407** (5MC in E. coli) gets no structural partner in any chain, including 2VQE. Worth a look, since it borders the aminoglycoside site.
- **12 of 7K00 A's contacted residues are its own unmapped positions** (1808, 1862, 2447 and others). These are second-copy or non-consensus residues and need the usual missing-vs-absent handling.

**First cross-kingdom result from 8 chains.** Every eukaryotic chain carries one decoding-site copy plus 2–5 **secondary** aminoglycoside copies (the full 5NDV entry has about 45). Of 17 secondary sites with mapped positions, two groups recur across species:
- **E. coli 1066, 1092, 1109, 1186:** human 18S (NMY) and T. brucei (GET).
- **E. coli 243, 281:** S. cerevisiae (PAR) and C. albicans (GET).
- A third pair is suggestive but below the clustering threshold: L. major and C. albicans copies share 1487 and 1488, next to the decoding site.
These secondary sites are invisible to the sequence-based studies and are exactly what a shared axis makes countable.

**Next.** Extend the bridge to the large subunit (RF02541 and RF02543) for the PTC and exit-tunnel drugs, where the eukaryotic chains are much larger; add ligand-quality flags and copy deduplication before any statistics.
