# Literature reference contacts: curation notes

Reference set for the RNA–ligand interaction-tool benchmark (PDBe-Arpeggio, PLIP, fingeRNAt). All 12 pilot ligand copies (10 entries) are curated. The contacts come from the primary structure papers: main text, figure legends, figure drawings and supplementary files. No contact was derived from the coordinates. I used the mmCIF files only to check that each residue number and residue identity named in a paper exists in the entry (102 named residue positions checked; all match), and to find the ligand copy (author chain and residue number).

- Output: `reference_contacts.json` has 107 records: 65 hbond, 20 stacking, 13 electrostatic, 4 metal_mediated, 1 water_mediated, 4 other.
- Evidence: 78 records are `text`, 28 are `figure` and 1 is `secondary`.
- Retrieval: 2026-10-23. Citations and DOIs come from the PDBe publications API (`/pdbe/api/pdb/entry/publications/{id}`), cross-checked against the RCSB entry primary citation.
- Access: full text came from Europe PMC/PMC (author manuscripts or OA) or from nature.com through institutional access. For the Nature/NSMB papers, figure drawings were read from the article images, and for 1FJG from the article PDF.

## Fields and conventions

| field | meaning |
|---|---|
| `pdb`, `ligand_chain`, `ligand_resnum`, `ligand_ccd` | the ligand copy (author chain and number), the same copies as `prepare.py` `PILOT` |
| `rna_chain`, `rna_resnum`, `rna_resname` | author numbering **of the PDB entry**. When the paper numbers differently, its label is in `paper_numbering` (3E5C only) |
| `rna_atom_or_group` | A bare PDB atom name (`N3`, `O2'`, `O6`…) is used only when the atom is identified. Otherwise the field holds a description, which `score.py` will not parse as an atom. Phosphate oxygens are written "phosphate O (OP1/OP2 not specified)" because no paper says which one |
| `ligand_atom_or_group` | CCD atom names where unambiguous (SAM/SAH: `N` = α-amino, `O`/`OXT` = carboxylate, `SD` = S, `CE` = S-methyl; adenine `N1`, `N6`, `N7`…). The aminoglycosides use conventional ring nomenclature |
| `type` | hbond, stacking, electrostatic, water_mediated, metal_mediated, other. `other` is used for van der Waals contacts and drawn C–H···O lines, and is not scored by `score.py` |
| `evidence` | `text` = stated in words (main text **or figure legend**). `figure` = the contact, or the atom pair, is only readable from a figure drawing. `secondary` = from another paper. `supplementary` was not needed: the supplementary files held no extra contacts |
| `atom_basis` | `stated` = atoms named in words. `figure` = atoms read from the drawing. `implied_by_pairing` = the authors state a pairing edge or type (e.g. Watson–Crick A–U) but not the atoms. `group_only` = no atom identified |
| `type_basis` | `stated`, or `curator_assigned` when the paper gives no type for an individual contact (1FJG streptomycin: "salt bridges and hydrogen bonds" collectively; 1FJG paromomycin ring IV dashed lines) |
| `unit` | Set when the authors describe stacking against a base triple or base pair as a whole, not against one base. Score these at unit level: any member counts |
| `qualifier`, `note` | Hedges from the paper ("likely", "protomers A and C only"), conflicts and ambiguities |
| `distance_A` | Only distances stated in the text or printed in a figure. Figure distances that I paired with a specific line from the drawing are flagged in `note`/`source` |

How `score.py` handles these records:
- Records with `rna_resnum: null` cannot be matched. There is one: the 1Y26 Mg5 record, where no RNA partner is named.
- Records with `unit` are listed per residue, so residue-level recall counts each member separately. Read the 1Y26 triple stacking (6 records) and the 6LAS pair stacking (G9, G7) with that in mind.
- The 4FE5 U22 contact appears twice: once at residue level (primary) and once with the atom pair (secondary). The residue set collapses the duplicate.

## Completeness (for the 'unsupported' counts in `score.py`)

A paper is treated as complete for a class only when the authors explicitly say that they enumerate all of the ligand's contacts of that kind. Otherwise it is partial, and a tool contact absent from the reference is not evidence of a false positive.

| copy | polar (hbond + electrostatic) | stacking | water | metal |
|---|---|---|---|---|
| 1Y26 ADE | **complete**: "all heteroatoms along the periphery … recognized" (5 heteroatoms, 4 residues) | complete at triple level (above and below) | partial: only Mg5–N6 is mentioned | partial |
| 4FE5 HPA | complete at residue level (U22, U47, U51, C74); partial at atom level | partial ("stacks on P1"; residues not named) | not described | not described |
| 2GIS SAM | partial (figure-level; "almost every functional group" recognized) | partial (C47 only) | not described | not described |
| 2QWY SAM | partial | not stated (adenine sits in a gap between U21 and G22; stacking not stated) | not described | not described |
| 3E5C SAM | **likely complete**: the adenine Watson–Crick and Hoogsteen faces, the ribose and the sulfonium are each assigned; the methionine tail is disordered and not recognized | complete (U72, G90) | partial (one water) | none |
| 6FZ0 SAM | **complete**: "hydrogen bonded … at both ends, with no direct contacts with the remaining methionyl chain" | complete (G21, A48) | none described | complete (one Mg) |
| 6LAS SAM | likely complete for the adenine, ribose and amino group | complete (C32/G33 explicit; G9-C32 and G7·G33 pairs) | none | none |
| 3NPQ SAH | **complete**: "All of the hydrogen bonding interactions … are mediated through the adenine ring and the main chain atoms" | complete (A29, C16) | none | none |
| 2GDI TPP | partial (HMP atoms not named; Mg2-bridged waters not curated) | complete (G42, A43) | partial (Mg-bound waters not curated) | partial (Mg1 only) |
| 1FJG PAR | partial (figure-level; ring III "weak contacts") | complete (G1491) | none | none |
| 1FJG SRY | partial (figure-level); also contacts protein S12 K45, which is not in the reference | none | none | none |
| 1FJG SCM | partial (figure-level) | none | none | none |

## Per complex

### 1. 1Y26 chain X, adenine (ADE X/90): add A-riboswitch
- **Citation**: Serganov et al. 2004 *Chem Biol* 11:1729, doi:10.1016/j.chembiol.2004.11.018.
- **Access**: PMC4692365 (author manuscript; the Europe PMC fullTextXML endpoint returned nothing, so I used the PMC HTML). Figure 7 image checked.
- **Curated**: all text-stated H-bonds.
  - Adenine N1 and N6 pair with the Watson–Crick edge of U74; the U74 atoms N3/O4 are implied by the stated Watson–Crick A–U pair.
  - N3 and N9-H pair with the Watson–Crick edge of U51 (U51 atoms not named).
  - N9-H to U47 O2.
  - N7 to the U22 2′-OH.
  - Stacking over the A73•(A52-U22) and C50•(U75-A21) triples.
  - Hydrated Mg5 within H-bond distance of N6 (no RNA partner named).
- **Numbering**: the paper uses xpt (*B. subtilis*) numbering for the add aptamer. The PDB uses the same numbering (13–83). Verified.
- **Ambiguities**: Fig. 7A draws the two U51 lines, but the U51 atoms are not labelled, so they stay group-level. Mg5 is not mapped to an MG residue number.
- **Not curated**: water molecules in the A73•(A52-U22) triple, which is RNA–RNA.

### 2. 4FE5 chain B, hypoxanthine (HPA B/101): xpt guanine riboswitch
- **Citation**: the PDB primary citation of 4FE5 is **Batey, Gilbert & Montange 2004 *Nature* 432:411**, doi:10.1038/nature03037.
  - 4FE5 is the re-refined successor of 1U8D (1U8D is no longer served). The coordinates were deposited in 2012 with Stoddard et al. 2013 *J Mol Biol* (doi:10.1016/j.jmb.2013.02.023, PMC3769107), which lists 4FE5 as "variant 1A:HX" and thanks Richardson/Murray for correcting backbone errors in the original 1U8D model.
  - Stoddard 2013 states that the ligand H-bond pattern is unaltered among its variants, but it does not enumerate the contacts. So the contacts come from Batey 2004.
- **Access**: nature.com full text (institutional) and supplementary .doc. Fig. 2b viewed: dashed lines are drawn but atoms are not labelled.
- **Curated**: H-bonds with U22, U47, U51 and C74, at residue level (text). C74 is a Watson–Crick partner, with the atoms implied.
- **Secondary**: U22 O2′–N7 is from Batey 2012 *Q Rev Biophys* 45:345 (same lab, PMC3760793, §2.1, Fig. 2c), which describes the purine-riboswitch ligand as recognized by the "N7 and the 2′-hydroxyl group of U22".
- **Numbering**: full-length mRNA numbering (paper Fig. 1a legend). It matches the PDB.
- **Explicit non-contacts**: hypoxanthine lacks the 2-NH2. The 2-NH2 of guanine would H-bond the O2 carbonyls of C74 and U51 (Batey 2004, text), so **HX makes no contact at its C2 position with C74 O2 or U51 O2**. Serganov 2004 says the same: two H-bonds are absent in the HX complex. A tool contact from HX C2-H to these oxygens would be a weak contact at most.
- **Not curated**:
  - Atom-level U47/U51 pairs. Gilbert et al. 2009 *Structure* (PMC2721690) discusses N3(U51) and O2(U51) for **xanthine and 2-fluoroadenine**, not HX, so it was not transferred.
  - Stacking: the text says the quadruple "stacks directly on the P1 helix" without naming residues.
  - Cobalt hexammine sites.

### 3. 2GIS chain A, SAM (SAM A/301): SAM-I (*T. tengcongensis* metF)
- **Citation**: Montange & Batey 2006 *Nature* 441:1172, doi:10.1038/nature04819.
- **Access**: nature.com (institutional) and supplementary .doc (which has no contact detail). Fig. 4b–d read at enlarged scale.
- **Curated**:
  - Adenine triple: U57 O4–N6 and U57 N3–N7 (2.9/2.9 Å); A45 N1–N6 (2.8 Å) and A45 2′-OH–N1 (3.0 Å). The atom pairs are from Fig. 4b; the text states the triple.
  - Stacking with C47.
  - SAM 2′-OH to "O4 of C47" (see below).
  - Carboxylate to G11 N1/N2 (text; 3.1/3.1 Å, Fig. 4c).
  - Carboxylate–G58 N2 (2.8 Å) and α-amino–G58 ribose OH (3.2 Å) (figure).
  - Sulfonium ~4 Å from the U7 and U88 O2 (electrostatic, text).
  - van der Waals packing with U7/U88 (`other`).
- **Ambiguities**:
  - "O4 of C47": cytidine has no O4; O4′ is probably meant, so the atom is left unparsed.
  - The G58 ribose oxygen is drawn but not labelled.
  - Atom naming conflict: Sun et al. 2019 (SAM-VI paper, citing other SAM-I work) call the SAM-I sulfonium partners O4 of U7/U88. The primary paper says O2 (minor groove), and the primary was followed.
- **Explicit non-contacts**:
  - The S-methyl is not directly recognized; it points into a solvent cavity.
  - The RNA "has evolved to not use the anionic phosphate backbone" for either positive group, so phosphate–sulfonium or phosphate–amine contacts are not reported.
  - The ribose makes only the one H-bond (to C47); P1 contacts are van der Waals only.
- **Not curated**: the intramolecular π–cation interaction (amino group over the adenine), which is ligand-internal.

### 4. 2QWY chain A, SAM (SAM A/100): SAM-II (Sargasso Sea metX)
- **Citation**: Gilbert, Rambo, Van Tyne & Batey 2008 *Nat Struct Mol Biol* 15:177, doi:10.1038/nsmb.1371.
- **Access**: nature.com (institutional) and supplementary PDF. Fig. 3b–d enlarged.
- **Curated**:
  - The adenine Hoogsteen face pairs with U44: U44 N3–N7 and U44 O2–N6 (figure atoms).
  - U10 O4 drawn to adenine C8, typed `other`. The text only places the adenine in a U10/U44 triple.
  - Sulfonium and S-methyl recognized by the U11 and U21 carbonyls (text; O4 from the figure; electrostatic).
  - A47 N1–α-amino and A47 N6–carboxylate (text and figure).
- **Ambiguities and conflicts**:
  - The text says A47 "N1 and N6 [bond] to the carbonyl and amino groups … respectively". That pairing is chemically reversed and contradicts Fig. 3d, which draws N1···amino and N6···carboxylate. The figure was used.
  - Methionine main-chain contacts differ among the three protomers. Molecules A and C bind along the A47 Watson–Crick face; in **molecule B the amino group does not H-bond**. Chain A is assumed to be "molecule A", and the paper's protomer labels are not mapped to chain IDs explicitly.
  - A45 is described as part of the "hydrogen-bonding network that recognizes the charged sulfur". In Fig. 3c, A45 N6 bonds U11 O4 and U21 O4, not SAM, so A45 is **not** recorded as a direct contact.
- **Explicit non-contacts and hedges**:
  - The adenine Watson–Crick face is solvent-exposed in the minor groove.
  - The methionine methylenes are poorly ordered and are "probably" held by weak van der Waals contacts with the U11/U21 carbonyls. This is speculative and not recorded.
  - Stacking of the adenine in the P2b gap between U21 and G22 is implied by the description of the "opening" but never stated as stacking, so it is not recorded.

### 5. 3E5C chain A, SAM (SAM A/216): SAM-III / S_MK box (*E. faecalis* metK)
- **Citation**: Lu et al. 2008 *Nat Struct Mol Biol* 15:1076, doi:10.1038/nsmb.1494.
- **Access**: PMC3467307 (author manuscript HTML). Figs. 1 and 3 viewed.
- **Numbering (important)**: the paper uses metK transcript numbering, while 3E5C chain A is numbered 1–53. The mapping was derived from the construct drawing (Fig. 1b: U21 label, G31 label, U65, A74, SD^89) and checked by base identity at every named position and by P1/P2/P4 pairing:
  - PDB 1–12 = paper 20–31 (PDB 1 is GTP = G20).
  - PDB 13–28 = engineered P3 and loop, with no native numbers.
  - PDB 29–41 = paper 64–76.
  - PDB 42–45 = the GAAA that replaces the hypervariable linker 77–87.
  - PDB 46–53 = paper 88–95.
  - Contacts are recorded with PDB numbers, and the paper label is in `paper_numbering`. G26 → 7, G71 → 36, U72 → 37, A73 → 38, A74 → 39, G89 → 47, G90 → 48 (also C25 → 6, U65 → 30, A64 → 29).
- **Curated**: all text-stated contacts.
  - G26 minor groove: SAM N1···G26 N2, SAM N6···G26 N3, SAM N6···G26 2′-OH (2.9/2.9/2.8 Å; the two 2.9 Å labels were assigned from the drawing).
  - A73 N6–N7 (2.9 Å).
  - G89 N7···SAM 2′-OH (π-H-bond, 3.1 Å).
  - G89···SAM 3′-OH (2.5 Å; conflict, see below).
  - Water-mediated SAM 2′-OH···A74 N6.
  - Sulfonium electrostatics with U72 O4 (4.1 Å) and the G71 2′-OH (3.9 Å).
  - Stacking between U72 and G90.
- **Conflict**: the main text says the SAM 3′-OH H-bonds the **2′-hydroxyl** of G89. The Fig. 3e legend and drawing say a **phosphoryl oxygen** of G89. The residue is certain and the atom is left unparsed.
- **Explicit non-contacts**:
  - No density beyond the sulfonium: the methionine tail is not recognized, and the modelled tail conformation is not evidence.
  - The methyl group is not specified.
  - The sulfonium–O4′ (3.0 Å) and sulfonium–N3 (3.8 Å) contacts are **intramolecular** and not recorded.
- **Not curated**: the 3E5E (SAH) and 3E5F (Se-SAM) complexes.

### 6. 6FZ0 chain A, SAM (SAM A/104): SAM-V (*Ca.* Pelagibacter ubique metY)
- **Citation**: Huang & Lilley 2018 *Nucleic Acids Res* 46:6869, doi:10.1093/nar/gky520. RCSB lists the authors as Huang L, Lilley DMJ.
- **Access**: Europe PMC OA full-text XML (PMC6061858). Fig. 6 checked.
- **Curated**: all text-stated contacts.
  - Trans Hoogsteen pair with U47: N6–U47 O2 and U47 N3–N7.
  - Stacking between G21 and A48.
  - A50 N6 → carboxylate O and A50 N1 ← amine.
  - Carboxylate–Mg²⁺ (inner-sphere waters)–U10 pro-R phosphate O, metal-mediated. The metal–O distance of 2.1 Å is why the authors call it Mg.
  - Sulfonium to U20 O4 (3.2 Å) and U9 O4 (3.3 Å), electrostatic. The authors say "there is no direct contact" but an electrostatic interaction.
- **Numbering**: same as the PDB. Nucleotides 33–36 (P3 loop) are disordered. Position 4 is 5-bromocytidine (phasing).
- **Explicit non-contacts**:
  - No direct contacts between the two ends of the extended SAM chain.
  - The A52 interaction with U17 is probably a construct artefact, and it is RNA–RNA.

### 7. 6LAS chain A, SAM (SAM A/101): SAM-VI (*B. angulatum* metK)
- **Citation**: Sun, Gasser, Li, …, Micura, Ren 2019 *Nat Commun* 10:5728, doi:10.1038/s41467-019-13600-9.
- **Access**: Europe PMC OA XML (PMC6914780) and the Fig. 2 PNG (panels f, g read).
- **Curated**:
  - The Hoogsteen edge of the adenine pairs with the U8 Watson–Crick edge: U8 N3–N7 and U8 O2–N6 (figure atoms).
  - U8 O4–α-amino (text).
  - A36 N6–SAM N1 and A37 N1–SAM N6 (text; the ligand atom for A37 follows from donor/acceptor logic and the drawing).
  - G33 non-bridging phosphate O–SAM 2′-OH.
  - Stacking: C32 and G33 explicitly; the G9-C32 and G7·G33 pairs at pair level.
  - Sulfonium–U6 O4, 4.8 Å, "likely assists" (electrostatic, hedged).
- **Ambiguities**:
  - A "3.3" label in Fig. 2f sits between the SAM N6···U8 O2 and the G38 N2···U8 O2 lines, so it was not assigned.
  - U8 O4 is 5.1 Å from the sulfonium (Fig. 2g, "nearest neighbours") but is not described as an interaction. It is not recorded.
- **Explicit non-contacts**:
  - G38 "does not form direct interactions with SAM" (it contacts U8 O2).
  - G7 "does not hydrogen-bond with SAM" (it stacks).
- **Other copies**: chain B also holds SAM 101; it is not curated separately because the paper describes one complex. U1A protein chains are present.

### 8. 3NPQ chain A, SAH (SAH A/55): SAH riboswitch (*R. solanacearum* ahcY)
- **Citation**: Edwards, Reyes, Héroux & Batey 2010 *RNA* 16:2144, doi:10.1261/rna.2341610.
- **Access**: PMC2957054 (HTML). Fig. 3 viewed; atoms are not labelled.
- **Curated**: all text-stated contacts.
  - Sheared G15·A(SAH) pair: two H-bonds drawn, atoms not named.
  - Adenine intercalated between A29 and C16.
  - α-amine near G30 N3 and a G30 non-bridging phosphate O ("charge–charge"), electrostatic.
  - Carboxylate to the G47 and G31 2′-OH groups.
  - S in van der Waals contact with A29 O4′ (3.2 Å, `other`).
- **Hedges**: the main-chain contacts differ slightly among the three protomers, and "likely … only a majority at any given time". The chain A choice follows `prepare.py`; the paper does not say which protomer Fig. 3A shows.
- **Numbering**: crystal-construct numbering equals the PDB; nucleotides 40–42 are not modelled (domain swap). The wild-type probing numbering differs (crystal 32–38 = WT 46–52, Fig. 2B), but that offset applies only to the probing panel.
- **Explicit non-contacts**: the 2′-OH, 3′-OH and adenine sugar edge of SAH project away from the RNA. SAM/SAH discrimination is steric (a modelled methyl clashes with the A29 C4′ and G31 N7), and that model is not a contact.

### 9. 2GDI chain X, TPP (TPP X/100): *E. coli* thiM TPP riboswitch
- **Citation**: Serganov, Polonskaia, Phan, Breaker & Patel 2006 *Nature* 441:1167, doi:10.1038/nature04740.
- **Access**: PMC4689313 (author manuscript HTML) and supplementary PDF from nature.com (Suppl. Figs. S1–S11 text and S1 image read).
- **Chains**: the entry has two complexes, chains X and Y, and the task's "chain A" does not exist. The paper's figures show "molecule 1", assumed to be X. Suppl. Fig. S9 says both molecules contain "virtually all described interactions", so the same contacts should hold for Y/100.
- **Curated**:
  - The HMP ring intercalates between G42 and A43 (stacking) and H-bonds the G40 base and the G19 2′-OH. The HMP atoms are not named.
  - Mg1 directly coordinates the O6 of G60 and G78 and the terminal phosphate (metal-mediated).
  - C77 and G78 H-bond the terminal phosphate oxygens (atoms not named).
- **Explicit non-contacts**: the central thiazole is "not recognized by the RNA". Suppl. Fig. S8: G72 N3 is the residue closest to the ligand's central ring (a modelled PTPP C–H···N), so a tool contact between G72 and the thiazole would be at most a weak contact that the authors do not claim for TPP.
- **Not curated**:
  - Mg2: it coordinates both phosphates and is positioned by direct and water-mediated contacts to RNA, but no residue is named in the text, and Fig. 2d's partners cannot be resolved.
  - Individual Mg-bound waters bridging to the RNA.
  - HMP atom identities: none are given in the paper or the supplement.

### 10. 1FJG 16S chain A: paromomycin (PAR A/1545), streptomycin (SRY A/1634), spectinomycin (SCM A/1633)
- **Citation**: Carter et al. 2000 *Nature* 407:340, doi:10.1038/35030019.
- **Access**: nature.com (institutional) and the article PDF. Panel b of Figs. 4–6 (chemical interaction diagrams) was read from the PDF page images at about 1500 px.
- **Numbering**: the paper states that it uses *E. coli* 16S numbering (Wimberly et al. 2000). The 1FJG chain A author numbering follows it, and every named residue matches by identity.
- **Paromomycin**:
  - Text: ring I stacks on G1491 and H-bonds A1408; ring I H-bonds the A1493 phosphate; ring IV contacts the backbone on both sides of H44; ring III makes only weak contacts.
  - Fig. 6b: ring I 6′-OH···A1408 N1; ring I 4′-OH···A1493 phosphate; ring III 5″-OH···G1491 N7; ring II N3···G1494 N7; ring II N1···U1495 carbonyl (O4/O2 not resolvable); ring IV 6‴-NH2···C1490 phosphate; ring IV 3‴-OH···G1405 phosphate.
  - The ring-position labels are my reading of the drawn substituent pattern and are noted as such.
  - No A1408–A1493 pair, unlike the NMR structure.
- **Streptomycin**:
  - Text: bound to the phosphate backbone of four regions through salt bridges and H-bonds, and also to S12 K45. K45 H-bonds two streptomycin OH groups; protein contacts are excluded from the reference and from the crops in `prepare.py`.
  - Fig. 5b: streptidine guanidine···C1490 phosphate; streptidine OH···G1491 phosphate; the other guanidine···two A914 ribose oxygens; N-methylglucosamine NHMe···A914 phosphate; its OH···U14 phosphate; streptose 3′-OH···G527 phosphate; formyl O···C526 phosphate and ···an A913 ribose O.
  - The text's list "13, 526, 915, 1490" is the biochemical (protection/crosslinking/mutation) list. The structure figure labels U14, A913 and A914, and the text adds that A915 makes no contacts with other bases. Per-contact types are curator-assigned (`type_basis`).
- **Spectinomycin**:
  - Text: minor groove at the end of H34, a single contact with a 2′-OH, most interactions with G1064 and C1192.
  - Fig. 4b: an actinamine OH···G1064 O6 and N7; NHMe···C1066 O2; the other NHMe···G1068 2′-OH; a bridgehead hemiketal OH···C1192 N3.
  - The text says a loop of S5 and part of H28 are within 5 Å but make **no direct contacts**, which is an explicit non-contact.
- **Caveats**: the structure is at 3.0 Å with three antibiotics soaked together, and the figure line endpoints were read from 2000-era drawings. Treat all 1FJG atom pairs as figure-level with moderate confidence. Ligand atoms are described chemically, not by CCD atom names.

## Summary of access routes

| entry | primary paper access | figures | supplementary |
|---|---|---|---|
| 1Y26 | PMC4692365 (manuscript) | PMC images | not needed |
| 4FE5 | nature.com (institutional) | Fig. 2 (atoms unlabelled) | .doc read (no contacts) |
| 2GIS | nature.com (institutional) | Fig. 4 | .doc read (no contacts) |
| 2QWY | nature.com (institutional) | Fig. 3 | PDF downloaded (not needed for contacts) |
| 3E5C | PMC3467307 (manuscript) | Figs. 1, 3 | not needed |
| 6FZ0 | Europe PMC OA XML | Fig. 6 | not needed |
| 6LAS | Europe PMC OA XML | Fig. 2 | not needed |
| 3NPQ | PMC2957054 | Fig. 3 | not needed |
| 2GDI | PMC4689313 (manuscript) | Fig. 2 | PDF read (no atom detail) |
| 1FJG | nature.com and article PDF | Figs. 4–6 | not used |

There is one secondary record, 4FE5 U22 O2′–N7, from Batey 2012 *Q Rev Biophys* (PMC3760793). No complex was left uncurated.
