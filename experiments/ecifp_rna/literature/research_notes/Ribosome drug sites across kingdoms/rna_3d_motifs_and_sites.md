# RNA 3D motifs and small-molecule binding sites (ribosome focus) and cross-kingdom rRNA comparison

Scope note (2026-10-24): about 15 tool calls. Items marked **[V]** were checked against a fetched page or search result in this session. Items marked **[K]** come from prior knowledge of the literature; the DOI is given but was **not re-verified** in this session and should be checked before citing. Where nothing was found, the item is listed under Gaps.

## Q1. Which studies link small-molecule binding to RNA 3D motifs or local structural modules, and do any include ribosomes?

### Takeaway
There is a substantial literature linking ligands to RNA *secondary-structure* motifs (Disney's Inforna, R-BIND 2.0, RNALigands, fpocketR's motif model) and to *base-pair graph* neighbourhoods (RNAmigos). None of the resources found here annotates binding sites with **BGSU RNA 3D Motif Atlas loop classes** or recurrent 3D modules (A-minor, kink-turn, T-loop and so on). Ribosomes are either excluded (R-BIND explicitly), outside the motif framework (Inforna is built from small-motif screens), or present only as large cryo-EM entries that are not motif-annotated (HARIBOSS).

### Cited Findings
- **Inforna** (Disney lab) matches RNA secondary-structure motifs (internal loops, hairpins, bulges) to small molecules, using data from two-dimensional combinatorial screening (2DCS) — Velagapudi, Gallo & Disney, *Nat Chem Biol* 2014, 10:291–297, https://doi.org/10.1038/nchembio.1452 **[K]**
- Disney's 2DCS probed more than 3,000,000 RNA motif–small molecule combinations to find "privileged" RNA motifs, and found hairpin loops to be the preferred motif space for small-molecule binding; typical Kd 1–10 µM — Tran & Disney, *Nat Commun* 2012, https://doi.org/10.1038/ncomms2119 (summary seen in search, [science.gov topic page](https://www.science.gov/topicpages/r/rna+structural+motifs.html)) **[V, secondary]**
- **R-BIND** was "the first database of non-ribosomal RNA-binding ligands to include bioactivity as a criterion for inclusion" (2017), so ribosomal ligands are **excluded by design**. R-BIND 2.0 added curated RNA target structural motifs (secondary and tertiary) and a structure search that returns small molecules that target similar motifs — [R-BIND 2.0 bioRxiv 10.1101/2022.03.14.484334](https://www.biorxiv.org/content/10.1101/2022.03.14.484334.full.pdf) **[V]**
- **RNALigands** is a database and web server that pairs RNA secondary-structure motifs with bound small molecules, drawn from PDB, R-BIND and Inforna — [RNA 28(2):115](https://rnajournal.cshlp.org/content/28/2/115.full) **[V]**
- **HARIBOSS** (Panei, Torchet, Menager, Gkeka, Bonomi 2022) is a curated set of PDB RNA–small molecule complexes (X-ray, NMR, cryo-EM) with ligand and pocket annotations. Ribosomes are **included**: "all the structures above 102 kDa are ribosomal RNA/protein complexes", and cryo-EM, "in particular of large ribosomal RNA/protein complexes", made up about 11% of the database (typical resolution 3.0 Å). Pockets are clustered by backbone RMSD (cutoff 1.5 Å) within 90%-identity sequence clusters. The preprint text contains no loop-type or 3D-motif classification of pockets. 267 unique RNA binders — [bioRxiv 10.1101/2022.05.17.492306](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full.pdf) (journal: *Bioinformatics* 2022, https://doi.org/10.1093/bioinformatics/btac483 **[K]**) **[V for content]**
- A 2025 benchmark built from HARIBOSS complexes added a single high-resolution ribosome (PDB 7K00), which suggests ribosomes are usually treated as a special case in ML benchmarks — [PMC11974927 (search summary)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11974927) **[V, snippet only]**
- **fpocketR** (Veenbaas, Koehn, Irving, Lama, Weeks; *PNAS* 2025, 10.1073/pnas.2422346122) adapts fpocket to RNA, because plain fpocket over-predicts pockets in duplex grooves. The paper notes that only about 50 drug-like RNA pockets had been visualized, reports "dozens" of new sites, and "supports a model for RNA secondary structural motifs able to bind quality ligands". Pockets were validated in cells with Frag-MaP — [Nature Index record](https://www.nature.com/nature-index/article/10.1073/pnas.2422346122); [bioRxiv 2025.03.13.643147](https://www.biorxiv.org/content/10.1101/2025.03.13.643147.full.pdf) **[V]**
- Hewitt, Calabrese & Schneekloth (*Bioorg Med Chem* 2019, 27:2253–2260) ran a protein pocket-detection method over all RNA structures in the PDB. RNA pocket buriedness, hydrophobicity and volume overlapped substantially with protein pockets, and many folded RNAs adopt plausibly ligandable conformations. The approach is not motif-based — [NCI Frederick record](https://ncifrederick.cancer.gov/scientificlibrary/publication/Approved.aspx?id=16402) **[V]**
- **RNAmigos** (Oliver, Mallet, Sarrazin Gendron et al., *NAR* 2020, https://doi.org/10.1093/nar/gkaa583 **[K]**) represents binding sites as "augmented base-pairing networks", i.e. graphs with Leontis–Westhof non-canonical edges. This is the closest ML analogue to an FR3D/motif representation. **RNAmigos2** (*Nat Commun* 2025, https://doi.org/10.1038/s41467-025-57852-0) does virtual screening with deep graph learning and large-scale docking data — [PMC11928640](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11928640/), [bioRxiv 2023.11.23.568394](https://www.biorxiv.org/content/10.1101/2023.11.23.568394.full.pdf) **[V for existence and DOI]**
- **RLBind**: Wang et al., *Brief Bioinform* 24(1) bbac486, 2023, https://doi.org/10.1093/bib/bbac486 — [issue page](https://academic.oup.com/bib/issue/24/1) **[V]**. **RNAsite**: Su, Peng & Yang, *Bioinformatics* 37:36–42 (2021), combining structure-based and sequence-profile predictors — [RNAsite server](https://yanglab.qd.sdu.edu.cn/RNAsite/) **[V]**. **RLsite** (2025, *Bioinformatics*, 10.1093/bioinformatics/btaf447) is a different tool that combines an RNA language model with graph attention — [PMC12417085](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12417085/) **[V]**
- The RNA 3D Motif Atlas classifies hairpin, internal and multi-helix junction loops from the BGSU representative sets, automatically and with regular updates. Recurrent modules include GNRA/UNCG, the anticodon and T-loop hairpin loops, sarcin–ricin / loop E, kink-turns, C-loops and the 11-nt receptor — Petrov, Zirbel & Leontis, *RNA* 2013, https://doi.org/10.1261/rna.039438.113 ([PMC3854523](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3854523/)); methods update *Methods* 2016 ([PMC4921307](https://pmc.ncbi.nlm.nih.gov/articles/PMC4921307/)) **[V]**

### Inferences
- The field's "motif" vocabulary for ligands is overwhelmingly **secondary-structural** (loop type and size, from 2DCS or from dot-bracket). RNAmigos uses base-pair graphs but does not map to Motif Atlas classes. Cross-referencing ligand sites against **Motif Atlas loop IDs and FR3D annotations** (base triples, A-minor, stacking) appears not to have been done systematically.
- Ribosomes are the largest single source of RNA–drug structures, yet motif-oriented resources leave them out (R-BIND) or include them without motif annotation (HARIBOSS). Ribosomal drug sites are therefore an under-annotated gap that a motif-on-consensus-position framework could fill.

### Gaps
- I could not confirm whether the RNAmigos, RNAmigos2, RNAsite or RLBind training sets include or exclude ribosomal pockets: the PMC full text was blocked by a browser check. Check the methods sections directly.
- I found no study that tabulates ligand-binding sites by Motif Atlas loop class (HL/IL/J3) across the PDB. The searches returned none, but absence of evidence after about 15 calls is not proof.
- I did not verify whether HARIBOSS's current web version adds secondary-structure or loop annotations to pockets. The 2022 preprint had none.

## Q2. Have ribosome antibiotic sites been described in terms of RNA 3D motifs? Is there a systematic motif-level annotation?

### Takeaway
Individual ribosomal drug sites are routinely described in structural terms: the h44 decoding-site internal loop (A1492/A1493/A1408), the PTC nucleotides (A2062, U2504, U2506, U2585, A2602), the exit-tunnel loops, and the GTPase centre. These descriptions are per site and per paper. I found **no systematic annotation that assigns every ribosome drug site to Motif Atlas loops and recurrent 3D modules**.

### Cited Findings
- A cryo-EM study of 17 compounds from six antibiotic classes found high structural conservation in how drugs that share a scaffold bind, including ordered waters in the pockets. This is a structural, not motif-level, comparison — [PMC10497419 "Structural conservation of antibiotic interaction with ribosomes"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10497419/) **[V]**
- Conserved PTC nucleotides (A2062, U2504, U2506, U2585, A2602) differ in orientation between bacterial and eukaryotic ribosomes. U2504 sits in mutation networks that weaken PTC antibiotic binding — S. aureus 50S structure document (patent, a weak source) [USPTO 10723770](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10723770) **[V, low-quality source]**
- Antibiotic target sites are highly conserved, and small structural differences between bacterial and eukaryotic cytoplasmic ribosomes usually explain selectivity. Inhibition of mitochondrial translation is a major side-effect problem — [Mol Microbiol 2006, 10.1111/j.1365-2958.2006.05063.x](https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2958.2006.05063.x) **[V]**
- Standard reviews of the antibiotic site catalogue: Wilson, *Nat Rev Microbiol* 2014, 12:35–48, https://doi.org/10.1038/nrmicro3155 **[K]**; Lin et al., *Annu Rev Biochem* 2018, https://doi.org/10.1146/annurev-biochem-062917-011942 **[K]**
- The PTC as an ancient core: in the accretion model, the PTC and decoding regions are the oldest common-core rRNA — Petrov et al., *PNAS* 2015, "History of the ribosome and the origin of translation", https://doi.org/10.1073/pnas.1509761112 (summary in search; [Williams lab PDFs](https://williams.chemistry.gatech.edu/publications/LDW_115.pdf)) **[V, summary]**

### Inferences
- The A-site internal loop (h44), the PTC central loop of domain V (a multi-helix junction in Motif Atlas terms), the h18 530 loop and the GTPase-associated hairpin and internal loops (the thiostrepton site) map naturally onto BGSU loop classes. Producing that mapping per kingdom on Rfam consensus positions would appear to be **new as a systematic resource**, even though every individual site is well characterised.
- A good framing for usefulness is the aminoglycoside story: the A-site internal loop is the canonical, independently crystallised "RNA motif that binds a drug" (Westhof lab model oligonucleotides, **[K]**), which shows that motif-level abstraction of a ribosomal site has precedent and value.

### Gaps
- I found no published database or table of "ribosome drug site → Motif Atlas loop ID / FR3D triples / A-minor". Two leads to check: the BGSU group's own ligand-site annotations in RNA 3D Hub, and RiboVision2's data layers.
- The Westhof A-site model papers (Vicens & Westhof; Francois et al.) were not re-verified this session.

## Q3. How have rRNA 3D structure and motifs been compared across kingdoms? Is there a common coordinate system?

### Takeaway
The Williams/Petrov (Georgia Tech) programme is the main cross-kingdom framework: the accretion model, the "common core", RiboVision/RiboVision2 and ProteoVision, built on structure-informed alignments across a phylogenetic sample. BGSU supplies structure-to-alignment mapping (R3D-2-MSA, Rfam column-to-unit files). There is **no single universally adopted rRNA coordinate system**: E. coli numbering is the de facto reference, and RiboVision alignments are a curated alternative.

### Cited Findings
- Petrov et al., *PNAS* 2014, 111:10251–10256, "Evolution of the ribosome at atomic resolution", https://doi.org/10.1073/pnas.1407205111. Bacterial and eukaryotic ribosomes were superimposed; eukaryote-specific rRNA inserts into the common core without perturbing local conformation; 64 putative branch-helix insertions involving 59 ancestral expansion segments — [Williams lab](https://williams.chemistry.gatech.edu/publications/LDW_115.pdf), [PMC4555209 etc. via search](https://pmc.ncbi.nlm.nih.gov/articles/PMC4555209) **[V, summary]**
- Critique: Caetano-Anollés & Caetano-Anollés (*Front Genet* 2015) argue that the Petrov–Williams deconstruction is subjective and point to a different origin — search summary [PMC4451634](https://pmc.ncbi.nlm.nih.gov/articles/PMC4451634) **[V, summary]**
- **RiboVision2** (McCann et al., *J Mol Biol* 2024, 436:168556) visualises rRNA in sequence, 2D and 3D. It uses 152 species from DESIRE in a phylogenetic browser, with preloaded alignments for 16S/18S, 23S/28S, 5.8S and 5S. Any wwPDB ribosome is mapped to the alignments by **MAFFT** (sequence alignment, not structural). The 2013 version handled only about a dozen hand-processed structures — [RiboVision2 data page](https://apollo2.chemistry.gatech.edu/AboutRiboVision2/ribovision_data), [about](https://apollo2.chemistry.gatech.edu/AboutRiboVision2/about), [Benasque 2024 talk](https://www.benasque.org/2024rna/talks_contr/0113_RiboVision2_Petrov.pdf) **[V]**. Original RiboVision: Bernier et al., *Faraday Discuss* 2014, https://doi.org/10.1039/C3FD00126A **[K]**
- **ProteoVision** is the companion server for ribosomal proteins — [ProteoVision navigation](https://apollo2.chemistry.gatech.edu/AboutProteoVision/basic_navigation/) **[V for existence]**
- **SEREB** (Sparse and Efficient Representation of Extant Biology): a small, phylogenetically balanced rRNA/rProtein sequence set used to define the universal core — Bernier, Petrov, Kovacs, Penev & Williams, *Mol Biol Evol* 2018, 35:2065–2076, https://doi.org/10.1093/molbev/msy101 **[K, not found by search this session]**
- Template-based secondary structures for rRNAs of all three domains (a cross-kingdom 2D reference): Petrov et al., *PLoS One* 2014, https://doi.org/10.1371/journal.pone.0088222 **[K]**
- **R3D-2-MSA** (BGSU and Gutell) links PDB nucleotides to columns of large rRNA multiple sequence alignments — Cannone et al., *NAR* 2015 (web server issue), https://doi.org/10.1093/nar/gkv543 **[K]**
- BGSU publishes, per Rfam family and release, PDB chains aligned to the Rfam seed and **column-to-unit-ID** maps (for example `RF00003_15.0_column_to_unit_id.txt`) — [rna.bgsu.edu/data/alignments/rfam](https://rna.bgsu.edu/data/alignments/rfam/) **[V]**
- DARC (database of aligned ribosomal complexes) superposes ribosome structures in a common frame (3D superposition, not a sequence coordinate system) — [CIPSM entry](https://www.cipsm.de/publications/research-area-c/the-darc-site-a-database-of-aligned-ribosomal-complexes/) **[V for existence]**

### Inferences
- The two existing cross-kingdom coordinate systems are E. coli numbering and RiboVision's curated alignments. Rfam consensus (CM match-state) positions are a third option. Their advantages are that they are reproducible, versioned, have published PDB mappings (Rfam/BGSU) and come from the same machinery for any family. Their weakness is that **each kingdom has its own family** (RF00177 bacterial SSU, RF01960 eukaryotic SSU, RF02541 bacterial LSU, RF02543 eukaryotic LSU, plus mitochondrial families), so cross-kingdom equivalence needs a family-to-family position map, for example through the RiboVision alignment or 3D superposition.
- RiboVision2 maps structures to alignments by MAFFT, so motif-level (FR3D/Motif Atlas) annotation of drug sites at cross-kingdom aligned positions does not appear to be a RiboVision feature. This was not confirmed: its data layers were not inspected.

### Gaps
- Whether RiboVision2 shows ligand/antibiotic contacts or Motif Atlas loops as data layers was not checked.
- I found no formal "universal rRNA position numbering" standard comparable to the GtRNAdb/Sprinzl tRNA numbering.

## Q4. Do Rfam rRNA families map PDB rRNA chains to consensus positions across kingdoms? Known issues?

### Takeaway
Yes, within a family: Rfam lists PDB structures per family (RF02541 shows 1730 structures), and BGSU provides column-to-unit maps. Rfam 15 used 3D structures to improve 65 families. Across kingdoms there is no shared consensus axis, because bacterial, archaeal, eukaryotic and mitochondrial rRNAs are separate families (grouped in clans).

### Cited Findings
- RF02541 (bacterial LSU) lists 1730 structures, 10816 species and 52036 sequences, in clan CL00112 with the other LSU/5.8S families. RF02540 (archaeal LSU) lists 115 structures — [Rfam RF02541](https://rfam.org/family/RF02541), [Rfam RF02540](https://rfam.org/family/RF02540) **[V]**
- RF02543 (eukaryotic LSU) is mapped to cryo-EM/X-ray entries in RNA 3D Hub, for example 6QZP — [BGSU NR release list](https://rna.bgsu.edu/rna3dhub/nrlist/release/3.92/all) **[V, snippet]**
- Rfam motif pages list RF01959 (archaeal SSU), RF01960 (eukaryotic SSU), RF02540, RF02541 and RF02542 (microsporidia SSU) among families with annotated 3D modules — [Rfam motif RM00003](https://rfam.org/motif/RM00003) **[V]**
- Rfam per-structure pages colour mapped families on the 3D structure — [rfam.org/structure/5fk3](https://rfam.org/structure/5fk3/av) **[V]**. The rfaRm R package exposes residue-to-CM-position correspondences per family — [rfaRm rfamPDBMapping](https://rdrr.io/github/LaraSellesVidal/rfaRm/man/rfamPDBMapping.html) **[V]**
- Rfam 15 (*NAR* 2025, 53:D258–D267, https://doi.org/10.1093/nar/gkae1023): 65 families enhanced using experimentally determined 3D structures, and R-scape refinement of 26 families. The current release is Rfam 15.1 (January 2026, 4227 families) — [PMC11701678](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11701678/), [rfam.org](https://rfam.org/) **[V]**
- Rfam's 2022 3D update talk reported that 126 families map to 3D structures. Most families were built before their structures were solved, so Rfam is realigning PDB sequences in seeds and reviewing SS_cons against 3D — [Benasque 2022 slides](https://www.benasque.org/2022rna/talks_contr/1610_Rfam_3D_and_R-scape_updates_-_Benasque_2022.pdf) **[V]**

### Inferences
- Known issues, reasoned from the above and from this project's own RF00177 experience rather than from literature:
  - rRNA CMs are huge, and PDB chains are often truncated or modelled incompletely, which gives local-end and insert-state residues.
  - Separate per-kingdom families mean drug-site positions need a cross-family correspondence.
  - Mitochondrial rRNAs are highly reduced and divergent, and may map poorly or only partly.
  - Whether the rRNA families are among the 65 3D-enhanced families in Rfam 15 is unconfirmed.

### Gaps
- I did not find which rRNA families were 3D-revised in Rfam 15, or any published assessment of Rfam rRNA CM mapping quality for PDB chains.
- Mitochondrial rRNA Rfam family accessions and their PDB coverage were not checked.

## Q5. Is there work on motif-level conservation or variation at drug sites that explains selectivity?

### Takeaway
Selectivity is explained in the literature at the **nucleotide or identity level**: A1408G in eukaryotes for aminoglycosides, the mitochondrial A1555G deafness mutation, and PTC geometry and solvent exposure in mitoribosomes. Recent cryo-EM shows *context* differences, not just sequence. I found no work that frames selectivity as conservation or variation of **recurrent 3D motifs** (for example, a Motif Atlas loop class retained or lost) across bacteria, mitochondria and eukaryotes.

### Cited Findings
- Chloramphenicol and linezolid on human mitoribosomes (*NAR* 2025, 53(3) gkaf046): the PTC is well conserved between bacterial and mitochondrial ribosomes, which explains susceptibility. The mitoribosomal PTC is more solvent-exposed than in E. coli, which changes how linezolid's tail sits, and inhibition is context-specific. The authors propose this could guide less mitotoxic designs — [NAR gkaf046](https://academic.oup.com/nar/article/53/3/gkaf046/8001965) **[V]**
- Genetic and structural basis of aminoglycoside ototoxicity via the mitoribosomal A-site: Hobbie et al., *PNAS* 2008, https://doi.org/10.1073/pnas.0811258106 **[K]**. Aminoglycosides on the eukaryotic 80S: Prokhorova et al., *PNAS* 2017, https://doi.org/10.1073/pnas.1715501114 **[K]**. Inhibitors on the yeast 80S, a cross-kingdom site comparison: Garreau de Loubresse et al., *Nature* 2014, https://doi.org/10.1038/nature13737 **[K]**
- Antibiotic target sites are among the least diverse ribosomal regions, so resistance may also involve less conserved regions — [PMC10039652](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10039652/) **[V]**

### Inferences
- The candidate novelty is real but narrow. The **sites themselves and their per-nucleotide selectivity determinants are well known**. What is not available is a uniform, reproducible annotation layer with these features:
  - drug-contacting nucleotides on family consensus positions;
  - the FR3D base pairs, triples and stacks and the Motif Atlas loops at those positions;
  - per-structure observation counts with numerators and denominators;
  - side-by-side comparison of bacterial, mitochondrial and cytosolic eukaryotic structures.
  This would be useful mainly as an *integrative resource* (a PDBe-KB-like view) and for hypothesis generation, for example finding where selectivity coincides with a change in a non-canonical pair or triple and not just a base identity. It is unlikely to be seen as a new biological discovery on its own.
- The main risks to novelty claims are RiboVision2 (cross-kingdom alignments, any wwPDB structure) and possible BGSU or RNA 3D Hub ligand layers. Both should be inspected before claiming "first".

### Gaps
- No source was found that compares Motif Atlas loop instances at antibiotic sites across kingdoms. This would need a targeted search of the RNA 3D Hub and the RiboVision2 feature lists, and of 2024–2026 bioRxiv.
- Lin et al. 2018, Wilson 2014, Hobbie 2008, Prokhorova 2017 and Garreau de Loubresse 2014 are cited from memory with DOIs and should be verified.
