# Databases, datasets and large-scale analyses of RNA–small molecule complexes and RNA binding pockets

Research date: 2026-10-23. Abstracts and full texts were read through Europe PMC and NCBI E-utilities. HARIBOSS methods come from the Bioinformatics article page and the live HARIBOSS site. "Reachability" below means an HTTP probe from our network on 2026-10-23. A probe that timed out is weak evidence, because it may be blocked by region or firewall, so it is not proof that a site is dead.

## 1. What resources exist and what do they contain?

### Takeaway
Only one resource is a structure-based (PDB) RNA–small-molecule database that is curated and still maintained: HARIBOSS / HARIBOSS++. It has about 2,131 PDB entries, 550 compounds and 5,002 pockets, is licensed CC0 and is updated regularly. Most other resources list ligands from the literature or from screens without 3D structures (R-BIND, SMMRNA, NALDB, Inforna, ROBIN, RNAInter), or are small, frozen ML benchmarks (RNAsite TR60/TE18, RLBind, RBind). BioLiP2 treats RNA as a *ligand* of proteins, not as a receptor. NAKB is the maintained nucleic-acid PDB portal that carries Rfam annotations.

### Cited Findings
**Structure-based (PDB) collections**
- **HARIBOSS** (Panei, Torchet, Ménager, Gkeka, Bonomi; *Bioinformatics* 38:4185, 2022; Institut Pasteur and Sanofi). This is a curated collection of RNA–small-molecule structures solved by X-ray, NMR and cryo-EM, with physico-chemical annotations of ligands and pockets — [Europe PMC 35799352](https://europepmc.org/abstract/MED/35799352); [doi:10.1093/bioinformatics/btac483](https://doi.org/10.1093/bioinformatics/btac483)
- Contents of the HARIBOSS 2022 release (May 2022): the redundant set had 716 structures, 1,226 pockets and 267 unique ligands; the non-redundant set had 484 structures and 676 pockets. Methods were about 83% X-ray, 11% cryo-EM and 6% NMR. Of the pockets, 67% sit on one chain, 29% on two chains and 4% on three or more — [Bioinformatics full text](https://academic.oup.com/bioinformatics/article/38/17/4185/6633309)
- **HARIBOSS++** (Marengo, Chapeaublanc, Ménager, Gkeka, Bonomi; bioRxiv 2026-09-13, doi 10.64898/2026.09.07.749797). It adds pockets where the ligand also contacts other biomolecules, analysis "by family", detailed 3D characterisation, ligand chemical space, experimental binding affinities and MD data — [bioRxiv / Europe PMC PPR1319252](https://doi.org/10.64898/2026.09.07.749797)
- The live HARIBOSS site on 2026-10-23 showed 2,131 PDBs, 550 compounds and 5,002 pockets. By environment: 624 RNA-only, 1,239 RNA/protein, 20 RNA/DNA and 248 RNA/DNA/protein complexes. 1,953 complexes have an "equivalence class", and there are MD simulations for 11 complexes (12 pockets). Data come from the "GitHub COSB Lab" repository (commit 7cf281b1, imported 2026-10-19) — [HARIBOSS About/General](https://hariboss.pasteur.cloud/about-general)
- HARIBOSS pocket statistics by "RNA family": rRNA 3,092, riboswitch 445, synthetic construct 376, unknown 333, aptamer 188, ribozyme 170, tRNA 146, multiple 120, RNA/DNA hybrid 61, TAR 41, mRNA 13, telomerase RNA 10, snRNA 5. 127 pockets have binding-affinity data. By environment: RNA 2,081, RNA/Protein 1,955, RNA/Other 712, RNA/DNA/Protein 234, RNA/DNA 20 — [HARIBOSS About/Pocket](https://hariboss.pasteur.cloud/about-pocket)

**Literature or screening ligand databases (mostly no 3D structures)**
- **R-BIND** (Hargrove lab, Duke). It is the first database of non-ribosomal RNA ligands to require bioactivity for inclusion (2017). The website rbind.chem.duke.edu has been online since 2019 and was updated every two years — [R-BIND 2019, ACS Chem Biol](https://doi.org/10.1021/acschembio.9b00631); [R-BIND 2.0, ACS Chem Biol 2022](https://doi.org/10.1021/acschembio.2c00224)
- R-BIND 2.0 covers reports published 2018-06-01 to 2019-12-31. The small-molecule set grew from 97 to 153. Ligands are now classed SM / MV (multivalent) / LM (large molecule). RNA target sequences and secondary-structure motifs were added, along with a motif search that works without knowing the secondary structure in advance — [R-BIND 2.0 full text, PMC9343015](https://pmc.ncbi.nlm.nih.gov/articles/PMC9343015)
- **SMMRNA** (Mehta et al., *NAR* 2014). About 770 unique ligands with Kd, Ki, IC50, ΔTm and descriptors, plus 2D/3D structures and RNA images, at smmrna.org — [Europe PMC 24163098](https://doi.org/10.1093/nar/gkt976)
- **NALDB** (Kumar Mishra & Kumar, *Database* 2016). More than 3,500 ligand entries for all nucleic-acid types (DNA and RNA), with Ki, Kd, IC50 and up to 255 conformers per ligand — [doi:10.1093/database/baw002](https://doi.org/10.1093/database/baw002)
- **Inforna 2.0** (Disney lab, *ACS Chem Biol* 2016). 1,936 RNA motif–small molecule interactions, 244 unique small molecules and 1,331 secondary-structure motifs. It is a sequence-based design tool, not a 3D resource — [doi:10.1021/acschembio.6b00001](https://doi.org/10.1021/acschembio.6b00001)
- **RNALigands** (Sun, Yang, Zhang; *RNA* 2022). It links RNA secondary-structure motifs to the small molecules that bind them. It derives 386 motif–ligand pairs from the PDB and adds R-BIND (104 pairs) and other sets. The web server takes a sequence, predicts the secondary structure and searches for similar motifs. Code is on GitHub (SaisaiSun/RNALigands) — [PMC8906548](https://doi.org/10.1261/rna.078889.121)
- **ROBIN** (Repository Of BInders to Nucleic acids; Yazdani et al., *Angew Chem* 2023, Schneekloth lab). The authors ran 36 small-molecule-microarray screens of nucleic acids against 24,572 compounds (1,627,072 interactions assayed). They found 2,003 RNA-binding small molecules, which they call the "largest fully public, experimentally derived library of its kind". Code and data are on GitHub (ky66/ROBIN) — [PMC9992102](https://doi.org/10.1002/anie.202211358)
- **RNAInter v4.0** (*NAR* 2022). An RNA interactome with more than 47 million entries from the literature and six databases. It includes RNA–compound interactions and uses a redefined confidence score. It is not structure-based — [doi:10.1093/nar/gkab997](https://doi.org/10.1093/nar/gkab997)

**General PDB-derived resources**
- **BioLiP2** (Zhang et al., *NAR* 2024). It holds biologically relevant **ligand–protein** interactions. Compared with BioLiP it covers far more nucleic acid–protein interactions, but RNA appears as a *ligand* ("RNA ligands" browse view) and the receptors are proteins. Datasets come redundant and non-redundant (receptors below 90% identity by CD-HIT), and the curation code is open — [PMC10767969](https://doi.org/10.1093/nar/gkad630)
- **NAKB** (Nucleic Acid Knowledgebase, *NAR* 2024). It replaces the NDB. It covers NA-containing structures from X-ray, NMR and EM and is "updated weekly". Annotations include BGSU RNA equivalence class, RNAcentral ID and **Rfam ID**, PDB nonpolymer ligands and PANTHER classes. Search can filter on "ligand present by MW cutoff". PDB metadata come from the RCSB API. The paper is CC BY 4.0 — [PMC10767938](https://doi.org/10.1093/nar/gkad957)
- **SPIRAL** (Lu & Wang, bioRxiv 2026). A DSSR-based curated database of 1,098 RNA–small-molecule PDB structures with 1,137 binding events in 6 functional RNA categories. It defines 6 binding modes by unsupervised clustering and a Composite Binding Quality Score. It has 275 entries with affinity data — [Europe PMC PPR1237231](https://doi.org/10.64898/2026.05.19.726393)

**Binding-site prediction and docking benchmarks**
- **RNAsite** (Su, Peng, Yang; *Bioinformatics* 2021). It defines the **TR60 / TE18** benchmark: 78 RNA chains (RB78), split into 60 training RNAs from 42 clusters and 18 test RNAs from 15 clusters. It also uses the older RB19 set. The web server is at yanglab.nankai.edu.cn/RNAsite — [PMC8034527](https://doi.org/10.1093/bioinformatics/btaa1092)
- **RBind** (Wang et al., *Bioinformatics* 2018). A network-based predictor of binding sites. Its average accuracy is 0.82 on RNA–ligand and 0.63 on RNA–protein test sets — [doi:10.1093/bioinformatics/bty345](https://doi.org/10.1093/bioinformatics/bty345)
- **RLBind** (Wang et al., *Brief Bioinform* 2023). A CNN that combines sequence and local-neighbour channels. Datasets and code are on GitHub (KailiWang1/RLBind) — [doi:10.1093/bib/bbac486](https://doi.org/10.1093/bib/bbac486)
- Later predictors still report results on TE18 and RB9 (for example RNetsite, MultiModRLBP). MultiModRLBP trained on 851 RNA–ligand interactions from RNAglib plus the RLBind training set — [CoBRA preprint 2025](https://www.biorxiv.org/content/10.1101/2025.09.03.673903.full.pdf); [PMC12417085](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12417085/)
- **Rsite** (Zeng et al., *Sci Rep* 2015). It predicts functional sites from extreme points of the distance curve in a 3D structure. It was shown on only 2 RNAs (tRNA-Lys and the Diels–Alder ribozyme) — [doi:10.1038/srep09179](https://doi.org/10.1038/srep09179)
- **RLDOCK** (Sun et al., *JCTC* 2020). Tested on 230 RNA–ligand bound structures. Within its top 3 poses it placed 27.8%, 58.3% and 69.6% of ligands within 1, 2 and 3 Å RMSD — [doi:10.1021/acs.jctc.0c00798](https://doi.org/10.1021/acs.jctc.0c00798)
- **rDock** (Ruiz-Carmona et al., *PLoS Comput Biol* 2014). It evolved from RiboDock and docks to proteins and nucleic acids. It is licensed LGPL, and its test sets are distributed with it — [doi:10.1371/journal.pcbi.1003571](https://doi.org/10.1371/journal.pcbi.1003571)
- **Docking pose benchmark** (Agarwal, Smith; *JCIM* 2023). 173 RNA–small-molecule crystal structures. With the search space defined by a known ligand, rDock succeeded in 48% and Vina in 63% of cases (the abstract wording is ambiguous about which number goes with which program). In blind (large-box) docking both were about 27%. AD4 did poorly — [doi:10.1021/acs.jcim.3c01533](https://doi.org/10.1021/acs.jcim.3c01533)
- **DrugPred_RNA** (Rekand & Brenk, *JCIM* 2021). A druggability predictor trained on protein pockets, using only descriptors that can also be computed for RNA. The RNA set started from 1,084 PDB structures; after removing entries whose only ligands were detergents, buffer salts or crystallisation components, 427 unique entries remained. There is also a separate ribosome set — [PMC8389535](https://doi.org/10.1021/acs.jcim.1c00155)
- **DRLiPS** (Krishnan et al., *NAR* 2025). Takes 819 binding sites from 399 RNA–ligand complexes, including 51 ribosomal A-site complexes, and clusters them **by Rfam annotation plus manual classification** into 56 non-redundant positive sites. Negative sites are sampled by backbone motif search, fpocket and blind docking. On the external test set: F1 0.70, precision 0.61, recall 0.73. Training data come from the PDB and HARIBOSS — [PMC11963762](https://doi.org/10.1093/nar/gkaf239); [DRLiPS details](https://web.iitm.ac.in/bioinfo2/DRLiPS/Details.html)
- **fpocketR** (Veenbaas et al., *PNAS* 2025; *ACS Chem Biol* 2026). fpocket 4.0.3 retuned for RNA on a curated subset of HARIBOSS. Starting from 538 non-redundant HARIBOSS structures, the authors kept RNAs of ≤200 nt at ≤3.5 Å resolution with ligands of QED ≥0.3, giving 48 complexes. Of these, 38 had deep pockets, split 13 training / 25 test. fpocketR found 100% of known sites with a minimum PPV of 78% — [PNAS PMC12054788](https://doi.org/10.1073/pnas.2422346122); [ACS Chem Biol 2026](https://doi.org/10.1021/acschembio.5c00805)
- **SHAMAN** (Bonomi lab, *Nat Commun* 2024). Finds small-molecule binding sites in RNA conformational ensembles using MD and fragment probes — [doi:10.1038/s41467-024-49638-7](https://doi.org/10.1038/s41467-024-49638-7)

### Inferences
- For our purposes the only reusable 3D source for cross-checking our PDBe Arpeggio ligand contacts is HARIBOSS: its data are CC0, and its data repository is on GitHub with a commit hash. NAKB is the reusable source of chain-level Rfam IDs.
- R-BIND, Inforna, ROBIN, SMMRNA and NALDB are chemistry and bioactivity resources. They could add "known bioactive binders" annotations on a family page, but they give no residue-level contacts.
- TR60/TE18 and the HARIBOSS-derived sets (fpocketR's 13/25, DRLiPS's 56 sites) are where riboswitch pockets would serve as external benchmarks for our fpocketR pocket tracks.

### Gaps
- **PDBbind nucleic acid–ligand subset**: I found no primary source for its current size or curation. Search results only gave PDBbind v2020 protein–nucleic-acid subsets (186 complexes in a derived set), which are a different thing — [arXiv SVSBI](https://arxiv.org/pdf/2212.13617). The PDBbind release notes need to be checked by hand.
- I could not confirm the exact composition of the RLBind training set (whether it reuses TR60/TE18) from primary text.
- The HARIBOSS++ full text (bioRxiv) could not be retrieved because of rate limiting. Its exact family-mapping method (Rfam accession or coarse category) and its updated curation rules were not checked against the paper.
- I did not cover RNAglib datasets, RNAmigos or RsiteDB in depth.

## 2. How do they curate (ligand filtering, redundancy, resolution, pocket definition, Rfam mapping)?

### Takeaway
Curation choices differ a lot. HARIBOSS keeps "drug-like" ligands by mass (160–1,000 Da) and element filters, uses a 6 Å interaction radius, and deliberately keeps buffers and crystallisation agents. BioLiP2 and RNAsite use the BioLiP artifact list (463 ligands). Binding-site benchmarks use a 4 Å nucleotide contact cut-off. Redundancy is usually cut by sequence clustering (CD-HIT 90–95%) or TM-score. True Rfam mapping is rare: DRLiPS clusters by Rfam, NAKB annotates Rfam IDs, and HARIBOSS uses coarse RNA categories.

### Cited Findings
- **HARIBOSS: retrieval and parsing.** Entries are fetched from the PDBe search API when they contain at least one RNA and one ligand. Molecules are typed as RNA, DNA, protein, ion or water. Modified residues are recognised by their covalent bonds, and anything left is a ligand — [HARIBOSS, Bioinformatics](https://academic.oup.com/bioinformatics/article/38/17/4185/6633309)
- **HARIBOSS: "drug-like" filter.** Mass 160–1,000 Da, at least one carbon, and only C, H, N, O, Br, Cl, F, P, Si, B, S or Se. Atoms within **6 Å** of the ligand count as interacting, and RNA chains with fewer than 10 atoms inside that radius are dropped. The 2022 version kept only ligands that contact RNA alone, so ligands that also touch protein or DNA were removed. **Buffers and crystallisation enhancers were deliberately not excluded**, because they may help fragment-based design. The paper gives no explicit additive list — [HARIBOSS, Bioinformatics](https://academic.oup.com/bioinformatics/article/38/17/4185/6633309)
- **HARIBOSS: redundancy.** CD-HIT-EST at 90% identity on RNA chains, with chains that share a ligand linked into one cluster. Within clusters, entries are subclustered by ligand and then by pocket-backbone RMSD with a 1.5 Å hierarchical cut-off. The representative is the entry nearest the cluster centre with the best resolution. HARIBOSS sets no resolution cut-off; typical values are 2.8 Å for X-ray and 3.0 Å for cryo-EM — [HARIBOSS, Bioinformatics](https://academic.oup.com/bioinformatics/article/38/17/4185/6633309)
- **HARIBOSS: pocket definition.** Two methods are used. Schrödinger SiteMap runs on a 12 Å region with PROPKA at pH 7.4 and OPLS3 minimisation, with a 6 Å sitebox; pockets are kept if they hold at least 20% of the ligand atoms. mkgridXf cavities use 1 and 8 Å probes and watershed sub-cavities, also kept at ≥20% ligand-atom occupancy — [HARIBOSS, Bioinformatics](https://academic.oup.com/bioinformatics/article/38/17/4185/6633309)
- **HARIBOSS++** extends the pocketome to sites where the ligand also contacts other biomolecules, and adds "enhanced analysis of RNA molecules by family" — [HARIBOSS++ preprint](https://doi.org/10.64898/2026.09.07.749797). The live site's "RNA families" are coarse categories (rRNA, riboswitch, aptamer, tRNA, ...), not Rfam accessions — [HARIBOSS pocket page](https://hariboss.pasteur.cloud/about-pocket)
- **BioLiP2 artifact rule.** The list has 463 artifact ligands: 353 from the first BioLiP and 110 added since. A listed ligand is discarded if it appears more than 15 times in the structure or if its site has only two consecutive residues. Otherwise the title and abstract of the primary citation (PubMed ID from mmCIF or SIFTS) are checked to judge whether it is biologically relevant — [BioLiP2, PMC10767969](https://doi.org/10.1093/nar/gkad630)
- **RNAsite curation.** Started from 1,673 RNA chains. Removed chains shorter than 20 or longer than 1,500 nt, structures whose only small molecules are crystallisation additives (BioLiP list), and chains with no interacting ligand. A nucleotide counts as binding if any atom is **< 4 Å** from the ligand, as in RBind. 712 chains remained. Redundancy was cut by RNA-align TM-score (>0.3 counts as redundant; for example the 22 TPP riboswitches collapse to one), leaving RB78. Clusters at 30% sequence similarity (cd-hit-est and BLASTclust) set the train/test split — [RNAsite, PMC8034527](https://doi.org/10.1093/bioinformatics/btaa1092)
- **RNALigands curation.** Kept entries with "biologically relevant" ligands (BioLiP), RNA–ligand minimum distance ≤4 Å and CD-HIT ≤95% identity, and removed identical ligand–binding-mode duplicates. Motifs from RNAView were assigned to ligands at <4 Å, giving 386 motif–ligand pairs — [RNALigands, PMC8906548](https://doi.org/10.1261/rna.078889.121)
- **R-BIND inclusion rules.** MW < 2,000 Da, binds a *non-ribosomal* RNA in vitro without covalent bonds, and shows activity in cell culture or animals. Ribosomal ligands are excluded by design — [R-BIND 2.0, PMC9343015](https://pmc.ncbi.nlm.nih.gov/articles/PMC9343015)
- **DrugPred_RNA.** Removes detergents, buffer salts and crystallisation components (the rejected three-letter codes are in its supplement). Binding sites are cut out as all residues within 15 Å of the ligand. Ribosome ligands were also checked by eye for buffer components — [PMC8389535](https://doi.org/10.1021/acs.jcim.1c00155)
- **Hewitt et al. 2019.** RNA chains of 12–700 residues with at most 2 protein chains (1,552 structures). "Liganded" means a small molecule with more than 7 heavy atoms, to leave out metal ions and solvents, giving 300 structures. Pockets were found with ICM PocketFinder — [PMC8283815](https://doi.org/10.1016/j.bmc.2019.04.010)
- **fpocketR training set.** Built from non-redundant HARIBOSS with ≤200 nt, ≤3.5 Å and ligand QED ≥0.3 — [PNAS 2025](https://doi.org/10.1073/pnas.2422346122)
- **DRLiPS.** Clusters binding sites by Rfam annotation and manual classification, and checks non-redundancy with RMscore — [PMC11963762](https://doi.org/10.1093/nar/gkaf239)

### Inferences
- The thresholds in use are 4 Å (benchmarks), 5 Å (RCSB neighbours, which we use as a fallback) and 6 Å (HARIBOSS). Our Arpeggio contacts are typed interactions, so they are stricter than any of these. Any cross-check against HARIBOSS should expect HARIBOSS to have larger residue sets.
- HARIBOSS keeps buffers and additives on purpose. If we use it as a ligand allow-list, we still need our own additive filter, and the BioLiP 463-ligand list is a citable, reusable one. The way R-BIND rules out ribosomal and covalent ligands is a model for flagging rather than deleting, which fits our rule to keep flagged observations.
- Collapsing all TPP riboswitches into one chain (RNAsite) shows that ML benchmarks remove exactly the within-family ensembles our project aggregates. Per-family conservation is something those benchmarks cannot supply.

### Gaps
- I did not obtain the BioLiP artifact list itself; it is in the BioLiP download area. HARIBOSS's exact ion and water definitions were not checked either.
- Whether HARIBOSS++ maps to Rfam accessions is unconfirmed (see section 1).

## 3. What did large-scale analyses find about RNA vs protein pockets and about ligand chemistry?

### Takeaway
Analyses agree that RNA pockets look like protein pockets in volume and buriedness but are clearly more polar or less hydrophobic. "Ligandable" pockets are rarer on RNA, and only about 50 pockets that bind drug-like ligands have been seen in structures. RNA-binding ligands tend to be more nitrogen-rich, more aromatic and more rod-shaped than protein binders, and R-BIND's "RNA-privileged chemical space" has held across updates.

### Cited Findings
- **Hewitt, Calabrese, Schneekloth (*Bioorg Med Chem* 2019)**, 1,552 RNA structures compared with all PDB proteins:
  - Pocket volume medians: proteins 175 Å³, RNA 187 Å³.
  - Buriedness medians: proteins 0.75, RNA 0.72.
  - Hydrophobicity medians: proteins **0.49**, RNA **0.35**.
  - DLID (drug-like density) medians: proteins −0.84, RNA −1.12, so "ligandable pockets are considerably more frequent on proteins".
  - For liganded pockets: volume 352 Å³ (protein) vs 280 Å³ (RNA); buriedness 0.81 vs 0.77.
  - 90% of the 300 liganded RNAs had several pockets.
  - [PMC8283815](https://doi.org/10.1016/j.bmc.2019.04.010)
- **HARIBOSS (2022), SiteMap.** RNA pockets have donor/acceptor character and hydrophilicity similar to protein benchmarks, but are **less hydrophobic and more solvent-exposed**. Overall, 67% are ligandable (SiteScore ≥0.8), while for druggability (Dscore) 35% are druggable, 24% difficult and 41% undruggable. Cavities under 100 Å³ were never druggable, and the largest druggable pocket was about 230 Å³ — [HARIBOSS, Bioinformatics](https://academic.oup.com/bioinformatics/article/38/17/4185/6633309)
- **HARIBOSS ligand properties.**
  - Mass ranges from 162 to 972 Da, and 82% are under 600 Da.
  - About half are above the Veber TPSA limit of 140 Å², and about 45% are above the H-bond donor plus acceptor limit.
  - 12.4% have more than 10 rotatable bonds.
  - About 52% have poor predicted Caco-2 permeability.
  - Frequent binders are antibiotics and polar endogenous metabolites.
  - [HARIBOSS, Bioinformatics](https://academic.oup.com/bioinformatics/article/38/17/4185/6633309)
- **fpocketR (PNAS 2025).** Reviewing earlier work, the authors summarise that RNA and protein pockets have "similar bulk properties of volume, buriedness, and solvent accessibility", while RNA pockets "have been proposed to be less hydrophobic and to bind more rod-like ligands". They say earlier studies were biased by a "preponderance of simple organic" ligands. Only "roughly 50" pockets that bind drug-like ligands have ever been seen in structures. fpocket's polarity filters (A, p) did not work on RNA — [PNAS PMC12054788](https://doi.org/10.1073/pnas.2422346122)
- **DrugPred_RNA.** Among RNA-set ligands with QED ≥0.67, 12 (67%) sat in pockets predicted druggable, compared with 37% of all metal-free sites. Drug-like ligands are therefore enriched in druggable pockets. Examples are the malachite green aptamer and the FMN riboswitch (ligand 6YG). The tool advises scoring several structures of the same site, because predictions change with conformation — [PMC8389535](https://doi.org/10.1021/acs.jcim.1c00155)
- **R-BIND chemical space.**
  - Compared with protein-targeting ligands, R-BIND small molecules have distinctive physicochemical, structural and shape properties ("RNA-privileged chemical space").
  - These differences held after a 50% growth of the set (97 to 153), and PMI analysis shows the set remains "enriched in rodlike shapes".
  - By target class: G-quadruplex binders average about 5 rings and more than 4 aromatic rings; dsRNA binders are more lipophilic; bulge binders have higher MW, more rotatable bonds and more surface area.
  - [R-BIND 2.0, PMC9343015](https://pmc.ncbi.nlm.nih.gov/articles/PMC9343015)
- **ROBIN.** RNA binders and protein binders have similar drug-like properties. Even so, descriptors for nitrogen-containing and aromatic rings, van der Waals surface area and topological charge separate the two groups. In TMAP space ROBIN binders spread across many scaffolds — [ROBIN, Angew Chem 2023](https://doi.org/10.1002/anie.202211358)
- **Inforna 2.0.** Chemotype analysis identifies features that are "privileged" for RNA binding — [doi:10.1021/acschembio.6b00001](https://doi.org/10.1021/acschembio.6b00001)
- **SPIRAL (2026 preprint).** Unsupervised clustering finds 6 binding modes, and their distribution depends strongly on the functional class of the RNA. Riboswitches score highest on binding quality, and C2′-endo pucker count and buried contact surface area best predict affinity across 275 entries — [PPR1237231](https://doi.org/10.64898/2026.05.19.726393)

### Inferences
- These results suggest which descriptors to show on our riboswitch pocket tracks: buriedness or SASA (we already compute SASA), polarity and pocket volume. They also suggest putting the context next to each other: RNA pockets are about as buried as protein pockets but less hydrophobic.
- The finding that conformation changes the predictions (DrugPred_RNA) and SHAMAN's ensemble approach both support our per-family aggregation over many instances. Pocket presence in our data should be reported as a fraction of chains (numerator over denominator), not taken from a single structure.
- Most of the "about 50" drug-like pockets are riboswitch or aptamer sites. So our riboswitch families (RF00167, RF00162, RF00059) are where those pocket definitions apply most directly.

### Gaps
- No study I found compared RNA and protein pockets with the same pocket-detection tool and a non-redundant set controlled by family. All of them pool rRNA with riboswitches.
- The exact R-BIND descriptor values (for example nitrogen counts, Fsp3) are in its supplementary tables, which I did not retrieve.

## 4. Update/maintenance status, licences, programmatic access

### Takeaway
The maintained resources are HARIBOSS (CC0, regular automated imports, data on GitHub), NAKB (updated weekly), BioLiP2 (open curation code), and DRLiPS and fpocketR (code). R-BIND, SMMRNA, NALDB, Inforna and RNAInter appear static since about 2016–2022, and several did not respond to our network probe.

### Cited Findings
- **HARIBOSS.** "Regularly updated"; the 2022 paper says monthly. The last data import shown on the site was 2026-10-19 (commit 7cf281b1 of the COSB Lab GitHub repository). Data files are under **CC0 1.0** — [HARIBOSS home](https://hariboss.pasteur.cloud/); [HARIBOSS policies](https://hariboss.pasteur.cloud/about-policies-logo); [HARIBOSS general](https://hariboss.pasteur.cloud/about-general)
- **NAKB** is updated weekly — [PMC10767938](https://doi.org/10.1093/nar/gkad957)
- **BioLiP2** provides redundant and non-redundant downloads plus curation and web source code — [PMC10767969](https://doi.org/10.1093/nar/gkad630). A BioLiP page lists 182,239 RNA-ligand entries, meaning RNA bound to protein — [zhanggroup BioLiP](https://www.zhanggroup.org/BioLiP/)
- **R-BIND.** Its last major paper is 2.0, from 2022, covering literature to the end of 2019 — [ACS Chem Biol 2022](https://doi.org/10.1021/acschembio.2c00224)
- **Inforna.** No database update since 2.0 (2016) was found. It is free for academic use under a licence agreement — [search summary of Disney lab pages](https://disney.scripps.ufl.edu/software/)
- **ROBIN**: data and models on GitHub — [Angew 2023](https://doi.org/10.1002/anie.202211358)
- **RLBind**: code and data on GitHub — [Brief Bioinform 2023](https://doi.org/10.1093/bib/bbac486)
- **RNALigands**: code on GitHub — [RNA 2022](https://doi.org/10.1261/rna.078889.121)
- **rDock** is LGPL — [PLoS Comput Biol 2014](https://doi.org/10.1371/journal.pcbi.1003571)
- **DRLiPS** is free online, and its paper is CC BY-NC 4.0 — [PMC11963762](https://doi.org/10.1093/nar/gkaf239)
- HTTP probe on 2026-10-23 (my own measurement; codes are HTTP status, and 000 means a timeout):

| Result | Sites |
|---|---|
| 200 (responded) | hariboss.pasteur.cloud, nakb.org, zhanggroup.org/BioLiP, DRLiPS, pdbbind-plus.org.cn, GitHub ROBIN and RLBind |
| 302 (redirect) | RNALigands |
| 404 | NALDB HOME.php |
| 000 (timed out) | rbind.chem.duke.edu, smmrna.org, rnainter.org, yanglab.nankai.edu.cn/RNAsite, pdbbind.org.cn |

### Inferences
- For a reusable, citable, maintained 3D benchmark under a compatible licence, HARIBOSS is the clear choice. Pin the GitHub commit in provenance, as our provenance rules require.
- The NALDB 404 and the timeouts (R-BIND, SMMRNA, RNAInter, RNAsite) should be confirmed from another network before we call these resources defunct.

### Gaps
- No programmatic API was documented for HARIBOSS apart from the web interface and the GitHub data repository; I did not inspect the repository's structure.
- The PDBbind and PDBbind+ licence and nucleic-acid subset were not checked from primary pages.
- I could not determine the current licences for R-BIND, SMMRNA or NALDB.
