# Existing resources and surveys of small molecules bound to ribosomes in the PDB

Compiled 2026-10-24 for the planned study: archive-wide, cross-kingdom ribosome–small-molecule binding sites placed on Rfam CM consensus positions (RF00177, RF02541, RF01960, RF02543) and linked to RNA 3D motifs. Our own inventory, for comparison: about 1,180 usable drug copies in 568 PDB entries with PDBe-Arpeggio contacts.

Access caveats: bioRxiv, Cell Press/ScienceDirect, HAL and PMC web pages were often blocked (HTTP 403/429). Most full texts were read via Europe PMC fullTextXML. Items read only as abstracts are marked as such. "UNVERIFIED" means I could not confirm the point from a primary source.

## Q1. Which resources cover ribosome ligands, and what do they do or not do?

### Takeaway
Several RNA–ligand databases now include ribosomal pockets, and in some they are the largest category: HARIBOSS/HARIBOSS++, SPIRAL, FURNA and a DrugPred_RNA ribosome set. All of them store binding sites in each structure's own chain numbering, or as per-entry descriptors. None places ribosomal contacts on a shared rRNA coordinate system across kingdoms. The ligand-centric chemistry databases (R-BIND, SMMRNA, NALDB, RNALigands) are not structural ribosome resources. The ribosome-specific structural portals (RiboXYZ, RiboVision2, ProteoVision) handle nomenclature, visualisation and phylogeny, not systematic ligand-contact aggregation.

### Cited findings

**HARIBOSS (2022) and HARIBOSS++ (2026 preprint), Bonomi lab, Institut Pasteur**
- HARIBOSS is a curated database of RNA–small-molecule structures (X-ray, NMR, cryo-EM) with drug-like ligands. Each entry has physicochemical annotations of ligands and pockets. Panei FP, Torchet R, Ménager H, Gkeka P, Bonomi M. *Bioinformatics* 38(17):4185–4193 (2022), doi:10.1093/bioinformatics/btac483 — [Europe PMC record](https://europepmc.org/article/MED/35799352); [preprint](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full)
- One downstream user reports that the HARIBOSS release it used had 746 complexes and 1,268 binding sites, with sites defined as residues within 6 Å of the ligand — [arXiv 2310.18985](https://arxiv.org/pdf/2310.18985)
- On the live website today: 550 compounds, 2,131 RNA–small-molecule complexes and 5,002 pockets. The homepage says nothing about ribosomes and gives no update date — [hariboss.pasteur.cloud](https://hariboss.pasteur.cloud/)
- **Ribosomes are included.** In HARIBOSS++, rRNA is the largest RNA family at 60.9% of pockets, which the authors attribute to "experimental bias toward ribosomal structures in the PDB". Other HARIBOSS++ details:
  - Pocket counts (redundant set): 4,718 pockets, split as RNA 2,016, RNA/Protein 1,810, RNA/Other 669, RNA/DNA/Protein 203, RNA/DNA 20.
  - Pocket = residues within 6.0 Å of the ligand, with at least one RNA molecule of at least 10 atoms.
  - Chains are assigned to 12 coarse families by keywords (macromolecule names, organism, title). Rfam is used only as a secondary chain-level cross-check.
  - New features: DSSR secondary structure of pocket nucleotides, SiteMap descriptors, affinities from PDBbind (93 values), MD links via MDDB (11 complexes), BGSU equivalence classes for apo–holo pairing, SMILES/SMARTS search.
  - Source: Marengo M, Chapeaublanc E, Ménager H, Gkeka P, Bonomi M, "HARIBOSS++: An Integrated Platform for RNA-Targeted Drug Design", bioRxiv 2026, doi:10.64898/2026.09.07.749797 — [bioRxiv full text](https://www.biorxiv.org/content/10.64898/2026.09.07.749797v1.full). The non-redundant pocket counts differ between sections of the preprint (2,583 / 2,093 / 2,073 / 2,908).
- What HARIBOSS does not do: it has no ribosome-specific nomenclature, no rRNA residue correspondence across species and no Rfam consensus-position mapping. Rfam serves only as a family cross-check — [HARIBOSS++ preprint](https://www.biorxiv.org/content/10.64898/2026.09.07.749797v1.full)

**SPIRAL (2026 preprint)**
- SPIRAL is a curated set of 1,098 RNA–small-molecule PDB structures with 1,137 binding events, covering PDB depositions through May 2026. Its pipeline is DSSR-based and records stacking, H-bond topology by RNA moiety (base/sugar/phosphate), groove engagement and tertiary-motif context of contacted nucleotides (hairpin loop, junction loop, pseudoknot, multiplet, non-canonical pair and others). Lu, Wang, bioRxiv 2026, doi:10.64898/2026.05.19.726393 — [Europe PMC PMC13228370](https://europepmc.org/article/PMC/PMC13228370)
- **Ribosomal RNA is the largest category: n = 543 of 1,137 binding events (47.8%).**
  - Unsupervised clustering of the interaction fingerprints separates ribosomal sub-sites.
  - Cluster R-C1 (n = 208) is 98% aminoglycosides at the 30S decoding A-site.
  - A second large ribosomal cluster is the PTC/exit tunnel: macrolides/ketolides 26%, lincosamides 14%.
  - Other ribosomal entries spread over clusters without a distinctive signature, for example tetracyclines.
  - Contacts are defined as within 4.8 Å (the Å values were partly garbled in the extracted text, so check the paper).
  - Categories were cross-checked against Rfam.
  - Source: [PMC13228370](https://europepmc.org/article/PMC/PMC13228370)
- SPIRAL is the closest existing "archive-wide quantitative interaction-fingerprint" analysis that includes ribosomes. It works per binding event with structure-local descriptors. I found no text indicating a common rRNA numbering, Rfam consensus positions or cross-kingdom comparison — [PMC13228370](https://europepmc.org/article/PMC/PMC13228370)

**FURNA (2024)**
- FURNA provides functional annotations for all PDB RNA chains of 10 nt or more: GO, EC, ligand-binding sites, Rfam families, RNAcentral, ATtRACT motifs. Ligand interactions are computed on the biological assembly. A nucleotide counts as ligand-binding with 2 or more intermolecular atomic contacts. Unlike BioLiP, FURNA does not filter "biologically irrelevant" ligands; only water is removed. Zhang C, Freddolino L, *PLoS Biol* 22(7):e3002476 (2024), doi:10.1371/journal.pbio.3002476 — [Europe PMC PMC11309384](https://europepmc.org/article/PMC/PMC11309384)
- rRNAs dominate: 58.0% of RNAs with GO molecular-function terms are rRNAs. Rfam families are assigned from the Rfam–PDB mapping, or by cmscan if a chain is absent from it — [PMC11309384](https://europepmc.org/article/PMC/PMC11309384)
- Inference from the full text: binding sites are listed per chain in that chain's own numbering. Rfam is a family label, not a coordinate system — [PMC11309384](https://europepmc.org/article/PMC/PMC11309384)

**BioLiP2 (2024)**
- BioLiP2 extracts biologically relevant ligand–**protein** interactions from mmCIF. Compared with BioLiP, it widens nucleic-acid coverage by including all DNA/RNA–protein interactions (nucleic acids as ligands of protein receptors). Zhang C, Zhang X, Freddolino L, Zhang Y, *NAR* 52(D1):D404–D412 (2024), doi:10.1093/nar/gkad630 — [Europe PMC PMC10767969](https://europepmc.org/article/PMC/PMC10767969)
- BioLiP2 is receptor-protein centric, so rRNA–drug contacts are not its subject. The RNA-receptor counterpart from the same group is FURNA — [PMC11309384](https://europepmc.org/article/PMC/PMC11309384) (FURNA says it was "inspired by BioLiP2")

**R-BIND / R-BIND 2.0 (Hargrove lab)**
- R-BIND is a database of bioactive RNA-targeting small molecules with analysis of their RNA secondary-structure motifs. Morgan BS et al., *ACS Chem Biol* 2019, doi:10.1021/acschembio.9b00631; Donlic A et al., *ACS Chem Biol* 17(6):1556–1566 (2022), doi:10.1021/acschembio.2c00224 — [Europe PMC search record](https://europepmc.org/article/MED/35594415)
- **UNVERIFIED (could not open the full text):** R-BIND's published inclusion criteria are widely understood to exclude ribosome-targeting antibiotics. Neither abstract mentions ribosomes — [Europe PMC](https://europepmc.org/article/MED/35594415)

**RNALigands (2022)**
- RNALigands is a database of RNA secondary-structure motif–ligand pairs: 386 pairs from the PDB, 67 from R-BIND and 388 from Inforna. It keeps nucleotides from a single RNA chain only and works on secondary structure, not 3D. Sun S, Yang J, Zhang Z, *RNA* 28(2):115–122 (2022), doi:10.1261/rna.078889.121 — [Europe PMC PMC8906548](https://europepmc.org/article/PMC/PMC8906548)
- Ribosomes appear in the text only as background examples. There is no ribosome-specific content — [PMC8906548](https://europepmc.org/article/PMC/PMC8906548)

**SMMRNA (2014) and NALDB (2016)**
- SMMRNA collects experimentally validated small-molecule modulators of RNA with 2D ligand coordinates. Its introduction cites rRNA as a validated target (aminoglycosides, macrolides, tetracyclines, oxazolidinones). It is a literature or chemistry database, not a structural contact resource. Mehta A et al., *NAR* 42:D132 (2014), doi:10.1093/nar/gkt976 — [PMC3965028](https://europepmc.org/article/PMC/PMC3965028)
- NALDB has more than 3,500 nucleic-acid-binding small-molecule entries, ligand-centric. Kumar Mishra S et al., *Database* 2016:baw002, doi:10.1093/database/baw002 — [PMC4761111](https://europepmc.org/article/PMC/PMC4761111)

**DrugPred_RNA (2021): a PDB ribosome binding-site dataset**
- The authors built a ribosome set by querying the PDB for "ribosome": 497 unique PDB IDs, 613 pockets (metal-free) and 546 (metal-containing), 217 unique ligands. Ribosomal pockets could include ribosomal proteins.
  - 31% of ribosomal pockets were predicted druggable, against 36% for RNA-only pockets.
  - Ribosomal ligands had higher median MW, HBA/HBD count, rotatable bonds, clogP and Fsp3 than RNA-only ligands.
  - Rekand IH, Brenk R, *J Chem Inf Model* 61:4068 (2021), doi:10.1021/acs.jcim.1c00155 — [PMC8389535](https://europepmc.org/article/PMC/PMC8389535)
- The ChemRxiv v1 preprint gave different figures: 24% druggable, 731/732 sites. Cite the JCIM version — [ChemRxiv v1](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.13303187.v1)

**Ribosome-specific structural portals**
- RiboXYZ (verified earlier). I found no 2024–2026 update paper in Europe PMC or web searches.
  - The GitHub repo is active (last push 2026-05-15).
  - Recent commits concern nascent-peptide exit-tunnel (NPET) extraction, PTC/constriction-site detection and tunnel meshes, not ligand-site features — [github.com/rtviii/riboxyz](https://github.com/rtviii/riboxyz)
- RiboVision2 is a web server for phylogenetic, structural and evolutionary properties of rRNAs at the 1D, 2D and 3D levels.
  - It works for any wwPDB ribosome, using Mol* viewers, an alignment viewer and a 2D viewer.
  - McCann HM, Meade CD, Banerjee B, Penev PI, … Petrov AS, Williams LD, *J Mol Biol* 436(17):168556 (2024), doi:10.1016/j.jmb.2024.168556 — [Europe PMC record](https://europepmc.org/search?query=RiboVision2); [Georgia Tech news](https://sites.gatech.edu/icool/2024/10/17/a-new-bioinformatics-tool-ribovision2); [Benasque 2024 talk](https://www.benasque.org/2024rna/talks_contr/0113_RiboVision2_Petrov.pdf)
  - Its predecessor RiboVision handled about a dozen manually processed structures. Bernier CR et al., *Faraday Discuss* 169:195 (2014), doi:10.1039/c3fd00126a — [RSC](https://pubs.rsc.org/en/content/articlehtml/2014/fd/c3fd00126a)
  - **UNVERIFIED:** whether RiboVision2 has a ligand/antibiotic contact layer. The abstract mentions none, and the PDF failed to download.
- ProteoVision visualises ribosomal **proteins** across phylogeny (1D/2D/3D). Its DESIRE mode uses RiboXYZ annotations and PDBe REST APIs. No ligand aggregation. Penev PI et al., *NAR* 49:W578 (2021), doi:10.1093/nar/gkab351 — [PMC8265156](https://europepmc.org/article/PMC/PMC8265156)

**PDBe-KB and RCSB**
- PDBe-KB aggregated views are protein (UniProt)-centric: per-protein ligand interaction sites and annotations. I found no evidence of RNA or Rfam-based aggregated views — [ELIXIR TeSS event description](https://dev.tess.elixir-europe.org/events/pdbe-kb-aggregated-views-presenting-pdb-data-in-a-different-context-131d7dc8-cb41-4b94-b565-eb902a25b347); [BioExcel webinar](https://bioexcel.eu/events/making-sense-of-pdb-data-with-pdbe-kb/)

**"Antibiotic binding site atlas"**
- I found no dedicated database or web atlas. The nearest is Paternoga H, Crowe-McAuliffe C, … Wilson DN, "Structural conservation of antibiotic interaction with ribosomes", *Nat Struct Mol Biol* 30:1380 (2023), doi:10.1038/s41594-023-01047-y — [PMC10497419](https://europepmc.org/article/PMC/PMC10497419)
  - This is a set of 17 compounds from 6 classes on *E. coli* 70S at 1.6–2.2 Å, with ordered waters. It reports high conservation of binding mode within a scaffold.
  - It also reports that earlier, lower-resolution models of the same drug classes differ in drug pose and surrounding rRNA conformation. This matters for any archive-wide contact statistics.
  - It is one species, new structures, and not an archive survey.

### Inferences
- The claim that "nobody has aggregated ribosome–ligand structures archive-wide" is false. HARIBOSS++ (60.9% rRNA pockets), SPIRAL (543 rRNA binding events), FURNA and the DrugPred_RNA set (497 entries) all contain ribosomal sites at a scale comparable to our 568 entries.
- The novel element must be the coordinate system and the cross-kingdom comparison, not the inventory.
- ChEMBL ribosome targets: not researched within budget. See Gaps.

### Gaps
- HARIBOSS (2022) original handling of ribosomes: the HAL and bioRxiv full texts were blocked. HARIBOSS++ confirms inclusion now.
- R-BIND's explicit ribosome exclusion: the full text was not accessible.
- RiboVision2 ligand features and alignment method (which alignments, whether mitochondria are included): the PDF was not readable.
- ChEMBL ribosome target records (e.g. "30S/50S ribosomal subunit" organism-level targets) and how bioactivity is linked to structures: not checked.
- RCSB ligand-validation and "binding affinity" annotations for ribosome entries: not checked.

## Q2. Have there been systematic, quantitative, archive-wide analyses of ribosome antibiotic binding sites?

### Takeaway
Yes, but these are either structure-local statistics without a shared coordinate frame, or sequence-variation studies built on a handful of representative structures. The quantitative structural surveys are David-Eden 2010 (65 complexes, 4 species), DrugPred_RNA 2021 (613 pockets), SPIRAL 2026 (543 rRNA binding events) and HARIBOSS++ 2026. The two Melnikov-lab papers of 2025 are the most relevant cross-species analyses of drug-binding residues. RNA–ligand ML benchmarks mostly derive from HARIBOSS or the PDB, and I found no clear statement of ribosome exclusion.

### Cited findings

**Foundational survey**
- David-Eden H, Mankin AS, Mandel-Gutfreund Y, "Structural signatures of antibiotic binding sites on the ribosome", *NAR* 38(18):5982–5994 (2010), doi:10.1093/nar/gkq411 — [NAR](https://academic.oup.com/nar/article/38/18/5982/1065889)
  - Data: 65 ribosome–antibiotic complexes (*T. thermophilus* and *E. coli* 30S, *H. marismortui* and *D. radiodurans* 50S). Sites are rRNA atoms within 6 Å, mapped onto apo structures, and compared with 17,610 background pockets from 3V.
  - Sites were enriched for unpaired bases, *syn* conformation, C2′-endo/C4′-exo puckers, Watson–Hoogsteen pairs and non-adjacent stacking, high network centrality and high conservation (ARB-SILVA). These form a "druggability index".
  - Predicted unexplored pockets: the sarcin–ricin loop (DI 0.68) and the A-site finger (DI 0.49).
  - Bacteria and one archaeon only; no eukaryotes or mitochondria.

**Pocket druggability**
- DrugPred_RNA (2021): 497 PDB ribosome entries, 613/546 pockets, 217 ligands, druggability predicted per pocket — [PMC8389535](https://europepmc.org/article/PMC/PMC8389535)

**Interaction-fingerprint clustering**
- SPIRAL (2026): 543 rRNA binding events with DSSR-derived fingerprints. Ribosomal sub-sites cluster by drug class (A-site aminoglycosides R-C1 n = 208; PTC/tunnel macrolides and lincosamides). Affinity predictors were computed across 275 affinity-characterised entries, all RNA categories pooled — [PMC13228370](https://europepmc.org/article/PMC/PMC13228370)

**Cross-species drug-binding residue evolution (Melnikov lab, Newcastle)**
- Ekemezie CL, Chan LI, Brown CR, Helena-Bueno K, Williams TA, Melnikov SV, "Extensive natural variation in bacterial ribosomal drug-binding sites", *Cell Rep* 44(7):115878 (2025), doi:10.1016/j.celrep.2025.115878; preprint doi:10.1101/2024.10.31.621234 — [Cell Reports](https://www.cell.com/cell-reports/fulltext/S2211-1247(25)00649-7); [Europe PMC search](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:PPR933343&format=json&resultType=core)
  - They reconstruct the evolutionary history of drug-binding residues across bacteria.
  - Several rRNA features thought universal to bacteria occur only in a subset. Divergent sites often arise from old rRNA polymorphisms.
  - They show intrinsic resistance in selected species that carry divergent sites.
  - Abstract only: the methods could not be opened (HTTP 429). The design probably mirrors the eukaryotic paper below, with SILVA and *E. coli* numbering (UNVERIFIED).
- Chan LI, Ekemezie CL, Brown CR, Helena-Bueno K, Williams TA, Melnikov SV, "Evolution of drug-binding residues in eukaryotic ribosomes", *Cell Rep* 44(9):116244 (2025), doi:10.1016/j.celrep.2025.116244; preprint doi:10.1101/2024.11.28.625670 — [bioRxiv full text](https://www.biorxiv.org/content/10.1101/2024.11.28.625670v1.full)
  - Methods from the preprint:
    - 29 PDB structures of eukaryotic ribosomes, each with one representative of a drug family (emetine, homoharringtonine, puromycin, blasticidin S, lactimidomycin, cycloheximide, paromomycin, edeine, cryptopleurine and others).
    - Binding residues are heavy atoms within 4.3 Å in ChimeraX, restricted to rRNA bases and r-protein side chains: 58 residues, 26 of them in 18S.
    - *E. coli* numbering throughout, with yeast and human correspondence in Table S2.
    - Sequences: SILVA NR99 v138.1 aligned with MAFFT. Organelle-encoded rRNA was **excluded**, so there is no mitochondrial analysis. Archaea were not analysed.
    - Final set: 8,563 eukaryotes from 5,148 genera. Parasites covered include *Plasmodium*, *Leishmania*, *Trypanosoma* and Microsporidia.
  - Findings:
    - About 23% of eukaryotes carry 1–4 fixed changes relative to human. There are more than 60 unique binding-residue combinations.
    - 15 of 26 18S residues vary; U1489 and A1491 variants occur in nearly 2,000 species each.
    - Only U2393 and G2422 vary in the LSU.
    - Some lineages differ more from human than human differs from bacteria.
  - Data are on Figshare at doi:10.6084/m9.figshare.27918987. No web resource.

**Variation within bacterial rRNA, older work**
- A study of species-specific susceptibility from phylogenetic rRNA variation at drug sites — *Antimicrob Agents Chemother* (2011), doi:10.1128/aac.01398-10 — [ASM](https://journals.asm.org/doi/10.1128/aac.01398-10). Abstract only; authors not verified here.

**ML datasets and benchmarks**
- RNAmigos2 trains on 1,740 RNA–ligand binding sites from the PDB (436 clusters by RMAlign, 264 unique ligands). The full text contains no explicit statement about including or excluding ribosomes. Carvajal-Patiño JG et al., *Nat Commun* 16:2799 (2025), doi:10.1038/s41467-025-57852-0 — [PMC11928640](https://europepmc.org/article/PMC/PMC11928640)
- HARIBOSS-derived benchmarks:
  - SMARTBind: HARIBOSS 10-fold and 5-fold sequence, structure, pocket, pair and ligand splits — [Zenodo](https://zenodo.org/records/17197893)
  - CoBRA: TR60 + HARIBOSS. Jang W, Shin WH, *Brief Bioinform* 27:bbaf713 (2026), doi:10.1093/bib/bbaf713 — [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=CoBRA%20RNA%20language%20model&format=json)
  - GerNA-Bind: splits by date — [Benasque 2026 talk](https://www.benasque.org/2026rna/talks_contr/099_Waldispuhl_Benasque_July2026.pdf)
- The fpocketR evaluation filtered HARIBOSS to remove redundancy and limit RNA size, then added a single ribosome (PDB 7K00). In practice, size filters remove most ribosomes from such benchmarks — [search summary citing fpocketR curation; see PMC11974927](https://pmc.ncbi.nlm.nih.gov/articles/PMC11974927)
- RLBind (Wang K et al., *Brief Bioinform* 24:bbac486, 2023, doi:10.1093/bib/bbac486), RBind (2018) and GATRsite (2025): the abstracts give no dataset composition or ribosome handling — [Europe PMC search](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE:RLBind&format=json)

### Inferences
- No existing survey quantifies, per drug class and across kingdoms, which consensus positions are contacted, how contact patterns vary across copies and species of the same drug, and how they relate to annotated 3D motifs.
- SPIRAL comes closest on fingerprints but stays structure-local. The Melnikov papers come closest on cross-species comparison, but they define binding residues from about 29 representative structures and infer variation from sequence alone.
- ML benchmarks that use HARIBOSS inherit a ribosome-heavy pocket distribution unless size or redundancy filters remove ribosomes (as in the fpocketR curation). Reporting how ribosomal pockets are treated is an open methodological point.

### Gaps
- Exact structure counts and methods of the bacterial Melnikov paper (blocked).
- Whether RLBind's, RNAsite's (TR60/TE18) or RNAmigos2's datasets explicitly exclude ribosomes: not stated in the abstracts or the RNAmigos2 full text.
- A review-level synthesis such as Lin J, Zhou D, Steitz TA, Polikanov YS, Gagnon MG, *Annu Rev Biochem* 2018 ("Ribosome-targeting antibiotics: modes of action, mechanisms of resistance, and implications for drug design") was not checked for quantitative content — [ResearchGate listing](https://www.researchgate.net/publication/323976028_Ribosome-Targeting_Antibiotics_Modes_of_Action_Mechanisms_of_Resistance_and_Implications_for_Drug_Design)

## Q3. Has anyone mapped ribosome ligand contacts onto a common alignment or consensus numbering across bacteria, archaea, eukaryotes and mitochondria?

### Takeaway
Only partially, and never archive-wide from structures. The prevailing practice is *E. coli* numbering via sequence alignment of a few representative binding sites: Melnikov lab (SILVA + MAFFT; bacteria and eukaryotes separately; organelle rRNA excluded) and David-Eden 2010 (conservation via ARB-SILVA). RiboXYZ transfers 5 Å sites pairwise between two chains with pairwise2.globalxx. I found no resource that places every structure's contacts on Rfam CM consensus positions, or on any shared cross-kingdom coordinate system that includes mitochondrial rRNA.

### Cited findings
- Melnikov eukaryotic paper: *E. coli* numbering with yeast and human correspondence tables. MAFFT alignment of SILVA NR99 sequences. Organelle-encoded rRNA excluded. Binding residues come from 29 structures, not from all PDB copies — [bioRxiv 2024.11.28.625670](https://www.biorxiv.org/content/10.1101/2024.11.28.625670v1.full)
- Melnikov bacterial paper: bacteria only — [Cell Reports](https://www.cell.com/cell-reports/fulltext/S2211-1247(25)00649-7) (abstract)
- HARIBOSS++ uses Rfam only as a chain-level family cross-check, not for positions — [HARIBOSS++](https://www.biorxiv.org/content/10.64898/2026.09.07.749797v1.full)
- FURNA assigns Rfam families per chain (Rfam–PDB mapping or cmscan) and builds Infernal CMs for its sequence search. Ligand-binding residues stay in chain numbering — [PMC11309384](https://europepmc.org/article/PMC/PMC11309384)
- SPIRAL cross-checks functional categories against Rfam only — [PMC13228370](https://europepmc.org/article/PMC/PMC13228370)
- Mitochondrial-specific structural work is per-study, not systematic. Examples:
  - Tigecycline occupies three sites on the T-cell mitoribosome, including an A-site and a PTC site that align with bacterial helices — [PMC12045974](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12045974/)
  - Mitoribosome profiling of antibiotic off-target activity — [eLife reviewed preprint 103699](https://elifesciences.org/reviewed-preprints/103699v1)
  - Older biochemical evidence of aminoglycoside binding to human mt-rRNA — [JBC](https://www.jbc.org/article/S0021-9258(20)42296-3/fulltext)
- RiboVision2 shows rRNA alignments and phylogenetic properties on 2D/3D for any PDB ribosome. I could not confirm whether it maps ligand contacts or covers mitochondrial rRNA — [GT news](https://sites.gatech.edu/icool/2024/10/17/a-new-bioinformatics-tool-ribovision2)

### Inferences
- A structure-derived, copy-level mapping of all PDB ribosome–drug contacts onto Rfam CM consensus positions appears to be unaddressed. It would cover the bacterial SSU/LSU families (RF00177/RF02541) and the eukaryotic families (RF01960/RF02543), with mitochondrial and parasite rRNAs mapped through whichever family their CM hits.
- So does linking those positions to RNA 3D motifs (BGSU loops) across kingdoms. The Melnikov work is the main prior art on cross-lineage residue variation and should be cited and compared directly. Their 58-residue eukaryotic set and *E. coli* numbering could serve as a benchmark for our consensus-position mapping.
- One caveat for our design: Rfam SSU/LSU CMs are domain-specific. Bacterial-to-eukaryotic correspondence will still need a bridge (structure superposition or a cross-family alignment). Mitochondrial rRNAs may score poorly against all four CMs. This is our own inference, not checked against Rfam documentation here.

### Gaps
- Whether RiboVision2 or the Petrov/Williams "common core" rRNA alignments (bacteria, archaea, eukaryotes) carry ligand annotations: unverified.
- Whether any 2025–2026 work (e.g. from the RiboXYZ/Dao Duc group) added cross-species ligand-site transfer beyond pairwise2: none found. The repo activity is tunnel-focused.

## Q4. What is clearly missing or unaddressed?

### Takeaway
No current resource or survey combines all of the following, which the planned study proposes:
- every PDB copy of each ribosome-bound drug, with explicit denominators per copy and per entry;
- contacts on a model-based, shared rRNA coordinate system (Rfam CM consensus positions) instead of per-chain numbering or *E. coli* numbering via ad hoc alignment;
- explicit cross-kingdom coverage including mitochondrial and parasite ribosomes from structures;
- linkage of contacted positions to annotated RNA 3D motifs;
- a quantitative comparison of contact patterns across drug classes and kingdoms.

### Cited findings
- Existing databases keep chain-local numbering and use Rfam only as a family label: HARIBOSS++ ([preprint](https://www.biorxiv.org/content/10.64898/2026.09.07.749797v1.full)), FURNA ([PMC11309384](https://europepmc.org/article/PMC/PMC11309384)), SPIRAL ([PMC13228370](https://europepmc.org/article/PMC/PMC13228370)).
- Cross-species residue studies use few representative structures and sequence alignments, and exclude organellar rRNA — [Chan et al. preprint](https://www.biorxiv.org/content/10.1101/2024.11.28.625670v1.full)
- Lower-resolution structures of the same drug class disagree in pose and in local rRNA conformation. Aggregation across copies needs quality-aware handling, which no resource provides — [Paternoga et al. 2023](https://europepmc.org/article/PMC/PMC10497419)
- RiboXYZ transfers sites only pairwise, through global sequence alignment (verified earlier; [NAR 2023](https://academic.oup.com/nar/article/51/D1/D509/6777803)).
- The motif context in SPIRAL is coarse DSSR categories (hairpin, junction, pseudoknot, multiplet and others), not curated recurrent 3D motif families such as BGSU/RNA 3D Motif Atlas — [PMC13228370](https://europepmc.org/article/PMC/PMC13228370)

### Inferences
- Novelty is strongest for:
  1. the consensus-position coordinate system applied to all copies;
  2. structure-based (not sequence-only) cross-kingdom comparison including mitoribosomes;
  3. per-position occupancy statistics with explicit denominators (structures where a position is observed and eligible);
  4. linking to Motif Atlas loops.
- Novelty is weak for:
  - the claim that the inventory is large: HARIBOSS++, SPIRAL and FURNA already hold hundreds of ribosomal sites;
  - interaction-fingerprint clustering by drug class (SPIRAL did this in 2026);
  - drug-site sequence variation across bacteria and eukaryotes (Melnikov 2025).
- The planned study should cite and position against SPIRAL, HARIBOSS++ and both Melnikov papers explicitly.

### Gaps
- No systematic check of ChEMBL ribosome targets, PDBbind/BindingDB coverage of ribosome affinities, or the RCSB ligand-validation reports for ribosome entries.
- Archaeal ribosome–drug structures (e.g. *H. marismortui*) are covered only in David-Eden 2010. No modern archaeal survey was found, though absence of evidence here reflects a limited search.
