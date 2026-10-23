# RNA-specific resources for RNA–small-molecule interaction annotation (as of 2026-10-23)

Scope: RNA-focused tools, services, databases and benchmark sets that annotate or precompute RNA–ligand contacts, as alternatives or cross-checks to PDBe-Arpeggio (PDBe `bound_ligand_interactions`). About 17 search/fetch calls were made. Several resources in the brief could not be verified in that budget; they are listed under Gaps and not described from memory. Unreachable sites are flagged.

## Q1. Which resources give atom-level, machine-readable RNA–ligand interaction typing (not just residue lists)?

### Takeaway
Only general interaction profilers extended to nucleic acids give atom-level typed contacts that can be run in bulk: fingeRNAt (Apache-2.0, nine interaction classes, built for nucleic-acid receptors) and PLIP 2021 (open source, nucleic acids supported since 2021). The curated RNA–ligand databases (HARIBOSS, R-BIND, ROBIN, NALDB, RNALigands, BioLiP2) give complexes, pockets, residue lists or activity labels, not typed atom-level contacts. DSSR reports atom-level H-bonds, but by default only within nucleic acids. Whether its free Basic licence covers ligand H-bonds was not confirmed.

### Cited Findings
**fingeRNAt (tool, recommended cross-check)**
- Python tool that detects non-covalent nucleic acid–ligand interactions and encodes them as Structural Interaction Fingerprints (SIFt). Published in PLoS Comput Biol (June 2022). Code at https://github.com/n-szulc/fingeRNAt — [PMC9197077](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9197077/); [DOAJ record](https://doaj.org/article/6cbd49b2c37344c7a5bbee1e99d33373)
- It detects nine interaction classes: hydrogen and halogen bonds, cation–anion, π–cation, π–anion, π-stacking, inorganic-ion-mediated, water-mediated and lipophilic contacts. It has a plugin system to add more, a PyMOL plugin, and a fingerprint-similarity alternative to RMSD. It was applied to cluster RNA–ligand complexes from the PDB — [PMC9197077](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9197077/)
- Licence: Apache-2.0. Zenodo releases exist (v0-1-1, December 2022). The search found no evidence of releases after 2022, so treat it as possibly unmaintained — [Zenodo 7486222](https://zenodo.org/records/7486222/latest); [science.ecosyste.ms](https://science.ecosyste.ms/projects/56315)

**PLIP 2021 (tool and web server)**
- Adasme et al., NAR 2021, 49(W1):W530–W534, doi 10.1093/nar/gkab294. "Expanding the scope … to DNA and RNA". It detects hydrogen bonds, hydrophobic contacts, π-stacking, π-cation, salt bridges, water bridges, metal complexes and halogen bonds, with adjustable thresholds. Open source on GitHub, plus a web server — [PMC8262720](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8262720/); [NAR](https://academic.oup.com/nar/article/49/W1/W530/6282136)
- The paper's demonstrations are a drug bound to DNA and the FUS zinc finger bound to RNA, so a protein receptor with an RNA ligand. Confirm RNA-as-receptor with a small-molecule ligand separately — [PMC8262720](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8262720/)

**x3dna-DSSR**
- Geometry-based H-bond detection (`--get-hbonds`) for nucleic acids and proteins. By default it reports only H-bonds within nucleic acids. Output can be JSON, with atom IDs, distance and a residue-pair type such as "nt:aa" — [x3dna.org: Get hydrogen bonds with DSSR](https://x3dna.org/highlights/get-hydrogen-bonds-with-dssr); [H-bonds reported by DSSR and SNAP](https://x3dna.org/articles/h-bonds-reported-by-dssr-and-snap)
- The most recent release mentioned is v2.7.2-2026jan12, which adds `--pair-wise`. The Columbia Technology Ventures (CTV) page lists a v2.8.1 download — [x3dna.org posts](https://x3dna.org/posts-list/); [CTV DSSR licences](https://inventions.techventures.columbia.edu/technologies/dssr-an-integrated-software--CU20391/licenses)
- Licensing: CTV has handled DSSR licensing since late 2020. Under NIH R24GM153869, DSSR Basic is free for academics via CTV Express Licensing. Basic is "as-is", has no support, and has fewer analysis and annotation features than the paid DSSR Pro. Commercial use needs a CTV licence — [x3dna.org: DSSR is funded and DSSR Basic academic is free](https://x3dna.org/highlights/x3dna-dssr-is-funded-and-dssr-basic-academic-is-free); [CTV Basic Academic licence](https://inventions.techventures.columbia.edu/licenses/dssr-basic-academic-1-seat-now--188); [DSSR-Pro Commercial](https://inventions.techventures.columbia.edu/licenses/dssr-pro-commercial--215)
- DSSR blog case studies analyse ligand complexes, for example the class V GTP aptamer–GTP complex — [x3dna.org](https://x3dna.org/highlights/dssr-analysis-of-the-class-v-gtp-aptamer-gtp-complex-structure)

**RNApolis annotator**
- Its documented scope is base pairs (Leontis–Westhof and Saenger), stacking, base–ribose and base–phosphate contacts. Output formats are dot-bracket, BPSEQ, CSV, JSON and DOT. An adapter converts FR3D, DSSR and RNAView output. No ligand handling is documented — [PyPI RNApolis](https://pypi.org/project/RNApolis)
- An RNApdbee 3.0 talk (Benasque 2026) lists base pairs, triples and pseudoknots, with no ligands — [Benasque 2026 slides](https://www.benasque.org/2026rna/talks_contr/088_RNApdbee3Benasque2026.pdf)

### Inferences
- fingeRNAt is the closest RNA-native analogue to Arpeggio or PLIP. It is the best independent, open-licence cross-check for atom-level typed contacts, including ion- and water-mediated ones. PLIP 2021 is a second independent open tool.
- DSSR suits an RNA-centric H-bond listing, but its licence (academic-only Basic, commercial Pro) limits redistribution in a public EMBL-EBI pipeline. Its ligand H-bond behaviour still needs checking.
- RNApolis, MC-Annotate and FR3D-type annotators are RNA–RNA tools. Their value here is the base-edge vocabulary, not ligand contacts. That is an inference: no source showed ligand support.

### Gaps
- DSSR: no source found confirming whether ligand (HETATM) H-bonds are reported, or under which option and licence tier. Needs the DSSR manual or forum.
- fingeRNAt: no maintenance status for 2023–2026 was found. Its exact distance and angle defaults were not extracted.
- MC-Annotate, ClaRNA, NucleicNet, RNAProSite and rsRNASP were not researched in this budget. rsRNASP is, by its name, a statistical potential, not an annotator, but this is not verified.

## Q2. Which resources explicitly annotate base edges, stacking, and ion- or water-mediated contacts?

### Takeaway
Ion- and water-mediated contacts plus π-stacking are typed by fingeRNAt and PLIP. No tool or database found annotates ligand contacts by nucleobase edge (Watson–Crick, Hoogsteen, sugar). Edge-level classification of base:ligand pairs exists only as a published, curated analysis of riboswitch crystal structures (RNA 2019).

### Cited Findings
- fingeRNAt has explicit "inorganic ion-mediated" and "water-mediated" classes, plus π-stacking — [PMC9197077](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9197077/)
- PLIP 2021 types water bridges, metal complexes and π-stacking — [PMC8262720](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8262720/)
- RNA journal 2019, vol. 25(10), p. 1274: "Pairing interactions between nucleobases and ligands in aptamer:ligand complexes of riboswitches". The study analysed 67 riboswitch aptamer crystal structures and found 80 base:ligand pairs with at least two H-bonds.
  - Contacts are classified by the interacting nucleobase edge. Watson–Crick-edge contacts are most common, followed by purine sugar-edge and uracil Hoogsteen-edge contacts.
  - Pairs are also classed by ligand chemistry: base:base, base:amino acid, base:sugar, base:phosphate and base:other.
  - In QM checks, 67 of the 80 pairs were stable near their crystal geometry.
  - Sources: [RNA 25:1274](https://rnajournal.cshlp.org/content/25/10/1274.full); [PMC6800475](https://pmc.ncbi.nlm.nih.gov/articles/PMC6800475); [bioRxiv 603159](https://www.biorxiv.org/content/10.1101/603159.full.pdf)
- HARIBOSS classifies ions and water separately from ligands and excludes them from the ligand set. Modified residues are assigned to RNA, DNA or protein using mmCIF covalent-bond information. It does not type interactions systematically: H-bonds, salt bridges and π-stacking appear only as figure examples — [HARIBOSS bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full)
- The RNApolis annotator gives base–phosphate and base–ribose classes for RNA–RNA contacts only — [PyPI RNApolis](https://pypi.org/project/RNApolis)

### Inferences
- An edge-level ligand-contact annotation would have to be derived: for example, map typed atom contacts from fingeRNAt, PLIP or Arpeggio onto base-edge atom sets (WC, Hoogsteen, sugar edge, ribose, phosphate). The RNA 2019 dataset of 80 pairs can validate that derivation for riboswitches.

### Gaps
- No machine-readable supplementary table for the RNA 2019 base:ligand pair set was confirmed (PDB IDs, residues, edges). Check the paper's supplement.
- Whether BGSU RNA 3D Hub or FR3D annotate any RNA–ligand interaction was not verified in this session.

## Q3. Which could provide an independent "reference" annotation for benchmarking a contact tool on RNA?

### Takeaway
No expert-curated, atom-level ground truth for RNA–ligand contacts was found. The usable references are:
- (a) HARIBOSS, a curated complex and pocket set with explicit selection rules;
- (b) binding-site benchmark sets (RNAsite TE18 and others, as used by RLBind and CoBRA), which label nucleotides by distance cutoffs that differ between papers (4 Å vs 6 Å);
- (c) the RNA 2019 riboswitch base:ligand pair catalogue, which is H-bond and edge-level and literature-derived;
- (d) BioLiP2 for nucleic-acid receptors, whose server was unreachable during this research.

Activity databases (R-BIND, ROBIN, NALDB, SMMRNA, Inforna) give binders and non-binders, not 3D contacts.

### Cited Findings
**HARIBOSS**
- 2,131 RNA–small-molecule complexes, 550 compounds and 5,002 pockets, from X-ray, NMR and cryo-EM. No last-update date, download format or licence is shown on the page (fetched 2026-10-23) — [hariboss.pasteur.cloud](https://hariboss.pasteur.cloud/)
- Panei, Torchet, Ménager, Gkeka, Bonomi, bioRxiv doi 10.1101/2022.05.17.492306. A HAL record suggests final publication in *Bioinformatics* (supplement prefix "btac"), not JCIM; verify the venue. The Institut Pasteur page says it is "regularly updated" — [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full.pdf); [HAL](https://hal-pasteur.archives-ouvertes.fr/hal-03834209v2); [Pasteur](https://research.pasteur.fr/en/b/qJN)
- Ligand filter: non-polymer, 160–1000 Da, at least one C, and only the elements C, H, N, O, Br, Cl, F, P, Si, B and Se.
  - Interacting chains are those within 6 Å of any ligand atom; an RNA chain needs at least 10 atoms within that radius.
  - The redundant set keeps only structures where the ligand contacts RNA exclusively.
  - Pockets come from mkgridXf with watershed sub-cavities; a sub-cavity is "occupied" if ≥20% of ligand atoms fall inside it. SiteMap is also used.
  - No RSCC or model-quality filter. The non-redundant representative is chosen by best resolution.
  - Source: [HARIBOSS bioRxiv full text](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full)
- A downstream benchmark (arXiv 2310.18985) used a HARIBOSS snapshot of 746 complexes and 1,268 sites, with the site cutoff at 6 Å — [arXiv 2310.18985](https://arxiv.org/pdf/2310.18985)

**Binding-site benchmark sets (RNAsite, RLBind, CoBRA, rnaglib)**
- RNAsite (Su, Peng, Yang, *Bioinformatics* 2021, 37:36–42) distributes its training and test sets as a single archive, `training_and_test.tar.br2`. The page gives no composition or cutoff — [RNAsite benchmark](https://yanglab.qd.sdu.edu.cn/RNAsite/benchmark); [RNAsite](https://yanglab.qd.sdu.edu.cn/RNAsite/)
- The rnaglib benchmark (arXiv 2503.21681) takes its datasets and splits from RNASite, which RLBind also uses. It labels residues as binding at < 6 Å — [arXiv 2503.21681](https://arxiv.org/html/2503.21681v3)
- The CoBRA preprint (bioRxiv 2025) evaluates on TE18, RB9, JL10 and TL12. It labels a nucleotide as binding if any atom is ≤ 4 Å from a ligand atom — [CoBRA bioRxiv](https://www.biorxiv.org/content/10.1101/2025.09.03.673903.full.pdf)
- An RLBind-related GAT and language-model paper is at [PMC12417085](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12417085/). MultiModRLBP uses "Test18" and "Test3" (likely TE18; not confirmed) — [PMC12417085](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12417085/)

**rnaglib**
- Python package with Zenodo-hosted, versioned datasets of annotated RNA 3D structures. You can choose all PDB RNAs at build time or the BGSU representative set. It includes "small molecule binding partners" with molecule-, base-pair- and nucleotide-level attributes. The code and data licences were not stated in search results — [GitHub cgoliver/rnaglib](https://github.com/cgoliver/rnaglib); [PyPI](https://pypi.org/project/rnaglib); [arXiv 2402.09330](https://arxiv.org/pdf/2402.09330)

**Riboswitch literature ground truth**
- 80 base:ligand pairs across 67 riboswitch aptamer structures, edge-classified — [RNA 25:1274](https://rnajournal.cshlp.org/content/25/10/1274.full)
- An older NMR study of purine riboswitches reports a conserved base pair between a uridine and the purine ligand's N3/N9 edge — [PMC547832](https://pmc.ncbi.nlm.nih.gov/articles/PMC547832)

**BioLiP2 (includes DNA and RNA ligand categories; protein receptors primarily)**
- NAR 2024 52(D1):D404, doi 10.1093/nar/gkad630.
  - Free downloads: `BioLiP.txt.gz`, `BioLiP_nr.txt.gz`, and FASTA files per receptor and ligand type.
  - 24,809 PDBbind-CN affinity records were removed in January 2025 over licensing.
  - The download readme was unreachable (ECONNREFUSED) on 2026-10-23.
  - Sources: [download page](https://seq2fun.dcmb.med.umich.edu/BioLiP2/download.html); [NAR](https://academic.oup.com/nar/article/52/D1/D404/7233921); [readme](https://seq2fun.dcmb.med.umich.edu/BioLiP2/download/readme.txt)
- Q-BioLiP (Yang lab, Shandong) is a related quaternary-structure resource with weekly updates. It explicitly covers "RNA–small molecule" interactions and offers mmCIF and PDB downloads — [Q-BioLiP download](https://yanglab.qd.sdu.edu.cn/BioLiP2/Download/)

**Activity and ligand databases (no 3D contacts)**
- R-BIND 2.0: ACS Chem Biol 2022, 17:1556–1566. The small-molecule count grew from 97 to 153. Version 2.1 lists 175 small and 48 large molecules. Site: rbind.chem.duke.edu — [PMC9343015](https://pmc.ncbi.nlm.nih.gov/articles/PMC9343015); [R-BIND updates](https://rbind.chem.duke.edu/updates.html)
- ROBIN: 36 small-molecule-microarray screens of nucleic acids against 24,572 compounds, giving 2,003 RNA binders; machine-learning models were trained on them — [PMC9992102](https://pmc.ncbi.nlm.nih.gov/articles/PMC9992102); [bioRxiv 2022.08.01.502065](https://www.biorxiv.org/content/10.1101/2022.08.01.502065v1)
- NALDB (2016) catalogues nucleic-acid-targeting small molecules. URL given in the paper: http://bsbe.iiti.ac.in/bsbe/naldb/HOME.php, current status unverified — [PMC4761111](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4761111/)
- SMMRNA: about 1,000 literature-curated RNA–small-molecule interactions. No working URL was found — [IITM related resources](https://web.iitm.ac.in/bioinfo2/R_SIM/Related_resources.html)
- RNALigands: RNA 2022, 28:115. Server at http://rnaligands.ccbr.utoronto.ca, code at https://github.com/SaisaiSun/RNALigands, current status unverified — [RNA 28:115](https://rnajournal.cshlp.org/content/28/2/115)

### Inferences
- For a contact-tool benchmark on RNA, the most defensible independent references are:
  1. a residue-level cross-check against HARIBOSS or RNAsite labels, with the 4 Å vs 6 Å cutoff stated and harmonised;
  2. an atom- and type-level cross-check against fingeRNAt and/or PLIP run locally on the same structures;
  3. an edge-level, literature-validated subset from the RNA 2019 riboswitch base:ligand pairs (purine, SAM, TPP, FMN, lysine, glmS, preQ1 are likely among its 67 structures, not verified per family).
- None of the curated 3D sets apply ligand-quality filters (RSCC or validation) as far as could be verified (HARIBOSS explicitly does not). Ligand-quality filtering would have to come from PDBe or wwPDB validation data.

### Gaps
- Per-family, literature-curated contact sets for purine, SAM, TPP, FMN, lysine, glmS and preQ1 riboswitches were not located as standalone datasets.
- Docking benchmark sets (rDock RNA set, RLDOCK, NLDock, AnnapuRNA test sets) were not researched in this budget.
- Exact sizes and complex lists of TE18, RB9, JL10 and TL12 were not extracted.

## Q4. Bulk access and licences suitable for reuse in a public EMBL-EBI context

### Takeaway
The cleanest open options are fingeRNAt (Apache-2.0) and PLIP (open source on GitHub; exact licence not verified here), run in the pipeline; BioLiP2 flat files (free download); and RNAsite and rnaglib datasets (downloadable; licence unstated). DSSR is unsuitable for open redistribution: Basic is free for academics only and Pro is commercial, both via Columbia Technology Ventures. HARIBOSS, NAKB and R-BIND show no explicit data licence on their landing pages.

### Cited Findings
- fingeRNAt: Apache-2.0 — [Zenodo](https://zenodo.org/records/7486222/latest); [GitHub](https://github.com/n-szulc/fingeRNAt)
- DSSR: free Basic academic licence, Pro commercial, via CTV — [x3dna.org](https://x3dna.org/highlights/x3dna-dssr-is-funded-and-dssr-basic-academic-is-free); [CTV](https://inventions.techventures.columbia.edu/licenses/dssr-basic-academic-1-seat-now--188)
- BioLiP2: full data free to download, but PDBbind-derived affinities were removed in January 2025 over licensing — [BioLiP2 download](https://seq2fun.dcmb.med.umich.edu/BioLiP2/download.html)
- NAKB: 21,987 nucleic-acid-containing 3D structures (homepage dated 2026-06-24). It advertises "NA, Protein, and Ligand Annotations" (annotations.html), a Modified Nucleotide Browser and a Download page. No licence beyond "© 2022-2026 NAKB" and NIH NIGMS support is shown.
  - Paper: Lawson et al., NAR 2024 52:D245–D254, doi 10.1093/nar/gkad957. It links Rfam, RNAcentral, DNATCO, 3DFootprint and DNAproDB.
  - A Hugging Face mirror exists (RosettaCommons/NAKB).
  - Sources: [nakb.org](https://www.nakb.org/); [RCSB news](https://bioinformatics.rcsb.org/news/645e4abc1d621d75127a7a49); [HF NAKB](https://huggingface.co/datasets/RosettaCommons/NAKB/blob/main/README.md)
- HARIBOSS: the website allows downloading query results; formats and licence are not stated — [hariboss.pasteur.cloud](https://hariboss.pasteur.cloud/); [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full)

### Inferences
- For an EMBL-EBI public context, the pattern that avoids licence risk is to run open tools (fingeRNAt, PLIP) locally on wwPDB data (CC0) and use HARIBOSS or RNAsite only as external comparison labels, cited rather than redistributed, until their licences are confirmed.

### Gaps
- NAKB's ligand-annotation content (annotations.html, download.html) and API could not be read: the PMC full text was blocked by a browser check and the homepage gives no detail.
- PLIP's exact licence (believed to be GPL-2.0; not verified here) and rnaglib's code and data licences.
- Not researched in this budget, with no claims made: Inforna, PDBsum/LigPlot for nucleic acids, PDBe-KB ligand pages, RCSB ligand interaction view, BGSU RNA 3D Hub ligand annotation, MC-Annotate, ClaRNA, NucleicNet, RNAProSite, rsRNASP, rDock RNA set, RLDOCK, NLDock, AnnapuRNA, HARIBOSS++ or any HARIBOSS update. Live status of NALDB, SMMRNA and RNALigands was not checked.
