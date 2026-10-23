# Comparing, aligning and clustering RNA small-molecule binding sites and binding modes

Scope: (a) grouping ligand copies into binding modes within a site, (b) deciding whether two binding sites in different RNA families are equivalent, (c) comparing pockets structurally. Project context: sites are sets of Rfam CM match positions; ligand copies are grouped by ECIFP Tanimoto; clan siblings are mapped to the main family's CM and pockets checked by backbone RMSD of the mapped pocket positions. Research date 2026-10-23. Licences were checked on the GitHub API on that date.

## Q1. RNA structural alignment tools usable for pocket/site correspondence, and their accuracy for local regions

### Takeaway
Global RNA 3D aligners (RNA-align/US-align, RMalign, ARTEMIS, STAR3D, SARA, ARTS, Rclick) are benchmarked on whole chains with TM-score, RMSD or structure overlap. None is benchmarked on ligand-pocket correspondence. For local regions, the most relevant are R3D-Align (local neighbourhood superpositions merged by maximum clique), LocalSTAR3D (multiple local alignments), and ARTEMIS/ARTEM, which are topology-independent and seeded by near-zero-RMSD residue matches. The project's approach (CM mapping, then RMSD of mapped pocket positions) is a sequence/CM-anchored local superposition. A structure aligner is the independent check on it, and the fallback for sites the CM does not align.

### Cited Findings
- **US-align** (Nature Methods 2022): universal aligner for proteins, RNA and DNA. On 637 monomeric RNA chains its TM-scores were 5.8% higher than RMalign, 27.5% higher than STAR3D, 34.5% higher than ARTS and 38.6% higher than Rclick — [US-align, Nat Methods 2022](https://www.nature.com/articles/s41592-022-01585-1); [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.04.18.488565v1.full). Code at [github.com/pylelab/USalign](https://github.com/pylelab/USalign) (GitHub reports the licence as "NOASSERTION", i.e. a custom free licence; active, last push Oct 2026).
- **US-align2** (iScience 2022) adds sequential and non-sequential (fNS, sNS) modes. For detecting same-Rfam-family RNA pairs, fNS and sNS had about 1.2% higher AUROC than the sequential mode, which in turn beat RMalign by 1.4% and ARTS by 12.5%. Caveat: 68.5% of same-Rfam pairs had ≥40% sequence identity, which compresses the differences — [Zhang et al. 2022, iScience](https://www.sciencedirect.com/science/article/pii/S2589004222014900).
- **RNA-align** (Bioinformatics 2019): open source, with a size-independent TM-score_RNA. It "significantly outperforms other state-of-the-art programs in both alignment accuracy and running speed" — [Gong et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6821192).
- **RMalign** (BMC Genomics 2019): size-independent RMscore, built to replace SARA's size-dependent score in the PRIME template pipeline. tRNA precision 0.88 vs 0.78 (SARA). Functional classification AUC 0.95, comparable to ESA-RNA (CLICK 0.67, SETTER 0.75) — [Zheng et al. 2019](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5631-3).
- **ARTEMIS** (NAR 2024, Bohdan, Bujnicki, Baulin): topology-independent superposition plus structure-based sequence alignment, developed from ARTEM (pairwise superposition of RNA 3D modules). It "significantly outperforms state-of-the-art tools in both sequentially-ordered and topology-independent RNA 3D structure superposition". The heuristic finds the optimal solution provided it includes at least one residue–residue match with near-zero RMSD. Using it, the authors found a helical packing motif preserved across different backbone topologies in ribozymes and riboswitches — [ARTEMIS, NAR 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11472068/). Code: [github.com/david-bogdan/ARTEMIS](https://github.com/david-bogdan/ARTEMIS) (Apache-2.0, active 2026).
- **R3D-Align** (Bioinformatics 2010; NAR web server 2013): nucleotide-to-nucleotide alignment of homologous structures from many local neighbourhood superpositions, merged by a maximum-clique algorithm. Neighbourhoods include non-WC and WC pairs plus tertiary pairs and stacking. On 5S/16S/23S rRNA it was the most accurate against hand alignments. The authors note that ARTS, which uses a single rigid transformation, leaves whole regions unaligned. Web server rna.bgsu.edu/r3dalign; MATLAB code on GitHub — [Rahrig et al. 2010](https://academic.oup.com/bioinformatics/article/26/21/2689/213314); [BGSU page](https://www.bgsu.edu/research/rna/web-applications/r3d-align.html); [web server paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3465099).
- **STAR3D / LocalSTAR3D / CircularSTAR3D**: anchor on conserved stacks. LocalSTAR3D (NAR 2020) reports multiple local alignments between two RNAs, filling the stated gap that there was "still a lack of local alignment tools". CircularSTAR3D (NAR 2023) finds matches under circular permutation — [LocalSTAR3D](https://pmc.ncbi.nlm.nih.gov/articles/PMC7367197); [CircularSTAR3D](https://doi.org/10.1093/nar/gkad222).
- **SARA** minimises RMSD between unit-vector representations. **Rclick** uses clique matching plus least-squares superposition. On 55 difficult pairs Rclick's mean structure overlap was 51.4%, against 32.0% for SARA and 33.4% for ARTS (self-reported) — [Rclick benchmark](https://mspc.bii.a-star.edu.sg/minhn/benchmark_rclick.html); [overview](https://pmc.ncbi.nlm.nih.gov/articles/PMC4787758).
- Shared limitation: aligner scores rely on base-pairing scores and RMSD, which depend on aligned length — [PMC4787758](https://pmc.ncbi.nlm.nih.gov/articles/PMC4787758).
- **Zurkowski, Antczak, Szachniuk 2023** (Bioinformatics, btad315): customisable heuristics for RNA 3D structure alignment — [PMC10199316](https://pmc.ncbi.nlm.nih.gov/articles/PMC10199316/).

### Inferences
- Within a clan, the CM mapping already supplies residue correspondence. A structure aligner adds value as an independent check on whether the mapped pocket positions superpose, and for finding equivalences the CM cannot express (permuted or rearranged scaffolds, e.g. SAM-I vs SAM-IV, below). ARTEMIS (topology-independent, Apache-2.0, Python) and US-align with the non-sequential option are the practical candidates.
- For a local question ("do these 10–20 pocket nucleotides correspond?"), anchor-based local superposition is closer to R3D-Align/ARTEM logic: superpose on the pocket plus a shell of neighbours, then report per-position deviation. This is less sensitive to hinge motions elsewhere than a global TM-score. It matches the base-anchored local superposition the Interactions tab already uses.
- Reporting both the CM-mapped RMSD and an independent aligner's correspondence (agreement fraction over pocket positions) would make "equivalent site" claims auditable.

### Gaps
- No published benchmark of these aligners on small-molecule pocket correspondence specifically, and none on accuracy over sub-regions of about 10–30 nt.
- LaRA (graph-based, sequence–structure alignment) and CLICK/ARTS licences and current availability not verified. LaRA is mainly a 2D/sequence-structure aligner rather than a 3D one (not verified in this session).
- US-align's exact licence terms not read; GitHub reports it as non-standard.

## Q2. RNA pocket comparison/similarity methods, and protein methods applied to RNA

### Takeaway
No established, benchmarked RNA-specific pocket alignment or similarity tool (an RNA analogue of SiteAlign, ProBiS, G-LoSA or APoc) was found. RNA work falls into four groups:
- pocket detection with descriptors: fpocketR, RPocket, HARIBOSS annotations;
- interaction-fingerprint similarity: fingeRNAt SIFt, ProLIF;
- graph learning that embeds base-pair-network pockets: RNAmigos/RNAmigos2;
- global shape descriptors: a 2026 inertia-eigenvalue study.

Protein site-comparison tools are mature but benchmarks warn of inflated performance.

### Cited Findings
- **fpocketR** (PNAS 2025, Weeks lab): detects pockets capable of binding drug-like ligands in RNA structures. "Remarkably few, roughly 50" such pockets have been visualised. Frag-MaP validates them in cells — [Veenbaas et al., PNAS 2025](https://doi.org/10.1073/pnas.2422346122) (PMC12054788). A companion preprint describes fpocketR as a platform for pocket identification, characterisation and visualisation, best at ≤3.5 Å resolution — [bioRxiv 2025.03.25.645323](https://www.biorxiv.org/content/10.1101/2025.03.25.645323v2.full). On a ~160-structure library (many binding the same ligands, e.g. TPP, SAM) it found 139 known pockets and 237 novel ones — [via HARIBOSS-related summary, PMC11974927](https://pmc.ncbi.nlm.nih.gov/articles/PMC11974927). Code: [github.com/Weeks-UNC/fpocketR](https://github.com/Weeks-UNC/fpocketR) (MIT).
- Earlier RNA pocket detection reused protein tools (V3 and mkgrid, both rolling-probe; PocketFinder, energy-based). These studies suggest RNA and protein pockets share bulk volume, buriedness and solvent accessibility — [fpocketR preprint](https://www.biorxiv.org/content/10.1101/2025.03.13.643147.full.pdf).
- **HARIBOSS** (Panei, Torchet, Menager, Gkeka, Bonomi, 2022): curated X-ray/NMR/cryo-EM RNA–small-molecule complexes annotated with physico-chemical properties of ligands and pockets; web interface hariboss.pasteur.cloud — [bioRxiv 2022.05.17.492306](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full.pdf).
- **RPocket** (BMC Bioinformatics 2021): 240 RNA pockets from 94 non-redundant complexes, with geometric size, centroid, shape and secondary-structure-element descriptors (RPDescriptor) — [PMC8424408](https://doi.org/10.1186/s12859-021-04349-4).
- **Shape descriptors** (2026): inertia-matrix eigenvalue ratios of pocket atoms classify RNA and protein pockets into sphere-, rod-, disk-like and strongly anisotropic archetypes, showing pocket shapes shared between RNA and protein — [CSBJ 2026](https://spj.science.org/doi/10.34133/csbj.0022).
- **RNAmigos** (NAR 2020): GNN over augmented base-pairing networks (including non-canonical pairs) of a binding site, predicting a ligand fingerprint. Non-canonical pairs were "the only representation able to uncover a significant signal". It assumes the input is already a binding site — [RNAmigos help](https://rnamigos.cs.mcgill.ca/help); [PMC7430648](https://pmc.ncbi.nlm.nih.gov/articles/PMC7430648). **RNAmigos2** (Nat Commun 2025): coarse-grained 3D, synthetic docking data, active compounds in the top 2.8% on structurally distinct test sets, robust to binding-site variation — [PMC11928640](https://doi.org/10.1038/s41467-025-57852-0). Code [github.com/cgoliver/rnamigos2](https://github.com/cgoliver/rnamigos2) (MIT).
- **fingeRNAt** (PLoS Comput Biol 2022): detects nine interaction types (H-bond, halogen, cation–anion, pi–cation, pi–anion, pi-stacking, ion-mediated, water-mediated, lipophilic) in nucleic acid–ligand complexes and encodes them as SIFt bit strings — [PMC9197077](https://doi.org/10.1371/journal.pcbi.1009783); [github.com/n-szulc/fingeRNAt](https://github.com/n-szulc/fingeRNAt) (GPL-3.0; last push 2022). The authors propose IFP similarity as an alternative to RMSD for finding complexes that share interactions despite different folds, and for grouping ligands that form similar binding networks — [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2021.12.23.474073.full.pdf). A follow-up reports SIFt plus ML models outperforming classic scoring functions in RNA virtual screening — [bioRxiv 2023.01.11.523582](https://www.biorxiv.org/content/10.1101/2023.01.11.523582.full.pdf).
- **ProLIF** supports ligand, protein, DNA or RNA complexes from MD, experimental structures or docking. Interaction definitions are user-extensible and output is a pandas DataFrame — [ProLIF docs via StreamD](https://streamd.readthedocs.io/en/latest/plif.html); [J Cheminf 2021](https://link.springer.com/doi/10.1186/s13321-021-00548-6); [github.com/chemosim-lab/ProLIF](https://github.com/chemosim-lab/ProLIF) (Apache-2.0, active).
- **Protein pocket-comparison tools**:
  - A benchmark of APoc, SiteEngine and G-LoSA found that including structurally similar proteins in datasets overestimates performance. Alignment-based methods cost more than alignment-free ones but give local alignments that inform binding modes — [Govindaraj & Brylinski 2018, BMC Bioinformatics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5845264/).
  - Other tools: PocketMatch (sorted distance lists) — [PMC2639437](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2639437/); a convolution kernel over atom clouds — [BMC Bioinf 2010](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-99); eMatchSite — [PLoS CB](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1003829); SiteMotif (graph-based, multiple-site frameworks) — [PLoS CB 2022](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1009901).
- Beware a name clash: "Structure alignment-based classification of RNA-binding pockets" (SAMO network clustering) concerns protein pockets that bind RNA, not RNA pockets that bind small molecules — [PMC5225598](https://pmc.ncbi.nlm.nih.gov/articles/PMC5225598).
- Recent ML binding-site predictors (GATRsite, JCIM 2025; BC-GNN 2026; CoCoBind 2026; SMARTPocket 2026 preprint) predict sites and do not compare them — [GATRsite](https://doi.org/10.1021/acs.jcim.5c00605); [SMARTPocket](https://doi.org/10.64898/2026.06.18.732920).

### Inferences
- In a CM-anchored framework, the cleanest pocket-similarity measures are position-based:
  - Jaccard or overlap of the Rfam positions contacted;
  - position-indexed interaction fingerprints, i.e. a SIFt or ECIFP polymer side keyed by Rfam position instead of residue number. Different families in a clan can then be compared bit-for-bit wherever the CM maps them.
  - Alignment-free shape or physico-chemical descriptors (fpocketR or RPocket-style volume, buriedness, eigenvalue ratios) are a complementary check where positions do not map.
- Protein tools (G-LoSA, APoc, ProBiS, SiteAlign) assume amino-acid chemistry and protein atom types. No evidence was found of validated use on RNA pockets, so adapting them would be a research project and not a drop-in.
- Equivalence across families should require agreement on several independent axes (CM positions, local superposition, and interaction pattern). Fingerprint similarity alone can link non-homologous sites (see the convergent triples in Q4).

### Gaps
- No RNA pocket comparison benchmark found. No evidence that ProBiS, SiteAlign, G-LoSA, APoc or PocketMatch has been validated on RNA. 3D Zernike descriptors were searched and no RNA-pocket application found.
- "Rsite/RNAsite" descriptor details were not retrieved. These are binding-site predictors, not comparison tools.
- The HARIBOSS current entry count and update status were not verified.

## Q3. Ligand pose comparison and binding-mode clustering (RMSD after site superposition, symmetry-corrected RMSD, IFP clustering), cut-offs

### Takeaway
Standard practice:
1. Superpose sites on pocket atoms.
2. Compute ligand heavy-atom RMSD with symmetry correction by graph isomorphism (spyrmsd, DockRMSD). The naive Hungarian assignment can underestimate RMSD badly.
3. Complement with interaction-fingerprint similarity (Tanimoto or alternatives) and cluster.

No RNA-specific cut-off standard was found. The 2 Å pose-RMSD convention comes from docking, and the 0.7–0.85 Tanimoto heuristic comes from chemical fingerprints, not IFPs.

### Cited Findings
- **spyrmsd** (J Cheminf 2020, Meli & Biggin): symmetry-corrected RMSD by graph isomorphism. On 1DRJ the Hungarian method gives an artificially low 1.00 Å against the correct 2.46 Å — [PMC7457508](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7457508/); [github.com/RMeli/spyrmsd](https://github.com/RMeli/spyrmsd) (MIT, active).
- **DockRMSD** (J Cheminf 2019, Bell & Zhang): symmetry correction as graph-isomorphism search; evaluated on Vina poses from CSAR Hi-Q — [PMC6556049](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6556049/).
- IFPs can be clustered to reduce candidate binding modes to a set small enough for visual inspection, and have also been applied to MD unbinding — [ProLIF paper, J Cheminf 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8466659).
- IFP bits encode (residue, interaction type). Tanimoto against a reference pose is the usual comparison. Rácz, Bajusz & Héberger (J Cheminf 2018) evaluated alternative similarity metrics for IFPs ("Life beyond the Tanimoto coefficient") — [PMC6755604](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6755604/).
- A common ad hoc chemical-fingerprint similarity cut-off is about 0.7–0.85, but it depends on fingerprint type and length. This is a generic heuristic, not an IFP result — [search summary of PMC4456712](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4456712/).
- Applied example (proteins): IFP clustering of 384 RdRP structures gave four binding modes with corresponding pockets — [J Proteome Res 2020](https://doi.org/10.1021/acs.jproteome.0c00623).
- fingeRNAt explicitly proposes IFP similarity to group RNA ligands with similar binding networks, as an RMSD alternative — [fingeRNAt preprint](https://www.biorxiv.org/content/10.1101/2021.12.23.474073.full.pdf).

### Inferences
- For the project's within-site grouping:
  - Keep ECIFP Tanimoto but report it alongside symmetry-corrected ligand RMSD after local pocket superposition (spyrmsd, MIT, pure Python). This applies only to copies of the same CCD, since RMSD needs identical topology.
  - When the two disagree (same contacts but flipped pose, or same pose with different water or ion-mediated contacts), flag the pair rather than silently choosing one.
- Cut-off choice:
  - Derive the cut-off empirically for each site, e.g. from the distribution of pairwise distances, a silhouette scan, or a gap in average-linkage dendrogram heights.
  - Report the chosen value and how stable cluster membership is across nearby values, as the Structures tab does for k.
  - Use complete or average linkage on 1 − Tanimoto rather than single linkage, which chains copies together.
- Polymer-side IFPs keyed by Rfam position make cross-structure fingerprints comparable. Residue-number keys would not.
- For multi-copy ligands in one asymmetric unit, the same ligand copy should not be double-counted with different chain pairings. Treat copy × chain contact sets as the unit, consistent with the project's numerator/denominator rule.

### Gaps
- No RNA-specific benchmark of binding-mode clustering thresholds (RMSD or IFP Tanimoto) found.
- The original Marcou & Rognan 2007 IFP paper and its thresholds were not retrieved.
- No source found on ECIFP specifically applied to RNA binding-mode clustering.

## Q4. Studies comparing the same ligand across different RNA families: how recognition was compared

### Takeaway
Published cross-class comparisons are mostly expert structural comparisons: superposition of the recognition elements (triples) with an RMSD, side-by-side contact maps, and mutational or affinity data. They are not systematic pocket-similarity scores. Two cases arise:
- **Convergent recognition with unrelated folds:** preQ1-II/III/SAM-II triples; SAM-I vs SAM-II vs SAM-III; c-di-GMP-I vs II.
- **Shared scaffold with tuned pockets:** purine vs 2′-dG-I/II/III; Guanine-I vs Guanine-II; SAM-I vs SAM-IV; ykkC guanidine-I vs ppGpp vs PRPP.

The second case is the clan-like one most similar to this project.

### Cited Findings
- **preQ1-I/II/III**: three classes with distinct folds — [Chem Biol 2014](https://www.cell.com/ccbio/fulltext/S1074-5521(14)00203-8).
  - preQ1-II (2.3 Å) binds at a three-way junction flanking a pseudoknotted RBS, a previously uncharacterised fold — [Liberman et al., Nat Chem Biol 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3661761); NMR structure — [Kang et al., PNAS 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC3926045).
  - preQ1-III (2.75 Å, KD,app 6.5 nM) binds between P1 and P2 of a three-way junction — [Liberman et al., PNAS 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4500280).
  - In that paper, preQ1-II and preQ1-III recognise preQ1 with a common constellation of bases on U•A-U triples, superimposing at 1.1 Å RMSD. SAM-II's major-groove triples that read SAM's adenine overlay the preQ1-III triples at 0.82 Å RMSD — [PNAS 2015 / OSTI copy](https://www.osti.gov/pages/servlets/purl/1349644) (figures paraphrased from search summary; full text not read).
  - Despite the shared motif, class II pocket mutations were highly disruptive while four equivalent mutations in class III kept tight binding — [NAR 2024](https://academic.oup.com/nar/article/52/21/13152/7824224).
- **SAM classes**:
  - SAM-I (four-way junction plus k-turn), SAM-III (three-way junction) and SAM-II (pseudoknot) are completely different from each other, and SAM-I and SAM-II bind SAM in different conformations — [Batey-lab review text, PMC3264916 / PMC6061858](https://pmc.ncbi.nlm.nih.gov/articles/PMC6061858).
  - SAM-II structure: 2.8 Å, PDB 2QWY — [Gilbert et al., NSMB 2008](https://www.rcsb.org/structure/2QWY).
  - SAM-IV shares conserved nucleotide positions with SAM-I "despite rearranged structures". Sequence analysis and molecular recognition experiments indicate similar binding sites on different scaffolds — [Weinberg et al., RNA 2008](https://rnajournal.cshlp.org/content/14/5/822.full).
  - SAM-V is predicted to form a pocket similar to SAM-II's but is phylogenetically distinct — [NAR 2018](https://academic.oup.com/nar/article/46/13/6869/5042024).
- **c-di-GMP-I vs II**: c-di-GMP-II (PDB 3Q3Z, ~2.5 Å) binds within a pseudoknot-formed triplex, with non-canonical guanine pairing and no RNA contact to the phosphodiester backbone. The authors describe this as quite different from c-di-GMP-I, i.e. at least two independent solutions. A designed analog binds class II selectively over class I — [Smith et al., PNAS 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3093532). c-di-GMP-I binds at a Y-shaped three-helix junction with a conserved C44/G83 pair bridging P2/P3 — [Smith & Strobel review](https://pmc.ncbi.nlm.nih.gov/articles/PMC3689584).
- **Purine scaffold vs 2′-dG**:
  - A few sequence changes switch a guanine riboswitch to 2′-dG and reshape the pocket to fit the sugar (PDB 3DS7, 1.85 Å) — [Edwards & Batey, JMB 2009 / RCSB 3DS7](https://www.rcsb.org/structure/3ds7).
  - 2′-dG-II keeps the purine-family global architecture but changes the core, with a two-nucleotide junction insertion. It needs only local pocket changes, unlike 2′-dG-I — [Matyjasik & Batey 2019, PMC6847200](https://pmc.ncbi.nlm.nih.gov/articles/PMC6847200).
  - 2′-dG-III structures (2025) compared against 2′-dG-I, 2′-dG-II and Guanine-I report "subtle yet significant" variations controlling affinity and specificity — [RCSB 9LKW](https://rcsb.org/structure/9LKW).
- **Guanine-II vs Guanine-I / Xanthine-II** (NAR 2025): same conserved scaffold, distinct ligand-binding characteristics. Local rearrangements and differences in pocket composition and peripheral architecture determine specificity. Shared recognition patterns with the THF-II riboswitch were used to find new binders — [PMC12445698](https://doi.org/10.1093/nar/gkaf884).
- **ykkC guanidine-I / ppGpp (subtype 2a) / PRPP (2b)**: share a conserved core and nearly overlapping scaffolds in the bound state. Selectivity comes from an extra helical element forming a tunnel for direct and Mg2+-mediated binding. PRPP and ppGpp bind in overlapping positions while guanidine binds higher in the stem. A G96A mutation switches specificity 40,000-fold towards ppGpp — [Knappenberger, Reiss & Strobel, eLife 2018](https://elifesciences.org/articles/36381); also Peselis & Serganov, Nat Chem Biol 2018 (PDB 6DLQ–6DNR; [RCSB 6DMD](https://rcsb.org/structure/6DMD)).

### Inferences
- "Equivalent site" across families needs at least two senses:
  - **Homologous site**: same CM-mapped positions and local superposition. Examples: purine/2′-dG, Guanine-I/II, ykkC subtypes, plausibly SAM-I/IV, although SAM-IV's rearranged scaffold may defeat a CM-to-CM mapping.
  - **Analogous recognition**: a convergent motif such as U•A-U triples reading a purine face, without positional homology. Examples: preQ1-II/III/SAM-II, SAM-I/II/III, c-di-GMP-I/II.

  The project's clan approach (CM mapping plus pocket RMSD) addresses only the first. An IFP or motif-level comparison (base-triple geometry keyed by ligand atom groups) would be needed for the second, and should be labelled differently.
- The literature's comparison recipe could be implemented as a reproducible output: superpose only the recognition nucleotides, report RMSD with the number of atoms, plus a contact table by ligand functional group. Labs already report numbers like 0.82 Å and 1.1 Å over triples this way.
- The ykkC and 2′-dG examples show that a ligand can bind in an overlapping but shifted location within a homologous scaffold. Comparing pockets requires position-level overlap statistics (shared/total contacted positions), not a binary same/different.

### Gaps
- No systematic computational survey found that compares all ligand-sharing riboswitch classes with a single pocket-similarity metric.
- SAM-III structure paper (Lu et al. 2008) and SAM-I/SAM-II pocket-level comparison text not retrieved directly. NAD+ (I/II), Mn2+, and THF-I/II cross-class comparisons not searched.
- Full texts of the preQ1-III paper (exact RMSD atom sets) and the 2′-dG-III paper (journal and DOI) were not read. Figures come from search summaries.
