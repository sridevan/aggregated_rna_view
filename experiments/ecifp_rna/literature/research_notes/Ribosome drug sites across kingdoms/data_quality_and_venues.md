# Data quality pitfalls, recent comparable studies and venues for an archive-wide ribosome–small-molecule survey

Retrieval date: 2026-10-24. Bibliographic details (authors, journal, volume, DOI) were checked against Europe PMC REST records unless a line is marked **[unverified]**. Full texts were mostly paywalled or blocked (PMC returned a browser check and nature.com a login redirect), so most findings rest on abstracts and indexed metadata.

## Q1. Ligand model quality in ribosome structures, and validation metrics usable at scale

### Takeaway
There are well-documented cases of antibiotic models in early ribosome crystal structures that disagree with each other and with later data. The clearest is chloramphenicol, and macrolides in eubacterial versus archaeal 50S. Bulkley et al. (2010) attribute these disagreements to how the electron density was interpreted, not to differences between species. An archive-wide survey therefore needs per-ligand density-fit filters: RSCC/RSR for X-ray, and Q-score or atom inclusion for cryo-EM, both from the wwPDB validation reports. It should also keep flagged copies in the data rather than drop them, and treat 2000–2010 crystal structures with caution.

### Cited Findings
- Bulkley D, Innis CA, Blaha G, Steitz TA. "Revisiting the structures of several antibiotics bound to the bacterial ribosome." PNAS 2010;107(40):17158–17163. doi:10.1073/pnas.1008685107. Their T. thermophilus 70S chloramphenicol model "differs substantially from the prevailing model". Their erythromycin, azithromycin and telithromycin structures match the archaeal (H. marismortui) complexes but differ from some published eubacterial (D. radiodurans) models. They argue the discrepancies come from how the electron density was interpreted, not from species, and that macrolide binding is "largely conserved across species". — [Europe PMC record](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE:%22revisiting%20the%20structures%20of%20several%20antibiotics%22&resultType=core&format=json); [PNAS](https://pnas.org/content/107/40/17158.full); [PMC2951403](https://pmc.ncbi.nlm.nih.gov/articles/PMC2951403)
- A later higher-resolution T. thermophilus 70S study of chloramphenicol and erythromycin reported "unambiguous placement" of the drugs. It still could not place the dichloroacetic moiety of chloramphenicol unambiguously, so parts of a ligand can stay uncertain even in better maps. — [OSTI full text](https://www.osti.gov/pages/servlets/purl/1543007) (seen only as a search snippet; title and authors **[unverified]**, probably Svetlov/Polikanov et al. 2019)
- Paternoga H, …, Wilson DN. "Structural conservation of antibiotic interaction with ribosomes." Nat Struct Mol Biol 2023;30(9):1380–1392. doi:10.1038/s41594-023-01047-y. Cryo-EM structures of 17 compounds from 6 antibiotic classes at 1.6–2.2 Å resolve water-mediated contact networks. Antibiotics that share a scaffold bind in highly similar ways, and some waters are pre-ordered in the site before the drug binds. This gives a modern high-resolution reference set against which older models could be checked. — [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE:%22structural%20conservation%20of%20antibiotic%20interaction%20with%20ribosomes%22&resultType=core&format=json); [Nature](https://www.nature.com/articles/s41594-023-01047-y)
- An in-cell cryo-EM study (M. pneumoniae 70S with chloramphenicol, 3.7 Å, multi-particle refinement in M) used a negative control: fitting the same model into a 5.6 Å map of untreated ribosomes shows no chloramphenicol density. — [ResearchGate record](https://www.researchgate.net/publication/341968396_Multi-particle_cryo-EM_refinement_with_M_visualizes_ribosome-antibiotic_complex_at_37_A_inside_cells) (Tegunov et al., Nat Methods 2021 **[citation details unverified]**)
- General critiques of ligand models (not specific to ribosomes):
  - Weichenberger, Pozharski and Rupp, "Visualizing ligand molecules in twilight electron density", Acta Cryst F 2013 ([PMC3564628](https://pmc.ncbi.nlm.nih.gov/articles/PMC3564628)), with a companion paper by Pozharski, Weichenberger and Rupp in Acta Cryst D 2013;69:150–167. Both use RSCC as the key metric and a cutoff of 0.6 to judge whether the experimental data support a ligand. A secondary summary reports that 7.6% of ligands in their sample had RSCC < 0.6 ([blopig summary](https://blopig.com/blog/?p=442); the figure is **[unverified against the primary paper]**). The Twilight GUI is discontinued, but its tables can still be downloaded ([ruppweb.org/twilight](https://ruppweb.org/twilight/default.htm)).
  - Deller MC, Rupp B. "Models of protein-ligand crystal structures: trust, but verify." J Comput Aided Mol Des 2015;29(9):817–836. doi:10.1007/s10822-015-9833-8. A review of overenthusiastic interpretation of ligand density. — [Europe PMC](https://www.ebi.ac.uk/europepmc/)
  - Wlodawer A, Dauter Z, Porebski PJ, Minor W, Stanfield R, Jaskolski M, Pozharski E, Weichenberger CX, Rupp B. "Detect, correct, retract: How to manage incorrect structural models." FEBS J 2018;285(3):444–466. doi:10.1111/febs.14320. — [Europe PMC](https://www.ebi.ac.uk/europepmc/)
  - Weichenberger, Pozharski, Rupp. "Twilight reloaded: the peptide experience." Acta Cryst D 2017;73:211–222. — [IUCr](https://journals.iucr.org/d/issues/2017/03/00/ba5251/index.html)
- Kowiel M, Brzezinski D, Porebski PJ, Shabalin IG, Jaskolski M, Minor W. "Automatic recognition of ligands in electron density by machine learning" (CheckMyBlob). Bioinformatics 2019;35(3):452–461. doi:10.1093/bioinformatics/bty626. Benchmarked on up to 219,931 binding sites for the 200 most common PDB ligands. It works on X-ray density only. — [Europe PMC](https://www.ebi.ac.uk/europepmc/)
- Pintilie G, Zhang K, Su Z, Li S, Schmid MF, Chiu W. "Measurement of atom resolvability in cryo-EM maps with Q-scores." Nat Methods 2020;17(3):328–334. doi:10.1038/s41592-020-0731-1. Averaged Q-scores track map resolution for both protein and RNA. — [Europe PMC](https://www.ebi.ac.uk/europepmc/)
- Lawson CL, Kryshtafovych A, … Chiu W. "Outcomes of the EMDataResource cryo-EM Ligand Modeling Challenge." Nat Methods 2024;21(7):1340–1348. doi:10.1038/s41592-024-02321-7.
  - Scope: 61 models from 17 groups, built into three maps at 1.9–2.5 Å.
  - Main conclusion: "A composite rather than a single score was needed to assess macromolecule+ligand model quality."
  - Recommended checks: Q-score, EMRinger and map-model FSC. For ligands also occupancy, missing atoms, surrounding ions and strain energy. The challenge also introduced LIVQ, which scores the density fit of the ligand and its surrounding residues.
  - [PMC10854310](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10854310/); [EMDataResource news](https://www.emdataresource.org/news/ligand-challenge-publication.html); [RCSB news](https://www.rcsb.org/news/6679793750bece1b815a1450)
- One deep-learning ligand-identification study reports that "most cryoEM ligands have Q-scores below 0.8, whereas most X-ray ligands have RSCC between 0.9 and 1.0". It argues that additional ligand validation metrics are needed for cryo-EM. — [PMC11383698](https://pmc.ncbi.nlm.nih.gov/articles/PMC11383698/) (seen as a search summary only; authors **[unverified]**)
- Current wwPDB EM validation reports include Q-scores (for the whole model and per chain) and atom inclusion at the recommended contour level. — example reports [8EEW (2026)](https://files.rcsb.org/validation/view/8eew_full_validation.pdf), [9LJJ (2025)](https://files.rcsb.org/validation/view/9ljj_full_validation.pdf). Whether these reports give per-ligand Q-scores for every ligand, and whether those values are available in the validation XML for bulk parsing, was not confirmed **[gap]**.
- van Zundert GCP et al. "Macromolecular refinement of X-ray and cryoelectron microscopy structures with Phenix/OPLS3e for improved structure and ligand quality." Structure 2021. doi:10.1016/j.str.2021.03.011. This is an example of re-refinement improving ligand quality at scale (scope unverified). — [Europe PMC listing](https://www.ebi.ac.uk/europepmc/)

### Inferences
- Use the validation report XML (RSCC/RSR/LLDF for X-ray; Q-score and atom inclusion for EM) as per-copy quality flags rather than as hard filters. Report results with and without the flags (for example RSCC ≥ 0.8 for X-ray, ligand Q-score ≥ 0.5–0.6 for EM). The cutoffs would need calibrating, because cryo-EM ligand Q-scores sit systematically lower than X-ray RSCC.
- Treat the D50S macrolide and chloramphenicol entries from 2001–2003 as a known test case for any quality flagging: the survey should flag them, or explain why it does not.
- When one compound has both an old and a new structure (for example chloramphenicol in a 2001 D50S crystal versus the 2023 1.6–2.2 Å cryo-EM), a survey can measure how consistently the binding pose is reproduced. That would add value beyond simply listing contacts.

### Gaps
- I found no systematic, archive-wide re-refinement (PDB-REDO-style) of ribosome–antibiotic crystal structures. PDB-REDO generally does not process the largest ribosome entries, but that is from memory and **[unverified]**.
- I could not confirm a published distribution of ligand Q-scores for cryo-EM ribosome entries specifically.
- No full text was obtained for Bulkley 2010 or Paternoga 2023 (blocked), so specific per-entry misplacements beyond chloramphenicol and macrolides were not extracted.

## Q2. Redundancy in ribosome–ligand structures and how surveys deduplicated

### Takeaway
None of the sources I found gives a ribosome-specific deduplication protocol. RNA–ligand resources either collapse ribosomes drastically, deduplicate by sequence clustering, or apply a drug-likeness filter and keep everything. Two recent resources also treat structure-, sequence- and pocket-level clusters as separate axes.

### Cited Findings
- DRLiPS (IIT Madras) handles ribosomes explicitly. Of 513 ribosomal RNA–ligand complexes, it retained only 51 (the 16S A-site) alongside the non-ribosomal ones. It defines a pocket as the residues within 6.5 Å of the ligand. Deduplication ran in two steps: 819 pockets remained after removing copies within homo-multimers, and 56 after Rfam-based clustering. — [DRLiPS workflow](https://web.iitm.ac.in/bioinfo2/DRLiPS/Details.html)
- RPocket clustered 269 structures with CD-HIT at 0.80 and 0.95 identity, giving 75 and 94 clusters respectively. The final set is 240 pockets from 94 non-redundant complexes. Zhou T et al., BMC Bioinformatics 2021, doi:10.1186/s12859-021-04349-4. — [PMC8424408](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8424408/)
- HARIBOSS (Panei FP, Torchet R, Ménager H, Gkeka P, Bonomi M. Bioinformatics 2022;38(17):4185–4193, doi:10.1093/bioinformatics/btac483) selects structures with a drug-likeness filter. The sources I found do not describe ribosome-specific handling ([bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full.pdf)). HARIBOSS++ (Marengo M et al., bioRxiv 2026, doi:10.64898/2026.09.07.749797) extends the pocketome to sites where the ligand also contacts other biomolecules, which matters for ribosomes where both rRNA and r-proteins line the pockets. It also adds RNA-family analysis and chemical-space analysis. — [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE:%22HARIBOSS%22&resultType=core&format=json)
- SMARTFlexDB assigns separate sequence-, structure- and pocket-level clusters so that dataset splits respect redundancy. — [SMARTFlexDB](https://aidd.rc.ufl.edu/app/smartflexdb/) (ribosome handling unknown)
- BGSU Representative Sets give equivalence classes for RNA chains at several resolution thresholds, updated weekly. This could serve as a chain-level deduplication layer, though it ignores ligands. — [BGSU](https://www.bgsu.edu/research/rna/databases/non-redundant-list.html); RNA-NRD deduplicates by sequence and 3D similarity ([PMC10132383](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10132383/))

### Inferences
- For ribosomes, redundancy runs along four axes:
  - two 70S copies per asymmetric unit in T. thermophilus crystals;
  - long series of near-identical entries from the same lab (for example the Polikanov, Steitz and Yonath series);
  - the same drug across species;
  - several sites for one drug.
  Collapsing hierarchically fits the project's Rfam-position coordinate system well: copy → entry → (drug, site defined on Rfam/rRNA positions, species) → (drug, site). Every count should then carry its numerator and denominator at the chosen level.
- The DRLiPS figure (513 → 51) shows that generic RNA–ligand resources mostly discard ribosomes. A survey that keeps ribosomes and treats them carefully fills a gap rather than duplicating existing work.

### Gaps
- I found no published survey that documents a rule for choosing one copy between the two ribosomes in the asymmetric unit (for example the copy with the better ligand RSCC, or keeping both as replicates).
- How HARIBOSS handles ribosomes and redundancy could not be confirmed from the abstracts.

## Q3. Coverage gaps in contact annotations, and ion/water modelling at cryo-EM resolution

### Takeaway
Ion identity in ribosome models is unreliable. Mg²⁺ is systematically over-assigned, and at least one of the two re-evaluated models also misses ions. Stereochemistry at better than ~2.4 Å is the only robust check for Mg²⁺/Na⁺/K⁺, so contacts and "water-mediated" interactions that involve ions should be reported with caution and stratified by resolution. Coverage of large cryo-EM entries by PDBe and RCSB contact services was not documented in the external literature I found. The project's own findings are the best evidence here.

### Cited Findings
- Auffinger P, Ennifar E, D'Ascenzo L. "Deflating the RNA Mg2+ bubble: stereochemistry to the rescue!" RNA 2021;27(3):243–252. doi:10.1261/rna.076067.120. They re-evaluated the ribosome models 6QNR and 6SJ6. "Numerous studies overestimate the number of Mg2+" in deposited structures. Their analysis substantially reduced the number of Mg²⁺ assigned in 6QNR and suggests K⁺ is the dominant ion in the first solvation shell. "stereochemistry coupled with high resolution structures (<2.4 Å) remains the best currently available option". — [RNA journal](https://rnajournal.cshlp.org/content/27/3/243)
- Leonarski F, Henning-Knechtel A, Kirmizialtin S, Ennifar E, Auffinger P. "Principles of ion binding to RNA inferred from the analysis of a 1.55 Å resolution bacterial ribosome structure – Part I: Mg2+." Nucleic Acids Res 2025;53(1):gkae1148. doi:10.1093/nar/gkae1148.
  - Data: re-modelling of the 1.55 Å E. coli 70S cryo-EM structure 8B0X (Fromm et al., Nat Commun 2023;14:1095, doi:10.1038/s41467-023-36742-3), deposited as 9Q87.
  - Method: a 2.9 Å distance between two waters bound cis to Mg²⁺ lets them assign every Mg²⁺ that has 2–4 non-water oxygen ligands.
  - Framing: incomplete ion assignments in the original model "blurred our understanding".
  - [9Q87 at RCSB](https://rcsb.org/experimental/9Q87); [8B0X at RCSB](https://rcsb.org/structure/8B0X)
- Leonarski, D'Ascenzo, Auffinger (RNA 2019) found that nucleobase carbonyls are poor inner-sphere Mg²⁺ binders but excellent binders of monovalent ions. They note "a substantial number of assignment issues" in the PDB. — [PMC6348993](https://pmc.ncbi.nlm.nih.gov/articles/PMC6348993/)
- Paternoga et al. 2023 (above) shows that structured water networks in antibiotic sites are only resolvable at about 2 Å or better.
- Internal project finding (not published literature): PDBe Arpeggio/bound-ligand interaction data cover only the ligand copies in the preferred assembly, while RCSB 5 Å neighbour lists cover all copies (project memory "ligand-contact-sources"). `ligands.json` is built from PDBe entity-level contacts, and RF00167 showed copy gaps, including a wrong bestChainId for 5SWD (project memory "ligands-tab-investigation"). See also `experiments/ligand_annotation_investigation/REPORT.md`.

### Inferences
- An archive-wide ribosome survey should compute contacts itself from the coordinates (for every copy, with explicit distance and geometry rules, as Arpeggio, PLIP or the project's ECIFP does). It should use PDBe and RCSB annotations only for cross-checking, and record which source and which ligand copy each contact comes from.
- Ion- and water-bridged contacts should be a separate class, restricted to entries better than about 2.0–2.4 Å. Metal identity (Mg versus K/Na) should not be treated as reliable in older or lower-resolution models.

### Gaps
- I found no external paper that quantifies missing PDBe Arpeggio or RCSB interaction annotations for large ribosome entries. The claim rests on project-internal checks.
- I found no paper that quantifies how often waters are modelled, or how reliably, across cryo-EM ribosome entries by resolution bin.

## Q4. Recent (2022–2026) large-scale or meta analyses relevant to ribosome-targeting compounds

### Takeaway
Recent work clusters in three areas:
- high-resolution structural series (Paternoga 2023);
- evolutionary variation of drug-binding residues across bacteria (Ekemezie 2025);
- general RNA–small-molecule databases and meta-analyses (HARIBOSS/HARIBOSS++, SPIRAL 2026, fpocketR 2026, RNALigands, R-SIM).

I found no recent archive-wide structural survey of ribosome-bound small molecules across kingdoms, and no ribosome-specific docking or ML benchmark. Both are plausible gaps.

### Cited Findings
- Ekemezie CL, Chan LI, Brown CR, Helena-Bueno K, Williams TA, Melnikov SV. "Extensive natural variation in bacterial ribosomal drug-binding sites." Cell Rep 2025;44(7):115878. doi:10.1016/j.celrep.2025.115878. Many rRNA drug-binding residues regarded as bacteria-specific are conserved only in subsets of bacteria. Divergent sites are widespread and trace back to ancient polymorphisms. They show intrinsic resistance experimentally in some of these species. This is the closest recent precedent for "site conservation across pathogens and kingdoms", but it is sequence- and evolution-based rather than an archive-wide structural survey. — [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%2210.1016/j.celrep.2025.115878%22&resultType=core&format=json)
- Paternoga et al. 2023 NSMB (above): a multi-compound high-resolution structural comparison with conserved binding modes and water networks.
- Nandi S et al. "Distant ribose 2'-O-methylation of 23S rRNA helix 69 pre-orders the capreomycin drug binding pocket…" Nucleic Acids Res 2025. doi:10.1093/nar/gkaf618. Shows rRNA modifications shaping a pocket, which matters for annotating modified nucleotides in drug sites. — [Europe PMC listing]
- Lu X, Wang Y. "Structural Pockets and Interacting RNA-Associated Ligands (SPIRAL): a DSSR-enabled meta-analysis of RNA-small molecule recognition." bioRxiv 2026. doi:10.64898/2026.05.19.726393. Covers 1,098 PDB structures and 1,137 binding events across six RNA functional categories. Fingerprints from DSSR (stacking, H-bond topology, groove, tertiary motif) cluster into six binding modes. A 7-metric "Composite Binding Quality Score" ranks riboswitches highest. Across 275 entries with measured affinities, the C2′-endo count and buried contact area were the dominant affinity predictors. Redundancy handling is not stated in the abstract. — [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%2210.64898/2026.05.19.726393%22&resultType=core&format=json)
- HARIBOSS (Bioinformatics 2022) and HARIBOSS++ (bioRxiv 2026): above.
- Veenbaas SD et al. "fpocketR: a platform for identification and analysis of ligand-binding pockets in RNA." ACS Chem Biol 2026. doi:10.1021/acschembio.5c00805 (preprint 2025, doi:10.1101/2025.03.25.645323).
- Padroni G, Patwardhan NN, Schapira M, Hargrove AE. "Systematic analysis of the interactions driving small molecule-RNA recognition." RSC Med Chem 2020;11(7):802–813. doi:10.1039/d0md00167h. Stacking and H-bonds dominate RNA recognition, cation–π interactions favour G and purines, and contacts to the backbone are scarce even for cationic ligands.
- Other RNA–ligand resources:
  - Sun S et al. "RNALigands." RNA 2022. doi:10.1261/rna.078889.121
  - Ramaswamy Krishnan S et al. "R-SIM." J Mol Biol 2023. doi:10.1016/j.jmb.2022.167914
  - Nithin C et al. "Comparative analysis of RNA 3D structure prediction methods: towards enhanced modeling of RNA-ligand interactions." Nucleic Acids Res 2024. doi:10.1093/nar/gkae541
  - Childs-Disney JL et al. "Targeting RNA structures with small molecules." Nat Rev Drug Discov 2022 (review)
  - [Europe PMC listing](https://www.ebi.ac.uk/europepmc/)
- Older baseline: "Structural signatures of antibiotic binding sites on the ribosome," Nucleic Acids Res 2010;38(18):5982. It analysed 65 antibiotic sites from ribosome–antibiotic PDB complexes. The most prominent features were unpaired bases, a strong bias towards syn bases and unusual sugar puckers. — [NAR](https://academic.oup.com/nar/article/38/18/5982/1065889) (authors David-Eden H, Mankin AS, Mandel-Gutfreund Y **[from memory, unverified]**)
- ML and docking: a 2025 review notes that RNA-targeting ML is fragmented and that benchmarks remain protein-centric (PDBbind, PLINDER, …) ([arXiv 2512.15645](https://arxiv.org/html/2512.15645v1)). RmsdXNA (Brief Bioinform 2024) rescores docking poses for nucleic acid–ligand complexes. An AlphaFold-enabled docking benchmark for antibiotics reported a mean auROC of 0.48, rising to 0.63 with ML rescoring, but against E. coli proteins rather than the ribosome ([ResearchGate](https://www.researchgate.net/publication/363306375_Benchmarking_AlphaFold-enabled_molecular_docking_predictions_for_antibiotic_discovery); Wong F et al., Mol Syst Biol 2022 **[unverified]**). A 2026 arXiv paper proposes a pretrained geometric deep-learning docking method for nucleic acids ([arXiv 2606.05198](https://arxiv.org/pdf/2606.05198); benchmark details unverified).

### Inferences
- The closest precedents are:
  - Ekemezie 2025 for conservation;
  - Paternoga 2023 for structural conservation of binding modes;
  - SPIRAL and HARIBOSS++ for RNA–ligand meta-analysis.
  A cross-kingdom structural survey that maps contacts onto a common rRNA coordinate system and reports conservation of contact residues by kingdom or pathogen would extend all three. A benchmark split with explicit redundancy and quality flags would add a further layer.
- SPIRAL and HARIBOSS++ (both 2026 preprints) are direct competitors for "RNA–ligand meta-analysis". Novelty would have to come from treating ribosomes properly (copies, assemblies, mixed protein and RNA pockets, cross-kingdom mapping) and from quality-aware statistics.

### Gaps
- I found no 2022–2026 systematic review that catalogues the whole structural landscape of ribosome–antibiotic complexes in the PDB (counts by class, site and species). Narrative reviews exist (for example "Mechanistic Insights into Clinically Relevant Ribosome-Targeting Antibiotics", 2024, [ResearchGate](https://www.researchgate.net/publication/384717662_Mechanistic_Insights_into_Clinically_Relevant_Ribosome-Targeting_Antibiotics), venue unverified).
- I found no cheminformatics analysis specific to ribosome inhibitors (chemical space versus PDB coverage) in the searches.
- I found no ribosome-pocket docking or ML benchmark.

## Q5. Venues and precedent for comparable structural-bioinformatics surveys

### Takeaway
Comparable RNA–ligand resources and surveys appear in the following venues, by type of contribution:

| Contribution | Venues with precedent |
|---|---|
| Database / resource | Bioinformatics (HARIBOSS), RNA (RNALigands), BMC Bioinformatics (RPocket), J Mol Biol (R-SIM), NAR (2010 ribosome site signatures; 2024 RNA–ligand modelling comparison) |
| Chemistry-facing systematic analysis | RSC Med Chem (Padroni 2020), ACS Chem Biol (fpocketR) |
| Biological insight | Cell Reports (Ekemezie 2025), NSMB (Paternoga 2023) |

What made these publishable was one of three things: a reusable curated resource with a web interface, a clear biological finding, or a quantitative benchmark or method.

### Cited Findings
- HARIBOSS: a curated resource plus web interface (hariboss.pasteur.cloud), published in Bioinformatics 2022. — [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full.pdf); [Pasteur](https://research.pasteur.fr/en/b/qJN)
- RNALigands: database and web server, RNA 2022. RPocket: database, BMC Bioinformatics 2021. R-SIM: affinity database, J Mol Biol 2023. — Europe PMC records above
- "Structural signatures of antibiotic binding sites on the ribosome": a structural-feature survey of 65 sites, NAR 2010. — [NAR](https://academic.oup.com/nar/article/38/18/5982/1065889)
- Padroni et al. 2020: interaction statistics for medicinal chemists, RSC Med Chem.
- Ekemezie et al. 2025: biological insight plus experimental validation of intrinsic resistance, Cell Reports.
- Leonarski et al. 2025: a principles paper that re-analyses one deposited ribosome structure, NAR. It shows that NAR accepts careful re-analysis of archive data when it yields general rules.
- Nithin et al. 2024: a benchmark-style comparative analysis focused on RNA–ligand interactions, NAR.

### Inferences
- Plausible venues by framing:
  - **Resource or database** (ribosome drug-site atlas on Rfam/rRNA positions, per-copy quality flags, web viewer): NAR Database issue (needs a maintained web resource), NAR Genomics & Bioinformatics, Bioinformatics, or Database.
  - **Biological insight** (site conservation across kingdoms and pathogens, how contacts map onto Ekemezie-style polymorphisms, eukaryotic and mitochondrial off-target liability): NAR main section, RNA, Structure, or PLoS Comput Biol.
  - **Benchmark or ML** (redundancy- and quality-aware ribosome pocket set for docking and ML): JCIM, J Chem Inf Model / J Cheminform, or Bioinformatics.
- Reviewers in this area are likely to expect three things:
  - explicit deduplication and quality filtering, with sensitivity analysis (the crystallographic critiques above make this an obvious review point);
  - a comparison against HARIBOSS/HARIBOSS++ and SPIRAL;
  - at least one non-obvious biological result or a reusable dataset with documented provenance.

### Gaps
- Journal acceptance criteria were not fetched. Statements about what reviewers expect are inferences, not citations.
- I found no ribosome-specific structural-bioinformatics survey from 2022–2026 in Structure, Protein Science or NAR Genomics & Bioinformatics, so there is no direct venue precedent for exactly this study.
