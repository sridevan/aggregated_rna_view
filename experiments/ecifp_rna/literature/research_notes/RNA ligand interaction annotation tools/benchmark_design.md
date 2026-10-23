# Benchmarking RNA–ligand contact annotation tools: evidence for benchmark design

Research date 2026-10-23. PMC article pages (pmc.ncbi.nlm.nih.gov) were blocked by a browser check, so PLIP 2015 was read from the OUP page, fingeRNAt from the PLOS page plus its GitHub `config.py`, and two arXiv PDFs were text-extracted locally. Cutoffs for PLIP and fingeRNAt that were already collected in the sibling note `../RNA small molecule 3D interactions/fingerprints_contact_typing.md` are reused here with their original URLs. Candidate tools in scope: PDBe-Arpeggio, PLIP, fingeRNAt, ProLIF, LUNA, BINANA, GetContacts.

## 1. Which published studies compare interaction-detection tools head to head, on what datasets, and with what results?

### Takeaway
No published, systematic head-to-head benchmark of interaction-detection tools exists for nucleic acid–ligand complexes. For proteins, the comparisons that exist are small and secondary. The best documented is Kokh et al. 2020: 40 complexes, consensus of PLIP, FLIP, LPC and MOE, hydrophobic contacts excluded. It found good agreement for H-bonds and aromatic contacts and the most disagreement for water bridges. Tool papers (PLIP, fingeRNAt) validate against literature-reported interactions or compare features only, without agreement statistics.

### Cited Findings
- **Kokh et al. 2020 (MD-IFP), *J Chem Phys* 153:125102.** This is the only quantitative multi-tool comparison found.
  - **Dataset:** 50 PDB entries originally used to validate FLIP, reduced to 40 after preparation failures. Hydrogens were added with Chimera 1.14 ("unspecified" protonation, H-bond-aware). Only altloc A was kept. Ligands were converted to mol2, with manual fixes for azole, amidine and urea groups.
  - **Rules:** an MD-IFP interaction counted as a false positive if none of PLIP, FLIP, LPC or MOE found it, and as a false negative if all four did. LPC cannot detect water bridges or halogen bonds, so for those the consensus was FLIP+PLIP+MOE. Hydrophobic contacts were excluded "since their definition is not as clear cut", and metal and cofactor interactions were also excluded.
  - **Results:** of 250 non-hydrophobic interactions, 5 were FP (2 H-bonds, 3 aromatic) and 6 were FN (5 H-bonds, 1 salt bridge). Aromatic contacts had 3 FP among 56. Of 69 water bridges, 23 were found by none of FLIP, PLIP or MOE, and 4 were found by all three but missed by MD-IFP. The authors attribute the water-bridge inconsistency to "differences in the approaches used and in the hydrogen bond parameters employed, as well the ambiguity of the assignment of hydrogen orientation (which depends on the procedure used for protonation)". Except for MOE, none of the benchmark tools distinguishes donated from accepted H-bonds.
  - **HSP90 sub-comparison** (3 complexes; OpenEye OEChem vs RCSB ligand-interaction tool vs PLIP vs MD-IFP): "Apart from HY [hydrophobic], there are only a few differences". Hydrophobic counts differ because OEChem types every carbon as hydrophobic, whereas PLIP and RDKit exclude carbons next to N/O. Some H-bonds were missed because the distance was "slightly longer than the H-bond detection threshold".
  - Sources: [arXiv 2006.11066 (PDF)](https://arxiv.org/pdf/2006.11066); [arXiv abstract with journal reference](https://arxiv.org/abs/2006.11066)
- **PLIP 2015 validation.** A test suite of 30 literature-validated protein–ligand complexes, resolution 1.2–3.3 Å, covering all PLIP interaction types. Each test checks that every interaction reported in the original paper is detected (pass/fail). The page gives no recovery percentage and no comparison with other tools. Thresholds were taken from the literature, "carefully adapted", and made somewhat more permissive for low-quality structures — [PLIP, NAR 2015](https://academic.oup.com/nar/article/43/W1/W443/2467865)
- **PLIP 2021 (DNA/RNA)** describes Arpeggio as "most closely related to PLIP offering comparable functionality except for flexible adjustment of thresholds". No quantitative cross-tool benchmark was found in search summaries — [PLIP 2021, NAR](https://academic.oup.com/nar/article/49/W1/W530/6266421)
- **fingeRNAt (Szulc et al., PLoS Comput Biol 2022).**
  - Compares Arpeggio, PLIP2021 and ProLIF on features only (S1 Table); it reports no quantitative agreement benchmark against them.
  - Validation is against literature statistics (de Freitas & Schapira; Padroni et al.; Kligun & Mandel-Gutfreund; Kondo & Westhof; Lu et al.), plus clustering (HIV-1 TAR actives/decoys) and redocking.
  - Across 207 structures it counted 5026 H-bonds, 3582 lipophilic, 899 cation–anion, 151 water-mediated, 146 π-stacking and 6 halogen contacts. Padroni et al. had reported 764 π-stacking and 762 H-bond interactions, a ranking that differs from fingeRNAt's (π-stacking about as frequent as H-bonds there, H-bonds dominant here).
  - S3 Table compares H-bond detection with the default heavy-atom method against detection using hydrogen positions added by OpenBabel, RDKit or PyMOL.
  - Source: [fingeRNAt, PLoS Comput Biol](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009783)
- **LigPlot/LigPlot+** show only H-bonds and hydrophobic contacts, so they cannot serve as a full-coverage comparator — [fingeRNAt, PLoS Comput Biol](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009783)
- **Errington et al. 2024 ("Assessing interaction recovery of predicted protein-ligand poses").**
  - Uses ProLIF interactions in the crystal pose as the reference, with a "PLIF recovery rate" metric: the percentage of crystal interactions reproduced in a predicted pose.
  - Thresholds: H-bond 3.7 Å, cation–π 5.5 Å, ionic 5 Å, other ProLIF v2.0.3 defaults. Hydrophobic and vdW contacts were excluded as non-directional and "highly correlated with RMSD".
  - The authors state that PLIF libraries are "very sensitive to the protonation state of both the protein and the ligand as it can decide whether an interaction gets labelled as ionic or hydrogen bond". They note that PLIP and Arpeggio add hydrogens but do not optimise the H-bond network by default.
  - Dataset: PoseBusters (protein–ligand).
  - Source: [arXiv 2409.20227](https://arxiv.org/pdf/2409.20227)
- A deep-learning study describes PLIP and ProLIF as the two most widely used open tools. In its setting, performance correlated with the level of detail the tools provide (worst for PLIP). This is downstream task performance, not detection accuracy — [bioRxiv 2023.11.13.566916](https://www.biorxiv.org/content/10.1101/2023.11.13.566916v1.full.pdf)
- **BINANA 2** assigns donors/acceptors and charged groups from residue and atom names, assuming standard PDB nomenclature. This makes it sensitive to ligand atom naming and to non-standard residues — [BINANA 2, JCIM](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01461)
- **LUNA** (Keiser Lab) was validated by an ML proxy task: its EIFP fingerprints reproduced DOCK3.7 scores for about 1M D4-docked molecules (EIFP-4096 R² = 0.61). This is not validation of contact detection — [LUNA docs](https://luna-toolkit.readthedocs.io/en/latest/overview.html); [bioRxiv 2022.05.25.493419](https://www.biorxiv.org/content/10.1101/2022.05.25.493419.full.pdf)
- A 2024 PCCP article, "Elucidation of protein–ligand interactions by multiple trajectory analysis methods", appears to compare several analysis methods on MD trajectories. The page returned 403 and its content was not verified — [RSC PCCP 26:6903](https://pubs.rsc.org/cp/article/26/8/6903/802169/Elucidation-of-protein-ligand-interactions-by)

### Inferences
- The design precedent to copy is the Kokh consensus scheme with per-type FP/FN, reported separately per interaction type, with hydrophobic/proximal contacts analysed separately or excluded. For RNA it should be extended with pairwise tool–tool agreement (Jaccard per type) and an independent reference, so that tool consensus is not the only "truth".
- An RNA-specific benchmark would be novel: no paper reports PLIP vs Arpeggio vs fingeRNAt vs ProLIF agreement on RNA–ligand complexes.
- The fingeRNAt vs Padroni discrepancy in π-stacking counts shows that stacking definitions alone can reorder interaction-type frequencies. Stacking should be a primary stratum in the benchmark.

### Gaps
- No published head-to-head study of PLIP vs Arpeggio vs LigPlot+/HBPLUS vs Schrödinger vs MOE with agreement metrics (kappa or Jaccard) was found, for proteins or nucleic acids. GetContacts has no validation paper (not retrieved).
- The fingeRNAt S1 (tool features) and S3 (H-bond method comparison) table contents were not retrieved.
- The PCCP 2024 multi-method comparison could not be read (403).

## 2. Metrics used (precision/recall vs reference, tool–tool agreement, per-type, cutoff sensitivity, hydrogen placement, protonation)

### Takeaway
Published practice is thin. It consists of pass/fail recovery of literature-reported interactions (PLIP), FP/FN against a multi-tool consensus (Kokh), and recovery rate of reference interactions (Errington). Cutoff sensitivity and protonation sensitivity are acknowledged qualitatively but rarely quantified. fingeRNAt's S3 Table (heavy-atom vs H-based H-bonds) is the only RNA-side sensitivity test found.

### Cited Findings
- Recovery rate = percentage of reference (crystal) interactions reproduced. Reported as thresholds of ≥50% and 100% recovery — [Errington et al., arXiv 2409.20227](https://arxiv.org/pdf/2409.20227)
- Consensus-based FP/FN per type (H-bond, aromatic, salt bridge, water bridge, halogen), with water bridges reported separately as the least consistent — [Kokh et al., arXiv 2006.11066](https://arxiv.org/pdf/2006.11066)
- Pass/fail detection of all literature-reported interactions per complex — [PLIP 2015](https://academic.oup.com/nar/article/43/W1/W443/2467865)
- Fingerprint similarity (Tanimoto/Tversky) vs ligand RMSD in redocking gave a weak correlation: mean R² 0.417 (Tanimoto) and 0.380 (Tversky) over 144 complexes — [fingeRNAt](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009783)
- ProLIF authors on cutoffs: an H-bond can be defined as H···A ≤ 3.0 Å, or D···A ≤ 3.5 Å or ≤ 4.1 Å, with variable angle constraints. "There is usually no consensus on the empirical thresholds", so ProLIF exposes them as parameters — [ProLIF, J Cheminform 2021](https://link.springer.com/article/10.1186/s13321-021-00548-6)

### Inferences
Recommended metric set for the RNA benchmark:
- **Units.** Compare at two granularities: atom-pair level (ligand atom, RNA atom, type) and residue level (nucleotide, type). Residue level is what PLIP/fingeRNAt/ProLIF fingerprints are built on; atom level is needed for LW-edge assignment.
- **Per-type precision, recall and F1** against (a) the expert-curated reference subset and (b) the multi-tool majority consensus. Report both, because consensus inherits shared biases, for example common OpenBabel typing in PLIP, Arpeggio and fingeRNAt.
- **Pairwise tool–tool Jaccard per type**, and Fleiss' kappa over shared detectable types. Use a harmonised type map (for example Arpeggio "ionic" vs PLIP "salt bridge" vs fingeRNAt "cation–anion").
- **Cutoff sensitivity.** Sweep each tool's main distance/angle threshold (for example H-bond D···A 3.2–4.1 Å), plot recall/precision curves, and report how many contacts flip type (H-bond ↔ ionic, stacking ↔ hydrophobic).
- **Hydrogen/protonation sensitivity.** Rerun the same structures with different hydrogenation pipelines and report the per-type Jaccard between runs (the non-determinism floor).

### Gaps
- No study found that reports Cohen's or Fleiss' kappa between interaction tools.
- No published threshold-sweep study for nucleic-acid contacts was found.

## 3. Sources of reference truth

### Takeaway
No ready-made gold standard of atom-level RNA–ligand contacts exists. Usable references are: (1) author-described contacts in the primary structure papers, as in PLIP's 30-complex suite but curated by hand for RNA; (2) the fingeRNAt 207-structure set and HARIBOSS as structure lists, not contact truth; (3) for ions, stereochemistry-validated sets (Zheng et al. 2015, 15,334 "reliable" Mg2+ sites; Auffinger criteria); (4) physics-based references such as SAPT0 fragment energies (Splinter) for the energetic plausibility of contact types, which are protein-oriented.

### Cited Findings
- **fingeRNAt dataset:** 207 non-redundant RNA–small-molecule structures (180 X-ray, 27 NMR), based on Philips et al. Where RNAs with >90% identity bind the same ligand, only the highest-resolution structure was kept. A separate redocking set has 144 complexes, taken from AnnapuRNA — [fingeRNAt](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009783)
- **HARIBOSS** (Panei, Torchet, Ménager, Gkeka, Bonomi; Bioinformatics 38:4185, 2022):
  - Curated RNA–small-molecule complexes from X-ray, NMR and cryo-EM, with a "drug-like" filter (160–1,000 Da), a 6 Å interaction radius, and buffers kept on purpose.
  - The 2022 release had 716 structures (non-redundant 484; about 83% X-ray, 11% cryo-EM, 6% NMR).
  - On 2026-10-19 the live site showed 2,131 PDB entries and 5,002 pockets; the data are CC0.
  - Sources: [HARIBOSS bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full.pdf); [HARIBOSS Bioinformatics](https://academic.oup.com/bioinformatics/article/38/17/4185/6633309); [HARIBOSS site](https://hariboss.pasteur.cloud/about-general)
- **Literature-validated interaction suite** (precedent: 30 protein–ligand complexes in which every paper-reported interaction must be detected) — [PLIP 2015](https://academic.oup.com/nar/article/43/W1/W443/2467865)
- **Mg2+ reference sites.** Zheng et al. 2015 validated Mg2+ sites in RNA crystal structures and filtered them to a benchmark of 15,334 "reliable" sites. They classified 41 inner-sphere and 95 outer-sphere coordination patterns, because many deposited Mg2+ are poorly modelled or misidentified — [Zheng et al., NAR 2015, PMC4402538](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4402538/)
- **Quantum-chemical reference.** The Splinter dataset (Spronk et al., Sci Data 2023) gives more than 1.6 million SAPT0 configurations for more than 9,000 protein/ligand fragment dimers. Each energy is split into electrostatics, exchange, induction and dispersion — [PMC10497680](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10497680/)
- First-order SAPT alone is not accurate enough to rank ligands correctly; the induction and dispersion terms are critical — [arXiv 2110.01589](https://arxiv.org/pdf/2110.01589)
- **CSD-derived criteria** are used as the basis of ion distance tables (Mg2+···O CSD mean 2.06 Å; Na+···O 2.41 Å) — [Auffinger et al., RNA 2021](https://rnajournal.cshlp.org/content/27/3/243.full); [Supplementary table](https://rnajournal.cshlp.org/content/suppl/2020/12/02/rna.076067.120.DC1/Supp_Figures_Table_Refs.docx)
- **Arpeggio's rules come from CREDO** (pdbe-arpeggio "calculates interatomic contacts ... based on the rules defined in CREDO") — [pdbe-arpeggio PyPI](https://pypi.org/project/pdbe-arpeggio)

### Inferences
- **Proposed tiered reference:**
  - Tier A: hand-curated atom-level contacts for about 20–40 riboswitch/aptamer complexes, taken from the primary papers' figures and text (author-reported H-bonds, stacking partners, metal-mediated contacts), with each claim citing the paper.
  - Tier B: multi-tool consensus on a larger stratified set (HARIBOSS non-redundant, or the fingeRNAt 207).
  - Tier C: ion sites checked against Auffinger/CheckMyMetal stereochemistry.
- Author annotations are themselves heterogeneous (different cutoffs, often no stated criteria). Tier A therefore tests "does the tool recover what experts considered meaningful" (recall-oriented). It cannot test precision for weak contacts.
- An NCI/QM check (for example SAPT0 on a small number of disputed contact geometries such as CH···O or long H-bonds) could settle borderline cases, but there are no RNA-specific SAPT reference sets.

### Gaps
- No RNA-specific SAPT or NCI-plot reference set for ligand contacts was found.
- PDBbind-NA and the NAKB ligand subsets were not checked from primary pages (the sibling note `databases_pockets.md` records the same gap).
- No existing compilation of author-annotated riboswitch–ligand contacts was found.

## 4. Best-supported geometric criteria (H-bonds, stacking, cation/anion–π, ion coordination)

### Takeaway
Tool defaults converge on permissive heavy-atom criteria: H-bond D···A 3.9–4.1 Å with angle ≥ 100°; stacking centroid ≤ 5.5 Å with offset ≤ 2.0 Å; π–ion ≤ 6.0 Å; ionic ≤ 5.5 Å. The best-supported ion criteria are the CSD/PDB-derived stereochemical distances of Auffinger and colleagues:
- Mg2+–O about 2.06–2.1 Å (1.9–2.4 Å), octahedral.
- Na+–O about 2.41 Å (2.3–2.6 Å).
- K+–O about 2.8 Å (2.6–3.2 Å), coordination number 7–8.
- 2.4–2.6 Å is ambiguous between Mg2+ and K+.

Mg2+ inner-sphere binding to nucleobase N7 or carbonyls is rare. Deposited cases are mostly misassigned Na+, K+, NH4+ or water.

### Cited Findings
- **fingeRNAt defaults:**
  - H-bond: D···A ≤ 3.9 Å, D–H···A angle 100–260°.
  - Halogen bond: ≤ 4.0 Å, angles 165°/120° ± 30°.
  - Cation–anion: ≤ 5.5 Å.
  - π–cation/π–anion: ≤ 6.0 Å, 90° ± 30°.
  - Stacking: centroid ≤ 5.5 Å, offset ≤ 2.0 Å, parallel ≤ 30°, T-shaped 90° ± 30°.
  - Ion-mediated: Mg ≤ 3.2 Å, Na ≤ 3.6 Å, K ≤ 3.9 Å, other ≤ 3.5 Å.
  - Water: ≤ 3.5 Å to the water O. Lipophilic: ≤ 4.0 Å.
  - Source: [fingeRNAt config.py](https://raw.githubusercontent.com/n-szulc/fingeRNAt/master/code/config.py)
- fingeRNAt literature basis for the criteria:
  - H-bonds: Torshin 2002; Panigrahi & Desiraju 2007.
  - Halogen bonds: Auffinger 2004.
  - Cation–π: Gallivan & Dougherty 1999.
  - π-stacking: McGaughey 1998.
  - Anion–π: Borozan 2016.
  - Metal bond lengths: Gagné & Hawthorne 2016.
  - Water-mediated: Lu 2007.
  - Lipophilic: Bissantz 2010.
  - Source: [fingeRNAt reference list](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009783)
- **PLIP defaults:**
  - Hydrophobic: ≤ 4.0 Å. H-bond: ≤ 4.1 Å, donor angle ≥ 100°.
  - π-stacking: ≤ 5.5 Å, ± 30°, offset ≤ 2.0 Å. π–cation: ≤ 6.0 Å. Salt bridge: ≤ 5.5 Å.
  - Water bridge: water–polar atom 2.5–4.1 Å, omega 71–140°, theta ≥ 100°.
  - Metal: ≤ 3.0 Å. Halogen: ≤ 4.0 Å.
  - Source: [PLIP config.py](https://raw.githubusercontent.com/pharmai/plip/master/plip/basic/config.py)
- **MD-IFP alternative criteria:**
  - H-bond: H···A ≤ 3.3 Å with angle ≥ 100°. A smaller angle is "usually chosen to make the definition less strict", which is needed for less accurate crystal structures.
  - Salt bridge: ≤ 5.0 Å, based on a protein–ligand contact analysis finding acceptor–H at 2.8–3.3 Å, which is stricter than PLIP's 5.5 Å.
  - Aromatic: distance only, with at least 5 heavy-atom contacts. This is less strict than PLIP's angle criterion, yet produced only 3 FP among 56.
  - Halogen–carboxyl: C–Hal···O > 170°.
  - Source: [Kokh et al., arXiv 2006.11066](https://arxiv.org/pdf/2006.11066)
- **H-bond definition variants:** H···A 3.0 Å, or D···A 3.5 or 4.1 Å, with no consensus — [ProLIF 2021](https://link.springer.com/article/10.1186/s13321-021-00548-6)
- **Arpeggio contact classes** are separated by atom type plus distance/angle constraints. They include H-bond, weak H-bond (CH···O class), "polar" (H-bond without the angle term), ionic, metal, aromatic, carbonyl, halogen and atom–plane classes (cation–π, donor–π, carbon–π) — [Arpeggio, J Mol Biol 2017 (PMC5282402)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5282402/); [GitHub PDBeurope/arpeggio](https://github.com/PDBeurope/arpeggio)
- **Ion distances (Auffinger, Ennifar, D'Ascenzo, RNA 2021, "Deflating the RNA Mg2+ bubble"):**
  - Mg2+ is about 2.1 Å from its ligands (favoured 2.06 Å, range 1.9–2.4 Å) and strictly octahedral.
  - K+ is about 2.8 Å (2.6–3.2 Å), with coordination number 7–8.
  - 2.4–2.6 Å fits neither ion. Mg2+ prefers phosphate/carboxylate; K+ prefers nucleobase carbonyls (for example the G•U major groove).
  - Water–O contacts < 2.4 Å flag misplaced Mg2+ hydration spheres.
  - Reliable Mg2+/Na+ assignment needs about ≤ 2.4 Å resolution. Do not assign 0.5-occupancy ions at medium/low resolution. For light ions "stereochemistry remains the only viable option".
  - In ribosome 6QNR, 853 of 1005 Mg2+ (about 85%) fail at least one criterion: 138 look like Na+, 361 like K+, and 144 have an empty first shell. About 400 are probably K+.
  - Source: [RNA 27:243 full text](https://rnajournal.cshlp.org/content/27/3/243.full)
- **Supplementary Table S1 (stringent vs relaxed):** Mg2+···O 1.9–2.3 Å stringent / 1.8–2.4 Å relaxed (CSD mean 2.06 Å); Na+···O 2.3–2.6 Å stringent / 2.4–2.6 Å relaxed (CSD mean 2.41 Å) — [Supp. Table](https://rnajournal.cshlp.org/content/suppl/2020/12/02/rna.076067.120.DC1/Supp_Figures_Table_Refs.docx). The relaxed Na+ lower bound as rendered in the search summary (2.4) is higher than the stringent one (2.3), which may be a transcription error. Verify against the table.
- **Mg2+ and nucleobase nitrogens/carbonyls.** Mg2+ assigned near imine N1/N3/N7 mostly comes from poor reading of density (Na+, K+, NH4+, water, spurious peaks), and Mg2+ does not bind N7 apart from a few exceptions (Leonarski, D'Ascenzo, Auffinger, "Mg2+ ions: do they bind to nucleobase nitrogens?", NAR). The follow-up RNA 2019 survey reports that nucleobase carbonyls are poor Mg2+ inner-sphere binders but excellent monovalent-ion binders, and proposes stereochemical validation guidelines — [RNA 25:173](https://rnajournal.cshlp.org/content/25/2/173); [PMC6348993](https://pmc.ncbi.nlm.nih.gov/articles/PMC6348993/); [IBMC publication list](https://ibmc.cnrs.fr/en/laboratoire/m3i/publications/?tgid=9002)
- **CheckMyMetal (Gucwa, Bijak, Zheng, Murzyn, Minor; IUCrJ 11:871, 2024):**
  - Uses valence-bond-theory rules. VALENCE, nVECSUM, gRMSD and VACANCY are each labelled ACCEPTABLE/BORDERLINE/DUBIOUS.
  - Neighbouring metal contacts are limited to 4 Å. Waters count in the first shell. Cryo-EM maps (MRC) are now accepted.
  - In ribosome 8B0X (361 Mg2+, 168 K+), 27 Mg2+ did not match the map or CMM score, and 8 sites modelled as Mg2+ in 6QNR/7K00 are K+ in 8B0X.
  - Cryo-EM CMM values are consistently lower than X-ray values at comparable resolution. About 40% of PDB structures contain metals, which are not always correctly identified.
  - Sources: [IUCrJ 2024](https://journals.iucr.org/m/issues/2024/05/00/be5298/index.html); [CMM 2017, PMC5349434](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5349434/)
- **MESPEUS:** a metal-site geometry database for proteins (2008; updated NAR 2023/24, 40 metal types, symmetry mates generated, monthly updates; July 2026 version with 84,772 PDB IDs). It is a source of empirical metal–donor distance distributions, though protein-centred — [MESPEUS](https://mespeus.nchu.edu.tw/); [Edinburgh record](https://www.research.ed.ac.uk/en/publications/mespeus-a-database-of-the-geometry-of-metal-sites-in-proteins/); [metal-distance overview, PMC11066882](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11066882/)

### Inferences
- **Inner-sphere ion criteria.** Use Auffinger's ion-specific windows (Mg2+ ≤ 2.4 Å, Na+ 2.3–2.6 Å, K+ 2.6–3.2 Å) instead of a single "metal ≤ 3.0 Å" (PLIP) or the per-ion fingeRNAt bridge cutoffs (Mg 3.2 Å etc.). Then classify the contact as inner-sphere vs outer-sphere (water-mediated, about 4–4.5 Å). fingeRNAt's 3.2 Å Mg cutoff will admit misassigned Na+/K+ sites.
- **Flag, do not drop.** Mg2+ contacts to N7 or base carbonyls, and any ion at 2.4–2.6 Å, should be kept as "ion identity questionable" observations. This fits the project's rule to keep flagged observations.
- **H-bonds.** Report heavy-atom D···A in two bands (≤ 3.5 Å "strong", 3.5–4.1 Å "permissive") so tools with different defaults can be compared in the overlapping range. Report CH···O (Arpeggio weak H-bond) as its own class; most other tools do not detect it.
- **Stacking.** Fix a common parallel/T-shaped definition (centroid ≤ 5.5 Å, offset ≤ 2.0 Å, interplanar angle ≤ 30° or 90° ± 30°) and, in sensitivity analysis, compare it with the distance-only (≥ 5 atoms) definition.

### Gaps
- No RNA-specific CSD survey of base–ligand stacking geometry or of nucleobase H-bond distance distributions was retrieved. Classical sources (for example Panigrahi & Desiraju 2007; McGaughey 1998) are cited by fingeRNAt but were not read.
- Anion–π criteria (Borozan 2016) were not read directly.
- The exact Leonarski N7 paper citation (NAR, 2016/2017) was not confirmed from a primary page: search results gave conflicting years.
- Arpeggio's numeric cutoffs (CREDO) were not retrieved.

## 5. Handling hydrogens, protonation and alternate locations

### Takeaway
Every candidate tool depends on hydrogens, whether added internally (PLIP and Arpeggio via OpenBabel; BINANA from atom names) or supplied by the user (ProLIF and the fingeRNAt receptor). No method exists for hydrogenating or protonating RNA–ligand complexes specifically. Kokh et al. trace much of the disagreement in water bridges and H-bonds to protonation and hydrogen orientation. A benchmark therefore has to fix one hydrogenation pipeline for all tools and measure sensitivity to it separately.

### Cited Findings
- **fingeRNAt.** The receptor must have explicit hydrogens. Ligand H are added by OpenBabel (default) or RDKit; ligand protonation and formal charges are the user's job. In the paper, ligand protonation was set from literature and Chemicalize pKa, and RNA and ligand were protonated in isolation ("an approximation"). OP1/OP2 are always anionic. Hydrogen counts differ between OpenBabel, RDKit, PyMOL and Chimera — [fingeRNAt GitHub](https://github.com/n-szulc/fingeRNAt); [fingeRNAt paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009783); [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9197077/fullTextXML)
- **PLIP** needs no structure preparation: it hydrogenates internally and uses OpenBabel for molecule representation — [PLIP 2015](https://academic.oup.com/nar/article/43/W1/W443/2467865)
- **ProLIF** requires explicit hydrogens, so X-ray PDB files need preprocessing — [bioRxiv 2023.11.13.566916](https://www.biorxiv.org/content/10.1101/2023.11.13.566916v1.full.pdf); [ProLIF 2021](https://link.springer.com/article/10.1186/s13321-021-00548-6)
- **Protonation pipeline precedent (Errington et al.):** PDB2PQR for protein H, RDKit for ligand H, then a short MMFF minimisation of the ligand hydrogens with heavy atoms fixed, to optimise the H-bond network consistently. They note that Reduce, the CSD API and Schrödinger PrepWizard are alternatives — [arXiv 2409.20227](https://arxiv.org/pdf/2409.20227)
- **Reduce** adds H in standardised geometry and optimises OH, NH3+ and amide flips. It supports nucleic acids and recognises backbone connections — [Reduce README](https://genomics-lab.fleming.gr/fleming/tools/molprobity/modules/reduce/README.usingReduce.txt). The cctbx `mmtbx.reduce2` / `hydrogenate` alternatives exist, but even the maintainer was unsure how they differ — [cctbxbb mailing list](https://phenix-online.org/archives/list/cctbxbb@phenix-online.org/message/2QL7BA47KFXZNN4NB4ZUQZS52N3JSTCH/)
- **Hydride** (Kunzmann, Anter, Hamacher, Algorithms Mol Biol 2022) places H by superimposing fragments from a library. It works for organic compounds and biomacromolecules, needs no molecule-specific parameters, and is open source — [PMC8966362](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8966362/); [hydride PyPI](https://pypi.org/project/hydride/1.2.1)
- **Altlocs.** Kokh et al. kept conformation A only — [arXiv 2006.11066](https://arxiv.org/pdf/2006.11066). pdbe-arpeggio does not handle altlocs or missing density and "may give anomalous results" — [GitHub PDBeurope/arpeggio](https://github.com/PDBeurope/arpeggio)
- **Charge-dependent typing.** Protonation decides whether an interaction is labelled ionic or H-bond — [Errington et al.](https://arxiv.org/pdf/2409.20227)

### Inferences
Recommended protocol:
- **(i) Common input.** Feed every tool the same hydrogenated mmCIF/PDB. Add RNA H with Reduce (or Hydride). Add ligand H from CCD ideal coordinates with the CCD formal charges, keeping the CCD state as the deposited reference. Use an optional pH 7 variant from a pKa tool for known titratable ligands (aminoglycosides, for example).
- **(ii) Sensitivity arm.** Rerun with each tool's native hydrogenation (PLIP/Arpeggio OpenBabel) and with Hydride vs Reduce. Report the per-type Jaccard against the common-input run.
- **(iii) Altlocs.** Evaluate altloc A (or highest occupancy) as the primary run, with all-altloc enumeration as a flagged secondary analysis. Record ligand altloc presence as a stratum.
- **(iv) Charged groups.** Treat RNA phosphate OP1/OP2 as anionic in all tools. Check that Arpeggio's protein-oriented typing does this; the sibling note `fingerprints_contact_typing.md` flags this as a needed spot check.
- **Keep H-bond and ionic together.** Because the H-bond vs ionic label depends on protonation, compare tools first on a merged "polar/electrostatic" class, then on the split classes.

### Gaps
- No study found that quantifies RNA-specific hydrogenation differences (Reduce vs Hydride vs OpenBabel vs PDB2PQR) on contact detection.
- No source validates CCD ideal charges against bound-state protonation for RNA ligands.

## 6. Sample sizes and stratification

### Takeaway
Published validation sets are small: 30 complexes (PLIP), 40 (Kokh, from FLIP's 50) and 3 for the HSP90 sub-comparison. RNA tool papers use the fingeRNAt 207 set, the AnnapuRNA 144 redocking set, and HARIBOSS (484 non-redundant structures / 676 pockets in 2022). No formal power analysis or stratification guidance for interaction-tool benchmarks was found. The evidence that does exist identifies resolution, experimental method (cryo-EM lowers metal-site reliability) and ion presence as the main confounders.

### Cited Findings
- PLIP's validation suite covered resolutions 1.2–3.3 Å — [PLIP 2015](https://academic.oup.com/nar/article/43/W1/W443/2467865)
- Kokh: 40 complexes, 250 non-hydrophobic interactions and 69 water bridges — [arXiv 2006.11066](https://arxiv.org/pdf/2006.11066)
- fingeRNAt: 207 structures (180 X-ray, 27 NMR; NMR structures have no waters); the authors also note that low resolution can hide interactions — [fingeRNAt](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009783); [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9197077/fullTextXML)
- HARIBOSS 2022: 716 structures (484 non-redundant), about 83% X-ray, 11% cryo-EM, 6% NMR. Redundancy was reduced by CD-HIT-EST 90% clustering plus pocket RMSD 1.5 Å; no resolution cutoff was applied (typical 2.8 Å X-ray, 3.0 Å cryo-EM) — [HARIBOSS Bioinformatics](https://academic.oup.com/bioinformatics/article/38/17/4185/6633309)
- Ion-assignment reliability depends on resolution: about ≤ 2.4 Å is needed for Mg2+/Na+ — [Auffinger et al. 2021](https://rnajournal.cshlp.org/content/27/3/243.full)
- The number of detectable metal sites drops with resolution, and cryo-EM metal-site scores are lower than X-ray at the same resolution — [CheckMyMetal 2024](https://journals.iucr.org/m/issues/2024/05/00/be5298/index.html)
- fpocketR selected ≤ 200 nt, ≤ 3.5 Å and ligand QED ≥ 0.3 from HARIBOSS (48 complexes), a precedent for quality-filtered RNA subsets — [fpocketR, PNAS 2025](https://doi.org/10.1073/pnas.2422346122) (via sibling note `databases_pockets.md`)

### Inferences
Suggested design (inference, not from a source):
- **Tier A:** about 30–40 hand-curated complexes, matching the PLIP and Kokh precedents. Cover riboswitch classes in the project (purine RF00167, SAM-I RF00162, TPP RF00059), aminoglycoside–rRNA A-site, and synthetic aptamers.
- **Tier B:** the HARIBOSS non-redundant set (about 500 structures).
- **Strata:**
  - Method: X-ray, cryo-EM, NMR.
  - Resolution bins: ≤ 2.0, 2.0–2.5, 2.5–3.0, > 3.0 Å. The 2.4 Å boundary matters for ions.
  - Ligand class: purine-like/planar, aminoglycoside/polycationic, phosphate- or carboxylate-bearing (Mg2+-mediated), and fragments or additives.
  - Presence of modelled ions and waters.
  - Ligand altlocs or partial occupancy.
- **Precision of estimates.** Report bootstrap CIs per stratum. With 30–40 complexes and roughly 5–10 polar contacts each, per-type recall CIs will be wide for rare types (halogen, anion–π). Pool rare types or state their n explicitly, following the project's numerator/denominator rule.
- **Cluster-aware sampling.** Cluster by Rfam family and ligand so that the many purine-riboswitch structures do not dominate. Report per-family as well as pooled results.

### Gaps
- No published recommendation on sample size or statistical power for interaction-tool benchmarks was found.
- No study stratifies tool agreement by resolution or by experimental method.
