# RNA–small molecule interaction fingerprints and contact classification for 3D structures

Context for the reader: the project already prototypes ECIFP (PDBe "Extended Connectivity Interaction Fingerprints", Kunnakkattu et al., bioRxiv 2026, doi 10.64898/2026.04.17.719144: PDBe Arpeggio contacts, Morgan radius-2 1024-bit fingerprints of the contacting atoms on the ligand and polymer sides, similarity = Tanimoto_ligand^0.4 × Tanimoto_polymer^0.6) and a Leontis–Westhof (LW) edge classification of contacts. These notes survey what other tools do. Research date 2026-10-23. PMC article pages were blocked by a browser check, so details for fingeRNAt came from Europe PMC full text and the GitHub repo, and details for PLIP from its repo and config. The ECIFP preprint itself did not come up in web search (see Gaps).

## Q1. Which RNA-capable interaction fingerprints exist, what do they encode, and how is similarity computed?

### Takeaway
fingeRNAt (2022) is the only fingerprint tool built for nucleic acids. It is per-residue: 12 binary interaction bits in FULL mode, or 3 bits for phosphate / base / sugar in PBS mode. It covers Mg2+/K+/Na+/other ion-mediated and water-mediated contacts, pi-cation, pi-anion and stacking, and offers Tanimoto, Tversky and six other metrics. PLIP (2021) and ProLIF accept DNA/RNA receptors, but their atom typing and thresholds are the protein-ligand ones. Arpeggio (and so ECIFP) is officially protein-only. RNAmigos uses graph representations of the binding site rather than contact fingerprints. Of everything found, only RLEC, an RNA adaptation of PLEC, is atom-environment-based like ECIFP.

### Cited Findings
**fingeRNAt (Szulc, Mackiewicz, Bujnicki, Stefaniak; PLoS Comput Biol 2022, 18(6):e1009783, doi 10.1371/journal.pcbi.1009783)**
- Detects "nine types of interactions" in nucleic acid–ligand complexes, encoded as a binary Structural Interaction Fingerprint (SIFt). It can be extended with a plugin system for custom interactions — [PMC9197077](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9197077/)
- FULL fingerprint: 12 values per residue. These are H-bond, halogen bond, cation–anion, pi–cation, pi–anion, pi-stacking (split into sandwich/displaced and T-shaped in debug output), ion-mediated as four bits (Mg2+, K+, Na+, other ions), water-mediated (only with `-h2o`) and lipophilic — [GitHub n-szulc/fingeRNAt](https://github.com/n-szulc/fingeRNAt)
- SIMPLE: 1 bit per residue for any residue–ligand heavy-atom distance within 4.0 Å. PBS: 3 bits per residue (Phosphate, Base, Sugar), same 4.0 Å cutoff, canonical residues only. The fingerprint is per residue. `-detail` output lists atom IDs, coordinates and distances for each interaction — [GitHub](https://github.com/n-szulc/fingeRNAt)
- Similarity tool fingerDISt offers Tanimoto, cosine, Manhattan, Euclidean, squared Euclidean, half-squared Euclidean, Soergel and Tversky. Tanimoto works on binary SIFts only. Tversky is asymmetric, with alpha/beta hard-coded — [GitHub](https://github.com/n-szulc/fingeRNAt)
- Input: receptor as a single-model PDB with hydrogens already added. Ligands as SDF with correct formal charges and protonation; missing ligand hydrogens are added with OpenBabel or RDKit. Waters and ions in the receptor file are allowed. Python 3, CLI plus GUI, PyMOL plugin; licence GPL-3.0 — [GitHub](https://github.com/n-szulc/fingeRNAt)
- Paper analysis: 207 PDB RNA–ligand complexes (27 NMR, 180 X-ray, 0.61–4.50 Å). Pi–anion contacts: about 82% involve nucleobase rings with RNA as the anion acceptor; in about 18% the RNA phosphate is the anion — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9197077/fullTextXML)

**PLIP 2021 (Adasme et al., NAR 49:W530–W534, 2021, doi 10.1093/nar/gkab294; current citation Schake, Bolz et al. NAR 2025, doi 10.1093/nar/gkaf361)**
- Extended to DNA and RNA. Interaction types: H-bonds, hydrophobic contacts, pi-stacking, pi-cation, salt bridges, water bridges, metal complexes, halogen bonds — [PLIP 2021, PMC8262720](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8262720/)
- `--dnareceptor` treats DNA/RNA as part of the receptor; by default nucleic acid is treated as the ligand. Nucleotides are recognised by residue name (`RNA = ['U','A','C','G']`, `DNA = ['DT','DA','DC','DG']`), so modified nucleotides are not recognised as RNA by name — [PLIP DOCUMENTATION.md](https://raw.githubusercontent.com/pharmai/plip/master/DOCUMENTATION.md); [plip/basic/config.py](https://raw.githubusercontent.com/pharmai/plip/master/plip/basic/config.py)
- Output is XML/text/images/PyMOL sessions; the README lists no fingerprint output. Licence GPL-2.0 — [GitHub pharmai/plip](https://github.com/pharmai/plip)

**Arpeggio / pdbe-arpeggio (Jubb et al., J Mol Biol 2017, 429:365–371, doi 10.1016/j.jmb.2016.12.004)**
- Contact types (pdbe-arpeggio). Atom–atom: clash, covalent, vdW clash, vdW, proximal (within 5 Å, beyond vdW), H-bond, weak H-bond, halogen bond, ionic, metal, aromatic, hydrophobic, carbonyl, polar, weak polar. Atom–plane: carbon–pi, cation–pi, donor–pi, halogen–pi, sulphur–pi. Plane–plane (Chakrabarti & Bhattacharyya 2007 classes). Group: amide–amide, amide–ring. Output is a JSON contact list with `type` atom-atom / atom-plane / plane-plane / group-group and `interacting_entities` (including water categories) — [GitHub PDBeurope/arpeggio](https://github.com/PDBeurope/arpeggio)
- Built on CREDO's SMARTS atom typing with modified definitions, for example any covalently bound oxygen may act as an H-bond acceptor — [Jubb 2017, PMC5282402](https://pmc.ncbi.nlm.nih.gov/articles/PMC5282402)
- The pdbe-arpeggio README says it supports "only protein structures in mmCIF format" and does not mention nucleic acids. Requires OpenBabel ≥ 3.0; alt-locs and missing density are not handled. GPL-3.0 — [GitHub PDBeurope/arpeggio](https://github.com/PDBeurope/arpeggio)

**ProLIF (Bouysset & Fiorucci, J Cheminform 2021, doi 10.1186/s13321-021-00548-6)**
- Generates interaction fingerprints for complexes of ligands, proteins, DNA or RNA from MD, docking and experimental structures. Interactions are SMARTS patterns with distance/angle/dihedral constraints, giving a binary vector per interacting residue pair (HBDonor, HBAcceptor, PiStacking, etc.) — [ProLIF docs v1.0.0](https://prolif.readthedocs.io/en/v1.0.0/); [PMC8466659](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8466659/)
- Older docs (v0.3.1) mentioned DNA–ligand support only; RNA appears in newer docs — [ProLIF v0.3.1 docs](https://prolif.readthedocs.io/en/v0.3.1/)

**ODDT PLEC/SPLIF and RLEC**
- An ODDT developer describes PLIF as per-residue. SPLIF has only a depth-1 environment, and PLEC adds depth 0. Nothing found on RNA receptors — [ODDT Google group](https://groups.google.com/g/oddt/c/Xsq8U-2xC3Y)
- RLEC (PyPI package `rlec`) adapts PLEC (Wójcikowski et al., Bioinformatics 2019) to RNA–ligand systems. It pairs Morgan-style environments of the RNA atom and the ligand atom across increasing depths and hashes the pairs into a count vector — [search-result summary of PyPI rlec](https://pypi.org/project/rlec/0.3.2/) (the page itself did not render; details unverified)

**RNAmigos / RNAmigos2 (Waldispühl lab)**
- RNAmigos takes coordinates or a base-pairing network. Nodes carry nucleotide identity (A/U/C/G/N); edges carry one of 13 labels from the Leontis–Westhof nomenclature — [RNAmigos help](https://rnamigos.cs.mcgill.ca/help)
- Original finding: only binding-site graphs augmented with non-canonical base pairs recovered a significant signal for ligand binding — [RNAmigos help / paper summary](https://rnamigos.cs.mcgill.ca/help); [bioRxiv 701326](https://www.biorxiv.org/content/10.1101/701326v3)
- RNAmigos2 ranks candidate ligands and passes the top ones to rDock. It is trained on 1,740 native PDB RNA–small-molecule complexes and pre-trained on about 200k binding-site-like subgraphs; code at github.com/cgoliver/rnamigos2 — [Waldispühl, Benasque 2026 slides](https://www.benasque.org/2026rna/talks_contr/099_Waldispuhl_Benasque_July2026.pdf)

**HARIBOSS (Panei, Torchet, Menager, Gkeka, Bonomi; bioRxiv 2022)**
- Curated database of RNA–small molecule structures (X-ray, NMR, cryo-EM) with physicochemical annotations of ligand and pocket and interaction details. Binding site = residues within 6 Å of ligand atoms — [HARIBOSS bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full.pdf); [arXiv 2310.18985](https://arxiv.org/pdf/2310.18985)
- One 2023 study used 746 HARIBOSS complexes with 1,268 binding sites — [arXiv 2310.18985 via search summary](https://arxiv.org/pdf/2310.18985)

### Inferences
- fingeRNAt is the only tool with RNA-specific bits: separate Mg2+/K+/Na+ bridges and PBS moiety bits. Its PBS mode is a coarse version of the project's ribose / phosphate / base split. It is not chemistry-aware, so two different ligands that hit the same residues get identical SIFts. ECIFP/RLEC-style atom-environment hashing captures what fingeRNAt does not, and vice versa.
- None of these tools uses Rfam or family coordinates. Per-residue SIFts are indexed by residue number within one structure. Mapping their bits onto Rfam consensus positions (which the project already does for Arpeggio contacts) would be a new contribution.
- No tool combines LW edge-level ligand contact typing with a fingerprint. The project's edge classification has no direct precedent among these tools.

### Gaps
- No tool called "SIFt-R" was found under that name. The RNA SIFt in the literature is fingeRNAt's.
- The ECIFP preprint (Kunnakkattu et al. 2026) was not indexed by the search engine used, so its RNA handling could not be checked independently.
- RLEC authorship, paper, atom typing and benchmarks could not be verified (PyPI page did not render).
- The exact ProLIF version that added RNA, and whether its SMARTS cover nucleotide-specific groups such as phosphate charge, were not confirmed.

## Q2. How do these tools define nucleotide atoms/edges for ligand contacts? Is Leontis–Westhof used for ligand recognition?

### Takeaway
The tools group nucleotide atoms by moiety at most. fingeRNAt PBS uses fixed atom-name lists for phosphate, base and sugar. PLIP, ProLIF and Arpeggio use generic SMARTS/OpenBabel typing. None classifies ligand contacts by LW edge. LW labels appear only as base-pair edges in RNAmigos graphs and in manual structural surveys: Seelam et al. 2019 found the Watson–Crick edge dominates riboswitch base:ligand H-bonded pairs.

### Cited Findings
- fingeRNAt groups (config.py): Phosphate = P, OP1, OP2, OP3. Base = C2, C4, C5, C6, C8, N1, N2, N3, N4, N6, N7, N9, O2, O4, O6, C7. Sugar = C1', C2', C3', C4', C5', O2', O3', O4', O5'. Note that O3' and O5' are in the sugar group, not the phosphate group — [fingeRNAt code/config.py](https://raw.githubusercontent.com/n-szulc/fingeRNAt/master/code/config.py)
- fingeRNAt modified residues: a non-canonical residue name with canonical atom names works with FULL and SIMPLE, but the ACUG/PuPy wrappers skip it. Non-canonical atom names break PBS. Duplicate residue numbers can make pi results unreliable — [GitHub](https://github.com/n-szulc/fingeRNAt)
- PLIP finds aromatic rings generically (OpenBabel SSSR plus aromaticity, planarity check as fallback), not nucleobase-specifically. The documentation lists phosphate as a negative charge for ligands, but does not say how receptor nucleic acid phosphates are charged — [PLIP DOCUMENTATION.md](https://raw.githubusercontent.com/pharmai/plip/master/DOCUMENTATION.md)
- RNAmigos edges use 13 LW labels for base–base pairs, which describe the binding-site context rather than the ligand contact — [RNAmigos help](https://rnamigos.cs.mcgill.ca/help)
- Seelam, Mitra, Sharma, RNA 2019, 25(10):1274–1290: 67 riboswitch aptamer crystal structures and 80 unique base:ligand pairs with at least 2 H-bonds, found by visual inspection. The WC edge is most common, then the purine sugar edge, then the uracil Hoogsteen edge. Ligand moieties fall into base:base, base:amino acid, base:sugar, base:phosphate and base:other. QM interaction energies range from about −3 (base:other) to about −48 kcal/mol (base:phosphate); 67/80 pairs were stable after QM optimisation — [RNA 2019](https://rnajournal.cshlp.org/content/25/10/1274)
- Kligun & Mandel-Gutfreund (RNA Biology 2013): ligand-binding pockets from 10 RNA groups are enriched in rare nucleobase conformations and sugar puckers compared with computed pockets across RNA structures — [RNA Biol 2013 via search summary](https://informahealthcare.com/doi/full/10.4161/rna.24682)

### Inferences
- The project's H-bond-based WC / Hoogsteen / sugar-edge assignment is consistent with Seelam et al.'s manual edge classification, which could serve as a qualitative check: WC-edge-dominant counts in riboswitches. Their per-edge counts were not available on the page.
- If fingeRNAt's PBS split is reused, note the different O3'/O5' assignment. The project's "phosphate" category may differ, so the definitions should be stated explicitly when comparing.

### Gaps
- No systematic PDB-wide quantification of ligand contacts per LW edge was found.
- Seelam et al.'s per-edge counts and H-bond criteria were not on the abstract page.

## Q3. Distance/angle criteria; known RNA issues (hydrogens, protonation, modified nucleotides)

### Takeaway
fingeRNAt and PLIP use very similar defaults: H-bond D···A 3.9 or 4.1 Å with angle ≥ 100°, stacking centroid ≤ 5.5 Å with offset ≤ 2.0 Å, pi–ion ≤ 6.0 Å, ionic ≤ 5.5 Å. fingeRNAt adds ion-specific bridge cutoffs: Mg 3.2, Na 3.6, K 3.9, other 3.5 Å. All rely on added hydrogens and assigned protonation. Both groups flag that no RNA-specific protonation or hydrogenation method exists, and that hydrogen placement makes results non-deterministic.

### Cited Findings
- fingeRNAt defaults (config.py): MAX_HB_DIST 3.9 Å; H-bond angle 100–260°; halogen ≤ 4.0 Å, angles 165° / 120° ± 30°; cation–anion ≤ 5.5 Å; pi–ion ≤ 6.0 Å at 90° ± 30°; ring–ring centroid ≤ 5.5 Å, offset ≤ 2.0 Å, parallel/displaced ≤ 30°, T-shaped 90° ± 30°; ion-mediated MAX_ION_DIST 3.9, Mg 3.2, K 3.9, Na 3.6, other 3.5 Å; water ≤ 3.5 Å to the water oxygen; lipophilic ≤ 4.0 Å; SIMPLE/PBS 4.0 Å; recognised cations MG, CA, ZN, NA, K, MN, FE, CO, NI, LI, CU, AL, PT, HG, CD, YB, SR, GD, HO, AU, RB, CS, EU, SM, PB, CE, BA — [fingeRNAt config.py](https://raw.githubusercontent.com/n-szulc/fingeRNAt/master/code/config.py)
- fingeRNAt H-bonds use the heavy-atom D···A distance by default. `-dha` adds the D–H···A angle and needs explicit hydrogens. Hydrogen counts varied between OpenBabel, RDKit, PyMOL and Chimera. The authors note that no hydrogen-adding method for RNA–ligand complexes exists, that RNA/ligand protonation was assigned in isolation, that low resolution can hide interactions, and that NMR structures lack waters — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9197077/fullTextXML)
- PLIP defaults (config.py): binding-site residues within 7.5 Å; hydrophobic ≤ 4.0 Å; H-bond ≤ 4.1 Å, donor angle ≥ 100°; pi-stacking ≤ 5.5 Å, ±30°, offset ≤ 2.0 Å; pi–cation ≤ 6.0 Å; salt bridge ≤ 5.5 Å; water bridge water–polar atom 2.5–4.1 Å, omega 71–140°, theta ≥ 100°; metal ≤ 3.0 Å; halogen ≤ 4.0 Å, 165° / 120° ± 30° — [PLIP config.py](https://raw.githubusercontent.com/pharmai/plip/master/plip/basic/config.py)
- PLIP: OpenBabel hydrogen placement can be non-deterministic, so repeated runs may give different interaction sets. The advice is to protonate once beforehand or use `--nohydro`. Only the first NMR model is used by default — [GitHub pharmai/plip](https://github.com/pharmai/plip)
- pdbe-arpeggio does nothing beyond BioPython defaults. Alt-locs and missing density are not handled and may give anomalous results — [GitHub PDBeurope/arpeggio](https://github.com/PDBeurope/arpeggio)

### Inferences
- The project's Mg2+ contacts would be most comparable to fingeRNAt if they adopted its per-ion cutoffs (Mg 3.2 Å), which approximate inner-sphere coordination.
- Arpeggio-based ECIFP inherits Arpeggio's protein-only atom typing. RNA atom types (for example phosphate OP1/OP2 as ionic, N7 as acceptor) should be spot-checked in the PDBe Arpeggio output for the project's families.
- PLIP's name-based nucleotide detection drops modified nucleotides. fingeRNAt keeps them only when atom names are canonical. The project's mmCIF / parent-CCD route is more robust than either.

### Gaps
- fingeRNAt's donor/acceptor SMARTS typing and its full geometric table (supplementary) were not retrieved.
- Arpeggio's numeric per-type cutoffs were not retrieved (they are in the paper/source, not the README).

## Q4. Benchmarks: do RNA fingerprints separate binding modes or ligand classes, or work for virtual screening?

### Takeaway
The evidence is mostly from docking and virtual screening, not from separating native PDB binding modes. fingeRNAt SIFts clustered docked HIV-1 TAR compounds with partial active enrichment, and correlated only weakly with RMSD. SIFt plus ML outperformed general scoring functions on a six-target benchmark (Brief Bioinform 2023). No published benchmark of fingerprint-based clustering of native ligands within an Rfam family was found.

### Cited Findings
- fingeRNAt clustering: PCA plus k-means on 1,508 docked compounds against HIV-1 TAR (30 active, 1,478 inactive). Fifteen clusters with all interactions gave silhouette 0.516, with some clusters enriched for actives and two containing none. Removing lipophilic interactions gave 9 clusters, silhouette 0.72 and clearer active/inactive separation — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9197077/fullTextXML)
- SIFt similarity vs RMSD: RNA-Puzzles round 23 R² 0.31–0.61. Across 144 redocked complexes the mean R² was 0.417 (Tanimoto) and 0.380 (Tversky), with Spearman about −0.52 and −0.49. In one example a lower-RMSD pose recovered only 6/17 reference interactions while a higher-RMSD pose recovered 13/17 — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9197077/fullTextXML)
- Szulc, Mackiewicz, Bujnicki, Stefaniak, "Structural interaction fingerprints and machine learning for predicting and explaining binding of small molecule ligands to RNA", Brief Bioinform 24(4) bbad187 (2023). SIFt-based ML models outperform classic general-purpose scoring functions in virtual screening. A benchmark of ligands for six RNA targets is released (Zenodo; GitHub filipsPL/fingernat-ml, CC0). Explainable AI is used to interpret recognition — [bioRxiv 2023.01.11.523582](https://www.biorxiv.org/content/10.1101/2023.01.11.523582.full.pdf); [Zenodo 7486222](https://zenodo.org/records/7486222); [DOI bbad187](https://unpaywall.org/10.1093%2FBIB%2FBBAD187)
- A separate nucleic-acid benchmark of ligand-based fingerprints, similarity measures and docking found performance depends strongly on descriptor, similarity measure and target; a consensus method performed best — [real.mtak.hu bbaf620](https://real.mtak.hu/247422/1/bbaf620.pdf) (search summary only; not read in full)
- RNAmigos: graph binding-site encodings with non-canonical LW edges were the only representation carrying a significant ligand-binding signal — [RNAmigos](https://rnamigos.cs.mcgill.ca/help)

### Inferences
- The fingeRNAt finding that removing lipophilic bits improves separation suggests that, in the project's prototype, low-specificity contacts (vdW/proximal in Arpeggio) may dilute similarity. Weighting or excluding them is worth testing.
- Because SIFt similarity and RMSD disagree, an interaction-based comparison adds information beyond the project's existing structural clustering (RMSD/TM-score). That supports showing ECIFP/edge similarity as a separate axis.
- Candidates to adopt: (1) per-residue moiety bits (P/B/S, or the project's WC/H/S edge + ribose + phosphate), keyed by Rfam consensus position, so fingerprints are comparable across chains in a family; (2) separate Mg2+/K+/Na+ ion- and water-mediated bits with fingeRNAt's cutoffs; (3) Tversky similarity, which is asymmetric (how many reference interactions does this structure recover); (4) checking within-family ligand class separation, for example guanine vs adenine vs 2-aminopurine in RF00167, as a validation that has no published precedent.

### Gaps
- The final numeric results (AUC/enrichment) and the six targets of the Brief Bioinform 2023 paper were not retrieved (bioRxiv returned HTTP 429).
- No study evaluating ECIFP-like (atom-environment) interaction fingerprints on RNA binding-mode separation was found. RLEC's benchmarks could not be verified.
