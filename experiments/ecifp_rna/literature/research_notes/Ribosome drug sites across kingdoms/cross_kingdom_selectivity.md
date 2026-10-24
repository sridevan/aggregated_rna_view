# Structural basis of species selectivity of ribosome-targeting small molecules across kingdoms

Scope note: research conducted 2026-10-24 via web search/fetch. Several primary full texts were
not retrievable (PMC browser-verification walls, bioRxiv 429, Nature IDP redirect), so some claims
below rest on abstracts, PDB entry pages, or secondary summaries — flagged where that is the case.
"Established" vs "speculative" is marked per item.

## Q1. Comparative structural studies of the same inhibitor on ribosomes from different kingdoms

### Takeaway
A respectable set of single-drug or single-organism structural studies now exists (yeast 80S
inhibitor panel, aminoglycosides on 80S, *Plasmodium*, *Leishmania*, *Trypanosoma*, *Candida*,
human mitoribosome), but genuinely *paired* cross-kingdom studies in which the same compound is
solved on ribosomes from two or more kingdoms are rare — the tigecycline study (human 55S/39S +
human 80S + yeast 80S in one paper) is the clearest example, and the only one found that spans
mitochondrial and two cytosolic ribosomes in a single work.

### Cited Findings
- **Yeast 80S inhibitor panel (established).** Garreau de Loubresse, Prokhorova, Holtkamp, Rodnina,
  Yusupova & Yusupov, "Structural basis for the inhibition of the eukaryotic ribosome", *Nature*
  513:517–522 (2014), DOI [10.1038/nature13737](https://doi.org/10.1038/nature13737). 16
  high-resolution X-ray structures of *S. cerevisiae* 80S with 12 eukaryote-specific and 4
  broad-spectrum inhibitors; all inhibitors bound mRNA/tRNA-binding sites. Example PDB entries:
  cycloheximide [4U3U](https://rcsb.org/structure/4u3u), lactimidomycin
  [4U4R](https://rcsb.org/structure/4u4r), anisomycin [4U3M](https://www.rcsb.org/structure/4U3M),
  edeine [4U4N](https://rcsb.org/structure/4u4n). The paper's own framing: structural information on
  ribosome-targeting antibiotics in bacteria is extensive whereas "our understanding of the binding
  mode of ribosome inhibitors in eukaryotes is currently limited" — [RCSB 4U4Q entry page](https://www.rcsb.org/structure/4U4Q).
- **Aminoglycosides on the eukaryotic ribosome (established).** Prokhorova et al., "Aminoglycoside
  interactions and impacts on the eukaryotic ribosome", *PNAS* 114(51):E10899–E10908 (2017), DOI
  [10.1073/pnas.1715501114](https://doi.org/10.1073/pnas.1715501114), PMID 29208708. X-ray
  structures of yeast 80S with paromomycin, geneticin (G418), gentamicin and TC007 at 3.3–3.7 Å,
  plus smFRET. Key result: **multiple** aminoglycoside-binding sites in both subunits (not just the
  decoding centre), and the ring-I 6′-hydroxyl is a determinant of binding at the canonical decoding
  centre. PDB: [5NDV](https://rcsb.org/structure/5ndv) (paromomycin, 3.30 Å),
  [5NDW](https://rcsb.org/structure/5NDW) (TC007, 3.70 Å), also 5NDG, 5NDJ, 5NDK,
  [5OBM](https://rcsb.org/structure/5obm). *Caveat:* PDBj and RCSB pages disagree on some
  resolutions (5OBM 3.4 Å vs other values; 5NDK 2.95 Å) — verify before quoting.
- **Aminoglycosides on the human 80S (established).** Tomono et al. (2024, *J. Biochem.*) cryo-EM
  structures of arbekacin- and dibekacin-bound **human** 80S ribosomes at ~2.6 Å:
  [8IFD](https://rcsb.org/structure/8IFD), [8IFE](https://rcsb.org/structure/8IFE). (Found via RCSB
  entry pages; full text not retrieved.)
- **Tigecycline — the clearest cross-kingdom paired study (established).** Li et al., "Structural
  basis for differential inhibition of eukaryotic ribosomes by tigecycline", *Nat Commun* 15 (2024),
  DOI [10.1038/s41467-024-49797-7](https://doi.org/10.1038/s41467-024-49797-7);
  [DOAJ record](https://doaj.org/article/adc1134f2314428ab0c1e0d0a7f3d6fb). Cryo-EM of
  tigecycline-bound **human 55S and 39S mitoribosomes, human cytoplasmic 80S and yeast cytoplasmic
  80S** in one paper. At clinically relevant concentrations it targets the 55S mitoribosome
  (hindering A-site tRNA accommodation, blocking the PTC) and "does not bind to human 80S ribosomes
  under physiological concentrations"; at high concentrations both human and yeast 80S bind it at a
  different, conserved site restricting L1-stalk movement. Related PDB entries from this work
  include [8XSX](https://rcsb.org/structure/8xsx), 8XT0, 8XT2, 8K2A, 8K2B, EMD-38639.
- **Tigecycline follow-up (established, 2025).** Shao et al., *Nat Commun* (2025): 5–10 µM
  tigecycline inhibits mitochondrial but not cytosolic translation in T cells; cryo-EM showed three
  binding sites on T-cell mitoribosomes — [PMC12045974](https://pmc.ncbi.nlm.nih.gov/articles/12045974).
- **Linezolid on the human mitoribosome (established, 2025).** Bibel et al., "Context-specific
  inhibition of mitochondrial ribosomes by phenicol and oxazolidinone antibiotics", *Nucleic Acids
  Research* (2025), [PMID 39907106](https://pubmed.ncbi.nlm.nih.gov/39907106/); preprint
  [bioRxiv 2024.08.21.609012](https://www.biorxiv.org/content/10.1101/2024.08.21.609012.full.pdf).
  High-resolution cryo-EM of linezolid at the PTC of the human mitoribosome, PDB
  [9CN3](https://rcsb.org/structure/9cn3) (human 39S + linezolid, 2.62 Å). Binding mode shows
  "extensive similarity to the mode of bacterial inhibition" and the authors say it "suggests
  potential avenues for altering selectivity". Mitoribosome profiling: chloramphenicol and linezolid
  stall mitoribosomes mainly when Ala/Ser/Thr is in the penultimate position of the nascent chain —
  i.e. context-specific, as in bacteria. The paper's introduction states that no phenicol–mitoribosome
  structure existed at the time (so no chloramphenicol–human mitoribosome structure was found).
- **Broad off-target screen on mitoribosomes (preprint, peer review open).** "Determining the
  off-target activity of antibiotics and novel translation initiation sites in mitochondria", eLife
  reviewed preprint [103699](https://elifesciences.org/reviewed-preprints/103699v2): profiled several
  bacterial translation inhibitors on mitoribosomes and found the mechanism of mitochondrial
  inhibition is the same as the bacterial one for all tested inhibitors **except telithromycin**.
- **Streptomycin off-target on the human mitoribosome (established, 2023).** eLife paper on the
  "structural basis of streptomycin off-target binding to human mitoribosome" (Amunts group); a
  group project summary states they "recently shown that some antibiotics, including streptomycin
  bind to human mitoribosomes" —
  [STIAS project page](https://www.stias.ac.za/publications/structural-basis-of-streptomycin-off-target-binding-to-human-mitoribosome/).
  Reported there: streptomycin hydrogen-bonds directly to **A1555 and C1556** of the mitoribosomal
  small subunit, rigidifying it and impairing translation (see Q2).
- **Macrolides on the eukaryotic ribosome (established, 2021).** "Context-specific action of
  macrolide antibiotics on the eukaryotic ribosome", *Nat Commun* 12 (2021),
  [s41467-021-23068-1](https://www.nature.com/articles/s41467-021-23068-1) — macrolide action is not
  purely bacteria-restricted; relevant to the usual "A2058 vs G2058 settles it" story.
- **Parasites — *Plasmodium* (established).** Wong et al., "Cryo-EM structure of the *Plasmodium
  falciparum* 80S ribosome bound to the anti-protozoan drug emetine", *eLife* 3:e03080 (2014),
  [eLife](https://elifesciences.org/articles/03080), [PMC4086275](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4086275/).
  3.2 Å; emetine binds the SSU E-site, overlapping the pactamycin site, blocking mRNA/tRNA
  translocation. Described as the first entire eukaryotic ribosome at atomic resolution by cryo-EM.
  Separately, Wong et al., "Mefloquine targets the *Plasmodium falciparum* 80S ribosome to inhibit
  protein synthesis", *Nat Microbiol* (2017), DOI
  [10.1038/nmicrobiol.2017.31](https://doi.org/10.1038/nmicrobiol.2017.31), PMID 28288098: 3.2 Å
  Pf80S with (+)-mefloquine at the GTPase-associated centre; mutating binding residues increased
  resistance.
- **Parasites — *Leishmania* (established).** Shalev-Benami et al., *Nat Commun* (Nov 2017), PMID
  29150609, DOI [10.1038/s41467-017-01664-4](https://doi.org/10.1038/s41467-017-01664-4);
  [DOAJ record](https://doaj.org/article/4b7b9a6e66c84f36a4b0e669a2a51d88). Cryo-EM of the
  *Leishmania* ribosome with paromomycin: SSU [6AZ1](https://www.rcsb.org/structure/6AZ1) at 2.7 Å,
  LSU [6AZ3](https://rcsb.org/structure/6AZ3) at 2.5 Å. Conclusions: paromomycin interferes with
  several aspects of **cytosolic** translation, "thus highlighting the cytosolic rather than the
  mitochondrial ribosome as the primary drug target"; the PAR-binding pocket contains "unique as well
  as conserved elements" proposed as hotspots for new therapeutics.
- **Parasites — *Trypanosoma* (ribosome structures yes, drug complexes not found).** Hashem et al.,
  *Nature* 494:385–389 (2013), DOI [10.1038/nature11872](https://doi.org/10.1038/nature11872):
  high-resolution cryo-EM of the *T. brucei* 80S, noting kinetoplastid-specific divergence is
  "extreme". Rajan et al., *Nat Commun* 14 (2023), DOI
  [10.1038/s41467-023-43263-6](https://doi.org/10.1038/s41467-023-43263-6): *T. brucei* 80S at
  ~2.47 Å ([EMD-17208](https://www.ebi.ac.uk/pdbe/entry/emdb/EMD-17208),
  [EMD-17250](https://www.ebi.ac.uk/emdb/EMD-17250)) studying a single rRNA pseudouridine. Liu et
  al., *PNAS* (2016): *T. cruzi* 60S at 2.5 Å, explicitly proposing trypanosome-specific interactions
  for structure-based drug design — [PMC5087005](https://pmc.ncbi.nlm.nih.gov/articles/PMC5087005).
  A *T. brucei* small mitoribosomal subunit with mt-IF-3 is deposited at ~3.3 Å
  ([6HIY](https://rcsb.org/experimental/6HIY)).
- **Fungal pathogens — *Candida albicans* (established).** 2022 *Science Advances* study: *C.
  albicans* 80S at 2.3 Å and inhibitor complexes better than 2.9 Å; a change in a conserved amino
  acid of ribosomal protein **eL42** explains cycloheximide resistance in *C. albicans*
  ("E-site drug specificity of the human pathogen *Candida albicans* ribosome",
  [HAL record](https://hal.archives-ouvertes.fr/hal-03813407);
  [press summary](https://www.genengnews.com/news/therapeutic-target-against-common-fungus-uncovered-by-cryo-electron-microscopy/)).
  Note this is a **protein**-borne selectivity determinant, not an rRNA one.
- **Fungal pathogens — *Candida auris* (preprint).** bioRxiv, June 2025, "Unveiling the Molecular
  Architecture of *Candida auris* Ribosome",
  [2025.06.26.661723](https://www.biorxiv.org/content/10.1101/2025.06.26.661723.full.pdf): cryo-EM of
  the vacant *C. auris* ribosome and complexes with **cycloheximide, blasticidin S and geneticin
  (G418)**, plus cell-free translation inhibition assays. Full text could not be fetched (HTTP 429),
  so binding-site details and any human comparison are unverified. An earlier *C. auris* ribosome
  architecture paper: Moscow Univ. Biol. Sci. Bull. (2023),
  [Springer](https://link.springer.com/article/10.3103/S0096392523700220), framed as groundwork for
  antifungals. A Groningen group conference abstract claims "significant structural differences in
  the ribosome of the drug-resistant pathogen *Candida auris*" —
  [RUG record](https://research.rug.nl/en/publications/structural-study-of-the-candida-auris-ribosome/);
  treat as preliminary.
- **Human 80S baseline structures (established).** Khatter, Myasnikov, Natchiar & Klaholz,
  "Structure of the human 80S ribosome", *Nature* 520:640–645 (2015), PMID
  [25901680](https://pubmed.ncbi.nlm.nih.gov/25901680/) — explicitly framed as paving the way for
  analysing antibiotic side effects. Natchiar et al., *Nature* (2017), 2.9 Å human 80S with
  chemical modifications, map [EMD-4245](https://www.ebi.ac.uk/emdb/EMD-4245) (with ligands HHT and
  HYG). Myasnikov et al., *Nat Commun* (2016), "Structure-function insights reveal the human
  ribosome as a cancer target for antibiotics".

### Inferences
- The field is organised **by organism**, not by drug: most papers solve one ribosome species with a
  small panel of ligands. The comparison across kingdoms is then done narratively, by superposing
  onto previously published structures from other papers — which means it is done ad hoc, with
  inconsistent residue numbering between 16S/18S/12S/23S/28S/16S-mito.
- Tigecycline (2024) demonstrates that a single-paper, multi-kingdom design yields a
  concentration-dependent selectivity picture (mito-specific at therapeutic dose, a second conserved
  site at high dose) that a one-ribosome-per-paper approach would miss. This is direct evidence that
  systematic cross-kingdom comparison changes conclusions.
- Two distinct classes of selectivity determinant appear: **rRNA identity** (aminoglycoside A site,
  macrolide 2058) and **ribosomal protein** identity (eL42 for cycloheximide/E-site drugs in
  *Candida*). An rRNA-consensus-position-only analysis would capture the first and miss the second —
  a scope caveat for the proposed project.

### Gaps
- No chloramphenicol–human mitoribosome structure found (the 2025 NAR paper states none existed).
- No *T. brucei* cytosolic ribosome–inhibitor complex structure found.
- No paper found that solves an aminoglycoside on bacterial, human cytosolic and human mitochondrial
  ribosomes within one study; that comparison exists only across separate papers.
- The *C. auris* preprint's actual selectivity conclusions are unverified (fetch failed).

## Q2. Known determinant nucleotides, and how comprehensively catalogued

### Takeaway
A small number of determinants are very well established — 16S **A1408** (A in bacteria/mitochondria/
chloroplasts, G in eukaryotic cytosol) with position **1491**, 23S **A2058** (A in eubacteria, G in
eukaryotes and archaea), and the mitochondrial 12S **A1555G / C1494T** ototoxicity alleles — but no
comprehensive, drug-class-by-drug-class catalogue of binding-site residue identity across kingdoms
was found; the knowledge is scattered across individual papers and no source presenting it as a
systematic table turned up.

### Cited Findings
- **A1408 / 1491 (aminoglycosides, established).** Residue 1408 is adenine in bacterial and
  mitochondrial ribosomes but guanine in cytosolic ones, and contributes to aminoglycoside
  selectivity together with residue 1491 — "4′-O-substitutions determine selectivity of
  aminoglycoside antibiotics", *Nat Commun* (2014),
  [ncomms4112](https://www.nature.com/articles/ncomms4112). (Full text not retrievable; this is from
  the search-indexed summary. Needs verification of the exact 1491 claim.)
- **A1408G confers resistance (established).** The identity of position 1408 (*E. coli* numbering) is
  adenosine in prokaryotic and guanosine in eukaryotic ribosomes; an A1408G substitution in *E. coli*
  confers resistance to the aminoglycoside subclass bearing a 6′-amino group on ring I, attributed
  by footprinting to lower drug affinity for the eukaryotic-type sequence —
  [PMC3390888](https://pmc.ncbi.nlm.nih.gov/articles/PMC3390888) (Matt et al., *PNAS* 2012,
  apramycin; PMID [22699498](https://pubmed.ncbi.nlm.nih.gov/22699498)).
- **Fitness cost / lethality at 1408 (established).** "Structural Bases for the Fitness Cost of the
  Antibiotic-Resistance and Lethal Mutations at Position 1408 of 16S rRNA" —
  [PMC6983231](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6983231/). Shows the 1408 position is
  not a free parameter: substitutions carry a functional penalty.
- **Hybrid-ribosome methodology (established).** Böttger's group replaced a central 34-nucleotide
  segment of bacterial 16S helix 44 with human-derived homologues, producing *E. coli* hybrid
  ribosomes carrying the A-site sequences of the human cytosolic and human mitochondrial ribosomes —
  [PMC3390888](https://pmc.ncbi.nlm.nih.gov/articles/PMC3390888). This is the functional counterpart
  to a structural cross-kingdom comparison.
- **Apramycin selectivity (established functionally; structural basis less direct).** Apramycin, a
  structurally unique veterinary aminoglycoside, "shows little activity toward eukaryotic ribosomes,
  including hybrid ribosomes which were genetically engineered to carry the mitoribosomal
  aminoglycoside-susceptibility A1555G allele", with little hair-cell damage or hearing loss in
  cochlear explants and guinea pigs while remaining a potent antibacterial — Matt et al., *PNAS*
  2012, [PMC3390888](https://pmc.ncbi.nlm.nih.gov/articles/PMC3390888);
  [Proteopedia 4AQY](https://proteopedia.org/w/4aqy) (apramycin–ribosome structure).
  *Important caveat:* no source found states explicitly that apramycin's eukaryotic selectivity is
  caused by position 1408 specifically — that link is an extrapolation from the general 6′-amino /
  1408 rule, and apramycin lacks the canonical ring-I 6′-NH3 arrangement.
- **A2058 (macrolides, established).** Macrolide selectivity is commonly attributed to 2058,
  adenosine in eubacteria and guanosine in eukaryotes and archaea. A2058G in *E. coli* raised
  macrolide resistance by ~2 orders of magnitude, and a guanine at the mitochondrial binding site
  accounts for tiamulin resistance — Zhang et al., "Antibiotic susceptibility of mammalian
  mitochondrial translation", *FEBS Letters* (2005), DOI
  [10.1016/j.febslet.2005.09.103](https://doi.org/10.1016/j.febslet.2005.09.103).
- **2058 is not a clean kingdom split (established, important caveat).** Higher eukaryotes carry
  guanine at 2058 in both cytoplasmic and mitochondrial rRNA; but the tapeworm *Echinococcus
  multilocularis* has **adenine at 2058 in its mitochondrial** rRNA and guanine in the cytoplasmic
  one, making its mitoribosome a clarithromycin target; and *Acanthamoeba castellanii* has adenine at
  mitochondrial 2058 and is macrolide-susceptible —
  [PMC1196245](https://pmc.ncbi.nlm.nih.gov/articles/PMC1196245/).
- **12S A1555G and C1494T (established).** Both fall in the A (aminoacyl-tRNA) site of the small
  mitoribosomal subunit — the same region bound by bacteria-targeting aminoglycosides — and are the
  primary genetic traits in aminoglycoside-induced deafness
  ([PMC2897308](https://pmc.ncbi.nlm.nih.gov/articles/PMC2897308)). Proposed mechanism: the mutant
  12S secondary structure more closely resembles the corresponding *E. coli* 16S region, creating a
  new G–C pair that can serve as an aminoglycoside binding site
  ([PMC2772318](https://pmc.ncbi.nlm.nih.gov/articles/PMC2772318) / related reviews). Hybrid-ribosome
  data: A1555G or C1494T increased aminoglycoside susceptibility of cells carrying the mitochondrial
  hybrid ribosomes **4- to 16-fold** ([PMC3390888](https://pmc.ncbi.nlm.nih.gov/articles/PMC3390888)).
- **Direct structural contact at A1555 (established, 2023).** Itoh et al. (cryo-EM) report that
  streptomycin binds directly to nucleotides **A1555 and C1556** of the mitoribosomal small subunit
  via hydrogen bonds, rigidifying the subunit and impairing translation — summarised in
  [eLife reviewed preprint 103699](https://elifesciences.org/reviewed-preprints/103699v2) figure 1
  and in the STIAS summary of the streptomycin off-target work
  ([STIAS](https://www.stias.ac.za/publications/structural-basis-of-streptomycin-off-target-binding-to-human-mitoribosome/)).
- **Dissent in the literature (noted).** A letter by Böttger to *Antimicrobial Agents and
  Chemotherapy* argues Qian & Guan's footprinting data show no clear link between aminoglycoside
  chemical structure and the observed modification changes, while agreeing with the general
  conclusion — surfaced via [PMC2634874](https://pmc.ncbi.nlm.nih.gov/articles/PMC2634874) /
  [PMC6011804](https://pmc.ncbi.nlm.nih.gov/articles/PMC6011804) result set. Exact citation not
  pinned down.
- **Protein-borne determinant (established).** eL42 residue change explains cycloheximide resistance
  in *C. albicans* ([HAL hal-03813407](https://hal.archives-ouvertes.fr/hal-03813407)).

### Inferences
- The catalogue of determinants is **drug-class anecdotal**: aminoglycosides (1408, 1491, 1555/1494),
  macrolides/pleuromutilins (2058), E-site drugs (eL42). For oxazolidinones, phenicols,
  tetracyclines/tigecycline, pactamycin/emetine-site drugs and the GTPase-associated-centre
  (mefloquine) site, no comparable published list of kingdom-discriminating positions was found.
- The *Echinococcus* and *Acanthamoeba* 2058 exceptions show single-position rules break down once
  more species are examined — exactly the failure mode a systematic, archive-wide,
  consensus-position-indexed survey would expose and quantify.
- Because determinants are reported in at least four numbering systems (*E. coli* 16S/23S, human 12S,
  yeast/human 18S/28S), cross-paper comparison is currently manual and error-prone. A single family
  coordinate system (Rfam CM match states for SSU/LSU) is a plausible, concrete contribution.

### Gaps
- No source found that tabulates binding-site nucleotide identity per drug class per kingdom
  comprehensively.
- Could not verify the *Nat Commun* 2014 (ncomms4112) claims about position 1491 from the primary
  text (Nature IDP redirect blocked the fetch).
- No determinant positions found for oxazolidinones/phenicols distinguishing human mitoribosome from
  bacteria, despite the 2025 linezolid structure; the paper says the mode is "extensively similar"
  to bacteria, which implies selectivity there is not explained by a single rRNA identity change.

## Q3. Systematic / archive-wide analyses of selectivity determinants

### Takeaway
No archive-wide analysis of ribosome drug-binding sites across kingdoms was found. The three
closest works are each limited: a 2010 NAR analysis of 65 binding sites using **bacterial and
archaeal structures only**; a 2023 NSMB study of 17 compounds on a single organism; and a 2025 Cell
Reports study of natural bacterial variation using AlphaFold3 predictions. This is the strongest
evidence for a genuine gap.

### Cited Findings
- **NAR 2010, "Structural signatures of antibiotic binding sites on the ribosome"**, *Nucleic Acids
  Research* 38(18):5982, [article](https://academic.oup.com/nar/article/38/18/5982/1065889). Method:
  65 experimentally confirmed antibiotic binding sites (ABSs) defined as rRNA atoms within **6 Å** of
  any drug atom, mapped onto apo structures, compared against ~17,600 computed pockets of similar
  size (3V solvent program). Structures used: *Thermus thermophilus* and *E. coli* 30S, *Haloarcula
  marismortui* and *Deinococcus radiodurans* 50S. Findings: ABSs show an RNA signature — more
  non-paired bases, preference for *syn* base conformation, unusual sugar pucker, high evolutionary
  conservation and high centrality in the RNA interaction network; a druggability index flagged new
  candidate pockets mainly in 23S rRNA. **Crucially: "No eukaryotic or mitochondrial ribosome
  structures were analysed as binding sites"** — eukaryotic 18S/28S sequences entered only the
  conservation scoring; mitochondrial ribosomes are not mentioned (verified by fetching the article
  text).
- **NSMB 2023, "Structural conservation of antibiotic interaction with ribosomes"**, *Nat Struct Mol
  Biol*, [article](https://www.nature.com/articles/s41594-023-01047-y). Cryo-EM structures of 17
  compounds from six antibiotic classes, compared with earlier structures. Findings: high structural
  conservation of binding mode between antibiotics sharing a scaffold, **including ordered water
  molecules**; but across classes "in many cases profound differences in the exact binding
  interactions". Single-organism (bacterial) scope.
- **Cell Reports 2025, "Extensive natural variation in bacterial ribosomal drug-binding sites"**,
  [Cell Reports](https://www.cell.com/cell-reports/fulltext/S2211-1247(25)00649-7). Surveys natural
  sequence variation in drug-binding sites; predicts structures of ribosomes with divergent binding
  sites using **AlphaFold3** and compares to wild-type *E. coli* antibiotic complexes; reports that
  even subtle changes at the direct ribosome–drug interface typically caused resistance. Scope is
  *within bacteria*, and the structural component is predicted rather than experimental.
- **PNAS 2010, "Structures of the *E. coli* ribosome with antibiotics bound near the peptidyl
  transferase center explain spectra of drug action"**,
  [PNAS](https://www.pnas.org/doi/10.1073/pnas.1007988107) — an early attempt to rationalise drug
  *spectra* structurally, but bacterial.
- **Teaching resource (not an analysis).** RCSB PDB-101 "Antibiotics and Ribosome Function"
  [page](https://pdb101.rcsb.org/learn/structural-biology-highlights/antibiotics-and-ribosome-function).

### Inferences
- The 2010 NAR paper is the only methodologically "systematic binding-site" study found, and it
  explicitly excludes eukaryotic and mitochondrial ribosomes as binding-site sources — the eukaryotic
  and mitoribosome structures that would now make it cross-kingdom (yeast 80S panel 2014, human 80S
  2015, mitoribosome complexes 2023–2025, parasite/fungal structures 2013–2025) were all published
  **after** it. A 2026 repeat of that analysis across kingdoms is therefore newly possible and has
  not been done, on the evidence found.
- The 2023 NSMB finding that scaffold determines binding mode down to ordered waters supports the
  premise that mapping contacts onto a common coordinate system is informative rather than noisy.
- An experimental-structure, cross-kingdom, consensus-position-indexed survey would be
  complementary to (not duplicative of) the AlphaFold3-based Cell Reports 2025 work.

### Gaps
- Searches found **no** review that tabulates binding-site residue identity across species per drug
  class. This is a negative result reported with moderate confidence (two targeted searches); it is
  possible such a table exists inside a review whose indexing did not surface it.
- No source found that counts ribosome–drug structures per kingdom (see Q5).

## Q4. Open questions highlighted by recent (2020–2026) reviews

### Takeaway
Recent sources repeatedly name the same unresolved problem — the near-identity of fungal/parasite and
human cytosolic ribosomes makes selectivity the bottleneck for antifungal/antiparasitic ribosome
drugs, and mitoribosome off-target binding is the bottleneck for antibacterial safety — and 2024–2025
structural papers state outright that selectivity-altering avenues are still only "suggested".
However, no single 2024–2026 review matching all the requested criteria was located.

### Cited Findings
- **Antifungal selectivity bottleneck (recent review).** A *Journal of Microbiology* review on
  ribosome-associated proteins in fungal ribosome homeostasis argues antifungal ribosome-targeting
  approaches "have lagged because of high sequence and structural conservation with human ribosomes",
  proposes ribosome-*associated* proteins (biogenesis, turnover, hibernation) as better-diverged
  targets, and notes **tavaborole** (a leucyl-tRNA synthetase inhibitor) is the only FDA-approved
  antifungal that directly targets protein synthesis —
  [J Microbiol](https://www.jmicrobiol.or.kr/journal/view.php?number=3061).
- **Eukaryotic inhibitor landscape review.** "A Quick Guide to Small-Molecule Inhibitors of
  Eukaryotic Protein Synthesis", [PMC7689648](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7689648/):
  covers ribosome-targeting compounds in anticancer, anthelmintic and antifungal use, plus
  nonsense-suppression (read-through) applications.
- **Cross-pathogen selectivity by binding-energy prediction (established study, not a review).**
  "Broad anti-pathogen potential of DEAD box RNA helicase eIF4A-targeting rocaglates",
  [PMC10250355](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10250355/): silvestrol sensitivity
  across many eukaryotic pathogens tracked with predicted binding energy and thermal shift, with
  predicted sensitivity in some parasites and resistance in others. This is the closest example found
  of a *systematic cross-species* selectivity analysis for a translation target — though the target
  is eIF4A, not rRNA.
- **Selectivity reassessment review.** "Blast from the Past: Reassessing Forgotten Translation
  Inhibitors, Antibiotic Selectivity, and Resistance Mechanisms to Aid Drug Development", *Molecular
  Cell*, [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S1097276515008102) — a
  selectivity-focused review, but 2015 rather than 2020+.
- **Classic reviews for baseline (pre-2020).** Wilson, "Ribosome-targeting antibiotics and mechanisms
  of bacterial resistance", *Nat Rev Microbiol* (2014) —
  [ResearchGate record](https://www.researchgate.net/publication/259321011_Ribosome-targeting_antibiotics_and_bacterial_resistance_mechanisms);
  Lin, Zhou, Steitz, Polikanov & Gagnon, "Ribosome-Targeting Antibiotics: Modes of Action, Mechanisms
  of Resistance, and Implications for Drug Design", *Annu Rev Biochem* (2018) —
  [ResearchGate record](https://www.researchgate.net/publication/323976028_Ribosome-Targeting_Antibiotics_Modes_of_Action_Mechanisms_of_Resistance_and_Implications_for_Drug_Design);
  Vázquez-Laslop/Mankin-era "ANTIBIOTICS TARGETING RIBOSOMES: Resistance, Selectivity, Synergism, and
  Cellular Regulation", *Annu Rev Biochem* (2005) —
  [Annual Reviews](https://www.annualreviews.org/content/journals/10.1146/annurev.biochem.74.082803.133130).
  *Caveat:* I verified these exist and are selectivity-relevant but did not read their full texts, so
  I cannot quote their stated open questions.
- **Open question stated by a 2025 primary paper.** The linezolid–mitoribosome structure's conclusion
  that the binding mode "shows extensive similarity to the mode of bacterial inhibition and also
  suggests potential avenues for altering selectivity" ([PMID 39907106](https://pubmed.ncbi.nlm.nih.gov/39907106/))
  is an explicit statement that mito-vs-bacterial selectivity for oxazolidinones is unsolved.
- **Off-target scope still being mapped (2024–2025).** The eLife reviewed preprint's finding that
  telithromycin alone among the tested inhibitors differs in mechanism between bacteria and
  mitochondria ([eLife 103699](https://elifesciences.org/reviewed-preprints/103699v2)) implies the
  per-drug mito off-target map is incomplete.

### Inferences
- The recurring framing — "too conserved to target selectively" — is asserted qualitatively in
  reviews but, as far as the sources found show, has not been **quantified** structure-wide: nobody
  appears to have measured, per binding site and per drug, how many contacting positions actually
  differ between kingdoms. That is a well-defined quantitative question the proposed project answers.
- The 2024 tigecycline result (mito-selective at therapeutic concentration, conserved second site at
  high concentration) and the telithromycin exception both suggest selectivity is often about
  *which of several sites* is occupied, not just one determinant nucleotide. A survey that catalogues
  **all** contact sets per structure (including secondary sites, as Prokhorova 2017 found for
  aminoglycosides) is better matched to that reality than a single-position lookup table.

### Gaps
- I did not locate a dedicated 2020–2026 review of eukaryotic ribosome inhibitors or of mitoribosome
  toxicity with an explicit "open questions" section; searches returned mostly primary papers and
  older reviews. A targeted search of *Nat Rev Microbiol* / *Annu Rev Biochem* / *Trends Biochem Sci*
  2022–2026 would be the next step.
- Polikanov/Mankin reviews from 2020+ were not surfaced; only the 2018 *Annu Rev Biochem* and 2005
  *Annu Rev Biochem* appeared.

## Q5. How many structures exist per kingdom for key drug classes

### Takeaway
No source was found that counts ribosome–drug structures per kingdom or per drug class. The
per-paper numbers below are the only quantitative anchors found; the count would have to be computed
directly from the PDB/EMDB.

### Cited Findings
- 16 yeast 80S inhibitor structures in one 2014 paper (12 eukaryote-specific + 4 broad-spectrum
  inhibitors) — [RCSB 4U4Q](https://www.rcsb.org/structure/4U4Q).
- 6 yeast 80S aminoglycoside-related PDB entries from Prokhorova 2017 (5NDG, 5NDJ, 5NDK, 5NDV, 5NDW,
  5OBM) — [RCSB 5NDV](https://rcsb.org/structure/5ndv).
- 2 human 80S aminoglycoside entries (arbekacin, dibekacin): [8IFD](https://rcsb.org/structure/8IFD),
  [8IFE](https://rcsb.org/structure/8IFE).
- 65 antibiotic binding sites catalogued in bacterial+archaeal structures as of 2010 —
  [NAR 38:5982](https://academic.oup.com/nar/article/38/18/5982/1065889).
- 17 compounds from 6 antibiotic classes in the 2023 NSMB cryo-EM series —
  [NSMB](https://www.nature.com/articles/s41594-023-01047-y).
- 4 ribosome species in the single tigecycline paper (human 55S, human 39S, human 80S, yeast 80S) —
  [Nat Commun 2024](https://doi.org/10.1038/s41467-024-49797-7).
- Mitochondrial translation scope context: mitoribosomes translate 13 proteins from 11 transcripts,
  all electron-transport-chain components, which is why their inhibition is linked to lactic acidosis
  and myelosuppression — [bioRxiv 2024.08.21.609012](https://www.biorxiv.org/content/10.1101/2024.08.21.609012.full.pdf).
- An RCSB/PDB-101 educational page collects ribosome–antibiotic structures but gives no counts —
  [PDB-101](https://pdb101.rcsb.org/learn/structural-biology-highlights/antibiotics-and-ribosome-function).
- A 2019-era EU project report stated that for mitoribosome–antimicrobial binding "no structural
  information is available, and the molecular interactions have not been shown" —
  [CORDIS 805230](https://cordis.europa.eu/project/id/805230/reporting/pl). The 2023–2025
  streptomycin, tigecycline and linezolid structures have since partly filled this.

### Inferences
- The absence of any per-kingdom census is itself a finding: it means the proposed project's first
  deliverable (a counted inventory of ribosome–ligand structures per kingdom per drug class, indexed
  on Rfam SSU/LSU consensus positions) would be new descriptive knowledge, independent of any
  mechanistic conclusion.
- Rough qualitative picture from what was found: bacterial structures dominate heavily; eukaryotic
  cytosolic is concentrated in yeast (2014/2017 panels) plus a smaller human set; parasite and fungal
  complexes are single-digit per organism; mitoribosome drug complexes are very recent (2023–2025)
  and few. This asymmetry is consistent with every source's framing but is not itself a cited count.

### Gaps
- No published count per kingdom or per drug class. Would need direct PDB/EMDB queries (e.g. by
  polymer entity taxonomy × nonpolymer ligand) — feasible in this repository's existing tooling but
  not available from the literature.
- Resolution values on some PDB entries disagree between RCSB and PDBj (noted in Q1); any census
  should take values from the primary mmCIF.
