# Data-quality pitfalls and validation practice for small-molecule ligands in RNA/nucleic-acid PDB structures

Scope: rules a pipeline aggregating PDBe Arpeggio RNA–ligand contacts per Rfam family (X-ray + cryo-EM incl. ribosomes) could adopt. Searches done 2026-10-23; ~20 tool calls; several primary PDFs were blocked (reCAPTCHA/429), so some thresholds are reported as gaps rather than guessed.

## 1. Ligand validation metrics and thresholds (X-ray, cryo-EM, geometry) and their use in curated datasets

### Takeaway
For X-ray, the wwPDB report's per-ligand RSCC (flag < 0.8) and RSR (flag > ~0.4 in ligand tables; > 0.3 combined with RSCC < 0.8 in the summary graphic) are the standard density-fit metrics since LLDF was retired (2018); Mogul bond/angle RMSZ > 1 and |Z| > 2 flag geometry. For cryo-EM, wwPDB XML carries per-ligand Q-score and atom-inclusion, but there is no wwPDB-endorsed ligand Q-score cutoff, and the EMDataResource Ligand Challenge concluded a composite score is needed. No RNA-specific ligand-quality thresholds were found; curated datasets reuse protein-ligand cutoffs.

### Cited Findings
- Smart et al. 2018 (Acta Cryst. D74:228–236, PDBe/Masaryk) describe the wwPDB ligand-validation pipeline and conclude "the LLDF score currently used to identify ligand electron-density fit outliers can give misleading results and that better ligand-validation metrics are required." — [Smart et al. 2018](https://doi.org/10.1107/S2059798318002541) (text extracted from [IUCr PDF](https://journals.iucr.org/d/issues/2018/03/00/ba5278/ba5278.pdf))
- LLDF = (RSR_ligand − ⟨RSR_site⟩) / σ(RSR_site), using the RSR of neighbouring polymer residues; flagged as outlier if LLDF > 2 (chosen by analogy with the RSRZ > 2 cutoff for polymer residues). RSRZ itself is computed only for standard amino acids and nucleotides (per residue type and resolution bin), so it does not exist for ligands. — [Smart et al. 2018](https://journals.iucr.org/d/issues/2018/03/00/ba5278/ba5278.pdf)
- Analysis of 589,965 X-ray ligands released up to 28 June 2017: just over 11% had RSCC < 0.8; 14% of ligands with RSCC > 0.95 were nevertheless LLDF "outliers"; at 0.8 Å resolution the median LLDF is ~2 (half of ligands flagged); about one in four ligands with poor RSCC were *not* flagged by LLDF. LLDF fails at high resolution and when the surrounding polymer also fits poorly. RSCC and RSR improve with resolution, i.e. they are resolution dependent. — [Smart et al. 2018](https://journals.iucr.org/d/issues/2018/03/00/ba5278/ba5278.pdf)
- RSCC "rule of thumb" bands used by PDBe: < 0.8 poor, 0.8–0.95 OK, > 0.95 very good. — [Smart et al. 2018 Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC5947763/table/table1)
- Smart et al. note ligand placement and binding-site waters are particularly ambiguous below 3 Å resolution, and that the Mogul |Z| threshold of 2 judges ligand geometry "far more strictly" than polymer geometry (the X-ray VTF uses Z = 5 for proteins/nucleic acids). — [Smart et al. 2018](https://journals.iucr.org/d/issues/2018/03/00/ba5278/ba5278.pdf)
- wwPDB X-ray validation report help (current): LLDF retired and replaced by RSR > 0.4 and RSCC < 0.8; RSCC ≥ 0.95 very good, ~0.90 generally OK, ≤ ~0.80 "data may not support the ligand placement"; RSR 0 = perfect, ≥ ~0.4 poor fit and/or low resolution. Overall-summary density mark when RSR > 0.3 and RSCC < 0.8 (note the two sections give different RSR cutoffs). RSR/RSCC may be unreliable for structures worse than 3 Å with residue B > 90 Å². — [wwPDB X-ray VR help](https://wwpdb-beta.rcsb.org/validation/2017/XrayValidationReportHelp)
- Same help page, geometry: Mogul bond/angle Z outside ±2 merit inspection; RMSZ expected 0–1, > 1 indicates problems; ≥ 20 CSD examples needed per bond/angle; acyclic torsion outlier if local density < 5% (≥ 15 examples); ring outlier if mean and min RMSD > 60°; geometry mark if > 40% of bonds/angles/torsions/chirals/rings are outliers. Occupancy column counts atoms with occupancy < 0.9; B-factors reported as min/median/95th/max per instance. Ligand images shown for depositor "Ligands of Interest" and any ligand > 250 Da with outliers. — [wwPDB X-ray VR help](https://wwpdb-beta.rcsb.org/validation/2017/XrayValidationReportHelp)
- EBI training: a ligand with RSCC < 0.8 and/or RSR > 0.4 is highlighted as a potential issue. — [EBI training: global quality summary](https://www.ebi.ac.uk/training/online/courses/analysing-evaluating-macromolecular-models/global-quality-assessment/summary)
- Pozharski, Weichenberger & Rupp 2013 (Acta Cryst. D69:150–167) found low RSCC for "a large number" of deposited small-molecule ligands and argue electron density is the primary evidence for a binding mode; the companion Twilight tool (Weichenberger, Pozharski & Rupp 2013, Acta Cryst. F) flagged ligands with RSCC < 0.6 among 2,815 complexes deposited by 15 Jan 2012. — [IUCr newsletter](https://www.iucr.org/news/newsletter/volume-21/number-1/ligand-electron-density-analysis); [Twilight](https://ruppweb.org/twilight/default.htm)
- Cryo-EM: Q-score (Pintilie et al. 2020, Nat. Methods) is computed per atom for proteins, nucleic acids, water, ligands and ions; wwPDB 3DEM validation reports include Q-scores since September 2022 (per-residue/per-chain averages in the mmCIF/XML validation files). Pintilie et al. 2025 (Acta Cryst. D; ~10,000 EMDB/PDB map–model pairs) treat Q-score as a reliability measure for small molecules too. — [Pintilie 2020](https://link.springer.com/article/10.1038/s41592-020-0731-1); [wwPDB Q-score announcement](https://lists.wwpdb.org/empathy/thread/7MURONPQROS43GV4TLDPKU6I5XW5FPJP); [Pintilie 2025](https://journals.iucr.org/d/issues/2025/08/00/ic5125/ic5125.pdf)
- EMDataResource Ligand Model Challenge (Lawson et al., Nat. Methods 2024): a composite rather than single score was needed to assess macromolecule+ligand models; recommends standard geometry checks plus map-fit measures. — [EMDataResource news](https://emdataresource.org/news/ligand-challenge-publication.html)
- Curated datasets:
  - PoseBusters Benchmark set (Buttenschoen, Morris & Deane, Chem. Sci. 2024): PDB entries released 1 Jan 2021–30 May 2023; ligand 100–900 Da, ≥ 3 heavy atoms, elements H/C/O/N/P/S/F/Cl only; resolution ≤ 2 Å; ligand RSR ≤ 0.2 (as summarised by search; supplementary PDF not retrievable, so other filters unverified); 428 initial cases cut to 308 after peer review found some ligands at crystal contacts. — [PoseBusters arXiv](https://arxiv.org/pdf/2308.05777); [Zenodo record](https://zenodo.org/records/8278563)
  - Iridium (Warren et al. 2012): HT/MT/NT/NA trust categories from global (DPI, Rfree, resolution) and local (ligand density coverage, active-site density coverage, ligand occupancy) criteria; NA when non-X-ray or no density. Only 17% of 728 docking-validation structures were acceptable. OpenEye notes the paper is not always explicit about numeric thresholds. — [OpenEye Iridium docs](https://docs.eyesopen.com/toolkits/cpp/sprucetk/theory/iridium.html); [Warren et al. record](https://pure.gonzaga.edu/p/en/research-outputs/essential-considerations-for-using-protein-ligand-structures-in-d)
  - Benchmark best-practice paper (2021) recommends checking ligand density visually or RSCC > 0.90. — [arXiv 2105.06222](https://arxiv.org/pdf/2105.06222)
  - BioLiP / BioLiP2 (Zhang lab; NAR 2013, NAR 2024): separate biologically relevant ligands from additives (e.g. glycerol, ethylene glycol) via an additive filter, geometric rules and literature validation. — [BioLiP2 NAR 2024](https://academic.oup.com/nar/article/52/D1/D404/7233921); [BioLiP 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531193)
  - PLINDER: 500+ annotations per system including quality; downstream users additionally filtered resolution ≤ 3.5 Å, R ≤ 0.40, Rfree ≤ 0.45, Rfree−R ≤ 0.075 (OMTRA) and excluded PLINDER-classified artifacts and single-atom ions (2603.27950). — [PLINDER PyPI](https://pypi.org/project/plinder/0.2.21/); [OMTRA](https://arxiv.org/pdf/2512.05080); [arXiv 2603.27950](https://arxiv.org/pdf/2603.27950)
  - HARIBOSS (Panei et al. 2022, RNA–small molecule): X-ray, NMR and cryo-EM PDB structures with "drug-like" ligands interacting with RNA; a secondary study using it defines ligands as < 1000 Da and binding sites as residues within 6 Å. — [HARIBOSS bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.17.492306.full.pdf); [secondary use, arXiv 2310.18985](https://arxiv.org/pdf/2310.18985)

### Inferences
- A defensible default for our X-ray contacts: flag (not drop) ligand instances with RSCC < 0.8 or RSR > 0.4, keep RSCC 0.8–0.9 as "moderate", and treat resolution > 3 Å with caution, mirroring the wwPDB report. Stricter dataset cutoffs (RSR ≤ 0.2 at ≤ 2 Å as in PoseBusters) would remove most riboswitch and nearly all ribosome entries and suit benchmarking, not aggregation.
- Cryo-EM ligands should be judged against resolution-dependent expectations for Q-score rather than a fixed cutoff, since no endorsed threshold was found (Pintilie 2020 gives an expected Q vs resolution relation; not verified here).
- Do not use LLDF anywhere (retired); old validation XMLs may still contain it.

### Gaps
- No RNA-specific ligand-quality thresholds were found in any source.
- Exact Iridium numeric thresholds, PLINDER's own ligand RSCC/RSR fields and cutoffs, and the full HARIBOSS inclusion rules (e.g. exact additive exclusion list) could not be retrieved (PDFs blocked/429).
- The suggested "Deller & Rupp 2015" and "Shao et al." ligand-quality papers were not located by search; citation unverified.
- No numerical Q-score cutoff for ligands found in wwPDB documentation.

## 2. Known artefacts: additives, buffers, ions vs water, poor density, 5'-GTP, crystal contacts

### Takeaway
Additives (glycerol, ethylene glycol, acetate, sulfate, etc.) are routinely modelled and must be filtered by a curated list; ion identity in nucleic-acid structures is frequently wrong (Mg²⁺ over-assigned where K⁺/Na⁺/water fits better), so ion "ligands" should not be aggregated as bound small molecules without coordination checks. Even in RF00167 entries, modelled additives can have poor density.

### Cited Findings
- Smart et al. 2018: deciding whether binding-site density is the soaked/co-crystallised ligand or water/buffer molecules "is often far from trivial", aggravated by partial occupancy and resolution worse than 3 Å. — [Smart et al. 2018](https://journals.iucr.org/d/issues/2018/03/00/ba5278/ba5278.pdf)
- BioLiP exists because small molecules used as crystallisation/purification additives (glycerol, ethylene glycol) are deposited as ligands; it uses an additive filter plus geometric rules. — [BioLiP2](https://academic.oup.com/nar/article/52/D1/D404/7233921); [BioLiP 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531193)
- Leonarski, D'Ascenzo & Auffinger 2019 (RNA): PDB survey of nucleic-acid X-ray structures ≤ 2.9 Å; the PDB "contains many ion assignment errors" that skew database analyses; nucleobase carbonyls are poor Mg²⁺ inner-sphere binders but excellent monovalent-ion binders. — [Leonarski et al. 2019, PMC6348993](https://pmc.ncbi.nlm.nih.gov/articles/PMC6348993/)
- Auffinger, Ennifar & D'Ascenzo 2021 ("Deflating the RNA Mg²⁺ bubble", RNA 27:243): numerous studies overestimate Mg²⁺ counts due to assignment errors; criteria are coordination distance (Mg²⁺–O ≈ 2.1 Å vs K⁺–O ≈ 2.8 Å) and strongly preferred octahedral geometry for Mg²⁺; reassessment of ribosome models 6QNR and 6SJ6 suggests K⁺, not Mg²⁺, is the prevalent ion in the first solvation shell; an ion in 4V6F assigned Mg²⁺ fits K⁺ better. For light ions without anomalous signal (Mg²⁺, Na⁺), stereochemistry plus high resolution (< 2.4 Å) is the best available option. — [Auffinger et al. 2021](https://rnajournal.cshlp.org/content/27/3/243); [PMC7901845](https://pmc.ncbi.nlm.nih.gov/articles/PMC7901845)
- Own check of the wwPDB validation XML for RF00167 entry 4FE5 (fetched 2026-10-23 from files.wwpdb.org): hypoxanthine HPA B101 RSCC 0.945, RSR 0.072, EDIAm 1.012, occupancy 1.0; two acetates in the same entry differ sharply, ACT B102 RSCC 0.962 vs ACT B103 RSCC 0.706, EDIAm 0.213, B 42 Å²; cobalt hexammine (NCO) RSCC 0.97–0.99. — [4fe5_validation.xml.gz](https://files.wwpdb.org/pub/pdb/validation_reports/fe/4fe5/4fe5_validation.xml.gz)
- Real triphosphate ligands exist (e.g. ZTP riboswitch ligands), so triphosphates cannot be excluded generically. — [search context: RNA journal supplement](https://rnajournal.cshlp.org/content/suppl/2023/07/07/rna.079427.122.DC1/Supplemental_Material.pdf)
- PoseBusters removed benchmark cases after finding ligands at crystal contacts (428 → 308). — [Zenodo](https://zenodo.org/records/8278563)

### Inferences
- Maintain an explicit, versioned additive/buffer/cryoprotectant CCD list (e.g. GOL, EDO, PEG fragments, ACT, SO4, PO4, MPD, DMS, TRS, cacodylate) and classify rather than delete; the project already has ligand classification in `annotation_config.json`.
- Treat metal ions as a separate category with their own caveats; do not present Mg²⁺ site counts as biological without flagging the Auffinger/Leonarski evidence, especially for ribosomes.
- A 5'-terminal GTP from transcription is normally part of the polymer (5'-triphosphate on the first nucleotide), so it would appear as a modified/terminal residue rather than a nonpolymer ligand; it should be checked case-by-case if a free GTP/GDP is listed. (Inference; not sourced.)

### Gaps
- No source found documenting 5'-GTP transcription artefacts being deposited as ligands in RNA entries.
- No published curated additive list (BioLiP's actual list not retrieved).
- No systematic study found of ligands at crystal contacts in RNA crystal structures specifically.

## 3. Structural pitfalls: altlocs, partial occupancy, copies vs assembly, symmetry, low-resolution cryo-EM, NMR models

### Takeaway
The wwPDB XML reports per-instance occupancy, altcode and fit metrics, so each ligand copy can be judged separately; antibiotic poses in ribosomes have been disputed (e.g. chloramphenicol, macrolides) and recent ≤ 2.2 Å cryo-EM resolves many. Low-resolution cryo-EM ligand sites often cannot determine pose without priors.

### Cited Findings
- wwPDB VR lists occupancy (< 0.9 counted) and B-factor spread (approaching 80 Å² suggests disorder) per ligand instance. — [wwPDB X-ray VR help](https://wwpdb-beta.rcsb.org/validation/2017/XrayValidationReportHelp)
- Smart et al. 2018 discuss zero-occupancy atoms and alternative conformations of ligands and the Ligand Validation Workshop (Adams et al. 2016) recommendation on handling them; partial occupancy complicates density interpretation (Pearce et al. 2017). — [Smart et al. 2018](https://journals.iucr.org/d/issues/2018/03/00/ba5278/ba5278.pdf)
- 4FE5 example: two copies of the same additive (ACT) in one entry have RSCC 0.962 vs 0.706, i.e. quality is per copy, not per entity. — [4fe5 validation XML](https://files.wwpdb.org/pub/pdb/validation_reports/fe/4fe5/4fe5_validation.xml.gz)
- Bulkley et al. 2010 (PNAS): their chloramphenicol–ribosome structure suggested a binding model "substantially" different from the prevailing one, and macrolide structures differed from published D. radiodurans models; attributed to species differences or experimental problems. — [RCSB 4V7X](https://www.rcsb.org/structure/4V7X)
- Paternoga et al. 2023 (NSMB): cryo-EM of 17 antibiotics (six classes) at 1.6–2.2 Å shows conserved binding modes within scaffolds, including ordered waters. — [mbexc summary](https://mbexc.de/structural-conservation-of-antibiotic-interaction-with-ribosomes/)
- GemSpot review: cryo-EM maps at ligand sites are often too low resolution to fix a pose without prior knowledge, and modelling can give the wrong pose. — [hypothes.is bioRxiv review](https://api.hypothes.is/users/bioRxiv_admin?page=44)
- Auffinger 2021 shows even ~2 Å ribosome models (6QNR) carry questionable ion identities. — [Auffinger et al. 2021](https://rnajournal.cshlp.org/content/27/3/243)

### Inferences
- Aggregate per ligand *instance* with its own fit metrics; the project memory notes PDBe Arpeggio covers only preferred-assembly ligand copies, so copies outside the preferred assembly or symmetry-generated contacts are absent from the source — report as "not annotated by source", not absence.
- For cryo-EM ribosomes, combine entry resolution, ligand Q-score and residue_inclusion, and flag ligands in maps worse than ~3.5 Å as pose-uncertain (threshold is a judgment, not sourced).
- NMR: no density metric exists; use model 1 or report per-model consistency; Iridium classes non-X-ray as "NA" (not assessable), a useful precedent for an explicit "not assessable" state.

### Gaps
- No source found quantifying symmetry-mate (crystal-contact) ligand contacts in RNA structures.
- No source found on best practice for multi-model NMR RNA–ligand contacts.
- Specific published reinterpretations of individual ribosome antibiotic cryo-EM poses (beyond Bulkley 2010) not located.

## 4. Programmatic access: PDBe API / wwPDB validation XML / Q-score

### Takeaway
The richest machine-readable source is the wwPDB validation XML (files.wwpdb.org), whose per-ligand `ModelledSubgroup` elements carry RSCC, RSR, EDIAm, OPIA, occupancy, OWAB and Mogul RMSZ (X-ray) or Q_score and residue_inclusion (cryo-EM). The PDBe API residuewise outlier endpoint returns outlier type labels, not ligand fit values.

### Cited Findings
- wwPDB validation XML for 4FE5 (X-ray), per-ligand attributes observed: `rscc`, `rsr`, `EDIAm`, `OPIA`, `owab`, `avgoccu`, `NatomsEDS`, `mogul_bonds_rmsz`, `mogul_angles_rmsz`, `altcode`, `said` (label_asym_id), `resnum`, `chain`. Schema v6.07 referenced in XML headers. — [4fe5_validation.xml.gz](https://files.wwpdb.org/pub/pdb/validation_reports/fe/4fe5/4fe5_validation.xml.gz)
- wwPDB validation XML for 7K00 (E. coli ribosome cryo-EM, 1.98 Å, EMD-22586): ligands carry `Q_score` and `residue_inclusion`, e.g. paromomycin PAR A1601 Q 0.717, inclusion 0.905, Mogul bond RMSZ 0.50; Mg²⁺ ions Q ≈ 0.80 with `mogul-ignore="yes"`. — [7k00_validation.xml.gz](https://files.wwpdb.org/pub/pdb/validation_reports/k0/7k00/7k00_validation.xml.gz)
- Q-scores in 3DEM validation reports since Sept 2022; per-residue/per-chain averages in validation mmCIF/XML. — [wwPDB announcement](https://lists.wwpdb.org/empathy/thread/7MURONPQROS43GV4TLDPKU6I5XW5FPJP)
- PDBe API `GET https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/{pdb}` returns per residue `outlier_types` (e.g. RSRZ, clashes, planes) keyed by entity/struct_asym_id/residue_number (live call on 4fe5, 2026-10-23); related `/validation/global_percentiles/entry/{pdb}` and `/validation/rama_sidechain_listing/entry/{pdb}`. — [ToolUniverse PDBe validation wrapper](https://zitniklab.hms.harvard.edu/ToolUniverse/_modules/tooluniverse/pdbe_validation_tool.html); [PDBe API](https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/4fe5)
- Validation maps are computed with Servalcat (refmacat) from deposited structure factors; RSR against 2mFo−DFc. — [wwPDB X-ray VR help](https://wwpdb-beta.rcsb.org/validation/2017/XrayValidationReportHelp)

### Inferences
- Implement a `ligand_quality` step that streams each entry's validation XML (small: 14 KB for 4FE5, 237 KB for 7K00), joins on `said` + `resnum` to Arpeggio ligand instances (label_asym_id), and records RSCC/RSR/EDIAm/occupancy or Q-score/inclusion plus XML revision date as provenance. This fits the PDBe-first rule (XML is wwPDB, mirrored at PDBe).
- Ligand-level quality from the PDBe API appears limited to outlier labels; whether the residuewise endpoint includes nonpolymer ligands was not confirmed.

### Gaps
- Not verified whether PDBe API exposes ligand RSCC/RSR/Q-score values directly (no official ligand-validation endpoint found).
- No documented ligand-specific Q-score threshold in wwPDB/EMDB reports.
