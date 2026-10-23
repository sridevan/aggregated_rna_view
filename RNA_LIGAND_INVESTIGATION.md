# RNA–ligand interactions: investigation summary and hand-over

*Status as of 2026-10-23. Exploratory work only: nothing here is in the pipeline, the data files or the site, and nothing is committed yet. This file is the starting point for picking the work up again.*

## 1. Where we are

**Question.** Can PDBe's ECIFP interaction fingerprints, and related ideas, add a binding-site and binding-mode view to the Ligands tab, with sites placed on Rfam consensus positions?

**Answer so far:** yes, for riboswitches and ribosomes. That scope was the user's decision on 2026-10-23: not every family, since most have few or no small-molecule ligands.

**What works**
- PDBe Arpeggio contacts exist for RNA–ligand complexes.
- ECIFP can be computed from them unchanged.
- Sites defined as sets of contacted Rfam positions come out clean.
- Contacts can be classified by nucleotide edge.

**What changed after the literature review.** Excluding ions and water hides the main recognition element for many RNA ligands. Per-copy ligand quality and ion identity also need handling before anything is built.

**Next decision.** Implement the riboswitch version (section 7), after adding the top literature improvements (section 6). Ribosomes need more data work first.

## 2. Files

| Path | What |
|---|---|
| `experiments/ecifp_rna/REPORT.md` | Detailed experiment log with every result and number, in sections: ECIFP basics, purine poses, 16S antibiotics, clan, base edges |
| `experiments/ecifp_rna/ecifp_rna.py` | ECIFP reproduction plus sites on Rfam positions. Functions `configure`, `collect(acc, fam_dir)`, `cluster(rows)`; CLI `run` per family |
| `experiments/ecifp_rna/purine_poses.py` | Pocket-superposed ligand pose comparison (RF00167) |
| `experiments/ecifp_rna/clan_cross.py` | Sibling-family chains mapped to the main family's CM, a 3D pocket check, and ECIFP across the clan |
| `experiments/ecifp_rna/usalign_check.py` | US-align (pylelab/USalign 23f14ff, built with `g++ -std=c++17 -O3 -ffast-math USalign.cpp`) structural correspondence vs CM correspondence, pairwise and `-mm 4`; output `out/clan/usalign_check.json` |
| `experiments/ecifp_rna/base_edges.py` | Leontis–Westhof edge / face / ribose / phosphate classification of contacts per Rfam position |
| `experiments/ecifp_rna/sam_classes.py` | All SAM/SAH riboswitch classes (9 families) mapped on their own CMs; entries without a PDBe Rfam mapping were found by RCSB text search and confirmed by `cmscan`. Output `out/sam/<ACC>/` (family-like directories) |
| `experiments/ecifp_rna/sam_recognition.py` | SAM recognition across classes: contacts by ligand moiety × RNA edge, SAM conformation (χ, N6–CA), sulfonium environment, water and ion bridges, cross-class ECIFP. Output `out/sam_recognition.json` |
| `experiments/ecifp_rna/literature/reports/RNA ligand interaction annotation tools.md` | Survey of RNA–ligand interaction tools and resources (Arpeggio, PLIP, fingeRNAt, ProLIF, Mol*, DSSR, HARIBOSS …) with a shortlist and benchmark design |
| `experiments/ecifp_rna/tool_benchmark/` | Pilot benchmark (12 ligand copies): `prepare.py` (crops, common hydrogens), `run_tools.py` (Arpeggio API/local, PLIP, fingeRNAt → one schema), `score.py`, literature reference in `reference/`, results in `REPORT.md` and `out/` |
| `experiments/ecifp_rna/out/` | Results: `<ACC>.json` (fingerprints, sites, modes, matrices), `<ACC>_edges.json`, `RF00167_poses.json`, `clan/` (cross mappings, `clan_RF00167.json`, `clan_RF00162.json`) |
| `experiments/ecifp_rna/literature/reports/RNA small molecule 3D interactions.md` | Literature review with citations and 11 prioritised improvements |
| `experiments/ecifp_rna/literature/research_notes/RNA small molecule 3D interactions/` | Five research-note files behind the review |

**Source.** ECIFP is from Kunnakkattu et al., bioRxiv 10.64898/2026.04.17.719144; code at github.com/PDBeurope/ecifp, commit a21de95, CC0. Only the abstract and the code were read; the full text was rate-limited.

**Re-running.**
- The scripts need RDKit, pdbeccdutils, gemmi, numpy, scipy and requests. The environment used here was a throwaway one: `uv venv -p 3.12`, then `uv pip install rdkit pdbeccdutils gemmi numpy scipy requests` (RDKit 2026.09.1, pdbeccdutils 2.0.0, gemmi 0.7.5).
- API caches were in the session scratchpad and are gone, so re-runs re-fetch from the PDBe API, PDBe CCD and the ModelServer.
- Run from the repository root:

      python -I experiments/ecifp_rna/ecifp_rna.py <repo> <cache> experiments/ecifp_rna/out RF00167 RF00059 RF00162 RF00177
      python -I experiments/ecifp_rna/purine_poses.py <repo> <cache> experiments/ecifp_rna/out/RF00167_poses.json
      python -I experiments/ecifp_rna/clan_cross.py <repo> <scratch> experiments/ecifp_rna/out/clan 4 RF00167:RF01510 RF00162:RF00634,RF01725
      python3 experiments/ecifp_rna/base_edges.py <cache> experiments/ecifp_rna/out RF00167 RF00059 RF00162 RF00177   # needs the cache from ecifp_rna.py

## 3. Method as prototyped

**Contacts.** PDBe API `bound_ligand_interactions/{pdb}/{chain}/{num}`, which are Arpeggio contacts on the preferred assembly.
- The contact types proximal, clash, vdw_clash and vdw are dropped, along with water and the other components of the same bound molecule.
- Only contacts with family chains are kept.
- Ligands in our `ligands.json` classes water, metal_ion and additive are skipped.

**ECIFP.** Exactly as in PDBe's code:
- Morgan fingerprints (radius 2, 1024 bits) of the CCD component; only the bits of the contacting atoms are kept, for the ligand and for each contacted residue type.
- Polymer bits are merged by residue type, so the fingerprint carries no positions.
- Similarity = Tanimoto(ligand)^0.4 × Tanimoto(polymer)^0.6.

**Sites.** The contacted residues are mapped to Rfam consensus positions through our validated mappings. Sites are grouped by the overlap coefficient with single linkage at a distance of 0.4. This is ECIFP's UniProt-based site rule, with Rfam positions in place of UniProt numbers.

**Binding modes.** Average linkage on 1 − ECIFP, cut at 0.5. This is a prototype choice; the paper uses minimax linkage with the cut at the largest gap between 0.3 and 0.6.

**Edges** (`base_edges.py`; the atom table is a draft for review):

| Class | Atoms |
|---|---|
| Watson-Crick (W) | A N1; G N1; C N3; U N3 |
| Hoogsteen (H) | A/G N7, C8; C/U C5, C6 |
| Sugar edge (S) | A/G N3; O2′ |
| W/H corner | A N6; G O6; C N4; U O4 |
| W/S corner | A C2; G N2; C/U O2 |
| Base face | Arpeggio ring-centroid (π) contacts |
| Ribose | C1′ C2′ C3′ C4′ O4′ C5′ |
| Phosphate | P OP1 OP2 OP3 O5′ O3′ |

**Rules agreed on 2026-10-23**
- Base edges, including O2′, count only from hydrogen-bond-type contacts (hbond, polar, ionic, metal_complex).
- Ribose and phosphate count from any contact that passes the ECIFP filter.
- A corner counts only if no single-edge atom of the same base was contacted.

## 4. Results

### Coverage
- **Riboswitches:** every non-ion ligand copy has contacts. The counts are RF00167 71/71, RF00059 28/28 and RF00162 27/37.
- **The 10 missing SAM-I copies:** these are assembly copies on chain `A_2`, which could be recovered by mapping `X_n` to X.
- **16S:** only 32 of 87 copies are usable.
  - 28 return no PDBe contacts (HTTP 404). These include all apramycin (7PJS, 7PJV, 7PJY), SY5 (7B5K) and some tetracyclines in 9PIH–9PIJ. The cause is unknown and is not a chain-case issue.
  - 17 contact only chains outside our selection: the second ribosome or 23S.

### Sites on Rfam positions
- **One main site per riboswitch:**
  - purine: 25, 26, 51, 55, 56, 80, 81;
  - TPP: 15, 33, 35, 36, 76, 81;
  - SAM-I: 7, 11, 46, 48, 69–71, 103, 104.
- **Purine 5′-GTP artefact sites:** at positions 19 and 20, apart from the main site. These are already flagged in `ligands.json`.
- **16S sites:**
  - tetracycline (1044, 1186–1189);
  - aminoglycoside A-site (1396–1399, 1482–1486);
  - spectinomycin (1054–1058, 1182–1185, 1377–1379);
  - kasugamycin;
  - several one-copy secondary aminoglycoside sites.

### Binding modes and poses
- **Purine.**
  - **ECIFP groups:**
    - adenine-like (U at position 80 in 14 of 19);
    - guanine-like (C at 80 in 19 of 19);
    - hypoxanthine/xanthine.
  - **The pose check (`purine_poses.py`) shows these are one pose:** the purine core sits within 0.5 Å of adenine in 38 of 55 copies, and no copy is flipped.
  - **Larger deviations of about 1 Å:** guanosine (3DS7, 6P2H), GMP, the OXG copies in 9NDD/9NBC, and 5SWE.
  - **So ECIFP groups record contact chemistry, not pose.**
- **TPP.** TPP-like compounds form one group (19 copies). The 7TZR–7TZU series is a separate group (6), which lacks most of TPP's contacts at 61, 76 and 77.
- **SAM-I.** One mode (27 copies).
- **16S.**
  - The tetracyclines (6 copies, 5 drugs) form one tight mode: median ECIFP 0.75, against 0.56 for the ligand alone.
  - Between classes the median is 0.10 (tetracycline vs aminoglycoside).
  - At a shared site, negamycin, hygromycin B and thermorubin each stay separate.
  - The 2-DOS aminoglycosides split into paromomycin + dibekacin + arbekacin and gentamicin + G418. That is not the textbook 4,5 vs 4,6 split, and the reason was not examined.

### Fingerprint caveats
- **Position-free polymer side.** The same chemistry in another region scores high: in 9V4Y, an ANG at the 5′ end grouped with the pocket ANG. So sites must be defined first, with ECIFP used only within a site.
- **Small RNA vocabulary.** Four nucleotides on one backbone push scores higher than for proteins. Hypoxanthine in TPP vs purine pockets scores 0.62, against 0.94 within the purine family.

### Clan comparisons (`clan_cross.py`)
- **SAM clan.** SAM-IV and SAM-I/IV align poorly to the SAM-I CM overall: SAM-IV scores −6.5 bits, with 55 of 119 nucleotides in local ends. Their pockets still correspond, with 10 of 14 pocket positions mapped, an RMSD of 0.57–0.62 Å and 1.24 Å for 6UES. Their SAM copies join the SAM-I site and mode, with ECIFP 0.90–1.00.
- **Purine clan, 2dG-I (RF01510).** It scores 10–12 bits on the purine CM, and the pocket RMSD is 2.2–2.7 Å. By region:
  - J3/1 and P1 (79–82) correspond;
  - J2/3 at 54–57 roughly corresponds;
  - J1/2 (26–28) is probably off by one nucleotide;
  - J2/3 at 51–53 has no structural equivalent.

  2′-dG and its nucleotides still share one mode with the purine family's GNG and GMP copies, at ECIFP 0.72–0.96.
- **Correction (US-align, `usalign_check.py`).** The "51–53 no equivalent" finding was an artefact: the superposition went through the CM-mapped, shifted J1/2 and P2 positions. US-align puts 51–57 and 79–82 in agreement with the CM, at 1.0–2.7 Å. Only J1/2 (26–28) and P2/L2 (30–48) are one nucleotide off in the CM mapping. At J1/2, US-align puts A30 at 26 at 0.5 Å.
- **SAM clan with US-align.** The pocket positions agree with the CM. US-align also recovers positions 102–105, which the CM left in skipped stretches: 0.2–0.7 Å for SAM-I/IV and 1.3–4.3 Å for SAM-IV. A single pairwise run on 6UES was one nucleotide off across P1. The multiple alignment (`-mm 4`) fixes this, so use the multiple alignment.
- **Conclusion.** Cross-family correspondence should come from a structure-based multiple alignment (US-align `-mm 4`), cross-checked against the CM mapping per position, with C3′ distances as the reliability measure. The whole-CM bit score is a poor guide.

### SAM recognition across all SAM riboswitch classes (`sam_classes.py`, `sam_recognition.py`)
Covers SAM-I, -II, -III, -IV, -V, -VI, SAM-I/IV and SAM/SAH, plus the SAH riboswitch. Full tables are in REPORT.md.
- **Coverage gap.** PDBe has no Rfam mapping for SAM-II, SAM-III, SAM-VI, or most SAM/SAH entries. These were found by text search and `cmscan`. Copies per class: 1–3, except SAM-I (27) and SAM/SAH (21).
- **Three strategies for the sulfonium.**
  - **SAM-I clan:** a folded, syn SAM (N6–CA 5.3–6.3 Å), with S⁺ and CH₃ against uracil O2 (W/S corner).
  - **SAM-II and SAM-V:** an extended, anti SAM (11.6–11.8 Å), with S⁺ and CH₃ in a cradle of three uracil O4 carbonyls (W/H corner) plus an adenine N6.
  - **SAM-III and SAM-VI:** no carbonyl cradle within 4 Å.
- **SAM/SAH riboswitch.** It reads only the adenosine: the amino acid is unmodelled in 14 of 21 copies, and adenosine, AMP and MTA bind the same site.
- **Adenine.** Every class contacts it through the Watson–Crick edge plus stacking.
- **ECIFP follows the strategy, not the Rfam family.** SAM-I vs SAM-I/IV scores 0.92, and SAM-II vs SAM-V 0.90; other class pairs score 0.5–0.78.
- **SAM-II and SAM-V share one pocket.** US-align pairs all 9 SAM-II pocket nucleotides with SAM-V pocket nucleotides (TM 0.43), but Rfam has no clan for the two. Cross-family sites therefore need structural correspondence beyond Rfam clans.

### Interaction tools: survey and pilot benchmark (`tool_benchmark/REPORT.md`)
- **The survey.** No surveyed tool has a curated nucleotide atom table or assigns Leontis–Westhof edges. Only PLIP, fingeRNAt and Mol* explicitly support RNA as a receptor.
- **The pilot.** 12 copies, scored against 107 contacts curated from the papers.
  - **Polar contacts:** the PDBe API (Arpeggio) found 48 of 59, PLIP 46–49, fingeRNAt 55 (with 42 contacts the papers don't list).
  - **Stacking:** the API found 14 of 20, PLIP 11, fingeRNAt 10.
  - **Metal-mediated contacts:** only fingeRNAt detects them (2 of 3).
  - **Water bridges:** only PLIP checks their angles.
- **No tool sees the SAM sulfonium–carbonyl contact,** the electrostatic contact that makes the riboswitch select SAM over SAH. A project-side rule is needed.
- **The local pip pdbe-arpeggio sees no aromatic rings in RNA,** so it reports no stacking, unlike PDBe production.
- **Recommendation:** keep the API contacts and add bridges and the cation–carbonyl rule from coordinates.

### Edges (`base_edges.py`)
Examples, as copies with the contact out of copies at the site:

| Family | Contacts |
|---|---|
| Purine | Position 80 by its Watson-Crick edge (58/66); 26 by sugar edge plus face (59/66); 51 by the W/S corner (37/66) |
| SAM-I | Watson-Crick edge at 69 and 11; sugar edge at 46; face plus ribose at 48; ribose at 7, 71 and 104 |
| TPP | Sugar edge at 15; face at 35 and 36; phosphate at 61, 62 and 77 (11–12/28) |
| 16S tetracycline site | Mostly phosphate and ribose |
| 16S aminoglycoside site | W/H corners plus face at 1482 |

## 5. Decisions made

- **Scope:** riboswitches and ribosomes only (2026-10-23).
- **Edge information:** encode it as Watson-Crick, Hoogsteen or sugar edge, plus face, ribose and phosphate (2026-10-23).
- **Backbone contacts** do not need the hydrogen-bond rule; base edges do (2026-10-23).
- **Insertions** stay out of the aggregation. Work is on consensus positions; decided earlier on 2026-10-23 in the mapping review.

- **Scope of the resource paper (user decision, 2026-10-23):** the aggregated RNA view's resource paper does not need this ligand-interaction layer. Showing ligand-contacting residues as a sequence track, as the Ligands tab already does, is sufficient. The interaction-tool benchmark is the candidate stand-alone publication, if pursued.

## 6. Literature: what to change

The full review, with citations, is in `experiments/ecifp_rna/literature/reports/RNA small molecule 3D interactions.md`. Its 11 improvements, in priority order:

1. **Ion- and water-bridged contacts.** fingeRNAt's cutoffs give a starting point: Mg 3.2, Na 3.6, K 3.9, other ions 3.5 and water 3.5 Å. Each structure would carry a "solvent modelled" eligibility flag and a resolution stratum. Bridged contacts are essential for TPP, glmS, glycine, fluoride, ZTP, NAD⁺, SAM-V and most ribosome antibiotics.
2. **Per-copy ligand quality** from the wwPDB validation XML: RSCC and RSR for X-ray, Q-score for cryo-EM. Copies would be flagged, not dropped, at RSCC < 0.8 or RSR > 0.4. For example, the two acetates in 4FE5 have RSCC 0.96 and 0.71.
3. **A versioned additive/buffer list**, seeded from BioLiP's 463 artefact ligands. Triphosphates cannot be excluded generically.
4. **Ions as their own class** with plausibility checks. Mg²⁺ is often mis-assigned: frequently K⁺, or density near N7.
5. **Sharper edge definitions.**
   - State the O3′/O5′ convention.
   - Record the Leontis–Westhof family of each ligand pseudo-pair.
   - Add cation–π and electrostatic labels.
   - Spot-check Arpeggio's typing of RNA atoms. Its documentation says it supports proteins only.
6. **Similarity.**
   - Weight van der Waals and lipophilic contacts down.
   - Add a Tversky "reference recovery" score.
   - Use symmetry-corrected ligand RMSD (spyrmsd) for copies of the same ligand.
7. **Cut-offs derived per site from the data.** Use average or complete linkage for modes, keeping single linkage only for site overlap.
8. **Validation benchmarks.**
   - Edge fractions against Seelam 2019 (Watson-Crick > sugar > Hoogsteen) and Kligun 2013 (65% base / 35% backbone).
   - Whether ECIFP separates guanine, adenine and 2-aminopurine in RF00167.
   - Residue sets against HARIBOSS (CC0, 2,131 entries).
   - Tetracycline, paromomycin and kasugamycin as 16S checks.
9. **Cross-family labels.** "Homologous site" (CM positions plus 3D check, optionally ARTEMIS or US-align) vs "analogous recognition" (a convergent motif). Report position overlap, not a yes/no.
10. **Apo vs holo.** Use ligand-bound chains as the denominator, and treat apo chains as a separate pocket state.
11. **Modified nucleotides at drug sites** (m8A2503, m6,6A2058): route them through parent-CCD handling, never drop them.

## 7. Proposed build (not started)

**Riboswitches first:** RF00167, RF00059, RF00162, then any other riboswitch families we add.

**Backend.** A `binding_modes` command writing `data/families/<ACC>/binding_modes.json`, with a schema and tests. It would hold:
- the sites on Rfam positions (member copies, core positions, counts);
- an edge profile per position;
- ECIFP groups within each site;
- per-copy quality flags;
- copies with no contacts, each with a reason;
- provenance (ECIFP commit, RDKit and pdbeccdutils versions, PDBe API).

New backend dependencies, RDKit and pdbeccdutils, still need the user's approval.

**Frontend.** A "Binding sites" section in the Ligands tab: per site, its ligands, the core positions with their edge profile on the tracks, and the groups.

**Before building**, settle on:
- improvements 1–4 from section 6;
- the edge atom table.

**Ribosomes, later.** They need:
- a ligand-focused selection: all antibiotic-bound 16S chains (one per entry) and 23S (RF02541);
- the 28 missing PDBe contact sets, either by asking PDBe why or by running Arpeggio locally.

## 8. Open questions

- The edge atom table, and whether corner atoms should stay as corners or be assigned to one edge.
- Whether to include ion- and water-bridged contacts, and how to define eligibility by resolution.
- Quality thresholds for cryo-EM ligands. No published cut-off exists for the Q-score.
- Why PDBe returns no contacts for some 16S ligand copies.
- The 2-DOS aminoglycoside split in 16S (paromomycin with dibekacin and arbekacin; gentamicin with G418).
- Whether ECIFP for RNA should eventually live at PDBe level, keyed by Rfam, rather than in this project.
- **Interaction tools.** Whether to scale the pilot benchmark to about 30–40 copies (with a pH-7 ligand-protonation arm and per-unit stacking), and whether to report the local Arpeggio ring-perception problem to the PDBe-Arpeggio maintainers.
- **SAM classes.** Whether to (a) report SAM-II/SAM-V as one structurally corresponding pocket although Rfam has no clan for them, and (b) add the families missing from PDBe's Rfam mappings (SAM-II, SAM-III, SAM-VI, the SAM/SAH crystal structures) to the project through `cmscan`-confirmed whole-chain segments, or report them to the Rfam team instead.
- Whether to build a clan-level CM from a structure-based seed (US-align `-mm 4` plus each family's seed sequences), or to keep RF00167 positions and store a per-sibling correspondence table. The latter is the current suggestion.

## 9. Related project items

- **Mapping review (2026-10-23), with three open items:**
  - the stage-1 `has_coordinates` check counts atoms from all models, while the analyses use model 1;
  - `docs/RFAM_PDB_MAPPING.md` is out of date (31 of 60 16S CM spans differ from Rfam; RF00005 and RF00011 are missing from its tables);
  - there is no check of each modified nucleotide's parent base against the CCD.
- **Draft slides on the mapping:** https://claude.ai/artifact/E1szKiUJ1o454akb3VuwfS
- **Existing ligand work:**
  - `LIGANDS_TAB_VALIDATION.md` and `experiments/ligand_annotation_investigation/REPORT.md`;
  - the fpocketR pockets;
  - `ligands.json` (PDBe binding sites, at entity level).
