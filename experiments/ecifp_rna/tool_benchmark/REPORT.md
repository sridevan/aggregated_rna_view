# RNA–ligand interaction tools: pilot benchmark (2026-10-23)

Exploratory. Survey behind it: `experiments/ecifp_rna/literature/reports/RNA ligand interaction annotation tools.md`.

**Question.** Do other tools annotate RNA–ligand contacts better than the PDBe-Arpeggio contacts we use (PDBe API `bound_ligand_interactions`)?

## Set-up
- **Ligand copies (12 from 10 entries):**
  - 1Y26 ADE and 4FE5 HPA (purine)
  - 2GIS (SAM-I), 2QWY (SAM-II), 3E5C (SAM-III), 6FZ0 (SAM-V) and 6LAS (SAM-VI), all SAM
  - 3NPQ SAH (SAH riboswitch)
  - 2GDI TPP (thiM)
  - 1FJG paromomycin, streptomycin and spectinomycin (30S)
- **Reference** (`reference/`): 107 contacts curated from the primary papers by a subagent, from text (78), figures (28) or a secondary source (1).
  - Built only from what the papers say; coordinates were used only to confirm that the residues exist.
  - `CURATION_NOTES.md` records numbering (3E5C differs from its paper; 1FJG uses E. coli numbering), conflicts, and which contact classes each paper claims to list completely.
  - The curation script is `reference/build_ref.py`.
- **Inputs** (`prepare.py`). Model 1, first alternate location, from the PDBe updated mmCIF (kept in scratch only).
  - **Crop:** the ligand plus nucleic-acid residues, waters and single-atom ions within 12 Å.
  - **Proteins left out** for every tool, because fingeRNAt takes nucleic-acid receptors only. This drops the S12 contact for streptomycin, which is not in the reference either.
  - **Common hydrogens:** hydride 1.2.3 on CCD bond orders and CCD formal charges. Backbone OP1 is set anionic and H on OP1/OP2/OP3 removed.
- **Tool arms** (`run_tools.py`, normalised records in `out/pilot_records.json`):

  | Arm | What |
  |---|---|
  | `arpeggio_api` | PDBe API contacts: production PDBe-Arpeggio on the full entry, which is what the project uses |
  | `arpeggio_local` | pdbe-arpeggio 1.4.4 (pip, openbabel-wheel) on the crop |
  | `plip` | PLIP 3.0.1 `--dnareceptor`, with its own OpenBabel hydrogens |
  | `plip_H` | PLIP `--nohydro` on the common hydrogens |
  | `fingernat_H` | fingeRNAt cc62eb8 FULL `-h2o -detail` on the common hydrogens; it requires hydrogens |

- **Scoring** (`score.py`, results in `out/pilot_scores.json`).
  - **Classes compared:** polar (hydrogen bonds plus ionic), stacking, water-mediated and metal-mediated.
  - **Recall:** residue-level against the reference, plus atom-level for polar contacts where the paper names the RNA atom (51).
  - **Tool contacts not in the reference:** counted as "unsupported", not as false positives, because most papers don't claim to be complete.
  - **Agreement:** residue-level Jaccard between tools.

## Results

| Arm | Polar | Stacking | Water | Metal | Polar atoms | Unsupported polar / stacking / water |
|---|---|---|---|---|---|---|
| arpeggio_api | 48/59 | **14/20** | 0/1 | 0/3 | **36/51** | 21 / 21 / 0 |
| arpeggio_local | 26/59 | 0/20 | 0/1 | 0/3 | 13/51 | 13 / 0 / 0 |
| plip | 46/59 | 11/20 | 1/1 | 0/3 | 31/51 | 24 / 3 / 12 |
| plip_H | 49/59 | 11/20 | 1/1 | 0/3 | 29/51 | 22 / 3 / 13 |
| fingernat_H | **55/59** | 10/20 | 0/1 | **2/3** | 35/51 | **42** / 3 / 1 |

**Residue-level agreement (Jaccard).**
- **Polar:** API vs PLIP 0.62–0.69, API vs fingeRNAt 0.69, PLIP vs fingeRNAt 0.65–0.68.
- **Stacking:** PLIP vs fingeRNAt 0.93; API vs either 0.37–0.40. The API counts every ring–ring geometry, including edge-to-face and offset, so it finds more stacked residues, some of which the papers don't describe.

**Findings**
1. **No tool recognises the SAM sulfonium–carbonyl contact.** This electrostatic interaction is what makes SAM riboswitches selective for SAM over SAH.
   - **The papers state it:** S⁺···O4 of U11/U21 (SAM-II), U9/U20 (SAM-V, 3.2–3.3 Å) and U6 (SAM-VI, 4.8 Å), and S⁺···O2 of U7 (SAM-I, 4.0 Å).
   - **What the tools miss:** Arpeggio (API) misses all six; PLIP misses all but SAM-I U7; fingeRNAt picks up only U7 and U21, probably as distance-only hydrogen bonds.
   - **No type exists for it:** none of the tools has a cation–lone-pair or sulfonium-carbonyl interaction type. The geometric sulfonium check in `sam_recognition.py` catches exactly these contacts, so the project needs its own rule for them.
2. **Metal-mediated contacts.**
   - **fingeRNAt is the only tool that detects them:** it finds both TPP pyrophosphate–Mg²⁺–G60/G78 bridges.
   - **Missed by every tool:** SAM-V's Mg²⁺, which binds through its coordinated waters (outer sphere), so a direct-distance cutoff can't see it.
   - **Arpeggio and PLIP** report no ligand–metal–RNA bridge at all.
3. **Water bridges.** The papers name only one, which PLIP finds. PLIP and fingeRNAt report 12–13 more and 1–2 more respectively, so the reference can't score water.
4. **fingeRNAt finds the most polar contacts and the most extras.** It recovers 55 of 59 polar contacts but has 42 unsupported ones. Its hydrogen bond is a donor–acceptor distance up to 3.9 Å with no angle check, so it over-calls.
5. **Stacking is detected best by the PDBe API Arpeggio** (14/20).
   - **1Y26:** most misses for every tool are the base triples above and below the adenine. The reference splits these per residue, so not every residue overlaps the ligand ring, and residue-level scoring is harsh here.
   - **6LAS:** the G9–C32 / G7·G33 pairs, same pattern.
6. **The local pdbe-arpeggio install perceives no aromatic rings in RNA:**
   - So it reports no stacking or cation–π, and its polar recall drops to 26/59.
   - The production PDBe contacts do include ring contacts, so this is an OpenBabel build or version effect in the open-source install. Anyone rerunning Arpeggio locally (for example to cover non-preferred-assembly copies) has to check ring perception first.
7. **Hydrogens matter less than typing.** PLIP with its own hydrogens and PLIP with the common hydride hydrogens agree at 0.76 Jaccard (polar) and 1.0 (stacking), and their recall differs by 3.
8. **Ligand charge states come from the CCD** (aminoglycoside amines neutral, TPP pyrophosphate protonated). This undercounts ionic contacts for 1FJG and TPP; a pH-7 ligand-protonation arm is still to do.
9. **fingeRNAt is fragile with its inputs.**
   - **Kekulé bonds needed in the ligand SDF:** with aromatic query bond codes it silently finds no stacking or charged contacts.
   - **Formal charges must be supplied:** deposited files often leave them empty.
   - **Duplicate rows in its detail output:** removed during normalisation.

## Reading
- **No single tool is better than the PDBe API contacts across the board.** The API has the best stacking and atom-level polar recall, and polar recall comparable to PLIP's.
- **What the other tools add:**
  - fingeRNAt adds ion-mediated bridges (direct coordination only).
  - PLIP adds angle-checked water bridges.
- **The SAM sulfonium contact needs a project-side rule** (cation to carbonyl or ring oxygen within about 4 Å). No tool provides it.
- **Recommendation for now:**
  - Keep the API contacts.
  - Add ion- and water-mediated contacts from coordinates, using fingeRNAt-style per-ion cutoffs plus PLIP-style water-bridge angles, or run fingeRNAt and PLIP for those classes only.
  - Add the cation–carbonyl rule.

## Caveats and next steps
- **Size:** 12 copies and 107 reference contacts. Water and metal have very few reference contacts.
- **Completeness:** "unsupported" is not "false". Papers are complete only for some classes (see CURATION_NOTES).
- **Next:**
  - score stacking per unit (base pair or triple) rather than per residue;
  - add a pH-7 ligand-protonation arm;
  - add Mol*, if a batch route exists;
  - expand to about 30–40 copies (TPP variants, more 16S drugs, preQ1, FMN, glmS) for per-class numbers with confidence intervals.
