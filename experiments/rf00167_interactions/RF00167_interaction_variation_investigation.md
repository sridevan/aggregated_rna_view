# RF00167 (purine riboswitch): RNA–RNA interaction conservation and variation

Date: 2026-10-18. This is an exploratory investigation only.

**Unchanged:**
- no application code, schema or user-facing dataset;
- the research worktree `~/Desktop/aggregated_rna_clustering`, which was only read.

**Conventions:**
- Rfam consensus positions (CM match states, 1–102) are the coordinate system.
- Residue labels such as "U22–A52" use the numbering of the *xpt* G-riboswitch (4FE5 B) or the *add* A-riboswitch (1Y26 X, same numbering) at that consensus position. Most papers use this numbering.

**Scripts and outputs (this directory):**

| File | What it contains |
|---|---|
| `map_interactions.py` | Stage 1: BGSU unit ID → chain instance → polymer position → Rfam state, verified against the mmCIF |
| `analyse.py` | Stage 2: eligibility, conservation, regions, structural groups, representatives |
| `vv_junction_table.py` | Per-chain table of the adenine-riboswitch junction states |
| `RF00167_mapped_interactions.csv` | 9,991 normalised intra-chain interactions with provenance |
| `RF00167_interaction_conservation.csv` | One row per (category, consensus pair): eligibility breakdown, frequency, entries, FR3D classes |
| `RF00167_cluster_interaction_comparison.csv` | Group vs rest, both metrics, \|Δf\| ≥ 0.2 |
| `RF00167_representative_comparison.csv` | Pairwise comparison of group representatives |
| `RF00167_vvulnificus_junction_states.csv` | The 13 identical *add* sequences plus 1Y26 |
| `instance_metadata.csv`, `base_completeness.csv`, `inter_chain_interactions.csv` | Supporting tables |
| `stage1_summary.json`, `stage2_summary.json` | Run summaries |

Reproduce from the repository root:

    python3 experiments/rf00167_interactions/map_interactions.py && python3 experiments/rf00167_interactions/analyse.py
    (cd experiments/rf00167_interactions && python3 vv_junction_table.py)

## Short answers

1. **Conserved interactions:**
   - All Rfam `SS_cons` helix pairs (P1–P3).
   - The L2–L3 kissing-loop Watson–Crick pairs: G37–C61 in 80 of 80 eligible chains (65 entries); G38–C60 in 75 of 80.
   - The junction pairs U22–A52 (73/80), A23–G46 (69/80) and C50–U75 (71/80).
   - Each class of these pairs is invariant.
2. **Variable interactions:**
   - Variation is concentrated in the three-way junction (J1/2, J2/3, J3/1), at the P3 base and in L3. Helix pairs are almost invariant.
   - 13 of the 14 junction base pairs with at least 10 eligible chains have frequency below 0.9. 18 of the 22 helix pairs are at 0.9 or above.
   - The most variable positions are G46, C53, A66, G72, G55 and A24.
3. **Associations:** almost all strong differences have an identifiable cause:
   - **Sequence subclass:** guanine-II and xanthine-II riboswitches (TM-score group 2).
   - **Ligand state, among identical sequences:** the apo and intermediate XFEL structures of the *add* riboswitch (5E54, 5SWD; TM-score group 3, geometric-discrepancy group 2).
   - **Construct:** the cryo-EM scaffold 9NDD (TM-score group 4).
   - **Point mutants:** for example the U25A mutants.
   - Among same-sequence guanine riboswitches bound to different purine analogues, no interaction difference consistently follows the ligand. The differences are scattered helix stacks.

The clearest result reproduces a published finding chain by chain (§F, §G): in the *add* adenine riboswitch, the bound-state junction network is replaced in the apo and intermediate chains by the alternative Watson–Crick pairs that Stagno et al. (2017) describe.

## A. Data availability and mapping validation

### Already available

| Data | Where |
|---|---|
| 66 PDB entries, 84 chain instances | `instances.json` |
| Validated per-chain residue correspondence (author number and insertion code → label_seq_id → match, insert or flank) | `data/families/RF00167/mappings/*.json` |
| Structural groups for both metrics, with representatives | `structure_clustering.json` |
| Ligand contacts per entry | `ligands.json` |
| Reference structural elements (P1–P3, J1/2, L2, J2/3, L3, J3/1, transferred from 4FEO) | `reference_structural_elements.json` |
| FR3D per-entry CSVs, base pairs and stacking, for all 66 entries | Fetched on 2026-10-18 by the research worktree's `bp_fetch.py`; `raw/bgsu` with `PROVENANCE.json` (URL, date, sha256). Read in place, not copied |
| PDBe updated mmCIF for all 66 entries | Local, git-ignored `cache/`. Read for verification only |

Nothing new was downloaded or recomputed.

**Missing:** BGSU's CSV export has no near base pairs (no `n`-prefixed base-pair classes). Near stacks (`ns35`, `ns53`, `ns33`, `ns55`) are present and are kept separate from observed stacks.

### Unit ID handling

- **Format:** BGSU unit IDs are `PDB|Model|Chain|ResName|ResNum|Atom|AltId|InsCode|SymOp`, with trailing fields optional.
- **What the 42,404 endpoint occurrences contain:**
  - 41,760 have 5 fields; 644 have 7 (an alternate location in the AltId field).
  - All are model 1.
  - None has an insertion code or a symmetry operator, so no `P_1` or `1_555` handling was needed here (it was needed for 5AFI in RF00177).
- **Join:**
  1. `(PDB, auth chain)` → instance. No chain was ambiguous.
  2. `(author number, insertion code)` → mapped residue → `polymer_seq_position` (label_seq_id) → Rfam state.
- **Check against the mmCIF:** each mapped endpoint was compared with `_pdbx_poly_seq_scheme` (label asym, seq_id, mon_id) and `_atom_site` (residue present in model 1).

### Endpoint outcomes

Counts are over all 42,404 endpoint occurrences. BGSU lists every interaction in both orders, so each endpoint is counted twice.

| Outcome | Endpoint occurrences | Example |
|---|---|---|
| Consensus (match-state) position | 37,382 | `1Y26\|1\|X\|U\|74` → 1y26_1_X_A label_seq 62 → consensus 80 |
| Insertion (between consensus positions) | 98 | `9IWF\|1\|A\|A\|61` → label_seq 60 → insertion 77^2 (second residue after 77; xanthine-II) |
| Mapped residue outside the aligned region (construct flank) | 4,130 | `1Y26\|1\|X\|C\|13` → label_seq 1, `not_submitted_flank` |
| Not a polymer residue: ligand, ion or other heteroatom in the RNA chain | 794 | `1Y26\|1\|X\|ADE\|90` (adenine, cWW with U74) |
| Other model, symmetry mate, other polymer chain, ambiguous | 0 | none |

- **Verification:** 41,604 of the 41,610 mapped endpoint occurrences agree with the mmCIF in label asym, seq_id and residue name.
- **The 6 that don't:** all are 1Y27 C81. The mmCIF has `CCC` (cytidine 2′,3′-cyclic phosphate) and BGSU writes `C`. Position and chain agree, so the mapping is correct. This is a residue-name difference for a modified nucleotide.
- **Modified polymer residues:** these map normally, for example 5′ `GTP`, `IU2` / `IU3` and `50N`.
- **Nothing was discarded.** Non-polymer and flank endpoints are counted above. Insertion and flank interactions are in the mapped table with their state, but are not aggregated by consensus pair.

## B. Interaction dataset

### Rows and normalisation

- **Rows:** 21,202 (5,668 base-pair rows, 15,534 stacking rows).
- **Both orders listed:** every row has its mirror row with the class mirrored (cWH ↔ cHW, s35 ↔ s53; cWW, s33 and s55 are symmetric).
- **Normalisation:** interactions are written 5′→3′ by label_seq_id, flipping the class when the row was written from the 3′ nucleotide. Alternate locations are then merged, so FR3D's orientation-dependent classes stay correct.
- **Grouping after normalisation:**
  - 20,080 intra-chain rows give **9,991 interactions** in 84 chains.
  - 332 rows pair two RF00167 chains: the 9NBH and 9NDD scaffold assemblies (328) and 6VWV (4). These are kept in `inter_chain_interactions.csv`, not aggregated.
  - 790 rows involve a ligand or ion. They are out of scope, and the ligand contacts already come from PDBe.

### The 9,991 intra-chain interactions

| Category | Consensus–consensus | With an insertion or flank | Near (stacking only) |
|---|---|---|---|
| Base pairs | 2,397 | 305 | n/a |
| Stacking | 5,817 | 701 | 740 consensus, 31 other |

### Aggregated pairs

- 199 distinct (category, consensus pair) combinations are observed.
- With at least 10 eligible chains: 45 base-pair and 119 stacking pairs.

### Eligibility (the denominator)

- **Eligible:** both consensus positions are match-state residues of the chain, both bases are complete (all FR3D standard heavy atoms in model 1; ring atoms for modified residues), and the chain has FR3D annotations. All 84 chains have them.
- **Breakdown:** the conservation table splits the rest into:
  - observed;
  - near only;
  - both present but not annotated;
  - position absent from the chain (deleted or outside the alignment);
  - unmodelled;
  - base incomplete;
  - annotation unavailable.
- **Interactions annotated on ineligible bases:** kept as `n_observed_but_ineligible`, not added to the frequency.
- **Units:** frequencies count chain instances. PDB-entry support is given next to them to show redundancy, for example 80 chains from 65 entries.

## C. Family-wide conservation

### Most conserved

Base pairs with at least 75 eligible chains and frequency of 0.9 or more, outside `SS_cons`:

| Pair (consensus) | Residues | Context | FR3D class | Chains | Entries |
|---|---|---|---|---|---|
| 41–66 | G37–C61 | L2–L3 kissing loop | cWW | 80/80 | 65/65 |
| 42–65 | G38–C60 | L2–L3 kissing loop | cWW | 75/80 | 63/65 |
| 26–56 | U22–A52 | J1/2 × J2/3 | cWW | 73/80 | 62/65 |
| 54–81 | C50–U75 | J2/3 × P1 | cWS | 71/79 | 61/64 |

- **`SS_cons` pairs:** all are observed, mostly at 0.95 or more. The exceptions are explained in §E: 59–77 (G55–C71) at 0.79, 30–48 (C26–A44) at 0.84 and 19–87 (G15–C81) at 0.86.
- **Stacks at 0.9 or more:** 43 consecutive stacks (i, i+1), plus long-range stacks such as A59/G68 (63–74) and G57/C70 (61–76) inside P3.

### Variable, present in some comparable structures

Explanations come from §E and §F.

| Pair | Residues | Context | f (chains) | Who has it / who lacks it |
|---|---|---|---|---|
| 53–82 | U49–A76 | J2/3 × P1 | 0.71 (56/79) | Absent in 5E54, 5SWD, 6VWT and 6VWV; 10/17 in guanine-II |
| 39–69 | A35–A64 | L2–L3 | 0.65 (49/75) | 39/51 in group 1, 3/17 in guanine-II |
| 42–71 | G38–A66 | L2–L3 | 0.66 (53/80) | 40/51 in group 1, 5/17 in guanine-II, 0/4 in 9NDD |
| 50–57 | G46–C53 | within J2/3 | 0.82 (66/80) | 4/17 in guanine-II |
| 57–78 | C53–G72 | J2/3 × P3 | 0.16 (13/80) | Only guanine-II and xanthine-II (13/17) |

### Class changes

FR3D classes are kept verbatim.

| Pair | Residues | Classes | Cause |
|---|---|---|---|
| 29–49 | U25–A45 | cWW 75, tHS 4 | The 4 tHS chains (4FEL, 4FEN, 4FEO, 4FEP) are U25A mutants, so the pair is A·A |
| 59–77 | G55–C71 | cWW 59, tWW 4 | The 4 tWW chains (9IWF, 9IWG, 9LKE, 9LKF) have A at that position: A41–C59 in guanine-II and xanthine-II numbering |
| 38–70 | U34–A65 | tWH 63, tWW 1 | 9V50, a single 3.25 Å chain |

Class changes in this family therefore follow nucleotide identity. A class-change view must show the nucleotides.

### Uncommon or restricted

These come from 1–4 chains of 1–2 entries each:
- U47–A52, U49–A73 and U48–U74 / U48–U75: only 5E54 and 5SWD;
- A33–G37 cSS: 3 chains;
- G38–A66 and C61/G62 stacks: only 9NDD.

None of these should be read biologically without the entry context in §F.

## D. Where variation occurs

Pairs with at least 10 eligible chains, in frequency bands.

**Base pairs:**

| Region | f ≥ 0.9 | 0.1–0.9 | < 0.1 |
|---|---|---|---|
| Helix pairs (P1–P3) | 18 | 4 | 0 |
| L2–L3 loop–loop | 2 | 4 | 0 |
| Involving J1/2, J2/3 or J3/1 | 1 | 5 | 8 |
| Involving a ligand-contact position* | 2 | 1 | 6 |

**Stacks:**

| Region | f ≥ 0.9 | 0.1–0.9 | < 0.1 |
|---|---|---|---|
| Adjacent (i, i+1) | 43 | 13 | 6 |
| Within helices | 5 | 10 | 1 |
| Involving a junction | 1 | 9 | 15 |

\* Ligand-contact positions are those contacted by the entry's PDBe-annotated ligand in at least 50% of ligand-bound entries: consensus 25, 26, 51, 55, 56, 80 and 81 (A21, U22, U47, U51, A52, C74/U74, U75).

**Patterns:**
- **Helices:** base pairs are invariant. Their variable stacks are mostly cross-strand stacks in P1 and P2 (e.g. 21–86, 22–85, 32–47), which depend on small helical-geometry differences. These need resolution-aware interpretation.
- **Junction:** most of the variation is here. The three-way junction around the ligand pocket holds most rare and variable base pairs and stacks.
- **Hotspots:** the positions in the most variable interactions are G46 (50), C53 (57), A66 (71), G72 (78), G55 (59) and A24 (28). They lie in J2/3, at the base of P3 and in L3.
- **Kissing loop:** the canonical pairs are invariant, while the noncanonical loop–loop pairs (A33–A66, U34–A65, A35–A64, G38–A66) vary. That variation is mostly between subclasses (§E).

## E. Structural groups

Groups come from the existing clustering, without recomputing anything:
- **Geometric discrepancy:** 72 + 4 chains.
- **TM-score:** 51 / 17 / 4 / 4 chains.

Each group is compared with the rest of the clustered chains using eligible denominators. Fisher exact p-values are in the CSV as a descriptive aid only: the groups are small and their chains are not independent.

**TM group 2 (17 chains):** guanine-II (*Paenibacillus*, *S. aureus*; 9LJN, 9LKC, 9LKE, 9LKF, 9LMC, 9V4U–9V51) and xanthine-II (*P. beijingensis*; 9IWF, 9IWG).
- **Has, unlike the other groups:** C53–G72 cWW (13/17) and the stacks A24/G55 (17/17), C53/A73 (13/17) and G72/A73 (17/17).
- **Mostly lacks:**
  - G46–C53 (4/17) and G55–C71 (4/17);
  - the L3 pair G62–U63 cSH (0/17) and the noncanonical kissing-loop pairs (A35–A64 3/17, G38–A66 5/17);
  - the P3 stacks G55/G72, C71/G72 and A66/U67 (0/17).
- **Interpretation:** J2/3 and the base of P3 pair in a different register. This fits these chains using the 64–72 consensus region (pair 64–72 in 16/17) and inserting after position 77 (the 77^n insertions). It is a subclass and sequence difference, not a conformational state of one RNA.

**TM group 3 = geometric-discrepancy group 2 (4 chains):** 5E54 A/B (apo) and 5SWD A/B (intermediate, adenine-bound). This is the only between-group difference that is not explained by sequence (§F).
- **Have only here:** U47–A52 cWW, U49–A73 cWW (4/4 vs 0/72), and the stacks G46/U47, U47/C50 and C50/A52 (4/4 vs 0/68–72).
- **Lack:** U22–A52 (0/4 vs 69/72), C50–U75 (0/4 vs 68/72), A23–G46 (0/4 vs 66/72), U49–A76 (0/4 vs 53/72), and the stacks C50/U51, U51/A52, U22/A23 and A23/U25.

**TM group 4 (4 chains, 9NDD A–D):** the cryo-EM RNA scaffold carrying an 8-oxoguanine aptamer.
- **Lacks:** the canonical kissing-loop pair G38–C60 (0/4), the P3 pairs C58–G68 and G55–C71 (0/4), and C50–U75 / U49–A76 (0/4).
- **Has:** G38/A66 and C61/G62 stacks seen nowhere else.
- **Interpretation:** 9NBC has the same sequence, yet its single chain has these pairs. The difference may come from the construct or the map rather than biology; it was not investigated further.

**Representatives** (`RF00167_representative_comparison.csv`):
- **6UC8 B (group 1) vs 5E54 B (group 3):** 86 shared interactions, 12 only in 6UC8 B and 19 only in 5E54 B where both chains are eligible.
  - 3 more differences are not comparable: two because a 5E54 B residue is unmodelled, one because a position is absent.
  - The eligible differences are the junction differences above, and the other three 5E54/5SWD chains carry the same differences.
- **Comparisons with 9NDD A:** dominated by positions absent from the other chain (`position_not_in_chain`; scaffold-flank consensus positions). The table marks such cases instead of counting them as gains.

## F. Ligand-associated variation

**Comparison design:** only identical sequences were compared (`sequence_key`), so organism, mutations and construct boundaries can't explain a difference. Five sequences have more than one ligand state.

### 1. *add* adenine riboswitch, *V. vulnificus*: 13 chains of one sequence

`RF00167_vvulnificus_junction_states.csv`, with 1Y26 (a different construct of the same riboswitch) added as a reference:

| Chain(s) | Ligand | U22–A52 | U49–A76 | C50–U75 | A23–G46 | U47–A52 | U49–A73 | C50 on U47–A52 | A23–U74 | U48–U75 | U48–U74 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1Y26, 4TZX, 4TZY, 4XNR, 5SWE, 5UZA | ADE | cWW | cWS | cWS | tWS | – | – | – | – | – | – |
| 6VWT A/B, 6VWV A/B | ADE | cWW | **–** | cWS | tWS | – | – | – | – | – | – |
| 5E54 B, 5SWD B | none / ADE | – | – | – | – | cWW | cWW | s55, s35 | tHW | cWS | – |
| 5E54 A | none | – | – | – | – | cWW | cWW | s55, s35 | – | – | cWW |
| 5SWD A | ADE | – | – | – | – | cWW | cWW | s55, s35 | – | – | – |

"–" means both bases are complete and annotated, but no interaction was found.

**Reading the table:**
- **Clean split:** the bound-state junction network is present in every final adenine-bound chain. The apo-type alternative pairs are present in every 5E54 and 5SWD chain, and in no other chain.
- **The split is not simply bound vs unbound:** 5SWD is adenine-bound but has apo-type junction pairing. That is the intermediate state of a time-resolved experiment, so the interaction pattern follows the reaction state, not the presence of the ligand.
- **6VWT and 6VWV:** these crystal phase-transition structures (Ramakrishnan et al. 2021) are bound but lack U49–A76. This is an observation only; the paper was not read at nucleotide level.

### 2. Guanine-binding sequences with many purine analogues

- ***xpt* (B. subtilis), 10 chains with 9 ligands:** none of the variable interactions follows ligand identity. They are P1–P3 stacks and C26–A44 (6/10).
- **Guanine-II (Paenibacillus), 13 chains with 11 ligands:** the same holds. 9LMC (Mn²⁺ soak) and 9V50 (3.25 Å) lack the most interactions, which suggests a resolution or soaking effect rather than a ligand effect.

### 3. 8-oxoguanine scaffold: 9NDD (4 chains) vs 9NBC (1 chain), both bound

See §E. Same ligand, but different construct and map; no ligand conclusion.

### 4. Other apo entries

- **3IVN:** a *pbuE* U65C mutant with no ligand. It has no bound counterpart of identical sequence here.
- **9NBH:** an aptamer-core scaffold. Its chains are not clustered.

**Conclusion:**
- **Evidence:** in this dataset, ligand-associated variation in RNA–RNA interactions is supported only by one identical-sequence comparison, the *add* XFEL series. There it matches the published description of the reaction states (§G).
- **Other ligand differences:** among guanine riboswitches bound to different purine analogues, no RNA–RNA interaction difference is associated with the ligand at the level FR3D annotates.
- **Not established:** the data don't show that the ligand causes any of these changes. Crystal form differs between 5E54/5SWD and the other entries.

## G. Literature comparison

| Publication | PDB in RF00167 | Reported observation | Reproduced here? |
|---|---|---|---|
| Serganov et al. 2004, *Chem Biol* 11:1729, doi:10.1016/j.chembiol.2004.11.018 | 1Y26 (A-switch + adenine), 1Y27 (G-switch + guanine) | Tuning-fork fold; ligand pairs with U74 (A) or C74 (G); loop–loop interaction | Yes (structure-level). Kissing-loop G37–C61 and G38–C60 cWW in 80/80 and 75/80 chains. Ligand–U74 cWW is in the BGSU rows (e.g. `1Y26\|1\|X\|ADE\|90` cWW U74); RNA–ligand pairs were not aggregated |
| Batey, Gilbert & Montange 2004, *Nature* 432:411, doi:10.1038/nature03037 | 4FE5 (*xpt* + hypoxanthine, re-deposition of 1U8D) | Binding pocket formed by the three-way junction, which nearly encloses the ligand | Partly. The junction pairs around the ligand-contact positions (U22–A52, C50–U75, U49–A76) are present in 4FE5 and in the large majority of bound chains |
| Stagno et al. 2017, *Nature* 541:242, doi:10.1038/nature20599 (full text read, PMC5502819) | 5E54 (apo1 and apo2), 5SWD (intermediate), 5SWE (bound) | In both apo forms, G46, U47 and U49 form Watson–Crick pairs with C53, A52 and A73, and C50 stacks on U47–A52; A73·(A52–U22), U49·(U20–A76), C50·(U75–A21) and the U51 triple are absent; apo1: A23–U74 Hoogsteen and U48–U75 wobble; apo2: U48–U74 wobble | **Yes, chain by chain** (§F table): U47–A52 cWW and U49–A73 cWW in all four 5E54/5SWD chains, none elsewhere; C50 stacks on U47 and A52; U22–A52, U49–A76 and C50–U75 absent. 5E54 B has A23–U74 tHW (Hoogsteen) and U48–U75, as apo1 does; 5E54 A has U48–U74, as apo2 does. G46–C53 is also present in the bound chains here, so it does not separate the states |
| Ramakrishnan et al. 2021, *Nat Commun* 12:1762, doi:10.1038/s41467-021-21838-5 | 6VWT, 6VWV (adenine riboswitch phase transition in crystals) | Ligand-triggered cooperative rearrangement in the crystal | Not tested. U49–A76 is absent in all four chains; the abstract does not describe this |
| Li et al. 2025, *Nucleic Acids Res* 53, doi:10.1093/nar/gkaf884 | 9LJN, 9LKC, 9LKE, 9LKF, 9LMC, 9V4U–9V51 (guanine-II) | Conserved scaffold shared with guanine-I; local rearrangements adapt to ligands | Partly. Overall architecture is conserved (helix and kissing-loop canonical pairs present). The distinct J2/3–P3 pairing (C53–G72; no G46–C53 or G55–C71) is a subclass difference. Ligand-specific rearrangements are not visible at FR3D pair level (§F2). Not checked against the paper at nucleotide level |
| Xu et al. 2025, *Sci China Life Sci* 68:2073, doi:10.1007/s11427-024-2800-0 | 9IWF, 9IWG (xanthine-II) | Guanine-riboswitch-like three-way junction; a U→G change and insertions create a xanthine-specific pocket | Consistent. The insertion residues map as Rfam insertions (77^n), and the chains group with guanine-II (TM group 2) |

**Not testable with this dataset:**
- adenine vs guanine riboswitches beyond position 74, because their sequences differ in many places;
- a guanine riboswitch with and without ligand from the same sequence. No apo xpt structure is present; the only other apo entries are 3IVN and 9NBH (§F4).

The research worktree's `literature_review.md` covers R3DCID, which compared 35 purine-riboswitch chains by eye. This analysis is a systematic version of that comparison.

## H. Implications for an Interactions tab

What the data support, smallest useful version first:

1. **Family-wide interaction table on consensus pairs.** Supported.
   - Columns: pair (both positions, with residue labels in the selected numbering), category, FR3D classes with counts, `k / n eligible`, entries, and why the rest are not counted (missing, unmodelled, incomplete base, unannotated).
   - Filters: region / context and `SS_cons` vs tertiary.
   - **Why:** it immediately separates the invariant core (helices, kissing loop) from the variable junction. Class changes must show the nucleotides (§C).
2. **Instance-vs-instance comparison** (shared, only in A, only in B, reclassified, with a "not comparable" state when the partner position is missing). Supported, and the best way to find biology.
   - 6UC8 B vs 5E54 B immediately shows the apo-junction rearrangement. The same comparison is how the Stagno observations were reproduced.
   - It must check eligibility before calling something "lost", as the representative table does.
3. **Group comparison** (frequency in group vs rest). Useful but easy to misread:
   - here every strong group difference is subclass, construct or reaction state;
   - group sizes are 4–17 chains from few entries.
   
   Show `k/n` and entry counts, never p-values alone, and show the group's sequence and ligand composition beside it.
4. **Mapping to 3D.** Feasible now. Every interaction keeps `seq_name`, label_seq_ids and the original unit IDs. Highlighting both nucleotides of a selected pair in the existing viewers uses the same identifiers as the current selection.
5. **Ligand context.** Useful only when comparisons are restricted to identical sequences. The table should be able to filter or split by ligand and by `sequence_key`; across sequences, ligand differences are confounded.

**Not supported yet:**
- automatic "enriched interaction" calls;
- near base pairs (not in the BGSU CSV export);
- aggregation of insertion-mediated interactions (only 49 interactions).

## I. Limitations and next steps

**Limitations:**
- **Redundancy:** chains are not independent. 84 chains come from 66 entries, many of the same sequence, and there are multi-chain entries (9NBH and 9NDD with 4 chains each). Instance frequencies with entry support are given; no weighting was applied.
- **Coverage:** frequencies are FR3D annotation frequencies. A missing interaction in a 3.0–3.25 Å chain (e.g. 9V50, 9LMC) may be a resolution effect. Resolution was not modelled.
- **Unchecked differences:** the 9NDD differences were not investigated (construct or map quality).
- **Not aggregated:**
  - insertion interactions (49) and flank interactions (988) are kept but not aggregated;
  - guanine-II's P3 register difference is partly in insertion columns, so it is understated by consensus-pair statistics.
- **Not done:**
  - base triples are visible only as pairs; no triple reconstruction was attempted;
  - RNA–ligand FR3D pairs (790 rows) were not analysed;
  - the ligand-contact positions come from PDBe annotations, not FR3D.
- **Data source:** BGSU data were fetched on 2026-10-18. Newer BGSU releases may differ.

**Next steps:**
- **Weighting:** add a per-entry (or per-sequence) frequency to show redundancy.
- **Resolution:** split variable interactions by resolution band.
- **Literature:** read Li et al. 2025 and Ramakrishnan et al. 2021 at nucleotide level and test their specific statements.
- **Scope:** decide whether insertion columns (guanine-II P3) should be aggregated by anchor.
- **Production:** if the tab goes ahead, build the observation table in the backend (`rna_explorer`) from these two stages, with schema and tests. This report is the method reference, not production data.
