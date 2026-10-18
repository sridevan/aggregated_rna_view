# RF00005 (tRNA, 60 curated chains): family page build

Date: 2026-10-18.

The selection is the research tRNA demo set (`data/selections/RF00005_60_chains.csv`), copied unchanged from the read-only research worktree; see `data/selections/README.md`. Its structural groups were chosen by the user on 2026-10-17 (`RF00005_DEMO_CLUSTERING_RESULTS.md` in the research worktree).

## Phase 1: mapping

**Build:** `python3 -m rna_explorer build RF00005 --cpu 4 -v --selection data/selections/RF00005_60_chains.csv`.

| Measure | Value |
|---|---|
| Chains / entries | 60 / 59 |
| Unique mapped sequences | 43 |
| Chains passing validation | 60 / 60 |
| Skipped | 0 |
| Model length | 71 (Rfam excludes the 3′ CCA) |
| Matched positions per chain | 62–71, all with coordinates in the chain with the fewest |
| Insertions | 51 chains; at most 22 nt (9QRP C: type II variable arm after consensus 44, plus D-loop insertions after 16 and 19) |
| PDBe 2D topology | 46 of 60 chains |

**Agreement with the research:** `compare_with_research.py` shows all 60 chains identical to the research build `data/rf00005_demo`, for every residue and every consensus column.

**Representative:** the automatic rule picks **1EHZ A** (yeast tRNA-Phe, free, X-ray 1.93 Å). This is also the research's reference chain, so no override is needed.

**Robustness fix:** the PDBe API returned intermittent HTTP 500s (1EHZ, 5CCB). The download cache now retries server errors (5xx) as well as transport errors, up to 3 attempts. Client errors (4xx) are not retried. Test: `tests/test_selection_build.py::test_cache_retries_server_errors_not_client_errors`.

## Phase 2: orientation, SASA, modified nucleotides, annotations

**Orientation** (whole-chain superposition on 1EHZ A): 17 good, 34 moderate, 9 poor.

The poor fits are the structurally unusual tRNAs from the research:
- 1J2B (λ-form);
- the mitoribosome central-protuberance tRNAs 7QI4 JA and 8OIN I;
- the mitochondrial TRMT10C–SDR5C1 tRNAs 7ONU, 8CBM and 9GCH;
- 3QSY, 5CCB and 9Q3R.

**SASA:** 60 of 60 chains measured; 4,481 nucleotides complete, 18 incomplete, 149 without coordinates. All 71 positions have a value.

**Nonstandard nucleotides:**
- 182 of 4,648 polymer residues, in 30 chains (29 entries) and 27 distinct components.
- 170 lie on the consensus axis, at 26 positions. The other 12 are in insertions (10) or flanks (2).
- The usual tRNA modification sites are covered, including 4SU8, D (H2U) in the D-loop, m1A58, T54 (5MU), Ψ55 (PSU), m2G/m2,2G and the anticodon-loop modifications (MIA, T6A, QUO, YYG).

**Annotations:**
- **Loops:** 295 (HL 210, IL 17, J3 1, J4 62, J5 5). BGSU labels them as anticodon loops, D-loops, T-loops and tRNA junctions.
- **Ligands:** 52 bound-molecule records (22 ligand, 22 metal ion, 8 additive). The ligand rows are SAH (2 entries), GAI, GTP, MQH (bottromycin A2) and SAM.

## Fix: junction loops beyond J4 (affects RF00059 and RF00177 as well)

**Problem:** the aggregated junction track (Annotations "JL", `tracks.json` `junction_loops_union`) counted only J3 and J4 loops. BGSU also reports higher-order junctions:

| Family | J5 loops |
|---|---|
| RF00005 | 5 |
| RF00059 | 1 |
| RF00177 | 178 |

They were in `loops.json` but missing from the track.

**Change:** the union now covers every junction type present (J3, J4, J5 and any higher order). It records `loop_types` and a per-type `by_type_n_instances`, and the schema is extended. The Annotations details list each type; the Family Overview names any Jn ("Five-way junction loop").

Tests:
- `tests/test_annotations.py::test_j5_junction_loops_are_in_the_junction_union`;
- `src/lib/loops.test.ts`.

**Effect:**

| Family | Change |
|---|---|
| RF00177 | 53 consensus positions (e.g. 36–39, 58–61, 99–102) previously showed no junction loop and now show junction loops in 60 of 60 chains. HL, IL and ligand tracks are unchanged. |
| RF00059 | Track unchanged; its J5 has no consensus-mapped nucleotides. |
| RF00167, RF00162 | No J5. |

All five families were re-annotated. RF00167, RF00162 and RF00059 loop and ligand records are identical to the committed versions. Their files gain only the new fields and fields from earlier approved changes (per-chain loop files; PUT listed as an additive).

## Fix: modified nucleotides broke the 3D cartoon

**Problem:** 3Dmol draws a continuous nucleic-acid cartoon only through residues named A, C, G or U. Modified nucleotides (PSU, 2MG, H2U and others) therefore broke the backbone tube in the Structures tab. tRNAs are heavily modified, so this was visible there, and it also showed in the homepage thumbnails.

**Change:** `drawableNucleotides` gives each modified nucleotide its parent base name, inferred from its base atoms, in the drawing copy only. Hover labels keep the deposited names.

Test: `src/components/structures/AlignmentViewer.test.ts`. All five thumbnails were regenerated.

## Pockets: not shown (user decision, 2026-10-18)

fpocketR ran on the 60 chains in 376 s:

| Result | Chains |
|---|---|
| With pockets | 23 (42 pockets) |
| No pockets | 34 |
| Failed | 3 (8CCS Bb, 8CDR Bb, 7ZRZ ZN1) |

The failures are author chain ids too long for the PDB format used as fpocketR input.

The most frequent pocket positions are the elbow (consensus 20 and 57), in 13 of 54 chains.

**Decision:** not shown, as for 16S. tRNAs have no small-molecule binding site, and most chains are taken out of protein or ribosome complexes. The outputs are kept in `experiments/fpocketr_full_RF00005/not_shown/`. The page shows pockets as not available.

## Structural groups (user decision, 2026-10-17)

Imported with `python3 -m rna_explorer clustering RF00005` from the research run T60 (60 chains, 43 shared positions). Both cuts are provisional, and the memberships are identical to the research `structural_groups_T60.json`.

| Metric | k | Groups |
|---|---|---|
| Discrepancy | 4 | 56 canonical; 7QI4 JA + 8OIN I (mitoribosome CP mt-tRNAs); 3QSY C (aIF2-bound initiator); 1J2B A (λ-form) |
| TM-score | 5 | the same, plus 5CCB C (TRMT6/61A-bound tRNA-Lys3) |

**Floor:** omitted, because there are only 2 duplicate pairs.

**Rationale:** recorded in `rna_explorer/structure_clustering_config.json`. The groups separate structurally unusual tRNAs from one canonical group; they are not binding-context classes.

## Front end and tests

**Front end:**
- RF00005 is added as "tRNA" to the family list, the home page (thumbnail of 1EHZ A, with a curated-selection note) and the family selector.

**Tests:**
- `tests/family-rf00005.spec.ts`:
  - Family Overview, including the 9QRP C variable-arm insertion marker;
  - Structures: groups, no floor, 1J2B alone;
  - Annotations: J5 in the junction details, modifications and ligands, no pockets.
- The home page test covers five cards.

## Representatives without backbone breaks (user request, 2026-10-18; all families)

**Rule:** a group's representative is now the medoid among the members whose chain has no internal coordinate gap, i.e. no unmodelled nucleotide between its first and last modelled nucleotide. Its 3D backbone in the Structures tab therefore has no breaks. Unmodelled terminal nucleotides do not count.

**Fallbacks and records:**
- If every member has gaps, the plain medoid is used, and the page lists its missing nucleotides.
- The plain medoid is always recorded as `medoid`.
- `representative_rule` and `representative_internal_gaps` are new schema fields.

**Unchanged:** group memberships, the cut and all statistics.

| Family | Metric, group | Plain medoid (internal gaps) | Representative now | Mean dissimilarity |
|---|---|---|---|---|
| RF00005 | discrepancy, 1 (56) | 9TNL Y (16, 17, 20, 47) | 6UFG B | 0.281 → 0.288 (+2.6 %) |
| RF00177 | discrepancy, 1 (50) | 8R8M GA (205–213, 841–845) | 4V9D AA | 0.0443 → 0.0457 (+3.2 %) |
| RF00059 | discrepancy, 1 (35) | 7TDC A (27–29) | 9C6K X | 0.130 → 0.144 (+11 %) |
| RF00059 | TM, 1 (35) | 7TDA A (27–29) | 4NYA A | 0.148 → 0.180 (+22 %) |
| RF00059 | both, apo group (2) | 8F4O A (both members have gaps) | unchanged | |

All other groups were unchanged; RF00167 and RF00162 had no change.

**Tests:**
- `tests/test_structure_clustering.py::test_internal_gaps_and_gap_free_medoid`;
- `tests/test_structure_clustering.py::test_rf00005_groups_and_gap_free_representatives`;
- `frontend/tests/family-rf00005.spec.ts`: the representative is 6UFG B and the medoid note names 9TNL Y.
