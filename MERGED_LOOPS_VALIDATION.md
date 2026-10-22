# Merged BGSU loops supersede their parts (2026-10-21)

## Decision

BGSU FR3D emits each loop "single", and additionally emits a **merged** loop where two loops share one canonical closing pair and are linked across it by a base pair or a stack. The best-known case is RF00167: J3_4FE5_001 = J4_4FE5_001 + HL_4FE5_002, joined across the isolated pair G46–C53. Until now the tracks counted both readings, so the J2/3 strand showed up as a hairpin in most chains.

The user decided (2026-10-21) that a **merged loop replaces its two parts**. An internal loop absorbed into a junction is acceptable; a junction region called a hairpin is not.

**Basis:**
- `~/Desktop/aggregated_rna_clustering/results/open_source_rna_motif_annotation_review/`;
- the local fr3d-python `-c loops` run, which labels the same loops `_merged`.

## Rule (`rna_explorer/annotations.py`: `merged_loop_relations`, `_mark_merged_loops`)

**Which loops count as merged:**
- The BGSU loop download and motifs page do not flag merged loops, so the relation is inferred within each PDB entry from nucleotide sets.
- Loop M is merged from A and B when A ∪ B = M and A ∩ B is exactly two nucleotides, i.e. the shared closing pair.
- Problematic loops take no part in the inference.
- In 4FE5 and 1Y27 this inference matches the local fr3d-python `_merged` labels.

**What each loop record stores:** `decomposition` = `{role: single | merged, merged_from: [[A, B]], superseded_by: [M], counted_in_tracks}`.
- Nested merges supersede every intermediate loop. For example, J3_9IWF_002 = IL_9IWF_001 + J3_9IWF_001, and J3_9IWF_001 = HL_9IWF_002 + J4_9IWF_001.

**Where it takes effect:**
- `tracks.json`: the HL, IL, J3 and J4 tracks and the junction union count only loops with `counted_in_tracks`. Each position lists the excluded loops in `superseded_loop_ids`.
- **Family Overview loop rows:** superseded loops are not drawn. The merged loop's tooltip reads "Merged from A + B".
- **Annotations tab, pinned details:** a line reads "{n} loops merged by BGSU into a larger loop not counted here."

**Nothing is removed:**
- `loops.json`, `loops_by_instance/` and `annotations/` keep every record.
- Residue-level loop lists carry `counted_in_tracks`.

**Schemas:** `loops`, `loops_instance`, `annotation_instance` and `tracks` were extended (new required fields).

## Effect per family

Regenerated with `python3 -m rna_explorer annotate <ACC> --offline` from the existing cache:
- Before the change, an offline rebuild of RF00167 reproduced the committed files exactly.
- After the change, `loops.json`, `ligands.json` and `reference_structural_elements.json` are identical to the previous outputs apart from the new fields and the generation timestamps.
- For RF00005 and RF00177, `ligands.json` was restored from the previous build. The offline run had only replaced the detail text of uncached PDBe 404s.

| Family | Loops | Merged | Superseded | Merges | Track cells whose instance count changed |
|---|---|---|---|---|---|
| RF00167 purine | 496 | 88 | 176 | J3 = HL + J4 (58), J3 = IL + J3 (30) | HL 8, IL 19, J4 15 |
| RF00162 SAM-I | 351 | 2 | 4 | J4 = IL + J4 (3IQR and its symmetry copy) | IL 5 |
| RF00059 TPP | 284 | 7 | 14 | IL = IL + IL | none (same positions, fewer IL observations) |
| RF00005 tRNA | 295 | 0 | 0 | — | none |
| RF00011 RNase P | 1,303 | 0 | 0 | — | none |
| RF00177 16S | 7,601 | 295 | 590 | J3 = IL + J3 (229), J5 = IL + J5 (56), IL = IL + IL (5), HL = HL + IL (3), IL = HL + J3 (2) | HL 8, IL 71, J3 17, junction union 17 |

### RF00167

**Loop composition per chain, from counted, non-problematic, single-chain loops:**

| Pattern | Chains |
|---|---|
| 2 HL + 1 J3 | 42 |
| 2 HL + 1 IL + 1 J3 | 30 (1Y27, 4FE5, 2B57, …) |
| 2 IL + 1 J3 | 7 (9NBH/9NDD scaffold constructs) |
| 3 HL + 1 J4 | 4 |
| 1 HL + 2 IL + 1 J3 | 1 |

- **The four 3 HL + 1 J4 chains** are apo 5E54 A/B and 5SWD A/B. They have no merged loop because 47–52 is a two-pair helix (U47–A52 stacked on G46–C53), which is a real conformational difference.
- **Position 53, loop records** (audit values unchanged): HL 66, IL 13, J3 76.
- **Position 53, tracks:**
  - HL 8 of 80 (apo 5E54/5SWD and 9NDD);
  - IL 0;
  - J3 76 of 80, now one J3 per chain (76 observations, previously 93);
  - superseded at this position: 58 HL, 13 IL and 17 J3.
- **Internal loops absorbed into the junction:** in the Guanine-II structures (9LJN, 9LKE…), 9IWF, 6P2H and 9V50, BGSU also merges a neighbouring internal loop into the J3. The rule counts it as part of the junction.

## Tests

- **pytest: 105 passed.**
  - New: `test_merged_loop_relations_inference` (nested merges; one or three shared nucleotides are not a merge).
  - New: `test_merged_loop_supersedes_parts_in_tracks` (synthetic family: records kept, tracks and union exclude the parts, schemas).
  - New: `test_merged_loop_decomposition` (RF00167 patterns above, 4FE5 relation, apo J4s).
  - Updated: the RF00167 position-53 and traceability tests.
- **vitest: 101 passed**, including new merged-loop cases in `lib/loops.test.ts`.
- **Playwright:**
  - The specs that hovered or selected now-superseded loops (HL_4FE5_002, J4_4FE5_001, IL_9IWF_001) use visible loops (HL_4FE5_001/003, IL_4FE5_001, J3_9IWF_002).
  - The specs that pinned the old position-53 hairpin tooltip ("66 of 80") now expect "8 of 80 (10%)".
  - New checks: superseded loops are not drawn, and the merged-loop tooltip line and the pinned-details note appear.
- **Playwright results:**
  - The first full run (191 tests) had 30 failures: the stale assumptions above plus timeouts in unrelated specs (structure clustering, single selection). That run took 7.1 h, likely because the machine slept.
  - Re-running all 64 affected tests, including those specs: all pass after the fixes.
  - The screenshot spec rewrites `docs/screenshots/`; those files were restored to their committed versions, so no screenshots were regenerated.

## Limitations

- **Inference, not BGSU's own label.** It agrees with fr3d-python where checked. A local fr3d-python run, or a `merged` column in the BGSU download, would replace it without changing the counts.
- **Problematic loops are excluded from the inference.** A merge whose part is problematic is not recognised, and both readings are then counted.
- **No pinned note where everything was superseded.** A position where every loop of a type was superseded shows no mark in that track, so the pinned-details note only appears where at least one counted loop remains.
