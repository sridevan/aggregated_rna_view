# Near base pairs and the Rfam consensus comparison: validation (2026-10-22)

**Scope.** This implements the design agreed in `experiments/sscons_vs_3d/NEAR_CHECK.md` §8–9:
- base pairs only, with near base pairs included;
- no helix information;
- an "Rfam consensus" column with Both / 2D / 3D in the Interactions table and arcs.

The method and environment are in `docs/NEAR_BASEPAIRS_METHOD.md`.

## What changed

| Layer | Change |
|---|---|
| Backend | A new command `rna_explorer near <ACC>` (`rna_explorer/near_basepairs.py`) writes `near_basepairs.json` (schema `near_basepairs`) from a pinned local fr3d-python run (commit 3feb9b9). It streams each mmCIF into a temporary directory and never stores it. |
| Backend | `interactions` (schema 1.4) merges the near rows, tagged `fr3d_local`. It adds per pair: `consensus`, `n_paired`, `frequency_paired`, `n_entries_paired`, `near_classes`; and per chain: `near_basepairs`. Every SS_cons pair now gets a base-pair row. |
| Frontend | A new "Rfam consensus" column (sortable) and a matching filter under More filters. Near base pairs count in Conservation and appear in the ▸ class breakdown. A new per-chain state "Near" (details, viewer caption, Compare two). "2D" arcs are grey and dashed, with a legend row. Helix coverage counts near pairs. |
| Docs | `docs/TOOLTIPS.md` (new and changed texts) and `docs/NEAR_BASEPAIRS_METHOD.md`. |

## Checks

**BGSU-derived numbers are unchanged in all six families.** Comparing every existing pair with the committed `interactions.json`:
- the following fields are identical: `n_observed`, `n_eligible`, `n_entries_observed`, `classes`, `n_missing`, `n_incomplete`,
  `n_unmapped` and `n_annotation_unavailable`;
- chain annotation availability and `annotation_gaps` are unchanged, including the 7PJS a gap;
- BGSU file statistics are unchanged; the local rows have their own `summary.near_basepairs_local`.

**Near endpoints all verify against the mmCIFs.** No endpoint mismatched in any family; RF00177 alone has 30,688 verified endpoints.

**Pooled support reproduces the investigation.** Pooled = all consensus base-pair observations, paired / eligible:

| Family | Consensus pairs: Both / 2D | Pooled true (any class) | Pooled paired (true or near) | Investigation (NEAR_CHECK §7) |
|---|---|---|---|---|
| RF00167 | 22 / 0 | 0.953 | 0.987 | 0.987 |
| RF00162 | 32 / 1 | 0.932 | 0.960 | 0.960 |
| RF00059 | 22 / 6 | 0.887 | 0.908 | 0.908 |
| RF00005 | 21 / 0 | 0.920 | 0.971 | 0.973* |
| RF00011 | 97 / 0 | 0.898 | 0.938 | 0.938 |
| RF00177 | 443 / 19 | 0.906 | 0.951 | 0.951 |

\* The investigation also counted 2 tRNA observations that only the local run annotates as true pairs. Production does not use local true
pairs, so the figure is 0.971.

## What "2D" contains

**RF00059:**
- 64–83: 0 of 23 eligible chains. Position 64 pairs with 82 instead, as a G·G cWH.
- 38–48 to 42–44: no chain can be evaluated, because these positions are unmapped in every chain. They are labelled "not evaluable in
  any chain".

**RF00162:** 86–92, with 1 eligible chain and 28 unmapped.

**RF00177 (E. coli only):**
- Paired in no eligible chain (59–60 eligible each): 13–21, 372–385, 596–622, 934–1228, 950–965, 951–964, 1015–1026, 1043–1196,
  1054–1183, 1149–1169, 1150–1168 and 1497–1520.
- Not evaluable in any chain (0 eligible): 178–201 to 184–195.

**Caveat on "Both".** "Both" means paired in at least one chain, so weakly supported consensus pairs carry "Both" and show their
disagreement only in the frequency. Examples:
- SAM-I 50–64: 3 of 28 chains;
- TPP 21–27: 9 of 34;
- RNase P 326–340: 16 of 50;
- 16S 572–748: 3 of 59;
- 16S 648–737: 1 of 59.

Sorting by Conservation within the Both rows brings these to the top. This is the trade-off accepted in the design discussion.

## Tests

**Backend.** All 114 pytest pass. New or changed tests in `tests/test_interactions_rf00167.py`:
- **`test_near_basepairs_file`:** checks the schema, the pinned commit, every entry annotated, near classes only, and family chains only.
- **`test_local_near_basepairs_in_observations`:** checks the source tags; local rows are near base pairs only; BGSU rows have no near base
  pairs; endpoints verified; provenance commit.
- **`test_consensus_status_and_paired`:**
  - every SS_cons pair has a row;
  - the Both / 2D / 3D definitions hold;
  - `n_paired` = `n_observed` + `n_near_only`;
  - 3D rows need a true pair;
  - RF00167 pooled support is 0.987.
- **The reproduction test:** stacking near counts are still compared exactly; base-pair near-only chains were "not annotated" before.

**Frontend unit tests (vitest).** All 119 pass. New:
- the near state, `isPaired`, and comparison of a near pair as "only";
- RF00167 consensus counts and `paired()`;
- the consensus filter.

Helix coverage expectations are updated: 9V50 A is now 19 of 20 (13 true + 6 near), and the family median is 100%.

**Playwright (`tests/interactions.spec.ts`):**
- Existing tests now use paired counts and colours.
- The coverage test now shows 9V50 A at 95% (no longer sparse) and flags tRNA 1J2B C (16 of 21, 76%).
- New: a near-base-pair test (counts, the ncWW class row, the "Near" state in 9V50 A, details, caption, provenance) and a consensus test
  (TPP 64–83 2D row and arc, the filter for each value, base-pair-only behaviour, the details row, sorting).
- Full suite: 211 passed (2026-10-22).

## Open points

1. **Licence.** fr3d-python's `license.txt` is MIT but its `pyproject.toml` declares Apache-2.0. Both are recorded in the provenance.
2. **7PJS a.** The local run annotates true pairs in the region BGSU leaves unannotated. These are not used; the gap stays "annotation
   unavailable".
3. **BGSU web interface.** It was unavailable (HTTP 503) on 2026-10-22, so it is not yet re-checked for published near base pairs.
