# RF00162 (SAM-I riboswitch) — Phase 1: Rfam mapping

Command: `python3 -m rna_explorer build RF00162 --cpu 4 -v`, run on 2026-10-17. Output is in `data/families/RF00162/` and is local only, not deployed. The model is Rfam RF00162 "SAM" (SAM riboswitch aptamer, S box leader), with CLEN 108.

## Result

| Measure | Value |
|---|---|
| PDB entries | 28 |
| RNA chains | 29 |
| Unique mapped sequences | 23 |
| Chains passing validation | 29 / 29 |
| Skipped | 0 |
| Log warnings | 0 |

- **Method:** all chains are X-ray, from 1.71 to 3.3 Å.
- **2D topology:** PDBe 2D is available for all 29 chains.
- **Rfam agreement:** for every chain, the cmalign span agrees with the Rfam structures endpoint and with the PDBe-mapped region.
- **Research worktree:** it reports the same 29 chains (T91 run, 86 shared positions).

## Composition

| Group | Chains | Mapping pattern |
|---|---|---|
| *T. tengcongensis* (Caldanaerobacter subterraneus subsp. tengcongensis) constructs and point mutants | 25 | 93 nt mapped, 91 match. Consensus 51–63 deleted, with a 2-nt local end (GA) at anchor 60; 41 and 86–88 also deleted. 4AOB A differs slightly: 93 match, no local end. |
| 3V7E C, D: *T. tengcongensis*, engineered helix P3 (YbxF complex) | 2 | 99 match. 26-nt local end at anchor 60, the engineered P3 module, plus 9 deletions. |
| 4KQY A: *B. subtilis* yitJ | 1 | All 108 positions matched, plus a 9-nt insert (GCGAAAGCC) after consensus 58. |
| 3IQN A | 1 | Same pattern as *T. tengcongensis*, but no species: neither the PDBe API nor the mmCIF gives one. |

## Points for the next phases

1. **Low diversity:** 27 of 29 chains are one organism's constructs and mutants, so the family statistics mostly describe *T. tengcongensis* SAM-I. The page should state this.
2. **Variable P3 region (consensus around 51–66):** it is deleted, replaced by a local end, or engineered in most chains. The structural comparison already excludes it (86 of 108 positions shared). The insertion markers will show local-end runs as dashed "+N" markers, which the front end already supports.
3. **3IQN A has no species:** it is kept as "not annotated by the source". It is not inferred from sequence identity unless a source states it.
4. **Representatives:** the display representative is 5FJC A (X-ray 1.71 Å, *T. tengcongensis*, 91 of 108 matched). The coverage representative is 4KQY A (all 108 positions, 3.02 Å, *B. subtilis*). Which one is the family reference for orientation is a choice for Phase 2.

# Phase 2: orientation, SASA, modified nucleotides (2026-10-17)

Commands: `orient RF00162`, `sasa RF00162`, `modifications RF00162`. All three exited cleanly; the outputs are in `data/families/RF00162/`.

**Orientation** (reference 5FJC A, the display representative, using the same rule as RF00167):

| Fit quality | Chains |
|---|---|
| good | 26 |
| moderate | 3 |

- Every chain is superposed on 83 C4′ pairs at PP `*`.
- The three moderate fits are 3V7E C at 2.44 Å, 4KQY A at 2.38 Å and 3V7E D at 2.35 Å. These are the engineered P3 construct and the only *B. subtilis* chain, as expected.
- The largest RMSD among the good fits is 1.81 Å.

**SASA:**
- All 29 of 29 chains were measured, with no failures.
- 2,807 nucleotides are complete and 8 incomplete; incomplete ones are left out of the median track.
- 5 alternate-location copies were dropped.
- All 108 positions have a median value.

**Modified nucleotides:** 2,815 polymer residues were checked, and 1 is non-standard. It is GTP at 5FJC A position 1 (consensus 1). It is flagged as 5′-terminal, it is the typical in-vitro transcription 5′ triphosphate, and PDBe agrees. No other source was unavailable.

# Phase 3: annotations and pockets (2026-10-17)

**Backend fix:** BGSU has no nucleotide-annotation file for RF00162, the named-helix reference; the family file returns 404. Its per-chain endpoint replies HTTP 200 with the plain text "Could not find annotation file RF00162_annotations.txt", which crashed the parser. The parser now skips lines without a separator. The per-chain fetch classifies that reply as `not_available` and any other unparseable reply as `parse_error`, and the builder records the chain as unavailable. RF00167 is unaffected. Three tests were added in `tests/test_annotations.py`; the full backend suite has 57 passing tests.

**Annotations** (`annotate RF00162`, live fetches):
- **Named-element reference track:** unavailable for all 29 chains, because the source has no file. The page must show this as "source unavailable", not as absence.
- **BGSU loops:** 351 loops, 269 of them with a motif group (17 distinct groups), and 9 flagged problematic. 0 identifiers are unresolved.
- **Bound molecules:** 81 records with 10 CCD ids, all classified by the existing generic config; no change was needed and no additives are reported.

| CCD | Class | Entries |
|---|---|---|
| SAM | ligand | 25 |
| SFG (sinefungin) | ligand | 1 |
| SAH | ligand | 1 |
| BA | metal ion | 19 |
| K | metal ion | 15 |
| NA | metal ion | 9 |
| MG | metal ion | 8 |
| MN | metal ion | 1 |
| IRI (iridium hexammine) | metal ion | 1 |
| NCO (cobalt hexammine) | metal ion | 1 |

27 of 28 entries contribute to the default ligand track.

**Chain attribution:** the 3 records of 3V7E (SAM, MG, NCO) are entity-level, because chains C and D are one entity. They are kept, with the attribution note.

**SAM-contact positions** (in at least 13 of 25 entries): 7, 11, 46, 48, 69, 70, 71, 103, 104.

**Pockets** (`pockets RF00162 --experiment-dir experiments/fpocketr_full_RF00162`):
- 29 chains in 208 s: 28 with pockets (55 in total) and 3GX5 A with none.
- **Top-ranked pocket per chain:** mainly consensus 24–25 (69–72 %), 76–80 and 100 (59 %). It covers the SAM-contact positions in only 0–21 % of chains.
- **Any pocket, score ≥ 0:** the SAM-site positions 48, 103 and 104 are covered in 52 % of chains.
- **Reading:** fpocketR's best-scored cavity is usually not the SAM site; the SAM site appears as a secondary pocket in about half the chains. This is reported as found and not tuned.

# Phase 4: structural groups (2026-10-17)

Cuts chosen by the user: geometric discrepancy k = 3 ("stable") and TM-score k = 4 ("provisional"). The rationales are in `rna_explorer/structure_clustering_config.json`. Candidates and evidence are in `experiments/rf00162_structural_groups/`. The output is `data/families/RF00162/structure_clustering.json`, which passes schema validation.

**Dataset:** 29 chains clustered, 0 excluded, 86 shared positions. Every chain has a match-state residue with coordinates at every shared position.

**Discrepancy, k = 3:**
- Gap 0.108, the largest in the tree; cophenetic correlation 0.972.
- Group sizes 26, 2 and 1. The pair is 3V7E C/D, the engineered P3; the singleton is 4KQY A, from *B. subtilis*.

**TM-score, k = 4:**
- Gap 0.057, the largest in the tree; cophenetic correlation 0.938.
- Group sizes 22, 4, 2 and 1. The 4-chain group is 2YGH, 4B5R, 5FJC and 5FK4.

**Floor:** null. The research has only two duplicate pairs for RF00162 and judges a floor from them not useful.

**Apo structure:** 3IQP A, the only ligand-free chain, stays in the core group at both cuts.

**Importer changes:**
- It now reads both research summary layouts: per-threshold for RF00167, per-run for RF00162. A missing floor is recorded with its reason.
- Rebuilding RF00167 in memory gives identical groups, cut, values, floor and evidence, so its file was left unchanged.

**Correction:** an earlier note said RF00162's `ufix_tm_user.csv` was missing. It is present; a truncated directory listing hid it. The candidate analysis read the TM values from `pairs_long.csv`, which equals that file for all 406 pairs, so the candidates stand. The temporary fallback was removed.

**Tests:** 2 new tests in `tests/test_structure_clustering.py`, one for the per-run evidence layout and one for the RF00162 groups. The backend suite passes with 59 tests.

# Phase 5: front end, tests and deployment (2026-10-17)

**Front end:**
- `src/lib/route.ts` now lists RF00162 as an available family, served at `/RF00162`.
- A family selector in the header loads the other family's URL.
- The 404 page links both families.
- No other front-end change was needed.

**Browser review** of `/RF00162`, every tab:
- No page errors. The only console message is the vendored 2D viewer's known double-mount message in development mode, which also appears with RF00167.
- Family Overview: 5FJC A is the representative. The local-end run GA shows as a dashed "+2" at 60|61, and 3IQN shows as "Unknown taxonomy".
- Structures: groups of 26, 2 and 1, the representative on top, and no floor line.
- Annotations: SAM, SAH and SFG ligand rows, the GTP row, loops, SASA and pockets.
- The named-helix labels are not displayed anywhere in the app, so their unavailability changes nothing on screen.

**Tests:** new `tests/family-rf00162.spec.ts` with 3 tests, Family Overview, Structures and Annotations, with expected values from the RF00162 JSON.
