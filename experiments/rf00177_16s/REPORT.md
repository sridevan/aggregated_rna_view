# RF00177 (bacterial 16S rRNA, 60 E. coli chains) — family page build

## Phase 1: mapping (2026-10-17)

**Build:** `python3 -m rna_explorer build RF00177 --cpu 4 -v --selection data/selections/RF00177_60_chains.csv`. The `--selection` option is new: it maps only the listed chains, reads each chain's coordinates from the PDBe ModelServer, and rebuilds the validation scheme from the PDBe residue listing. The selection is the research's approved 60-chain set, copied unchanged; see `data/selections/README.md`.

| Measure | Value |
|---|---|
| Chains / entries | 60 / 58 |
| Unique mapped sequences | 15 |
| Chains passing validation | 60 / 60 |
| Skipped | 0 |
| Model length | 1,533 |
| Matched positions per chain | 1,506–1,516, of which 1,481–1,514 have coordinates |
| Insertions | 26 per chain; 27 for 9QQQ, 29 for 9G06 |
| PDBe 2D topology | 48 of 60 chains |

**Agreement with the research:** `compare_with_research.py` shows all 60 chains identical to the research build, for every residue (state, consensus position, insertion anchor and offset, coordinates) and every consensus column.

**Representative:** 4V9D Y, chosen by the user. It is recorded in `rna_explorer/family_overrides.json`; the automatic pick, 4V9O PG with its head rotated 10.7°, is kept in `automatic_choice`.

**Other code changes:**
- The download cache retries transient transport errors (3 attempts).
- Orientation, SASA and pockets read the per-chain files for selection-built families.
- Modified nucleotides need `_pdbx_poly_seq_scheme`, so they stream each whole-entry mmCIF: download, parse and delete, recording URL, size and MD5.

## Phases 2–3 (2026-10-17)

**Orientation** (whole-chain C4′ superposition on 4V9D Y): 28 good, 20 moderate, 12 poor. The poor and moderate fits are the rotated-head chains, as expected for a whole-chain superposition. A body-core frame, as used by the research, would separate head rotation; that is a possible later refinement.

**SASA:** 60 of 60 chains measured; 91,735 nucleotides complete, 7 incomplete, 619 without coordinates. 1,515 positions have a value.

**Modified nucleotides:** recorded at consensus positions including 506, 517, 956, 957, 1198, 1393, 1398, 1399, 1489, 1507, 1509 and 1510 (the m2G, m5C, m4Cm, m3U and m6,2A sites).

**Annotations:**

**Loops:**
- A new rule keeps loops with no nucleotide in a family chain only as compact records, because 70S entries list 23S, 5S and tRNA loops.
- 7,601 loops touch 16S and have full records; 20,326 other loops of the same entries are compact records.
- `loops.json` is 62 MB, down from 210 MB.
- RF00167, SAM-I and TPP have no loops outside their chains, so their data is unaffected.

**Named-element reference:**
- BGSU's family file uses 5J7L AA, which is not among the 60 chains, so its 1,534 identifiers are unresolved.
- The per-chain BGSU labels were evaluated for 60 of 60 chains.
- These labels are not displayed in the app.

**Ligands:** 107 bound-molecule records: 62 metal ions, 10 additives and 35 ligands. The ligands are 22 antibiotics plus putrescine (PUT, 6 entries).

**Pockets:** a one-chain probe (4V9D Y) ran fpocket in about 30 s, then fpocketR's analysis crashed inside its own code (`AttributeError: 'NoneType' object has no attribute 'getResnums'` in `analyze.add_basic_characteristics`). Not run for the family.

**Decisions (user, 2026-10-17):**
- **PUT reclassified:** putrescine is now an additive (`annotation_config.json` 1.1, changelog in the file). The RF00177 rerun has 29 ligand, 16 additive and 62 metal-ion records. No other family contains PUT.
- **No pockets for RF00177:** an isolated 16S without the 50S subunit and proteins gives mostly artificial surfaces, and fpocketR's analysis crashes on the 16S chain. The page shows pockets as not available for this family.

## Structural groups (user decision, 2026-10-17)

Both metrics use k = 2, marked "provisional". Groups are 50 and 10, the 10 being the most head-rotated chains (computed swivel 12.7–20.7°). The 10 chains are the same under both metrics.

| Metric | Gap at the cut | Silhouette | Cophenetic correlation |
|---|---|---|---|
| Discrepancy | 0.027 | 0.538 | 0.893 |
| TM-score | 0.019 | 0.727 | 0.858 |

**Inputs:** the research T1481 run of the 60-chain set (1,469 shared positions, none excluded).

**Floor:** not shown; the research computed no duplicate-pair floor for this run.

**Why provisional:** head rotation is a continuum, with gaps of similar size, and no bootstrap or jackknife was run. The groups describe head rotation, not functional state or antibiotic. Rationale: `rna_explorer/structure_clustering_config.json`. Candidates: `structural_groups/candidates.py`.

## Front end and tests

**Front end:**
- RF00177 is added as "Bacterial 16S rRNA" in the family selector and the 404 page.
- **Per-chain loop files:** `annotate` now writes `loops_by_instance/<seq_name>.json` (1.2–1.4 MB each for 16S), and the Family Overview loads only the selected chain's file. Families annotated before this change fall back to `loops.json`. `loops.json` stays the download.
- **Annotations axis:** the label step now follows the model length (10 for riboswitches, 100 for 16S); previously it drew 153 overlapping labels.
- **Download availability:** an item counts as available only when the response is JSON. Vite's dev server answers missing files with HTML.
- **Deploy copy:** `work/`, the per-chain `annotations/` files and the CSV tables are no longer deployed, because the app does not use them.

**Browser review** of every tab, no page errors:

| Item | Time |
|---|---|
| Page shell | 0.8 s |
| Mol\* (whole 70S entry) | 3.9 s |
| 2D topology | 3.9 s |
| Selected chain's loops | 0.3 s |
| Structures tab | 1.9 s |

**Tests:**
- `tests/family-rf00177.spec.ts`: 3 browser tests.
- `tests/test_structure_clustering.py::test_rf00177_head_rotation_groups`.
- `tests/test_annotations.py::test_loops_by_instance_files`.
- `tests/test_selection_build.py`: selection, scheme, retries, override, chain files, streaming.
