# Family Overview Phase 3: insertion and deletion audit (RF00167)

Run `python3 experiments/overview_insertions_phase3/audit_insertions.py`. The script reads the per-chain mappings and changes nothing. It writes `audit_insertions.json`.

## How the existing mappings represent each state

Source files are `data/families/RF00167/mappings/<seq_name>.json`, written by `rna_explorer/mapping.py` from Infernal cmalign output. The mapping pipeline was not changed.

| Concept | Field |
|---|---|
| Consensus match position | A residue with `alignment_state = "match"` and `rfam_consensus_pos` from 1 to 102. The column also appears in `consensus_columns` with state `match`. |
| Rfam insertion | A residue with `alignment_state = "insert"` and `rfam_consensus_pos = null`. |
| Insertion placement | `insertion_anchor` is the consensus position before the insert column; 0 means before position 1, and the model length means after the last position. `insertion_offset` runs 1..n within the run. |
| Local end | `alignment_state = "local_end"`, with the same anchor and offset fields. |
| Deletion | A `consensus_columns` entry with state `deletion`, inside the chain's aligned span. |
| Uncovered | A `consensus_columns` entry with state `uncovered`, outside the aligned span. |
| Outside the mapped region | `alignment_state = "not_submitted_flank"`, with no anchor and no consensus position. These are terminal and construct residues outside the PDBe/Rfam-mapped region; they are not insertions. |
| Missing coordinates | `has_coordinates = false` on any residue. |
| Identifiers | `polymer_seq_position` (the `label_seq_id`), `author_residue_number`, `author_insertion_code`, `residue_name`, `nucleotide_parent`, plus the chain's `label_asym_id` and `auth_asym_id`. |

Every field the feature needs is already present. Nothing comes from gaps in author numbering or from author insertion codes.

## RF00167 results

84 chain mappings.

| Residue state | Count |
|---|---|
| match | 5,932 |
| insert | 34 |
| not_submitted_flank | 635 |
| local_end | 0 |

| Column state | Count |
|---|---|
| match | 5,932 |
| deletion | 163 |
| uncovered | 2,473 |

**Insertion runs:** 17 runs, one per instance, all the same event. Each is a 2-nt insert `CA` between consensus 77 and 78. The instances are:
- 9IWF A and B, and 9IWG A and B: label_seq_id 59–60, author 60–61.
- 9LKE A: label_seq_id 59–60, author 60–61.
- 9LJN A, 9LKC A, 9LKF A, 9LMC A and B, and 9V4U, 9V4V, 9V4X, 9V4Y, 9V4Z, 9V50, 9V51 A: label_seq_id 60–61, author 60–61.

All 34 inserted residues have coordinates, and all have PP character `6`. The following match is always consensus 78.

**Terminal insert states:** none. The 635 terminal residues are flanks, 562 of them with coordinates. They are reported as instance-only residues and get no marker.

**Consistency checks:**
- Run counts equal `summary.n_insertion_sites`, and run lengths match `summary.longest_insertion` for all 84 chains.
- Offsets are 1..n within each run.
- No inserted residue has a consensus position.

## Run rules (frontend `lib/insertions.ts`)

- **One run:** consecutive residues (`label_seq_id` + 1) with the same state and anchor and offsets 1, 2, …
- **Separate runs:** a change of state, a numbering gap or an offset restart starts a new run.
- **Markers:** one per boundary; several runs at the same boundary stay separate inside the marker.
- **Collisions:** markers whose labels would overlap at the grid scale are aggregated into one marker. This never happens in RF00167.
- **Boundary wording:** `left` = anchor (null when the anchor is 0) and `right` = anchor + 1 (null at the model end).
- **Synthetic tests only:** terminal, local-end, split and aggregated cases do not occur in RF00167.
