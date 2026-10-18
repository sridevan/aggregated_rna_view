# Family Overview Phase 2: data audit for per-instance loop tracks

`audit_loops.py` is a read-only audit of `data/families/RF00167/loops.json` against the validated mappings and the per-instance annotation files. Its output is `audit_loops.json`.

## Source

`loops.json` holds the BGSU RNA 3D Motif Atlas loops for every RF00167 entry: 496 loops (234 HL, 82 IL, 114 J3, 66 J4) with 5,655 nucleotides. Each loop carries:

- the original `loop_id` and `loop_type`;
- `strands`, the ordered segments, with start and end unit ids and `nucleotide_indices` into the member list;
- the full `nucleotides` list, each with `seq_name`, `label_asym_id`, `auth_asym_id`, `label_seq_id`, author number and insertion code, BGSU `unit_id`, `alignment_state`, `rfam_consensus_pos`, `has_coordinates` and `mapping_status`.

## Results (no issues found)

**Membership is fully recoverable.** Every member resolves to a mapped instance and agrees with the mapping on `label_seq_id`, author number and insertion code, residue name, alignment state, consensus position and coordinate availability.

**Strands partition each loop's nucleotide list,** and `label_seq_id` has no gaps within a strand.

**The per-instance files agree with `loops.json`.** The annotation records' per-residue loop lists match the membership for every instance.

**Off-grid members are kept as members:**

| Kind | Count |
|---|---|
| Outside the mapped region (flanks) | 168 |
| Insertions | 68 |
| Loops with no consensus position at all (e.g. 9NBH A) | 20 |

**Every loop nucleotide has coordinates.** The missing-coordinate path is therefore covered by unit tests only.

**Four loops are inter-chain internal loops** in the dimeric constructs 9NBH and 9NDD (A–C and B–D, 4 + 4 nucleotides).

**Overlaps within one track occur only in the junction track:** J3 and J4 loops share nucleotides, giving 1,079 position-level overlaps. HL and IL have none.

## Residue correspondence used by the frontend

| Purpose | Identifier |
|---|---|
| BGSU loop membership | `unit_id` (PDB, model, author chain, residue, author number, insertion code) |
| Instance identity | `seq_name` / `label_asym_id` + `label_seq_id` |
| 2D viewer | `label_seq_id` (PDBe RNA Viewer, author chain) |
| Mol* | `label_asym_id` + `label_seq_id` |
| Track columns | `rfam_consensus_pos`, match state only |

Residue numbers alone are never used.
