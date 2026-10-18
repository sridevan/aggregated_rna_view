# Overlapping RNA loop counts at Rfam consensus position 53 (RF00167)

Investigation only (2026-10-16); no backend or frontend change. Script: `audit.py` (reads `data/families/RF00167/{tracks,loops,instances}.json`, the mapping files and the raw BGSU loop CSVs in `annotation_assessment/cache/bgsu/`); full numbers in `audit_pos53.json`.

## 1. Where 66 and 76 come from
The tooltip prints `tracks.json → positions[52].hairpin_loops.n_instances` (66) and `positions[52].junction_loops_union.n_instances` (76). Both are **counts of distinct RNA chain instances** (`seq_name`) that have at least one nucleotide of a loop of that class mapped to consensus 53. They are not loop observations (HL has 66 observations for 66 instances; J3 has 93 observations for 76 instances) and not track totals (the HL track's maximum elsewhere is 80 at L2/L3 positions). The chain of computation: BGSU loop membership (unit ids) → `resolve_unit` (author chain + author number → `label_seq_id`) → `rfam_consensus_pos` from the validated mapping → per position, set of `seq_name` per loop class → `n_instances`. Recomputing the sets directly from `loops.json` reproduces 66 and 76 exactly.

## 2. Audit at consensus 53 (deduplicated by chain instance)

| Quantity | Instances |
|---|---|
| HL | 66 |
| IL | 13 |
| J3 | 76 |
| J4 | 0 |
| JL = J3 ∪ J4 | 76 (J3 + J4 = 76; no overlap because J4 never covers 53) |
| HL ∩ JL | **62** |
| HL only | 4 (5E54 A, 5E54 B, 5SWD A, 5SWD B) |
| JL only | 14 (6P2H, 9LJN, 9LKC, 9LKE, 9LKF, 9LMC A/B, 9V4U, 9V4V, 9V4X, 9V4Y, 9V4Z, 9V50, 9V51) |
| Any loop class | 80 |
| Evaluable instances (match residue with coordinates, loop source evaluated) | 80 of 84 (the 4 9NBH chains lack coordinates at 53) |
| Evaluable without any loop | 0 |

So 66 + 76 = 142 > 84 is explained entirely by 62 instances carrying both an HL and a J3 loop at the same nucleotide. The denominators stored in `tracks.json` (`n_observable_and_evaluated = 80`) are consistent, and no positive instance is non-observable.

## 3. The overlap is in the BGSU annotations themselves
Representative instance 1Y26 chain X (consensus 53 = label_seq_id 37 = author U49):
- `HL_1Y26_002`, motif group HL_56046.1 "Purine riboswitch": raw BGSU units `1Y26|1|X|G|46 … C|53` (8 nt, one continuous segment) → consensus 50–57. This is the J2/3 strand of the three-way junction; FR3D treats it as a hairpin because the strand is closed by the non-canonical G46·C53 pair. The motif group's own name ("Purine riboswitch") confirms it is the riboswitch core strand, not a classical stem-loop.
- `J3_1Y26_001`, motif group J3_44163.4: raw units `A21–U25, A45–C54, G72–U75` (19 nt, three segments) → consensus 25–29, 49–58, 78–81. Its middle strand contains the same nucleotides 46–53.
The same pair (HL_56046.1 + J3_44163.4, 8 + 19 nt, identical consensus spans) occurs in 4FE5 B, 3RKF A–D, 9NDD A (as HL_9NDD_002 / J3_9NDD_002, problematic, one deleted position) and 48 further chains. `loops.json` membership equals the raw CSV membership for every loop (0 mismatches, 0 duplicated units, 0 units mapped to two consensus positions), so nothing is duplicated or expanded during mapping or aggregation. The overlap is a property of the FR3D/Motif Atlas classification, where hairpin and junction loops are not mutually exclusive.

The 4 "HL only" chains (5E54, 5SWD) have no J3 loop at 53 because their junction is annotated as J4 only at other positions, with `J4` excluding the 50–57 strand (3B finding). The 14 "JL only" chains are the 2025 entries (9IWF/9LJN/9LKC/9V4x series, 6P2H) whose J2/3 strand is not closed by the non-canonical pair in the deposited model, so FR3D emits no hairpin there; several of them carry **two J3 loops** that both include position 53 (e.g. 9IWF A: `J3_9IWF_001` J3_44163.4 with strands 7–11 / 31–40 / 62–65 and `J3_9IWF_002` J3_61367.1 with strands 7–11 / 31–42 / 58–65, i.e. alternative junction delimitations). That is why J3 has 93 loop observations for 76 instances; the instance count is deduplicated correctly.

## 4. Mapping correctness
- BGSU unit ids are matched to the chain by PDB id + author chain id (`1Y26|1|X|…` → instance 1y26 chain X, label asym A) and to the residue by author number + insertion code; all 5,655 loop nucleotides resolved (0 unresolved).
- Author number → label_seq_id follows the mapping files (1Y26: auth 49 → 37, offset 12; 4FE5: auth 49 → 35, offset 14), consistent with consensus 50 ↔ 1Y26 G34/auth 46 used by the viewer tests.
- Consensus 53 is a U in 82 and a C in 2 of the 84 chains (all 84 have the residue in sequence; 80 have coordinates); the BGSU reference-element transfer labels these nucleotides "J2/3" (e.g. `1Y26|1|X|U|49 → J2/3`, `4FE5|1|B|U|49 → J2/3`), i.e. the intended J2/3 junction nucleotide in every instance; `tracks.json` gives the same reference element at 53.

## Conclusion and options (for review, no change made)
The counts are correct as defined: HL and JL are overlapping FR3D classes at the J2/3 strand (consensus 50–57), and the tracks count distinct chain instances per class. If the display should not suggest mutually exclusive classes, the choices are: (a) keep both tracks but say "classes may overlap" in the ⓘ text and show the per-position "any loop" count; (b) add an "any loop (HL ∪ IL ∪ JL)" track, which at 53 reads 80 of 80; (c) flag HL_56046.1 as "junction-strand hairpin" in the HL track; or (d) resolve overlaps by precedence (junction over hairpin) in a derived view while keeping the raw classes in the data. None of these changes the underlying statistics.
