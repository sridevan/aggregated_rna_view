# UX and scientific-interpretation refinements (2026-10-18, frontend only)

Follow-up to the UX review notes (`Aggregated RNA 3D View UX notes`). Scope: Family Overview layout and selection guidance, annotation denominators and missing-data wording, the representative-comparison summary, and one explicit cross-tab action.

**Unchanged in this first round:** no backend, schema or data file was changed, and no clustering, aggregation or mapping rule was touched. The follow-up round (section 6) corrected the ligand denominator in the backend (issue A) and simplified the Structures header.

## Already present before this work (kept as is)

- **Position numbering:** shared by the Family Overview and Annotations tabs (`NumberingProvider` above the tabs).
- **Selection:**
  - selection of a nucleotide, whole loop or insertion, with toggling, Clear selection, camera focus and reset;
  - horizontal zoom and scroll.
- **Long-RNA 2D view:**
  - Clear (or a click on the selection again) returns to the fitted whole topology;
  - covered by `rna2d-focus.spec.ts`;
  - the 2D viewer internals were not modified.
- **Annotations:**
  - concise tooltips;
  - an "Evidence details" section in the pinned panel;
  - the per-category "Not counted" evidence lines.
- **Structures:**
  - representative comparison (overlay, dendrogram path, group colours);
  - the poor-fit flag;
  - the 3D alignment ⓘ, which already says the overlay differs from the metric's superposition.

## 1. Family Overview: compact family summary

`FamilyOverview.tsx` (`FamilySummary`).

**Layout:** one card replaces the About card and the six statistic cards.
- **Description:** the Rfam `comment` is clamped to three lines. A Show more / Show less button (`aria-expanded`) appears only when the text overflows. The full text and the Rfam source link are unchanged.
- **Primary counts (one row):** PDB entries, RNA chain instances ("· N validated"), Rfam consensus positions.
- **Secondary row:** a lighter row with unique sequences, mean sequence occupancy and mean structural occupancy. One ⓘ holds the definitions, with the same values and wording as before.

**Effect:** at 1400 px the molecular view now starts within 260 px of the summary card's top, compared with about 500 px before. The page width is unchanged.

## 2. Family Overview: selection guidance

- **Hint:** while nothing is selected, the toolbar shows "Click a nucleotide or loop to highlight it in 2D and 3D." (`selection-hint`). The selection summary replaces it.
- **Unchanged:** there is still no text input, and all selection behaviour is as before.

## 3. Annotations: audit of units, denominators and missing vs absent

Each pinned details panel now shows structured rows (`annotation-summary-*`): the measured quantity, the result with numerator, denominator and percentage, what the absence of a mark means, what was not counted and why, and the counting unit. The longer provenance stays under Evidence details. Tooltips are unchanged.

### Per-category definitions

| Category | Unit | Numerator / denominator (unchanged fields) | Empty mark means | Not counted |
|---|---|---|---|---|
| RNA loops (HL, IL, JL) | RNA chain instance | `n_positive_observable_and_evaluated` / `n_observable_and_evaluated` | Not annotated by BGSU (not a confirmed absence) | No observed nucleotide (`n_instances − n_observable_instances`); source unavailable (`n_observable_instances − n_observable_and_evaluated`) |
| Ligand binding sites | PDB entry | Entries with this CCD in contact at the position / entries with bound CCD (`n_entries_in_default_track`) | No PDBe contact recorded (not a confirmed absence) | See issue A |
| Solvent accessibility | RNA chain instance; Å² | Median, IQR, `n_contributing_instances` (no frequency) | n/a (every observed position has a value or is grey) | No observed nucleotide; incomplete; unsupported |
| Predicted pockets | RNA chain instance | `n_supporting_instances` / `n_evaluable_instances` | The prediction ran and no pocket was found (a prediction, not an absence of binding) | `unevaluable_reasons` (e.g. no coordinates) |
| Nonstandard nucleotides | RNA chain instance | Instances with the CCD / instances with a determined identity (`denominator`) | Standard nucleotide in the deposited sequence: a determined identity, so a real negative | `n_no_residue` (deleted / outside the aligned region), `n_missing_identity`, `n_source_unavailable` |

### Wording corrections (frontend)

- **Legend:** "No annotation" is now "Not observed". Its tooltip says that, for loops, ligands and pockets, this is not a confirmed absence. "No data" now names its cause: no observed nucleotide, or source unavailable.
- **Missing categories:** a category without data for the family is named ("Not available for RF00005: Predicted binding pockets.") instead of disappearing silently.
- **No nonstandard nucleotides:** when modification data exist but show none on the axis, the page now says so.

### Issue A: ligand denominator not adjusted for unobserved entries (resolved in section 6)

**Problem:**
- The per-ligand denominator is every PDB entry with that bound ligand.
- An entry whose RNA has no observed nucleotide at the position is still counted, so it acts as a negative.
- `tracks.json` has observability only for all ligands together (`positions[].observable.n_entries`), not per ligand. The frontend cannot correct this without loading every chain mapping.

**Size of the effect:**
- Measured with the mapping files (`experiments/ux_refinements/ligand_denominator_audit.py`, run from `data/families`).
- RF00167: 3 of 227 ligand marks are affected, all for GNG. At positions 25 and 81 it shows 2 of 3 entries where the adjusted value is 2 of 2; at position 82, 1 of 3 where it is 1 of 2.
- RF00162, RF00059, RF00005 and RF00177: 0 marks affected.

**Frontend now:** when any family entry lacks an observed nucleotide at the position, the panel says that the frequency may be underestimated. It shows the count of such entries and invents no per-ligand number.

**Smallest backend change proposed:**
- In `annotations.py`, for each default-track CCD and each consensus position, record the number of entries with that bound ligand that have an observed nucleotide at the position. Either field would do: `ligand_binding_sites.by_ccd_n_entries_observable` in `tracks.json`, or `consensus_positions_n_entries_observable` in `ligands.json` `by_ccd`.
- Use that number as the denominator, and list the unobserved entries as "not counted".
- This needs a schema field and a test. The default track's existing counts would be unchanged except the 3 GNG marks.

**Other checks:**
- No category treats unavailable source data as a negative.
- Pockets already exclude failed predictions and missing residues.
- Nonstandard nucleotides already exclude chains with no residue or an undetermined identity.

## 4. Structures: representative comparison summary (simplified in section 6)

In comparison mode, a two-line summary sits under the 3D caption (`compare-summary`). It uses stored data only.
- **Line 1:** "6UC8 B (Group 1) vs 5E54 B (Group 2) · Geometric discrepancy 0.583 · 57 shared Rfam positions".
  - The value is read from `structure_clustering.json` (`values_upper`) for the active metric.
  - TM-score is shown with 1 − TM-score.
  - The shared positions are `dataset.n_shared_positions`, the fixed correspondence used for every pair.
- **Line 2 (muted):** "3D overlay: each chain superposed on 4FE5 B using C4′ atoms (6UC8 B good, 0.27 Å · 5E54 B poor, 4.63 Å); this is not the superposition behind the geometric discrepancy". It comes from `orientation.json`.

The summary is absent in single-structure mode. Groups, representatives, colours, overlays and the camera are unchanged.

Code: `comparisonSummary` in `StructuresTab.tsx`, with a vitest file.

## 5. Cross-tab navigation (explicit only)

### Annotations → Family Overview

- **Action:** the pinned panel has a "View in Family Overview" button. It selects that Rfam position (`selectPosition`), keeps the current family, RNA instance and position numbering, and opens the tab.
- **Missing residues:** the existing selection status reports deleted, uncovered and no-coordinate positions, e.g. "Rfam position 64 is deleted in 1Y26 X: nothing highlighted".
- **Origin note:** a one-line note reads "Opened from Annotations (HPA at Rfam position 55). Family-level evidence: not necessarily observed in 1Y26 X." It shows only while that transferred position is the selection.
- **No continuous synchronisation:** a click in Annotations still changes nothing elsewhere.
- **Implementation:** a small context (`state/navigation.tsx`) provided by the page shell, which owns the tab. There is no event bus.
- **Stale state:** families are separate page loads (`/{acc}`), so selections cannot carry over between families.

### Structures → future Interactions tab (preparation only; nothing implemented)

**Identifiers already available:**
- `useStructures()`: `metric`, `groupByMetric[metric]` (the selected group) and `compareGroup`.
- `structure_clustering.json` `metrics[metric].clustering.groups[]`: `group_id`, `representative`, `members`.
- `useSelection().instance`: the app-wide selected chain.
- All chains are identified by `seq_name` (`pdb_entity_auth_label`), which carries `label_asym_id`.

**Simplest hand-over when the tab exists:**
1. Add `openInteractions(target)` to `NavigationApi`, with `target = { metric, groups: number[], representatives: string[], instance: string | null }`.
2. Implement it in the shell next to `viewInOverview`: store the target and set the tab.
3. Read the target in the new tab.

**Constraints:**
- Only identifiers are passed, never derived data.
- Group IDs are only meaningful with their metric, so the metric must always travel with them.

## Tests

- **New `tests/ux-refinements.spec.ts` (6 tests):**
  - the clamped description, its toggle and the full text, with all six statistics and their values (RF00167), and no toggle on RF00177;
  - the selection hint, plus nucleotide and loop selection, toggling and Clear;
  - the pinned panel rows for each category against the JSON (loops, ligands with the GNG caution, SASA units and IQR, pockets, nonstandard nucleotides); tooltips of 3 lines; evidence details kept; colour values unchanged;
  - View in Family Overview: deletion (1Y26 X, position 64) and a mapped position (55), with instance, numbering and URL kept, the origin note dropped on a new selection, and no pinned panel or selection change caused by the Annotations tab;
  - RF00177 transfer with 2D focus, and Clear back to the whole 2D view;
  - comparison summary values from the stored matrix for discrepancy 1–2 and TM-score 1–3 and 1–4, absent in single mode.
- **New `src/components/StructuresTab.test.ts`:** vitest for `comparisonSummary`.
- **Updated tests (also `family-rf00177.spec.ts`: pockets named as not available):**
  - `overview-layout.spec.ts`: compact summary instead of six cards;
  - `annotations-ui.spec.ts`: legend "Not observed";
  - `family-rf00005.spec.ts`: pockets named as not available.

Screenshots: the final UI is shown by the `m19_*` screenshots listed in section 6.

## 6. Follow-up: ligand denominator correction and compact comparison header

### 6.1 Ligand denominator (issue A, backend)

**Field:** `tracks.json` `positions[].ligand_binding_sites.by_ccd_entries`, keyed by the bound CCD and present for each default-track CCD with a recorded contact at the position (the same keys as `by_ccd_n_entry_entities`).
- `n_positive`: PDB entries with a recorded contact here (the existing numerator, unchanged).
- `n_observable_and_evaluated`: the new denominator. These are PDB entries with the bound CCD in which a chain of a ligand-source-evaluated entity has an observed nucleotide at the position.
- `n_positive_observable_and_evaluated`: the numerator restricted to those entries.
- `n_not_observable`: entries with no observed nucleotide here in any chain.
- `n_source_unavailable`: entries observed only in chains whose ligand annotation is unavailable.
- `frequency`: the observed-contact frequency among coordinate-observable entries, or null when no entry is eligible. It is not a binding probability.

**Counting rules:**
- The unit is the PDB entry, and the chains of one entry count once.
- The three buckets add up to `ligands.json` `by_ccd[ccd].n_entries_in_default_track`.
- Contacts, ligand identity and provenance are unchanged.

**Why `tracks.json`:** it already holds the per-position denominators, and one object avoids parallel maps. `ligands.json` is unchanged. The schema was extended (`tracks.schema.json`) and the field is required.

**Regeneration:**
- All five families' `tracks.json` were regenerated offline from the cached sources.
- The unmodified pipeline first reproduced the existing `tracks.json`, `ligands.json`, `loops.json` and reference files exactly, apart from timestamps. In RF00005 `ligands.json` the only difference was the error text of uncached 404 sources.
- After the change, every value other than the new field is identical. Only `tracks.json` was replaced.

**Results:**
- The contacts in the pipeline output and the audit script (`experiments/ux_refinements/ligand_denominator_audit.py`) agree.

| Family | Ligand marks | Marks changed | Contacts without an observed nucleotide |
|---|---|---|---|
| RF00167 | 227 | 3 | 0 |
| RF00162 | 32 | 0 | 0 |
| RF00059 | 108 | 0 | 0 |
| RF00005 | 16 | 0 | 0 |
| RF00177 | 276 | 0 | 0 |

- RF00167 GNG: the GNG entries are 3DS7 (two chains, counted once), 6P2H and 9LKC. 6P2H A has no observed nucleotide at 25, 81 and 82.

| Rfam position | Previous | Corrected | Not counted |
|---|---|---|---|
| 25 | 2/3 | 2/2 | 1 (6P2H) |
| 81 | 2/3 | 2/2 | 1 (6P2H) |
| 82 | 1/3 | 1/2 | 1 (6P2H) |

**Frontend:**
- Mark value: `by_ccd_entries[ccd].frequency`. Zero eligible entries gives the grey "No data" mark.
- Tooltip: "Recorded contact in k of n observable PDB entries (x%)".
- Pinned rows:
  - Recorded contact: "k of n observable PDB entries with bound CCD".
  - No contact: "… with an observed nucleotide: no contact recorded by PDBe (not a confirmed absence)".
  - Not counted: e.g. "1 PDB entry: no observed nucleotide at this position".
  - Counting unit.
- Evidence details: the supporting entries and the counting rule, with "observed-contact frequency, not a binding probability". The earlier "Not adjusted … may be underestimated" caution was removed.

**Tests:**
- `tests/test_annotations.py`, two new tests:
  - an observed entry with and without a contact;
  - an entry without an observed nucleotide;
  - a two-chain entry counted once;
  - an entry with unavailable ligand data excluded from the denominator;
  - zero eligible entries, giving a null frequency;
  - the buckets partitioning the entries, with the numerator never above the denominator.
- `ux-refinements.spec.ts` recomputes the three GNG values from the mapping files' coordinates in the browser and checks the panel, the mark value and the tooltip.

### 6.2 Structures: compact comparison header

**In comparison mode:**
- The subtitle "Representatives of Structural group 1 and Structural group 2" and the chain/group line are removed. The comparison controls and the structure chips already name them.
- Under "3D alignment view ⓘ" there is one line, e.g. "Geometric discrepancy 0.583 · 57 shared Rfam positions", or "TM-score X (1 − TM-score Y) · 57 shared Rfam positions".
- Values are read from the stored matrix.

**Outside comparison mode:** the existing caption is kept, and there is no comparison line.

**The ⓘ (`alignmentInfo`, values from `orientation.json`) now lists:**
- the reference chain;
- the C4′ superposition;
- each displayed chain's fit quality and C4′ RMSD;
- the sentence that this visual superposition may differ from the alignment used for the geometric discrepancy or TM-score.

**Kept:** the poor-fit flags on the chips, and the viewer, groups, colours, controls and camera.

**Size:** the heading and the metric line take 59 px above the comparison controls at 1400 px. Before, the subtitle and a two-line summary sat there.

**Tests:**
- `ux-refinements.spec.ts` covers:
  - single mode;
  - groups 1–2 in discrepancy mode, with a good fit (6UC8 B) and a poor fit (5E54 B);
  - TM-score 1–3 and 1–4;
  - the ⓘ content;
  - a 760 px window.
- The vitest file covers `comparisonSummary` and `alignmentInfo`.

Screenshots: `docs/screenshots/m19_annotations_ligand_gng.png`, `m19_structures_compare_header.png`, `m19_structures_alignment_info.png`, `m19_structures_compare_header_narrow.png`.
