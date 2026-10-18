# Milestone 3F — Nonstandard nucleotides (final Annotations category)

Date: 2026-10-16. Family: RF00167 (Rfam 15.1, CM 1292453307, CLEN 102; 66 entries, 84 RNA instances).
Basis: `experiments/modifications_RF00167/REPORT.md`, whose methodology is implemented unchanged.

## Backend

Command:

    python3 -m rna_explorer modifications RF00167 [--pdbe cache|off] [--nakb cache|off] [--seed-dir experiments/modifications_RF00167/raw] [--offline]

Module: `rna_explorer/modifications.py`. Schemas: `rna_explorer/schemas/modifications.schema.json` and `modification_tracks.schema.json`. The RF00167 build runs in about 3 s.

**Identity**

- The primary source is `_pdbx_poly_seq_scheme.mon_id` at `label_asym_id` + `label_seq_id` in the cached PDBe updated mmCIF.
- It is checked against `_entity_poly_seq`, against `_atom_site` (model 1, where coordinates exist) and against the validated mapping's `residue_name`.
- Polymer membership comes from the sequence scheme, never from the CCD type. GTP has CCD type "non-polymer" even as a 5′ polymer residue, and non-polymer GTP ligands in four other entries are never read.

**Status of each deposited position**

| Status | Meaning |
|---|---|
| positive | a component other than A, C, G or U |
| negative | A, C, G or U; coordinates are not required |
| missing | identity cannot be established |

A position is missing for one of these reasons: microheterogeneity, an unknown component (N, UNK or X), a conflict between sequence records, a conflict with the coordinates, a conflict with the mapping, or no row in the mapping. Components are not classified.

**Consensus positions and aggregation**

- Consensus positions come only from the validated mapping, in match state. Insert, local-end and outside-region residues keep a null position and stay in the detailed data.
- Occurrence at position p is the number of instances with a positive at p divided by the instances with a positive or negative at p. Each instance counts once.
- Missing identity, no residue, and mmCIF unavailable are counted separately and excluded from the denominator.

**Optional sources**

- **PDBe `modified_AA_or_NA`** is used for cross-validation only. Each observation records reported, not_reported or not_checked, and a missing PDBe record is never treated as a negative.
- **CCD parent** (`mon_nstd_parent_comp_id`) is taken first. When it is absent, the NAKB `standard_base` is used, labelled as NAKB. Otherwise the parent is unknown. The mapping's `nucleotide_parent`, which is "X" for 50N and 50L, is not used.
- **NAKB change locations** are retained without further classification.
- **Caching.** PDBe and NAKB responses are stored as cached envelopes holding status, retrieval time, Last-Modified, ETag and body. The `--seed-dir` option seeds them from the investigation's raw files, so no new requests were needed.
- **Missing sources.** If a source is unavailable it is recorded as such, and the core output is still written.

**Outputs**

- **`modifications.json` (150 KB):**
  - `counts`; `components`, keyed by CCD, with name, type, parents and sources, NAKB change locations and PDBe agreement.
  - `observations`: one per nonstandard polymer residue, all 13 including off-axis ones.
  - `instances`: all 84, with deposited components per `label_seq_id`, residues without coordinates, missing-identity records and the mmCIF source with retrieval date.
  - `provenance`: Rfam, CM, INFERNAL, mapping and pipeline versions; mmCIF, PDBe, NAKB and CCD sources with dates, ETag and sha256.
- **`modification_tracks.json` (45 KB):**
  - Per consensus position: n_positive, n_negative, n_missing_identity with reasons, n_no_residue, n_source_unavailable, denominator, frequency, n_supporting_entries and n_denominator_without_coordinates.
  - Also per position: the components with instance and entry counts, and the supporting observation ids.
  - A summary that includes the off-axis observations, plus components, rules, method and provenance.

## RF00167 results (identical to the investigation)

| Item | Value |
|---|---|
| Entries / instances / polymer residues | 66 / 84 / 6,601 |
| Standard / nonstandard / missing identity | 6,588 / 13 / 0 |
| Residues without coordinates | 422 |
| Instances / entries with nonstandard components | 9 / 8 |
| Distinct CCD | 7 (GTP 4, CCC 2, 23G 2, IU 2, GDP 1, 50N 1, 50L 1) |
| On the consensus axis / off-axis | 6 / 7 (all `not_submitted_flank`, 5′ or 3′ terminal) |
| PDBe agreement | 11 of 13 reported and agreeing; 50N and 50L not reported; no PDBe record on a family chain without an observation |

| Rfam position | Component | Positive / denominator | Frequency |
|---|---|---|---|
| 22 | 50N | 1 / 83 | 1.2 % |
| 32 | IU | 1 / 84 | 1.2 % |
| 35 | IU | 1 / 84 | 1.2 % |
| 84 | 50L | 1 / 83 | 1.2 % |
| 87 | CCC | 2 / 77 | 2.6 % |

The other 97 positions are 0 of N, with N from 9 to 84. All 102 positions agree with `experiments/modifications_RF00167/audit/provisional_occurrence.csv` on positives, denominators, no-residue counts and residues without coordinates.

23G in 3IVN, at label_seq_id 70 (author 74) in chains A and B, is kept in `modifications.json`. It is 3′ terminal, has no coordinates, lies outside the mapped region, and has no consensus position.

**Parents**

| Source | Components |
|---|---|
| CCD | GTP, GDP and 23G (G); CCC (C); IU (U) |
| NAKB | 50N (C), 50L (G) |

## Frontend

The Annotations tab gains category 5, "Nonstandard nucleotides", with one track, "Occurrence". The data load from the optional `modification_tracks.json`: when the file is missing, only this category is hidden, and its two downloads show as unavailable.

**Rendering**

- `ScaleDef` gains a `presence` kind: one fixed colour, `#c2185b`, wherever the track has a mark. This is a deliberate exception for sparse data, because all positives are below 3 %.
- The legend reads "Nonstandard nucleotide: ■ observed".
- Positions with a known zero occurrence are left empty. A position where no instance has a determined component would be grey; RF00167 has none.
- `data-value` still carries the exact frequency.
- Nothing in the renderer is specific to RF00167 positions.

**Tooltip and pinned details**

- The tooltip reads "Rfam residue 87 / Nonstandard nucleotide: CCC / Observed in 2 of 77 RNA instances (2.6%)". When a position has several components, it lists them as "Nonstandard nucleotides: X (n), Y (m)".
- Pinned details show the occurrence with its PDB entries, then one line per component with name, parent and its source, instances and entries.
- They also show how many denominator instances have no coordinates, any missing identities, the denominator rule, and the number of off-axis residues held in the downloads.

**Popover and downloads**

- The category popover uses the requested text. Popover text now keeps paragraph breaks (`whitespace-pre-line`), and the SASA and pockets popovers use the two-paragraph text originally requested.
- Two new download items are "Nonstandard nucleotides — detailed (JSON)" for `modifications.json` and "Nonstandard nucleotides — aggregated (JSON)" for `modification_tracks.json`, making 9 items.

**Final review fixes**

- **Tooltip clipping.** The residue hover tooltip, shared by all tracks, now flips above the pointer and is clamped horizontally. Before, a tooltip on the bottom track could extend below the viewport.
- **Truncated label.** The new track's sub-label was shortened to "84 RNA instances", because it had been cut off.

## Tests

- **Backend `tests/test_modifications.py`** uses a synthetic mmCIF in which entity 1 is shared by chains A and B and a non-polymer GTP ligand is present. It covers:
  - standard and nonstandard components, and polymer membership against a same-CCD ligand;
  - the asym/entity check, and ignoring model 2;
  - shared entities, and missing coordinates with a known identity;
  - insert state, outside-region residues, an unknown component and microheterogeneity;
  - all conflict reasons, and positive/negative/missing/no-residue/unavailable counts;
  - denominators, and deduplication per instance and across duplicate PDBe records;
  - builds without PDBe and NAKB, and with seeded ones;
  - schema-valid JSON, and reconstruction from the supporting observations.
- **Backend `tests/test_modifications_rf00167.py`** checks every count above and the 5 regression frequencies with 97 zeros. It also checks 23G off-axis, 50N and 50L included without PDBe, and parents with their sources. Finally, it reconstructs every position independently from `deposited_components` and the mapping files.
- **Frontend `tests/modification-track.spec.ts`** checks:
  - one category with one track, in fifth position, with 5 fixed-colour marks at 22, 32, 35, 84 and 87 and no others;
  - the same geometry as other tracks, and axis alignment;
  - all five tooltips, and the tooltip staying in the viewport with the flip-above case exercised;
  - pinned details with CCD and NAKB parents;
  - the popover text, no toggling, Escape, keyboard focus and collapse;
  - both downloads, and the other categories unchanged;
  - a missing dataset, which hides the category and disables both downloads.
- **Updated existing specs:** category and icon counts (5), the download count (9), and category-order and popover loops that now include the new category.
- **Vitest:** the presence scale uses one colour, and the frequency scale is unchanged.

Screenshots:

- `docs/screenshots/m3f_annotations_all_categories.png`
- `docs/screenshots/m3f_nonstandard_pinned_details.png`
- `docs/screenshots/m3f_popover_nonstandard.png`

## Limitations

- Only RF00167 has been built. Families with natural modifications, such as tRNA or rRNA, may hit microheterogeneity, UNK or conflicts. These are handled as missing identity, but the rules have only been tested on synthetic data.
- Most RF00167 observations are terminal components, and 7 of 13 are off-axis. The track therefore shows only 6 of the 13 residues; the rest are in `modifications.json`.
- Denominators include residues without coordinates, because identity comes from the deposited sequence. They can therefore differ from the loop, pocket and SASA denominators at the same position.
- The PDBe inclusion criterion and the meaning of NAKB `count` are undocumented. Neither affects the output.


## Revision: one row per component (user request, 2026-10-16)

The single "Occurrence" track was replaced by one row per CCD component, mirroring the per-ligand rows of "Ligand binding sites". The change is frontend-only, and `modification_tracks.json` is unchanged.

- **Which rows.** A component gets a row only when it occurs on the consensus axis, as ligand rows do for ligands with mapped contacts. RF00167 has four rows, in order CCC (2 instances), 50L, 50N and IU, sorted by instances and then id.
- **Off-axis components.** GTP, GDP and 23G occur only outside the mapped region. They have no row and stay in `modifications.json`. A row's pinned details say how many further residues of that component lie off-axis.
- **Row value.** Each row shows instances with that component at p, out of instances with a determined component at p. The denominator is the same as before. Marks keep the fixed presence colour.
- **Row label and tooltip.** The label shows the CCD id, with the component name and instance count beneath. The tooltip reads "Rfam residue 87 / CCC — Cytidine-5'-phosphate-2',3'-cyclic phosphate / Observed in 2 of 77 RNA instances (2.6%)".
- **Pinned details.** These show the name, parent with its source, instances and entries, and residues without coordinates. When several components share a position, a combined line is added.
- **Empty category.** The category is hidden when a family has no on-axis components.
- **Tests and screenshots.** `tests/modification-track.spec.ts` now checks four rows, their order, absent rows for GTP, GDP and 23G, marks per row, tooltips, the flip-above tooltip on the bottom row, and pinned details. The screenshots were regenerated.

## Revision: category order (user request, 2026-10-16)

The categories and their legend entries now appear in this order:

1. Solvent accessibility
2. RNA loop annotations
3. Nonstandard nucleotides
4. Ligand binding sites
5. Predicted binding pockets

Only the order changed; the contents of each category are as before. The tests that check the order were updated. All 39 Playwright tests pass, and the screenshots were regenerated.

## Final Annotations UI simplification (frontend only, 2026-10-16)

**Header**

- The title is now "Annotations", and the redundant eyebrow chip was removed.
- The text below it is generated from the family counts: "Residue-level annotations across 84 RNA instances from 66 PDB entries."
- The long instructions were removed, and the Download dropdown stays at the top right.
- At 1400×1000 the consensus axis now starts about 90 px higher.

**Colour scales**

- RNA loops, ligand binding sites and predicted pockets share one teal scale, 0–100% in five bins. A given frequency gets the same colour in all three categories.
- The ligand blue and pocket purple scales were removed. Numerators, denominators and tooltips are unchanged, and each category's popover still describes its own source.
- SASA keeps its orange quantity scale, with the domain taken from `sasa_tracks.json`.
- Nonstandard nucleotides keep solid magenta. It is labelled "Nonstandard nucleotide observed" and does not encode frequency; exact frequencies stay in tooltips and pinned details.

**Legend**

- One compact legend replaces the five separate ones, on one row on desktop and wrapping on narrow screens. It shows occurrence frequency, median SASA, the magenta swatch, "No data" and "No annotation".
- Hovering "No data" explains: status cannot be determined or the data are unavailable; for SASA, no qualifying measurement.
- Hovering "No annotation" explains: assessed, and the annotation is not present.
- Swatches have a faint ring so the lightest teal stays visible on white.

**Unchanged**

Calculations, JSON files, tracks, track and category order, the axis, mark geometry, popovers, tooltips, pinned details and downloads are all unchanged.

**Tests**

- New: `tests/annotations-ui.spec.ts` checks the header and dynamic counts, the Download position and its 9 items, and the reduced space above the axis.
- It also checks that every loop, ligand and pocket mark is filled with the shared bin colour for its value, and that equal bins match across categories. SASA must stay orange and nonstandard marks magenta. Backend values and tooltips must be unchanged.
- It checks the legend contents and order, category order and track counts, axis alignment, and the narrow screen at 390 px.
- A final test confirms "no data" stays distinct from "no annotation", using a pocket file altered only in the browser for the test.
- Updated: the legend IDs in the pocket and modification specs.
- New screenshot: `docs/screenshots/m3g_annotations_final.png`. The earlier screenshots were regenerated.
