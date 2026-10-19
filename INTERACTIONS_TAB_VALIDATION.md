# Interactions tab: implementation and validation (RF00167)

> **Revision 2 (2026-10-19)** adds base–ribose (Br) and base–phosphate (Bph), an interaction-type dropdown, a continuous conservation
> scale, "Multiple interactions" grouping, new default sorting, a 300-position cutoff for circular diagrams and a default comparison pair.
> It also shortens the tooltips. See **§11**, which supersedes the earlier sections where they differ: filters (§1), the five colour bins
> (§4), table columns (§1) and the test counts (§9).
>
> **Revision 3 (2026-10-19)** adds per-category colour families and line patterns, a two-dimension legend, Select all / Reset to Bp and a
> unified table search that replaces the toolbar lookup. The diagram is also larger. See **§12**, which supersedes §11.5 (the single blue
> scale) and the lookup parts of §5 and §11.8.
>
> **Revision 4 (2026-10-19)** adds table pagination, the split data format (schema 1.2) and the **16S (RF00177)** dataset. See **§13**.

Date: 2026-10-19. Scope: the pasted specification "Implement the Interactions Tab — RNA 3D Family Explorer". The scientific basis is the
RF00167 investigation in `experiments/rf00167_interactions/` (`RF00167_interaction_variation_investigation.md`).

## 1. What was implemented

There is a new **Interactions** tab (fourth tab). It shows RNA–RNA base pairs and base stacking, as annotated by BGSU FR3D, on Rfam
consensus positions. The tab has two states, chosen with a compact selector: *Family conservation* | *Compare two*.

**Family conservation (default)**
- One circular diagram of every consensus-mapped interaction observed anywhere in the family. Each arc is coloured and sized by the
  fraction of eligible chains in which it is observed.
- A details card for the selected interaction.
- A searchable, sortable table with these columns: Interaction, Type, FR3D class, Conservation (k/n, %), PDB entries, Region, and
  **In displayed structure**.
- The PDBe Mol* viewer, showing one *displayed structure*. The selected interaction's nucleotides are highlighted and focused.
- When the interaction is absent from, or not comparable in, the displayed structure, the details card offers an **Inspect in another
  structure** selector. It lists the actual supporting chains.

**Compare two**
- Structure A and Structure B selectors, a *Back to family conservation* button, and an *All interactions / Emphasise differences*
  switch.
- Two circular diagrams on the identical layout. Arcs are coloured by comparison status: grey = shared, vermilion = only in this
  structure, reddish purple = reclassified, dotted = not comparable in the other structure.
- A comparison table with status filter chips (Shared, Only in A, Only in B, Reclassified, Not comparable). It shows the FR3D class in A
  and in B, and the region.
- A two-structure 3D view with an interaction-local superposition.

**Shared controls**
- Numbering reference: default reference, any other chain, or Rfam consensus.
- Residue and interaction lookup.
- Base pairs (default on) and Stacking checkboxes.
- A *More filters* popover for base-pair type (canonical / noncanonical), context (secondary / tertiary) and structural region.

## 2. Files

| File | Change |
|---|---|
| `rna_explorer/interactions.py` | **new** – `interactions` pipeline step (production version of the investigation's `map_interactions.py` + `analyse.py`) |
| `rna_explorer/schemas/interactions.schema.json` | **new** – schema; `build_interactions` validates its output |
| `rna_explorer/cli.py` | `interactions` sub-command |
| `data/families/RF00167/interactions.json` | **new** generated dataset (964 KB uncompressed; loaded only when the tab opens) |
| `experiments/rf00167_interactions/raw/` | **new** – the 132 BGSU FR3D CSVs (gzip, 540 KB) + `PROVENANCE.json` (URL, retrieval date, sha256), copied read-only from the research worktree so the build is reproducible offline |
| `tests/test_interactions_rf00167.py` | **new** – backend tests |
| `frontend/src/types/backend.ts`, `frontend/src/api/data.ts` | `InteractionsDataset` types, `loadInteractions` (memoised, optional) |
| `frontend/src/lib/interactions.ts` (+ `.test.ts`) | **new** – per-instance status, comparison, filters, numbering labels, lookup, circular geometry, base atom correspondence, Kabsch fit |
| `frontend/src/components/interactions/CircularDiagram.tsx` | **new** – reusable SVG circular diagram |
| `frontend/src/components/interactions/InteractionPairViewer.tsx` | **new** – two-chain 3D view with local superposition (uses the Structures tab's 3Dmol helpers) |
| `frontend/src/components/InteractionsTab.tsx` | **new** – the tab |
| `frontend/src/components/Tabs.tsx`, `frontend/src/App.tsx` | fourth tab (the tab bar scrolls horizontally on narrow screens) |
| `frontend/src/components/structures/AlignmentViewer.tsx` | helpers only: `Atom.seq` (label_seq_id) is now parsed, and `applyRT` is exported. The viewer's behaviour is unchanged. |
| `frontend/tests/interactions.spec.ts` | **new** – 8 Playwright tests |

## 3. Data processing

```
python3 -m rna_explorer interactions RF00167 --seed-dir experiments/rf00167_interactions/raw [--offline]
```

**Inputs**
- BGSU RNA 3D Hub FR3D CSVs per entry: `https://rna.bgsu.edu/rna3dhub/pdb/{PDB}/interactions/fr3d/{basepairs|stacking}/csv`. They
  are fetched through the cache; a seed directory provides them offline.
- The validated per-chain mappings.
- The coordinate files the pipeline already uses: PDBe updated mmCIF, or ModelServer chain files for selection builds.

**Method.** This is the investigation's method, unchanged:
- **Unit-id join:** (PDB, author chain) → instance; (author number, insertion code) → mapping residue → `label_seq_id` → Rfam state.
- **Endpoint check:** every mapped endpoint is checked against the mmCIF (author chain, number, insertion code and residue name at that
  `label_seq_id`, model 1).
- **Normalisation and de-duplication:** interactions are normalised 5′→3′ by `label_seq_id`, with the FR3D class flipped where needed
  (cWH↔cHW, s35↔s53). Duplicates over the two row orders and over alternate locations are removed. FR3D classes are never merged.
- **Eligibility:**
  - the chain is FR3D-annotated;
  - both positions are match residues with coordinates;
  - both bases are complete (FR3D standard base heavy atoms; for modified residues, the six pyrimidine-ring atoms).
  - Per instance and position this is stored as a code string: `E` eligible, `I` incomplete base, `M` not modelled, `-` no match
    residue.
- **Near interactions:** these ('n' classes) are counted as `near_only`, never as observed.

**`interactions.json` keeps**
- Per instance: every intra-chain interaction, including insertion and flank endpoints. Each record has the original BGSU unit ids,
  `label_seq_id`s, the Rfam position or `anchor^offset`, the FR3D class, the mmCIF check result and the altlocs.
- Per instance: the eligibility and nucleotide strings, and annotation availability with a reason.
- Family-wide pairs: id `category:i-j`; observed / eligible chains; observed / eligible entries; the class and nucleotide-pair
  distributions; canonical count; and every non-eligibility count (no residue, not modelled, incomplete, not annotated, near only,
  observed-but-ineligible).
- Regions from `reference_structural_elements.json`; `structure_class` (secondary = Rfam SS_cons pair; tertiary = pair between two
  different validated elements; otherwise unclassified).
- The default numbering reference.
- The between-chain interactions, which are excluded from aggregation.
- Per-file sources (URL, retrieval date, sha256) and pipeline provenance.

**Family independence.** Nothing RF00167-specific is in the code; residue numbers, frequencies and PDB ids all come from the data. As a
dry run (not written to `data/`), the step was run on RF00059 and RF00162:

| Family | Interactions | Consensus pairs | Chains annotated | Endpoints verified against mmCIF | Default numbering reference |
|---|---|---|---|---|---|
| RF00059 | 5,018 | 166 | 38/38 | 20,176/20,176 | 2GDI X |
| RF00162 | 4,424 | 199 | 29/29 | 17,750/17,750 | 5FJC A |

Those families' datasets are **not** generated in `data/` yet, so their tabs show "No RNA–RNA interaction data have been prepared".

## 4. Circular diagrams

**Layout**
- Consensus positions run clockwise from just right of 12 o'clock, with a fixed 10° gap before the 5′ end.
- The span is the family's first to last position that any chain covers (RF00167: 1–102). Angles depend only on (span, size), so every
  diagram of a family shares them. This was checked in Playwright: an arc shared by two diagrams has an identical path, and position
  markers have identical coordinates.

**Arcs**
- Arcs are circular arcs orthogonal to the circle, the R3DCID convention (Fonga, Abubaker, Thaman & Zirbel, *RNA* 2026,
  doi:10.1261/rna.081031.126).
- They are drawn back to front by importance: low conservation first, base pairs over stacking.
- Stacking is dashed.

**Family encoding:** one blue ramp in 5 bins (< 25, 25–50, 50–75, 75–95, ≥ 95 %). Line width rises with the bin, so the encoding never
relies on colour alone. A legend explains it.

**Labels and markers**
- The selected numbering's round numbers (every 10).
- A grey region band with element names (P1, J1/2, P2, L2, J2/3, P3, L3, J3/1).
- The selected or looked-up nucleotides, with their full labels and an amber marker.

**Interaction**
- **Selection:** the selected arc is redrawn on top with a dark halo, and the other arcs are not dimmed.
- **Hover and click:** both are handled once on the SVG through invisible 8 px hit paths. Hovering shows residues in the selected
  numbering, Rfam positions, type, the FR3D class distribution, observed/eligible (%), entries and the displayed-structure status.

**Compare mode:** each diagram draws only its own observed interactions. Positions that are unmodelled, incomplete or unmapped in that
chain are drawn as hollow markers.

**Performance**
- Arcs, rings and hit paths are memoised, so selection, hover and numbering changes redraw only the overlay; numbering redraws only the
  labels.
- Precomputed data and display state are kept separate.
- RF00167 has 199 pairs, 49 of them base pairs by default.
- For a 16S-sized family (about 1,500 positions and several thousand pairs) the plain SVG remains workable, but labels would need
  thinning (already rule-based) and zooming. That has not been tested here.

## 5. Numbering and lookup

**Default reference: 4FE5 chain B.** This is the family display representative (`family.json display_representative`): *B. subtilis*
*xpt* guanine riboswitch at 1.32 Å, the highest resolution in the family. Its numbering is the one used in the literature (G37–C61,
U22–A52, C74).
- It has 67 eligible match positions.
- It has no residue at Rfam 18, 64, 72 and 88. Those are the terminal P1 pair and two L3 columns.
- These positions fall back to "Rfam 18" and so on, and are never given an invented author number.
- 153 of the 199 family pairs can be evaluated in it.

**Alternatives**
- 4TZX X (2.01 Å) covers 69 positions.
- The 9NDD scaffold covers 97, but it is a synthetic construct.
- Neither improves on 4FE5 B for user-facing numbering, so the rule stays: use the family display representative. It is recorded in
  `default_numbering_reference`.

**Behaviour**
- The numbering selector offers the default reference, any other chain, or Rfam consensus.
- It is separate from the **displayed structure** selector, which drives the 3D view and the status column.
- Changing the numbering changes labels, tooltips and lookup only. Ids, frequencies, arc geometry and the selection stay the same
  (Playwright test).

**Lookup** accepts `G37`, `37`, `G37-C61`, `C61 G37`, `Rfam 41` and `41-66`:
- **Number resolution:** numbers are read in the selected numbering, as author number plus insertion code, or as Rfam positions under
  Rfam consensus or with the `Rfam` prefix.
- **Residue name check:** a residue name must match. "A37" gives "37 of 4FE5 B is G37, not A37".
- **Insertions and ambiguity:** a residue that is an insertion or a flank of the reference is reported as such, and an author number
  matching several positions is reported as ambiguous. Neither is guessed.
- **A pair:** selects the interaction (base pair preferred). If no interaction exists, both positions are marked and a message says so.
- **A single residue:** marks it and limits the table to its interactions.

## 6. Present / Absent / Not comparable

`statusIn` (frontend) applies the backend's eligibility definition to the displayed chain:
- **Not comparable:** the chain is not annotated, or either position is not `E`. The reason is given per position, e.g. "Rfam 18: no
  residue at this Rfam position".
- **Present:** a non-near observation of that category exists. Its classes are shown.
- **Absent:** otherwise. Any other annotation between the same two nucleotides (another category, or a near class) is listed in the
  details.

The column is headed **In displayed structure** and follows the displayed chain. *Inspect in another structure* lists exactly the chains
where the interaction is present, with the numbering reference first and then by resolution. Choosing one:
- changes only the displayed chain;
- keeps the numbering, the frequencies and the family mode.

## 7. 3D highlighting and local superposition

**Family mode** reuses `PdbeMolstarViewer`:
- Mol* structure switching never re-initialises the canvas.
- The selection is the `label_seq_id`s of the selected positions that have coordinates in the displayed chain.
- Mol* draws them as ball-and-stick on the slate cartoon and centres on them. Its focus representation has no contact lines, so an
  absent interaction shows two nucleotides without implying a contact.
- A not-comparable selection highlights only the nucleotides that are modelled, and the caption explains the limitation.
- A selection change never reloads the structure.

**Compare mode** uses `InteractionPairViewer`, built on the Structures tab's 3D helpers: the same 3Dmol build, ModelServer chain fetch,
`orientation.json` family-frame transforms, principal-axes view and tube.
- **Without a selection:** whole-RNA superposition in the family frame, as in the Structures tab.
- **With a selection:** chain B is moved by a least-squares fit (Horn quaternion) of the **base heavy atoms of one nucleotide** (the
  interaction's first; the second if the first is not modelled in both) onto A's. Its partner keeps its own position relative to that
  base, so pair or stacking geometry differences remain visible.
- **Atom correspondence** follows nucleotide identity:
  - the same parent nucleotide uses the ring and exocyclic atoms present in both;
  - two different purines or two different pyrimidines use the ring atoms;
  - a purine against a pyrimidine uses N9/N1, C4/C2, C8/C6 and C1′;
  - modified residues use the atoms they have.
- **Fit failures:** with fewer than 3 atoms, or a missing residue, the previous valid superposition is kept and the caption says why.
- **Display:** both chains' nucleotides are drawn as sticks over faded tubes, also when the interaction is absent in one. Missing
  nucleotides are simply not drawn.
- **No refetching:** chain coordinates are cached; a new selection rebuilds the scene from the cache. Playwright checks that the
  fetch counter is unchanged.
- **Geometry check:** the C1′–C1′ distance of the selected pair is identical before and after the fit in both chains (Playwright,
  < 1e-6 Å).

The Structures tab viewer and its whole-RNA alignment are unchanged; the Structures Playwright specs still pass (see §9).

## 8. Scientific validation (RF00167)

1. **Frequencies reproduce the investigation.** All 199 consensus pairs match `RF00167_interaction_conservation.csv` exactly (pytest).
   This covers observed and eligible chains, observed and eligible entries, near-only, no residue, not modelled, incomplete base and
   the class distribution. The summary matches stage 1 exactly:
   - 9,991 interactions;
   - 41,604/41,610 endpoints verified; the 6 mismatches are 1Y27 C81, written `CCC` in the mmCIF;
   - endpoint status: 37,382 consensus, 98 insertion, 4,130 outside the alignment, 794 not in the RNA polymer.
2. **Denominators.**
   - For every pair, the eligible + missing + incomplete + unmapped + unannotated chains sum to 84.
   - Observed + near-only + not-annotated sums to the eligible count.
   - A recount from the per-instance observations equals the stored count (pytest).
3. **Kissing loop.** G37–C61 (Rfam 41–66) is 80/80 eligible chains (65/65 entries), cWW, tertiary, L2 × L3. The other L2–L3 pairs are
   also listed: G38–C60 75/80, A33–A66 68/80 and others.
4. **Junction variation.** With the region filter set to J2/3, most junction base pairs are below 90 % (Playwright). Low-frequency
   junction pairs such as U47–U48 and U47–U51 (1/75) appear alongside the conserved ones.
5. **5E54 B (apo) vs 5SWE (adenine-bound, final state)**, compared through the data alone (vitest and Playwright). Residues use *V. vulnificus* *add* numbering (5E54 B as the numbering reference):
   - **only in 5SWE:** U22–A52 (cWW), A23–G46 (tWS), U49–A76, C50–U75;
   - **only in 5E54 B:** U47–A52, U49–A73, A23–U74 (tHW), U48–U75;
   - **shared:** U20–A76 and the kissing loop.
   This is the junction rearrangement described by Stagno et al. 2017 (*Nature*, doi:10.1038/nature20599). No observation is
   hard-coded in the components.
6. **Not comparable is never absent.** In 4FE5 B, the 18–88 pair (no residue at Rfam 18) is Not comparable with its reason. Changing
   the displayed chain updates every status to the independently computed expectation (Playwright, all 49 base-pair rows).

## 9. Tests

- **pytest:** 97 passed, including the 5 new `test_interactions_rf00167.py` tests (helpers, schema and summary, reproduction of the
  investigation, denominators, known interactions).
- **vitest:** 88 passed, including 9 new `lib/interactions.test.ts` tests:
  - status and reasons; comparison categories, including reclassified; filters;
  - RF00167 frequencies independent of numbering; lookup and its errors; the Stagno comparison; supporting chains;
  - layout and arc geometry; Kabsch; base correspondence.
- **tsc:** clean.
- **Playwright:** `tests/interactions.spec.ts`, 8 tests:
  - table and diagram completeness with stored k/n and status per row;
  - arc ↔ row ↔ details ↔ 3D residue synchronisation, and numbering invariance;
  - the absent-interaction workflow and supporting-chain inspection;
  - hover content; lookup; region and tertiary filters;
  - the 5E54 B vs 5SWE comparison with the local fit, geometry preservation and no refetch;
  - the no-data notice.
- **Full Playwright suite** (167 tests): 162 passed and 5 failed on the first run.
  - Four were narrow-screen overflow checks: the fourth tab made the tab bar too wide. The tab bar now scrolls horizontally
    (`Tabs.tsx`), and all four pass.
  - One was a screenshot test whose Structures 3D viewer reported a load error. It passed on rerun, together with all screenshot and
    Interactions tests.
  - Historical screenshots regenerated by the suite were restored; two new ones were added:
    `docs/screenshots/m23_*` (replaced in revision 2, §11.14).

## 10. Limitations and suggested refinements

- **Interactions not compared.** Insertion and flank interactions are kept in the data and counted in the comparison note, but they are
  neither drawn nor compared. Between-chain interactions (332 in RF00167) are excluded.
- **No near base pairs.** BGSU's CSV export has none (near stacks only), so near base pairs cannot be shown.
- **Base triples** appear only as their component pairs.
- **Chain weighting.** Chains are not independent: frequencies count chains, with entries shown alongside, and are not weighted.
- **Low resolution.** Low-resolution chains (e.g. 9V50, 9LMC) may lack annotations for reasons of resolution.
- **Structure classes.** Secondary/tertiary is classified only for base pairs, and only from SS_cons and the validated elements. Pairs
  inside one element that are outside SS_cons are "unclassified".
- **Tooling.** The 3D superposition in compare mode uses the 3Dmol helpers of the Structures tab, as PDBe Mol* has no two-chain
  superposition in this app. Family mode uses PDBe Mol*.
- **Other families.** RF00059 and RF00162 can be built with the same command, but their datasets are not generated or reviewed yet.
- **Not tested at large scale.** Sizes like 16S (thousands of arcs) have not been tested; zoom/pan and adaptive label density would be
  the next step.
- **Possible refinements:**
  - an R3DCID-style "show / dim / hide" per category;
  - linking an arc to the BGSU interaction page;
  - a download of the filtered table;
  - near-stacking display.


## 11. Revision 2 — four interaction types, conservation scale and UX

### 11.1 Files

| File | Change |
|---|---|
| `rna_explorer/interactions.py` | base-ribose and base-phosphate added to the same step (no second processor); directional handling; category-specific eligibility; per-category summary; schema 1.1 |
| `rna_explorer/schemas/interactions.schema.json` | four categories; BR/BPh class patterns; `directional`; `sugar_eligibility`, `phosphate_eligibility`, `unavailable_categories` |
| `experiments/rf00167_interactions/raw/` | 132 further BGSU files (`{PDB}_baseribose.csv.gz`, `{PDB}_basephosphate.csv.gz`) with URL, retrieval date and sha256 in `PROVENANCE.json` |
| `data/families/RF00167/interactions.json` | rebuilt (1.3 MB uncompressed) |
| `tests/test_interactions_rf00167.py` | 7 tests (2 new) |
| `frontend/src/lib/interactions.ts` (+ test) | categories and codes; category-specific `statusIn`; `conservationColor`; `groupByEndpoints` and `groupOccurrence`; `defaultFamilyOrder`; `defaultComparePair`; `CIRCULAR_MAX_CLEN` and `showCircular`; `categoryAvailability` |
| `frontend/src/components/InteractionsTab.tsx` | the UI changes below |
| `frontend/src/components/interactions/CircularDiagram.tsx` | double-line arcs for multiple interactions |
| `frontend/src/types/backend.ts` | new fields |
| `frontend/tests/interactions.spec.ts` | 14 tests (6 new) |

`InteractionPairViewer`, the local fit, `PdbeMolstarViewer` and the Structures tab are unchanged.

### 11.2 Base–ribose and base–phosphate: retrieval and processing

**Source.** BGSU RNA 3D Hub FR3D exports, one file per entry and category:
- `https://rna.bgsu.edu/rna3dhub/pdb/{PDB}/interactions/fr3d/baseribose/csv`
- `https://rna.bgsu.edu/rna3dhub/pdb/{PDB}/interactions/fr3d/basephosphate/csv`

They have the same 3-column format as base pairs: unit id, class, unit id.

**Classes** are numbered 0–9 (`0BR`…`9BR`, `0BPh`…`9BPh`, after Zirbel et al. 2009). An `n` prefix marks a near interaction, which is
counted as near-only and never as observed.

**Directionality.** Each row is written from the base nucleotide to the nucleotide whose ribose or phosphate it contacts.
- RF00167 has 3,018 rows. None is self-referential, and no row reappears with its endpoints swapped.
- 14 rows have a partner with the same class and the roles exchanged. These are two distinct contacts (each base touches the other's
  backbone).
- So, unlike base pairs and stacking, Br and Bph are **not** re-ordered by sequence and their classes are **not** flipped.
- The pair id `category:pos_1-pos_2` has `pos_1` = base and `pos_2` = contacted backbone; `pos_1` may be greater than `pos_2`.
- De-duplication is over alternate locations only.
- Observations keep the BGSU order of unit ids.

**Join and checks** are the same as before: unit id → instance → mapping residue → `label_seq_id` → Rfam position, and every endpoint is
checked against the mmCIF.
- Br: 4,764 of 4,764 endpoints verified.
- Bph: 1,176 of 1,176 endpoints verified.
- Rows involving a ligand (e.g. HPA or ADE contacting a ribose) are counted as unmapped and are outside RNA–RNA scope.

**Comparison with the raw files.** For 4FE5 B, 5E54 B and 1Y26 X, every RNA–RNA Br/Bph row of the raw BGSU file is present in its
original order and class (pytest). Examples:
- 4FE5 `U49 9BPh U48`;
- 1Y26 `U22 9BR U47`.

**Br and Bph do not change Bp or Bs.**
- All 199 Bp/Bs pairs and every per-instance Bp/Bs observation are identical to the previous build, as are the eligibility strings and
  annotation states.
- The pytest that reproduces the investigation CSV still passes.
- Chain "annotated" status is still defined from base pairs and stacking only.

**RF00167 after the rebuild**

| Category | Interactions | Consensus pairs |
|---|---|---|
| Bp | 2,702 | 49 |
| Bs | 7,289 | 150 |
| Br | 2,284 | 41 |
| Bph | 555 | 22 |

The four categories give 209 nucleotide pairs, 48 of them with more than one annotation.

### 11.3 Category-specific eligibility

| Category | Nucleotide 1 (pos_1) | Nucleotide 2 (pos_2) |
|---|---|---|
| Bp, Bs | complete base (FR3D standard heavy atoms; six ring atoms for modified residues) | complete base |
| Br | complete base | complete ribose: C1′ C2′ C3′ C4′ O4′ O2′ |
| Bph | complete base | complete phosphate: P OP1 OP2 O5′ (or O1P/O2P) |

- **Shared conditions:** the chain must be FR3D-annotated (from base pairs and stacking), and the category's file must have been
  retrieved for the entry (`unavailable_categories`). Coordinates come from model 1.
- **Storage:** completeness is stored per position (`sugar_eligibility`, `phosphate_eligibility`; codes `E`/`I`/`M`/`-`).
- **Effect in RF00167:** 31 chain positions have an incomplete or absent phosphate, typically a 5′ end without P, and 4 have an
  incomplete ribose. A Bph or Br contact to such a nucleotide is Not comparable ("Rfam 1: incomplete or absent phosphate"), never absent.

### 11.4 Filtering

- **Interaction types:** the Base pairs / Stacking checkboxes are replaced by an **Interaction types** dropdown. It has one checkbox per
  category, Base pairing (Bp), Base stacking (Bs), Base–ribose (Br) and Base–phosphate (Bph), and the button shows the active codes.
  Only Bp is on at first.
- **Where it applies:** the family diagram and table, both comparison diagrams, the comparison table and its counts. Filtering never
  recomputes the dataset.
- **Missing categories:** a category with no source files is shown greyed out as "no data" and cannot be enabled. A category with data
  but no consensus interaction shows "none observed".
- **Base-pair-only filters:** base-pair type (canonical / noncanonical) and context (secondary / tertiary) are in *More filters*, with a
  note that they apply to base pairs only; other categories are hidden while they are set. Region applies to every category.

### 11.5 Continuous conservation scale

- **Scale:** arc colour is interpolated continuously between seven ColorBrewer *Blues* stops, from `#c6dbef` at 0% to `#08306b` at 100%.
  It encodes observed / eligible chains. Darkness strictly increases with frequency (vitest).
- **Encoding:** arc width is constant (1.7). Stacking stays dashed.
- **Selection:** the selected arc keeps its dark halo.
- **Legend:** a compact 0–100% gradient bar.
- **Exact values** appear in the tooltip, the table (k/n and %, with a colour chip whose colour is checked against the frequency in
  Playwright) and the details panel, so the colour is never the only cue.
- **Comparison diagrams** keep their status colours and do not use the gradient.

### 11.6 Multiple interactions

**Grouping.** Interactions are grouped for display by their unordered pair of nucleotides. Each member keeps its id (`category:i-j`), its
classes and its own frequency; nothing is merged in the data.

**Different classes within one category** across chains (e.g. U34–A65 tWH/tWW) are a separate phenomenon. They are shown as that
interaction's class distribution and are not a "multiple" group.

**When a group has more than one member** under the active filters:
- **Diagram:** one arc, drawn as a double line (widened stroke with a white core).
- **Table:** a group row with the type "Multiple (Bs+Br)", each member's code and class, and a collapsible list of member rows, each with
  its own conservation, entries and status.
- **Details panel:** lists every member with its own k/n and status. Members are selectable, and the selected member is highlighted.

**Group occurrence** ("any") is computed explicitly and never as a sum or average of member frequencies:
- denominator: chains in which at least one member can be evaluated;
- numerator: those chains in which at least one member is present.

It colours the grouped arc and is labelled "any" wherever it is shown.

### 11.7 Sorting, terminology and tooltips

**Family table**
- **Default order:** frequency descending, then eligible chains descending, then Rfam positions, all on numeric values. So 80/80 comes
  before 4/4, and 73/80 before 13/80. Playwright checks the exact order.
- **Column sorting:** every column header sorts, and **Default order** restores the default. Sorting keeps the selection and never
  reloads 3D.
- **Columns:** Interaction, Type (Bp/Bs/Br/Bph or Multiple), **Interaction class** (renamed from "FR3D class", with an info icon
  explaining BGSU/FR3D codes per category), Conservation (info icon on how the frequency is computed), Entries, In displayed structure.
- **Region** moved from a column to the row tooltip and the details panel (the region filter remains), so the table no longer scrolls
  horizontally.

**Comparison table:** differences first (only in A, only in B, reclassified, not comparable), then shared, then positions. Column sorting
and the status filter chips remain.

**Tooltips**
- The diagram info uses the requested wording.
- Arc tooltips show the residue pair, Rfam positions, one line per annotation (code, class, k/n, %, status in the displayed structure)
  and "any" for groups.
- Methodology moved into the info icons.

**Details panel**
- With nothing selected, it is a one-line prompt.
- With a selection, it shows the reference residues, Rfam positions (an arrow for directional types), type, class distribution,
  conservation, PDB entries, region, chains not evaluated, and the Present / Absent / Not comparable box with **Inspect in another
  structure**.

### 11.8 Residue numbering

- The control is labelled **Residue numbering**. Its options are "4FE5 B (default reference)", "Rfam consensus positions" and every
  other chain.
- An info icon carries the requested explanation, and a screen-reader text states the active labels.
- The displayed-structure selector has its own short info.
- Numbering logic, fallbacks and lookup are unchanged (existing tests pass). Lookup now also finds Br/Bph by their two residues.

### 11.9 300-position cutoff for circular diagrams

- **Rule:** `CIRCULAR_MAX_CLEN = 300` in `lib/interactions.ts`; `showCircular(clen)` uses the dataset's Rfam consensus model length.
- **Above the cutoff:**
  - neither the family diagram nor the comparison diagrams are rendered;
  - a one-line note is shown;
  - the table, details, displayed-structure 3D view, comparison table, pair viewer, numbering, lookup and filters all remain;
  - no data are filtered or truncated.
- **Tests:**
  - vitest: 299, 300, 301;
  - Playwright: the RF00167 dataset served with `clen` 299, 300 and 301. It checks that the diagrams are present, present and absent, and
    that the full table, the 3D viewer and the comparison are there every time.

### 11.10 Default comparison pair

`defaultComparePair` (vitest + Playwright) chooses the pair as follows:
- **Source:** the stored geometric discrepancy D (`structure_clustering.json`), never recomputed.
- **Candidates:** only the clustered chains, which is the Structures tab's existing inclusion and coverage rule (76 chains, 57 shared
  positions).
- **Exclusions:** missing or non-finite values; pairs where either chain has no FR3D annotation; pairs where no interaction is
  comparable in both. No new cutoff is introduced.
- **Ties:** broken by the lexicographically smallest (A, B) seq_names, with A the smaller id.

**RF00167:** 5E54 B vs 9V50 A, D = 0.698 (the matrix maximum). The note "Showing the pair with the highest geometric discrepancy among
eligible comparisons." is shown while that automatic pair is active.

**Manual choice:** once the user picks a chain, the pair is kept for the page session, across Family / Compare switches and tab switches.
"Emphasise differences" is the default emphasis.

### 11.11 3D captions

The pair viewer caption is now short: "Locally aligned on G32 to compare the relative position of A33." The method, anchor, number of
fitted atoms, correspondence rule and RMSD are in an info icon. The viewers and the alignment algorithm are unchanged.

### 11.12 Tests (revision 2)

- **pytest:** 99 passed.
  - `test_interactions_rf00167.py` checks: Bp+Bs summary and pairs equal to the investigation; Br/Bph endpoints all verified with files
    for all 66 entries; four source kinds; raw BGSU Br/Bph rows kept in order and class; directional ids; phosphate-specific
    eligibility.
- **vitest:** 95 passed, including 16 interaction tests:
  - category-specific eligibility and reasons; both directions distinct;
  - category availability; grouping and explicit group occurrence; the default order;
  - cutoff boundaries; default pair (null skipped, unannotated excluded, ties, A < B);
  - the RF00167 default pair equals the matrix maximum; Bp/Bs unchanged;
  - gradient monotonicity, plus the earlier tests.
- **tsc:** clean.
- **Playwright:** `interactions.spec.ts`, 14 passed. Includes default Bp only, multiple categories and grouping, colour vs frequency,
  default and header sorting with reset, the default and manual pair, and the cutoff at 299, 300 and 301. All earlier workflows still
  pass: Present/Absent/Not comparable, supporting structure, lookup, the Stagno 5E54 B vs 5SWE comparison, the local fit and no refetch.
- **Full Playwright suite:** 173 passed, 0 failed (18.2 min); see §11.14.

### 11.13 Remaining limitations

- **Interpretation of the new classes.** Br/Bph class numbers are reported as annotated. The tab does not explain each class number
  (0–9) beyond the info icon; a link to the FR3D documentation could be added.
- **Low-resolution default pair.** The default pair is purely the maximum discrepancy. For RF00167 it involves 9V50 A, a 3.25 Å cryo-EM
  chain; no resolution criterion was added (none exists in the Structures tab rule). A family without `structure_clustering.json` gets
  no automatic pair: A is the displayed chain and B is left to the user.
- **Comparison grouping.** Status colours use the most different member's status. The individual members are listed in the expandable
  row and in the tooltip.
- **Base-pair-only filters.** While canonical or context filters are active, other categories are hidden by design.
- **Self contacts and between-chain contacts.** Intra-nucleotide Br/Bph rows (none in RF00167) are kept in the observations but not
  aggregated. Between-chain contacts are excluded, as before.
- The earlier limitations in §10 still apply.

### 11.14 Full suite and screenshots

- The full Playwright suite (173 tests, all tabs) passed with no failures or retries.
- The historical screenshots that the suite regenerates were restored.
- The Interactions screenshots are now:
  - `docs/screenshots/m23_interactions_family_default.png`;
  - `m23_interactions_four_types_multiple.png`;
  - `m23_interactions_compare_default_pair.png`.


## 12. Revision 3 — visual encoding, filtering and search

The data, frequencies, eligibility, statuses, numbering, supporting chains, default pair, local superposition and viewers are unchanged.
No dataset was regenerated.

### 12.1 Files

| File | Change |
|---|---|
| `frontend/src/lib/interactions.ts` | `CATEGORY_STYLE` (colour family, ramp, dash, pattern), `categoryColor`, `MULTI_RAMP` / `multiColor`, `parseSearch`, `matchesPositions`; the single-hue `conservationColor` was removed |
| `frontend/src/components/InteractionsTab.tsx` | category arcs; status-plus-pattern arcs in compare mode; legend in the side column; dropdown actions and empty state; unified `SearchBox` in both tables; toolbar lookup removed |
| `frontend/src/components/interactions/CircularDiagram.tsx` | `data-category`, `data-status` and `data-dash` on arcs; diagram up to 780 px wide |
| `frontend/src/lib/interactions.test.ts`, `frontend/tests/interactions.spec.ts` | new and updated tests |

The circle radius in `InteractionsTab.tsx` is now 204 in the 560 viewBox (it was 186), which leaves less blank margin.

### 12.2 Category conventions

| Type | Colour family (rare → conserved) | Pattern (SVG dash, round caps) |
|---|---|---|
| Bp base pairing | blue `#6baed6` → `#08306b` | solid |
| Bs base stacking | green `#74c476` → `#00441b` | dashed `6 3.5` |
| Br base–ribose | orange `#fd8d3c` → `#7f2704` | dotted `0.1 3.4` |
| Bph base–phosphate | purple `#9e9ac8` → `#3f007d` | dash–dot `7 2.6 0.1 2.6` |
| Multiple interactions | neutral grey `#a8a29e` → `#292524`, double line | solid |

**How the scales work**
- Each type has a 5-step ColorBrewer ramp, interpolated continuously by observed / eligible chains.
- Arc width is constant (1.8).
- The lightest step of every ramp has at least 2:1 contrast on white, so rare interactions stay visible (vitest).
- Darkness rises strictly with frequency (vitest).

**Palette check.** The four mid-tones (`#2171b5`, `#238b45`, `#d94801`, `#6a51a3`) pass the dataviz palette validator. Its worst
colour-vision-deficiency separation is 8.4 ΔE (orange/green, deutan), and the line pattern provides a second cue.

**Other markers.** The selected arc keeps its dark halo, the endpoint markers are unchanged, and the arc geometry is unchanged.

**Compare two:** colour still means comparison status (shared grey, only-here vermilion, reclassified pink, not comparable khaki
`#c9b98f`). Category is carried only by the line pattern, and multiple interactions by the double line. Not comparable was previously
dotted; it is now a distinct colour, so dotted means Br in every view. The legend states "Colour = status" and "Pattern = type".
Playwright checks every arc's colour against its status and its dash against its category.

### 12.3 Legend

The legend sits in the side column above the details.
- One line sample per type, drawn in its pattern and with its own light-to-dark ramp.
- The multiple-interactions double line and the selected-nucleotide marker.
- Each type (and Multiple interactions) has its own small 0% → 100% intensity bar on its legend row, so every scale is shown, not only
  the base-pair one. A short note explains that intensity is observed / eligible chain instances within each type, and that Multiple
  interactions use the "any" occurrence. (This replaced the single illustrative blue bar after review.)

### 12.4 Multiple interactions

- **Grouping:** unchanged in the data. Each member keeps its id, classes and frequency, and alternative classes within one type are not a
  "multiple" group.
- **Arc:** a double line in neutral grey, coloured by the explicit "any" occurrence. It is never drawn in one category's colour (Playwright
  checks the colour is not in any category ramp).
- **Inspection:** clicking the grouped arc, or its row, selects a member. The details panel lists every constituent with its type, class,
  own k/n and status, and the expandable table row shows the same.

### 12.5 Interaction-type dropdown

- **Contents:** four checkboxes, each with its line sample, plus **Select all** (every type that has data) and **Reset to Bp**. There is
  no "All interactions" checkbox.
- **Summary label:** `Bp`, `Bp, Bs`, …, `All types` when all four are on, or `none selected` (amber).
- **No types selected:** the diagram and the tables show "No interaction types selected" instead of suggesting that no interactions
  exist.

### 12.6 Unified search

- **Toolbar:** the lookup field and Find button were removed. One search field sits above each table, "Search interactions or residues".
  The query is shared by the two modes, and the lookup utilities are reused unchanged.
- **Residue queries:**
  - `G37`, `37`, `G37-C61` in the selected numbering, or Rfam positions under Rfam numbering;
  - `Rfam 41` and `Rfam 41-66` always as Rfam positions.
- **Text queries:** class codes (`cWW`, `s35`, `9BR`, `0BPh`), region names (`J2/3`, `L2`, `P1`) and type names.
- **Matching:** a single residue matches every interaction involving it. A pair matches its group in either direction, so a pair with
  several types shows all of them.
- **What the search does:** it only filters the table. The diagram, the selection and the 3D view are untouched while typing.
  **Enter**, or clicking a row, is the explicit selection; it updates the table, diagrams, details and 3D view.
- **Hints and errors:**
  - a hint states the resolved residues and Rfam positions;
  - "n more in types not shown" appears when the pair has annotations in disabled types;
  - invalid or ambiguous residues give an explicit message ("37 of 4FE5 B is G37, not A37", "No residue 9999 in 4FE5 B") and an empty
    table;
  - type filters and numbering are preserved.
- **Large families:** the search works the same way above 300 positions, where no diagram is shown.

### 12.7 Layout

- The diagram uses more of its card: a larger radius in the same viewBox, up to 780 px wide.
- The side column holds the legend and a compact one-line prompt when nothing is selected, so there is no large empty panel.
- The details layout is unchanged when an interaction is selected.
- The 300-position cutoff is unchanged.

### 12.7a Later wording changes (same day)

- The Interactions tab now comes before Annotations: Family Overview · Structures · Interactions · Annotations.
- In Compare two, the "Emphasise differences" button is now **Show differences** (still the default).
- The status chips show counts in parentheses, e.g. "Only in A (16)".
- The status "Reclassified" is renamed **Different class**: same interaction type between the same nucleotides, different FR3D class
  in A and B (e.g. cWW vs tWW). The internal key (`reclassified`), the logic and the counts are unchanged. The chips and legend carry a
  short tooltip for each status. "Multiple interactions" (different types between two nucleotides) remains a separate notion.
- In both modes, the interaction table and the 3D viewer share one row grid (CSS subgrid: header, a 460 px body, caption). They have the
  same width (equal columns), height and top edge. The compare viewer gets a "3D superposition" heading so the header rows match.
  Playwright checks this.

### 12.8 Tests (revision 3)

- **pytest:** 99 passed (no backend change).
- **vitest:** 98 passed, including 19 interaction tests:
  - per-category family and pattern; distinct dashes; endpoints and monotonic intensity per ramp;
  - at least 2:1 contrast at 0% for every ramp and for the multiple-interactions ramp;
  - search parsing (residues, pairs, Rfam, class codes and regions as text, errors); position matching.
- **tsc:** clean.
- **Playwright:** `interactions.spec.ts`, 15 passed. New or updated:
  - category encoding of arcs and legend; selected halo; Select all / Reset to Bp / none-selected empty state;
  - the grouped arc's neutral colour and its constituents;
  - the toolbar lookup is absent;
  - unified search: pair, residue, Rfam, class code, errors, hidden-type note, Rfam numbering; typing does not select, Enter does,
    with 3D residues checked;
  - compare arcs: colour = status, dash = category.
  - Every earlier regression still passes: the Stagno comparison, local fit, absent/supporting workflow, default pair and the cutoff at
    299, 300 and 301.
- **Full Playwright suite:** see §12.10.

### 12.9 Remaining limitations

- **Orange in two meanings.** Br orange in family mode and the vermilion "only in this structure" in compare mode are both orange-ish.
  They never appear in the same view, and each view's legend states its meaning.
- **Small dash patterns.** The dotted and dash–dot patterns are less distinct on very short arcs between adjacent nucleotides. The type is
  always in the tooltip, the table and the details.
- **Shared search.** The search query is shared by both modes, so a search typed in Family conservation also filters the comparison
  table. The hint under the field shows the active query.
- Earlier limitations (§10, §11.13) still apply.


## 13. Revision 4 — pagination, split data files and 16S

### 13.1 Table pagination

Both tables (family and comparison) are paged.
- **Page size:** 15 top-level rows by default, chosen to fill the 460 px box beside the 3D viewer. "Rows per page" offers 15, 50 or All,
  and the choice is shared by the two tables.
- **Controls:** « First · Previous · Page [n] of N · Next · Last ». The page number is an input: type a page and press Enter, and it is
  clamped to the range.
- **Selection and paging:**
  - Paging never changes the selection.
  - A selection made elsewhere (a diagram arc, or Enter in the search) shows its page.
  - A new sort, search, filter or page size starts at **page 1**, so the new order is read from the top; the selection is kept
    (corrected after review: it previously jumped to the selected row's page, which made a sort look wrong). The sort arrow stays on the
    same line as the column name.
- **Groups:** an expanded multiple-interaction group keeps its members on its page.
- **Filling the box:** on a full page the rows stretch to fill the box; a short last page keeps natural row heights.
- **Tests:** Playwright covers first, previous, next, last, a typed page (including clamping), selection persistence, jump-to-selection,
  the page-size choice and the paged comparison table.

### 13.2 Split data files (schema 1.2)

The main `interactions.json` no longer embeds every raw observation.
- **Main file, per chain:** the compact index the tab needs:
  - `observed`: indices into `pairs`, by FR3D class;
  - `near`: near annotations on consensus positions;
  - `n_off_consensus`.
- **Per-chain files:** the full records (BGSU unit ids, `label_seq_id`s, Rfam location, mmCIF check, altlocs), in
  `interactions_by_instance/<seq_name>.json` (schema `interactions_instance`), referenced by `observations_files`.
- **Checks:**
  - pytest asserts that the index equals the per-chain records;
  - RF00167 `pairs` and every observation are identical to schema 1.1;
  - the RF00167 main file went from 1.3 MB to 414 KB.
- **Compatibility:** the frontend still reads 1.0/1.1 files that have inline observations.

### 13.3 16S rRNA (RF00177)

```
python3 -m rna_explorer interactions RF00177          # 35 s with a warm cache; BGSU files fetched once into cache/
```

The step is the same, with no family-specific code. It uses the existing per-chain ModelServer coordinates of the selection build and 4V9D
chain Y (BA) as the default numbering reference (the family display representative, with all 1,513 matched positions).

| Category | Interactions | Consensus pairs | Endpoints verified against mmCIF |
|---|---|---|---|
| basepair | 36,346 | 794 | 145,980 / 145,980 |
| stacking | 106,755 | 1,910 | 427,552 / 427,552 |
| baseribose | 41,114 | 808 | 83,074 / 83,074 |
| basephosphate | 15,982 | 407 | 32,326 / 32,326 |

**Coverage.** All 58 entries have all four BGSU files, and all 60 chains are FR3D-annotated. No mapped endpoint disagrees with the
mmCIF.

**Size and git.** The main file is 4.6 MB uncompressed. The per-chain files total 17 MB and are git-ignored with the other
regenerable RF00177 outputs over 10 MB (`.gitignore`); they are needed only for provenance.

**In the tab**
- **Layout:** 1,533 consensus positions is above the 300 cutoff, so there are no circular diagrams. The family table, details, Mol* view,
  comparison table, pair viewer, numbering, search and filters all work.
- **Measured locally:**
  - the tab is ready in about 2.5 s;
  - base pairs alone are 794 rows on 53 pages, and all four types give 3,196 rows;
  - search responds in about 0.2 s, and switching on all types takes about 0.3 s.
- **Default pair:** 8BF7 chain 2 (ZA) vs 8UZ3 chain 2 (B), the maximum stored discrepancy among annotated clustered chains.
- **Regions:** RF00177 has no validated structural-region track, so Region shows "–" and the region filter is empty.
- **Tests:** a Playwright test covers the RF00177 tab (no diagram, paging, numbering-based search, details, comparison and two-chain
  viewer). A pytest checks the RF00177 schema, endpoint verification, the default reference and the index bounds.

### 13.4 Other changes in this revision

- The legend shows one 0% → 100% bar per type (and for multiple interactions) instead of a single blue bar.
- "Reclassified" is now **Different class** (see §12.7a).

### 13.5 Limitations

- **Default 16S pair.** It is chosen purely by maximum discrepancy, like RF00167's. It may pair chains from different organisms, which
  gives many differences (e.g. 530 only in A, 319 only in B for Bp+Bs+Br+Bph).
- **First load.** The 16S main file is 4.6 MB uncompressed; it is served compressed and loaded only when the tab opens.
- **RF00059 and RF00162** can be built the same way but have not been generated yet.


### 13.6 Annotation coverage and the default pair (review of 5E54 B vs 9V50 A)

**What the review found.** The earlier default pair, 5E54 B vs 9V50 A, shared only 14 base pairs. The count is correct for the BGSU data,
but it understates the shared structure.
- **Sparse source annotation.** BGSU FR3D annotates only 19 base pairs in 9V50 chain A (17 cWW), a 71-nt guanine-II aptamer at 3.25 Å.
  5E54 B has 32.
- **Paired but unannotated.** Some unannotated pairs are paired in the coordinates (hydrogen-bond distances in the PDBe mmCIF):
  - G3–U70: 2.3 and 3.2 Å
  - G16–U28: 2.5 and 3.0 Å

  Even the annotated G2–C71 has compressed contacts (2.3–2.6 Å). The model's geometry is outside FR3D's cutoffs, not unpaired.
- **No near base pairs to fall back on.** BGSU exports no near base pairs, so these pairs appear as plain "absent" / "only in A".
- **"Only in A" therefore mixes two things:** real subclass differences (guanine-II register, e.g. C53–G72 only in B) and annotation gaps.

**Changes, decided with the user**
1. **Annotation coverage in Compare two** (`helixCoverage`, `medianHelixCoverage`):
   - **Definition:** for each chain, the Rfam consensus (SS_cons) helix base pairs FR3D annotates, out of those both of whose bases are
     resolved and complete.
   - **Display:** shown for A and B with the family median. A chain below 80% is highlighted with "sparse annotation; 'only in'
     differences may be annotation gaps". The 80% is a display cue only and changes no classification.
   - **RF00167:** the median is 95%. Only 9V50 A (13 of 20, 65%) and 9NBC A (13 of 19, 68%) fall below 80%. Chains with no evaluable
     helix pair (e.g. the 9NBH scaffold chains) have no value.
2. **Default pair resolution limit:** both chains at ≤ 3.0 Å (`DEFAULT_PAIR_MAX_RESOLUTION`). Chains without a resolution are excluded,
   and the note says "(both chains ≤ 3.0 Å)".
   - **RF00167:** the default is now **5E54 B vs 9LJN A** (D = 0.694; coverage 100% and 95%). It was 5E54 B vs 9V50 A, D = 0.698.
   - **16S:** recomputed with the same rule.
   - **Tests:** vitest checks the resolution rule, its boundary (3.0 Å included, 3.25 Å excluded), and the unrestricted result.
     Playwright checks the coverage note and the flag for 9V50 A.

**Near base pairs remain a possible next step.** They could be computed locally with BGSU's fr3d-python, so "near" can be told apart
from "absent".

### 13.7 Default table order changed to PDB entries

The user changed the family table's default order to **supporting PDB entries, descending**. Ties fall back to conservation frequency,
then eligible chains, then Rfam positions, all on numeric values. A group is ordered by its largest member entry count. This supersedes
the frequency-first default in §11.7. Column sorting and "Default order" behave as before, and vitest and Playwright check the new order.

### 13.8 Per-class breakdown

When an interaction's chains carry more than one FR3D class (e.g. U25–A45: cWW and tHS), its class cell has a ▸ toggle. The toggle does
not change the selection. Expanding it shows one sub-row per class:
- the chains with that class out of the interaction's eligible chains;
- that class's PDB entries;
- whether the displayed chain has this class ("this class", "other class (…)", Absent or Not comparable).

The breakdown is computed from the per-chain index (`classBreakdown`). It also works for members of a "Multiple interactions" group.
Playwright checks the toggle and the per-class counts.

### 13.9 All families built

`python3 -m rna_explorer interactions <ACC>` was run for every family in the app. The code is the same for all five, and every BGSU file was
retrieved for every entry.

| Family | Consensus length | Chains / entries | Interactions | Consensus pairs | Endpoints verified | Default numbering | Circular |
|---|---|---|---|---|---|---|---|
| RF00167 Purine riboswitch | 102 | 84 / 66 | 12,830 | 262 | 47,544 / 47,550 | 4fe5_1_B_A | yes |
| RF00059 TPP riboswitch | 105 | 38 / 28 | 6,816 | 226 | 23,808 / 23,808 | 2gdi_1_X_A | yes |
| RF00162 SAM-I riboswitch | 108 | 29 / 28 | 5,684 | 242 | 20,328 / 20,328 | 5fjc_1_A_A | yes |
| RF00005 tRNA | 71 | 60 / 59 | 8,579 | 355 | 30,802 / 30,802 | 1ehz_1_A_A | yes |
| RF00177 16S rRNA | 1533 | 60 / 58 | 200,197 | 3,919 | 688,932 / 688,932 | 4v9d_1_BA_Y | no |

The only endpoint mismatches are RF00167's 6 (1Y27 C81, `CCC` in the mmCIF, see §8).

### 13.10 Incomplete source annotation (schema 1.3)

The scientific review (`experiments/interactions_review/REVIEW.md`) found that BGSU annotates RF00177 7PJS chain a only up to about
nucleotide 650. The chain is fully modelled, so 2,343 interactions had been counted absent.
- **Rule:** the build now stores `annotation_gaps` per chain: [start, end] Rfam runs of ≥ 15, and ≥ 8% of the chain's, consecutive
  eligible positions without any annotation.
- **Effect:** interactions touching a gap are annotation-unavailable for that chain. In the frontend they are not comparable, with the
  reason "BGSU annotation incomplete for this chain (no annotation at Rfam a–b)".
- **Scope:** only 7PJS chain a (Rfam 652–1531) has a gap. The RF00167, RF00059, RF00162 and RF00005 pairs are unchanged.
- **RF00177:** 2,344 pairs were corrected, and 7PJS helix coverage went from 38% to 95.5% (169/177).
- **Tests:** pytest and vitest.
