# Ligands tab (2026-10-22)

A dedicated **Ligands** tab between Interactions and Annotations. It holds:
- a PDBe-KB-style gallery of RNA-contacting ligands;
- the per-ligand contact-residue tracks;
- the fpocketR predicted-pocket track, moved out of the Annotations tab.

The change is frontend only. No backend code, data file, schema or calculation was changed. Basis: `experiments/ligand_annotation_investigation/REPORT.md`.

## Layout

**Tab order:** Family Overview · Structures · Interactions · **Ligands** · Annotations.

### A. RNA-contacting ligands (`LigandsTab.tsx`)

**Which ligands appear.** The default ligand set of `ligands.json`: class `ligand`, not flagged. This is the same set as the residue tracks, in the same order (most PDB entries first, then CCD id; `defaultLigandCcds`).

**Each card shows:**
- the 2D depiction;
- CCD id and PDBe name;
- "Found in *n* PDB entries" (`by_ccd.n_entries_in_default_track`); the hover lists the entries from the default-set `bound_molecules`;
- a "PDBe ligand page ↗" link.

No ligand-copy counts are shown.

**Interaction:**
- The first 12 cards are shown, then "Show all".
- Search matches the CCD id or the name, case-insensitively. A search always shows every match and filters the gallery only, never the tracks.
- Selecting a card (image or CCD id) highlights that ligand's track row, re-opens its category if collapsed, and scrolls the row into view. Selecting it again clears the highlight. Values are never changed.

### B. Ligand-contact residues and C. Predicted binding pockets

- One `TrackBrowser` with two titled categories on one axis, so contacts and pockets stay aligned and share the position numbering and zoom.
- Their marks, colours, tooltips, pinned details, missing-data handling and frequencies are the existing ones.
- "View in Family Overview" now records the source tab ("Opened from Ligands …").

### About these annotations

A collapsed section listing: source, contact identification (Arpeggio; computed contacts, not proof of binding), annotation level (RNA entity), coverage (preferred assembly), mapping, alternative numbering, frequency unit (distinct PDB entries), limitation, and pockets.

### Annotations tab

- It now shows Solvent accessibility, RNA loop annotations and Nonstandard nucleotides.
- Its summary points to the Ligands tab.
- Its download menu keeps loops, tracks, SASA and nonstandard-nucleotide files. The ligand and pocket files moved to the Ligands tab's menu.

## Implementation

**Shared track viewer.** `AnnotationsTab.tsx` exports `TrackBrowser`: the former body of the Annotations tab, moved unchanged. It contains:
- zoom and the horizontal window;
- the axis with position numbering;
- the legend built from the scales in use;
- collapsible categories;
- the hover tooltip and pinned details.

**Track builders.** The per-category builders are exported: `loopTracks`, `ligandTracks` (former `buildTracks` ligand part, unchanged), `pocketTrack`, `sasaTrack`, `modificationTracks`.

**TrackBrowser additions:**
- `highlight` / `revealKey` for the gallery link;
- optional ⓘ texts next to the numbering selector and the legend;
- the `origin` tab for View in Family Overview.

**Navigation.** `NavigationOrigin.from` is `annotations | ligands`.

**Depictions.**
- Source: PDBe `https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/{CCD}_200.svg`, the images PDBe-KB uses, loaded as SVG through `<img>` (no CORS needed).
- Fallback: RCSB `https://cdn.rcsb.org/images/ccd/unlabeled/{first char}/{CCD}.svg`.
- If both fail: a "No depiction available" placeholder, with a skeleton while loading.
- Checked on 2026-10-22: all 72 default-set CCDs of the six families return `200 image/svg+xml` from PDBe. RCSB also served ADE, HPA, GUN, 6AP, A1JAI and SAM. `{CCD}.svg` without `_200` returns 404 at PDBe.
- No image is generated locally.

**Fixed while testing.** The hover tooltip was rendered inside a `space-y` container. Its margin changed the page height under the pointer near the bottom of the page, which made the tooltip flicker. It is now rendered outside that container, as the Annotations tab did before.

## Texts

The user-supplied texts (§8.1–8.4 and §7) are used verbatim:
- §8.1, gallery ⓘ;
- §8.2, ligand-contact residues ⓘ;
- §7, pocket ⓘ;
- §8.3, legend ⓘ, plus one sentence on pocket colours;
- §6 + §8.4, numbering ⓘ.

All are listed in `docs/TOOLTIPS.md` § Ligands.

## Tests

- **vitest: 104 passed.** New: default set and order, distinct-entry counts, search, depiction and ligand-page URLs.
- **Playwright, new `tests/ligands-tab.spec.ts` (9 tests):**
  - tab order;
  - gallery = default set in track order, names, entry counts and entry lists, no copy counts;
  - PDBe depiction actually loaded;
  - RCSB fallback and placeholder (network routes aborted);
  - tracks in gallery order with unchanged `data-value`s;
  - HPA tooltip and pinned details (3 of 12 at consensus 79, supporting 4FE5/4FEL/4FEN);
  - pocket tooltip;
  - every ⓘ text;
  - the methodology section;
  - selection highlight, category reveal and scroll with values unchanged;
  - search;
  - numbering changing labels only (all marks and values identical, shared with Annotations);
  - "Opened from Ligands";
  - independent category scroll and the four downloads;
  - RF00177 and RF00011 (no pockets; empty ligand set);
  - 390 px width.
- **Updated specs:** `annotation-tracks`, `annotations-ui`, `annotations-zoom`, `pocket-track`, `modification-track`, `sasa-track`, `info-popover`, `terminology`, `ux-refinements`, `structures-selection`, `family-rf00005`, `family-rf00177` and `screenshots`. Ligand and pocket checks now run in the Ligands tab; Annotations checks expect three categories and six downloads.
- **Also updated:** `annotations-download` (one test per tab menu; each transient popover is dismissed before the next hover, because it covers the next row's icon), `family-rf00011`, `family-rf00059`, `family-rf00162`.
- **Full Playwright suite: 201 passed** (2026-10-22). The screenshot spec rewrites `docs/screenshots/`; those files were restored to the committed versions, so no screenshots were regenerated.

## Limitations (unchanged data)

- **Contacts are entity-level PDBe annotations of preferred-assembly ligand copies** (see the investigation report). Selecting an author numbering never filters to that chain.
- **The terminus flag still excludes tRNA synthetase substrate analogues** from the default set, so they are not in the gallery either.
- **Names come from PDBe** and use the existing title-casing, e.g. "9h-purine-2,6-diamine" for 9H-PURINE-2,6-DIAMINE.
