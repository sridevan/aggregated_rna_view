# Mol* 3D viewer: molecular context and "Show ligands"

Date: 2026-10-19. Scope: the PDBe Mol* viewer (`frontend/src/components/PdbeMolstarViewer.tsx`), used by the Family Overview and by the
Interactions tab's family mode. Coordinates, Rfam mapping, numbering, annotations, interactions and the Structures tab's 3Dmol viewer are
unchanged.

## 1. What was there before

- **Rfam region (default):** one cartoon of the chain's Rfam-mapped region (`label_asym_id`, `label_seq_id` from `map_start` to
  `map_end`), slate `#64748B`. Every other component is hidden.
- **Complete structure:** PDBe Mol*'s default components are shown: Polymer (all chains, RNA and protein together), Ligand, Non-standard,
  Ion (Water hidden).
  - The whole Polymer component was painted slate, so in ribosomes and RNA–protein complexes the selected RNA could not be told apart
    (before: `docs/screenshots/m24_before_16s_4v9d_complete.png`, `m24_before_trna_1c0a_complete.png`).
  - Switching modes re-fitted the camera to every visible atom, so it zoomed out to the whole entry.
- **Ligands:** loaded by PDBe Mol* into its own "Ligand" component (not excluded), visible only in Complete structure. For very large
  entries PDBe Mol* may not create the component at all (4V9D: Polymer only).

## 2. Complete structure: visual hierarchy

Three cartoon layers replace the single slate Polymer component. They are built once per loaded structure and then only shown or hidden:
- **Selected RNA chain** (`label_asym_id`): slate `#64748B`, fully opaque. This is the same colour as Rfam region.
- **Other nucleic-acid chains** (Mol* `nucleic` query minus the selected chain): grey `#a8a29e` at 25% opacity.
- **Protein chains** (Mol* `protein` query minus the selected chain): lighter grey `#d6d3d1` at 20% opacity.

Other components:
- Non-standard residues and Ions keep PDBe Mol*'s atom-coloured representation; ions are small spheres.
- Water stays hidden.
- Small-molecule ligands are governed by Show ligands (§3).
- The same rules apply to every family. The opacities were lowered from 45%/35% after review on 4V9D, which has two 70S ribosomes per
  model.

**Camera**
- A mode switch keeps the camera exactly (verified: identical camera state before and after, in every test case).
- The initial view and Reset orientation / Fit now fit the **selected RNA chain** in Complete structure, not the whole deposited model.

**Selection**
- Nucleotide and whole-loop selection (violet focus representation) are unchanged, and both are kept across mode switches.
- A selection made before the layers existed is re-applied to them.

**Instance switch:** the layers are rebuilt for the new chain, and no stale layer remains (tested 4V9D → 7K00).

**Performance**
- The structure is loaded once and only visibility changes afterwards.
- Building the layers on the first switch to Complete structure takes 1.3 s for 4V9D (33,015 + 166,356 + 90,756 atoms in the three
  layers), 0.7 s for 1C0A, and under 0.1 s for riboswitches.
- No surfaces, and no additional Mol* instance.

## 3. "Show ligands"

**The toggle**
- A checkbox next to Rfam region / Complete structure, **off by default**.
- It works in both modes (independent of the mode), and its state is kept across mode switches.
- It is reset to off and re-evaluated on another RNA instance.

**Which ligands**
- *Why not the existing annotations:* the PDBe binding-site annotations behind `ligands.json` are entity-level. They give bound-molecule
  ids (e.g. `bm6`), contacted residue indices and a "best chain", but not the ligand instance (its chain and residue) or the RNA chain
  that makes the contact. They therefore cannot say which ligand copy belongs to the selected chain.
- *Spatial rule instead:* the toggle uses a documented spatial association, computed per instance in the loaded coordinates. A
  small-molecule residue qualifies when it is a non-polymer entity that is not water, ion, saccharide or lipid, and has any atom within
  **4.0 Å** of any atom of the selected RNA chain (Mol* `filter.within`, one atom group per residue).
- *No identity-based association:* the CCD identity is never used to associate a ligand.
- *Exclusions:* CCDs that the family's own classification (`ligands.json` `by_ccd`, from `annotation_config.json`) marks as **additive**
  (crystallisation agents, e.g. GOL, PEG, SO4, PUT, SPD) or **metal_ion** are excluded. That classification is loaded only when the toggle
  is first switched on.

**What is drawn**
- Ligands are sticks (equal atom and bond radii) with carbons in the RNA chain's colour and the other atoms in CPK colours (N blue, O red,
  P orange, S yellow). This was changed after review; previously they were ball-and-stick with grey carbons.
- In Rfam region, only the ligands near the selected chain are shown.
- In Complete structure, the other ligands of the model are shown too.
- PDBe's default Ligand component is retired in favour of these two layers, so with the toggle off Complete structure shows no
  small-molecule ligands.

**Ions** are not part of the toggle. They keep their previous behaviour: small spheres in Complete structure, hidden in Rfam region.

**Unavailable:** when no candidate lies within 4 Å, the checkbox is disabled and reads "(none nearby)". The hover explains: "No
small-molecule ligand within 4 Å of this RNA chain in the loaded coordinates". When candidates exist but all are ions or additives, the
label reads "(ions/additives only)" (renamed when the controls moved to the toolbar).

**Labels:** the checkbox label stays "Show ligands" (after review the shown CCD ids moved from the label to its tooltip, "Shown: …"); the instance identifiers ("CCD auth_chain auth_seq") are kept in the control's data and
the test hook.

### Results

| Structure | Candidates within 4 Å | Shown (instance) | Notes |
|---|---|---|---|
| Purine riboswitch 4FE5 B | 1 | HPA B101 | the hypoxanthine of chain B only; chain A's copy is not associated |
| Purine riboswitch 1Y26 X | 1 | ADE X | |
| TPP riboswitch 2GDI X | 1 | TPP X100 | chain Y's TPP not associated |
| TPP riboswitch in RNA scaffold 9C6K X | 1 | TPP X203 | |
| 16S rRNA 7K00 A | 3 | PAR A1601 (paromomycin) | 2 polyamine additives excluded |
| 16S rRNA 4V9D Y | 0 | – (disabled) | PDBe Mol* loads no ligand component for this entry |
| tRNA 1EHZ A | 0 | – (disabled) | |

- Toggling takes 60–90 ms.
- The camera never moves.
- Selection highlighting is unaffected.

## 4. Limitations

- **Spatial association is a proxy for binding.** A 4 Å contact is a geometric criterion, not a curated binding assignment. Ligands
  bridging two chains are associated with both.
- **Unclassified CCDs.** CCDs absent from the family's PDBe-derived classification (not in any PDBe binding site of the family) are not
  filtered. An unlisted crystallisation agent within 4 Å would therefore be shown.
- **What PDBe Mol* loads.** For very large entries PDBe Mol* may load no ligand or ion components (4V9D), so the toggle reports none
  nearby.
- **Ions** are not covered by the toggle; a separate ion control was not added, to keep the controls minimal.

## 5. Tests

**`tests/overview-complete-structure.spec.ts`** (4 tests)
- 1C0A tRNA–synthetase:
  - Rfam region is the default;
  - layer colours and opacity; non-standard residues visible; water, Polymer and Ligand hidden;
  - the camera unchanged on a mode switch, and the selection and its focus unchanged;
  - Reset fits the selected chain;
  - back in Rfam region, the layers are hidden.
- 4V9D 16S: subdued rRNA and protein layers; switching to 7K00 rebuilds the layers.
- Show ligands on 4FE5:
  - off by default; HPA B101 only; the camera kept;
  - kept in Complete structure; the default Ligand component hidden;
  - off again on 1Y26, and ADE shown when enabled.
- Show ligands elsewhere: 2GDI TPP X100; 7K00 PAR only (additives excluded); 1EHZ disabled with "none nearby".

**Updated tests**
- `overview-structure-colour.spec.ts`: the selected-chain layer counts as the slate RNA.
- `overview-loop-refinements.spec.ts`: Polymer and Ligand are replaced in Complete structure.

**Full Playwright suite:** 184 passed, 0 failed (21.1 min).

**Screenshots** (`docs/screenshots/m24_*`)
- Before / after, Complete structure: 16S 4V9D and tRNA 1C0A.
- Selected nucleotide in 16S Complete structure; selected tRNA loop.
- TPP scaffold 9C6K.
- Show ligands off / on: purine 4FE5 in Rfam region; 4FE5 in Complete structure; TPP 2GDI; 16S 7K00; tRNA 1EHZ (none nearby).
