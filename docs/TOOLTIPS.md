# Tooltip text: Aggregated RNA 3D View

All tooltip and information text in the frontend as of 2026-10-19, extracted from `frontend/src/`.

The wording follows the editorial review of 19 October 2026 (39 replacements), which is now applied in the frontend. Algorithms,
thresholds, counting units and data-derived placeholders are unchanged. Method details removed from two shortened tooltips are kept in
collapsed details blocks, listed with those tooltips below.

**Kinds of tooltip**
- **ⓘ** is an information icon (`InfoTip`): it opens a popover on hover or focus, and a click pins it.
- **Hover** is a hover card or a native `title` tooltip.

**How the text is written here**
- Fixed text is quoted verbatim.
- Data-driven text is written as a template, with `{placeholders}` standing for values filled in from the data.
- Plain button labels that only repeat the visible text (e.g. "Zoom in") are listed briefly.

---

## Shared

| Where | Kind | Text |
|---|---|---|
| Position numbering selector (Family Overview, Annotations) | Hover | Changes residue labels only; alignment, aggregation and selections use Rfam consensus positions. |
| Zoom buttons (Family Overview tracks, Annotations) | Hover | Zoom out · Zoom in · Fit all Rfam positions |
| 3D superposition viewers (Structures, Interactions compare) | Hover on a residue | `{PDB}.{author chain} {residue name}{author number}{insertion code} (author)`, e.g. 5J7L.AA G230 (author) |

---

## Family Overview

### Statistics and cards

**Representative structure, ⓘ "How the representative is chosen"**
> Selected using {rule}. A coverage-based representative is retained as fallback.

The default `{rule}` is: "coverage representative: clan-best instance with most consensus-matched residues with coordinates, then best
resolution".

**Rfam consensus positions (stat), hover**
> Match-state positions in the Rfam covariance model.

**Statistics, ⓘ "How these statistics are defined"**
> RNA instances: {n} mapped chains ({n validated} validated); chains from the same PDB entry count separately.
>
> Unique sequences: {n} distinct mapped-region sequences.
>
> Mean sequence occupancy: average fraction of {n} unique sequences with a residue across {clen} Rfam positions.
>
> Mean structural occupancy: average fraction of {n} RNA instances with coordinates across {clen} Rfam positions.

**Experimental methods, ⓘ "How experimental methods are counted"**
> Counts {n} mapped RNA chains, including multiple chains from one PDB entry. Percentages use RNA instances; sequence statistics use {n} distinct mapped sequences.

**Species distribution, ⓘ "About the species statistics"**
> Species are assigned from source-organism NCBI Taxonomy records, not inferred from RNA sequences. Strains and subspecies are grouped by species; synthetic, unresolved and unknown sources remain separate. Counts use PDB entries (once per species) or RNA instances; only the top 10 identified species are listed.

**Species row, hover**
> {species}: {value} {unit} = {percent} % of {total} ({n chain instances | n PDB entries}). Source names: {names}

For the "Other" row, the last part reads "Species: {names}".

### Consensus tracks

**Sequence coverage bar, hover**
> Sequence coverage · Rfam {pos}
> {n}/{N} unique sequences ({percent})
> (when relevant) Deleted in {n} · not covered by the alignment in {n}

**Structural coverage bar, hover**
> Structural coverage · Rfam {pos}
> {n}/{N} RNA instances ({percent})
> (when relevant) {n} instance(s) contain the nucleotide in sequence

**Selected RNA chain track, hover**
> Selected RNA chain · {PDB chain}
> - Rfam position: {pos}
> - {chain} author residue: {residue} *(or "Author residue" when the chain is the numbering reference)*
> - {reference} reference residue: {residue} *(only when another chain is the numbering reference)*
> - Status: Coordinates available | Coordinates unavailable | Deletion | Outside aligned region | Not available in the mapping

The default chain tooltip, when no chain is selected, shows "Selected chain · Position {pos}" with one of:
- Coordinates available
- Coordinates unavailable, "Nucleotide present in the deposited sequence"
- Deleted in chain, "No nucleotide aligned to this position"
- Outside aligned region
- Not available in the mapping

With a numbering reference, the title reads "Rfam position {pos}".

**Selected RNA chain track label, hover card**
> **Selected RNA chain**
> Shows nucleotide availability at Rfam match positions. +N markers indicate insertions between positions.
> Symbols are explained in Legend & help.

**Insertion marker, hover**
> {Insertion | Local-end residues} · +{n} nt
> {boundary text}
> Sequence: {sequence}
> {chain} author residues: {residues}

**RNA loop annotations, ⓘ**
> BGSU hairpin, internal and junction loops in this RNA chain. Select a block to highlight the full loop, including discontinuous segments.

**Loop block, hover**
> {loop type} ({HL | IL | J3 | J4 | Jn})
> {BGSU loop id}
> {n} nucleotides · {k} segments *(· {m} in another chain)*
> Merged from {loop id} + {loop id} *(BGSU merged loops only)*
> Motif: {group} ({BGSU annotation}) · Motif Atlas {release} · {loop is a group member | homologous/geometric assignment via {loop id}} *(only when BGSU assigns a group)*

Loops that a BGSU merged loop supersedes are not drawn; the merged loop names them.

**Reference and alignment details, ⓘ**
> Rfam consensus sequence, reference secondary structure and alignment confidence.

**Alignment confidence bar, hover**
> Alignment confidence · Rfam {pos}
> Mean cmalign PP {value} ({high | medium | low})
> {n} aligned unique sequences

**Consensus sequence and secondary structure rows, hover**
> Consensus sequence (RF) · Position {pos}: {character}
> Secondary structure (SS_cons) · Position {pos}: {character}

### 3D viewer (PDBe Mol\*)

**Show ligands (checkbox), hover**
- With candidates: "Show ligand instances within 4 Å of the selected RNA as atom-coloured sticks. Complete structure also shows other model ligands. Excludes water, ions and crystallisation additives."
- Without candidates (disabled, greyed out): "No small-molecule ligand within 4 Å of this RNA chain in the loaded coordinates". The label
  stays "Show ligands"; this note is in the hover only, so the toolbar stays on one line.
- While ligands are shown, the hover adds a line "Shown: ADE" (the CCD ids shown; they are not listed next to the checkbox).
- When ligands are shown but only ions or additives are nearby, the hover adds a line "Only ions or crystallisation additives are within 4 Å."
- This control appears in Family Overview and in the Interactions tab's family-mode viewer. The Structures tab and Compare two have their own
  version (see 3D structure / 3D superposition).

**Reset view (button), hover:** "Reset orientation and fit the view"

There are no tooltips on the viewer itself.
- **Rfam region:** shows only the chain's Rfam-mapped region.
- **Complete structure:** shows the whole PDB entry.

*Viewer details* (a panel, not a tooltip) lists the 3D alignment:
- reference structure;
- superposition (superposed on the reference, or shown as deposited when the fit is insufficient);
- RMSD;
- matched atom pairs;
- alignment quality.

---

## Structures

### Table (RNA instances)

| Column | Kind | Text |
|---|---|---|
| PDB | Hover (link) | Open PDB entry {PDB} at PDBe |
| Chain | ⓘ | Deposited author chain identifier (auth_asym_id). The label asym ID identifies the chain internally. |
| Chain | Hover (link) | Open molecule (entity {n}) of {PDB} at PDBe |
| Structural group | ⓘ | Groups RNA instances by 3D similarity under the selected metric; groups do not necessarily represent biological states. Select a badge to explore a group. |
| Structural group "not clustered" badge | Hover | Not clustered: {reason} |
| "medoid" badge | Hover | Representative (medoid) of its structural group |
| "not in matrix" tag | Hover | Not clustered: no pairwise geometric discrepancy, so not in the matrix and without a seriation position; listed after the ordered instances |
| Title, Organism, Method | Hover | the full title, organism or method |
| Res. (Å) | ⓘ | Reported experimental resolution (Å), where available. Lower values generally indicate finer structural detail. |
| Length (nt) | ⓘ | Length in nucleotides of the PDBe/Rfam-mapped RNA region (hover a value for its polymer label_seq_id range). |
| Length value | Hover | RNA region: {start}–{end} |
| Seq. cov. % | ⓘ | Consensus positions of the Rfam model to which this chain contributes a match-state residue, divided by the model length ({clen}). Insertions and local-end residues do not count. |
| Struct. cov. % | ⓘ | Match-state residues that have experimental coordinates in this chain, divided by the model length ({clen}). Lower than sequence coverage when parts of the chain are unmodelled. |
| Avg PP | ⓘ | Mean cmalign posterior probability of residue placement in the family alignment. Higher values indicate greater alignment confidence. |
| Rfam bits | ⓘ | cmsearch bit score against the Rfam covariance model. Higher scores indicate stronger sequence–family match support; scores depend on alignment length. |

### Structural Dissimilarity

**Panel, ⓘ "About the structural dissimilarity matrix"**
> Pairwise geometric discrepancy for {n} RNA instances at {n} shared Rfam positions; lower means more similar. Rows and columns are ordered by tree-penalized path length seriation (TPPL), independently of structural groups.

*Ordering method* (collapsed details below the matrix, not a tooltip):
> Tree-penalized path length seriation (TPPL; Aliyev & Zirbel 2023), {repetitions} repetitions, seed {seed}. Similar structures are placed
> next to each other; the order uses the dissimilarity matrix only, not the structural groups.

**Heatmap cell, hover**
> {A} × {B}
> Geometric discrepancy {value} *(or "Geometric discrepancy: not measured")*
> {n} shared Rfam positions compared

**Heatmap diagonal cell, hover:** the instance details, i.e. PDB and chain, title, method, resolution and organism.

### Clustering

**Coverage, ⓘ "About clustering coverage"**
> {n} of {N} RNA instances were clustered. Inclusion requires at least {threshold} core positions with complete base, C4′ and C3′
> coordinates. Comparisons use {n} positions shared by all included instances.

**Metric, ⓘ "About the comparison metric"**
- *Geometric discrepancy:* "Geometric discrepancy compares nucleotide positions and orientations at shared Rfam positions. Lower values indicate greater similarity (0 = identical)."
- *TM-score:* "TM-score measures 3D similarity at shared Rfam positions (1 = identical). Clustering uses 1 − TM-score, where lower values indicate greater similarity."

**"provisional" badge, hover:** {clustering rationale}

**Pairwise heatmap, hover**
- Row label: "{seq_name} (group representative, medoid). Click to select this instance."
- Cell: "{A} × {B}: geometric discrepancy {value} ({n} positions)", or for TM-score "{A} × {B}: TM-score {value} (1 − TM {value}; {n}
  positions)".
- Diagonal: "{A}: itself".

**Dendrogram, hover**
- Floor line: "Duplicate floor {value}: differences below it are within the spread of duplicate structures"
- Cut line: "Cut at {value}"
- Group: "{group}: {n} RNA instances. Click to select."

**Structural groups card**
- ⓘ About structural groups: "Groups of RNA instances with similar 3D structure under the selected metric. They are descriptive and do not
  necessarily correspond to biological conformational states."
- "largest" badge: "Largest refers only to the number of RNA instances in the group; it does not imply biological importance."
- "small" or "single" badge: "Fewer than {n} RNA instances: possibly an outlier rather than a population; its statistics rest on few
  pairs."
- ⓘ About the representative: "Group member with the lowest mean dissimilarity to the others, selected among chains without backbone breaks."
- ⓘ About within-group dissimilarity: "Mean pairwise {metric} within this group; lower indicates greater similarity."
- ⓘ About across-group dissimilarity: "Mean pairwise {metric} between this group and all other clustered instances."

**Associated features**
- ⓘ About associated features: "Features are summarised across clustered RNA instances. These are descriptive comparisons, not
  statistical enrichment tests. RNA instances from the same PDB entry are not independent observations."
- Ligand associations ⓘ: "Ligand associations are derived from PDBe entity-level annotations and do not necessarily indicate direct
  contact with each RNA chain. A count of 0 means not annotated by PDBe, not confirmed absence."
- Species ⓘ: "Counts indicate RNA instances associated with each organism."
- Nonstandard nucleotides (Rfam-mapped region) ⓘ: "Counts indicate RNA instances containing the specified nonstandard nucleotide
  component within the Rfam-mapped region. Components are not necessarily post-transcriptional modifications."
- Nonstandard nucleotides (outside the mapped region) ⓘ: "Counts indicate RNA instances containing the specified nonstandard component
  outside the Rfam-mapped region (for example terminal nucleotides of the construct)."
- Experimental method ⓘ: "Counts indicate RNA instances determined using each experimental method."
- Continuous properties ⓘ: "Medians are calculated from RNA instances with available measurements."
- Count cell hover: "{percent} of {n} RNA instances". Median cell hover: "Median of {n} RNA instances with a value".

### 3D structure / 3D superposition

**ⓘ "About the 3D alignment"**
> Superposition: each chain is fitted onto {reference} using C4′ atoms at shared Rfam positions.
>
> Fit to the reference (C4′ RMSD):
> • {chain}: {quality} fit, RMSD {value} Å *(one line per displayed chain; "no alignment" when missing)*
>
> This visual superposition may differ from the alignment used to calculate the {geometric discrepancy | TM-score}. Poor fits are
> flagged on the structure chips.

**Other hover texts**
- Poor-fit warning on a structure chip: "The superposition onto the family frame is poor for this chain (C4′ RMSD > 4 Å) or unavailable;
  its placement may look offset."
- Structure chip: {seq_name}
- Compare button: "Superpose the two group representatives"
- Reset button: "Show the representative alone"

**Show ligands (checkbox in the viewer toolbar), hover** (Structures tab and Interactions Compare two)
> Show ligand instances within 4 Å of displayed RNAs as sticks. Excludes water, ions and crystallisation additives.
- While ligands are shown, the hover adds a line "Shown: {CCD ids}". The label adds "(loading…)" while they are fetched and
  "(none nearby)" when no displayed chain has one.

---

## Interactions

### Controls

**Residue numbering, ⓘ**
> Residue labels use the selected reference structure. Equivalent nucleotides across RNA structures are identified using Rfam consensus
> positions. Changing the numbering reference changes labels only.

**Interaction types dropdown:** each option is labelled with its type, e.g. Base pairing (Bp), with "no data" or "none observed" where it
applies; there are no tooltips.

**More filters:** an inline note, not a tooltip: "Type and context apply to base pairs only; other categories are hidden while they are
set."

**Search field, hover**
> G37, 37, G37-C61, Rfam 41, Rfam 41-66, cWW, s35, 9BR, J2/3 … Enter selects the first result

### Family conservation

**Interaction conservation (diagram), ⓘ**
> Interactions between corresponding Rfam positions. Arc colour shows the fraction of evaluable RNA instances with a recorded interaction.

**Arc, hover**
> {residue pair}{ · Multiple interactions}
> Rfam {i}–{j}
> {code} {class} · {k}/{n} ({percent}) · {Present | Absent | Not comparable} in {displayed chain} *(one line per interaction)*
> Any: {k}/{n} ({percent}) *(multiple interactions only)*

**Legend:** "Multiple interactions" has a hover, "More than one interaction type between the same two nucleotides". The legend note
(visible text) reads: "Each type has its own scale: colour intensity = observed / eligible chain instances (lighter = rarer, darker = more
conserved). Multiple interactions use the 'any' occurrence."

**Table header icons**
- Interaction class ⓘ:
  > Interaction classifications follow BGSU/FR3D annotations. Codes describe the interaction geometry or contact type.
  > Bp: Leontis–Westhof family (cis/trans, edges W/H/S).
  > Bs: stacking face orientation (s35, s53, s33, s55).
  > Br, Bph: base–ribose / base–phosphate class 0–9 (position of the base donor).
  > n prefix: near interaction (not counted).
- Conservation ⓘ: "Frequency is the fraction of evaluable RNA instances with a recorded interaction."

**Table cell hovers**
- Interaction: "Rfam {i}–{j} · {region}" (→ instead of – for base–ribose and base–phosphate).
- Type: the full type name; for a group, "Multiple interactions: {types}".
- Interaction class: the class distribution, e.g. "cWW 76 · tHS 4".
- In displayed structure:
  - *Present*: "Annotated: {classes}"
  - *Absent*: "Resolved and eligible; not annotated." followed by " Other: {other annotations}" when there are any
  - *Not comparable*: the reason, one of:
    - "Rfam {pos}: nucleotide not modelled | incomplete base | incomplete ribose | incomplete or absent phosphate | no residue at this
      Rfam position (deletion or outside the aligned region)"
    - "FR3D annotation unavailable ({note})"
    - "{type} annotations unavailable for this chain"
    - "BGSU annotation incomplete for this chain (no annotation at Rfam {a}–{b})"

**Details panel**
- Multiple interactions ⓘ: "More than one interaction annotation between the same two nucleotides. Each keeps its own class and
  frequency; 'any' counts chains with at least one of them among chains where at least one can be evaluated."
- "Not evaluated" row, hover: "Chains where this interaction cannot be evaluated; never counted as absent"

**Displayed structure, ⓘ**
> The structure shown in 3D. The "In displayed structure" column refers to it; it does not change the numbering or the family-wide
> frequencies.

### Compare two

**Annotation coverage, ⓘ**
> Of the Rfam consensus (SS_cons) helix base pairs that can be evaluated in the chain (both bases resolved and complete), how many
> BGSU/FR3D annotates. Distorted or moderate-resolution models can have paired nucleotides that FR3D does not annotate; BGSU exports no
> near base pairs.

**Sequence relation, ⓘ**
> Sequence relationship across the mapped chains: identical; same RNA at shared Rfam positions (allowing limited terminal construct differences); near-identical (≤1% difference); or different. Exact classification rules are in the methods.

*Sequence classification rules* (collapsed details below the sequence note, not a tooltip):
> Identical: the same chain sequence. Same RNA at shared Rfam positions: different deposited chain sequences that agree on every shared
> Rfam position, or differ only at up to 3 positions within 5 positions of a chain end (e.g. different constructs). Near-identical: at most
> 1% of shared positions differ (e.g. rRNA operon variants). Different: other sequence differences.

**Status chips and legend, hover**

| Status | Text |
|---|---|
| Shared | Same interaction type and FR3D class in both structures |
| Only in A | Annotated in A; both nucleotides eligible in B but not annotated there |
| Only in B | Annotated in B; both nucleotides eligible in A but not annotated there |
| Different class | Same interaction type, different FR3D class (e.g. cWW in A, tWW in B) |
| Not comparable | Cannot be evaluated in one structure (missing, incomplete or unmapped nucleotide, or no annotation) |

**Arc, hover**
> {residue pair}{ · Multiple interactions}
> Rfam {i}–{j}
> {code} **{status}** · A: {classes | state} · B: {classes | state} *(one line per interaction)*

**Comparison table hovers**
- Interaction: "Rfam {i}–{j}".
- "—" or "n/c" cells: the not-comparable reason; otherwise "Other: {annotations}" or "not annotated".
- Class in A ⓘ: the same as Interaction class above.

**Local superposition, ⓘ "Local superposition details"**
> Structure B is aligned to A using {n} corresponding base heavy atoms of {residue} (Rfam {pos}); RMSD {value} Å. The partner nucleotide is not fitted, preserving its relative position.

**Diagram titles** (accessible labels): "Family-wide interaction conservation diagram"; "Interactions of {chain}".

**Pager buttons:** "First page", "Last page".

---

## Ligands

The Ligands tab uses the Annotations track browser (axis, numbering, zoom, legend, track hover and pinned details as listed under
Annotations below). Texts specific to this tab:

**RNA-contacting ligands ⓘ (gallery heading)**
> Ligands shown here have contacts reported by PDBe for RNA entities belonging to this Rfam family. Contacts are derived from Arpeggio
> and may not represent every RNA chain or biological assembly in an entry.

**Gallery description:** "Ligands with PDBe-reported contacts to RNA chains of this family, with the number of PDB entries in which a contact
is recorded. Metal ions, crystallisation additives and flagged artefacts are not shown. Select a ligand to locate its residue track."

**Gallery card:** CCD id · name · "Found in {n} PDB entr(y|ies)" (hover: "RNA contact recorded by PDBe in: {entries}") · "PDBe ligand page ↗".
Depiction hover/alt: "2D structure of {CCD}"; when no image loads: "No depiction available". Image hover: "Locate the {CCD} residue track"
(selected: "Clear the highlighted track"). Search placeholder: "Name or CCD id, e.g. adenine, HPA"; no match: "No ligand matches “{query}”."
Empty family: "No ligands with PDBe-reported RNA contacts in the default ligand set for {ACC}."

**Category headers, ⓘ**

| Category | Text |
|---|---|
| Ligand-contact residues | Ligand-contact residues are obtained from PDBe binding-site annotations and mapped to Rfam consensus positions. The annotations are aggregated at RNA entity level rather than individual chain level. |
| Predicted binding pockets | Computationally predicted RNA binding pockets identified using fpocketR. These annotations are independent of the PDBe ligand-contact data. |

**Legend ⓘ (frequency)**
> Contact frequency is the proportion of eligible PDB entries containing a ligand in which PDBe reports a contact at this position. Entries
> with unobserved residues are excluded from the denominator. Pocket colours show the proportion of RNA instances with a predicted pocket.

**Position numbering ⓘ**
> Residue numbering is taken from the selected RNA chain. Ligand-contact annotations remain aggregated at RNA entity level and are not
> necessarily specific to the selected chain. Switching numbering systems changes residue labels only; the underlying ligand contacts and
> frequencies remain unchanged.

**About these annotations (collapsed section):** Source · Contact identification · Annotation level · Coverage · Coordinate mapping ·
Alternative numbering · Frequency unit · Limitation · Predicted pockets (one sentence each; see `frontend/src/components/LigandsTab.tsx`).

**Family Overview note after "View in Family Overview":** "Opened from Ligands ({label}). Family-level evidence: not necessarily observed in {chain}."

**Download menu:** ligands.json, tracks.json, pockets.json, pocket_tracks.json (texts as in the Annotations table below).

---

## Annotations

The Annotations tab shows Solvent accessibility, RNA loop annotations and Nonstandard nucleotides; its summary adds "Ligand contacts and
predicted binding pockets are in the Ligands tab." Track hover, pinned-detail and download texts below also apply to the Ligands tab.

**Category headers, ⓘ "About {category}"**

| Category | Text |
|---|---|
| Solvent accessibility | Per-nucleotide solvent-accessible surface area (SASA) of isolated RNA chains, summarised by Rfam position (median Å²). |
| RNA loop annotations | BGSU hairpin (HL), internal (IL) and junction (JL) loops mapped to Rfam positions. Annotations may overlap or vary between structures. |
| Nonstandard nucleotides | Nonstandard nucleotide components in deposited RNA sequences, mapped to Rfam positions. These are not necessarily biological RNA modifications. |

**Legend swatches, hover**
- No data: "No assessable RNA instance at this position: residues or source data unavailable (SASA: no qualifying measurement)."
- Not observed: "No annotation recorded at this position. Interpretation depends on the annotation type; for loops and ligands, non-detection does not establish absence."

**Track label, hover:** {track label} — {sub-label}

**Track cell, hover** (title "Rfam residue {pos}"). With a numbering reference, each line is followed by "{reference} reference residue:
{residue}".
- Loops: "{loop type} ({HL | IL | JL})" · "Observed in {k} of {n} RNA instances ({percent})"

**Track cell, click (pinned details), loops:** loop observations, instances and entries; "Not annotated as {HL | IL | JL} by the source in {k} of
{n} RNA instances; source unavailable for {m}"; per-type junction counts (JL); motif groups; "{n} loop(s) merged by BGSU into a larger loop not
counted here." *(only where such loops exist)*; "Frequencies are calculated from instances where annotation presence or absence can be
determined. BGSU loop classifications can overlap."
- SASA: "Median SASA: {value} Å²" · "Based on {n} RNA instance(s)"
- Pockets: "Pocket occurrence: {k} of {n} RNA instances ({percent})"
- Nonstandard nucleotide: "{CCD} — {name}" · "Observed in {k} of {n} RNA instances ({percent})"
- Ligand: "{CCD} — {name}" · "Recorded contact in {k} of {n} observable PDB entries ({percent})"

**"View in Family Overview" button, hover**
> Select this Rfam position for the RNA instance shown in the Family Overview

**Download menu, ⓘ per file**

| File | Text |
|---|---|
| RNA loop annotations (loops.json) | BGSU RNA 3D Motif Atlas loops (HL, IL, J3, J4) with motif-group assignments, nucleotide and segment membership and Rfam consensus mapping. |
| Ligand binding annotations (ligands.json) | Bound molecules (CCD id, name, class), binding-site nucleotides, supporting entries and chains, and provenance. |
| Aggregated annotation tracks (tracks.json) | Per-position annotation frequencies with numerators, denominators and missing-data counts. |
| Predicted pockets · Per structure (pockets.json) | fpocketR pockets per RNA chain with scores, volumes, descriptors and lining nucleotides mapped to Rfam positions. |
| Predicted pockets · Aggregated (pocket_tracks.json) | Per-position pocket frequencies for each of the three pocket filters. |
| Solvent accessibility · Per nucleotide (sasa.json) | Per-chain nucleotide solvent accessibility (Å²) with atom completeness and calculation details. |
| Solvent accessibility · Aggregated (sasa_tracks.json) | Per-position solvent accessibility summaries (median, quartiles) and contributing-chain counts. |
| Nonstandard nucleotides · Detailed (modifications.json) | Nonstandard nucleotide identities, locations, parent nucleotides, Rfam mappings and source information. |
| Nonstandard nucleotides · Aggregated (modification_tracks.json) | Per-position nonstandard-nucleotide frequencies, counts and contributing components. |

The download link hover shows the file name, "{ACC}_{file}".
