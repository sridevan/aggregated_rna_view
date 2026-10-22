# Investigation: the existing RNA–ligand annotation pipeline

**Scope:** an investigation of 2026-10-22, at commit e1e93c2. No application code, data, pipeline or UI was changed.

**Evidence labels:**
- **[V] verified here:** from code read in this session, generated data, or the coordinate check in §6.
- **[P] verified in an earlier investigation:** `experiments/chain_specific_ligands/`, which ran live API calls on 1Y26, 3RKF and 9NDD on 2026-10-15.
- **[D] documented by the provider, not verified here.**

**Coordinate check:** `scripts/inventory.py` (read-only) lists every non-water non-polymer copy in model 1 of the cached PDBe updated mmCIF. For each copy it computes the minimum heavy-atom distance to each family RNA chain and the RNA residues within 4 Å and 5 Å. Outputs are in `data/<ACC>_inventory.json`.
- Run for RF00167, RF00162 and RF00059.
- Not run for tRNA, 16S or RNase P: their mmCIFs are not cached, and nothing was downloaded.

## 1. Executive summary

`ligands.json` and the Annotations-tab "Ligand binding sites" tracks are **PDBe binding-site annotations at the level of the RNA *entity***.
- One record is one (PDB entry, RNA entity, CCD).
- Its residues are the union of the contacting residues of the ligand copies PDBe processed. These are the copies in PDBe's preferred assembly, with contacts derived from Arpeggio interactions [P].

**What the data does not contain:**
- **Ligand instances.** There are no ligand chain or residue identifiers, only PDBe bound-molecule ids such as `bm1`.
- **Distances.**
- **Chain-level contact attribution when an entity has several chains.**
- **Ligands that do not contact the RNA entity.**

The tracks give, per CCD and Rfam position, the **fraction of PDB entries containing that CCD in which PDBe records a contact at that position**. The denominator only includes entries that have an observed nucleotide there. It is an entry-level observed-contact frequency, not a chain-level or instance-level one, and not a binding probability.

**For RF00167 the data are complete at the (entry, CCD) level.**
- Every non-water component present in the coordinates (166 entry/CCD pairs) has a record.
- Every copy lies within 5 Å of a family RNA chain.

**At the instance level the data are incomplete.**
- In multi-chain entries PDBe processed one ligand copy per entry, so 11 ADE copies are represented by 9 bound molecules.
- `bestChainId` can be wrong: in 5SWD it names chain B, but the only ADE copy is on chain A.

A separate mechanism exists:
- Mol*'s **"Show ligands"** computes **chain-specific 4 Å proximity in the browser** from the loaded coordinates.
- It is a display rule only and feeds no statistics.

## 2. Data provenance

| Step | Location | What it does |
|---|---|---|
| Fetch | `rna_explorer/annotation_sources.py:249-264` `fetch_binding_sites` | GET `https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/binding_sites/{pdb}/{entity_id}` (base URL at l.30), cached under `cache/www_ebi_ac_uk/`. It is seeded from `annotation_assessment/cache` with `--raw-dir`. HTTP errors become `status: unavailable`. **[V]** |
| Build records | `rna_explorer/annotations.py:339-436` `build_ligands` | Called once per (entry, RNA entity) of the family (`idx.by_entity`). Each `site` becomes one record. Residues `startIndex..endIndex` (label_seq_id) are deduplicated per residue and keep their `boundMoleculeId`s. **[V]** |
| Classify | `annotations.py:331` `_classify`; `rna_explorer/annotation_config.json` | Matches CCD ids exactly against the `water`, `metal_ion` and `additive` lists; anything else is `ligand`. Two artefact rules apply (§3). **[V]** |
| Map to Rfam | `annotations.py:361, 375-381`; `annotation_mapping.py:147` `resolve_seq_index` | Each label_seq_id is resolved on one chain of the entity (`map_inst`, the `bestChainId` chain if it belongs to the entity, else the first chain), giving `rfam_consensus_pos`, `mapping_status` and PP. `chains_with_coordinates` lists the entity's chains in which that residue is modelled. **It is not contact evidence.** **[V]** |
| Tracks | `annotations.py:553-588` (in `build_tracks`) | Per position: entity-level totals and per-CCD `by_ccd_entries` (§5). **[V]** |
| Outputs | `data/families/<ACC>/ligands.json`, `tracks.json` `positions[].ligand_binding_sites`, `annotations/<seq_name>.json` `residues[].bound_molecules` | Static JSON; no runtime API calls for these tracks. **[V]** |
| CLI | `python3 -m rna_explorer annotate <ACC>` (`cli.py:22-25, 96-99`) | **[V]** |

**Consumers [V]:**

| Consumer | Location | What it reads |
|---|---|---|
| Annotations tab | `frontend/src/components/AnnotationsTab.tsx:208-243` `buildTracks` | One track per CCD with `by_ccd[ccd].n_entries_in_default_track > 0`, sorted by entry count. Colour and tooltip come from `tracks.json` `by_ccd_entries`; supporting entries from `ligands.json` `bound_molecules`. |
| Category header and downloads | `AnnotationsTab.tsx:76, 395, 436` | |
| Structures tab, "Ligand associations" | `frontend/src/components/structures/AssociatedFeatures.tsx:56-79` | Entity-level records, assigned to every chain of the entity; the tooltip says so. |
| Mol* "Show ligands" | `PdbeMolstarViewer.tsx:65, 164-169`; `molstar/MolstarSuperposition.tsx:30, 313-327`; `lib/molstarScene.ts:45-52`; `api/data.ts:61` | Independent of `ligands.json` contacts. It only uses `by_ccd[].ligand_class` to hide metals and additives. |

## 3. Ligand identification rules

| Question | Answer |
|---|---|
| 1. All non-polymers of the entry retrieved? | **No.** Only molecules PDBe lists as binding sites **of the RNA entity**. A molecule that contacts only a protein or another RNA does not appear. [V, by construction of the endpoint call] |
| 2. Only molecules close to the RNA? | Yes: binding-site partners of the entity, computed by PDBe. [V for the outputs; mechanism P/D] |
| 3. Distance cutoff? | **The pipeline applies none.** PDBe's residues are Arpeggio interaction partners, which use per-interaction geometric criteria, not one cutoff. For single-copy 1Y26 the v2 residues equalled the Arpeggio partner set. [P] RCSB `rcsb_ligand_neighbors` uses 5 Å (empirical maximum 4.85 Å) but is **not used** in production. [P] Mol* "Show ligands" uses 4.0 Å (`PdbeMolstarViewer.tsx:65`). [V] |
| 4. Any atom to any atom? | Not applicable to production (typed Arpeggio contacts). Mol* "Show ligands": any atom within 4 Å (`filter.within`, `max-radius`). [V] |
| 5. Hydrogens? | PDBe's pipeline places hydrogens before Arpeggio. [P/D, not re-verified] The Mol* rule uses the model as loaded. This report's check uses heavy atoms only. |
| 6. Per chain or per entry? | **Per entity.** Residues from all processed copies are pooled on the entity sequence (`observation_level: "entity"`). `chain_attribution` is `resolved` only when the entity has one chain (`annotations.py:392`). [V] |
| 7. All copies? | **No.** Only copies in PDBe's preferred assembly are processed [P]. RF00167 evidence is in §6. |
| 8. Symmetry? | PDBe computes on the preferred assembly, which may include symmetry-generated copies; no record exposes the operator [P]. The pipeline does nothing with symmetry. [V] |
| 9. Assembly or asymmetric unit? | PDBe: preferred assembly [P]. Mol* "Show ligands": ModelServer `residueSurroundings` on `model_nums=1` of the deposited (asymmetric-unit) model (`molstarScene.ts:46-47`). [V] |

**Chemical interactions or proximity?**
- PDBe's residue lists derive from typed interactions (hydrogen bonds, stacking, van der Waals, …), but **the types are discarded**. `ligands.json` keeps only "residue i is a binding-site residue of CCD X in entity E". [V]
- "Show ligands" is pure proximity. [V]

**Classification (`annotation_config.json` v1.1) [V]:**
- **Categories by CCD id:**
  - water: 4 ids;
  - metal_ion: 38 ids;
  - additive: 54 ids (crystallisation agents, buffers, polyamines SPD, SPM, PUT; SO4, PO4 and CL also count as additives);
  - everything else: `ligand`.
- **Default track** = class `ligand` and no excluding flag.
- **Flag `ccd_is_polymer_residue`:** the CCD is also a residue name of the polymer.
- **Flag `contacts_only_polymer_termini`:** all contacts lie within 2 residues of a terminus.

**What each type of molecule becomes:**

| Molecule | How it is handled |
|---|---|
| Waters | PDBe does not report them as bound molecules; the class exists for safety. |
| Modified nucleotides | Polymer residues, never bound molecules here; they have their own track. |
| Covalently attached groups | Reported by PDBe when they are bound molecules. Terminal ones are flagged and excluded from the default track. |
| Branched (saccharide) entities | Treated by the pipeline like any other CCD. |

**Possible unintended exclusions [V]:**
- **The terminus rule (tRNA).** 16 tRNA records are flagged and excluded:
  - aminoacyl groups (VAL, PHE, MET, CYS, HIS) attached to the 3′ end, which is a reasonable exclusion;
  - **non-covalent synthetase substrate analogues bound at A76**: aminoacyl-adenylate or sulfamoyl analogues QSI, A5A, GAP, AMO, GOM, plus AMP and CDP.
  - These are biologically meaningful ligands that the rule hides from the default track.
- **The terminus rule (RF00167):** it excludes GTP (4 entries), the 5′-terminal transcription nucleotide. This is intended.
- **Ions and additives:** SO4, PO4 and polyamines can be functional (for example, polyamines in ribosomes). They are kept in `ligands.json` but never shown as tracks.
- **Missing PDBe records (TPP):** three ions present within 5 Å have no PDBe record (2CKY OS ×12, 2HOP CA, 7TZR NA). That is PDBe coverage, not our filter.

**Ribosomes (RF00177) [V]:**
- 107 records in 52 evaluated entries (6 entity sources unavailable).
- The `ligand` class holds 22 antibiotics and inhibitors; MG has 2,870 bound-molecule instances.
- The class scheme works, but:
  - (a) contacts with the 16S entity only are recorded, which is correct;
  - (b) 6 records are entity-unresolved because two ribosomes are in the asymmetric unit;
  - (c) polyamines, K and MG are hidden. That is a scientific choice for ribosomes, where they are structural.

## 4. Contact mapping: what is preserved

| Item | In `ligands.json`? |
|---|---|
| PDB entry id | yes (`pdb_id`) |
| RNA chain | **only for single-chain entities.** Otherwise `entity_chains` plus `best_chain_label_asym_id`, which is metadata and wrong in 5SWD and 3RKF [V, P]. |
| CCD id, name | yes (`ccd_id`, `name` from PDBe) |
| Ligand instance id | PDBe `bound_molecule_ids` only (bm1…). There is no ligand `label_asym_id` or author chain and number. |
| Ligand chain and residue | **no** |
| Contacting RNA residues | yes, as label_seq_id on the entity (`contacts[].label_seq_id`, `residue_name`) |
| Distance | **no** |
| Rfam position | yes (`contacts[].rfam_consensus_pos`, `mapping_status`, PP), resolved on `map_inst` |

**How contacts reach Rfam positions.** Before mapping, contacts are not assigned to a chain or ligand instance:
- Entity residues are mapped through one representative chain.
- Every chain of the entity has the same sequence, and the same consensus mapping in all RF00167 cases checked.

**Aggregation across chains** counts **entries**, not chains, so an entry with several chains or copies contributes once (§5).

## 5. Occurrence frequencies (`annotations.py:553-588`, displayed at `AnnotationsTab.tsx:219-241`)

For CCD *c* at position *p* (`by_ccd_entries[c]`):
- `E_c` = PDB entries with a default-track record for *c*.
- **eligible** = entries in `E_c` with at least one chain that belongs to a ligand-source-evaluated entity **and** has a modelled match residue at *p* (`idx.observable`, `annotation_mapping.py:109`).
- **positive** = entries in `E_c` whose record for *c* contains *p*.
- **frequency = |positive ∩ eligible| / |eligible|**.
- `n_not_observable` = entries in `E_c` with no chain observed at *p*.
- `n_source_unavailable` = entries in `E_c` minus eligible minus not observable.

**Answers [V]:**

| Question | Answer |
|---|---|
| Unit | PDB entry |
| Several copies or chains in one entry | Counted once. An RNA instance cannot contribute twice. |
| Restricted to entries containing *c*? | Yes. |
| Missing residues | Excluded from the denominator, never counted as negatives. |
| Entry with *c* but no contact at *p* | Counted as "no contact recorded by PDBe", not as an absence. |

`ligand_binding_sites.frequency` is a separate value: it is all-ligand and counted per entity. It is only shown in the pinned-details text.

**Worked example: HPA, 12 entries** (2EES, 2EET, 2EEU, 2EEV, 2EEW, 4FE5, 4FEJ, 4FEL, 4FEN, 5C7U, 5C7W, 9LKE), all single-chain entities. Recomputed from `ligands.json` and the mappings, matching `tracks.json` exactly:

| Position | Positive | Eligible | Frequency |
|---|---|---|---|
| 25 | 12 | 12 | 1.0 |
| 51 | 12 | 12 | 1.0 |
| 79 | 3 (4FE5, 4FEL, 4FEN) | 12 | 0.25 |
| 54 | none | — | no `by_ccd_entries` entry, so no mark ("Not observed" style) |

**ADE, 9 entries:**
- Positions 52 and 53: 1 of 9. The only entry is 5SWD, whose single ADE sits in an unusual pose on chain A (minimum distance 0.47 Å to the RNA, i.e. overlapping alternate conformations; not investigated further).
- Position 79: 2 of 9.

## 6. RF00167 evidence

The table compares each CCD's coordinates with what PDBe recorded:

| CCD | Entries | Copies in coordinates | Copies ≤ 4 Å of a family chain | PDBe bound molecules | Multi-chain entities |
|---|---|---|---|---|---|
| HPA | 12 | 12 | 12 | 12 | 0 |
| ADE | 9 | **11** | 11 | **9** | 3 (5SWD, 6VWT, 6VWV) |
| GUN | 4 | **5** | 5 | **4** | 1 (9LMC) |
| 6AP | 3 | 3 | 3 | 3 | 0 |
| 6GU | 3 | 3 | 3 | 3 | 0 |
| GNG | 3 | **4** | 4 | **3** | 1 (3DS7) |
| OXG | 3 | 6 | 6 | 6 | 1 (9NDD, all 4 processed) |

**Per-copy check** (`data/RF00167_inventory.json`; chain-level residues ≤ 4 Å mapped through each chain's own mapping):
- **Unprocessed copies.** In 6VWT, 6VWV, 9LMC and 3DS7 the chain-B copy has the same pocket as chain A (consensus 25, 26, 51, 55, 56, 79, 80, 81), but PDBe processed only `bm1`. The entity record therefore cannot say which chain it describes.
- **The missing copies are exactly the copies outside one assembly.** PDBe annotates interactions only in the preferred assembly. `_pdbx_struct_assembly_gen` in the cached mmCIFs [V] shows:

  | Entry | Assembly 1 | Assembly 2 | Unprocessed copy |
  |---|---|---|---|
  | 6VWT | A + ADE C | B + ADE H | H |
  | 6VWV | A + ADE C | B + ADE F | F |
  | 9LMC | A + GUN C | B + GUN J | J |
  | 3DS7 | A + GNG I | B + GNG Q | Q |
  | 9NDD | one assembly with all four chains | — | none (4 of 4 processed) |

  - In each two-assembly entry, the unprocessed copy is the one in assembly 2.
  - The mmCIF has no "preferred" flag. The data imply that assembly 1 is preferred in these four entries; confirming this needs PDBe `api/pdb/entry/summary`.
  - This is a systematic, predictable gap: in any entry whose asymmetric unit holds several author-defined assemblies, only one assembly's ligand copies, and only one chain's contacts, are annotated.
- **Wrong `bestChainId`.** 5SWD's `bestChainId` is **B**, but its only ADE copy is on **chain A**, in assembly 1 (A, C, D, E, G) and within 5 Å of A only. The field is wrong whichever assembly is preferred.
- **Consensus 79.** Proximity includes it for almost every copy (stacking or packing with A/C74), while PDBe's typed contacts include it in only a minority of entries (HPA 3/12, ADE 2/9). This confirms the earlier finding that proximity sets are supersets by 1–2 residues. [P]
- **No mismatches.** For all 166 records (67 ligand-class, 99 ion or additive), PDBe's entity-level positions are a subset of the union of the ≤ 5 Å proximity sets of that CCD's copies in the entry.
  - Ion records pool many copies, so the comparison is against the union, not a single copy.
  - Chains of the same entity have identical consensus mappings in every RF00167 multi-chain entity, so mapping through one representative chain loses no position.

**UI consistency [V]:**
- 29 default-track CCDs, matching the 29 tracks.
- Track sub-labels use `n_entries_in_default_track`, e.g. HPA "12 entries".
- Tooltips reproduce `by_ccd_entries` exactly.
- No numerical discrepancy was found between the data and the display.
- The terminology discrepancy: "Ligand binding sites" and "in contact" are presented without the entity-level and preferred-assembly caveat. Only the category ⓘ mentions "Entity-level annotations may not apply to every RNA chain".

## 7. Are non-contacting ligands available?

1. **Only contacting ligands are in the dataset.** `ligands.json` contains only ligands that PDBe reports as binding sites of the RNA entity. [V]
2. **Non-contacting ligands can't be identified from current data.**
   - In RF00167, RF00162 and RF00059 this barely matters: every non-water copy is within 5 Å of a family chain, except 2 ions each in SAM-I and TPP.
   - In ribosomes, most entry ligands do not contact 16S, and they are absent.
3. **Sources for "all ligands in the entry"** (none used now):
   - PDBe `api/pdb/entry/ligand_monomers/{pdb}`: every instance, label and author ids, `bm_id` null if unprocessed [P];
   - RCSB `core/entry` and `nonpolymer_entity_instance` [P];
   - or the mmCIF itself, as `scripts/inventory.py` does.
4. **Cost.**
   - Low for small RNAs: one call or file per entry, ≈ 3 s for all 66 RF00167 entries from cached mmCIF.
   - Ribosomes need large mmCIFs or API calls, and produce hundreds of ions per entry.
5. **Ambiguity.** In multi-chain or RNA–protein entries, "present in the entry" says nothing about the selected RNA (e.g. a 70S antibiotic bound to 23S).
   - An "all ligands" filter would need a per-ligand relation: contacts this RNA chain / another chain / no polymer.

**For the first version, use "RNA-contacting ligands only"**, with the entity-level caveat. It is what the current data mean.

**Open question [V, cause unresolved]:**
- PDBe binding_sites returned 404 for 27 of 60 tRNA entities (e.g. 1ASY, 7K00 entity 54) and 6 of 60 16S entities.
- The pipeline records them as "source unavailable".
- Whether a 404 means "no processed bound molecule" or "not processed" was not determined. It decides whether those entities are true negatives.

## 8. Terminology

The data are **computationally identified contacts** (PDBe Arpeggio), pooled per entity, from processed copies. They are not experimentally characterised binding, and not plain proximity (production) — except "Show ligands", which is proximity.

| Term | Use |
|---|---|
| Tab name | "Ligands" is fine for the tab. |
| Ligand set | **"RNA-contacting ligands"**, with the qualifier "(PDBe)". |
| Residue tracks | **"Ligand-contact residues"** or keep "Ligand binding sites" with the method line "PDBe-computed contacts (Arpeggio), entity level". |
| Tooltips | Keep "recorded contact". |
| "Directly interacting ligands" | Avoid; the types are discarded and it implies validated chemistry. |
| "RNA–ligand proximity sites" | Only for the Mol* 4 Å rule. |

## 9. Ligands-tab readiness (no change to calculations)

| Feature | Status |
|---|---|
| Ligand gallery | **Ready.** CCD ids, PDBe names and class are in `by_ccd`. |
| 2D chemical diagrams | **Possible from CCD id.** PDBe and RCSB serve CCD images by id [D, URL patterns not verified here]. No local data is needed. |
| PDB entry counts | **Ready and correct** (`by_ccd.n_entries`, `n_entries_in_default_track`). Instance counts are **not** reliable (`n_bound_molecule_instances` = processed copies only). |
| RNA-contacting filter | **Default meaning** of the current data. |
| All-ligands filter | **Not available.** It needs a new source (§7). |
| Ligand-binding residue tracks | **Reusable as is** (`buildTracks` per CCD, `by_ccd_entries`). |
| Instance-specific numbering | **Reusable** (existing correspondence and mapping). For multi-chain entities, chain-specific highlighting is **not** supported by the data: residues are entity-level. |
| Predicted pockets (fpocketR) | Movable unchanged (`pocket_tracks.json` is independent of `ligands.json`). |

**Limitations to address before or during implementation:**
- **(a) Chain and instance level:** no chain-level contacts for multi-chain entities, and unprocessed copies are invisible.
- **(b) `bestChainId` is unreliable:** do not use it to pick a chain.
- **(c) No distances or interaction types.**
- **(d) The terminus flag hides synthetase substrate analogues** in tRNA.
- **(e) PDBe 404 semantics** are unresolved.
- **(f) The Mol* 4 Å rule and the PDBe tracks can disagree** on the same screen. For example, consensus 79 is near the ligand but absent from most PDBe records.

## 10. Recommendations (minimal)

1. **Ship version 1 as "RNA-contacting ligands (PDBe)"** using `ligands.json` and the existing per-CCD tracks. Label counts as **PDB entries**, never instances.
2. **State the method in the tab:**
   - PDBe-computed contacts;
   - entity level;
   - preferred-assembly copies only;
   - not proof of binding.
   - Use the existing 4 Å proximity wording only for "Show ligands".
3. **Do not display `n_bound_molecule_instances` or `bestChainId`** as facts.
4. **If chain-specific views or instance counts are wanted,** implement the chain-level layer already designed in `experiments/chain_specific_ligands/RECOMMENDATIONS.md`:
   - RCSB 5 Å neighbours, or local proximity as in `scripts/inventory.py`, per ligand instance;
   - joined to PDBe `ligand_monomers`;
   - as a separate counting unit.
5. **Before an "all ligands" filter:** add an entry ligand inventory with a per-instance relation to the selected chain. Defer it for ribosomes.
6. **Review the `contacts_only_polymer_termini` rule** for tRNA synthetase complexes, e.g. exempt non-covalent nucleotide analogues.
7. **Resolve the meaning of PDBe 404s** for tRNA and 16S entities before treating them as unavailable or negative.

## Files

- `scripts/inventory.py`: read-only coordinate inventory (gemmi 0.6.7, cached PDBe updated mmCIF, model 1, heavy atoms).
- `data/RF00167_inventory.json`, `data/RF00162_inventory.json`, `data/RF00059_inventory.json`: per-copy distances and 4 Å / 5 Å residues.
