# Modified (nonstandard) nucleotide components in RF00167 — investigation report

Family: RF00167 (purine riboswitch), Rfam 15.1, CM checksum 1292453307, CLEN 102; 66 PDB entries, 84 mapped RNA chain instances.
Scope: investigation only. No production extraction, aggregation, JSON outputs, downloads or frontend category were created. No existing data or code outside this directory was modified.
Date of investigation: 2026-10-16. Pipeline mappings used as-is: `data/families/RF00167/mappings/*.json`, `residue_mapping.csv` (pipeline 0.1.0, generated 2026-10-15T03:48:01Z).

Labels used below: **[verified]** = checked directly against source data in this directory; **[inference]** = interpretation not proven by the data.

---

## 1. Executive summary

- **The data are small.** 13 nonstandard polymer residues occur in 9 of 84 RNA instances from 8 of 66 PDB entries, using 7 distinct CCD components. All other 6,588 deposited polymer residues are A, C, G or U. [verified]
- **Only 6 of the 13 map to Rfam consensus match positions**, at 5 distinct positions: 22, 32, 35, 84 and 87. The other 7 are 5′ or 3′ terminal residues outside the PDBe/Rfam mapped region, so they have no consensus position. None falls in an insert state. [verified]
- **The deposited mmCIF finds every observation; PDBe misses two.** The PDBe `modified_AA_or_NA` endpoint reports 11 of the 13. It omits 50N and 50L in 5C7W. Where both sources report a residue, they agree exactly on entity, `label_asym_id`, `label_seq_id` and CCD. [verified]
- **PDBe cannot establish negatives.** An entry with no modified residues and a nonexistent entry return the same HTTP 404 response. [verified]
- **NAKB/BGSU covers all 7 components.** It provides a parent standard nucleotide and atom-level changes by location: base, ribose or phosphate. It is chemistry reference data, not per-chain evidence. [verified]
- **23G in 3IVN** is a polymer residue: entity 1, label_seq_id 70, author 74, the 3′-terminal residue of chains A and B. It has no coordinates and lies outside the mapped region, so it has **no Rfam position**. PDBe, mmCIF and NAKB agree on its identity. [verified]
- **Recommendation.** Use the deposited mmCIF polymer sequence as the primary source, with identity read from `_pdbx_poly_seq_scheme` at `label_seq_id`. Use PDBe only for cross-validation and NAKB only for optional enrichment. The data support the user's initial preference.
- **Provisional occurrence.** Frequencies are very low: 1 of 83 to 2 of 84 instances at four positions, and 2 of 77 at position 87. With the existing five-bin frequency scale, every mark would fall in the lightest bin. One track is still appropriate, but the colour encoding needs a decision (§11).

## 2. PDBe modifications API assessment

Endpoint: `POST https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/modified_AA_or_NA`, body = JSON string of comma-separated PDB ids. A `GET .../modified_AA_or_NA/<pdb>` form returns the same records. Retrieved 2026-10-16 (UTC 17:27 onwards). Raw responses: `raw/pdbe/`.

**Record structure** [verified on 3IVN, 1Y27, 5UZA, 7K00, 1EHZ, 122D]. Records are grouped under the lower-case PDB id key, one record per residue:

| Field | Meaning (verified against mmCIF) |
|---|---|
| `entity_id` | polymer entity id |
| `struct_asym_id` | `label_asym_id` |
| `residue_number` | `label_seq_id` (11/11 RF00167 records) |
| `chain_id` | `auth_asym_id` |
| `author_residue_number`, `author_insertion_code` | `auth_seq_id`, insertion code |
| `chem_comp_id`, `chem_comp_name` | CCD id and CCD name |
| `alternate_conformers` | always 0 in all records seen |
| `description` | `_entity.pdbx_description` |
| `weight` | `_entity.formula_weight` (3IVN 22574.332, 1Y27 21994.916, 5UZA 22939.178) |

**Findings**

- **Chain-specific, not entity-level.** 3IVN gives one record per chain: A/70 and B/70. [verified]
- **RNA and DNA are both covered.** 122D, a DNA duplex, returns 5CM. [verified]
- **Protein modifications also appear.** 7K00 returns IAS, D2T, 4D4 and MS6 on ribosomal proteins. Records carry **no polymer-type field**, so polymer type must come from the mmCIF entity. [verified]
- **Polymer residues only.** No RF00167 record referred to a non-polymer ligand. GTP is reported for 5′ polymer residues but not for the non-polymer GTP ligands in 9IWF, 9IWG, 9LKE and 9LMC. [verified for RF00167]
- **Empty and unavailable look the same.** `"3rkf"`, `"1y26"`, the invalid `"0zzz"` and `"3rkf,0zzz"` all return HTTP 404 with `{"message":"Requested endpoint does not contain any data"}`. Entries without records are simply omitted from a batch response. [verified]
- **Batch equals per-entry.** A batch of all 66 entries returned 7 entries and 11 records, identical to 66 single-entry calls. The other 59 single-entry calls returned 404. No duplicate records were seen. [verified]
- **Missing components.** 50N (5C7W label_seq_id 4) and 50L (5C7W label_seq_id 64) are polymer residues of the mapped chain but have no PDBe record. [verified] Both lack a CCD parent (`mon_nstd_parent_comp_id` is empty). The endpoint may only report components with a declared parent, but this is not documented and was not proven. [inference]

## 3. NAKB reference-data assessment

- **Page:** https://www.nakb.org/modifiednt.html, retrieved 2026-10-16. It is a browser of non-standard residues found in PDB nucleic acids with the closest parent nucleotide, stated "as of February 2025". It filters by parent, modification location (base, sugar, phosphate), chirality changes and racemates, and offers atom-mapping details. Saved as `raw/nakb/modifiednt.html`.
- **Data:** https://rna.bgsu.edu/modified/modified_to_change_data.json, linked from the page as the raw data source. HTTP 200, 1,315,013 bytes, `Last-Modified: Fri, 17 Apr 2026 22:31:29 GMT`, ETag `"1410c5-64faf8295f240"`, retrieved 2026-10-16T17:30Z. Saved under `raw/nakb/`. No version number is published, so Last-Modified and ETag are the version identifiers.
- **Format** [verified]: one object keyed by CCD id, 874 components. Each entry has:
  - `standard_base`: a list such as `["G"]`; DNA parents use DA/DC/DG/DT; 52 entries have an empty list.
  - `changes`: atom-level records with `change_type` (addition, removal, replacement, chirality, over_under, added_bond, removed_bond and reversals) and `change_location` (base, ribose, phosphate), plus parent and modified atom names.
  - `pdb`: CCD name, `chem_comp_type`, `par_comp_id`, `one_letter_code`, `pdbx_initial_date`, formula and SMILES.
  - Also `atom_count`, `count`, `base_rmsd` (804 entries), `modomics` (140 entries) and `error` (52 entries).
- **Meaning of `count`:** not documented on the page. It may be the number of PDB entries or occurrences. [inference] It should not be used.
- **Coverage:** all 7 RF00167 components are present. The standard nucleotides A, C, G and U are absent by design. [verified]
- **Parents differ from the CCD for two components.** For 50N and 50L, NAKB `standard_base` is C and G, while CCD `par_comp_id` is null. The NAKB parent is a mapping derived by NAKB, not a declared parent. [verified data; derivation per page text]
- **Scope of use.** NAKB describes component chemistry. It is not evidence that a component occurs in a particular chain. Its descriptors (parent, change locations) can be retained as-is. No classification was built from them.

## 4. mmCIF polymer component assessment

Source: the PDBe updated mmCIF already cached by the pipeline (`https://www.ebi.ac.uk/pdbe/entry-files/download/<pdb>_updated.cif`, cached 2026-10-14), all 66 entries. Script: `scripts/mmcif_inventory.py`. Output: `audit/polymer_positions.csv`, one row for each of the 6,601 polymer positions of the 84 instances.

**Categories present** [verified]:

- `_entity_poly`: the one-letter code with nonstandard components in parentheses, for example `(GDP)GAUC…GGU(CCC)`. In the canonical code, components without a CCD one-letter code become `X` (23G, 50N, 50L).
- `_entity_poly_seq`: one component per position. No microheterogeneity anywhere: `hetero` is never `y` and no position has two `mon_id` values.
- `_pdbx_poly_seq_scheme`: the authoritative join from `label_asym_id` and `label_seq_id` to `mon_id` and author numbering.
- `_atom_site`: coordinates, model 1. Nonstandard residues are `HETATM` records inside the polymer chain.
- `_chem_comp`: component type and name.
- `_pdbx_struct_mod_residue`: present in **only 1 of 66 entries** (1Y27 GDP, parent G). It cannot be used as a completeness source.
- `_pdbx_unobs_or_zero_occ_residues`: confirms unobserved polymer residues, for example 3IVN 23G with polymer_flag Y.

**Consistency checks** (`audit/checks.json`), all with zero discrepancies over 6,601 positions [verified]:

- `_pdbx_poly_seq_scheme` against `_entity_poly_seq`.
- `_pdbx_poly_seq_scheme` against `_atom_site`, where coordinates exist.
- Coordinate presence in mmCIF against the mapping's `has_coordinates`.
- `_pdbx_poly_seq_scheme` against the mapping's `residue_name`, which the pipeline built from the PDBe residue_listing API.
- Nonstandard status against the mapping's `is_modified`.
- Author number and insertion code against the mapping.

**Component type is not membership** [verified]. GTP has `_chem_comp.type` "non-polymer" even when it is the 5′ polymer residue. In 9IWF, 9IWG, 9LKE and 9LMC, GTP is a separate non-polymer ligand recorded in `_pdbx_nonpoly_scheme`. Polymer membership must therefore come from `_pdbx_poly_seq_scheme`, never from the CCD type.

**Missing coordinates** [verified]. 422 of 6,601 polymer residues have no coordinates, yet each still has a deposited component identity. Identity is therefore known for every residue, resolved or not.

## 5. PDBe–NAKB–mmCIF coverage comparison

| Question | mmCIF | PDBe endpoint | NAKB |
|---|---|---|---|
| Finds all 13 RF00167 observations | yes | 11/13 (no 50N, 50L) | n/a |
| Polymer membership explicit | yes (`_pdbx_poly_seq_scheme`) | implicit; no polymer type | no |
| Chain-specific | yes | yes | no |
| Standard (negative) identities | yes, every position | no | no |
| Distinguishes empty from unavailable | yes (file present) | no (404 in both cases) | n/a |
| Unresolved residues | yes (sequence) | yes (23G reported) | n/a |
| Parent nucleotide | `_pdbx_struct_mod_residue` in 1/66 entries; CCD `mon_nstd_parent_comp_id` | no | yes, all 7 |
| Chemical change description | no | no | yes (atom-level, location) |
| Network at run time | none (pipeline cache) | yes | one file |

## 6. RF00167 modification inventory

Counting units are kept separate. Sources: `audit/inventory.json`, `audit/nonstandard_observations.csv`.

| CCD ID | Chemical name (CCD) | Parent: CCD / NAKB | RNA instances | PDB entries | Residues | Rfam positions | Source agreement | Notes |
|---|---|---|---|---|---|---|---|---|
| GTP | GUANOSINE-5'-TRIPHOSPHATE | G / G | 4 | 4 (6P2H, 9V4V, 9V4Z, 9V51) | 4 | none | PDBe ✓ mmCIF ✓ NAKB ✓ | label_seq_id 1, 5′ terminal, outside mapped region. CCD type "non-polymer". Also a non-polymer ligand in 4 other entries, not counted. |
| CCC | CYTIDINE-5'-PHOSPHATE-2',3'-CYCLIC PHOSPHATE | C / C | 2 | 2 (1Y27, 9V51) | 2 | 87 (match) | PDBe ✓ mmCIF ✓ NAKB ✓ | 3′ terminal, with coordinates; NAKB change location ribose. |
| 23G | GUANOSINE-5'-PHOSPHATE-2',3'-CYCLIC PHOSPHATE | G / G | 2 | 1 (3IVN) | 2 | none | PDBe ✓ mmCIF ✓ NAKB ✓ | 3′ terminal, no coordinates, outside mapped region (§7). |
| IU | 5-IODOURIDINE-5'-MONOPHOSPHATE | U / U | 1 | 1 (5UZA) | 2 | 32, 35 (match) | PDBe ✓ mmCIF ✓ NAKB ✓ | internal; NAKB change location base. |
| GDP | GUANOSINE-5'-DIPHOSPHATE | G / G | 1 | 1 (1Y27) | 1 | none | PDBe ✓ mmCIF ✓ NAKB ✓ | 5′ terminal, outside mapped region; the only `_pdbx_struct_mod_residue` row (parent G). |
| 50N | (1S)-1-(6-amino-2-hydroxy-5-nitropyridin-3-yl)-1,4-anhydro-5-O-phosphono-D-ribitol | none / C | 1 | 1 (5C7W) | 1 | 22 (match) | PDBe ✗ mmCIF ✓ NAKB ✓ | internal; mapping parent "X", aligned as N with posterior "*". |
| 50L | (5S)-2-amino-8-(5-O-phosphono-beta-D-ribofuranosyl)imidazo[1,2-a][1,3,5]triazin-4(8H)-one | none / G | 1 | 1 (5C7W) | 1 | 84 (match) | PDBe ✗ mmCIF ✓ NAKB ✓ | internal; mapping parent "X", aligned as N with posterior "*". |

### B7 counts

| # | Item | Value |
|---|---|---|
| 1 | PDB entries examined | 66 |
| 2 | RNA chain instances examined | 84 (6,601 polymer residues; 6,179 with coordinates) |
| 3 | Instances with nonstandard polymer components | 9 (in 8 entries) |
| 4 | Distinct CCD ids | 7 |
| 5 | Nonstandard residue observations | 13 |
| 6 | On consensus match positions | 6 |
| 7 | In insert states | 0 |
| 8 | Outside the mapped Rfam region | 7 (all 5′ or 3′ terminal residues; `not_submitted_flank`) |
| 9 | Distinct consensus positions | 5 (22, 32, 35, 84, 87) |
| 10 | CCD ids reported by PDBe | 23G, CCC, GDP, GTP, IU (11 records) |
| 11 | CCD ids in NAKB | all 7 |
| 12 | CCD ids found only in mmCIF | 50N, 50L |
| 13 | Unresolved residue or entity mappings | none. All 13 resolve to entity 1 (polyribonucleotide), the instance's label_asym_id and label_seq_id. Every PDBe record lands on a mapped instance. |

Other entities: every one of the 66 entries has exactly one polymer entity, a polyribonucleotide, and no protein or DNA polymer. The risk of picking up residues from unrelated polymers is therefore nil in RF00167, but the polymer-entity check is still needed in general. The only same-CCD confusion is non-polymer GTP in 9IWF, 9IWG, 9LKE and 9LMC. [verified, `audit/entry_components.json`]

Shared entities: 3IVN chains A and B share entity 1 and are counted as 2 instances in 1 entry. Each chain was checked separately, since identity is per `label_asym_id`.

## 7. Verification of 23G in 3IVN

Evidence: `audit/23G_3ivn_evidence.json`. All points [verified]:

- **It is a polymer residue, not a ligand.** It appears in `_entity_poly_seq` and `_pdbx_poly_seq_scheme` of entity 1 (polyribonucleotide, "A-riboswitch") at label_seq_id 70, for asym A (author A 74) and asym B (author B 74). It is absent from `_pdbx_nonpoly_scheme`. The one-letter code ends `…CUGGUC(23G)`.
- **Position.** It is the last residue of the 70-nt deposited sequence, so it is 3′ terminal. Its CCD type is "RNA OH 3 prime terminus".
- **No coordinates.** It is listed in `_pdbx_unobs_or_zero_occ_residues` (polymer_flag Y) for both chains, and `_atom_site` has no atoms for it.
- **No Rfam position.** The mapped region of both 3IVN instances is label_seq_id 1–69, with the last mapped residue C69 at consensus 88. Residue 70 is `not_submitted_flank`. Residue number 70 is **not** Rfam position 70, which in both chains is A at label_seq_id 52.
- **NAKB coverage.** Present, with standard_base G, par_comp_id G and changes on the ribose only (the 2′,3′-cyclic phosphate).
- **Agreement.** PDBe (struct_asym_id A/B, residue_number 70, author 74), mmCIF and the pipeline mapping (`residue_name` 23G, `nucleotide_parent` G, `is_modified` true) all agree.
- **Conclusion.** 23G is a nonstandard nucleotide component of the mapped RNA polymer. It should be kept in a residue-level dataset, with its name and terminal, unobserved, outside-mapped-region status. It cannot contribute to any consensus-position track without forcing it onto a consensus position, which the project rules forbid. No biological significance is assigned.

## 8. Residue-to-Rfam mapping validation

For every observation, the chain was verified in this order [verified]:

1. PDB entry
2. Entity 1 (polyribonucleotide)
3. The instance's `label_asym_id`, against the instance list and the mapping file
4. `label_seq_id`, from `_pdbx_poly_seq_scheme`
5. CCD id, in both the sequence and `_atom_site`
6. Alignment state and consensus position, read from the validated mapping (the `residues` list), and independently from the mapping's `consensus_columns` and from `residue_mapping.csv`

No position was inferred from author numbering or from PDBe `residue_number`.

The independent reconstruction (`scripts/independent_check.py`) uses gemmi's structure model rather than CIF tables, the CSV mapping and the per-entry PDBe responses. It reproduces every count in §6 exactly (`audit/independent_check.json`).

**Ambiguity to note.** 50N and 50L are aligned as `N` because PDBe's one-letter sequence gives them parent `X`. They sit at match columns 22 and 84 with posterior `*`, between standard flanking residues A–A and U–U. The positions come from the validated mapping, as the rules require. The aligner placed these residues without base identity. [verified data; the interpretation that placement is driven by flanks is inference]

Positions with no residue: an instance with a deletion at p, or with p outside its mapped span, has no nucleotide at p. For example, 7 instances lack position 87.

## 9. Recommended primary data source

**Approach B, mmCIF-primary,** is the simplest defensible choice:

- **Coverage.** mmCIF finds 13 of 13 observations; PDBe finds 11 of 13.
- **Negatives.** Only the deposited sequence identifies standard residues, which is what makes a negative observation possible.
- **Reproducibility.** It needs no run-time API and reuses the cached mmCIF files the pipeline already depends on.
- **Missing coordinates.** These are handled naturally, since identity comes from the sequence.
- **Maintenance.** It needs one table read (`_pdbx_poly_seq_scheme`) plus `_chem_comp` names.

The pipeline mappings already carry `residue_name` and `is_modified`, identical to mmCIF at all 6,601 positions. A future implementation could read identity from the mapping and verify it against mmCIF. It must **not** use the mapping's `nucleotide_parent`, which is `X` for 50N and 50L.

**Approach A, PDBe-primary,** is rejected. It misses 50N and 50L, cannot separate "no modifications" from "no data", and depends on a live API.

**Roles of the other sources:**

- **PDBe:** cross-validation only. Store the response and record per residue whether PDBe agrees.
- **NAKB:** optional enrichment with `standard_base` and change locations, labelled as NAKB-derived, with the file's Last-Modified date and ETag as version.
- **Parent precedence:** CCD `mon_nstd_parent_comp_id` first, then NAKB `standard_base` (labelled as NAKB), otherwise none.

## 10. Proposed occurrence numerator and denominator

For a consensus position p, using only residues in match state at p in the validated mapping:

- **Positive:** the deposited component at p is not A, C, G or U.
- **Negative:** the deposited component at p is A, C, G or U. Coordinates are **not** required.
- **Missing:** identity cannot be established. Examples are an `UNK` or `N` component, microheterogeneity, or a sequence that conflicts with the coordinates. **None occur in RF00167.**
- **Not counted:** instances with no residue at p (deletion or uncovered).

The definition is F(p) = positives / (positives + negatives). Each instance counts at most once per position, and CCD-level counts are kept per position.

Provisional result (`audit/provisional_occurrence.csv`, not production):

| Rfam position | Positive / denominator | Frequency | CCD | Entries | Denominator residues without coordinates |
|---|---|---|---|---|---|
| 22 | 1 / 83 | 1.2 % | 50N | 1 | 4 |
| 32 | 1 / 84 | 1.2 % | IU | 1 | 4 |
| 35 | 1 / 84 | 1.2 % | IU | 1 | 4 |
| 84 | 1 / 83 | 1.2 % | 50L | 1 | 6 |
| 87 | 2 / 77 | 2.6 % | CCC | 2 | 14 |

All 102 positions have a determinable denominator, ranging from 9 to 84 with a median of 83. 97 positions are 0 of N.

**Points for review**

- **Denominators differ from other tracks.** Because negatives do not need coordinates, these denominators can be larger than those of the loop, pocket and SASA tracks, which require observed coordinates. 98 positions include at least one unresolved residue. This is consistent with "known nucleotide identity", but users comparing tracks may see different N at the same position. Showing the with-coordinates subcount in the pinned details would address it.
- **PDBe silence is not a negative.** An empty PDBe response alone never establishes a negative; negatives come only from the deposited sequence.
- **Constructs, not natural RNAs.** This measures occurrence of nonstandard components in deposited RNA constructs. It says nothing about modification prevalence in natural RNAs. Most observations here are 5′ triphosphate or diphosphate and 3′ cyclic-phosphate termini, iodouridine, and two components without a CCD parent.

## 11. Recommendation for one aggregated frontend track

- **One category with one track is appropriate.** Per-CCD tracks are not justified: there are 7 components, most occur once, and each mapped position carries a single CCD.
- **The colour scale is the open issue.** On the existing five-bin frequency scale (<20 %, …, ≥80 %), all 5 marks fall in the lightest bin and will be hard to see. Options, for a decision before implementation:
  - (a) Keep the shared frequency scale for consistency and accept faint marks.
  - (b) Use a dedicated scale with finer low-end bins.
  - (c) Use a presence mark of fixed colour, with the frequency in the tooltip.

  Option (c) best matches the sparse data.
- **Tooltip:** `Rfam residue 87` / `Modification occurrence: 2 of 77 RNA instances (2.6%)`.
- **Pinned details:**
  - CCD id with its name and parent, labelled CCD or NAKB.
  - Supporting instances and entries per CCD.
  - The count of denominator residues without coordinates.
  - PDBe agreement.
- **Off-axis observations.** The 7 terminal observations outside the mapped region would not appear on the axis. They should stay in `modifications.json`, and a one-line note in the details or popover should say they exist.
- **Popover wording.** The proposed text is accurate. Suggested, slightly tighter:

  > Nonstandard nucleotide components in deposited RNA polymer sequences, mapped to Rfam residue positions.
  >
  > Colours indicate occurrence across RNA instances. These components do not necessarily represent naturally occurring RNA modifications.

- **Naming.** Most RF00167 observations are terminal phosphate states or synthetic analogues, so "Nonstandard nucleotides" may describe the content more faithfully than "Modified nucleotides". This is the user's decision.

## 12. Open questions and limitations

1. Colour encoding for very low frequencies (§11).
2. Category name: "Modified nucleotides" or "Nonstandard nucleotides".
3. Should 5′ and 3′ terminal components that do map, such as CCC at 87, be shown alongside internal ones, or flagged as terminal in the details? This report recommends a factual "terminal" flag rather than exclusion.
4. Denominators may include unresolved residues (§10). The alternative is to restrict to residues with coordinates to match the other tracks. The data support either; the "known identity" definition is more faithful to B8.
5. PDBe's inclusion criterion is undocumented; the parent-based explanation is inference. NAKB `count` is undocumented.
6. 50N and 50L were aligned without base identity (`N`). Their positions follow the validated mapping and were not re-examined by re-alignment.
7. Only RF00167 was examined. Families with natural modifications, such as tRNA or rRNA, will have far more records, `_pdbx_struct_mod_residue` may be more common, and microheterogeneity may occur. The "missing" rules are untested on real cases here.

## 13. Minimal future implementation plan (not started; needs approval)

1. **Command** `python3 -m rna_explorer modifications RF00167`. For each instance, read `_pdbx_poly_seq_scheme` and `_chem_comp` from the cached updated mmCIF, and join to the mapping by `label_asym_id` and `label_seq_id`. Assert mmCIF and mapping `residue_name` agree.
2. **`modifications.json`**, one record per nonstandard residue, kept even when unmapped:
   - Instance identity: `seq_name`, pdb_id, entity_id, label_asym_id, label_seq_id, auth_asym_id, auth_seq_id, ins_code.
   - Component: ccd_id, ccd_name, ccd_type, ccd_parent, nakb_standard_base, nakb_change_locations.
   - Position: alignment_state, rfam_consensus_pos (or null), is_5prime_terminal, is_3prime_terminal, has_coordinates.
   - Cross-check: pdbe_reported (true/false), plus flags.
   - Instance-level status for all 84 instances: the count of standard residues assessed, and missing reasons.
   - Provenance: mmCIF source and date, PDBe URL, retrieval date and raw-response hash, NAKB URL with Last-Modified and ETag, Rfam/CM releases, pipeline and config versions.
3. **`modification_tracks.json`**, per consensus position:
   - n_positive, n_negative, n_missing, n_no_residue, denominator, frequency.
   - n_denominator_without_coordinates, ccd_counts, n_supporting_entries, pdbe_agreement counts.
   - Summary: off-axis observations, by state.
4. **Schemas and tests.** JSON schemas plus pytest:
   - Synthetic unit tests: microheterogeneity, UNK, a non-polymer copy of the same CCD, an insert-state residue, an unresolved residue.
   - RF00167 regression: 13 observations, 6 on match positions, 5 positions, 23G off-axis, 50N/50L not in PDBe.
5. **Frontend**, after a decision on §12 items 1–4: one category and one track, a popover, tooltip and pinned details, two download items, and Playwright and vitest coverage.

## Files in this directory

- `scripts/fetch_pdbe.py`: read-only PDBe fetch (batch and per entry).
- `scripts/mmcif_inventory.py`: per-position mmCIF audit and consistency checks.
- `scripts/occurrence_and_inventory.py`: inventory and provisional counts.
- `scripts/independent_check.py`: independent reconstruction.
- `raw/pdbe/`: all API responses, including probes (3IVN/3RKF, invalid id, 7K00, 1EHZ, 122D, protein entries).
- `raw/nakb/`: page, JSON, headers and retrieval time.
- `raw/ccd/`: CCD definitions of the 7 components (copied from the pipeline cache, wwPDB CCD via files.rcsb.org).
- `audit/`: `polymer_positions.csv` (6,601 rows), `nonstandard_observations.*`, `entry_components.json`, `checks.json`, `inventory.json`, `provisional_occurrence.csv`, `independent_check.json`, `23G_3ivn_evidence.json`.
