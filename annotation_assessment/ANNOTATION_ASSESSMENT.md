# Milestone 3A: RNA annotation data assessment for RF00167

Date: 2026-10-15. Data investigation only; no backend, frontend or orientation code changed. Inputs: the validated `data/families/RF00167/` mappings (66 entries, 84 chain instances, residue → Rfam consensus with PP and coordinate flags). Raw responses are cached in `cache/bgsu/` and `cache/pdbe/`; `scripts/assess.py` regenerates every table.

## 1. BGSU helix / structural-element annotations

| Item | Finding |
|---|---|
| Family file | `RF00167_annotations.txt`: one curated reference chain, **4FEO chain B** (B. subtilis xpt-pbuX guanine riboswitch), 67 nucleotides, element names following the cited paper (PMID via sciencedirect S0022283613001253): Helix P1, J1/2, Helix P2, L2, J2/3, Helix P3, L3, J3/1. Unit ids `PDB|model|chain|comp|author number`, one label per nucleotide. |
| Per-chain endpoint | `nucleotide_annotation?chain=PDB|1|CHAIN` transfers the reference labels through BGSU's Rfam-based correspondence. Fetched for all 84 instances: **80 chains annotated (65 entries)**; the four 9NBH chains return empty responses (their aptamer core has no coordinates). 65–67 rows per chain, 5,283 nucleotides in total. |
| Multiple assignments | None observed in RF00167 (one reference chain, one label per nucleotide). The parser keeps the full label string and takes the first entry as the preferred assignment, as BGSU recommends, should other families carry several references. |
| Identifier mapping | 100 % of 5,283 unit ids resolved to `(pdb_id, label_asym_id, label_seq_id)` via author number + insertion code; all are match-state residues with coordinates (status `mapped`); 0 unresolved, 0 insert/local-end (BGSU only annotates nucleotides that align to the reference). PP confidence: 5,134 high, 102 medium, 47 low. |
| Coverage | 67 of 102 consensus positions annotated (the aptamer core 19–87). Element → consensus span: P1 19–25/81–87, J1/2 26–28, P2 29–35/43–49, L2 36–42, J2/3 50–57, P3 58–64/72–78, L3 65–71, J3/1 79–80. Every consensus position receives exactly one element label across all instances (no disagreement), which confirms the labels are reference-transferred rather than per-structure helix observations. |
| Provenance | reference chain, family file URL, per-chain endpoint and retrieval date are recorded per row. |

Caveat: these are alignment-transferred labels (P1 in every chain because P1 in 4FEO), not structure-specific base-pair detection. A per-structure helix track would need FR3D/DSSR pairing, which is out of scope here; the BGSU PDBe basepair JSON used by the 2D viewer could provide it later.

## 2. BGSU RNA 3D Motif Atlas loops

Sources: `rna3dhub/loops/download/{PDB}` (loop id → nucleotide unit ids, complete membership including discontinuous strands) and `rna3dhub/pdb/{PDB}/motifs` (loop → motif-group id and annotation text, parsed from the HTML table; no JSON endpoint exists). All 66 entries returned data.

| Loop type | Loops | With motif group | Without | Distinct groups | PDB entries | Instances | Consensus positions |
|---|---|---|---|---|---|---|---|
| HL | 234 | 211 | 23 | 10 | 66 | 84 | 42 |
| IL | 82 | 76 | 6 | 7 | 46 | 50 | 33 |
| J3 | 114 | 109 | 5 | 9 | 64 | 80 | 42 |
| J4 | 66 | 60 | 6 | 3 | 52 | 66 | 15 |
| total | 496 | 456 | 40 | 29 | 66 | 84 | |

- 5,655 loop nucleotide memberships; 0 unresolved identifiers; 5,419 mapped to consensus positions (5,058 high, 242 medium, 119 low PP); 68 insert-state residues (IL/J3 in the 2025 entries 9IWF, 9IWG, 9LMC, 9LJN, 9LKC: residues inserted relative to the model); 168 nucleotides outside the Rfam-mapped region (the scaffold of 9NBH/9NDD and the 6P2H flanks). All loops lie on RF00167 chains; 4 loops span two chains (crystal contacts in multi-copy entries).
- Motif groups: HL_09122.3 (L3, 72 loops), HL_59138.3 and HL_35704.1 (L2, 52 + 21), HL_56046.1 "Purine riboswitch" (53), HL_90436.7 "Mini UNCG" (4), HL_34789.7 "GNRA" (1, consensus 99–102 of a scaffold); J4_07491.4 "Purine riboswitch junction" (54); J3_44163.4 (53) and eight minor J3 groups; IL_19025.1 / IL_01003.7 "Isolated non-canonical cWW pair", IL_48076.11 "Major groove platform", three others.
- Scientific ambiguities found: (a) FR3D classifies the J2/3 strand (consensus 50–57, closed by a non-canonical G·C pair) as a *hairpin* loop HL_56046.1, so "hairpin" counts at 50–57 are junction residues biologically; (b) 62 of 84 chains carry **both** a J3 and a J4 loop for the same three-way junction (J3 includes the 50–57 strand; J4 excludes it), so J3 and J4 memberships overlap at consensus 25–29, 49–50, 57–58, 78–81 and must not be summed; (c) loops without a motif group (40) are usually loops the Atlas release excluded for quality reasons, not loops without structure.

## 3. PDBe bound molecules

Source: `api/v2/pdb/entry/binding_sites/{pdb}/{entity}` (entity 1 for all RF00167 entries; all 66 returned 200). Response = one record per ligand CCD id with `residues` ranges (`startIndex`/`endIndex`, `indexType: "PDB"`) and `boundMoleculeId` per range; `additionalData.bestChainId` (a **label_asym_id**) and `entityId`.

- Numbering: the indices are **polymer sequence positions (label_seq_id)**, verified for 1Y26 (ADE contacts at indices 9–10, 35, 39–40, 62–63 whose sequence codes match the entity sequence; author numbering of 1Y26 is offset by 12). Mapping to canonical ids is therefore direct.
- Observation level: **entity-level**. Residues of all bound-molecule instances of a ligand are pooled; `boundMoleculeId` distinguishes instances but is not tied to a chain in this endpoint. For the 12 multi-chain entities (3RKF, 9NBH, 9NDD with 4 chains; nine with 2) instance-specific attribution is not determinable without Arpeggio/bound-molecule endpoints. The assessment attributes contacts to `bestChainId` (instance counts are therefore lower bounds; entry counts are exact).
- 37 distinct CCD ids: 30 ligands, 5 metal ions (MG 30 entries, NCO 34, K 5, MN 1, SR 1), 2 additives (ACT 27 entries, BR 1); no water records in this endpoint. 63 entries have a ligand-class molecule; 3IVN, 5E54 and 9NBH have none (apo/metal-only). 1,492 binding-residue rows: ligand class 497 (493 mapped, 4 outside the mapped region, 6P2H flank), metal ions 849 (834 mapped, 15 outside: 6P2H flanks and 5UZA MG at the 3′ end), additives 146 (all mapped).
- The ligand pocket is remarkably consistent. Entries (of 63 ligand-bound) contacting each consensus position: 25 and 26: 62; 55: 62; 80: 61; 56 and 81: 60; 51: 54; 79: 24; 54: 7; 53: 4; 52 and 19: 3; 20: 2; 27, 28, 82, 87: 1 (bulky analogues such as CMG and ANG; GTP at 19–20 is a 5′-terminal nucleotide artefact, not a pocket ligand). HPA (12 entries), ADE (9), GTP and GUN (4), 6GU, 6AP, GNG, OXG, XAN (3), 6GO and ANG (2) are the most frequent; the other 19 analogues occur once.

## 4. Identifier validation summary (`identifier_mapping_issues.csv`)

| Source | Unit ids | Mapped (with coords) | Mapped, no coords | Insert | Local end | Outside mapped region | Unresolved |
|---|---|---|---|---|---|---|---|
| BGSU helix | 5,283 | 5,283 | 0 | 0 | 0 | 0 | 0 |
| BGSU loops | 5,655 | 5,419 | 0 | 68 | 0 | 168 | 0 |
| PDBe binding sites | 1,492 | 1,473 | 0 | 0 | 0 | 19 | 0 |

No annotation is dropped: every non-mapped nucleotide is listed with its reason. Low-confidence consensus assignments (PP ≤ 0.75) are 47 / 119 / 1 rows respectively and are flagged per row.

## 5. Coverage summary
See `annotation_coverage.csv` / `.json` (one row per source with entries, instances, nucleotides, unique consensus positions, unresolved and missing counts; loop categories with motif-group counts; ligands with entry/instance counts and consensus positions).

## 6. Feasibility of the five tracks (per consensus position)

| Track | Occurrence counts | Entries / instances | Distinct motif groups | Provenance / confidence | Notes |
|---|---|---|---|---|---|
| Helices (P1–P3 + J/L elements) | yes (reference-transferred) | yes | n/a | reference chain + PP per residue | denominator = 80 annotated instances / 65 entries, not 84; absence at a position = position not aligned to the reference, not "no helix" |
| Hairpin loops | yes | yes | yes (10 HL groups) | loop id, motif group, PP | exclude or flag HL_56046.1 (J2/3 strand) when presenting "hairpins" |
| Internal loops | yes | yes (46 entries) | yes (7) | idem | IL present only in some constructs; absence in the other 20 entries is an experimentally supported negative (loops were computed for every entry) |
| Junction loops | yes, per type | yes | yes (9 J3, 3 J4) | idem | count J3 and J4 separately; never add them |
| Ligand-binding sites | yes (entry level) | entries exact, instances lower bound | n/a (by CCD id / class) | PDBe binding_sites, bestChainId | default display excludes metals/additives; apo entries are true negatives |

Denominators: helices 80 instances / 65 entries; loops and ligands 84 instances / 66 entries where the source returned a record for the entry (all did), with per-position occupancy also limited to residues present with coordinates in each instance (positions deleted or unmodelled in an instance are "not observable", distinct from "observed without annotation").

## 7. Answers

1. Completeness: yes for a demonstration. Helix labels cover 80/84 instances; loops cover all 66 entries and 84 instances with 92 % motif-group assignment; binding sites exist for all 66 entries and 63 carry a ligand. Gaps are explainable (9NBH unmodelled core; 2025 insertions).
2. Identifier mapping: yes. 0 unresolved identifiers across 12,430 nucleotide references; author-number+insertion-code (BGSU) and label_seq_id (PDBe) both map deterministically through the existing mapping files; non-mapped cases are insertions, scaffold residues outside the Rfam region, or flanks, all retained with reasons.
3. Distinct Motif Atlas groups per consensus position: yes, from loop membership × motif group (e.g. L2 positions 36–42 carry HL_59138.3 and HL_35704.1; junction 25–29/78–81 carry J3_44163.4 and J4_07491.4).
4. Ligand aggregation without Arpeggio: yes at PDB-entry level with the binding_sites endpoint (CCD id, name, bound-molecule id, label_seq_id ranges). Chain-specific attribution in multi-copy entities and interaction types need the bound-molecule/Arpeggio endpoints and are deferred. GTP at positions 19–20 shows that the endpoint also reports covalently attached 5′ nucleotides as bound molecules, so a small curated exclusion list is needed.
5. To resolve before implementation: (a) J3/J4 double annotation of the same junction and FR3D "hairpin" at J2/3: choose per-category counting and label the J2/3 hairpin group explicitly; (b) reference-transferred helices vs observed base pairs: present the helix track as "reference element assignment (4FEO)" and plan an observed-pairing track from FR3D data later; (c) entity-level ligands: report instances as lower bounds or use `bestChainId` only; (d) metal/additive classification lists are curated here (`scripts/assess.py`) and should be reviewed; (e) motif-group ids are release-specific (e.g. HL_09122.3): store the Atlas release and treat ids as versioned; (f) scaffold residues outside the mapped region (9NBH/9NDD) carry loops that must stay out of family tracks but could be shown per instance.
