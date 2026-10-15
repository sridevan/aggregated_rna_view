# RNA Structure Explorer: Milestone 1 design

Status 2026-10-14: implemented as described below; RF00167 built (66 entries, 84 instances, 29 unique sequences, 84/84 validation ok) and regression-tested against the POC (434/434 residue rows identical). Deviations from the proposal: none in structure; `get_molecules` falls back to per-entry GET when the PDBe batch endpoint refuses a request.

## 1. What exists and what is reused

| POC code | Reused as | Notes |
|---|---|---|
| `poc_rf00177/scripts/common.py` (`cif_category`, `atom_site_observed`, `entity_source`) | `rna_explorer/mmcif.py` | unchanged logic, typed, plus `poly_seq_scheme()` helper |
| `poc_rf00177/scripts/03_parse_stockholm.py` (RF column typing: letter = match, `.` = insert, `~` = local end; anchor/offset; deletion vs uncovered) | `rna_explorer/stockholm.py` + `rna_explorer/mapping.py` | validated against RF00167, RF00177 and BGSU; logic kept, split into parse (columns) and map (residues) |
| `poc_rf00177/scripts/02_build_inputs.py` (mapped-region extraction, parent-base normalisation, species from mmCIF) | `rna_explorer/sequences.py` | generalised to all instances of a family, deduplicated by mapped sequence |
| `poc_rf00177/scripts/04_validate.py` (mmCIF identifier checks, coordinates cross-check, continuity, insertion bookkeeping) | `rna_explorer/validation.py` | unchanged checks, returns typed results |
| `poc_rf00177/scripts/01_instances.py`, `01b_clan_competition.py` (PDBe mapping list, organism via search API, Rfam structures endpoint, clan competition) | `rna_explorer/pdbe.py`, `rna_explorer/rfam.py`, `rna_explorer/instances.py` | same endpoints; cached |
| cmalign invocation + `--sfile` parsing | `rna_explorer/infernal.py` | local mode (default), `--cpu`, no `-g` (RF00177 showed global mode collapses divergent chains) |
| PP handling, window flags (`06_reference_check.py`) | `rna_explorer/confidence.py` (PP only in Milestone 1) | window flags/reference checks stay in the POC; not part of the export |

Not reused: BGSU scripts (validation only), insertion-profile plotting, reference-motif checks.

## 2. Package architecture (minimal)

```
rna_explorer/            Python package (no framework; stdlib + requests, jsonschema for tests)
  config.py              paths (cache dir, output dir), cmalign discovery (env RNA_EXPLORER_CMALIGN, tools/infernal-1.1.5/bin, PATH)
  cache.py               fetch(url) -> cached file; JSON/text/binary; cache/<host>/<sha1 of url>  (+ readable name)
  rfam.py                family info, CM download + header (CLEN, CKSUM, DATE), structures endpoint, SSU/any clan membership, release
  pdbe.py                nucleic_mappings, molecules (POST batch), residue_listing, search API (organism/method/resolution), updated mmCIF
  mmcif.py               category reader, poly_seq_scheme, atom_site observed set, entity source
  instances.py           Instance records (one per mapped chain) with organism, method, resolution, Rfam bits, clan-best flag
  sequences.py           mapped-region extraction, parent-base normalisation, dedup -> UniqueSequence (md5 key) + FASTA
  infernal.py            cmalign runner (subprocess), scores parser, version string
  stockholm.py           Stockholm parser, column typing (match/insert/local_end), consensus numbering
  mapping.py             ResidueMapping / ConsensusColumn dataclasses; build per-chain mapping from unique-sequence alignment + residue listing + coordinates
  validation.py          mmCIF identifier agreement, coordinate cross-check, continuity, insertion bookkeeping, span vs Rfam
  confidence.py          PP bins, per-residue PP, per-consensus statistics (sequence-level deduplicated vs structure-level)
  export.py              writes frontend JSON + CSV, validates against JSON Schemas in rna_explorer/schemas/
  pipeline.py            build_family(acc): fetch -> instances -> unique sequences -> cmalign -> mapping -> validate -> confidence -> export
  cli.py                 python -m rna_explorer build RF00167 [--out data] [--cache cache]
  schemas/*.schema.json  JSON Schema (draft 2020-12) for family.json, instances.json, consensus_stats.json, mapping json
tests/                   pytest; offline regression tests use POC Stockholm/chains/mmCIF files already on disk
tools/infernal-1.1.5/    cmalign binary built from source (persisted from the POC session)
data/families/<ACC>/     pipeline output consumed by the frontend (static JSON)
cache/                   downloaded Rfam/PDBe data (gitignored)
```

Alignment unit = unique mapped sequence (entity-level dedup, md5 of the parent-base mapped region). Every chain instance with that sequence receives the same alignment; chain-specific information (author numbering, insertion codes, coordinates, modified residue names) is joined per chain from the residue listing and mmCIF.

## 3. Schemas

### Residue mapping (`mappings/<pdb>_<entity>_<auth_chain>.json`, one per chain instance)
```
{ "pdb_id", "entity_id", "label_asym_id", "auth_asym_id", "rfam_acc", "sequence_key",
  "mapped_region": {"start", "end"}, "alignment": {"bit_score", "avg_pp", "trunc", "cm_from", "cm_to"},
  "residues": [ { "polymer_seq_position", "author_residue_number", "author_insertion_code", "residue_name",
                  "nucleotide_parent", "is_modified", "alignment_state": match|insert|local_end|not_submitted_flank,
                  "rfam_consensus_pos" (int|null), "insertion_anchor" (int|null), "insertion_offset" (int|null),
                  "posterior_prob_char", "posterior_prob" (float|null, bin midpoint), "has_coordinates" (bool) } ],
  "consensus_columns": [ { "rfam_consensus_pos", "state": match|deletion|uncovered, "polymer_seq_position" (int|null) } ],
  "provenance": { "rfam_release", "cm_checksum", "cm_clen", "infernal_version", "cmalign_options", "generated" } }
```
CSV twins (`residue_mapping.csv`, `consensus_columns.csv`) with the same columns as the POC for regression.

### Confidence statistics (`consensus_stats.json`)
Per consensus position 1..CLEN, arrays of equal length CLEN:
```
{ "rfam_acc", "clen", "rf": [...], "ss_cons": [...],
  "sequence_level": { "n_unique_sequences", "n_match": [], "n_deletion": [], "n_uncovered": [], "coverage": [],
                      "pp_mean": [], "pp_n": [], "pp_min": [], "pp_high": [], "pp_medium": [], "pp_low": [] },
  "structure_level": { "n_instances", "n_match": [], "n_match_with_coordinates": [], "coverage": [], "structural_coverage": [],
                       "pp_mean": [], "pp_n": [], "pp_high": [], "pp_medium": [], "pp_low": [] },
  "pp_bins": {"*": 0.975, "9": 0.9, ...}, "pp_display_thresholds": {"high": 0.95, "medium": 0.80}, "provenance": {...} }
```
No composite score. Categories: high = bin midpoint >= 0.95 (`*`), medium = 0.80 <= midpoint < 0.95 (`9`, `8`), low otherwise; provisional, for display only.

### Instances (`instances.json`) and family (`family.json`)
Instances: one row per mapped chain: identifiers, organism/tax id (search API, mmCIF fallback), experimental method, resolution, release year, title, mapped region, Rfam bit score/E-value/model span, clan-best flag, sequence_key, n_match, n_insert, n_local_end, n_match_with_coordinates, model coverage, avg PP, mapping file path.
Family: accession, id, description, clen, consensus (RF) sequence, SS_cons, counts (entries, instances, unique sequences, methods), representative instance (documented rule: clan-best, highest n_match_with_coordinates, then best resolution), coverage summary, provenance.

## 4. How the frontend JSON is generated

`python -m rna_explorer build RF00167 --out data --cache cache` runs `pipeline.build_family`: all downloads go through `cache.fetch` (re-runs are offline); cmalign runs once on the deduplicated FASTA; mapping/validation/confidence are computed in memory; `export.write_family` writes JSON (validated against the schemas) and CSVs under `data/families/RF00167/`. The frontend (Milestone 2) reads those static files; no API and no Infernal at runtime.

## 5. Tests

- `test_stockholm.py`: synthetic alignment covering `.`/`~` columns, anchor 0, trailing insertion, deletion vs uncovered.
- `test_mapping_regression.py`: rebuild mappings from the POC Stockholm + chains files for RF00167 (5 chains) and RF00177 (8 chains) and compare row by row with the POC CSVs (offline).
- `test_mmcif.py`: category reader on a POC mmCIF (loop + key-value), observed set.
- `test_confidence.py`: PP bins, dedup, coverage arithmetic on a toy family.
- `test_pipeline_rf00167.py`: end-to-end on cached RF00167 data when cmalign is available (skipped otherwise); compares the 5 POC chains with `poc_rf00167/results/mapped_residue_mapping.csv`.
