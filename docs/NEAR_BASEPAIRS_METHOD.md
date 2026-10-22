# Near base pairs (local fr3d-python) and the Rfam consensus comparison

## Why a local run

BGSU RNA 3D Hub's per-entry base-pair CSVs, the production source of true base pairs, contain **no near base pairs**:

- None of the 286 cached base-pair files has an `n`-prefixed class. The stacking, base–ribose and base–phosphate files do have them.
- The PDBe per-chain base-pair JSON has none either.
- fr3d-python, BGSU's own annotation code, writes near base pairs only when the category list includes `near`. See
  `NA_pairwise_interactions.py` at commit 3feb9b9, line 4233.

So near base pairs are computed locally with the same code. This was validated in `experiments/sscons_vs_3d/NEAR_CHECK.md`: the local true
base pairs are identical to BGSU's in 298 of 321 chains, and differ by 1–2 pairs in 22 chains. The remaining chain is 7PJS a, where the
local run also annotates the region BGSU leaves unannotated.

## Environment

The checkout and the interpreter are kept under `tools/`, which is git-ignored:

    git clone https://github.com/BGSU-RNA/fr3d-python.git tools/fr3d-python
    git -C tools/fr3d-python checkout 3feb9b93e179187f45e1b79f3842d1a0ca19d797
    python3.11 -m venv tools/fr3d-env
    tools/fr3d-env/bin/python -m pip install numpy==2.4.6 scipy==1.17.1 biopython==1.88 mmcif-pdbx==2.1.0 mmcif==1.2.1

**Provenance.** `rna_explorer near` refuses a checkout with local modifications. It records in `near_basepairs.json`:

- the commit and commit date;
- the interpreter and package versions;
- the licence note.

**Licence.** fr3d-python's `license.txt` is MIT (Copyright 2026 BGSU RNA Group), but its `pyproject.toml` declares Apache-2.0. Both are
permissive, and both are recorded.

## Commands

    python3 -m rna_explorer near RF00167          # near_basepairs.json (default --fr3d tools/fr3d-python --python tools/fr3d-env/bin/python)
    python3 -m rna_explorer interactions RF00167 --seed-dir experiments/rf00167_interactions/raw   # merges them into interactions.json

**What `near` does, per entry:**
1. Streams the PDBe updated mmCIF into a temporary directory (`--tmp`).
2. Runs `NA_pairwise_interactions.py -c basepair,near --chain <family author chains>`.
3. Keeps only near base pairs whose two residues are both in family chains.
4. Deletes the mmCIF. Full mmCIFs are never stored.

Entries already annotated with the same commit and chain set are reused; `--force` re-annotates every entry. A failed entry is recorded as
`failed`, and its chains are `near_basepairs: unavailable` in `interactions.json`.

| Family | Entries | Near base-pair rows |
|---|---|---|
| RF00167 | 66 | 684 |
| RF00162 | 28 | 530 |
| RF00059 | 28 | 324 |
| RF00005 | 59 | 544 |
| RF00011 | 50 | 3,096 |
| RF00177 | 58 | 15,344 |

All runs are from 2026-10-22, with no failures.

## How `interactions` uses them (schema 1.4)

**Joining.** Near rows are joined to Rfam positions exactly like the BGSU rows: the author chain, then the author number plus insertion
code, then the validated mapping. Each endpoint is checked against the mmCIF. The rows are tagged `fr3d_local` in the last observation
field; BGSU rows are tagged `bgsu`.

**What near rows never change:**
- BGSU file statistics;
- chain annotation availability;
- annotation gaps (the 7PJS a gap is unchanged);
- `n_observed` or `frequency`;
- 3D-only rows. A non-consensus pair needs a true base pair to get a row, so near-only non-consensus pairs have none.

**Base pairs**
- `n_paired` = `n_observed` + `n_near_only`, which counts a true or near FR3D base pair of any class.
- `frequency_paired` and `n_entries_paired` follow from it.
- `near_classes` gives the near-only chains per near class.

**`consensus`**, comparing each base pair with the Rfam SS_cons:

| Value | Meaning |
|---|---|
| `both` | SS_cons pair paired in at least one eligible chain |
| `2d` | SS_cons pair paired in no eligible chain. It gets a row even then, including when no chain can evaluate it (`n_eligible` 0, labelled "not evaluable in any chain"). |
| `3d` | Observed base pair not in SS_cons |

**Other categories:** `n_paired` = `n_observed`, and `consensus` is null. Their near interactions come from BGSU and are still not
counted.

## Frontend

**Family table and arcs**
- Conservation counts paired chains (true or near).
- Near classes appear in the ▸ class breakdown.
- A new **Rfam consensus** column shows Both, 2D or 3D, and can be sorted.
- **More filters › Rfam consensus** sets the same filter for the table and the arcs. Like the other base-pair-only filters, it hides other
  categories while it is set.
- "2D" arcs are grey and dashed (`#a8a29e`, `2 3`).

**Per chain**
- A near pair has the state **Near**, shown with its class.
- In Compare two it counts as paired: a near pair against a true pair is "Different class", and a near pair against no pair is "Only in".
- Helix coverage counts near pairs.
