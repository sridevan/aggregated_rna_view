# Near base pairs: check (2026-10-22)

**Question.** Are near base pairs really unavailable? What do they show for the SS_cons pairs that have no true FR3D base pair?

**Status.** Investigation only. The app, backend and data are unchanged.

**How to run.** Script `near_check.py`; outputs in `followup/near/`.

    python3 experiments/sscons_vs_3d/near_check.py --fr3d <fr3d-python checkout> --python <interpreter> --workers 4

## 1. Where near base pairs are (and are not)

| Source | Near base pairs? | Evidence |
|---|---|---|
| BGSU per-entry CSV `interactions/fr3d/basepairs/csv` (our production source) | **No** | 0 of 286 cached files contain an `n` class. Stacking, base-ribose and base-phosphate files all do. |
| BGSU web pages / other endpoints | Not checked | rna.bgsu.edu returned HTTP 503 for every URL on 2026-10-22. |
| PDBe `static/entry/{pdb}_{chain}_basepair.json` | **No** | 1EHZ A: 30 pairs, none near. It does list G18–Ψ55 tWS, which BGSU's CSV lacks. |
| **fr3d-python** (BGSU's annotation code), commit 3feb9b9 | **Yes, on request** | Near pairs are written only when the category list includes `near`: `-c basepair,near`. The code skips `n` classes for plain `basepair` (`NA_pairwise_interactions.py`, line 4233). This explains why the published base-pair files have none. |

**Conclusion.** Near base pairs are not missing because nothing was nearly paired. The published base-pair export simply leaves them out, so they can only be obtained by running fr3d-python locally.

## 2. Local run and validation against BGSU

**Setup.**
- fr3d-python 3feb9b9, `-c basepair,near`, restricted to the family chains (`--chain`).
- Coordinates: PDBe updated mmCIF, model 1. Each file was downloaded to a temporary scratch folder and deleted after the run.
- All 289 entries, 321 chains; no failures.
- Endpoints were joined to Rfam positions the same way as production: author chain, author number + insertion code, validated mapping, match states.

**True (non-near) base pairs, local vs BGSU,** on unordered consensus position pairs + class:

| Family | Chains | Local | BGSU | Shared | Chains identical |
|---|---|---|---|---|---|
| RF00167 | 84 | 2,398 | 2,397 | 2,397 | 83 |
| RF00162 | 29 | 1,091 | 1,091 | 1,091 | 29 |
| RF00059 | 38 | 1,134 | 1,135 | 1,134 | 37 |
| RF00005 | 60 | 1,538 | 1,537 | 1,536 | 57 |
| RF00011 | 50 | 6,954 | 6,954 | 6,954 | 50 |
| RF00177 | 60 | 36,135 | 35,733 | 35,729 | 42 |

**Differences.**
- **7PJS a:** 387 pairs are local-only. BGSU's annotation of this chain stops near position 650; this is the known 652–1531 annotation gap, and the local run fills it.
- **Everywhere else**, differences are 1–2 cWW pairs per chain, in 22 chains. They are probably small version or threshold differences between BGSU's production run and the current code.

**So the local run reproduces the production source**, and its near pairs can be used alongside the BGSU true pairs.

## 3. Near pairs on SS_cons pairs without a true BGSU base pair

| Family | Observations without a BGSU pair | Near pair between i and j | No pair at all | True pair in the local run only | Near classes |
|---|---|---|---|---|---|
| RF00167 | 75 | 54 | 20 | 1 | ncWW 49, ntWW 5 |
| RF00162 | 56 | 23 | 33 | 0 | ncWW 21, ntSH 2 |
| RF00059 | 90 | 17 | 73 | 0 | ncWW 16, ntSH 1 |
| RF00005 | 100 | 64 | 34 | 2 | ncWW 55, others 9 |
| RF00011 | 494 | 193 | 301 | 0 | ncWW 121, ntSH 29, ncWS 26, others 14 |
| RF00177 | 2,550 | 1,218 | 1,315 | 17 | ncWW 1,069, others 149 |

**Agreement with the geometry check** (`followup/unpaired_geometry.csv`; five families, 16S not measured):

| Geometry | Near pair | No pair | Local true pair |
|---|---|---|---|
| WC-like (N–N ≤ 3.6 Å, C1′–C1′ 8.5–12.5 Å) | 146 | 27 | 3 |
| Distorted (N–N ≤ 5 Å) | 143 | 105 | 0 |
| Apart | 62 | 329 | 0 |

**Reading the geometry table.**
- **Near pairs are the FR3D-native form of the "paired but not annotated" state.**
  - 83% of WC-like observations are near pairs.
  - Most "apart" observations have no pair at all.
- **The 27 WC-like observations without a near pair** are mostly RNase P (hairpin 135–155 and 328–338 regions). These are likely cases where the base planes are not coplanar enough even for a near call.
- **The geometry measure is not needed** where the near annotation exists. Near is the better choice, because it uses FR3D's own criteria.

## 4. Effect on support (near shown separately, never counted as formed)

| Family | cWW | cWW or ncWW | Pairs ≥ 0.9 (cWW → cWW or ncWW) | Pairs < 0.5 (cWW → cWW or ncWW) |
|---|---|---|---|---|
| RF00167 | 0.948 | 0.979 | 18 → 21 | 0 → 0 |
| RF00162 | 0.900 | 0.926 | 24 → 28 | 4 → 3 |
| RF00059 | 0.876 | 0.897 | 19 → 19 | 3 → 3 |
| RF00005 | 0.916 | 0.959 | 17 → 20 | 0 → 0 |
| RF00011 | 0.865 | 0.890 | 72 → 75 | 8 → 6 |
| RF00177 | 0.898 | 0.938 | 375 → 418 | 32 → 22 |

66 SS_cons pairs change band.

## 5. The validation cases revisited

| Case | Before (cWW, no near) | Near pairs (chains without a true pair) | Updated reading |
|---|---|---|---|
| RNase P 110–115 (U·G) | 19 of 50 | ncWW in **31 of 31** | Formed in every chain: 19 cWW + 31 near cWW. |
| RNase P 7–361 (G·U) | 0 cWW; cWS 27 | ncWS 21, none 2 | Consistently a Watson–Crick / Sugar-edge contact, not a cWW pair. |
| RNase P 324–342 | tSH 37 | ntSH 6, ncSW 2, none 5 | Sheared. |
| RNase P 325–341 / 326–340 | none 40 / 43 | ntSH 19 / ntHS 8 + ncWS 1 | Partly near-sheared. A G·A tandem-like geometry is present in about half the chains, but never as a true pair. |
| SAM-I 46–69 | tSH 26 | ntSH 2, none 1 | Sheared in 28 of 29. |
| SAM-I 50–64 | 0 of 25 | none in 25 | Unchanged: alignment-placement problem. |
| TPP 64–83, 21–27 | 0 / 25 | none | Unchanged: no pair. |
| TPP 20–28 | 10 cWW | ncWW 7, ntSH 1 | The E. coli stem is near-paired in some chains. |
| tRNA 30–38 (A31–Ψ39) | 51 cWW | near in 7 of 9 (ncWW 3, ntSH, ncWH, ntHW, ncSW) | The Ψ cases (1EHZ, 1B23, 1TTT, 8UU5) are near pairs, not absent. |
| 16S 438–483 | 4 cWW | ncWW 34 | Mostly near-paired. |
| 16S 667–703 | 19 cWW | ncWW 39 | Mostly near-paired. |
| 16S 75–86 / 76–85 / 77–84 / 78–83 | 23 / 26 / 14 / 12 cWW | ncWW 9 / 6 / 7 / 14 | Partly near-paired; still often no pair. Local variation remains. |
| 16S 372–385 | 0 | none (60) | Unchanged: bulge-position case. |
| 16S 648–737 | 0 | ncWW 1 | Unchanged. |
| 3QSY C | 12 SS_cons pairs without a pair | 6 ncWW, 6 none (none includes C1·A72 and the open 2–69 and 27–41) | As in the follow-up report. |
| 6AHU M | 11 without a pair | 9 ncWW, 1 local true, 1 none | Stems present: near or true pairs. |
| 1J2B C | 5 without a pair (D stem) | none | Genuine rearrangement. |
| 9V50 A / 9NBC A | 7 / 6 without a pair | 6 / 4 ncWW | Mostly near-paired. |

## 6. Conclusions and recommendation

1. **Near base pairs exist and are informative.**
   - They account for about half of the SS_cons observations previously called "unpaired": 1,569 of 3,365.
   - Most of them are ncWW.
   - The status "near only" should be a first-class state in any consensus-structure layer, with its class. It should never be merged into "formed".
2. **The local fr3d-python run reproduces BGSU's true pairs** (identical in 298 of 321 chains). It also fills the 7PJS gap.
3. **Recommended source for the layer:**
   - keep BGSU true pairs (production, unchanged);
   - add near pairs from a pinned local fr3d-python run, with commit and run date in provenance;
   - this replaces the ad-hoc geometry measure.

   Implementing it would add an annotation step to the backend: fr3d-python as an optional dependency, a schema field and tests. The same run could serve the planned near-pair work in the Interactions tab.
4. **Open points:**
   - Re-check the BGSU web interface when it is back up, in case near pairs are published elsewhere.
   - fr3d-python's licence: `license.txt` says MIT but `pyproject.toml` says Apache-2.0, so this is unresolved.
   - Whether 7PJS should use the local true pairs.

## 7. Base–ribose, base–phosphate and stacking contacts (added 2026-10-22)

**Question.** For SS_cons observations with no true base pair, is there a base–ribose (BR), base–phosphate (BPh) or stacking contact between exactly i and j?

**Source.** These contacts come from the production BGSU observations, whose BR and BPh files include near classes. Script: `backbone_check.py`. Outputs: `followup/backbone/`.

**Observations are counted at the strongest available evidence.** The order is: true base pair, then near base pair, then true BR/BPh, then near BR/BPh, then stacking only. BR/BPh and stacking are not base pairs, so these would be separate statuses, never counted as "paired".

| Family | Any true base pair | + near base pair | + true BR/BPh | + near BR/BPh | Still nothing |
|---|---|---|---|---|---|
| RF00167 | 0.953 | 0.987 | 0.987 | 0.987 | 20 |
| RF00162 | 0.932 | 0.960 | 0.960 | 0.960 | 33 |
| RF00059 | 0.887 | 0.908 | 0.908 | 0.908 | 73 |
| RF00005 | 0.920 | 0.973 | 0.973 | 0.973 | 33 (+1 stacking only) |
| RF00011 | 0.898 | 0.938 | 0.941 | 0.948 | 237 (+14 stacking only) |
| RF00177 | 0.906 | 0.951 | 0.956 | 0.957 | 1,095 (+71 stacking only) |

**Coverage barely changes.**
- Base–backbone contacts add 0 for four families, and at most 0.01 for RNase P and 16S.
- The "no pair" cases (SAM-I 50–64, TPP 64–83 and 21–27, 16S 372–385) have no such contact either.

**A few pairs are explained consistently.**

| Pair | Nucleotides | Contact between i and j |
|---|---|---|
| 16S 572–748 | C·C | 6BPh in 56 of 59 chains (near base pair in 3); no base pair in any chain. |
| 16S 596–622 | G·U | 1BPh in 56 of 60 chains, n1BPh in 4; no base pair in any chain. |
| 16S 37–537 | U·A | Near base pair in 36 of 52, n6BR only in 14. |
| RNase P 326–340 | A·G | 6BR in 5 of 43, n6BR only in 16, near base pair in 9, nothing in 13. |
| RNase P 325–341 | G·A | Near base pair in 19 of 40, BR/BPh in 9, nothing in 12. |

For 572–748 and 596–622 the two residues are in contact in every chain, but always as a base–phosphate contact, never as a base pair. These are reproducible within E. coli. Whether the cause is the consensus or the alignment column placement is not established.

## 8. Decision (user, 2026-10-22)

**Only base pairs are considered, with near base pairs included.** Base–ribose, base–phosphate and stacking contacts are left out.

- **Paired** means a true or near FR3D base pair of any class between i and j. Each observation keeps its class and a near flag, so the breakdown (cWW, other classes, near) stays visible.
- **Not paired** means no base pair, true or near.
- **Not evaluable** means a residue is missing or incomplete, or falls in an annotation gap.
- **Sources:** BGSU true pairs, plus near pairs from the pinned local fr3d-python run.

**Support with this definition (any class, true or near):**

| Family | Support |
|---|---|
| RF00167 | 0.987 |
| RF00162 | 0.960 |
| RF00059 | 0.908 |
| RF00005 | 0.973 |
| RF00011 | 0.938 |
| RF00177 | 0.951 |

## 9. Agreed display design (user, 2026-10-22; not built)

**No helix information for now.** The consensus comparison is shown in the Interactions tab: in the table rows and in the arcs of the circular diagram.

**A new column (header "Rfam consensus") with three values:**

| Value | Meaning |
|---|---|
| **Both** | In the Rfam 2D structure (SS_cons) and paired in 3D, as a true or near FR3D base pair of any class, in at least one eligible chain. |
| **2D** | In the Rfam 2D structure; not paired in any eligible chain. These pairs get a row (e.g. SAM-I 50–64, 0 / 25). |
| **3D** | Paired in 3D; not in the Rfam 2D structure. |

**Table and frequency:**
- The tooltip defines "2D" as the Rfam consensus secondary structure, not the PDBe 2D viewer.
- The frequency column (paired / eligible) counts near pairs.
- The ▸ breakdown lists near classes separately.

**Arcs:** the same three states are shown in the circular diagram, which exists only for families with Rfam clen ≤ 300.

**Displayed structure and Compare two:** each pair gets the same values per chain, plus "not evaluable" where the pair can't be checked.

**Build order:**
1. A fr3d-python near base-pair step, pinned, with provenance, schema and tests. The licence question is settled first.
2. The consensus flag and the "2D" rows in `interactions.json`.
3. The frontend column, filter and arcs, with vitest and Playwright tests.
