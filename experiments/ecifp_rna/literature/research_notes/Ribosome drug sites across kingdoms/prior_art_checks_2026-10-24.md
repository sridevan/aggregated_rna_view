# Prior-art checks before the ribosome study (2026-10-24, done directly)

| Resource | Ligand or drug-site layer? | Shared positions across kingdoms? | Source checked |
|---|---|---|---|
| RiboVision2 | No. Its layers are nucleotide frequency, Shannon entropy, TwinCons, helices and expansion segments, RNA–protein contacts (3.5 Å), modifications from the CIF, and custom data. | Alignment of 152 DESIRE species; mitochondria not mentioned | https://apollo2.chemistry.gatech.edu/AboutRiboVision2/ribovision_data |
| BGSU RNA 3D Hub | No. It offers the Structure Atlas, representative sets, Motif Atlas, WebFR3D, JAR3D, R3DMCS, R3DSVS and pair, triple and modification catalogues. | No | http://rna.bgsu.edu/rna3dhub/ |
| SPIRAL (Lu & Wang, bioRxiv 2026, PMC13228370) | Yes (DSSR, 4.8 Å). **It excludes any ligand that also contacts protein.** RNAs are grouped into 6 categories by keyword, with an Rfam cross-check only. | No: no shared numbering or alignment; no eukaryotic, mitochondrial, Motif Atlas or FR3D terms in the full text | Europe PMC full-text XML |
| RiboXYZ | Ligand catalogue; 5 Å sites transferred between two chains with pairwise2.globalxx | No: names only | paper and GitHub code (see memory and the report) |
| Melnikov, bacterial (Ekemezie et al., Cell Rep 2025, 115878) | Abstract only: sequence-based evolution of drug-binding residues, plus intrinsic resistance | Not verified: full text blocked (bioRxiv 429, Cell 403, Newcastle ePrints 500). The preprint is 10.1101/2024.10.31.621234 | Europe PMC abstract |

**Conclusion:** none of these covers the gap identified in the report: every copy on one rRNA axis, compared across kingdoms (mitochondria included) from structures, and linked to the Motif Atlas.
