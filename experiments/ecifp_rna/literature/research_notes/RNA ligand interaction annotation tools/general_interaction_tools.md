# General-purpose, locally runnable ligand–nucleic-acid (RNA) non-covalent interaction annotators

Scope: candidates to replace or complement PDBe-Arpeggio (behind PDBe `bound_ligand_interactions`). Retrieved 2026-10-23. Repository metadata (licence SPDX, last commit, latest release) came from the public GitHub REST API (`https://api.github.com/repos/<owner>/<repo>`, `/releases`, `/commits`) and PyPI JSON (`https://pypi.org/pypi/<pkg>/json`) on 2026-10-23; the per-repo URLs are given inline. Items marked **[unverified]** come from background knowledge and were not checked against a primary source in this session.

## Inventory (one block per tool)

### Snapshot table (verified items only)

| Tool | Latest version (date) | Licence | Last commit | Install | RNA receptor explicit? | Input |
|---|---|---|---|---|---|---|
| PDBe-Arpeggio | 1.4.4 (2024-01-24) | GPL-3.0 | 2024-01-24 | pip `pdbe-arpeggio` + conda openbabel/gemmi/biopython | No (README says protein only); nucleotides get generic SMARTS typing | mmCIF only |
| Arpeggio (original, harryjubb) | no releases | GPL-3.0 | 2026-03-05 (README says not under active development) | Docker / source, Python 2 era | No | PDB |
| PLIP | 3.0.1 (GitHub 2026-07-03; PyPI 2026-07-31) | GPL-2.0 | 2026-07-03 | pip `plip` + OpenBabel 3.2.0, Docker | Yes, `--dnareceptor` | PDB format |
| fingeRNAt | v3.0 (2021-12-22) | GPL-3.0 | 2022-06-03 | conda env / apt / Singularity, source script | Yes (purpose-built for NA) | receptor PDB (with H), ligand SDF |
| ProLIF | 2.2.3 (2026-10-09) | Apache-2.0 | 2026-10-09 | pip/conda `prolif` (RDKit + MDAnalysis) | README says "ligands, protein, DNA or RNA"; no NA tutorial | anything MDAnalysis/RDKit reads; templates from SMILES or CIF |
| LUNA | 0.14.0 (2025-09-18) | MIT | 2025-09-18 | pip `luna` (RDKit + OpenBabel) | Not found in README/source grep | PDB/mol files [details unverified] |
| BINANA | 2.1 (2023-07-03) | Apache-2.0 | 2025-06-30 | source (Python, also JavaScript transpile) | No (protein residue/atom-name tables) | PDBQT (preferred) or PDB |
| ODDT | 0.7 (2019-12-24 release; PyPI 2019-12-23) | BSD-3-Clause | 2026-10-04 | pip/conda `oddt` | Generic (toolkit-perceived atom features) | via OpenBabel/RDKit |
| GetContacts | no releases | Apache-2.0 | 2023-03-22 (push 2026-06-28) | source + vmd-python | Not stated | any VMD-readable topology/trajectory |
| PyContact | v1.0.2 (2017-11-03) | GPL-3.0 | 2023-03-12 | source (MDAnalysis) | Not stated | MD trajectories |
| MDAnalysis (HydrogenBondAnalysis) | 2.11.0 (2026-10-21) | NOASSERTION on GitHub (LGPL/GPL mix [unverified]) | 2026-10-21 | pip/conda | Generic selection-based | many formats incl. PDB/mmCIF [unverified for mmCIF] |
| Mol* (computed interactions) | v5.13.0 (2026-10-04) | MIT | 2026-10-23 | npm (TypeScript) | Yes, nucleotide phosphates typed as negative charge | mmCIF/BinaryCIF native |

Sources for the table: [PDBeurope/arpeggio](https://github.com/PDBeurope/arpeggio), [PyPI pdbe-arpeggio](https://pypi.org/project/pdbe-arpeggio/), [harryjubb/arpeggio](https://github.com/harryjubb/arpeggio), [pharmai/plip](https://github.com/pharmai/plip), [PyPI plip](https://pypi.org/project/plip/), [n-szulc/fingeRNAt](https://github.com/n-szulc/fingeRNAt), [chemosim-lab/ProLIF](https://github.com/chemosim-lab/ProLIF), [keiserlab/LUNA](https://github.com/keiserlab/LUNA), [durrantlab/binana](https://github.com/durrantlab/binana), [oddt/oddt](https://github.com/oddt/oddt), [getcontacts/getcontacts](https://github.com/getcontacts/getcontacts), [maxscheurer/pycontact](https://github.com/maxscheurer/pycontact), [MDAnalysis/mdanalysis](https://github.com/MDAnalysis/mdanalysis), [molstar/molstar](https://github.com/molstar/molstar).

## Which tools explicitly support RNA as the receptor and with what atom typing?

### Takeaway
Three tools support a nucleic-acid receptor explicitly: PLIP (`--dnareceptor`, OpenBabel perception plus a hard-coded A/C/G/U name list), fingeRNAt (built for nucleic acids; OpenBabel/RDKit perception plus a rule that OP1/OP2 are always anions), and Mol* (NGL-derived feature typing; nucleotide phosphates are negative charges). ProLIF claims DNA/RNA support in its README but documents none. PDBe-Arpeggio types nucleotides only through generic OpenBabel SMARTS on geometry-perceived bonds. Its curated residue-name override table covers standard amino acids only.

### Cited Findings
**PDBe-Arpeggio**
- The README says it calculates "interatomic contacts in a protein" and that "This library only supports protein structures in mmCIF format". — [PDBe-Arpeggio README](https://github.com/PDBeurope/arpeggio/blob/master/README.md)
- Atom typing: SMARTS-based `ATOM_TYPES` are matched on the OpenBabel molecule. A curated `PROT_ATOM_TYPES` table then overrides the types, but only for residues in `config.STD_RES` (the 20 amino acids). The loop is commented "OVERRIDE PROTEIN ATOM TYPING FROM DICTIONARY". — [arpeggio/core/interactions.py](https://github.com/PDBeurope/arpeggio/blob/master/arpeggio/core/interactions.py) (around line 1965), [arpeggio/core/config.py](https://github.com/PDBeurope/arpeggio/blob/master/arpeggio/core/config.py)
- `config.py` defines `STANDARD_NUCLEOTIDES = set(['A','C','G','I','U','DA','DC','DG','DI','DT','DU','N'])` and maps `_chem_comp.type` values such as 'RNA LINKING' to a `ComponentType` (R = Polyribonucleotide). No nucleotide atom-type dictionary was found. — [config.py](https://github.com/PDBeurope/arpeggio/blob/master/arpeggio/core/config.py)
- Bonds come from `OBResidueData.AssignBonds`, then `ConnectTheDots()` and `PerceiveBondOrders()` (geometry-based), with extra bonds from `_struct_conn`. CCD bond orders are not used. — [arpeggio/core/protein_reader.py](https://github.com/PDBeurope/arpeggio/blob/master/arpeggio/core/protein_reader.py)
- Arpeggio's README says it "may be applied to DNA macromolecules". This is the fingeRNAt authors' characterisation. — [fingeRNAt paper, PLoS Comput Biol 2022](https://doi.org/10.1371/journal.pcbi.1009783)

**PLIP**
- `--dnareceptor`: "Treat nucleic acids as part of the receptor structure (together with any present protein) instead of as a ligand." — [plip/plipcmd.py](https://github.com/pharmai/plip/blob/master/plip/plipcmd.py)
- From the documentation: "A special mode allows to switch from treating DNA/RNA molecules as ligands to treating them as part of the receptor … start PLIP with the option `--dnareceptor`." — [PLIP DOCUMENTATION.md](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md)
- Config hard-codes `RNA = ['U','A','C','G']` and `DNA = ['DT','DA','DC','DG']`, so modified nucleotide residue names are not in these lists. — [plip/basic/config.py](https://github.com/pharmai/plip/blob/master/plip/basic/config.py)
- Donors and acceptors come from OpenBabel. Aromatic rings come from OpenBabel SSSR, with a planarity fallback (`AROMATIC_PLANARITY = 5.0`°). Charged groups are defined for protein residues (Arg, His, Lys, Asp, Glu) and for ligand groups. The docs say "The detection of charged groups is only exhaustive for the binding site, not the ligands". No nucleotide phosphate charge rule is documented. — [DOCUMENTATION.md](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md)
- PLIP 2021 paper: "the PLIP web tool has been adapted to recognise nucleic acid molecules as receptors". — [Adasme et al., NAR 2021, gkab294](https://doi.org/10.1093/nar/gkab294) (full text via [Europe PMC PMC8262720](https://europepmc.org/article/PMC/PMC8262720))
- Changelog 2.1.9: "fixes issues with DNA receptor handling". — [CHANGES.txt](https://github.com/pharmai/plip/blob/master/CHANGES.txt)
- PLIP 2025 adds protein–protein interactions and keeps "eight types of non-covalent interactions". It reports no new nucleic-acid features. — [Schake et al., NAR 2025, gkaf361](https://doi.org/10.1093/nar/gkaf361)

**fingeRNAt**
- "a software tool for detecting non-covalent interactions formed within complexes of nucleic acids with ligands". — [fingeRNAt README](https://github.com/n-szulc/fingeRNAt/blob/master/README.md)
- "fingeRNAt always treats OP1 and OP2 atoms as negatively charged anions". Receptor residues with fewer than 4 atoms are ignored. — [README](https://github.com/n-szulc/fingeRNAt/blob/master/README.md)
- Aromatic rings are detected with OpenBabel (receptor) and RDKit (ligand). Detection is "Inspired by PLIP implementation". Users can add their own SMARTS-defined interactions through YAML. — [implementation.md](https://github.com/n-szulc/fingeRNAt/blob/master/implementation.md)
- The `PBS` fingerprint splits each nucleotide into Phosphate, Base and Sugar groups using a 4.0 Å distance threshold. — [README](https://github.com/n-szulc/fingeRNAt/blob/master/README.md)

**ProLIF**
- README: "generate interaction fingerprints for complexes made of ligands, protein, DNA or RNA molecules extracted from molecular dynamics trajectories, docking simulations and experimental structures". — [ProLIF README](https://github.com/chemosim-lab/ProLIF/blob/master/README.rst)
- Changelog 0.3.1: "`ResidueId.from_string` method now supports 1-letter and 2-letter codes for RNA/DNA". — [CHANGELOG.md](https://github.com/chemosim-lab/ProLIF/blob/master/CHANGELOG.md)
- The PDB tutorial does not mention nucleic acids. Non-standard residues need SMILES or CIF templates; standard (protein) residue states are covered by default. — [ProLIF PDB tutorial](https://prolif.readthedocs.io/en/stable/notebooks/pdb.html)
- Interactions are SMARTS on RDKit molecules. Classes: Anionic, Cationic, CationPi, PiCation, EdgeToFace, FaceToFace, PiStacking, HBAcceptor, HBDonor, ImplicitHBAcceptor, ImplicitHBDonor, Hydrophobic, MetalAcceptor, MetalDonor, VdWContact, XBAcceptor, XBDonor. — [ProLIF PDB tutorial](https://prolif.readthedocs.io/en/stable/notebooks/pdb.html)

**Mol\***
- `charged.ts` (comment "based in part on NGL") has a branch for `BaseNames` components: when `isPhosphate(...)` is true, it records `FeatureType.NegativeCharge, FeatureGroup.Phosphate`. — [packages/model/src/props/computed/interactions/charged.ts](https://github.com/molstar/molstar/blob/master/packages/model/src/props/computed/interactions/charged.ts)
- Interaction types (enum): Ionic, CationPi, PiStacking, HydrogenBond, HalogenBond, Hydrophobic, MetalCoordination, WeakHydrogenBond, WaterBridge. — [common.ts](https://github.com/molstar/molstar/blob/master/packages/model/src/props/computed/interactions/common.ts)

**BINANA 2**
- "In the case of protein receptors, BINANA 2 now identifies hydrogen bond donors/acceptors and charged groups based on the residue and atom names (assuming standardized PDB nomenclature). In the case of the ligand, it makes rudimentary assumptions about protonation based on the geometry of the heavy atoms." The paper mentions no nucleic-acid receptor. — [Durrant et al., JCIM 2022 "BINANA 2"](https://doi.org/10.1021/acs.jcim.1c01461) (full text via [Europe PMC PMC8889568](https://europepmc.org/article/PMC/PMC8889568))

**LUNA**
- Toolkit "to analyze very large data sets of 3D molecular structures and complexes" with EIFP/FIFP/HIFP fingerprints. The README does not mention nucleic acids. — [LUNA README.rst](https://github.com/keiserlab/LUNA/blob/master/README.rst)
- Feature typing uses RDKit `.fdef` feature definitions (`luna/data/LUNA.fdef`, `BaseFeatures.fdef`). — [LUNA repo tree](https://github.com/keiserlab/LUNA/tree/master/luna/data)

**ODDT**
- The `oddt.interactions` functions (hbonds, halogenbonds, pi_stacking, salt_bridges, hydrophobic_contacts, pi_cation, acceptor_metal, pi_metal) work on toolkit-perceived atom dictionaries and are not receptor-type specific. — [oddt/interactions.py](https://github.com/oddt/oddt/blob/master/oddt/interactions.py)

**GetContacts**
- The hydrogen-bond subtypes are protein-specific: backbone/sidechain (hbbb, hbsb, hbss) and ligand–backbone/sidechain (hblb, hbls). The README does not mention nucleic acids. — [GetContacts README](https://github.com/getcontacts/getcontacts/blob/master/README.md)

### Inferences
- For RNA receptors, the explicit nucleic-acid rules in these tools are thin. They amount to phosphate O1P/O2P set as anions (fingeRNAt, Mol*) plus generic aromatic-ring and donor/acceptor perception. None ships a curated per-atom nucleotide table comparable to Arpeggio's `PROT_ATOM_TYPES`. A custom table for A/C/G/U atoms (donor/acceptor and edge) would be needed under any tool for reproducible typing.
- Arpeggio's typing of RNA atoms depends on OpenBabel bond-order perception from geometry. That is fragile for modified nucleotides and for poorly resolved bases, so contacts from PDBe's API for RNA should be treated as "generic SMARTS-typed".

### Gaps
- I did not read PLIP's code path that types nucleotide atoms in `--dnareceptor` mode (for example, whether phosphate O are set as negative charges for salt bridges). Check `plip/structure/preparation.py`.
- Could not confirm whether LUNA handles nucleic-acid residues (grep for "nucleic|dna|rna" in the downloaded core files found nothing). The GitHub code-search API requires authentication.
- The fingeRNAt S1 Table compares fingeRNAt, Arpeggio, PLIP2021 and ProLIF but is a supplementary PDF that was not retrieved. — [fingeRNAt paper](https://doi.org/10.1371/journal.pcbi.1009783)

## Which detect water-mediated and metal-ion-mediated contacts, with what geometric criteria?

### Takeaway
Water bridges: PLIP (one water, distance plus two angles), fingeRNAt (distance only, ≤ 3.5 Å each side, `-h2o`), LUNA, GetContacts (wb, wb2 and ligand-specific lwb, lwb2), ProLIF (WaterBridge since 2.1.0) and Mol* (WaterBridge type). Arpeggio reports only atom–water contacts (SELECTION_WATER), not bridges. Ion-mediated ligand–nucleic acid bridges with ion-specific cutoffs exist only in fingeRNAt (Mg 3.2, K 3.9, Na 3.6, other 3.5 Å). PLIP, PDBe-Arpeggio, BINANA, LUNA, ODDT and Mol* detect metal coordination by the ligand or receptor, not metal-mediated bridges.

### Cited Findings
- **PLIP water bridges**: one water only. Water O within `WATER_BRIDGE_MINDIST` 2.5 Å and `WATER_BRIDGE_MAXDIST` 4.1 Å of the polar atoms. ω (acceptor–O–donor H) 71–140°, θ (O–H–donor) ≥ 100°. At most two donor H-bonds per water (keeps the pair closest to 110°). Hydrogens are needed (donor-H angle). — [config.py](https://github.com/pharmai/plip/blob/master/plip/basic/config.py); [DOCUMENTATION.md](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md)
- **PLIP metal complexes**: more than 50 metal species; `METAL_DIST_MAX = 3.0` Å (Harding 2001). The geometry is fitted to ideal coordination geometries (linear to octahedral). The documented receptor target groups are protein side chains and main-chain O; ligand groups are alcohols, phenolates, carboxylates, phosphoryls, thiolates, imidazoles and pyrroles. — [DOCUMENTATION.md](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md)
- **fingeRNAt ion-mediated**: ligand N/O/S–ion ≤ X and nucleic acid N/O–ion ≤ X, where X = 3.2 Å Mg, 3.9 Å K, 3.6 Å Na, 3.5 Å other (Zheng et al. 2008/2015). Reported separately as Mg_mediated, K_mediated, Na_mediated and Other_mediated. — [implementation.md](https://github.com/n-szulc/fingeRNAt/blob/master/implementation.md); [README](https://github.com/n-szulc/fingeRNAt/blob/master/README.md)
- **fingeRNAt water-mediated**: ligand donor/acceptor–water O ≤ 3.5 Å and nucleic acid donor/acceptor–water O ≤ 3.5 Å. Only computed when the flag is passed. The README text names it both `-h2o` and `-hoh`, which is inconsistent. — [implementation.md](https://github.com/n-szulc/fingeRNAt/blob/master/implementation.md); [README](https://github.com/n-szulc/fingeRNAt/blob/master/README.md)
- **fingeRNAt ions as ligands**: without `-l`, all inorganic ions in the receptor file are treated as ligands, and only SIMPLE or PBS (contact) fingerprints are allowed. — [README](https://github.com/n-szulc/fingeRNAt/blob/master/README.md)
- **PDBe-Arpeggio**: `metal` ("part of a metal complex") is an atom–atom contact type. Water contacts are classified as SELECTION_WATER, NON_SELECTION_WATER or WATER_WATER; all water atoms are typed as H-bond donor and acceptor. — [README](https://github.com/PDBeurope/arpeggio/blob/master/README.md); [interactions.py](https://github.com/PDBeurope/arpeggio/blob/master/arpeggio/core/interactions.py)
- **LUNA**: interaction types include "Water-bridged hydrogen bond" and "Metal coordination" (`METAL_COMPLEX_DIST = 2.8`). — [luna/util/default_values.py](https://github.com/keiserlab/LUNA/blob/master/luna/util/default_values.py)
- **ProLIF**: "Support for `WaterBridge` interactions (PR #229)" in 2.1.0 (2026-01-20). MetalAcceptor and MetalDonor classes exist. — [CHANGELOG.md](https://github.com/chemosim-lab/ProLIF/blob/master/CHANGELOG.md); [PDB tutorial](https://prolif.readthedocs.io/en/stable/notebooks/pdb.html)
- **GetContacts**: wb, wb2 (extended, two waters), lwb and lwb2 (ligand water-mediated, extended), with extra columns naming the waters. — [README](https://github.com/getcontacts/getcontacts/blob/master/README.md)
- **BINANA 2**: "now also detects metal-coordination bonds whenever select metal cations come within a user-defined distance of nitrogen, oxygen, sulfur, iodine, bromine, chlorine, and fluorine atoms … defaults to 3.5 Å". — [BINANA 2 paper](https://doi.org/10.1021/acs.jcim.1c01461)
- **ODDT**: `acceptor_metal(cutoff=4, tolerance=30)` and `pi_metal(cutoff=5)`. — [oddt/interactions.py](https://github.com/oddt/oddt/blob/master/oddt/interactions.py)
- **Mol\***: MetalCoordination and WaterBridge interaction types exist. The H-bond code has an `includeWater` option. — [common.ts](https://github.com/molstar/molstar/blob/master/packages/model/src/props/computed/interactions/common.ts), [hydrogen-bonds.ts](https://github.com/molstar/molstar/blob/master/packages/model/src/props/computed/interactions/hydrogen-bonds.ts)

### Inferences
- Mg²⁺-mediated ligand–RNA bridges are common in riboswitches (for example TPP and fluoride). fingeRNAt is the only tool found with ion-specific bridge cutoffs and per-ion classes, but it computes them per residue in a fingerprint, not as atom triplets. Atom-level triplets can be rebuilt from its `-detail` output (not checked) or with a short custom script on Arpeggio or gemmi neighbour output.

### Gaps
- Mol* WaterBridge criteria (distance and angle) were not read; `water.ts` does not exist at the path tried.
- ProLIF WaterBridge parameters (maximum order, distances) were not read.

## Which handle mmCIF directly (large entries, e.g. ribosomes), altlocs and modified nucleotides?

### Takeaway
Among the open, locally runnable tools, only PDBe-Arpeggio (gemmi/Biopython reader) and Mol* read mmCIF natively. PLIP and fingeRNAt are PDB-format tools (fingeRNAt receptor PDB, ligand SDF), and BINANA takes PDBQT or PDB. ProLIF can read via MDAnalysis, and it takes CIF templates for ligand chemistry. Altloc handling varies: Arpeggio keeps every altloc atom; PLIP drops altlocs unless `--altlocation`; Arpeggio and PLIP both take only the first model by default. No tool has explicit modified-nucleotide typing. All rely on generic perception, and PLIP's RNA list covers only A/C/G/U.

### Cited Findings
- PDBe-Arpeggio refactor changes include "Support for mmCIF format to process also larger structures" and "Results in JSON format". The README warns: "Alternate locations and missing density are not explicitly accounted for and may result in anomalous results." — [README](https://github.com/PDBeurope/arpeggio/blob/master/README.md)
- PDBe-Arpeggio reader: "Trim all the other models but the first one from the mmcif structure". Atom identifiers use `pdbe_label_seq_id` if present, otherwise `auth_seq_id`. — [protein_reader.py](https://github.com/PDBeurope/arpeggio/blob/master/arpeggio/core/protein_reader.py)
- PDBe-Arpeggio keeps all altloc children (`for altloc, alt_atom in atom.child_dict.items(): self.s_atoms.append(alt_atom)`). — [interactions.py](https://github.com/PDBeurope/arpeggio/blob/master/arpeggio/core/interactions.py)
- PLIP: "Polar hydrogens are added to the structure and alternative conformations/models/positions removed". Options: `--altlocation` ("Do not discard alternate locations") and `--model`. By default it uses the first model. — [DOCUMENTATION.md](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md); [README](https://github.com/pharmai/plip/blob/master/README.md)
- PLIP: "Modified amino acids are identified and excluded using MODRES entries of the PDB files." BioLiP artefact ligands appearing 15 or more times are removed. — [DOCUMENTATION.md](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md)
- PLIP describes itself as analysing "protein-ligand interactions in PDB files". The preparation code parses PDB lines and fixes PDBQT atom types. — [DOCUMENTATION.md](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md); [preparation.py](https://github.com/pharmai/plip/blob/master/plip/structure/preparation.py)
- PLIP 2.2.0 changelog: "improved handling of multi-model structures". 2.1.5: "fixes a bug in alt-location handling". 3.0.0: "fixed handling of modified residues in protein and peptide ligands". — [CHANGES.txt](https://github.com/pharmai/plip/blob/master/CHANGES.txt)
- fingeRNAt: "currently operates on strictly defined file formats (pdb for nucleic acid receptors and sdf for ligands)". "The numbering of residues in the receptor structure must be unique". — [fingeRNAt paper, Known limitations](https://doi.org/10.1371/journal.pcbi.1009783)
- BINANA: "accepts receptor and ligand files in the PDBQT (preferred) or PDB formats". — [BINANA README](https://github.com/durrantlab/binana/blob/master/README.md)
- ProLIF 2.2.0: `MoleculeStandardizer` "to parse and correct the connectivity of your topology easily with SMILES or CIF templates". 2.2.3 added `FragmentedResidueError` "for when multiple residues share the same ResidueId". — [CHANGELOG.md](https://github.com/chemosim-lab/ProLIF/blob/master/CHANGELOG.md)
- GetContacts is "compatible with all topology and reimaged trajectory file formats readable by VMD". — [README](https://github.com/getcontacts/getcontacts/blob/master/README.md)
- Background: legacy PDB format cannot represent very large entries; PDBx/mmCIF can. — [PDBe "Improved representation of large structures"](https://www.ebi.ac.uk/pdbe/node/1270)

### Inferences
- For ribosome-scale entries (for example the RF00177 16S set), PLIP, fingeRNAt and BINANA would need a PDB-format extract around each ligand (for example a gemmi-cropped sphere), with residue renumbering for fingeRNAt's uniqueness rule. PDBe-Arpeggio can take the full mmCIF but loads the whole structure into OpenBabel, so cropping is still advisable for speed (speed not measured).
- Arpeggio keeping all altlocs can produce duplicate or contradictory contacts for disordered ligands or nucleotides. That is consistent with the README's "anomalous results" warning.

### Gaps
- No PLIP issue or documentation confirming or denying mmCIF input was found. The search returned nothing specific.
- MDAnalysis mmCIF reader status was not verified.
- Speed and scalability numbers were not found for any tool.

## Which produce atom-level, machine-readable output suitable for mapping to residue edges (Watson–Crick/Hoogsteen/sugar)?

### Takeaway
Atom-level, machine-readable output: PDBe-Arpeggio (JSON, with auth_asym_id, auth_seq_id, auth_atom_id, label_comp_id, contact types and distance), PLIP (XML with donor and acceptor atom indices, types, coordinates and angles), GetContacts (TSV with chain:resname:resid:atom per partner), BINANA 2 (JSON/CSV), Mol* (in-memory feature and atom indices; no file exporter documented), and ChimeraX `hbonds saveFile` (text with atom specs and distances). fingeRNAt's main output is per-residue fingerprints (TSV); `-detail` gives a separate per-interaction file. Edge assignment (WC/Hoogsteen/sugar) is not done by any tool. It needs a nucleotide atom→edge lookup on top of the atom names.

### Cited Findings
- The PDBe-Arpeggio JSON example has `bgn` and `end` objects with `auth_asym_id`, `auth_atom_id`, `auth_seq_id`, `label_comp_id`, `label_comp_type`, `pdbx_PDB_ins_code`, plus `contact` (a list such as `["proximal","hydrophobic"]`), `distance` and `interacting_entities`. There are also atom–plane, plane–plane and group–group records. — [README](https://github.com/PDBeurope/arpeggio/blob/master/README.md)
- PDBe-Arpeggio contact types: clash, covalent, vdw_clash, vdw, proximal, hbond, weak_hbond, xbond, ionic, metal, aromatic, hydrophobic, carbonyl, polar, weak_polar. Atom–plane: CARBONPI, CATIONPI, DONORPI, HALOGENPI, METSULPHURPI. Plane–plane follows Chakrabarti & Bhattacharyya 2007. Group: AMIDEAMIDE, AMIDERING. — [README](https://github.com/PDBeurope/arpeggio/blob/master/README.md)
- PLIP XML gives donor and acceptor atom IDs and types, ligand and protein coordinates, water_idx, metal_idx, target_idx and coordination number. `ligtype` can be SMALLMOLECULE, POLYMER, DNA, RNA or ION. Changelog 2.2.0: "increased detail level of report (individual atoms on protein side)". — [DOCUMENTATION.md](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md); [CHANGES.txt](https://github.com/pharmai/plip/blob/master/CHANGES.txt)
- PLIP types (8): hydrophobic, H-bond, water bridge, salt bridge, π-stacking (parallel and perpendicular), π-cation, halogen bond, metal complex. Thresholds: HBOND_DIST_MAX 4.1 Å, HBOND_DON_ANGLE_MIN 100°, PISTACK_DIST_MAX 5.5 Å, PISTACK_ANG_DEV 30°, PISTACK_OFFSET_MAX 2.0 Å, PICATION_DIST_MAX 6.0 Å, SALTBRIDGE_DIST_MAX 5.5 Å, HALOGEN_DIST_MAX 4.0 Å, HYDROPH_DIST_MAX 4.0 Å. — [config.py](https://github.com/pharmai/plip/blob/master/plip/basic/config.py); [PLIP 2025](https://doi.org/10.1093/nar/gkaf361)
- fingeRNAt FULL fingerprint covers 9 types per residue: HB, HAL, cation–anion, π–cation, π–anion (both directions, including RNA phosphate to ligand ring), π-stacking (sandwich, parallel displaced, T-shaped; centroid < 5.5 Å, offset < 2.0 Å), ion-mediated, water-mediated and lipophilic. Outputs are TSV. `-detail` writes per-interaction data for the PyMOL plugin. — [README](https://github.com/n-szulc/fingeRNAt/blob/master/README.md); [implementation.md](https://github.com/n-szulc/fingeRNAt/blob/master/implementation.md)
- fingeRNAt default H-bond is distance-only (D···A ≤ 3.9 Å); `-dha` adds a D-H-A angle. "The number of detected interactions for methods that consider the position of hydrogen is almost three times lower (ranging from 1735 to 1822) than for the default method (5026)" on 207 RNA–ligand complexes, and depends on the H-placement program. — [fingeRNAt paper](https://doi.org/10.1371/journal.pcbi.1009783)
- In fingeRNAt's RNA–ligand dataset, nucleobases form 61% of H-bonds, phosphate 23% and ribose 15%. — [fingeRNAt paper](https://doi.org/10.1371/journal.pcbi.1009783)
- GetContacts TSV rows look like `7 lwb2 A:EJ4:1219:NAN A:LYS:214:NZ A:HOH:1316:O A:HOH:1315:O`. Stacking and π-cation are "denoted by the CG atom", a protein-specific simplification that loses ring identity for nucleobases. — [README](https://github.com/getcontacts/getcontacts/blob/master/README.md)
- BINANA 2 "added JSON- and CSV-formatted output" and runs in Python and transpiled JavaScript. — [BINANA 2 paper](https://doi.org/10.1021/acs.jcim.1c01461)
- ChimeraX `hbonds saveFile` writes donor and acceptor atom specs and D–A distances, plus H specs and H–A distances when hydrogens are present. Explicit H are not required (inferred from atom type). Mills & Dean 1996 criteria; `relax` defaults to false (0.4 Å and 20° slop when enabled); `saltOnly` option. — [ChimeraX hbonds docs](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/hbonds.html)
- The Mol* Interactions extension "enables computing or providing custom interactions between multiple selections/structures". The example readme describes computing on demand, per-ligand computation from an SDF, and manual custom interactions. It is run through the dev server (`npm run dev -- -e interactions`). No headless CLI is documented. — [Mol* interactions extension doc](https://github.com/molstar/molstar/blob/master/apps/docs/src/content/docs/v5/extensions/interactions.md); [examples/interactions readme](https://github.com/molstar/molstar/blob/master/examples/interactions/src/readme.md)
- The Mol* `cli/` folder contains utilities such as cif2bcif, cifschema and chem-comp-dict. No interactions CLI was found in the repo tree listing. — [molstar repo tree](https://github.com/molstar/molstar/tree/master/cli)

### Inferences
- PDBe-Arpeggio JSON (or PLIP XML) already carries `auth_atom_id`. Mapping to edges needs only a static table (for example A: N6/N1 WC, N7/N6 Hoogsteen, O2'/N3 sugar; G: N1/N2/O6 WC, O6/N7 Hoogsteen, N2/N3/O2' sugar), the same convention as FR3D/Leontis–Westhof. Plane-based records (π-stacking, cation–π) map to "base face" rather than an edge.
- Running Mol* interactions headless is technically plausible in Node, since the model package is TypeScript and separate from graphics. It is not a supported or documented route, and results would need custom serialisation (not verified).

### Gaps
- Exact fingeRNAt `-detail` columns (whether they include atom names) were not read.
- Mol* H-bond geometric defaults were not extracted (parameter definitions not located in the grep).

## Licences compatible with a public academic web resource (EMBL-EBI/PDBe context)

### Takeaway
Permissive (MIT, Apache-2.0, BSD-3): Mol*, LUNA, ProLIF, BINANA, GetContacts and ODDT. Copyleft GPL: PDBe-Arpeggio and Arpeggio (GPL-3.0), PLIP (GPL-2.0), fingeRNAt (GPL-3.0) and PyContact (GPL-3.0). Running GPL tools server-side to produce data is generally fine; redistributing modified code would trigger GPL obligations. PDBe already runs GPL-3 Arpeggio. LigPlot+, HBPLUS and IChem have non-OSI or unclear terms. GetContacts depends on vmd-python, whose licence was not checked.

### Cited Findings
- Licence SPDX values from the GitHub API, 2026-10-23: PDBe-Arpeggio GPL-3.0; Arpeggio GPL-3.0; PLIP GPL-2.0; fingeRNAt GPL-3.0; ProLIF Apache-2.0; LUNA MIT; BINANA Apache-2.0; ODDT BSD-3-Clause; GetContacts Apache-2.0; PyContact GPL-3.0; Mol* MIT; MDAnalysis NOASSERTION. — repo URLs in the snapshot table above.
- PLIP README: "PLIP is published under the GNU GPLv2 … Using PLIP in your commercial or non-commercial project is generally possible when giving a proper reference". — [PLIP README](https://github.com/pharmai/plip/blob/master/README.md)
- The ProLIF README links the docs under a "Creative Commons Attribution 4.0 International License" (documentation), while the repository licence is Apache-2.0. — [ProLIF README](https://github.com/chemosim-lab/ProLIF/blob/master/README.rst)
- The LigPlot+ page states "LigPlot+ v.2.3 will shortly be available from UCLB's XIP website". It gives no licence terms and describes the tool as "a Java application for the automatic generation and onscreen editing of 2D ligand-protein interaction diagrams". — [LigPlot+ page](https://www.ebi.ac.uk/thornton-srv/software/LigPlus/)
- The fingeRNAt authors note that LigPlot/LigPlot+ "display only hydrogen bonds and lipophilic interactions" for nucleic acid–ligand complexes. — [fingeRNAt paper](https://doi.org/10.1371/journal.pcbi.1009783)
- IChem is listed as free with a blank licence field (directory listing, secondary source). It is a protein–ligand toolkit (C++). — [bioinformaticshome IChem](https://bioinformaticshome.com/db/tool/IChem); [IChem paper page](https://lab.drugdesign.unistra.fr/publications/2018_da-silva_1/)

### Inferences
- If PDBe/EMBL-EBI wanted to swap or extend the contact engine, the closest drop-in with the same JSON contract is PDBe-Arpeggio itself plus a nucleotide atom-type dictionary (GPL-3, already in production). A permissive alternative for the web front end is Mol* (MIT), which already runs in PDBe/RCSB viewers but is not packaged for batch use.

### Gaps
- HBPLUS licence (historically a free academic licence by request from UCL) is **[unverified]**.
- MDAnalysis licence (GitHub reports NOASSERTION; historically GPL-2+ moving to LGPL) is **[unverified]**.
- The vmd-python licence (affects GetContacts) is **[unverified]**.

## Other tools named in the brief (not verified this session)

### Takeaway
These were not verified against primary sources here and are listed for completeness. Treat every claim as **[unverified]**.

### Cited Findings
- OpenStructure 2.11 `mol.alg` lists `ligand_scoring`, `contact_score` and `scoring` submodules but no general ligand interaction detector. The page banner says 2.12 is the latest. — [OpenStructure mol.alg docs](https://openstructure.org/docs/2.11/mol/alg/molalg/)
- PLIP 2025 is "available as a web server, source code with containers, and Jupyter notebook". — [PLIP 2025](https://doi.org/10.1093/nar/gkaf361)

### Inferences
- **[unverified]** MDAnalysis `HydrogenBondAnalysis`: selection-based D–H···A geometry that needs explicit hydrogens and bonds (or guessing). Generic, so usable for RNA, but H-bonds only.
- **[unverified]** NCIPLOT: electron-density/RDG-based NCI surfaces. Not a contact-list annotator, and expensive for large systems.
- **[unverified]** PoseView/ProteinsPlus: 2D diagrams on a web server, protein-focused. Not locally runnable for free.
- **[unverified]** Schrödinger (Maestro interaction fingerprints / SIFt), MOE (PLIF) and BIOVIA Discovery Studio: commercial, protein-centred, not suitable for a public open pipeline.
- **[unverified]** IChem IFP: protein–ligand cavity-based typing (Rognan lab); per the PLIP 2021 text, binary fingerprint tools "do not cover DNA and RNA" (refs 35, 36 in that paper; identities not resolved). — [PLIP 2021](https://doi.org/10.1093/nar/gkab294)

### Gaps
- Not checked this session: OpenStructure's `ligand_scoring` contact definitions; NCIPLOT, ChimeraX `contacts` and IChem versions and licences; HBPLUS/DIMPLOT distribution and nucleic-acid support; PoseView availability; BINANA's nucleic-acid behaviour beyond the paper's protein-only typing statement.
- No head-to-head benchmark of these tools on RNA–ligand complexes was found, apart from the fingeRNAt S1 feature table (not retrieved).
