"""Build reference_contacts.json from hand-curated literature contacts (curation by reading papers; no coordinate analysis)."""
import json, sys

OUT = sys.argv[1]

SRC = {
    "1Y26": "Serganov A, Yuan YR, Pikovskaya O, Polonskaia A, Malinina L, Phan AT, Hobartner C, Micura R, Breaker RR, Patel DJ (2004) Chem Biol 11:1729-1741, doi:10.1016/j.chembiol.2004.11.018 (PMC4692365)",
    "4FE5": "Batey RT, Gilbert SD, Montange RK (2004) Nature 432:411-415, doi:10.1038/nature03037 (primary citation of 4FE5; re-refined model deposited with Stoddard CD et al. (2013) J Mol Biol, doi:10.1016/j.jmb.2013.02.023)",
    "4FE5_sec": "Batey RT (2012) Q Rev Biophys 45:345-381 'Structure and mechanism of purine-binding riboswitches', doi:10.1017/S0033583512000078 (PMC3760793) [secondary, same lab]",
    "2GIS": "Montange RK, Batey RT (2006) Nature 441:1172-1175, doi:10.1038/nature04819",
    "2QWY": "Gilbert SD, Rambo RP, Van Tyne D, Batey RT (2008) Nat Struct Mol Biol 15:177-182, doi:10.1038/nsmb.1371",
    "3E5C": "Lu C, Smith AM, Fuchs RT, Ding F, Rajashankar K, Henkin TM, Ke A (2008) Nat Struct Mol Biol 15:1076-1083, doi:10.1038/nsmb.1494 (PMC3467307)",
    "6FZ0": "Huang L, Lilley DMJ (2018) Nucleic Acids Res 46:6869-6879, doi:10.1093/nar/gky520 (PMC6061858)",
    "6LAS": "Sun A, Gasser C, Li F, Chen H, Mair S, Krasheninina O, Micura R, Ren A (2019) Nat Commun 10:5728, doi:10.1038/s41467-019-13600-9 (PMC6914780)",
    "3NPQ": "Edwards AL, Reyes FE, Heroux A, Batey RT (2010) RNA 16:2144-2155, doi:10.1261/rna.2341610 (PMC2957054)",
    "2GDI": "Serganov A, Polonskaia A, Phan AT, Breaker RR, Patel DJ (2006) Nature 441:1167-1171, doi:10.1038/nature04740 (PMC4689313)",
    "1FJG": "Carter AP, Clemons WM Jr, Brodersen DE, Morgan-Warren RJ, Wimberly BT, Ramakrishnan V (2000) Nature 407:340-348, doi:10.1038/35030019",
}

COPY = {  # pdb -> (ligand_chain, ligand_resnum, rna_chain)
    "1Y26": ("X", 90, "X"), "4FE5": ("B", 101, "B"), "2GIS": ("A", 301, "A"), "2QWY": ("A", 100, "A"),
    "3E5C": ("A", 216, "A"), "6FZ0": ("A", 104, "A"), "6LAS": ("A", 101, "A"), "3NPQ": ("A", 55, "A"),
    "2GDI": ("X", 100, "X"),
}
CCD = {"1Y26": "ADE", "4FE5": "HPA", "2GIS": "SAM", "2QWY": "SAM", "3E5C": "SAM", "6FZ0": "SAM", "6LAS": "SAM", "3NPQ": "SAH", "2GDI": "TPP"}
FJG = {"PAR": 1545, "SRY": 1634, "SCM": 1633}

rows = []


def c(pdb, resnum, resname, ratom, latom, typ, ev, loc, *, ccd=None, mediator=None, dist=None, atom_basis="stated",
      type_basis="stated", paper_num=None, unit=None, qualifier=None, note=None, src=None):
    if pdb == "1FJG":
        lch, lnum, rch = "A", FJG[ccd], "A"
    else:
        lch, lnum, rch = COPY[pdb]
        ccd = ccd or CCD[pdb]
    rows.append(dict(
        pdb=pdb, ligand_chain=lch, ligand_resnum=lnum, ligand_ccd=ccd, rna_chain=rch,
        rna_resnum=resnum, rna_resname=resname, paper_numbering=paper_num,
        rna_atom_or_group=ratom, ligand_atom_or_group=latom, type=typ, mediator=mediator, distance_A=dist,
        evidence=ev, atom_basis=atom_basis, type_basis=type_basis, unit=unit, qualifier=qualifier,
        source=f"{SRC[src or pdb]}; {loc}", note=note))


# ---------------- 1Y26 add A-riboswitch / adenine ----------------
L = "Results, 'Adenine Recognition Specificity'; Fig. 7A"
c("1Y26", 74, "U", "N3", "N1", "hbond", "text", L, atom_basis="implied_by_pairing",
  note="Text: adenine N1 and N6H2 form a pair of H-bonds with the Watson-Crick edge of U74 (Watson-Crick A-U74 pair). U74 atom (N3) implied by the stated Watson-Crick pairing, not named.")
c("1Y26", 74, "U", "O4", "N6", "hbond", "text", L, atom_basis="implied_by_pairing",
  note="U74 atom (O4) implied by the stated Watson-Crick pairing, not named.")
c("1Y26", 51, "U", "Watson-Crick edge (atoms not named)", "N3", "hbond", "text", L, atom_basis="group_only",
  note="Text: N3 and N9H of adenine form a pair of H-bonds with the Watson-Crick edge of U51; the U51 atom for each H-bond is not named.")
c("1Y26", 51, "U", "Watson-Crick edge (atoms not named)", "N9 (N9-H)", "hbond", "text", L, atom_basis="group_only")
c("1Y26", 47, "U", "O2", "N9 (N9-H)", "hbond", "text", L)
c("1Y26", 22, "U", "O2'", "N7", "hbond", "text", L)
L2 = "Discussion, 'Purine Binding Pocket' (bound purine stacks over both flanking triples); Fig. 5C, 5D"
for rn, nm in ((73, "A"), (52, "A"), (22, "U")):
    c("1Y26", rn, nm, "base", "adenine ring", "stacking", "text", L2, atom_basis="group_only", unit="A73*(A52-U22) triple",
      note="Authors state stacking of the ligand over the whole base triple above it; which of the three bases overlaps the ligand is not stated. Score at triple (unit) level.")
for rn, nm in ((50, "C"), (75, "U"), (21, "A")):
    c("1Y26", rn, nm, "base", "adenine ring", "stacking", "text", L2, atom_basis="group_only", unit="C50*(U75-A21) triple",
      note="Authors state stacking of the ligand over the whole base triple below it; individual stacked base not stated. Score at triple (unit) level.")
c("1Y26", None, None, "not stated (RNA ligands of Mg5 not described)", "N6", "metal_mediated",
  "text", "Results, 'Divalent Cation Binding Sites'; Fig. 5C, 5D", mediator="hydrated Mg2+ (Mg5)", atom_basis="group_only",
  note="Hydrated Mg5 lies within hydrogen-bonding distance of adenine N6 (presumably via a coordinated water; not stated). No RNA partner named, so this record has no rna_resnum and is not residue-scorable.")

# ---------------- 4FE5 xpt guanine riboswitch / hypoxanthine ----------------
L = "main text (HX bound through H-bonds with U22, U47, U51 and C74, forming a base quadruple that stacks on P1); Fig. 2b"
c("4FE5", 22, "U", "residue (atom not stated in primary)", "hypoxanthine (atom not stated in primary)", "hbond", "text", L,
  atom_basis="group_only",
  note="Atom detail from secondary source (same lab, purine riboswitch class description): N7 of the purine ... 2'-hydroxyl of U22 (Batey 2012 review, PMC3760793, section 2.1, Fig. 2c). Serganov 2004 (1Y26/xpt-G, Chem Biol) states the same U22 2'-OH...N7 contact for adenine and guanine.")
c("4FE5", 47, "U", "residue (atom not stated in primary)", "hypoxanthine (atom not stated)", "hbond", "text", L, atom_basis="group_only")
c("4FE5", 51, "U", "residue (atom not stated in primary)", "hypoxanthine (atom not stated)", "hbond", "text", L, atom_basis="group_only")
c("4FE5", 74, "C", "Watson-Crick edge (atoms not named)", "Watson-Crick edge (N1/O6; not named)", "hbond", "text", L, atom_basis="implied_by_pairing",
  note="Primary text: specificity via Watson-Crick pairing preference at nucleotide 74 (C74 pairs with the purine); with hypoxanthine (no 2-NH2) this implies HX N1-H...C74 N3 and HX O6...C74 N4-H, but atoms are not named for HX.")
c("4FE5", 22, "U", "O2'", "N7", "hbond", "secondary", "section 2.1; Fig. 2c", src="4FE5_sec",
  note="Secondary record giving the atom pair for the U22 contact above (do not double count at residue level).")

# ---------------- 2GIS SAM-I / SAM ----------------
F4b = "Fig. 4b (distances in A as labelled); text: adenine ring of SAM is the central base of a triple between A45 and U57"
c("2GIS", 57, "U", "O4", "N6", "hbond", "figure", F4b, dist=2.9, atom_basis="figure")
c("2GIS", 57, "U", "N3", "N7", "hbond", "figure", F4b, dist=2.9, atom_basis="figure")
c("2GIS", 45, "A", "N1", "N6", "hbond", "figure", F4b, dist=2.8, atom_basis="figure")
c("2GIS", 45, "A", "O2'", "N1", "hbond", "figure", F4b, dist=3.0, atom_basis="figure")
c("2GIS", 47, "C", "base", "adenine ring", "stacking", "text", "main text (placement of adenine ring further stabilized by stacking with C47)", atom_basis="group_only")
c("2GIS", 47, "C", "O4' (printed as 'O4 of C47')", "O2'", "hbond", "text", "main text (ribose 2'-OH forms the only direct ribose H-bond, to 'O4 of C47')",
  atom_basis="group_only", note="Paper says 'O4 of C47'; cytidine has no O4, so the ribose O4' is the likely meaning. rna_atom left unparsed on purpose.")
F4c = "Fig. 4c (distances in A); text: carboxylate recognized by the Watson-Crick face of G11 (N1 and N2); main-chain atoms recognized by the G11*C44-G58 triple"
c("2GIS", 11, "G", "N1", "carboxylate O (O/OXT)", "hbond", "text", F4c, dist=3.1, atom_basis="stated",
  note="RNA atoms N1 and N2 stated in text; which carboxylate O pairs with which, and the distances, from Fig. 4c.")
c("2GIS", 11, "G", "N2", "carboxylate O (O/OXT)", "hbond", "text", F4c, dist=3.1)
c("2GIS", 58, "G", "N2", "carboxylate O (O/OXT)", "hbond", "figure", F4c, dist=2.8, atom_basis="figure",
  note="G58 named in text only as part of the G11*C44-G58 triple; atom from Fig. 4c.")
c("2GIS", 58, "G", "ribose hydroxyl (O2' in drawing)", "alpha-amino N", "hbond", "figure", F4c, dist=3.2, atom_basis="figure",
  note="Drawn to a G58 ribose oxygen; most likely O2' but not labelled, so left unparsed.")
F4d = "main text (sulphur ~4 A from the O2 carbonyls of U7 and U88; favourable electrostatic interaction); Fig. 4d"
c("2GIS", 7, "U", "O2", "SD (sulfonium S+)", "electrostatic", "text", F4d, dist=4.0, qualifier="approximate distance ('about 4 A')")
c("2GIS", 88, "U", "O2", "SD (sulfonium S+)", "electrostatic", "text", F4d, dist=4.0, qualifier="approximate distance ('about 4 A')",
  note="Sun et al. 2019 (SAM-VI paper, Suppl. Fig. 4a) describe SAM-I as O4 of U7/U88 at 4.3/4.0 A, apparently for another SAM-I entry; the primary paper says O2 (minor groove). Primary followed.")
for rn in (7, 88):
    c("2GIS", rn, "U", "nucleotide (minor-groove face)", "ribose/sulfonium face of SAM", "other", "text",
      "main text (interactions with helix P1 restricted to van der Waals packing with U7 and U88)", atom_basis="group_only",
      type_basis="stated", note="van der Waals packing; not scored (type other).")

# ---------------- 2QWY SAM-II / SAM ----------------
F3b = "Fig. 3b (left); text: adenine uses its Hoogsteen face to pair with U44 in a U10/U44 base triple"
c("2QWY", 44, "U", "N3", "N7", "hbond", "figure", F3b, atom_basis="figure")
c("2QWY", 44, "U", "O2", "N6", "hbond", "figure", F3b, atom_basis="figure")
c("2QWY", 10, "U", "O4", "C8", "other", "figure", F3b, atom_basis="figure", type_basis="curator_assigned",
  note="Dashed line drawn from U10 O4 to the adenine C8 position (C-H...O); text only says the adenine is in a triple with U10 and U44. Type not stated; not scored.")
F3c = "main text (sulfonium and activated methyl recognized by the carbonyl oxygens of U11 and U21); Fig. 3c (left)"
c("2QWY", 11, "U", "O4", "SD (sulfonium S+) and CE (S-methyl)", "electrostatic", "text", F3c, atom_basis="figure",
  note="Carbonyl identity (O4) from Fig. 3c; the text says 'carbonyl oxygen groups'. Legend: hydrogen-bonding and electrostatic interactions.")
c("2QWY", 21, "U", "O4", "SD (sulfonium S+) and CE (S-methyl)", "electrostatic", "text", F3c, atom_basis="figure")
F3d = "main text (A47 forms two H-bonds through N1 and N6 to the amino-acid groups, in protomers A and C); Fig. 3d (left)"
c("2QWY", 47, "A", "N1", "alpha-amino N", "hbond", "text", F3d, qualifier="present in protomers A and C, not B",
  note="Text says N1 and N6 bond 'to the carbonyl and amino groups ... respectively', which is chemically reversed (N1 acceptor, N6 donor) and contradicts Fig. 3d, which draws N1...amino and N6...carboxylate. Figure assignment used.")
c("2QWY", 47, "A", "N6", "carboxylate O (O/OXT)", "hbond", "text", F3d, qualifier="present in protomers A and C, not B")

# ---------------- 3E5C SAM-III (SMK box) / SAM ----------------
# PDB 3E5C chain A is numbered 1-53; paper uses E. faecalis metK transcript numbering. Mapping (sequence + helix consistency
# with Fig. 1b): PDB 1-12 = paper 20-31; 13-28 = engineered P3 (no paper numbers); 29-41 = paper 64-76; 42-45 = replaced
# linker; 46-53 = paper 88-95.
L3d = "main text (three minor-groove contacts N1-amino, amino-N3, amino-2'-OH with G26); Fig. 3d (distances in A, assignment of the two 2.9 A labels read from the figure)"
c("3E5C", 7, "G", "N2", "N1", "hbond", "text", L3d, dist=2.9, paper_num="G26",
  note="Text lists 'N1-amino' (SAM N1 ... G26 amino N2), 'amino-N3' (SAM N6 ... G26 N3), 'amino-2'-OH' (SAM N6 ... G26 2'-OH).")
c("3E5C", 7, "G", "N3", "N6", "hbond", "text", L3d, dist=2.9, paper_num="G26")
c("3E5C", 7, "G", "O2'", "N6", "hbond", "text", L3d, dist=2.8, paper_num="G26")
c("3E5C", 38, "A", "N6", "N7", "hbond", "text", "main text (N7 accepts an H-bond from the N6 amine of A73); Fig. 3d", dist=2.9, paper_num="A73")
L3e = "main text and Fig. 3e (distances in A)"
c("3E5C", 47, "G", "N7", "O2'", "hbond", "text", L3e, dist=3.1, paper_num="G89", qualifier="pi-hydrogen bond (stated)")
c("3E5C", 47, "G", "phosphoryl O (non-bridging; legend) / 2'-OH (main text)", "O3'", "hbond", "text", L3e, dist=2.5, paper_num="G89",
  atom_basis="group_only",
  note="Conflict in the paper: main text says SAM 3'-OH H-bonds the 2'-hydroxyl of G89; Fig. 3e legend and drawing show a phosphoryl oxygen of G89 (2.5 A). Residue-level contact is unambiguous; atom left unparsed.")
c("3E5C", 39, "A", "N6", "O2'", "water_mediated", "text", "main text (2'-hydroxyl of SAM donates an additional water-mediated H-bond to N6 of A74)",
  mediator="water", paper_num="A74")
L3f = "main text (sulfonium recognized by electrostatic interactions with O4 of U72 and 2'-hydroxyl of G71); Fig. 3f (distances in A)"
c("3E5C", 37, "U", "O4", "SD (sulfonium S+)", "electrostatic", "text", L3f, dist=4.1, paper_num="U72")
c("3E5C", 36, "G", "O2'", "SD (sulfonium S+)", "electrostatic", "text", L3f, dist=3.9, paper_num="G71")
c("3E5C", 48, "G", "base", "adenine ring", "stacking", "text", "Fig. 3a legend (adenosine stacks between U72 and G90); main text (G90 stacks directly underneath SAM)",
  atom_basis="group_only", paper_num="G90")
c("3E5C", 37, "U", "base", "adenine ring", "stacking", "text", "Fig. 3a legend (adenosine of SAM base-stacks between U72 and G90); main text: sheared U72*A64 pair sits immediately on top of the SAM intercalation plane",
  atom_basis="group_only", paper_num="U72")

# ---------------- 6FZ0 SAM-V / SAM ----------------
L6a = "Results, 'SAM binding by the riboswitch' (trans Hoogsteen pair with U47: N6 to U47 O2, U47 N3 to N7); Fig. 6A"
c("6FZ0", 47, "U", "O2", "N6", "hbond", "text", L6a)
c("6FZ0", 47, "U", "N3", "N7", "hbond", "text", L6a)
for rn, nm in ((21, "G"), (48, "A")):
    c("6FZ0", rn, nm, "base", "adenine ring", "stacking", "text", "Results (adenine intercalated into P2, stacked between G21 and A48)", atom_basis="group_only")
L6b = "Results (aminoacyl group H-bonded to A50: N6 donates to carboxyl O, N1 accepts from the amine); Fig. 6B"
c("6FZ0", 50, "A", "N6", "carboxylate O (O/OXT)", "hbond", "text", L6b)
c("6FZ0", 50, "A", "N1", "alpha-amino N", "hbond", "text", L6b)
c("6FZ0", 10, "U", "pro-R non-bridging phosphate O (OP1/OP2 not specified)", "carboxylate O (both O and OXT)", "metal_mediated", "text",
  "Results (metal-O 2.1 A, most likely Mg2+; metal bonded directly to the pro-R phosphate O of U10; both carboxylate O H-bonded to inner-sphere waters); Fig. 6C",
  mediator="Mg2+ (assigned from 2.1 A metal-O distance) with inner-sphere waters", atom_basis="group_only",
  note="RNA side: inner-sphere Mg-OP bond; ligand side: outer-sphere (via coordinated waters). 2.1 A is the metal-oxygen distance, not RNA-ligand.")
L6d = "Results (sulfur 3.2 and 3.3 A from the O4 carbonyls of U20 and U9; 'no direct contact' but electrostatic stabilization); Fig. 6D"
c("6FZ0", 20, "U", "O4", "SD (sulfonium S+)", "electrostatic", "text", L6d, dist=3.2)
c("6FZ0", 9, "U", "O4", "SD (sulfonium S+)", "electrostatic", "text", L6d, dist=3.3)

# ---------------- 6LAS SAM-VI / SAM ----------------
L7f = "Results, 'Nucleoside alignments in the binding pocket' (Hoogsteen edge of SAM-adenine pairs with the Watson-Crick edge of U8); Fig. 2f"
c("6LAS", 8, "U", "N3", "N7", "hbond", "figure", L7f, atom_basis="figure", note="Pairing edges stated in text; atom pairs from Fig. 2f drawing.")
c("6LAS", 8, "U", "O2", "N6", "hbond", "figure", L7f, atom_basis="figure",
  note="Drawn from SAM N6 to U8 O2 (the same O2 also accepts from G38 N2 per text). A '3.3' label lies between this line and the G38 line; not assigned.")
c("6LAS", 8, "U", "O4", "alpha-amino N", "hbond", "text", "Results (methionine tail folds back; H-bond between the alpha-amino group and O4 of U8); Fig. 2f")
c("6LAS", 36, "A", "N6", "N1", "hbond", "text", "Results (6-NH2 of A36 and N1 of A37 each H-bond the Watson-Crick edge of SAM adenine); Discussion (A36 interacts with N1 of SAM adenine); Fig. 2f")
c("6LAS", 37, "A", "N1", "N6", "hbond", "text", "Results (N1 of A37 H-bonds the Watson-Crick edge of SAM adenine); Fig. 2f",
  atom_basis="stated", note="Ligand atom N6 is implied (A37 N1 is an acceptor; the only donor on the SAM Watson-Crick edge is N6); drawing agrees.")
c("6LAS", 33, "G", "non-bridging phosphate O (OP1/OP2 not specified)", "O2'", "hbond", "text", "Results (non-bridging phosphate oxygen of G33 H-bonds the 2'-OH of SAM); Discussion; Fig. 2f",
  atom_basis="group_only")
for rn, nm, u in ((32, "C", None), (33, "G", None)):
    c("6LAS", rn, nm, "base", "adenine ring", "stacking", "text", "Results (C32 and G33 sandwich the SAM-adenine base and stack on it from both sides); Fig. 2c, 2e", atom_basis="group_only")
c("6LAS", 9, "G", "base", "adenine ring", "stacking", "text", "Results (terminal pair G9-C32 of P2 stacks on one side of the SAM adenine); Fig. 2e", atom_basis="group_only",
  unit="G9-C32 pair", note="Stated at base-pair level.")
c("6LAS", 7, "G", "base", "adenine ring", "stacking", "text", "Results (trans WC pair G7*G33 stacks on the other side); Discussion (G7 stacks with the SAM adenine); Fig. 2e", atom_basis="group_only",
  unit="G7*G33 pair")
c("6LAS", 6, "U", "O4", "SD (sulfonium S+)", "electrostatic", "text", "Results (sulfonium near O4 of U6, 4.8 A; O4 likely assists recognition of the positive charge); Fig. 2g",
  dist=4.8, qualifier="authors: 'likely assists'; long distance")

# ---------------- 3NPQ SAH riboswitch / SAH ----------------
L8 = "Results, SAH recognition paragraph; Fig. 3A, 3B"
c("3NPQ", 15, "G", "sugar edge (sheared G*A; atoms not named)", "adenine Hoogsteen edge (sheared G*A; atoms not named)", "hbond", "text", L8,
  atom_basis="group_only", note="Sheared G15*A(SAH) pair stated; Fig. 3B draws two H-bonds; atoms not labelled.")
for rn, nm in ((29, "A"), (16, "C")):
    c("3NPQ", rn, nm, "base", "adenine ring", "stacking", "text", "Results (adenine ring intercalates between the nucleobases of A29 and C16); Fig. 1C", atom_basis="group_only")
L8b = "Results (alpha-amine close to N3 and a non-bridging phosphate oxygen of G30: charge-charge interaction); Fig. 3A"
c("3NPQ", 30, "G", "N3", "alpha-amino N", "electrostatic", "text", L8b)
c("3NPQ", 30, "G", "non-bridging phosphate O (OP1/OP2 not specified)", "alpha-amino N", "electrostatic", "text", L8b, atom_basis="group_only")
L8c = "Results (carboxylate oxygens H-bond the 2'-hydroxyls of G47 and G31); Fig. 3A"
c("3NPQ", 47, "G", "O2'", "carboxylate O (O/OXT)", "hbond", "text", L8c, qualifier="main-chain contacts differ slightly among the 3 protomers")
c("3NPQ", 31, "G", "O2'", "carboxylate O (O/OXT)", "hbond", "text", L8c, qualifier="main-chain contacts differ slightly among the 3 protomers")
c("3NPQ", 29, "A", "O4'", "SD (thioether S)", "other", "text", "Results (sulfur in van der Waals contact with O4' of A29, 3.2 A); Fig. 3A", dist=3.2,
  note="van der Waals contact; not scored (type other).")

# ---------------- 2GDI TPP riboswitch / TPP ----------------
L9 = "main text (HMP ring intercalates between G42 and A43 and H-bonds the G40 base and the 2'-OH of G19); Fig. 2b, 2c"
c("2GDI", 40, "G", "base (atoms not named)", "HMP (4-amino-5-hydroxymethyl-2-methylpyrimidine) polar groups (atoms not named)", "hbond", "text", L9, atom_basis="group_only")
c("2GDI", 19, "G", "O2'", "HMP polar group (atom not named)", "hbond", "text", L9)
for rn, nm in ((42, "G"), (43, "A")):
    c("2GDI", rn, nm, "base", "HMP (pyrimidine) ring", "stacking", "text", L9 + "; continuous stacking of A41, G42, HMP, A43 and G21", atom_basis="group_only")
L9b = "main text (Mg1 directly coordinated to O6 of G60 and G78; terminal phosphate coordinated to Mg1 and Mg2); Fig. 2b, 2d"
c("2GDI", 60, "G", "O6", "terminal (beta) phosphate O (PB oxygens)", "metal_mediated", "text", L9b, mediator="Mg2+ (Mg1, hexa-coordinated)",
  note="Inner-sphere on both sides: TPP terminal phosphate and G60 O6 both coordinate Mg1.")
c("2GDI", 78, "G", "O6", "terminal (beta) phosphate O (PB oxygens)", "metal_mediated", "text", L9b, mediator="Mg2+ (Mg1, hexa-coordinated)")
L9c = "main text (conserved C77 and G78 form H-bonds to the oxygen atoms of the terminal phosphate)"
c("2GDI", 77, "C", "residue (atom not named)", "terminal (beta) phosphate O", "hbond", "text", L9c, atom_basis="group_only")
c("2GDI", 78, "G", "residue (atom not named)", "terminal (beta) phosphate O", "hbond", "text", L9c, atom_basis="group_only")

# ---------------- 1FJG 30S / paromomycin, streptomycin, spectinomycin (E. coli numbering = PDB author numbering) ----------------
P = "main text (paromomycin section); Fig. 6b chemical diagram"
c("1FJG", 1491, "G", "base", "ring I", "stacking", "text", P, ccd="PAR", atom_basis="group_only")
c("1FJG", 1408, "A", "N1", "ring I O6' (6'-OH)", "hbond", "text", P, ccd="PAR", atom_basis="figure",
  note="Text: ring I H-bonds with A1408; atoms from Fig. 6b.")
c("1FJG", 1493, "A", "phosphate O (OP1/OP2 not specified)", "ring I O4' (4'-OH, per drawing)", "hbond", "text", P, ccd="PAR", atom_basis="figure",
  note="Text: tight H-bond of ring I with the phosphate backbone of A1493; ring I hydroxyl identity read from Fig. 6b.")
c("1FJG", 1491, "G", "N7", "ring III O5'' (5''-OH)", "hbond", "figure", P, ccd="PAR", atom_basis="figure",
  note="Text: ring III makes only weak contacts with the RNA.")
c("1FJG", 1494, "G", "N7", "ring II N3 (amine)", "hbond", "figure", P, ccd="PAR", atom_basis="figure",
  note="Text: ring II forms tight interactions with bases and backbone; atoms from Fig. 6b (ring II amine numbering inferred from the drawn ring substitution).")
c("1FJG", 1495, "U", "base carbonyl O (O4 or O2; not resolvable)", "ring II N1 (amine)", "hbond", "figure", P, ccd="PAR", atom_basis="figure")
c("1FJG", 1490, "C", "phosphate O (OP1/OP2 not specified)", "ring IV N6''' (6'''-NH2)", "hbond", "figure", P, ccd="PAR", atom_basis="figure",
  type_basis="curator_assigned", note="Text: ring IV contacts the backbone of both sides of H44; contact type not stated per contact (dashed line).")
c("1FJG", 1405, "G", "phosphate O (OP1/OP2 not specified)", "ring IV O3''' (3'''-OH)", "hbond", "figure", P, ccd="PAR", atom_basis="figure",
  type_basis="curator_assigned")
S = "main text (streptomycin bound to the phosphate backbone from four parts of 16S RNA through salt bridges and hydrogen bonds); Fig. 5b chemical diagram"
NS = "Per-contact type not stated (text: 'salt bridges and hydrogen bonds'); curator assigned electrostatic for cationic N...phosphate and hbond otherwise."
c("1FJG", 1490, "C", "phosphate O (OP1/OP2 not specified)", "streptidine guanidino NH2 (C1 guanidine, per drawing)", "electrostatic", "figure", S, ccd="SRY",
  atom_basis="figure", type_basis="curator_assigned", note=NS)
c("1FJG", 1491, "G", "phosphate O (OP1/OP2 not specified)", "streptidine C6-OH (per drawing)", "hbond", "figure", S, ccd="SRY",
  atom_basis="figure", type_basis="curator_assigned", note=NS)
c("1FJG", 914, "A", "ribose O (two lines; O2'/O3'/O4' not resolvable)", "streptidine guanidino NH2 (C3 guanidine, per drawing)", "hbond", "figure", S, ccd="SRY",
  atom_basis="figure", type_basis="curator_assigned", note=NS)
c("1FJG", 914, "A", "phosphate O (OP1/OP2 not specified)", "N-methyl-L-glucosamine methylamino N (C2'')", "electrostatic", "figure", S, ccd="SRY",
  atom_basis="figure", type_basis="curator_assigned", note=NS)
c("1FJG", 14, "U", "phosphate O (OP1/OP2 not specified)", "N-methyl-L-glucosamine hydroxyl (O4'' per drawing)", "hbond", "figure", S, ccd="SRY",
  atom_basis="figure", type_basis="curator_assigned", note=NS + " Text cites nucleotide 13 from biochemical data; the structure figure labels U14.")
c("1FJG", 527, "G", "phosphate O (OP1/OP2 not specified)", "streptose C3'-OH", "hbond", "figure", S, ccd="SRY",
  atom_basis="figure", type_basis="curator_assigned", note=NS)
c("1FJG", 526, "C", "phosphate O (OP1/OP2 not specified)", "streptose formyl (aldehyde) O", "hbond", "figure", S, ccd="SRY",
  atom_basis="figure", type_basis="curator_assigned", note=NS)
c("1FJG", 913, "A", "ribose O (O2' or O3'; not resolvable)", "streptose formyl (aldehyde) O", "hbond", "figure", S, ccd="SRY",
  atom_basis="figure", type_basis="curator_assigned", note=NS + " Text cites 915 from biochemical data; the structure figure labels A913 and A914 (and says A915 makes no contacts with other bases).")
SP = "main text (spectinomycin: minor groove at the end of H34, single contact with a 2'-OH, H-bonds to several bases, most with G1064 and C1192); Fig. 4b chemical diagram"
c("1FJG", 1064, "G", "O6", "actinamine hydroxyl flanked by an N-methylamino group and the dioxane O (drawn at top of ring A)", "hbond", "figure", SP, ccd="SCM", atom_basis="figure")
c("1FJG", 1064, "G", "N7", "same actinamine hydroxyl", "hbond", "figure", SP, ccd="SCM", atom_basis="figure")
c("1FJG", 1066, "C", "O2", "actinamine N-methylamino N (the one adjacent to the G1064-bonded hydroxyl)", "hbond", "figure", SP, ccd="SCM", atom_basis="figure")
c("1FJG", 1068, "G", "O2'", "actinamine N-methylamino N (the other one)", "hbond", "text", SP, ccd="SCM", atom_basis="figure",
  note="Text: single contact with a 2'-OH; Fig. 4b identifies G1068 and the methylamino group.")
c("1FJG", 1192, "C", "N3", "hemiketal/hydrated-ketone hydroxyl on the ring B/C bridgehead carbon (adjacent to the ketone)", "hbond", "figure", SP, ccd="SCM", atom_basis="figure",
  note="Text names C1192 among the most-contacting bases; atom pair from Fig. 4b.")

json.dump(rows, open(OUT, "w"), indent=1, ensure_ascii=False)
print(len(rows), "records")
import collections
print(collections.Counter((r["pdb"], r["ligand_ccd"]) for r in rows))
print(collections.Counter(r["type"] for r in rows), collections.Counter(r["evidence"] for r in rows))
