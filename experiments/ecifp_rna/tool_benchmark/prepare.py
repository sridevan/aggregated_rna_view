"""Pilot benchmark of RNA-ligand interaction tools: common inputs per ligand copy.

For every pilot ligand copy:
1. the PDBe updated mmCIF (scratch only, never stored in the repository), model 1, first alternate conformer only;
2. a crop: the ligand plus every nucleic-acid residue, water and single-atom ion with any atom within 12 A of the ligand (whole
   residues). Proteins are left out for all tools alike (fingeRNAt accepts nucleic-acid receptors only), which matters for 1FJG
   streptomycin (S12 contact) and is noted in the report;
3. written heavy-atom only (crop.cif, crop.pdb) and with a common hydrogenation (hydride 1.2.3 on CCD bond orders and the deposited
   formal charges from the CCD; hydrogens on phosphate OP1/OP2/OP3 removed and OP1 given a -1 charge, so the backbone phosphate is anionic):
   crop_H.cif, crop_H.pdb, rna_H.pdb (receptor: nucleic acid, water, ions) and ligand_H.sdf (for fingeRNAt).

usage: prepare.py <work_dir>
"""
import sys, json
from pathlib import Path
import numpy as np, requests, gemmi
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as bpdb
import biotite.structure.io.mol as bmol
import hydride

WORK = Path(sys.argv[1])
CROP = 12.0
PILOT = [  # pdb, ligand auth chain, ligand auth seq id, CCD, family
    ("1y26", "X", 90, "ADE", "RF00167"), ("4fe5", "B", 101, "HPA", "RF00167"),
    ("2gis", "A", 301, "SAM", "RF00162"), ("2qwy", "A", 100, "SAM", "RF00521"), ("3e5c", "A", 216, "SAM", "RF01767"),
    ("6fz0", "A", 104, "SAM", "RF01826"), ("6las", "A", 101, "SAM", "RF02885"), ("3npq", "A", 55, "SAH", "RF01057"),
    ("2gdi", "X", 100, "TPP", "RF00059"),
    ("1fjg", "A", 1545, "PAR", "RF00177"), ("1fjg", "A", 1634, "SRY", "RF00177"), ("1fjg", "A", 1633, "SCM", "RF00177"),
]


def entry(pdb):
    p = WORK / "entries" / f"{pdb}_updated.cif"
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}_updated.cif", timeout=300); r.raise_for_status()
        p.write_bytes(r.content)
    return p


def find_copy(st, ch, num, ccd):
    for c in st[0]:
        if c.name != ch:
            continue
        for r in c:
            if r.seqid.num == num and r.name == ccd:
                return r
    return None


def crop(pdb, ch, num, ccd, out):
    st = gemmi.read_structure(str(entry(pdb)))
    st.setup_entities()
    for m in list(st)[1:]:
        del st[1]
    st.remove_alternative_conformations(); st.remove_hydrogens()
    lig = find_copy(st, ch, num, ccd)
    if lig is None:
        raise ValueError(f"{pdb} {ch} {num} {ccd} not found")
    lc = np.array([a.pos.tolist() for a in lig])
    keep = {}
    for c in st[0]:
        for i, r in enumerate(c):
            info = gemmi.find_tabulated_residue(r.name)
            na = info is not None and info.is_nucleic_acid()
            ion = len(r) == 1 and r[0].element.is_metal
            wat = r.is_water()
            if not (r is lig or (r.name == ccd and r.seqid.num == num and c.name == ch) or na or ion or wat or
                    (r.het_flag == "A" and info is not None and not info.is_amino_acid())):
                continue
            xyz = np.array([a.pos.tolist() for a in r])
            if (c.name == ch and r.seqid.num == num and r.name == ccd) or np.min(np.linalg.norm(xyz[:, None] - lc[None], axis=2)) <= CROP:
                keep.setdefault(c.name, set()).add(i)
    for c in st[0]:
        for i in reversed(range(len(c))):
            if i not in keep.get(c.name, set()):
                del c[i]
    st.remove_empty_chains()
    out.mkdir(parents=True, exist_ok=True)
    doc = st.make_mmcif_document(); doc.write_file(str(out / "crop.cif"))
    st.shorten_chain_names() if any(len(c.name) > 1 for c in st[0]) else None
    st.write_pdb(str(out / "crop.pdb"), gemmi.PdbWriteOptions(minimal=True))
    return dict(residues=sum(len(c) for c in st[0]), atoms=sum(len(r) for c in st[0] for r in c))


def hydrogenate(out, ch, num, ccd):
    f = pdbx.CIFFile.read(str(out / "crop.cif"))
    a = pdbx.get_structure(f, model=1, include_bonds=True, extra_fields=["charge"], use_author_fields=True)
    a = a[a.element != "H"]
    # formal charges from the CCD (deposited files often leave pdbx_formal_charge empty, e.g. the SAM sulfonium in 2GIS)
    import biotite.structure.info as binfo
    tmpl = {}
    for rn in set(a.res_name):
        try:
            t = binfo.residue(rn); tmpl[rn] = dict(zip(t.atom_name, t.charge))
        except Exception:
            tmpl[rn] = {}
    a.charge = np.array([tmpl.get(rn, {}).get(an, c) for rn, an, c in zip(a.res_name, a.atom_name, a.charge)], dtype=int)
    # backbone phosphate anionic: OP1 -1 (hydride then adds no H to it)
    a.charge[(a.atom_name == "OP1") & np.isin(a.res_name, ["A", "C", "G", "U"])] = -1
    h, mask = hydride.add_hydrogen(a)
    h.coord = hydride.relax_hydrogen(h)
    # drop any H on OP1/OP2/OP3 of nucleotides
    bonds = h.bonds.as_array()
    po = np.isin(h.atom_name, ["OP1", "OP2", "OP3"]) & np.isin(h.res_name, ["A", "C", "G", "U"])
    drop = set()
    for i, j, _ in bonds:
        if h.element[i] == "H" and po[j]: drop.add(i)
        if h.element[j] == "H" and po[i]: drop.add(j)
    h = h[np.array([k not in drop for k in range(len(h))])]
    cf = pdbx.CIFFile(); pdbx.set_structure(cf, h, include_bonds=True); cf.write(str(out / "crop_H.cif"))
    pf = bpdb.PDBFile(); pf.set_structure(h); pf.write(str(out / "crop_H.pdb"))
    is_lig = (h.chain_id == ch) & (h.res_id == num) & (h.res_name == ccd)
    rec = h[~is_lig & (struc.filter_nucleotides(h) | (h.res_name == "HOH") | struc.filter_monoatomic_ions(h))]
    pf = bpdb.PDBFile(); pf.set_structure(rec); pf.write(str(out / "rna_H.pdb"))
    lig = h[is_lig].copy()
    lig.bonds.remove_aromaticity()          # V2000 has no aromatic bond order outside query files: write the Kekule form
    sd = bmol.SDFile(); rec_ = bmol.SDRecord(); rec_.set_structure(lig); sd[f"{ccd}_{ch}_{num}"] = rec_; sd.write(str(out / "ligand_H.sdf"))
    return dict(atoms_H=len(h), ligand_atoms_H=int(is_lig.sum()), removed_phosphate_H=len(drop))


if __name__ == "__main__":
    summary = {}
    for pdb, ch, num, ccd, fam in PILOT:
        out = WORK / "inputs" / f"{pdb}_{ch}_{num}_{ccd}"
        try:
            s = crop(pdb, ch, num, ccd, out); s.update(hydrogenate(out, ch, num, ccd)); s["family"] = fam
        except Exception as e:
            s = dict(error=f"{type(e).__name__}: {e}")
        summary[out.name] = s
        print(out.name, s)
    json.dump(summary, open(WORK / "inputs" / "summary.json", "w"), indent=1)
