"""Pilot benchmark: run the interaction tools on the prepared inputs and normalise every contact to one schema.

Arms (per ligand copy, inputs from prepare.py):
- arpeggio_api    PDBe API bound_ligand_interactions (production PDBe-Arpeggio on the full entry: the contacts this project uses)
- arpeggio_local  pdbe-arpeggio 1.4.4 on crop.cif (OpenBabel hydrogens). NB: with the pip OpenBabel build no ring is perceived as
                  aromatic, so this arm has no ring (stacking / cation-pi) contacts at all; kept to document that
- plip            PLIP 3.0.1 --dnareceptor on crop.pdb (PLIP's own OpenBabel polar hydrogens)
- plip_H          PLIP --dnareceptor --nohydro on crop_H.pdb (common hydride hydrogens)
- fingernat_H     fingeRNAt (cc62eb8) FULL -h2o -detail on rna_H.pdb + ligand_H.sdf (it requires hydrogens; common hydride ones)

Normalised types: hbond, weak_hbond, ionic, stacking, cation_pi, anion_pi, other_pi, hydrophobic, halogen, metal (ligand-metal
coordination), metal_mediated (ligand-metal-RNA), water_mediated. Every record keeps the tool's own type label.

usage: run_tools.py <work_dir> <fingernat_dir> <out_json>
"""
import json, sys, subprocess, collections, shutil, os
import xml.etree.ElementTree as ET
from pathlib import Path
import requests, gemmi

WORK, FR, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
BIN = Path(sys.executable).parent
API = "https://www.ebi.ac.uk/pdbe/api/v2/pdb/bound_ligand_interactions/{pdb}/{ch}/{num}"
ARP = {"hbond": "hbond", "polar": "hbond", "weak_hbond": "weak_hbond", "weak_polar": "weak_hbond", "ionic": "ionic",
       "hydrophobic": "hydrophobic", "halogen": "halogen", "xbond": "halogen", "metal_complex": "metal",
       "aromatic": "stacking", "FF": "stacking", "OF": "stacking", "EE": "stacking", "FT": "stacking", "OT": "stacking",
       "ET": "stacking", "FE": "stacking", "OE": "stacking", "EF": "stacking", "TT": "stacking", "TF": "stacking", "TE": "stacking",
       "TO": "stacking", "EO": "stacking", "FO": "stacking", "OO": "stacking",
       "CATIONPI": "cation_pi", "CARBONPI": "other_pi", "DONORPI": "other_pi", "HALOGENPI": "other_pi", "METSULPHURPI": "other_pi",
       "carbonyl": "other", "AMIDEAMIDE": "other", "AMIDERING": "other"}
SKIP = {"proximal", "vdw", "vdw_clash", "clash", "covalent"}
NUC = {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}


def rec(tool, copy, ch, num, name, ratoms, latoms, typ, label, mediator=None):
    return dict(tool=tool, copy=copy, rna=f"{ch}:{num}", rna_resname=name, rna_atoms=sorted(set(ratoms)), lig_atoms=sorted(set(latoms)),
                type=typ, label=label, mediator=mediator)


def rna_residues(inp):
    st = gemmi.read_structure(str(inp / "crop.cif"))
    out = set()
    for c in st[0]:
        for r in c:
            info = gemmi.find_tabulated_residue(r.name)
            if info is not None and info.is_nucleic_acid():
                out.add(f"{c.name}:{r.seqid.num}")
    return out


def arpeggio_records(js, tool, copy, ch, num, ccd, rna):
    out = []
    for x in js:
        if x.get("interacting_entities") != "INTER":
            continue
        b, e = x["bgn"], x["end"]
        if not (b["auth_asym_id"] == ch and b["auth_seq_id"] == num):
            b, e = e, b
        k = f"{e['auth_asym_id']}:{e['auth_seq_id']}"
        if k not in rna:
            continue
        for t in x["contact"]:
            if t in SKIP:
                continue
            out.append(rec(tool, copy, e["auth_asym_id"], e["auth_seq_id"], e["label_comp_id"], [e["auth_atom_id"]] if "auth_atom_id" in e else e.get("auth_atom_ids", []),
                           [b["auth_atom_id"]] if "auth_atom_id" in b else b.get("auth_atom_ids", []), ARP.get(t, "other"), t))
    return out


def run_arpeggio_api(pdb, ch, num, copy, rna, cache):
    p = cache / f"bli_{pdb}_{ch}_{num}.json"
    if not p.exists():
        r = requests.get(API.format(pdb=pdb, ch=ch, num=num), timeout=120)
        p.write_text(json.dumps(r.json() if r.status_code == 200 else None))
    d = json.loads(p.read_text())
    if not d:
        return None
    out = []
    for x in d[pdb][0]["interactions"]:
        e = x["end"]
        k = f"{e['chain_id']}:{e['author_residue_number']}"
        if k not in rna:
            continue
        for t in x["interaction_details"]:
            if t in SKIP:
                continue
            out.append(rec("arpeggio_api", copy, e["chain_id"], e["author_residue_number"], e["chem_comp_id"], e["atom_names"], x["ligand_atoms"], ARP.get(t, "other"), t))
    return out


def run_arpeggio_local(inp, ch, num, copy, rna, od):
    od.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(BIN / "pdbe-arpeggio"), str(inp / "crop.cif"), "-s", f"/{ch}/{num}/", "-o", str(od)], capture_output=True, check=True)
    js = json.load(open(od / "crop.json"))
    for x in js:   # ring contacts carry a list of atoms
        for side in ("bgn", "end"):
            if "auth_atom_id" not in x[side] and "auth_atom_ids" not in x[side]:
                x[side]["auth_atom_ids"] = []
    return arpeggio_records(js, "arpeggio_local", copy, ch, num, None, rna)


def pdb_index(p):
    idx = {}
    for line in open(p):
        if line.startswith(("ATOM", "HETATM")):
            idx[int(line[6:11])] = (line[21], int(line[22:26]), line[17:20].strip(), line[12:16].strip())
    return idx


def run_plip(inp, ch, num, ccd, copy, rna, od, hydro):
    od.mkdir(parents=True, exist_ok=True)
    f = inp / ("crop_H.pdb" if hydro else "crop.pdb")
    cmd = [str(BIN / "plip"), "-f", str(f), "--dnareceptor", "-x", "-o", str(od)] + (["--nohydro"] if hydro else [])
    subprocess.run(cmd, capture_output=True, check=True)
    xml = next(od.glob("*_report.xml"))
    prot = next(od.glob("*_protonated.pdb"), f)
    idx = pdb_index(prot)
    tool = "plip_H" if hydro else "plip"
    out = []
    root = ET.parse(xml).getroot()
    for bs in root.iter("bindingsite"):
        ident = bs.find("identifiers")
        if not (ident.find("hetid").text == ccd and ident.find("chain").text == ch and int(ident.find("position").text) == num):
            continue
        I = bs.find("interactions")
        def res(x):
            return x.find("reschain").text, int(x.find("resnr").text), x.find("restype").text
        def names(ids):
            return [idx[i][3] for i in ids if i in idx]
        def lig_or_rec(ids):
            r, l = [], []
            for i in ids:
                if i in idx:
                    (l if (idx[i][0] == ch and idx[i][1] == num) else r).append(idx[i][3])
            return r, l
        for x in I.find("hydrogen_bonds"):
            c, n, t = res(x); r, l = lig_or_rec([int(x.find("donoridx").text), int(x.find("acceptoridx").text)])
            out.append(rec(tool, copy, c, n, t, r, l, "hbond", "hydrogen_bond"))
        for x in I.find("salt_bridges"):
            c, n, t = res(x)
            r, l = lig_or_rec([int(i.text) for i in x.find("prot_idx_list")] + [int(i.text) for i in x.find("lig_idx_list")])
            out.append(rec(tool, copy, c, n, t, r, l, "ionic", "salt_bridge"))
        for x in I.find("pi_stacks"):
            c, n, t = res(x)
            r, l = lig_or_rec([int(i.text) for i in x.find("prot_idx_list")] + [int(i.text) for i in x.find("lig_idx_list")])
            out.append(rec(tool, copy, c, n, t, r, l, "stacking", f"pi_stack_{x.find('type').text}"))
        for x in I.find("pi_cation_interactions"):
            c, n, t = res(x)
            r, l = lig_or_rec([int(i.text) for i in x.find("prot_idx_list")] + [int(i.text) for i in x.find("lig_idx_list")])
            out.append(rec(tool, copy, c, n, t, r, l, "cation_pi", "pi_cation"))
        for x in I.find("hydrophobic_interactions"):
            c, n, t = res(x); r, l = lig_or_rec([int(x.find("ligcarbonidx").text), int(x.find("protcarbonidx").text)])
            out.append(rec(tool, copy, c, n, t, r, l, "hydrophobic", "hydrophobic"))
        for x in I.find("halogen_bonds"):
            c, n, t = res(x); r, l = lig_or_rec([int(x.find("don_idx").text), int(x.find("acc_idx").text)])
            out.append(rec(tool, copy, c, n, t, r, l, "halogen", "halogen_bond"))
        for x in I.find("water_bridges"):
            c, n, t = res(x); r, l = lig_or_rec([int(x.find("donor_idx").text), int(x.find("acceptor_idx").text)])
            w = idx.get(int(x.find("water_idx").text))
            out.append(rec(tool, copy, c, n, t, r, l, "water_mediated", "water_bridge", mediator=f"HOH{w[1]}" if w else "HOH"))
        for x in I.find("metal_complexes"):
            c, n, t = res(x)
            m = idx.get(int(x.find("metal_idx").text)); tg = idx.get(int(x.find("target_idx").text))
            if tg and not (tg[0] == ch and tg[1] == num):
                out.append(rec(tool, copy, c, n, t, [tg[3]], [], "metal", "metal_complex_receptor", mediator=m[2] if m else None))
            elif tg:
                out.append(rec(tool, copy, ch, num, ccd, [], [tg[3]], "metal_ligand", "metal_complex_ligand", mediator=f"{m[2]}{m[1]}" if m else None))
    # metal-mediated = a metal coordinated by the ligand (metal_ligand) and by an RNA atom (receptor metal complex of the same metal)
    lig_metals = {r["mediator"] for r in out if r["type"] == "metal_ligand"}
    for bs in root.iter("bindingsite"):
        ident = bs.find("identifiers")
        tag = f"{ident.find('hetid').text}{ident.find('position').text}"
        if tag in lig_metals:
            for x in bs.find("interactions").find("metal_complexes"):
                tg = idx.get(int(x.find("target_idx").text))
                if tg and f"{tg[0]}:{tg[1]}" in rna:
                    out.append(rec(tool, copy, tg[0], tg[1], tg[2], [tg[3]], [], "metal_mediated", "metal_complex_bridge", mediator=tag))
    return out


FRTYPE = {"HB": "hbond", "HAL": "halogen", "CA": "ionic", "Pi_Cation": "cation_pi", "Pi_Anion": "anion_pi", "Anion_Pi": "anion_pi",
          "Pi_Stacking": "stacking", "Sandwich": "stacking", "Displaced": "stacking", "T-shaped": "stacking", "Lipophilic": "hydrophobic",
          "Water-mediated": "water_mediated", "Mg_mediated": "metal_mediated", "K_mediated": "metal_mediated", "Na_mediated": "metal_mediated",
          "Other_mediated": "metal_mediated", "Ion-mediated": "metal_mediated"}


def run_fingernat(inp, ch, num, ccd, copy, rna, od):
    if od.exists():
        shutil.rmtree(od)
    od.mkdir(parents=True)
    p = subprocess.run([sys.executable, "code/fingeRNAt.py", "-r", str(inp / "rna_H.pdb"), "-l", str(inp / "ligand_H.sdf"), "-detail", "-h2o",
                        "-o", str(od)], cwd=FR, capture_output=True, text=True)
    det = next(od.glob("DETAIL_*.tsv"), None)
    if det is None:
        raise RuntimeError(p.stdout[-800:] + p.stderr[-800:])
    import biotite.structure.io.pdbx as pdbx
    h = pdbx.get_structure(pdbx.CIFFile.read(str(inp / "crop_H.cif")), model=1, use_author_fields=True)
    lig_names = list(h.atom_name[(h.chain_id == ch) & (h.res_id == num) & (h.res_name == ccd)])
    out, seen = [], set()
    rows = [l.rstrip("\n").split("\t") for l in open(det)]
    head = rows[0]
    for r in rows[1:]:
        d = dict(zip(head, r))
        key = tuple(d.get(k, "").strip() for k in head[4:])
        if key in seen:
            continue
        seen.add(key)
        if f'{d["Receptor_Chain"].strip()}:{int(float(d["Receptor_Number"]))}' not in rna:
            continue                                          # first row of a mediated pair is ligand-ion/water; keep the RNA side
        it = d["Interaction"]
        typ = FRTYPE.get(it, "metal_mediated" if "mediated" in it.lower() and "water" not in it.lower() else "other")
        la = d.get("Ligand_Atom", "")
        lname = [lig_names[int(a) - 1] for a in la.replace("[", "").replace("]", "").replace(",", " ").split() if a.isdigit() and int(a) - 1 < len(lig_names)]
        ra = d.get("Receptor_Atom", "").strip()
        out.append(rec("fingernat_H", copy, d["Receptor_Chain"].strip(), int(float(d["Receptor_Number"])), d["Receptor_Residue_Name"].strip(),
                       [a.strip() for a in ra.replace("[", "").replace("]", "").replace("'", "").split(",") if a.strip()], lname, typ, it,
                       mediator=(f'{d["Ligand_name"]}:{d["Ligand_Atom"]}' if typ in ("water_mediated", "metal_mediated") else None)))
    return out, head


if __name__ == "__main__":
    summary = json.load(open(WORK / "inputs" / "summary.json"))
    cache = WORK / "api"; cache.mkdir(exist_ok=True)
    allrec, status = [], {}
    for copy, s in summary.items():
        if "error" in s:
            continue
        pdb, ch, num, ccd = copy.split("_"); num = int(num)
        inp = WORK / "inputs" / copy
        rna = rna_residues(inp)
        st = {}
        for name, fn in (("arpeggio_api", lambda: run_arpeggio_api(pdb, ch, num, copy, rna, cache)),
                         ("arpeggio_local", lambda: run_arpeggio_local(inp, ch, num, copy, rna, WORK / "runs" / copy / "arpeggio")),
                         ("plip", lambda: run_plip(inp, ch, num, ccd, copy, rna, WORK / "runs" / copy / "plip", False)),
                         ("plip_H", lambda: run_plip(inp, ch, num, ccd, copy, rna, WORK / "runs" / copy / "plip_H", True)),
                         ("fingernat_H", lambda: run_fingernat(inp, ch, num, ccd, copy, rna, WORK / "runs" / copy / "fingernat")[0])):
            try:
                r = fn()
                if r is None:
                    st[name] = "no_data"; continue
                allrec += r
                st[name] = dict(collections.Counter(x["type"] for x in r))
            except Exception as e:
                st[name] = f"error: {type(e).__name__}: {str(e)[:300]}"
        status[copy] = st
        print(copy); [print(f"   {k}: {v}") for k, v in st.items()]
    json.dump(dict(status=status, records=allrec), open(OUT, "w"), indent=1)
