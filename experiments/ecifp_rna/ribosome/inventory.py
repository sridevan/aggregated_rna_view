"""Ribosome ligand inventory: how many usable rRNA-drug complexes exist per drug class, rRNA, species group, method and resolution?

Sources (PDBe, 2026-10-23):
- PDBe search (Solr) `rfam_accession:<ACC>` for RF00177 (bacterial SSU, incl. mitochondrial 12S), RF02541 (bacterial LSU, incl.
  mitochondrial 16S), RF01960 (eukaryotic SSU) and RF02543 (eukaryotic LSU): entries, rRNA chains, organism, method, resolution,
  component list. (The PDBe API nucleic_mappings endpoint did not answer within 5 min for these families.)
- PDBe compound summary for names and types; drug classes assigned by hand below (CCD -> class), everything else left as 'other'.
- PDBe API bound_molecules and bound_ligand_interactions per ligand copy of a classed drug: contacts present (or 404), contacted
  chains, and whether any contacted chain is an rRNA chain of the four families.

usage: inventory.py <cache_dir> <out_dir>
"""
import json, sys, os, time, collections, urllib.request, urllib.parse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

CACHE, OUT = Path(sys.argv[1]), Path(sys.argv[2])
CACHE.mkdir(parents=True, exist_ok=True); OUT.mkdir(parents=True, exist_ok=True)
FAMS = {"RF00177": "SSU", "RF02541": "LSU", "RF01960": "SSU", "RF02543": "LSU"}
KINGDOM = {"RF00177": "bacterial", "RF02541": "bacterial", "RF01960": "eukaryotic", "RF02543": "eukaryotic"}

CLASSES = {
    "aminoglycoside (2-DOS)": "PAR NMY GET SIS 84G LLL KAN TOY 84D RPO ON0 M5Z LUJ B6M HKO HJR HJO AKN 3TS EDS AB9 AM2",
    "streptomycin": "SRY", "hygromycin B": "HYG HY0", "spectinomycin": "SCM", "kasugamycin": "KSG", "negamycin": "NEG",
    "tetracycline": "TAC T1C DXT YQM V7A MIY U3B P8F 80P D2C 8UZ",
    "macrolide/ketolide": "ERY TEL ZIT DI0 CTY ROX TAO TYK JOS LMA MT9 MIV QU2 QTZ 773 ERN EM1 T17 4M2 1F2 1F3 1F4 RU6 UHP UH6 UI6 A1AEZ A1AE0 A1AE1 UIU AQJ",
    "lincosamide": "CLY 3QB 6IF WC9 A1A1F A1A1J A1A1K A1A1L A1I09",
    "phenicol": "CLM EZP EZM EZG A1AIX YXM",
    "oxazolidinone": "ZLD RD8 U7V U7Y ZC0 G6V G6M 917 A1BWO JJH UI9 UH0 UIF",
    "pleuromutilin": "MUL G34 62B G19 G80 A1CQ5 A1CQ7 A1CQ8 A1CQ9 A1CRD A1CRE",
    "streptogramin": "VIR VIF H8T DOL O7S O7V O7Y O8D O8J O8P O8S O8V UCX UEC UI0 UE6 A1C7K",
    "orthosomycin": "6O1 EVN 6NO 6UQ",
    "PTC/tunnel other (bacterial)": "SPS HGR BLS UAM M2D 3V6 LC2 3H3 MQH A1IC4 6EK KKL",
    "SSU other (bacterial)": "PCY EDE T8B FSD OIY OI9",
    "eukaryotic inhibitors": "ANM 3HE 3K5 YAT 34G K16 MQ6 HMT 3KD 3KF 3K8 3J2 ZBA 3J6 VDU 7AL YMZ WDP A1L32 SJH SJE 3L2 ZVM",
}
CCD2CLASS = {c: k for k, v in CLASSES.items() for c in v.split()}


def get(url, name):
    p = CACHE / name
    if p.exists():
        return json.loads(p.read_text())
    for k in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                d = json.loads(r.read())
            break
        except urllib.error.HTTPError as e:
            d = None if e.code == 404 else "retry"
            if d is None:
                break
        except Exception:
            d = "retry"
        time.sleep(2 + 3 * k)
    if d == "retry":
        return "error"
    p.write_text(json.dumps(d))
    return d


def search(acc):
    q = urllib.parse.urlencode({"q": f"rfam_accession:{acc}", "rows": 20000, "wt": "json",
                                "fl": "pdb_id,chain_id,organism_scientific_name,experimental_method,resolution,release_date,compound_id,molecule_name"})
    return get("https://www.ebi.ac.uk/pdbe/search/pdb/select?" + q, f"search_{acc}.json")["response"]["docs"]


def main():
    entries = {}
    for acc in FAMS:
        for x in search(acc):
            e = entries.setdefault(x["pdb_id"], dict(rrna={}, method=(x.get("experimental_method") or [None])[0], resolution=x.get("resolution"),
                                                     release=x.get("release_date", "")[:10], comps=set()))
            orgs = x.get("organism_scientific_name") or []
            for ch in x.get("chain_id") or []:
                e["rrna"][ch] = dict(acc=acc, organism=orgs[0] if orgs else None, molecule=(x.get("molecule_name") or [None])[0])
            e["comps"] |= set(x.get("compound_id") or [])
    drug_entries = {p: e for p, e in entries.items() if e["comps"] & set(CCD2CLASS)}
    print(f"rRNA entries {len(entries)}; with a classed drug {len(drug_entries)}")

    def per_entry(pdb):
        bm = get(f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/bound_molecules/{pdb}", f"bm_{pdb}.json")
        rows = []
        if not isinstance(bm, dict):
            return [dict(pdb=pdb, status="bound_molecules_" + str(bm))]
        for b in bm.get(pdb, []):
            for l in b["composition"]["ligands"]:
                ccd = l["chem_comp_id"]
                if ccd not in CCD2CLASS:
                    continue
                it = get(f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/bound_ligand_interactions/{pdb}/{l['chain_id']}/{l['author_residue_number']}",
                         f"bli_{pdb}_{l['chain_id']}_{l['author_residue_number']}.json")
                row = dict(pdb=pdb, ccd=ccd, cls=CCD2CLASS[ccd], chain=l["chain_id"], num=l["author_residue_number"])
                if not isinstance(it, dict) or pdb not in it:
                    row["status"] = "no_contacts" if it is None else "error"
                else:
                    ints = it[pdb][0]["interactions"]
                    good = [x for x in ints if set(x["interaction_details"]) - {"proximal", "vdw", "vdw_clash", "clash"}]
                    chains = collections.Counter(x["end"]["chain_id"] for x in good)
                    rr = {c: n for c, n in chains.items() if c in drug_entries[pdb]["rrna"]}
                    row.update(status="ok" if rr else "no_rrna_contacts", rrna_chains=rr, other_chains={c: n for c, n in chains.items() if c not in rr},
                               rrna_residues=sorted({(x["end"]["chain_id"], x["end"]["author_residue_number"]) for x in good if x["end"]["chain_id"] in rr}))
                    if rr:
                        top = max(rr, key=rr.get)
                        info = drug_entries[pdb]["rrna"][top]
                        row.update(rrna_acc=info["acc"], organism=info["organism"])
                rows.append(row)
        return rows

    rows = []
    with ThreadPoolExecutor(6) as ex:
        for r in ex.map(per_entry, sorted(drug_entries)):
            rows += r
    for r in rows:
        e = drug_entries[r["pdb"]]
        r.update(method=e["method"], resolution=e["resolution"], release=e["release"])
    json.dump(dict(n_entries=len(entries), n_drug_entries=len(drug_entries), classes=CLASSES, rows=rows,
                   entries={p: dict(e, comps=sorted(e["comps"])) for p, e in drug_entries.items()}),
              open(OUT / "inventory.json", "w"), indent=1)
    print(collections.Counter(r.get("status") for r in rows))


main()
