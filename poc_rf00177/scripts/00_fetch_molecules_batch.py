"""Fetch PDBe molecules + residue-listing-free entity metadata for all RF00177 entries via POST batch; fall back to GET."""
import json, subprocess, sys, time
from pathlib import Path
DATA = Path(sys.argv[1])
ids = sorted(json.load(open(DATA/"api/rf00177_nucleic_mappings.json"))["RF00177"]["PDB"].keys())
out = {}
def post(chunk):
    r = subprocess.run(["curl","-sS","-X","POST","https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/","-d",",".join(chunk)],capture_output=True,text=True)
    try: return json.loads(r.stdout)
    except Exception: return {}
CH=40
for i in range(0,len(ids),CH):
    chunk=ids[i:i+CH]; got=post(chunk)
    missing=[p for p in chunk if p not in got]
    for p in missing:
        r=subprocess.run(["curl","-sS",f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{p}"],capture_output=True,text=True)
        try: got.update(json.loads(r.stdout))
        except Exception: print("FAILED",p,file=sys.stderr)
    out.update(got); print(f"{i+len(chunk)}/{len(ids)} fetched ({len(missing)} via GET fallback)",flush=True)
json.dump(out,open(DATA/"api/rf00177_molecules_all.json","w"))
print("entries with molecules:",len(out))
