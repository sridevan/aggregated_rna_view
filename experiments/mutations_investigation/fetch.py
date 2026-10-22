"""Fetch PDBe mutated_AA_or_NA for every PDB entry of the six families (batches of 40), keep the raw responses."""
import json, datetime, requests
from pathlib import Path
OUT = Path(__file__).parent / "raw"; URL = "https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/mutated_AA_or_NA"
FAMS = ["RF00167", "RF00162", "RF00059", "RF00005", "RF00011", "RF00177"]
for acc in FAMS:
    ids = sorted({i["pdb_id"].lower() for i in json.load(open(f"data/families/{acc}/instances.json"))})
    body = {}
    for k in range(0, len(ids), 40):
        r = requests.post(URL, data=json.dumps(",".join(ids[k:k + 40])), headers={"Content-Type": "application/json"}, timeout=120)
        if r.status_code == 200: body.update(r.json())
        elif r.status_code != 404: raise SystemExit(f"{acc} HTTP {r.status_code}")
    json.dump({"url": URL, "method": "POST", "request_ids": ids, "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), "body": body}, open(OUT / f"{acc}.json", "w"), indent=1)
    print(acc, len(ids), "entries;", len(body), "with records;", sum(len(v) for v in body.values()), "records")
