"""Read-only fetch of PDBe modified_AA_or_NA for every RF00167 entry: one batched POST plus one POST per entry (to check batch completeness).
Raw responses (body + HTTP status + UTC time) are saved under raw/pdbe/. Usage: python3 scripts/fetch_pdbe.py"""
import json, subprocess, datetime, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[3]; OUT = pathlib.Path(__file__).resolve().parents[1] / "raw" / "pdbe"
URL = "https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/modified_AA_or_NA"
inst = json.load(open(ROOT / "data/families/RF00167/instances.json"))
inst = inst["instances"] if isinstance(inst, dict) else inst
entries = sorted({i["pdb_id"].lower() for i in inst})
def post(ids):
    r = subprocess.run(["curl", "-sS", "-X", "POST", URL, "-H", "accept: */*", "-H", "Content-Type: application/json", "-d", json.dumps(",".join(ids)), "-w", "\n%{http_code}"], capture_output=True, text=True, timeout=120)
    body, code = r.stdout.rsplit("\n", 1)
    return {"request_ids": ids, "http_status": int(code), "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), "body": json.loads(body) if body.strip() else None}
batch = post(entries); json.dump(batch, open(OUT / "batch_all_entries.json", "w"), indent=1)
single = {e: post([e]) for e in entries}; json.dump(single, open(OUT / "single_entries.json", "w"), indent=1)
print(len(entries), "entries; batch HTTP", batch["http_status"], "entries in batch body:", sorted(batch["body"]) if isinstance(batch["body"], dict) else batch["body"])
print("single-entry status counts:", {c: sum(1 for v in single.values() if v["http_status"] == c) for c in {v["http_status"] for v in single.values()}})
print("single 200:", sorted(e for e, v in single.items() if v["http_status"] == 200))
