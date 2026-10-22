"""Classify every PDBe mutated_AA_or_NA record: family RNA chain (mapped through the validated family mapping) or other chain."""
import json
from collections import Counter, defaultdict
from pathlib import Path
D = Path(__file__).parent
FAMS = ["RF00167", "RF00162", "RF00059", "RF00005", "RF00011", "RF00177"]
NA = {"A", "C", "G", "U", "DA", "DC", "DG", "DT", "N"}
rows = []
for acc in FAMS:
    raw = json.load(open(D / "raw" / f"{acc}.json"))["body"]
    inst = {(i["pdb_id"].lower(), i["label_asym_id"]): i for i in json.load(open(f"data/families/{acc}/instances.json"))}
    maps = {}
    for e, recs in raw.items():
        for x in recs:
            k = (e.lower(), x["struct_asym_id"]); i = inst.get(k)
            row = {"acc": acc, "pdb_id": e.lower(), "entity_id": x["entity_id"], "auth_asym_id": x["chain_id"], "label_asym_id": x["struct_asym_id"],
                   "author_residue_number": x["author_residue_number"], "author_insertion_code": x["author_insertion_code"] or "", "residue_number": x["residue_number"],
                   "chem_comp_id": x["chem_comp_id"], "from": x["mutation_details"]["from"], "to": x["mutation_details"]["to"], "type": x["mutation_details"]["type"]}
            if not i:
                row["chain_class"] = "family_entity_other_chain" if any(j["pdb_id"].lower() == e.lower() and int(j["entity_id"]) == x["entity_id"] for j in inst.values()) else "not_family_chain"
                row["looks_nucleic"] = x["chem_comp_id"] in NA
            else:
                if i["seq_name"] not in maps: maps[i["seq_name"]] = {r["polymer_seq_position"]: r for r in json.load(open(f"data/families/{acc}/" + i["mapping_file"]))["residues"]}
                r = maps[i["seq_name"]].get(int(x["residue_number"]))
                row.update({"chain_class": "family_rna_chain", "seq_name": i["seq_name"]})
                if r is None: row.update({"mapping_status": "residue_not_in_mapping"})
                else:
                    row.update({"map_residue_name": r["residue_name"], "map_author_number": r["author_residue_number"], "map_ins_code": r["author_insertion_code"] or "",
                                "alignment_state": r["alignment_state"], "rfam_consensus_pos": r["rfam_consensus_pos"], "has_coordinates": r["has_coordinates"],
                                "author_id_agrees": str(r["author_residue_number"]) == str(x["author_residue_number"]) and (r["author_insertion_code"] or "") == (x["author_insertion_code"] or ""),
                                "identity_agrees": r["residue_name"] == x["chem_comp_id"]})
                    row["mapping_status"] = "mapped" if r["alignment_state"] == "match" else f"no_consensus_position:{r['alignment_state']}"
            rows.append(row)
json.dump(rows, open(D / "records_classified.json", "w"), indent=1)
fam = [r for r in rows if r["chain_class"] == "family_rna_chain"]
print("records:", len(rows), Counter(r["chain_class"] for r in rows))
print("non-family records that look nucleic:", [(r["acc"], r["pdb_id"], r["label_asym_id"], r["chem_comp_id"], r["type"]) for r in rows if r["chain_class"] != "family_rna_chain" and r.get("looks_nucleic")][:20])
print("types (all):", Counter(r["type"] for r in rows)); print("types (family RNA):", Counter((r["acc"], r["type"]) for r in fam))
print("family RNA mapping:", Counter(r["mapping_status"] for r in fam), "author ids agree:", Counter(r.get("author_id_agrees") for r in fam), "identity agrees:", Counter(r.get("identity_agrees") for r in fam))
print("from/to values (family RNA):", Counter((r["from"], r["to"]) for r in fam).most_common(30))
