"""_struct_ref_seq_dif rows of the family RNA entities in the cached PDBe updated mmCIF (RF00167, RF00162, RF00059), mapped through the
validated family mappings; compared with the PDBe mutated_AA_or_NA response for the same entries."""
import glob, json, gemmi
from collections import Counter
from pathlib import Path
D = Path(__file__).parent
COLS = ["align_id", "pdbx_pdb_id_code", "mon_id", "pdbx_pdb_strand_id", "seq_num", "pdbx_pdb_ins_code", "pdbx_seq_db_name", "pdbx_seq_db_accession_code", "db_mon_id", "pdbx_seq_db_seq_num", "details", "pdbx_auth_seq_num", "pdbx_ordinal"]
rows = []
for acc in ["RF00167", "RF00162", "RF00059"]:
    inst = json.load(open(f"data/families/{acc}/instances.json"))
    by_pdb = {}
    for i in inst: by_pdb.setdefault(i["pdb_id"].lower(), []).append(i)
    pdbe = json.load(open(D / "raw" / f"{acc}.json"))["body"]
    for pdb, insts in sorted(by_pdb.items()):
        b = gemmi.cif.read(glob.glob(f"cache/www_ebi_ac_uk/pdbe_entry-files_download_{pdb}_updated.cif-*.cif")[0]).sole_block()
        if not b.find_mmcif_category("_struct_ref_seq_dif."): continue
        t = b.find("_struct_ref_seq_dif.", [c if c in ("align_id", "mon_id", "details") else "?" + c for c in COLS])
        # align_id -> entity via struct_ref_seq.ref_id -> struct_ref.entity_id
        rs = {r[0]: r[1] for r in b.find("_struct_ref_seq.", ["align_id", "ref_id"])}
        re_ = {r[0]: (r[1], r[2], r[3]) for r in b.find("_struct_ref.", ["id", "entity_id", "db_name", "db_code"])}
        # auth strand -> label asym of the same entity
        asym = {}
        for r in b.find("_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "pdb_strand_id"]): asym[(r[1], r[2])] = r[0]
        fam_entities = {str(i["entity_id"]) for i in insts}
        for r in t:
            v = {c: (gemmi.cif.as_string(r[k]) if r.has(k) else None) for k, c in enumerate(COLS)}
            ent, dbn, dbc = re_.get(rs.get(v["align_id"]), (None, None, None))
            if ent not in fam_entities: continue
            lab = asym.get((ent, v["pdbx_pdb_strand_id"]))
            i = next((x for x in insts if x["label_asym_id"] == lab), None)
            row = {"acc": acc, "pdb_id": pdb, "entity_id": ent, "ref_db": dbn, "ref_code": dbc, "auth_asym_id": v["pdbx_pdb_strand_id"], "label_asym_id": lab, "seq_num": v["seq_num"],
                   "auth_seq_num": v["pdbx_auth_seq_num"], "ins_code": v["pdbx_pdb_ins_code"], "mon_id": v["mon_id"], "db_mon_id": v["db_mon_id"], "details": v["details"],
                   "other_chains_of_entity": [x["label_asym_id"] for x in insts if str(x["entity_id"]) == ent and x["label_asym_id"] != lab],
                   "in_pdbe_api": any(x["struct_asym_id"] == lab and str(x["residue_number"]) == str(v["seq_num"]) for x in pdbe.get(pdb, []))}
            if i is None: row["mapping_status"] = "chain_not_resolved"
            elif v["seq_num"] in (None, "", "?", "."): row["mapping_status"] = "no_seq_num (residue absent from the deposited sequence)"
            else:
                m = {x["polymer_seq_position"]: x for x in json.load(open(f"data/families/{acc}/" + i["mapping_file"]))["residues"]}.get(int(v["seq_num"]))
                if m is None: row["mapping_status"] = "residue_not_in_mapping"
                else:
                    row.update({"seq_name": i["seq_name"], "map_residue": m["residue_name"], "map_author": f'{m["author_residue_number"]}{m["author_insertion_code"] or ""}',
                                "alignment_state": m["alignment_state"], "rfam_consensus_pos": m["rfam_consensus_pos"], "has_coordinates": m["has_coordinates"],
                                "identity_agrees": m["residue_name"] == v["mon_id"], "author_agrees": str(m["author_residue_number"]) == str(v["pdbx_auth_seq_num"]) and (m["author_insertion_code"] or "") == (v["pdbx_pdb_ins_code"] or "")})
                    row["mapping_status"] = "mapped" if m["alignment_state"] == "match" else f"no_consensus_position ({m['alignment_state']})"
            rows.append(row)
json.dump(rows, open(D / "mmcif_seqdif_family_rna.json", "w"), indent=1)
print(len(rows), "family-RNA seq_dif rows in", len({(r['acc'], r['pdb_id']) for r in rows}), "entries;", Counter(r["acc"] for r in rows))
print("details:", Counter(r["details"] for r in rows)); print("ref db:", Counter(r["ref_db"] for r in rows))
print("mapping:", Counter(r["mapping_status"] for r in rows)); print("identity agrees:", Counter(r.get("identity_agrees") for r in rows), "author agrees:", Counter(r.get("author_agrees") for r in rows))
print("from(db_mon)->to(mon):", Counter((r["db_mon_id"], r["mon_id"]) for r in rows).most_common(20))
print("in PDBe API:", Counter(r["in_pdbe_api"] for r in rows))
print("entries with >1 RNA chain of the entity:", sorted({r["pdb_id"] for r in rows if r["other_chains_of_entity"]}))
