"""Independent reconstruction of the Part B counts with different code paths and inputs:
  * residue identity from gemmi.read_structure (entity full_sequence + subchain residues), not from the raw CIF tables
  * mapping from data/families/RF00167/residue_mapping.csv and each mapping JSON's consensus_columns, not the JSON residues list
  * PDBe from the per-entry responses (raw/pdbe/single_entries.json), not the batch
Prints the reconstructed numbers and asserts they equal audit/inventory.json."""
import csv, json, glob, collections, pathlib, gemmi
HERE = pathlib.Path(__file__).resolve().parents[1]; ROOT = HERE.parents[1]; FAM = ROOT / "data/families/RF00167"
STD = {"A", "C", "G", "U"}
inst = json.load(open(FAM / "instances.json"))
csvrows = collections.defaultdict(dict)
for r in csv.DictReader(open(FAM / "residue_mapping.csv")): csvrows[r["seq_name"]][int(r["polymer_seq_position"])] = r
obs = []
for i in inst:
    st = gemmi.read_structure(glob.glob(str(ROOT / f"cache/www_ebi_ac_uk/pdbe_entry-files_download_{i['pdb_id']}_updated.cif-*.cif"))[0])
    st.setup_entities()
    ent = next(e for e in st.entities if i["label_asym_id"] in e.subchains)
    seq = [gemmi.Entity.first_mon(m) for m in ent.full_sequence]
    assert len(seq) == len(csvrows[i["seq_name"]]), i["seq_name"]
    resolved = {r.label_seq: r.name for r in st[0].get_subchain(i["label_asym_id"])}
    cols = {c["polymer_seq_position"]: c["rfam_consensus_pos"] for c in json.load(open(FAM / i["mapping_file"]))["consensus_columns"] if c["state"] == "match"}
    for k, mon in enumerate(seq, start=1):
        r = csvrows[i["seq_name"]][k]
        assert r["residue_name"] == mon and (k in resolved) == (r["has_coordinates"] == "1") and (k not in resolved or resolved[k] == mon), (i["seq_name"], k)
        if r["alignment_state"] == "match": assert cols.get(k) == int(r["rfam_consensus_pos"]), (i["seq_name"], k)
        if mon not in STD: obs.append((i["pdb_id"], i["seq_name"], k, mon, r["alignment_state"], r["rfam_consensus_pos"], k in resolved))
single = json.load(open(HERE / "raw/pdbe/single_entries.json"))
pdbe = {(e, x["struct_asym_id"], x["residue_number"], x["chem_comp_id"]) for e, v in single.items() if v["http_status"] == 200 for x in v["body"][e]}
lab = {(i["seq_name"]): i["label_asym_id"] for i in inst}
mine = {(o[0], lab[o[1]], o[2], o[3]) for o in obs}
rec = dict(n_nonstandard_residue_observations=len(obs), n_instances_with_nonstandard=len({o[1] for o in obs}), n_entries_with_nonstandard=len({o[0] for o in obs}),
           n_distinct_ccd=len({o[3] for o in obs}), n_on_match=sum(o[4] == "match" for o in obs), n_on_insert=sum(o[4] == "insert" for o in obs),
           n_outside_mapped_region=sum(o[4] == "not_submitted_flank" for o in obs), distinct_positions=sorted({int(o[5]) for o in obs if o[4] == "match"}),
           n_pdbe_records=len(pdbe), pdbe_minus_mmcif=sorted(pdbe - mine), mmcif_minus_pdbe=sorted(mine - pdbe))
print(json.dumps(rec, indent=1))
ref = json.load(open(HERE / "audit/inventory.json"))["summary"]
for k in ("n_nonstandard_residue_observations", "n_instances_with_nonstandard", "n_entries_with_nonstandard", "n_distinct_ccd", "n_on_match", "n_on_insert", "n_outside_mapped_region", "distinct_positions", "n_pdbe_records"):
    assert rec[k] == ref[k], (k, rec[k], ref[k])
json.dump(rec, open(HERE / "audit/independent_check.json", "w"), indent=1)
print("independent reconstruction matches audit/inventory.json")
