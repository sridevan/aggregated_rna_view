"""Read-only mmCIF inventory of polymer component identities for every RF00167 RNA instance (Part B investigation; not a production pipeline).

Inputs (never modified):
  * data/families/RF00167/instances.json and mappings/<seq_name>.json   (validated Rfam mapping; consensus positions come only from here)
  * cache/www_ebi_ac_uk/pdbe_entry-files_download_<pdb>_updated.cif-*.cif  (PDBe updated mmCIF already used by the pipeline)
  * raw/pdbe/batch_all_entries.json                                        (PDBe modified_AA_or_NA, fetched by fetch_pdbe.py)
Outputs (audit/):
  * polymer_positions.csv   one row per (instance, label_seq_id): deposited sequence, coordinates, mod_residue, mapping, PDBe record
  * nonstandard_observations.csv / .json   rows whose deposited component is not A/C/G/U
  * entry_components.json   every nonstandard component per entry and where it occurs (other polymers, non-polymer ligands)
  * checks.json             consistency checks and counts
Identity is read at label_seq_id from _pdbx_poly_seq_scheme (mon_id) and checked against _entity_poly_seq, _atom_site (model 1) and the mapping.
"""
import csv, glob, json, collections, pathlib, gemmi

HERE = pathlib.Path(__file__).resolve().parents[1]; ROOT = HERE.parents[1]
FAM = ROOT / "data/families/RF00167"; AUD = HERE / "audit"; AUD.mkdir(exist_ok=True)
STD = {"A", "C", "G", "U"}

def cif_path(pdb):
    hits = glob.glob(str(ROOT / f"cache/www_ebi_ac_uk/pdbe_entry-files_download_{pdb}_updated.cif-*.cif"))
    assert len(hits) == 1, (pdb, hits); return pathlib.Path(hits[0])

def table(block, cat, cols):
    t = block.find(cat, cols) if block.find_mmcif_category(cat) else None
    return [[gemmi.cif.as_string(r[i]) for i in range(len(cols))] for r in t] if t is not None else []

def has_cat(block, cat): return bool(block.find_mmcif_category(cat))

pdbe = json.load(open(HERE / "raw/pdbe/batch_all_entries.json"))["body"]
pdbe_idx = {}   # (pdb, struct_asym_id, residue_number) -> list of records (duplicates kept to detect them)
for e, recs in pdbe.items():
    for r in recs: pdbe_idx.setdefault((e, r["struct_asym_id"], int(r["residue_number"])), []).append(r)

instances = json.load(open(FAM / "instances.json"))
rows, entry_comp, checks = [], {}, collections.defaultdict(list)
blocks = {}
for inst in instances:
    pdb, asym, ent = inst["pdb_id"], inst["label_asym_id"], str(inst["entity_id"])
    if pdb not in blocks: blocks[pdb] = gemmi.cif.read(str(cif_path(pdb))).sole_block()
    b = blocks[pdb]
    mp = json.load(open(FAM / inst["mapping_file"]))
    assert mp["label_asym_id"] == asym and str(mp["entity_id"]) == ent
    mres = {r["polymer_seq_position"]: r for r in mp["residues"]}
    ep = {r[0]: r for r in table(b, "_entity_poly.", ["entity_id", "type", "pdbx_seq_one_letter_code", "pdbx_seq_one_letter_code_can"])}
    eps = collections.defaultdict(list)
    for e_id, num, mon, het in table(b, "_entity_poly_seq.", ["entity_id", "num", "mon_id", "hetero"]):
        if e_id == ent: eps[int(num)].append((mon, het))
    pss = {}
    for a, e_id, sid, mon, pnum, ins, pmon, amon, het in table(b, "_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "seq_id", "mon_id", "pdb_seq_num", "pdb_ins_code", "pdb_mon_id", "auth_mon_id", "hetero"]):
        if a == asym: assert e_id == ent; pss.setdefault(int(sid), []).append(dict(mon=mon, pdb_seq_num=pnum, ins="" if ins in ".?" else ins, pdb_mon=pmon, auth_mon=amon, het=het))
    modres = {}
    for a, sid, comp, parent, det in table(b, "_pdbx_struct_mod_residue.", ["label_asym_id", "label_seq_id", "label_comp_id", "parent_comp_id", "details"]):
        if a == asym: modres[int(sid)] = (comp, parent, det)
    atoms = collections.defaultdict(set); grp = collections.defaultdict(set)
    first_model = None
    for a, sid, comp, g, model in table(b, "_atom_site.", ["label_asym_id", "label_seq_id", "label_comp_id", "group_PDB", "pdbx_PDB_model_num"]):
        first_model = first_model or model
        if a == asym and model == first_model and sid not in ".?": atoms[int(sid)].add(comp); grp[int(sid)].add(g)
    checks["entity_poly_type"].append((inst["seq_name"], ep[ent][1]))
    assert set(eps) == set(pss) == set(mres), (inst["seq_name"], len(eps), len(pss), len(mres))
    for sid in sorted(pss):
        s = pss[sid]; m = mres[sid]
        mons_scheme = sorted({x["mon"] for x in s}); mons_eps = sorted({x[0] for x in eps[sid]})
        hetero = any(x["het"] in ("y", "Y") for x in s) or len(mons_eps) > 1
        coord_comp = sorted(atoms.get(sid, set()))
        comp = mons_scheme[0] if len(mons_scheme) == 1 else "|".join(mons_scheme)
        pr = pdbe_idx.get((pdb, asym, sid), [])
        row = dict(seq_name=inst["seq_name"], pdb_id=pdb, entity_id=ent, label_asym_id=asym, auth_asym_id=inst["auth_asym_id"], label_seq_id=sid,
                   auth_seq_id=s[0]["pdb_seq_num"], ins_code=s[0]["ins"], comp_poly_seq_scheme=comp, comp_entity_poly_seq="|".join(mons_eps), hetero=hetero,
                   comp_atom_site="|".join(coord_comp), atom_site_group="|".join(sorted(grp.get(sid, set()))), has_coordinates_mmcif=bool(coord_comp),
                   mod_residue_comp=modres.get(sid, ("", "", ""))[0], mod_residue_parent=modres.get(sid, ("", "", ""))[1],
                   mapping_residue_name=m["residue_name"], mapping_parent=m["nucleotide_parent"], mapping_is_modified=m["is_modified"], mapping_has_coordinates=m["has_coordinates"],
                   alignment_state=m["alignment_state"], rfam_consensus_pos=m["rfam_consensus_pos"], insertion_anchor=m["insertion_anchor"],
                   pdbe_comp="|".join(r["chem_comp_id"] for r in pr), pdbe_n_records=len(pr), pdbe_author_residue_number="|".join(str(r["author_residue_number"]) for r in pr),
                   pdbe_chain_id="|".join(r["chain_id"] for r in pr), pdbe_entity_id="|".join(str(r["entity_id"]) for r in pr))
        row["nonstandard"] = any(c not in STD for c in mons_scheme)
        row["is_terminal"] = "5prime" if sid == min(pss) else ("3prime" if sid == max(pss) else "")
        rows.append(row)
        # consistency checks
        if mons_scheme != mons_eps: checks["scheme_vs_entity_poly_seq"].append((inst["seq_name"], sid, mons_scheme, mons_eps))
        if coord_comp and coord_comp != mons_scheme: checks["scheme_vs_atom_site"].append((inst["seq_name"], sid, mons_scheme, coord_comp))
        if bool(coord_comp) != bool(m["has_coordinates"]): checks["coords_mmcif_vs_mapping"].append((inst["seq_name"], sid, bool(coord_comp), m["has_coordinates"]))
        if comp != m["residue_name"]: checks["scheme_vs_mapping_residue_name"].append((inst["seq_name"], sid, comp, m["residue_name"]))
        if row["nonstandard"] != bool(m["is_modified"]): checks["nonstandard_vs_mapping_is_modified"].append((inst["seq_name"], sid, comp, m["is_modified"]))
        if s[0]["pdb_seq_num"] != str(m["author_residue_number"]) or s[0]["ins"] != (m["author_insertion_code"] or ""): checks["author_number_vs_mapping"].append((inst["seq_name"], sid, s[0]["pdb_seq_num"], m["author_residue_number"]))
        if len(pr) > 1: checks["pdbe_duplicate_records"].append((inst["seq_name"], sid, len(pr)))
        if pr and pr[0]["chem_comp_id"] != comp: checks["pdbe_vs_mmcif_comp"].append((inst["seq_name"], sid, pr[0]["chem_comp_id"], comp))

# every nonstandard component per entry, wherever it occurs (polymers of any type, non-polymer ligands), to exclude unrelated occurrences
for pdb, b in blocks.items():
    etype = {r[0]: r[1] for r in table(b, "_entity_poly.", ["entity_id", "type"])}
    ent_desc = {r[0]: (r[1], r[2]) for r in table(b, "_entity.", ["id", "type", "pdbx_description"])}
    cc = {r[0]: dict(type=r[1], name=r[2]) for r in table(b, "_chem_comp.", ["id", "type", "name"])}
    comp_where = collections.defaultdict(lambda: collections.Counter())
    for a, e_id, sid, mon in table(b, "_pdbx_poly_seq_scheme.", ["asym_id", "entity_id", "seq_id", "mon_id"]):
        if mon not in STD: comp_where[mon][f"polymer entity {e_id} ({etype.get(e_id, '?')}) asym {a}"] += 1
    for a, e_id, mon in table(b, "_pdbx_nonpoly_scheme.", ["asym_id", "entity_id", "mon_id"]):
        comp_where[mon][f"non-polymer entity {e_id} asym {a}"] += 1
    entry_comp[pdb] = {"entities": {k: {"type": v[0], "description": v[1], "polymer_type": etype.get(k)} for k, v in ent_desc.items()},
                       "chem_comp": cc, "nonstandard_or_nonpolymer_components": {k: dict(v) for k, v in comp_where.items()},
                       "has_pdbx_struct_mod_residue": has_cat(b, "_pdbx_struct_mod_residue.")}

cols = list(rows[0])
with open(AUD / "polymer_positions.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
ns = [r for r in rows if r["nonstandard"]]
with open(AUD / "nonstandard_observations.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(ns)
json.dump(ns, open(AUD / "nonstandard_observations.json", "w"), indent=1)
json.dump(entry_comp, open(AUD / "entry_components.json", "w"), indent=1)
# PDBe records that did not match any mapped-instance position
matched = {(r["pdb_id"], r["label_asym_id"], r["label_seq_id"]) for r in rows if r["pdbe_n_records"]}
checks["pdbe_records_not_on_mapped_instance"] = [(e, k[1], k[2], [r["chem_comp_id"] for r in v]) for (e, *_), v in pdbe_idx.items() for k in [(e, *_)] if (e, k[1], k[2]) not in matched]
summary = {k: (len(v), v[:20]) for k, v in checks.items() if k != "entity_poly_type"}
summary["entity_poly_types"] = dict(collections.Counter(t for _, t in checks["entity_poly_type"]))
summary["n_instances"] = len(instances); summary["n_entries"] = len(blocks); summary["n_polymer_positions"] = len(rows)
summary["mmcif_source"] = "PDBe updated mmCIF (https://www.ebi.ac.uk/pdbe/entry-files/download/<pdb>_updated.cif), pipeline cache"
json.dump(summary, open(AUD / "checks.json", "w"), indent=1, default=str)
for k, v in summary.items(): print(k, v if not isinstance(v, tuple) else (v[0], v[1][:6]))
