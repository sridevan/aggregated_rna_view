"""Step 1: build rf00177_instances.csv and summary from PDBe nucleic_mappings + molecules."""
import json, csv, sys, collections, statistics, hashlib
from pathlib import Path
DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
maps = json.load(open(DATA/"api/rf00177_nucleic_mappings.json"))["RF00177"]["PDB"]
mols = json.load(open(DATA/"api/rf00177_molecules_all.json"))
solr = {(d["pdb_id"], d["entity_id"]): d for d in json.load(open(DATA/"api/rf00177_solr.json"))["response"]["docs"]}
rows = []
for pdb, segs in maps.items():
    ents = {e["entity_id"]: e for e in mols.get(pdb, [])}
    for s in segs:
        e = ents.get(s["entity_id"], {})
        src = (e.get("source") or [{}])[0]
        sd = solr.get((pdb, s["entity_id"]), {})
        org = (sd.get("organism_scientific_name") or [src.get("organism_scientific_name")])[0]
        tax = (sd.get("tax_id") or [src.get("tax_id")])[0]
        seq = e.get("sequence", "")
        L = s["end"]["residue_number"] - s["start"]["residue_number"] + 1
        rows.append({"pdb_id": pdb, "entity_id": s["entity_id"], "auth_asym_id": s["chain_id"], "label_asym_id": s["struct_asym_id"],
                     "map_start": s["start"]["residue_number"], "map_end": s["end"]["residue_number"],
                     "map_auth_start": s["start"]["author_residue_number"], "map_auth_end": s["end"]["author_residue_number"],
                     "map_auth_start_ins": s["start"]["author_insertion_code"] or "", "map_auth_end_ins": s["end"]["author_insertion_code"] or "",
                     "mapped_length": L, "entity_length": e.get("length"), "molecule_name": "; ".join(e.get("molecule_name") or []),
                     "organism": org, "tax_id": tax, "organism_source": "pdbe_search" if sd.get("organism_scientific_name") else ("pdbe_molecules" if src.get("organism_scientific_name") else ""), "exp_method": ";".join(sd.get("experimental_method", [])), "resolution": sd.get("resolution"), "release_year": sd.get("release_year"), "title": sd.get("title", ""), "expression_host": src.get("expression_host_scientific_name"),
                     "n_chains_of_entity": len(e.get("in_chains") or []), "entity_chains": ",".join(e.get("in_chains") or []),
                     "sequence_md5": hashlib.md5(seq.encode()).hexdigest() if seq else "", "n_modified_positions": len(e.get("pdb_sequence_indices_with_multiple_residues") or {}),
                     "mapped_sequence": seq[s["start"]["residue_number"]-1:s["end"]["residue_number"]] if seq else ""})
rows.sort(key=lambda r: (r["pdb_id"], r["entity_id"], r["auth_asym_id"]))
with open(RES/"rf00177_instances.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
S = {}
S["total_pdb_entries"] = len(maps); S["total_mapped_instances"] = len(rows)
S["unique_pdb_entity_pairs"] = len({(r["pdb_id"], r["entity_id"]) for r in rows})
S["unique_entity_sequences_md5"] = len({r["sequence_md5"] for r in rows if r["sequence_md5"]})
S["unique_mapped_region_sequences"] = len({r["mapped_sequence"] for r in rows if r["mapped_sequence"]})
Ls = [r["mapped_length"] for r in rows]
S["mapped_length"] = {"min": min(Ls), "q1": statistics.quantiles(Ls, n=4)[0], "median": statistics.median(Ls), "q3": statistics.quantiles(Ls, n=4)[2], "max": max(Ls)}
S["mapped_length_histogram_100nt_bins"] = dict(sorted(collections.Counter(f"{(l//100)*100}-{(l//100)*100+99}" for l in Ls).items(), key=lambda kv: int(kv[0].split("-")[0])))
org = collections.Counter((r["organism"], r["tax_id"]) for r in rows)
S["n_distinct_organisms"] = len(org)
S["top_organisms_by_instances"] = [{"organism": o, "tax_id": t, "instances": n} for (o, t), n in org.most_common(25)]
S["instances_missing_organism"] = sum(1 for r in rows if not r["organism"])
multi = collections.Counter((r["pdb_id"], r["entity_id"]) for r in rows)
S["entities_with_multiple_mapped_chains"] = sum(1 for v in multi.values() if v > 1)
S["entries_with_multiple_mapped_instances"] = sum(1 for v in collections.Counter(r["pdb_id"] for r in rows).values() if v > 1)
S["instances_where_entity_has_more_chains_than_mapped"] = sum(1 for (p, e), v in multi.items() if v < max(r["n_chains_of_entity"] for r in rows if (r["pdb_id"], r["entity_id"]) == (p, e)))
S["instances_with_modified_nucleotides_in_entity"] = sum(1 for r in rows if r["n_modified_positions"])
S["instances_with_author_insertion_codes_at_boundaries"] = sum(1 for r in rows if r["map_auth_start_ins"] or r["map_auth_end_ins"])
json.dump(S, open(RES/"rf00177_instances_summary.json", "w"), indent=1)
print(json.dumps(S, indent=1))
# organelle-ish keyword scan for selection
print("\n--- instances with mito/chloro/plastid keywords or non-bacterial organism, length-sorted sample")
kw = [r for r in rows if any(k in (r["molecule_name"] + " " + str(r["organism"])).lower() for k in ["mito", "chloro", "plastid", "12s", "15s", "homo sapiens", "saccharomyces", "spinacia", "arabidopsis", "trypanosoma", "leishmania", "tetrahymena", "neurospora", "chlamydomonas", "plasmodium", "toxoplasma", "polytomella", "sus scrofa"])]
print(len(kw))
for r in sorted(kw, key=lambda r: (r["organism"] or "", r["mapped_length"])):
    print(f"  {r['pdb_id']} ent{r['entity_id']} {r['auth_asym_id']}/{r['label_asym_id']} map {r['map_start']}-{r['map_end']} ({r['mapped_length']}) entlen {r['entity_length']} | {r['organism']} | {r['molecule_name'][:50]}")
