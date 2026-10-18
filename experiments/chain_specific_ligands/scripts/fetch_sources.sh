#!/usr/bin/env bash
# Reproducible retrieval of every response used in this investigation (retrieved 2026-10-15). Files land in ../cache/.
set -u; cd "$(dirname "$0")/../cache"
get() { local u="$1"; local n; n=$(echo "$u" | sed 's#https://##; s#[/?=&|]#_#g'); [ -f "$n.json" ] || curl -s -m 90 -o "$n.json" -w "%{http_code} $u\n" "$u"; }
for p in 1y26 3rkf 9ndd 9nbh; do P=$(echo $p | tr a-z A-Z)
  get "https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecules/$p"
  get "https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/$p"
  get "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/$p"
  get "https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/$p"
  get "https://data.rcsb.org/rest/v1/core/entry/$P"
  for bm in $(python3 -c "import json;print(' '.join(b['bm_id'] for b in json.load(open('www.ebi.ac.uk_pdbe_graph-api_pdb_bound_molecules_$p.json'))['$p']))"); do
    get "https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecule_interactions/$p/$bm"; done
  for c in $(python3 -c "import json;print(' '.join(sorted({r['label_asym_id'] for r in json.load(open('../../../data/families/RF00167/instances.json')) if r['pdb_id']=='$p'})))"); do
    get "https://data.rcsb.org/rest/v1/core/polymer_entity_instance/$P/$c"; done
  for a in $(python3 -c "import json;print(' '.join(m['struct_asym_id'] for m in json.load(open('www.ebi.ac.uk_pdbe_api_pdb_entry_ligand_monomers_$p.json'))['$p']))"); do
    get "https://data.rcsb.org/rest/v1/core/nonpolymer_entity_instance/$P/$a"; done
done
for q in 1y26/X/90 3rkf/A/91 3rkf/B/91 3rkf/C/91 3rkf/D/91 3rkf/A/101 9ndd/A/201 9ndd/C/201 9nbh/A/201; do get "https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_ligand_interactions/$q"; done
# endpoints that returned 404 on 2026-10-15 (kept for the record): graph-api/pdb/bound_ligand_interactions/{pdb}/{bm}, api/pdb/entry/binding_sites/{pdb},
#   api/v2/pdb/entry/bound_molecules/{pdb}, graph-api/pdb/ligand_sites/{pdb}, graph-api/residue_mapping/{pdb}, graph-api/pdb/bound_excluded_molecules/{pdb}
curl -s -m 60 -o pdbe_graph_api_openapi.json https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/openapi.json
curl -s -m 60 -o rcsb_schema_polymer_entity_instance.json https://data.rcsb.org/rest/v1/schema/polymer_entity_instance
curl -s -m 60 -o rcsb_schema_nonpolymer_entity_instance.json https://data.rcsb.org/rest/v1/schema/nonpolymer_entity_instance
curl -s -m 60 -o models_rcsb_1y26_ADE_r5.cif "https://models.rcsb.org/v1/1y26/ligand?label_comp_id=ADE&radius=5&encoding=cif"
