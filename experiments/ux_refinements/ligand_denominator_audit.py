"""Audit: does the per-ligand denominator (entries with bound CCD) include PDB entries whose RNA has no observed nucleotide at the position?
Usage (from data/families): python3 ligand_denominator_audit.py RF00167. Reads instances.json, ligands.json and mappings/; writes nothing."""
import json, sys, os, collections
acc = sys.argv[1]; d = acc + '/'
inst = json.load(open(d + 'instances.json')); inst = inst['instances'] if isinstance(inst, dict) else inst
lig = json.load(open(d + 'ligands.json'))
byentry = collections.defaultdict(list)
for i in inst: byentry[i['pdb_id']].append(i['seq_name'])
obs = {}
def observable(sn):
    if sn in obs: return obs[sn]
    m = json.load(open(d + 'mappings/' + sn + '.json'))
    res = {r['polymer_seq_position']: r for r in m['residues']}
    s = set(c['rfam_consensus_pos'] for c in m['consensus_columns'] if c['state'] == 'match' and res.get(c['polymer_seq_position'], {}).get('has_coordinates'))
    obs[sn] = s; return s
ccd_entries = collections.defaultdict(set); contact = collections.defaultdict(lambda: collections.defaultdict(set))
for bm in lig['bound_molecules']:
    if not bm['in_default_track']: continue
    ccd_entries[bm['ccd_id']].add(bm['pdb_id'])
    for p in bm['consensus_positions']: contact[bm['ccd_id']][p].add(bm['pdb_id'])
n_marks = n_aff = 0; examples = []
for ccd, pos_map in contact.items():
    N = len(ccd_entries[ccd])
    for p, ents in pos_map.items():
        n_marks += 1
        unobs = [e for e in ccd_entries[ccd] if not any(p in observable(sn) for sn in byentry[e])]
        if unobs:
            n_aff += 1
            if len(examples) < 6: examples.append(f'{ccd}@{p}: {len(ents)}/{N} shown, {len(unobs)} entries unobservable -> {len(ents)}/{N-len(unobs)}')
print(acc, 'ligand marks', n_marks, 'with unobservable entries in denominator', n_aff); [print('  ', e) for e in examples]
