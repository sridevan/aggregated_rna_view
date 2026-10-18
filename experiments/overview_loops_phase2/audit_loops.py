"""Read-only audit of BGSU loop membership for Family Overview Phase 2 (data/families/RF00167/loops.json + mappings + annotations).
Checks that every loop's complete nucleotide membership is recoverable and maps to the instance identifiers used by the 2D and 3D viewers."""
import json, collections, glob
F = 'data/families/RF00167/'
L = json.load(open(F + 'loops.json'))
loops = L['loops']
maps = {f.split('/')[-1][:-5]: json.load(open(f)) for f in glob.glob(F + 'mappings/*.json')}
res = {sn: {r['polymer_seq_position']: r for r in m['residues']} for sn, m in maps.items()}
c = collections.Counter(); issues = []
types = collections.Counter(l['loop_type'] for l in loops)
for l in loops:
    nts = l['nucleotides']
    c['loops'] += 1; c['nts'] += len(nts)
    if len(nts) != l['n_nucleotides']: issues.append(('n_nucleotides mismatch', l['loop_id'], len(nts), l['n_nucleotides']))
    if sum(s['n_nucleotides'] for s in l['strands']) != len(nts): issues.append(('strand total mismatch', l['loop_id']))
    # strands partition the nucleotide list in order
    for s in l['strands']:
        a, b = s['nucleotide_indices']
        if b - a + 1 != s['n_nucleotides']: issues.append(('strand index span', l['loop_id']))
    if l['multi_chain']: c['multi_chain'] += 1
    for n in nts:
        c['status:' + str(n.get('mapping_status'))] += 1; c['state:' + str(n.get('alignment_state'))] += 1
        c['coords:' + str(n.get('has_coordinates'))] += 1
        sn = n.get('seq_name')
        if sn in res:
            r = res[sn].get(n['label_seq_id'])
            if not r: issues.append(('label_seq_id not in mapping', l['loop_id'], sn, n['label_seq_id'])); continue
            if str(r['author_residue_number']) != str(n['author_residue_number']) or (r['author_insertion_code'] or '') != (n['insertion_code'] or ''): issues.append(('author id mismatch', l['loop_id'], sn, n['label_seq_id']))
            if r['residue_name'] != n['residue_name']: issues.append(('residue name mismatch', l['loop_id'], sn))
            if r['alignment_state'] != n['alignment_state'] or (r['rfam_consensus_pos'] if r['alignment_state'] == 'match' else None) != n['rfam_consensus_pos']: issues.append(('mapping mismatch', l['loop_id'], sn, n['label_seq_id']))
            if bool(r['has_coordinates']) != bool(n['has_coordinates']): issues.append(('coords mismatch', l['loop_id'], sn))
        else:
            c['nts_outside_family_chains'] += 1
    # consecutive label_seq_id within each strand (per chain)?
    for s in l['strands']:
        a, b = s['nucleotide_indices']; ids = [n['label_seq_id'] for n in nts[a:b + 1]]
        if any(y - x != 1 for x, y in zip(ids, ids[1:]) if x is not None and y is not None): c['strands_with_label_gaps'] += 1
# per-instance cross-check: annotation residues' loop ids == loops.json membership
mism = 0
for f in glob.glob(F + 'annotations/*.json'):
    a = json.load(open(f)); sn = a['seq_name']
    per_res = {(r['label_seq_id'], (x['loop_id'] if isinstance(x, dict) else x)) for r in a['residues'] for x in (r.get('loops') or [])}
    from_loops = {(n['label_seq_id'], l['loop_id']) for l in loops for n in l['nucleotides'] if n.get('seq_name') == sn}
    if per_res != from_loops: mism += 1
# overlaps within the same instance and track (HL / IL / JL = J3+J4)
cat = lambda t: 'JL' if t in ('J3', 'J4') else t
ov = collections.Counter()
by_inst = collections.defaultdict(lambda: collections.defaultdict(set))
for l in loops:
    for n in l['nucleotides']:
        if n.get('seq_name') and n['alignment_state'] == 'match': by_inst[(n['seq_name'], cat(l['loop_type']))][n['rfam_consensus_pos']].add(l['loop_id'])
examples = []
for (sn, t), d in by_inst.items():
    for p, ids in d.items():
        if len(ids) > 1: ov[t] += 1; examples.append((sn, t, p, sorted(ids))) if len(examples) < 6 else None
out = {'counts': dict(c), 'loop_types': dict(types), 'issues': issues[:20], 'n_issues': len(issues), 'annotation_vs_loops_mismatched_instances': mism,
       'same_track_overlaps_by_type': dict(ov), 'overlap_examples': examples}
print(json.dumps(out, indent=1))
json.dump(out, open('experiments/overview_loops_phase2/audit_loops.json', 'w'), indent=1)
