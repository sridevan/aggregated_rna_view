"""Shared helpers for the RF00177 POC (minimal mmCIF reading; no external deps)."""
import shlex

def _tokens(line):
    try: return shlex.split(line, posix=True)
    except ValueError: return line.split()

def cif_category(path, category):
    """Return list of dicts for an mmCIF category written as loop_ or as key-value pairs.
    Handles quoted tokens and ;-delimited multi-line values. Streams the file; non-target lines are skipped fast."""
    prefix = category + "."
    rows, keys, pending, cur = [], [], [], {}
    in_loop, target_loop, collecting = False, False, False
    with open(path) as f:
        it = iter(f)
        for line in it:
            if line.startswith("loop_"):
                in_loop, target_loop, collecting, keys, pending = True, False, True, [], []; continue
            if line.startswith("#"):
                if (target_loop and rows) or cur: break
                in_loop, target_loop, collecting = False, False, False; continue
            if line.startswith("_"):
                name = line.split(None, 1)[0]
                if in_loop and collecting:
                    if name.startswith(prefix): keys.append(name[len(prefix):]); target_loop = True
                    continue
                if not in_loop and name.startswith(prefix):
                    key = name[len(prefix):]; rest = line[len(name):].strip()
                    if rest == "":
                        nxt = next(it)
                        if nxt.startswith(";"):
                            val = nxt[1:].rstrip("\n")
                            for l2 in it:
                                if l2.startswith(";"): break
                                val += "\n" + l2.rstrip("\n")
                        else: val = nxt.strip()
                    else:
                        t = _tokens(rest); val = t[0] if t else rest
                    cur[key] = val
                continue
            # data line
            if in_loop and target_loop:
                collecting = False
                if line.startswith(";"):
                    val = line[1:].rstrip("\n")
                    for l2 in it:
                        if l2.startswith(";"): break
                        val += "\n" + l2.rstrip("\n")
                    pending.append(val)
                else:
                    pending.extend(_tokens(line.strip()))
                while keys and len(pending) >= len(keys):
                    rows.append(dict(zip(keys, pending[:len(keys)]))); pending = pending[len(keys):]
            elif in_loop and collecting:
                collecting = False   # data of a non-target loop; skip quickly
    if cur: rows.append(cur)
    return rows

def atom_site_observed(path, label_asym_id, entity_id):
    """Set of label_seq_id (int) with at least one atom for the given polymer asym (streams atom_site)."""
    keys, seen, a = [], set(), None
    with open(path) as f:
        for line in f:
            if line.startswith("_atom_site."): keys.append(line.split()[0].split(".",1)[1]); continue
            if keys and (line.startswith("ATOM") or line.startswith("HETATM")):
                if a is None:
                    a = keys.index("label_asym_id"); s = keys.index("label_seq_id"); e = keys.index("label_entity_id")
                t = line.split()
                if t[a] == label_asym_id and t[e] == str(entity_id) and t[s] not in (".", "?"): seen.add(int(t[s]))
            elif keys and a is not None and line.startswith("#"): break
    return seen

def entity_source(path, entity_id):
    """Organism for an entity from the mmCIF source categories."""
    for cat, org, tax in [("_entity_src_nat", "pdbx_organism_scientific", "pdbx_ncbi_taxonomy_id"),
                          ("_entity_src_gen", "pdbx_gene_src_scientific_name", "pdbx_gene_src_ncbi_taxonomy_id"),
                          ("_pdbx_entity_src_syn", "organism_scientific", "ncbi_taxonomy_id")]:
        for r in cif_category(path, cat):
            if r.get("entity_id") == str(entity_id):
                return {"cif_src_category": cat, "cif_organism": r.get(org), "cif_tax_id": r.get(tax)}
    return {"cif_src_category": None, "cif_organism": None, "cif_tax_id": None}
