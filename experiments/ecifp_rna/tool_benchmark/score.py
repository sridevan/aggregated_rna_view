"""Pilot benchmark scoring: tool contacts (run_tools.py) against the literature reference (reference/reference_contacts.json) and against
each other.

Classes compared (tool type -> class): hbond, ionic -> polar; stacking -> stacking; cation_pi, anion_pi -> charge_pi;
water_mediated -> water; metal_mediated -> metal. Reference types: hbond, electrostatic -> polar; stacking; water_mediated -> water;
metal_mediated -> metal. weak_hbond, hydrophobic, other_pi and other are reported but not scored against the reference (papers rarely
list them).

Per tool and class:
- residue-level recall: reference (copy, residue, class) found by the tool with the same class at that residue;
- atom-level recall (polar only, where the reference names an RNA atom): the tool has a polar contact involving that RNA atom;
- 'unsupported': tool (copy, residue, class) not in the reference, for complexes whose reference is marked complete for that class
  (CURATION_NOTES); reported as counts, not as precision, because papers do not claim to list every contact;
- agreement between tools: residue-level Jaccard per class over all copies.

usage: score.py <records.json> <reference.json> <out_json>
"""
import json, sys, collections, itertools, re

REC, REF, OUT = sys.argv[1:4]
CLASS = {"hbond": "polar", "ionic": "polar", "stacking": "stacking", "cation_pi": "charge_pi", "anion_pi": "charge_pi",
         "water_mediated": "water", "metal_mediated": "metal"}
REFCLASS = {"hbond": "polar", "electrostatic": "polar", "stacking": "stacking", "water_mediated": "water", "metal_mediated": "metal"}
TOOLS = ["arpeggio_api", "arpeggio_local", "plip", "plip_H", "fingernat_H"]


def atom_name(s):
    if not s:
        return None
    s = s.strip().replace("*", "'")
    m = re.fullmatch(r"(?:2'-?OH|2'OH|O2')", s, re.I)
    if m:
        return "O2'"
    m = re.fullmatch(r"(OP1|OP2|O1P|O2P|N[1-7]|O[2-6]|O[2-5]'|C[1-8]'?)", s)
    return {"O1P": "OP1", "O2P": "OP2"}.get(s, s) if m else None


def main():
    recs = json.load(open(REC))["records"]
    ref = json.load(open(REF))
    copies = sorted({r["copy"] for r in recs})
    key = lambda pdb, ch, num, ccd: next((c for c in copies if c.lower().startswith(f"{pdb.lower()}_") and c.endswith(f"_{ccd}")
                                          and (ch is None or c.split("_")[1] == ch)), None)
    tool_sets = collections.defaultdict(set)      # (tool, class) -> {(copy, res)}
    tool_atoms = collections.defaultdict(set)     # (tool) -> {(copy, res, atom)} polar
    for r in recs:
        c = CLASS.get(r["type"])
        if c:
            tool_sets[(r["tool"], c)].add((r["copy"], r["rna"]))
            if c == "polar":
                for a in r["rna_atoms"]:
                    tool_atoms[r["tool"]].add((r["copy"], r["rna"], a))
    refset = collections.defaultdict(set); refatoms = set(); unmatched = []
    for x in ref:
        c = REFCLASS.get(x["type"])
        cp = key(x["pdb"], x.get("ligand_chain"), x.get("rna_resnum"), x["ligand_ccd"])
        if not c or cp is None or x.get("rna_resnum") in (None, ""):
            unmatched.append(x); continue
        res = f"{x['rna_chain']}:{int(x['rna_resnum'])}"
        refset[c].add((cp, res))
        a = atom_name(x.get("rna_atom_or_group"))
        if c == "polar" and a:
            refatoms.add((cp, res, a))
    out = dict(copies=copies, reference_counts={c: len(v) for c, v in refset.items()}, reference_atoms=len(refatoms),
               reference_unmatched=len(unmatched), recall={}, atom_recall={}, unsupported={}, jaccard={}, per_copy={})
    for t in TOOLS:
        out["recall"][t] = {}
        for c, R in refset.items():
            T = tool_sets[(t, c)]
            # a tool with no output for a copy (no_data) is scored only on copies it has
            have = {cp for (tt, _), s in tool_sets.items() if tt == t for cp, _ in s}
            Rc = {x for x in R if x[0] in have}
            out["recall"][t][c] = dict(found=len(Rc & T), of=len(Rc), missed=sorted(Rc - T))
            out["unsupported"].setdefault(t, {})[c] = len(T - R)
        A = tool_atoms[t]
        out["atom_recall"][t] = dict(found=len(refatoms & A), of=len(refatoms), missed=sorted(refatoms - A))
    for c in ("polar", "stacking", "water", "metal", "charge_pi"):
        out["jaccard"][c] = {f"{a}|{b}": round(len(tool_sets[(a, c)] & tool_sets[(b, c)]) / max(1, len(tool_sets[(a, c)] | tool_sets[(b, c)])), 2)
                             for a, b in itertools.combinations(TOOLS, 2)}
    for cp in copies:
        out["per_copy"][cp] = {t: {c: sorted(r for (k, r) in tool_sets[(t, c)] if k == cp) for c in ("polar", "stacking", "water", "metal", "charge_pi")}
                               for t in TOOLS}
        out["per_copy"][cp]["reference"] = {c: sorted(r for (k, r) in refset[c] if k == cp) for c in refset}
    json.dump(out, open(OUT, "w"), indent=1)
    print("reference:", out["reference_counts"], "atoms", len(refatoms), "unmatched", len(unmatched))
    for t in TOOLS:
        print(f"{t:<15} " + "  ".join(f"{c} {v['found']}/{v['of']}" for c, v in out["recall"][t].items())
              + f"  | atoms {out['atom_recall'][t]['found']}/{out['atom_recall'][t]['of']}  | unsupported {out['unsupported'][t]}")
    for c, d in out["jaccard"].items():
        print(c, d)


main()
