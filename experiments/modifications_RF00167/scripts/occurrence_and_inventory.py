"""Provisional (NOT production) modification-occurrence counts and the RF00167 component inventory, from audit/polymer_positions.csv.
Policy investigated (B8): per consensus match position p,
  positive = instance whose deposited polymer residue at p (match state) is not A/C/G/U
  negative = instance whose deposited polymer residue at p is A/C/G/U (coordinates not required: identity comes from the deposited sequence)
  missing  = residue at p whose identity cannot be established (ambiguous/UNK/N, microheterogeneity, sequence-vs-coordinate conflict)
  no residue at p (deleted/uncovered in the validated mapping) = not in the denominator (no nucleotide to assess)
Each instance counts at most once per position. Writes audit/provisional_occurrence.csv, audit/inventory.json."""
import csv, json, collections, pathlib
HERE = pathlib.Path(__file__).resolve().parents[1]; ROOT = HERE.parents[1]
STD = {"A", "C", "G", "U"}; AMBIG = {"N", "UNK", "X"}
rows = list(csv.DictReader(open(HERE / "audit/polymer_positions.csv")))
nakb = json.load(open(HERE / "raw/nakb/modified_to_change_data.json"))
pdbe = json.load(open(HERE / "raw/pdbe/batch_all_entries.json"))["body"]
clen = 102
inst_all = sorted({r["seq_name"] for r in rows})

def status(r):
    c = r["comp_poly_seq_scheme"]
    if "|" in c or c in AMBIG or r["hetero"] == "True": return "missing"
    if r["comp_atom_site"] and r["comp_atom_site"] != c: return "missing"
    return "negative" if c in STD else "positive"

per = {p: collections.defaultdict(dict) for p in range(1, clen + 1)}
for r in rows:
    if r["alignment_state"] != "match": continue
    p = int(r["rfam_consensus_pos"]); per[p][r["seq_name"]][r["label_seq_id"]] = r
out = []
for p in range(1, clen + 1):
    pos, neg, mis, unres, ccd, entries = 0, 0, 0, 0, collections.Counter(), set()
    for sn, res in per[p].items():
        assert len(res) == 1, (p, sn)   # one match residue per instance per position
        r = next(iter(res.values())); s = status(r)
        if s == "positive": pos += 1; ccd[r["comp_poly_seq_scheme"]] += 1; entries.add(r["pdb_id"])
        elif s == "negative": neg += 1
        else: mis += 1
        if r["has_coordinates_mmcif"] == "False": unres += 1
    den = pos + neg
    out.append(dict(rfam_consensus_pos=p, n_positive=pos, n_negative=neg, n_missing=mis, denominator=den, frequency=round(pos / den, 4) if den else None,
                    n_instances_no_residue_at_position=len(inst_all) - len(per[p]), n_in_denominator_without_coordinates=unres,
                    n_positive_entries=len(entries), ccd_counts=json.dumps(dict(ccd)) if ccd else ""))
with open(HERE / "audit/provisional_occurrence.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)

ns = [r for r in rows if r["comp_poly_seq_scheme"] not in STD]
inv = {}
for c in sorted({r["comp_poly_seq_scheme"] for r in ns}):
    rs = [r for r in ns if r["comp_poly_seq_scheme"] == c]; nk = nakb.get(c)
    inv[c] = dict(n_residue_observations=len(rs), n_rna_instances=len({r["seq_name"] for r in rs}), n_pdb_entries=len({r["pdb_id"] for r in rs}),
                  instances=sorted({r["seq_name"] for r in rs}),
                  states=dict(collections.Counter(r["alignment_state"] for r in rs)),
                  rfam_positions=sorted({int(r["rfam_consensus_pos"]) for r in rs if r["alignment_state"] == "match"}),
                  terminal=dict(collections.Counter(r["is_terminal"] or "internal" for r in rs)),
                  with_coordinates=sum(r["has_coordinates_mmcif"] == "True" for r in rs),
                  in_pdbe=sum(r["pdbe_n_records"] != "0" for r in rs),
                  in_nakb=nk is not None, nakb_standard_base=nk and nk["standard_base"], nakb_par_comp_id=nk and nk["pdb"]["par_comp_id"],
                  nakb_change_locations=nk and sorted({x.get("change_location") for x in nk["changes"]} | {x.get("change_location_2") for x in nk["changes"] if x.get("change_location_2")}),
                  nakb_name=nk and nk["pdb"]["name"], nakb_chem_comp_type=nk and nk["pdb"]["chem_comp_type"])
summ = dict(
    n_entries=len({r["pdb_id"] for r in rows}), n_instances=len(inst_all), n_polymer_residues=len(rows),
    n_residues_with_coordinates=sum(r["has_coordinates_mmcif"] == "True" for r in rows),
    n_instances_with_nonstandard=len({r["seq_name"] for r in ns}), n_entries_with_nonstandard=len({r["pdb_id"] for r in ns}),
    n_distinct_ccd=len(inv), n_nonstandard_residue_observations=len(ns),
    n_on_match=sum(r["alignment_state"] == "match" for r in ns), n_on_insert=sum(r["alignment_state"] == "insert" for r in ns),
    n_local_end=sum(r["alignment_state"] == "local_end" for r in ns), n_outside_mapped_region=sum(r["alignment_state"] == "not_submitted_flank" for r in ns),
    distinct_positions=sorted({int(r["rfam_consensus_pos"]) for r in ns if r["alignment_state"] == "match"}),
    ccd_reported_by_pdbe=sorted({x["chem_comp_id"] for v in pdbe.values() for x in v}), n_pdbe_records=sum(len(v) for v in pdbe.values()),
    ccd_in_nakb=sorted(c for c in inv if inv[c]["in_nakb"]), ccd_mmcif_only=sorted(c for c in inv if not inv[c]["in_pdbe"]),
    match_positions_status_counts=dict(collections.Counter(status(r) for r in rows if r["alignment_state"] == "match")),
    all_positions_status_counts=dict(collections.Counter(status(r) for r in rows)),
    alignment_state_counts=dict(collections.Counter(r["alignment_state"] for r in rows)))
json.dump(dict(summary=summ, inventory=inv), open(HERE / "audit/inventory.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
for c, v in inv.items(): print(c, {k: v[k] for k in ("n_residue_observations", "n_rna_instances", "n_pdb_entries", "rfam_positions", "states", "terminal", "with_coordinates", "in_pdbe", "nakb_standard_base", "nakb_par_comp_id", "nakb_change_locations")})
print([o for o in out if o["n_positive"]])
dens = [o["denominator"] for o in out if o["denominator"]]
print("denominator range over covered positions", min(dens), max(dens), "positions with any residue", len(dens), "positions with unresolved residues counted in denominator", sum(1 for o in out if o["n_in_denominator_without_coordinates"]))
