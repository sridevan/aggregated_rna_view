"""Predicted RNA ligand-binding pockets (fpocketR 1.3.6 / fpocket 4.0.3) for every chain instance of a family, mapped to Rfam consensus positions.

Protocol (validated in experiments/fpocketr_pilot):
  * Input per instance: single ligand-free RNA chain written from the cached PDBe updated mmCIF (polymer residues of the instance's
    label_asym_id, whole deposited construct, heavy atoms, highest-occupancy altloc, model 1, modified residues kept as ATOM records,
    author chain id and numbering). A residue map (pdb_chain, pdb_resnum, icode) -> (label_asym_id, label_seq_id) is kept per input.
  * Run: `python -m fpocketR -pdb <file> -c <chain> -m 3.0 -M 5.7 -i 42 -D 1.65 -A 3 -p 0.0 -dpi 10 -y -o <relative dir>` in the isolated
    x86-64 environment (tools/fpocketr-env via tools/micromamba); the output directory must be relative (fpocketR 1.3.6 strips a leading '/').
  * Pocket nucleotides = residues of fpocket's pockets/pocketN_atm.pdb (atoms contacted by the pocket's Voronoi vertices); no extra cut-off.
  * Mapping through the validated family mapping (InstanceIndex.resolve_unit on author ids); insert/flank residues kept with consensus null.
  * Aggregation: per consensus position, supporting instances (nucleotide in >= 1 qualifying pocket) / evaluable instances (run succeeded and
    the match residue has coordinates); one count per instance; three filters: all Score >= 0 (primary candidate), top-ranked only, Score >= 0.4.
"""
from __future__ import annotations
import csv
import datetime as _dt
import json
import logging
import os
import re
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional
from . import __version__
from .annotation_mapping import InstanceIndex
from .cache import Cache
from .config import PROJECT_ROOT
from .export import validate_json
from .pdbe import family_uses_chain_files, instance_cif

log = logging.getLogger("rna_explorer.pockets")
SCHEMA_VERSION = "1.0"
FPOCKET_PARAMS: dict[str, Any] = {"m": 3.0, "M": 5.7, "i": 42, "D": 1.65, "A": 3, "p": 0.0}
FILTERS = {"all_score_ge_0": "all pockets with fpocket Score >= 0 (fpocketR default qualityfilter 0.0)", "top1": "highest-ranked pocket only", "score_ge_0.4": "pockets with fpocket Score >= 0.4"}
DEFAULT_MICROMAMBA = PROJECT_ROOT / "tools/micromamba/bin/micromamba"
DEFAULT_ENV = PROJECT_ROOT / "tools/fpocketr-env"
TERMINAL_WINDOW = 2


# ----------------------------------------------------------------------------- input preparation (pilot policy)
def _pick_atoms(res):
    by: dict[str, list] = {}
    for a in res:
        if a.element.name.upper() in ("H", "D"):
            continue
        by.setdefault(a.name, []).append(a)
    return [sorted(v, key=lambda a: (-a.occ, a.altloc))[0] for v in by.values()], sum(len(v) - 1 for v in by.values())


def write_chain_pdb(cif_path: Path, inst: dict[str, Any], out_pdb: Path) -> dict[str, Any]:
    """Ligand-free single-chain PDB + residue map; validates the written file by re-reading it."""
    import gemmi
    st = gemmi.read_structure(str(cif_path)); st.setup_entities(); model = st[0]
    sub = model.get_subchain(inst["label_asym_id"])
    if len(sub) == 0:
        raise ValueError(f"subchain {inst['label_asym_id']} not found in {cif_path.name}")
    out = gemmi.Structure(); out.name = inst["pdb_id"].upper(); m = gemmi.Model("1"); ch = gemmi.Chain(inst["auth_asym_id"])
    rows = []; dropped = 0; n_atoms = 0
    for res in sub:
        atoms, d = _pick_atoms(res); dropped += d
        nr = gemmi.Residue(); nr.name = res.name; nr.seqid = res.seqid; nr.het_flag = "A"
        for a in atoms:
            na = gemmi.Atom(); na.name = a.name; na.element = a.element; na.pos = a.pos; na.occ = a.occ; na.b_iso = a.b_iso; na.altloc = "\0"; nr.add_atom(na); n_atoms += 1
        ch.add_residue(nr)
        rows.append({"pdb_chain": inst["auth_asym_id"], "pdb_resnum": res.seqid.num, "pdb_icode": res.seqid.icode.strip(), "resname": res.name, "label_asym_id": inst["label_asym_id"], "label_seq_id": res.label_seq,
                     "n_heavy_atoms": len(atoms), "altloc_copies_dropped": d})
    m.add_chain(ch); out.add_model(m); out.setup_entities()
    out_pdb.parent.mkdir(parents=True, exist_ok=True); out.write_pdb(str(out_pdb))
    chk = gemmi.read_structure(str(out_pdb)); back = [(r.seqid.num, r.seqid.icode.strip(), r.name, len(r)) for r in chk[0][inst["auth_asym_id"]]]
    ok = [(x["pdb_resnum"], x["pdb_icode"], x["resname"], x["n_heavy_atoms"]) for x in rows] == back
    return {"residues": rows, "n_residues": len(rows), "n_heavy_atoms": n_atoms, "altloc_copies_dropped": dropped, "roundtrip_ok": ok,
            "modified_residues": sorted({x["resname"] for x in rows if x["resname"] not in ("A", "C", "G", "U")}), "multi_char_chain_id": len(inst["auth_asym_id"]) > 1,
            "has_insertion_codes": any(x["pdb_icode"] for x in rows)}


# ----------------------------------------------------------------------------- execution
def run_fpocketr(pdb_path: Path, chain: str, work_dir: Path, out_dir_name: str, micromamba: Path = DEFAULT_MICROMAMBA, env: Path = DEFAULT_ENV, params: dict[str, Any] = FPOCKET_PARAMS, dpi: int = 10) -> dict[str, Any]:
    """Runs fpocketR in `work_dir` with a RELATIVE output directory (fpocketR 1.3.6 absolute-path bug). Returns status, timing and the output dir."""
    work_dir.mkdir(parents=True, exist_ok=True)
    local = work_dir / pdb_path.name
    if local.resolve() != pdb_path.resolve():
        shutil.copy2(pdb_path, local)
    clean = work_dir / f"{pdb_path.stem}_clean.pdb"
    if clean.exists():
        clean.unlink()
    cmd = [str(micromamba), "run", "-p", str(env), "python", "-W", "ignore", "-m", "fpocketR", "-pdb", local.name, "-c", chain,
           "-m", str(params["m"]), "-M", str(params["M"]), "-i", str(params["i"]), "-D", str(params["D"]), "-A", str(params["A"]), "-p", str(params["p"]), "-dpi", str(dpi), "-y", "-o", out_dir_name]
    envv = dict(os.environ); envv["MAMBA_ROOT_PREFIX"] = str(micromamba.parent.parent / "root")
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(work_dir), capture_output=True, text=True, env=envv)
    secs = time.time() - t0
    out_dir = work_dir / out_dir_name / f"{pdb_path.stem}_clean_out"
    (work_dir / out_dir_name).mkdir(exist_ok=True)
    (work_dir / out_dir_name / f"{pdb_path.stem}.fpocketR-log.txt").write_text(proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    csv_path = out_dir / f"{pdb_path.stem}_out_pocket_characteristics.csv"
    status = "ok" if proc.returncode == 0 and csv_path.exists() else "failed"
    err = None if status == "ok" else (proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else f"exit {proc.returncode}, no characteristics file")
    return {"status": status, "returncode": proc.returncode, "seconds": round(secs, 2), "out_dir": out_dir, "command": " ".join(cmd[5:]), "error": err}


# ----------------------------------------------------------------------------- parsing fpocket / fpocketR output
def read_info(path: Path) -> dict[int, dict[str, Any]]:
    out: dict[int, dict[str, Any]] = {}; cur = None
    if not path.exists():
        return out
    for line in open(path):
        m = re.match(r"Pocket (\d+) :", line)
        if m:
            cur = int(m.group(1)); out[cur] = {}; continue
        if cur and ":" in line:
            k, v = line.split(":", 1); v = v.strip()
            out[cur][k.strip()] = float(v) if re.fullmatch(r"-?[\d.]+(e[+-]?\d+)?", v) else v
    return out


def lining_atoms(atm_pdb: Path) -> dict[tuple, list[str]]:
    res: dict[tuple, list[str]] = defaultdict(list)
    for line in open(atm_pdb):
        if line.startswith(("ATOM", "HETATM")):
            res[(line[21].strip(), int(line[22:26]), line[26].strip(), line[17:20].strip())].append(line[12:16].strip())
    return res


def parse_run(out_dir: Path, name: str) -> list[dict[str, Any]]:
    """Pockets of one fpocketR run: characteristics CSV + fpocket info + lining atoms (residue keys are author ids)."""
    csv_path = out_dir / f"{name}_out_pocket_characteristics.csv"
    info = read_info(out_dir / f"{name}_clean_info.txt")
    pockets = []
    for r in csv.DictReader(open(csv_path)):
        n = int(r["Pocket"]); atm = out_dir / "pockets" / f"pocket{n}_atm.pdb"
        la = lining_atoms(atm) if atm.exists() else {}
        f = lambda k: float(r[k]) if r.get(k) not in (None, "", "nan") else None
        pockets.append({"rank": n, "score": f("Score"), "drug_score": f("Drug_score"), "volume_A3": f("Volume"), "n_alpha_spheres": int(float(r["a-sphere"])) if r.get("a-sphere") else None, "pocket_sasa_A2": f("SASA"),
                        "hydrophobicity_score": f("Hydrophobicity_score"), "polarity_score": f("Polarity_score"), "hydrophobic_density": f("Hydrophobic_density"), "apolar_alpha_sphere_proportion": f("Apolar_a-sphere_proportion"),
                        "npr1": f("Pocket_NPR1"), "npr2": f("Pocket_NPR2"), "shape": r.get("Pocket_shape") or None, "filter_pass": r.get("Filter") == "Pass",
                        "pocket_nt_author_csv": json.loads(r["PocketNT"]) if r.get("PocketNT") else [], "fpocket_info": info.get(n, {}),
                        "lining": [{"pdb_chain": k[0], "pdb_resnum": k[1], "pdb_icode": k[2], "resname": k[3], "atoms": sorted(v)} for k, v in sorted(la.items(), key=lambda x: (x[0][0], x[0][1], x[0][2]))]})
    return pockets


# ----------------------------------------------------------------------------- mapping
def map_pockets(pockets: list[dict[str, Any]], residue_rows: list[dict[str, Any]], idx: InstanceIndex, inst: dict[str, Any], incomplete_nts: set[int], n_residues: int) -> list[dict[str, Any]]:
    rmap = {(x["pdb_chain"], int(x["pdb_resnum"]), x["pdb_icode"]): x for x in residue_rows}
    seq_positions = {int(x["label_seq_id"]) for x in residue_rows}
    lo, hi = (min(seq_positions), max(seq_positions)) if seq_positions else (0, 0)
    out = []
    for p in pockets:
        nts = []
        for L in p["lining"]:
            m = rmap.get((L["pdb_chain"], L["pdb_resnum"], L["pdb_icode"]))
            ref = idx.resolve_unit(f"{inst['pdb_id'].upper()}|1|{L['pdb_chain']}|{L['resname']}|{L['pdb_resnum']}" + (f"|||{L['pdb_icode']}" if L["pdb_icode"] else ""))
            label_seq = int(m["label_seq_id"]) if m else None
            nts.append({"pdb_chain": L["pdb_chain"], "pdb_resnum": L["pdb_resnum"], "pdb_icode": L["pdb_icode"], "resname": L["resname"], "label_asym_id": inst["label_asym_id"], "label_seq_id": label_seq,
                        "resolver_label_seq_id": ref.label_seq_id, "rfam_consensus_pos": ref.rfam_consensus_pos, "alignment_state": ref.alignment_state, "mapping_status": ref.mapping_status,
                        "n_lining_atoms": len(L["atoms"]), "lining_atoms": L["atoms"], "incomplete_nucleotide": label_seq in incomplete_nts if label_seq is not None else None,
                        "terminal": label_seq is not None and (label_seq - lo < TERMINAL_WINDOW or hi - label_seq < TERMINAL_WINDOW)})
        cons = sorted({x["rfam_consensus_pos"] for x in nts if x["rfam_consensus_pos"] is not None})
        q = {k: v for k, v in p.items() if k != "lining"}
        q.update({"pocket_id": f"{inst['seq_name']}:pocket{p['rank']}", "nucleotides": nts, "n_nucleotides": len(nts), "consensus_positions": cons,
                  "n_non_consensus_nucleotides": sum(1 for x in nts if x["rfam_consensus_pos"] is None), "nt_match_csv": sorted(p["pocket_nt_author_csv"]) == sorted({x["pdb_resnum"] for x in nts}),
                  "map_consistent": all(x["label_seq_id"] == x["resolver_label_seq_id"] for x in nts), "unresolved_nucleotides": sum(1 for x in nts if x["mapping_status"] == "unresolved_identifier"),
                  "flags": sorted({f for f in ["includes_incomplete_nucleotide" if any(x["incomplete_nucleotide"] for x in nts) else None, "touches_chain_terminus" if any(x["terminal"] for x in nts) else None,
                                                "includes_non_consensus_nucleotide" if any(x["rfam_consensus_pos"] is None for x in nts) else None] if f})})
        out.append(q)
    return out


def qualifying(pockets: list[dict[str, Any]], flt: str) -> list[dict[str, Any]]:
    ps = [p for p in pockets if p.get("score") is not None and p["score"] >= 0]
    if flt == "top1":
        return ps[:1]
    if flt == "score_ge_0.4":
        return [p for p in ps if p["score"] >= 0.4]
    return ps


# ----------------------------------------------------------------------------- aggregation
def aggregate(instances: list[dict[str, Any]], idx: InstanceIndex, clen: int) -> dict[str, list[dict[str, Any]]]:
    """Per filter, per consensus position: supporting / evaluable instances, entries, distinct pockets, score summary. One count per instance."""
    out: dict[str, list[dict[str, Any]]] = {}
    for flt in FILTERS:
        positions = []
        for pos in range(1, clen + 1):
            supporting: set[str] = set(); evaluable: set[str] = set(); uneval = Counter(); entries: set[str] = set(); pockets: list[tuple[str, float]] = []
            for rec in instances:
                sn = rec["seq_name"]
                has_res = idx.match_residue(sn, pos) is not None
                if not has_res:
                    continue   # position deleted/uncovered in this chain: not part of the denominator
                if rec["status"] == "failed":
                    uneval["failed_run"] += 1; continue
                if not idx.observable(sn, pos):
                    uneval["no_coordinates"] += 1; continue
                evaluable.add(sn)
                for p in qualifying(rec["pockets"], flt):
                    if pos in p["consensus_positions"]:
                        supporting.add(sn); entries.add(rec["pdb_id"]); pockets.append((p["pocket_id"], p["score"]))
            scores = sorted(s for _, s in pockets)
            positions.append({"rfam_consensus_pos": pos, "n_supporting_instances": len(supporting), "n_evaluable_instances": len(evaluable), "frequency": round(len(supporting) / len(evaluable), 4) if evaluable else None,
                              "n_unevaluable_instances": sum(uneval.values()), "unevaluable_reasons": dict(uneval), "n_supporting_entries": len(entries), "n_distinct_pockets": len(pockets),
                              "pocket_ids": sorted(p for p, _ in pockets), "score_summary": {"min": scores[0], "median": scores[len(scores) // 2], "max": scores[-1]} if scores else None})
        out[flt] = positions
    return out


# ----------------------------------------------------------------------------- build
def build_pockets(acc: str, cache_dir: Path, out_dir: Path, experiment_dir: Path, instances: Optional[list[str]] = None, micromamba: Path = DEFAULT_MICROMAMBA, env: Path = DEFAULT_ENV, write: bool = True,
                  reuse_existing_runs: bool = True) -> dict[str, Any]:
    fam_dir = out_dir / "families" / acc
    idx = InstanceIndex(fam_dir); cache = Cache(cache_dir)
    chain_files = family_uses_chain_files(fam_dir)
    rows = [r for r in idx.instances if not instances or r["seq_name"] in set(instances)]
    inputs_dir = experiment_dir / "inputs"; work_dir = experiment_dir / "work"; runs_name = "runs"
    sasa_status: dict[str, set[int]] = {}
    sasa_file = fam_dir / "sasa.json"
    if sasa_file.exists():
        for m in json.load(open(sasa_file))["instances"]:
            sasa_status[m["seq_name"]] = {n["label_seq_id"] for n in m["nucleotides"] if n["status"] == "incomplete"}
    generated = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    records = []
    t_all = time.time()
    for r in rows:
        sn = r["seq_name"]; rec: dict[str, Any] = {"seq_name": sn, "pdb_id": r["pdb_id"], "entity_id": r["entity_id"], "label_asym_id": r["label_asym_id"], "auth_asym_id": r["auth_asym_id"], "warnings": []}
        try:
            pdb_path = inputs_dir / f"{sn}.pdb"
            prep = write_chain_pdb(instance_cif(cache, r, chain_files), r, pdb_path)
            with open(inputs_dir / f"{sn}_residue_map.csv", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(prep["residues"][0].keys())); w.writeheader(); w.writerows(prep["residues"])
            rec["input"] = {k: v for k, v in prep.items() if k != "residues"}
            if not prep["roundtrip_ok"]:
                rec["warnings"].append("PDB round-trip mismatch")
            if prep["multi_char_chain_id"]:
                rec["warnings"].append("multi-character author chain id truncated in PDB")
            out_existing = work_dir / runs_name / f"{sn}_clean_out" / f"{sn}_out_pocket_characteristics.csv"
            if reuse_existing_runs and out_existing.exists():
                run = {"status": "ok", "returncode": 0, "seconds": None, "out_dir": out_existing.parent, "command": "(reused existing run)", "error": None}
            else:
                run = run_fpocketr(pdb_path, r["auth_asym_id"], work_dir, runs_name, micromamba, env)
            rec["run"] = {k: (str(v.relative_to(experiment_dir)) if isinstance(v, Path) else v) for k, v in run.items()}
            if run["status"] != "ok":
                rec.update({"status": "failed", "failure_reason": f"fpocketR execution: {run['error']}", "pockets": []})
            else:
                pockets = map_pockets(parse_run(run["out_dir"], sn), prep["residues"], idx, r, sasa_status.get(sn, set()), prep["n_residues"])
                if any(not p["nt_match_csv"] for p in pockets):
                    rec["warnings"].append("PocketNT differs from lining-atom residues")
                if any(not p["map_consistent"] for p in pockets):
                    rec["warnings"].append("residue map and resolver disagree")
                rec.update({"status": "pockets" if pockets else "no_pockets", "pockets": pockets})
        except Exception as e:   # input preparation / parsing failure
            log.error("%s failed: %s", sn, e)
            rec.update({"status": "failed", "failure_reason": f"{type(e).__name__}: {e}", "pockets": []})
        rec["n_pockets"] = len(rec["pockets"]); rec["n_pockets_by_filter"] = {flt: len(qualifying(rec["pockets"], flt)) for flt in FILTERS}
        log.info("%s %s %d pockets", sn, rec["status"], rec["n_pockets"])
        records.append(rec)
    total_secs = round(time.time() - t_all, 1)
    tracks = aggregate(records, idx, idx.clen)
    mp = idx.mapping(idx.instances[0]["seq_name"])["provenance"]
    method = {"software": {"fpocketR": "1.3.6", "fpocket": "4.0.3 (conda-forge h8a118fb_2, x86_64 under Rosetta 2)", "python": "3.11.17 (x86_64)", "pymol-open-source": "2.5.0", "rna_explorer": __version__},
              "environment": str(env), "fpocket_parameters": FPOCKET_PARAMS, "fpocketR_options": "-dpi 10 -y, relative --out (fpocketR 1.3.6 strips a leading '/')",
              "input": "ligand-free single RNA chain (polymer residues of label_asym_id, whole deposited construct, heavy atoms, highest-occupancy altloc, model 1, modified residues kept) from the PDBe updated mmCIF",
              "pocket_nucleotides": "residues of fpocket pockets/pocketN_atm.pdb (atoms contacted by the pocket's Voronoi vertices); no additional distance cut-off",
              "mapping": "author (chain, number, insertion code) -> label_seq_id via the input residue map, verified against InstanceIndex.resolve_unit; consensus positions from the validated family mapping",
              "filters": FILTERS, "terminal_window_residues": TERMINAL_WINDOW, "incomplete_nucleotides_source": "sasa.json status 'incomplete'" if sasa_status else None,
              "volume_note": "fpocket pocket volume is a Monte-Carlo estimate (about +/-3 % between runs); scores, alpha-sphere counts and lining atoms are deterministic"}
    prov = {"rfam_acc": acc, "rfam_release": mp.get("rfam_release"), "cm_checksum": mp.get("cm_checksum"), "cm_clen": idx.clen, "mapping_pipeline_version": mp.get("pipeline_version"), "mapping_generated": mp.get("generated"),
            "generated": generated, "experiment_dir": str(experiment_dir), "subset": sorted(instances) if instances else None}
    status_counts = Counter(rec["status"] for rec in records)
    pockets_json = {"schema_version": SCHEMA_VERSION, "rfam_acc": acc, "clen": idx.clen, "method": method,
                    "counts": {"n_instances": len(records), "status_counts": dict(status_counts), "n_pockets_total": sum(rec["n_pockets"] for rec in records),
                               "n_pockets_by_filter": {flt: sum(rec["n_pockets_by_filter"][flt] for rec in records) for flt in FILTERS}, "total_seconds": total_secs},
                    "instances": records, "provenance": prov}
    tracks_json = {"schema_version": SCHEMA_VERSION, "rfam_acc": acc, "clen": idx.clen, "track_id": "pocket_occurrence", "category": "Predicted binding pockets",
                   "measurement": {"name": "Pocket occurrence", "definition": "fraction of evaluable RNA chain instances whose nucleotide at the consensus position lines at least one qualifying predicted pocket",
                                   "evaluable": "fpocketR run succeeded and the match residue has coordinates; deleted/uncovered positions and failed runs are outside the denominator"},
                   "filters": FILTERS, "n_instances": len(records), "status_counts": dict(status_counts), "tracks": tracks, "method": method, "provenance": prov}
    if write:
        validate_json(pockets_json, "pockets"); validate_json(tracks_json, "pocket_tracks")
        with open(fam_dir / "pockets.json", "w") as f: json.dump(pockets_json, f, indent=1)
        with open(fam_dir / "pocket_tracks.json", "w") as f: json.dump(tracks_json, f, indent=1)
    return {"acc": acc, "n_instances": len(records), "status_counts": dict(status_counts), "n_pockets_total": pockets_json["counts"]["n_pockets_total"], "total_seconds": total_secs, "pockets": pockets_json, "tracks": tracks_json}
