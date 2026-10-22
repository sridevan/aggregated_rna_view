"""Near base pairs for the SS_cons comparison (investigation only; see NEAR_CHECK.md).

BGSU's published per-entry base-pair CSVs contain no near ('n'-prefixed) classes. fr3d-python (the BGSU annotation code) writes near
pairs only when the category 'near' is requested, so they are regenerated locally here:

  fr3d-python NA_pairwise_interactions.py -c basepair,near --chain <family chains>   (pinned commit given with --fr3d)

on the PDBe updated mmCIF (model 1). Each mmCIF is downloaded into a temporary scratch directory, annotated and deleted (never stored).
Endpoints are joined to Rfam positions exactly like the production pipeline: (pdb, author chain) -> instance; (author number, insertion
code) -> validated mapping residue -> match-state consensus position (insertions / flanks are kept as None).

Outputs in followup/near/:
  local_pairs/<seq_name>.tsv     every local base pair (true and near) of the chain on consensus positions
  concordance.csv                per chain: local true base pairs vs BGSU true base pairs (on unordered consensus position pairs + class)
  ss_cons_near.csv               every evaluable SS_cons pair observation without a true BGSU base pair: local near class(es) between i and j,
                                 near partners of i and j elsewhere, joined with the geometry category of followup/unpaired_geometry.csv
  near_summary.json
Usage: python3 experiments/sscons_vs_3d/near_check.py --fr3d <fr3d-python checkout> --python <interpreter> [--workers 4] [ACC ...]
"""
import argparse, csv, json, shutil, subprocess, sys, tempfile, urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).parent; OUT = HERE / "followup" / "near"
sys.path.insert(0, str(HERE))
from followup import load, evaluable  # noqa: E402
from analyse import ss_pairs  # noqa: E402


def flip(cls):
    import re
    m = re.fullmatch(r"(n?)([ct])([WHS])([WHS])(.*)", cls)
    return f"{m[1]}{m[2]}{m[4]}{m[3]}{m[5]}" if m else cls


def run_entry(pdb, chains, fr3d, python, tmp_root):
    """annotate one entry (only the family chains) locally; returns rows of the basepair output or raises"""
    d = Path(tempfile.mkdtemp(dir=tmp_root))
    try:
        cif = d / f"{pdb.upper()}.cif"
        with urllib.request.urlopen(f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb.lower()}_updated.cif", timeout=600) as r, open(cif, "wb") as f:
            shutil.copyfileobj(r, f)
        cmd = [python, str(Path(fr3d) / "fr3d/classifiers/NA_pairwise_interactions.py"), "-i", str(d), "-o", str(d), "-c", "basepair,near",
               "--chain", ",".join(sorted(chains)), cif.name]
        p = subprocess.run(cmd, cwd=fr3d, env={"PYTHONPATH": str(fr3d), "PATH": "/usr/bin:/bin"}, capture_output=True, text=True, timeout=3600)
        outs = list(d.glob("*basepair*.txt"))
        if not outs:
            raise RuntimeError(f"no output: {p.stderr[-300:]}")
        return [line.rstrip("\n").split("\t") for f in outs for line in open(f) if line.strip()]
    finally:
        shutil.rmtree(d, ignore_errors=True)


def unit(u):
    f = u.split("|") + [""] * 9
    return {"pdb": f[0].upper(), "model": f[1], "chain": f[2], "name": f[3], "num": int(f[4]), "alt": f[6], "ins": f[7], "sym": f[8]}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--fr3d", required=True); ap.add_argument("--python", required=True)
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--tmp", default=None); ap.add_argument("accs", nargs="*")
    a = ap.parse_args()
    accs = a.accs or ["RF00167", "RF00162", "RF00059", "RF00005", "RF00011", "RF00177"]
    commit = subprocess.run(["git", "-C", a.fr3d, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (OUT / "local_pairs").mkdir(parents=True, exist_ok=True)
    geo = {(r["acc"], r["seq_name"], int(r["pos_1"]), int(r["pos_2"])): r["geometry"] for r in csv.DictReader(open(HERE / "followup/unpaired_geometry.csv")) if r["kind"] == "unpaired"}
    conc_rows, near_rows, summary, failures = [], [], {"fr3d_python_commit": commit, "coordinates": "PDBe updated mmCIF, model 1", "families": {}}, []
    for acc in accs:
        fam, it, chains = load(acc)
        by_entry = defaultdict(dict)
        for sn, ch in chains.items():
            by_entry[ch["meta"]["pdb_id"].lower()][ch["meta"]["auth_asym_id"]] = sn
        local = {}  # seq_name -> {(p1,p2): set(classes)} on consensus positions, near and true
        def job(pdb):
            try:
                return pdb, run_entry(pdb, by_entry[pdb].keys(), a.fr3d, a.python, a.tmp), None
            except Exception as e:  # noqa: BLE001
                return pdb, None, str(e)[:300]
        with ThreadPoolExecutor(a.workers) as ex:
            for pdb, rows, err in ex.map(job, sorted(by_entry)):
                if err:
                    failures.append({"acc": acc, "pdb": pdb, "error": err}); print(acc, pdb, "FAILED", err, file=sys.stderr); continue
                for auth, sn in by_entry[pdb].items():
                    ch = chains[sn]
                    res = {(r["author_residue_number"], r["author_insertion_code"] or ""): r for r in ch["by_label"].values()}
                    d = defaultdict(set); out = []
                    for row in rows:
                        u1, cls, u2 = unit(row[0]), row[1], unit(row[2])
                        if u1["chain"] != auth or u2["chain"] != auth or u1["model"] != "1" or u2["model"] != "1" or (u1["sym"] not in ("", "1_555")) or (u2["sym"] not in ("", "1_555")):
                            continue
                        r1, r2 = res.get((u1["num"], u1["ins"])), res.get((u2["num"], u2["ins"]))
                        if not r1 or not r2: continue
                        l1, l2 = r1["polymer_seq_position"], r2["polymer_seq_position"]
                        if l1 > l2: l1, l2, r1, r2, cls = l2, l1, r2, r1, flip(cls)
                        p1 = r1["rfam_consensus_pos"] if r1["alignment_state"] == "match" else None
                        p2 = r2["rfam_consensus_pos"] if r2["alignment_state"] == "match" else None
                        out.append((l1, l2, p1, p2, cls))
                        if isinstance(p1, int) and isinstance(p2, int): d[(p1, p2)].add(cls)
                    out = sorted(set(out), key=lambda x: (x[0], x[1], x[4]))
                    with open(OUT / "local_pairs" / f"{sn}.tsv", "w") as f:
                        f.write("label_seq_id_1\tlabel_seq_id_2\trfam_1\trfam_2\tclass\n")
                        for o in out: f.write("\t".join("" if v is None else str(v) for v in o) + "\n")
                    local[sn] = d
                print(acc, pdb, file=sys.stderr)
        # concordance of true pairs on consensus positions
        for sn, d in local.items():
            ch = chains[sn]
            loc = {(k, c) for k, cs in d.items() for c in cs if not c.startswith("n")}
            bgsu = {(k, c) for k, cs in ch["bp"].items() for c in cs}
            conc_rows.append({"acc": acc, "seq_name": sn, "n_local_true": len(loc), "n_bgsu_true": len(bgsu), "shared": len(loc & bgsu),
                              "local_only": " ".join(f"{k[0]}-{k[1]}:{c}" for k, c in sorted(loc - bgsu)), "bgsu_only": " ".join(f"{k[0]}-{k[1]}:{c}" for k, c in sorted(bgsu - loc)),
                              "same_pairs_ignoring_class": len({k for k, _ in loc} & {k for k, _ in bgsu})})
        # SS_cons pair observations without a true BGSU base pair
        st = Counter(); bygeo = defaultdict(Counter)
        for i, j in ss_pairs(fam["ss_cons"]):
            for sn, ch in chains.items():
                if sn not in local or not evaluable(ch, i, j) or ch["bp"].get((i, j)): continue
                d = local[sn]; nears = sorted(c for c in d.get((i, j), set()) if c.startswith("n")); trues = sorted(c for c in d.get((i, j), set()) if not c.startswith("n"))
                elsewhere = lambda p: sorted(f"{q[0] if q[1] == p else q[1]}:{c}" for q, cs in d.items() if p in q and q != (i, j) for c in cs if c.startswith("n"))
                status = "local_true_pair" if trues else "near_pair" if nears else "no_pair"
                g = geo.get((acc, sn, i, j), "not_measured")
                st[status] += 1; bygeo[g][status] += 1
                near_rows.append({"acc": acc, "pos_1": i, "pos_2": j, "seq_name": sn, "nts": (ch["inst"]["nucleotides"][i - 1] + ch["inst"]["nucleotides"][j - 1]).upper(),
                                  "status": status, "near_classes": ",".join(nears), "local_true_classes": ",".join(trues), "geometry": g,
                                  "near_partners_i_elsewhere": " ".join(elsewhere(i)), "near_partners_j_elsewhere": " ".join(elsewhere(j))})
        cr = [r for r in conc_rows if r["acc"] == acc]
        summary["families"][acc] = {"chains_annotated": len(local), "chains_total": len(chains),
                                    "true_pairs_local": sum(r["n_local_true"] for r in cr), "true_pairs_bgsu": sum(r["n_bgsu_true"] for r in cr), "true_pairs_shared": sum(r["shared"] for r in cr),
                                    "chains_identical_true_pairs": sum(1 for r in cr if not r["local_only"] and not r["bgsu_only"]),
                                    "ss_obs_without_bgsu_pair": dict(st), "by_geometry": {g: dict(c) for g, c in bygeo.items()}}
        print(json.dumps({acc: summary["families"][acc]}), file=sys.stderr)
    summary["failures"] = failures
    for name, rows in (("concordance.csv", conc_rows), ("ss_cons_near.csv", near_rows)):
        if rows:
            with open(OUT / name, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    json.dump(summary, open(OUT / "near_summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
