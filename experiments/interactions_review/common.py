"""Shared helpers for the scientific review of the Interactions datasets (read-only use of data/families/<ACC>/).

Status of an interaction in a chain follows the production definitions (rna_explorer/interactions.py, frontend lib/interactions.ts):
present = non-near FR3D annotation of that category between the two consensus positions; absent = both endpoints eligible for the
category but not annotated; not comparable otherwise (chain not annotated, category file missing, unmapped / unmodelled / incomplete).
"""
from __future__ import annotations
import collections, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class Family:
    def __init__(self, acc: str):
        self.acc = acc
        d = ROOT / "data" / "families" / acc
        self.dir = d
        self.ds = json.load(open(d / "interactions.json"))
        self.family = json.load(open(d / "family.json"))
        self.inst = {i["seq_name"]: i for i in json.load(open(d / "instances.json"))}
        self.ix = {i["seq_name"]: i for i in self.ds["instances"]}
        self.pairs = self.ds["pairs"]
        self.by_id = {p["id"]: p for p in self.pairs}
        sc = d / "structure_clustering.json"
        self.groups = json.load(open(sc))["metrics"]["discrepancy"]["clustering"]["assignments"] if sc.exists() else {}
        # present: seq_name -> {pair_id: [classes]}
        self.present: dict[str, dict[str, list[str]]] = {}
        for sn, i in self.ix.items():
            m: dict[str, list[str]] = collections.defaultdict(list)
            for cls, ks in i["observed"].items():
                for k in ks:
                    m[self.pairs[k]["id"]].append(cls)
            self.present[sn] = dict(m)
        self.ss = self._ss_pairs()
        self.helix = self._helix_labels()

    # ------------------------------------------------------------------ structure labels
    def _ss_pairs(self) -> dict[int, int]:
        st, out = collections.defaultdict(list), {}
        for i, ch in enumerate(self.family["ss_cons"], 1):
            for o, c in ("()", "<>", "[]", "{}"):
                if ch == o:
                    st[o].append(i)
                elif ch == c:
                    j = st[o].pop(); out[i], out[j] = j, i
        return out

    def _helix_labels(self) -> dict[int, str]:
        """stem index (S1, S2, ... in 5' order of their opening strand) for every paired consensus position; contiguous stacked pairs
        (i, j), (i+1, j-1) with at most a 1-nt bulge on either side form one stem"""
        lab, n = {}, 0
        opens = sorted(i for i, j in self.ss.items() if i < j)
        prev = None
        for i in opens:
            j = self.ss[i]
            if prev and i - prev[0] <= 2 and prev[1] - j <= 2 and j < prev[1]:
                pass
            else:
                n += 1
            lab[i] = lab[j] = f"S{n}"
            prev = (i, j)
        return lab

    # ------------------------------------------------------------------ status
    def status(self, sn: str, p: dict) -> str:
        i = self.ix[sn]
        if i["annotation"] != "available" or p["category"] in i.get("unavailable_categories", []):
            return "nc"
        if any(a <= q <= b for a, b in i.get("annotation_gaps", []) for q in (p["pos_1"], p["pos_2"])):
            return "nc"
        second = {"baseribose": i["sugar_eligibility"], "basephosphate": i["phosphate_eligibility"]}.get(p["category"], i["eligibility"])
        if i["eligibility"][p["pos_1"] - 1] != "E" or second[p["pos_2"] - 1] != "E":
            return "nc"
        return "present" if p["id"] in self.present[sn] else "absent"

    def freq(self, p: dict, chains) -> tuple[int, int]:
        k = n = 0
        for sn in chains:
            s = self.status(sn, p)
            if s != "nc":
                n += 1; k += s == "present"
        return k, n

    # ------------------------------------------------------------------ numbering
    def numbering(self, sn: str) -> dict[int, str]:
        m = json.load(open(self.dir / self.inst[sn]["mapping_file"]))
        return {r["rfam_consensus_pos"]: f"{r['residue_name']}{r['author_residue_number']}{r['author_insertion_code'] or ''}"
                for r in m["residues"] if r["alignment_state"] == "match"}

    def rfam_of(self, sn: str) -> dict[str, int]:
        """author number (+ins) -> Rfam position, for one chain"""
        m = json.load(open(self.dir / self.inst[sn]["mapping_file"]))
        return {f"{r['author_residue_number']}{r['author_insertion_code'] or ''}": r["rfam_consensus_pos"] for r in m["residues"] if r["alignment_state"] == "match"}

    def coverage(self, sn: str) -> tuple[int, int]:
        """consensus helix base pairs annotated / evaluable (as helixCoverage in the frontend)"""
        k = n = 0
        for p in self.pairs:
            if p["category"] == "basepair" and p["in_ss_cons"]:
                s = self.status(sn, p)
                if s != "nc":
                    n += 1; k += s == "present"
        return k, n
