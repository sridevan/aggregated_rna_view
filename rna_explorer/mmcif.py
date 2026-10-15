"""Minimal mmCIF reading (validated in poc_rf00177/scripts/common.py): category reader, polymer scheme, observed residues, entity source."""
from __future__ import annotations
import shlex
from pathlib import Path
from typing import Iterator, Optional


def _tokens(line: str) -> list[str]:
    try:
        return shlex.split(line, posix=True)
    except ValueError:
        return line.split()


def _multiline(first: str, it: Iterator[str]) -> str:
    val = first[1:].rstrip("\n")
    for l2 in it:
        if l2.startswith(";"):
            break
        val += "\n" + l2.rstrip("\n")
    return val


def cif_category(path: Path | str, category: str) -> list[dict[str, str]]:
    """Rows of an mmCIF category written as loop_ or key-value pairs; quoted and ;-delimited values handled; streams the file."""
    prefix = category + "."
    rows: list[dict[str, str]] = []
    keys: list[str] = []
    pending: list[str] = []
    cur: dict[str, str] = {}
    in_loop = target_loop = collecting = False
    with open(path) as f:
        it = iter(f)
        for line in it:
            if line.startswith("loop_"):
                in_loop, target_loop, collecting, keys, pending = True, False, True, [], []
                continue
            if line.startswith("#"):
                if (target_loop and rows) or cur:
                    break
                in_loop = target_loop = collecting = False
                continue
            if line.startswith("_"):
                name = line.split(None, 1)[0]
                if in_loop and collecting:
                    if name.startswith(prefix):
                        keys.append(name[len(prefix):])
                        target_loop = True
                    continue
                if not in_loop and name.startswith(prefix):
                    key = name[len(prefix):]
                    rest = line[len(name):].strip()
                    if rest == "":
                        nxt = next(it)
                        val = _multiline(nxt, it) if nxt.startswith(";") else nxt.strip()
                    else:
                        t = _tokens(rest)
                        val = t[0] if t else rest
                    cur[key] = val
                continue
            if in_loop and target_loop:
                collecting = False
                if line.startswith(";"):
                    pending.append(_multiline(line, it))
                else:
                    pending.extend(_tokens(line.strip()))
                while keys and len(pending) >= len(keys):
                    rows.append(dict(zip(keys, pending[:len(keys)])))
                    pending = pending[len(keys):]
            elif in_loop and collecting:
                collecting = False
    if cur:
        rows.append(cur)
    return rows


def poly_seq_scheme(path: Path | str, label_asym_id: str, entity_id: int | str) -> dict[int, dict[str, str]]:
    """_pdbx_poly_seq_scheme rows for one polymer asym, keyed by seq_id (label_seq_id)."""
    return {int(r["seq_id"]): r for r in cif_category(path, "_pdbx_poly_seq_scheme")
            if r["asym_id"] == label_asym_id and r["entity_id"] == str(entity_id)}


def atom_site_observed(path: Path | str, label_asym_id: str, entity_id: int | str) -> set[int]:
    """label_seq_id values with at least one atom for the polymer asym (streams atom_site)."""
    keys: list[str] = []
    seen: set[int] = set()
    idx: Optional[tuple[int, int, int]] = None
    with open(path) as f:
        for line in f:
            if line.startswith("_atom_site."):
                keys.append(line.split()[0].split(".", 1)[1])
                continue
            if keys and (line.startswith("ATOM") or line.startswith("HETATM")):
                if idx is None:
                    idx = (keys.index("label_asym_id"), keys.index("label_seq_id"), keys.index("label_entity_id"))
                t = line.split()
                if t[idx[0]] == label_asym_id and t[idx[2]] == str(entity_id) and t[idx[1]] not in (".", "?"):
                    seen.add(int(t[idx[1]]))
            elif keys and idx is not None and line.startswith("#"):
                break
    return seen


def entity_source(path: Path | str, entity_id: int | str) -> dict[str, Optional[str]]:
    for cat, org, tax in (("_entity_src_nat", "pdbx_organism_scientific", "pdbx_ncbi_taxonomy_id"),
                          ("_entity_src_gen", "pdbx_gene_src_scientific_name", "pdbx_gene_src_ncbi_taxonomy_id"),
                          ("_pdbx_entity_src_syn", "organism_scientific", "ncbi_taxonomy_id")):
        for r in cif_category(path, cat):
            if r.get("entity_id") == str(entity_id):
                o, t = r.get(org), r.get(tax)
                return {"cif_src_category": cat, "cif_organism": None if o in (None, "?", ".") else o, "cif_tax_id": None if t in (None, "?", ".") else t}
    return {"cif_src_category": None, "cif_organism": None, "cif_tax_id": None}


def atom_site_coords(path: Path | str, label_asym_id: str, entity_id: int | str, atoms: tuple[str, ...] = ("C4'", "P", "C1'")) -> dict[int, dict[str, tuple[float, float, float]]]:
    """Coordinates of selected atoms per label_seq_id for one polymer asym: first model (pdbx_PDB_model_num) and first alternate location only."""
    keys: list[str] = []
    out: dict[int, dict[str, tuple[float, float, float]]] = {}
    idx: Optional[dict[str, int]] = None
    first_model: Optional[str] = None
    wanted = set(atoms)
    with open(path) as f:
        for line in f:
            if line.startswith("_atom_site."):
                keys.append(line.split()[0].split(".", 1)[1]); continue
            if keys and (line.startswith("ATOM") or line.startswith("HETATM")):
                if idx is None:
                    idx = {k: keys.index(k) for k in ("label_asym_id", "label_seq_id", "label_entity_id", "label_atom_id", "label_alt_id", "pdbx_PDB_model_num", "Cartn_x", "Cartn_y", "Cartn_z") if k in keys}
                t = line.split()
                if t[idx["label_asym_id"]] != label_asym_id or t[idx["label_entity_id"]] != str(entity_id):
                    continue
                model = t[idx["pdbx_PDB_model_num"]] if "pdbx_PDB_model_num" in idx else "1"
                if first_model is None:
                    first_model = model
                if model != first_model:
                    continue
                name = t[idx["label_atom_id"]].strip('"')
                if name not in wanted or t[idx["label_seq_id"]] in (".", "?"):
                    continue
                seq = int(t[idx["label_seq_id"]])
                res = out.setdefault(seq, {})
                if name in res:   # keep the first alternate location
                    continue
                res[name] = (float(t[idx["Cartn_x"]]), float(t[idx["Cartn_y"]]), float(t[idx["Cartn_z"]]))
            elif keys and idx is not None and line.startswith("#"):
                break
    return out
