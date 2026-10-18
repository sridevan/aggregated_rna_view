"""Annotation backend (Milestone 3B): reference structural elements, Motif Atlas loops, PDBe bound molecules, per-instance
annotation records and consensus-position tracks for one built family.

Counting rules (see docs in annotation_assessment/PROPOSED_ANNOTATION_SCHEMA.md):
  * reference_structural_elements: one categorical family-level track transferred from the BGSU reference chain through the
    validated Rfam mapping. Instance labels (BGSU per-chain transfer) are retained for highlighting, never counted as
    independent structural observations.
  * loops: one FR3D loop = one observation, whatever the number of chains it spans. Counting units at a consensus position are
    distinct loop ids, distinct instances (chains) and distinct PDB entries. J3 and J4 stay separate; the junction union is a
    set union, never a sum.
  * bound molecules: PDBe binding_sites are entity-level. Counting unit = (pdb_id, entity_id); chain attribution is only
    'resolved' when the entity has a single chain among the family instances.
  * denominators: for every position and source the numbers of observable instances (match residue with coordinates),
    source-evaluated instances, positives, 'not annotated by the source' and 'source unavailable' are stored explicitly.
"""
from __future__ import annotations
import datetime as _dt
import json
import re
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional
from . import __version__
from .annotation_mapping import InstanceIndex, ResidueRef, UnitId
from .annotation_sources import (LOOP_TYPES, AnnotationFetcher, fetch_binding_sites, fetch_chain_annotations, fetch_entry_loops,
                                 fetch_entry_motifs, fetch_family_reference, fetch_motif_atlas_release)
from .export import validate_json

log = logging.getLogger("rna_explorer.annotations")
SCHEMA_VERSION = "1.0"
CONFIG_PATH = Path(__file__).parent / "annotation_config.json"
BGSU_GROUP_URL = "https://rna.bgsu.edu/rna3dhub/motif/view/{group}"
BGSU_LOOP_URL = "https://rna.bgsu.edu/rna3dhub/loops/view/{loop_id}"


def load_config(path: Path = CONFIG_PATH) -> dict[str, Any]:
    return json.load(open(path))


def _frac(num: int, den: int) -> Optional[float]:
    return round(num / den, 4) if den else None


def _dump(path: Path, obj: Any) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)
    return path.stat().st_size


# ============================================================================ builder
class AnnotationBuilder:
    def __init__(self, fam_dir: Path, fetcher: AnnotationFetcher, config: Optional[dict[str, Any]] = None) -> None:
        self.fam_dir = Path(fam_dir)
        self.fetcher = fetcher
        self.config = config or load_config()
        self.idx = InstanceIndex(self.fam_dir)
        self.acc: str = self.idx.family["rfam_acc"]
        self.clen = self.idx.clen
        self.generated = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
        # filled by build steps
        self.reference: dict[str, Any] = {}
        self.inst_elements: dict[str, dict[int, dict[str, Any]]] = {}      # seq_name -> label_seq_id -> {labels, unit_id, status}
        self.inst_element_source: dict[str, dict[str, Any]] = {}           # seq_name -> source meta/status
        self.loops: list[dict[str, Any]] = []
        self.loops_outside: list[dict[str, Any]] = []     # loops of the same entries with no nucleotide in a family chain (compact records)
        self.loop_source: dict[str, dict[str, Any]] = {}                   # pdb_id -> source status
        self.atlas_release: dict[str, Any] = {}
        self.ligands: list[dict[str, Any]] = []
        self.ligand_source: dict[str, dict[str, Any]] = {}                 # "pdb:entity" -> status
        self.unresolved: list[dict[str, Any]] = []

    # ------------------------------------------------------------------ provenance
    def _provenance(self, **extra: Any) -> dict[str, Any]:
        mp = self.idx.mapping(self.idx.instances[0]["seq_name"])["provenance"]
        return {"rfam_acc": self.acc, "rfam_release": mp.get("rfam_release"), "cm_checksum": mp.get("cm_checksum"), "cm_clen": self.clen,
                "mapping_pipeline_version": mp.get("pipeline_version"), "mapping_generated": mp.get("generated"),
                "annotation_pipeline_version": __version__, "annotation_config_version": self.config.get("config_version"), "generated": self.generated, **extra}

    def _note_unresolved(self, source: str, ref: ResidueRef, context: str) -> None:
        if ref.mapping_status == "unresolved_identifier":
            self.unresolved.append({"source": source, "identifier": ref.source_identifier, "pdb_id": ref.pdb_id, "auth_asym_id": ref.auth_asym_id,
                                    "reason": "chain not among family instances" if ref.seq_name is None else "residue not found in chain mapping", "context": context})

    # ------------------------------------------------------------------ 1. reference structural elements
    def build_reference(self) -> dict[str, Any]:
        rows, comments, meta = fetch_family_reference(self.fetcher, self.acc)
        track: list[Optional[str]] = [None] * self.clen
        positions: list[dict[str, Any]] = []
        not_on_consensus: list[dict[str, Any]] = []
        conflicts: list[dict[str, Any]] = []
        ref_chains = sorted({"|".join(r.unit_id.split("|")[:3]) for r in rows})
        ref_pdb = ref_auth = ref_inst = None
        if ref_chains:
            u = UnitId.parse(rows[0].unit_id)
            ref_pdb, ref_auth = u.pdb_id, u.auth_asym_id
            ref_inst = self.idx.by_auth.get((ref_pdb, ref_auth))
        for r in rows:
            ref = self.idx.resolve_unit(r.unit_id)
            self._note_unresolved("bgsu_reference_annotation", ref, "family reference file")
            if ref.rfam_consensus_pos is not None:
                pos = ref.rfam_consensus_pos
                if track[pos - 1] is not None and track[pos - 1] != r.labels[0]:
                    conflicts.append({"rfam_consensus_pos": pos, "labels": [track[pos - 1], r.labels[0]], "unit_id": r.unit_id})
                    continue
                track[pos - 1] = r.labels[0]
                positions.append({"rfam_consensus_pos": pos, "element": r.labels[0], "all_labels": r.labels, "reference_unit_id": r.unit_id,
                                  "reference_label_seq_id": ref.label_seq_id, "reference_author_residue_number": ref.author_residue_number,
                                  "reference_nucleotide": ref.residue_name, "reference_pp_char": ref.pp_char, "reference_pp_category": ref.pp_category,
                                  "reference_has_coordinates": ref.has_coordinates})
            else:
                not_on_consensus.append({"unit_id": r.unit_id, "labels": r.labels, "mapping_status": ref.mapping_status, "label_seq_id": ref.label_seq_id})
        positions.sort(key=lambda p: p["rfam_consensus_pos"])
        # element catalogue with consensus spans (contiguous runs)
        elements: list[dict[str, Any]] = []
        for p in positions:
            if elements and elements[-1]["element"] == p["element"] and elements[-1]["spans"][-1][1] == p["rfam_consensus_pos"] - 1:
                elements[-1]["spans"][-1][1] = p["rfam_consensus_pos"]
            elif elements and any(e["element"] == p["element"] for e in elements):
                e = next(e for e in elements if e["element"] == p["element"])
                if e["spans"][-1][1] == p["rfam_consensus_pos"] - 1:
                    e["spans"][-1][1] = p["rfam_consensus_pos"]
                else:
                    e["spans"].append([p["rfam_consensus_pos"], p["rfam_consensus_pos"]])
            else:
                elements.append({"element": p["element"], "spans": [[p["rfam_consensus_pos"], p["rfam_consensus_pos"]]]})
        for e in elements:
            e["n_positions"] = sum(b - a + 1 for a, b in e["spans"])
        # per-instance transfer (BGSU per-chain endpoint), retained for highlighting only
        evaluated, unavailable, status_counts, agree, disagree = [], [], Counter(), 0, []
        n_labels = 0
        for inst in self.idx.instances:
            sn = inst["seq_name"]
            crow, cmeta = fetch_chain_annotations(self.fetcher, inst["pdb_id"], inst["auth_asym_id"])
            if cmeta.status != "ok" or not crow:
                unavailable.append(sn)
                self.inst_element_source[sn] = {"status": "unavailable", "detail": cmeta.status, "url": cmeta.url, "retrieved": cmeta.retrieved}
                self.inst_elements[sn] = {}
                continue
            evaluated.append(sn)
            self.inst_element_source[sn] = {"status": "evaluated", "url": cmeta.url, "retrieved": cmeta.retrieved, "n_labelled_nucleotides": len(crow)}
            d: dict[int, dict[str, Any]] = {}
            for r in crow:
                ref = self.idx.resolve_unit(r.unit_id)
                self._note_unresolved("bgsu_chain_annotation", ref, sn)
                status_counts[ref.mapping_status] += 1
                n_labels += 1
                if ref.label_seq_id is not None:
                    d[ref.label_seq_id] = {"labels": r.labels, "element": r.labels[0], "unit_id": r.unit_id, "mapping_status": ref.mapping_status,
                                           "rfam_consensus_pos": ref.rfam_consensus_pos}
                if ref.rfam_consensus_pos is not None:
                    if track[ref.rfam_consensus_pos - 1] == r.labels[0]:
                        agree += 1
                    else:
                        disagree.append({"seq_name": sn, "unit_id": r.unit_id, "instance_label": r.labels[0], "family_track_label": track[ref.rfam_consensus_pos - 1],
                                         "rfam_consensus_pos": ref.rfam_consensus_pos})
            self.inst_elements[sn] = d
        self.reference = {
            "schema_version": SCHEMA_VERSION, "rfam_acc": self.acc, "clen": self.clen, "track_name": "reference_structural_elements",
            "track_kind": "categorical_reference_transfer",
            "interpretation": "Element labels of one curated reference chain transferred to Rfam consensus positions through the validated "
                              "residue mapping. They are not independent observations of helices or loops in each RNA instance.",
            "reference": {"pdb_id": ref_pdb, "auth_asym_id": ref_auth, "label_asym_id": ref_inst["label_asym_id"] if ref_inst else None,
                          "entity_id": ref_inst["entity_id"] if ref_inst else None, "seq_name": ref_inst["seq_name"] if ref_inst else None,
                          "model": 1, "reference_chains_in_file": ref_chains, "n_annotated_nucleotides": len(rows), "citation": comments},
            "source": {"name": "BGSU nucleotide annotations (family file)", "url": meta.url, "retrieved": meta.retrieved, "status": meta.status,
                       "annotation_identifier": f"{self.acc}_annotations.txt"},
            "transfer_method": "reference unit id -> (pdb_id, label_asym_id, label_seq_id) by author number + insertion code -> rfam_consensus_pos of the "
                               "reference chain mapping (rna_explorer mapping pipeline)",
            "elements": elements, "positions": positions, "track": track,
            "n_consensus_positions_annotated": len(positions), "reference_nucleotides_not_on_consensus": not_on_consensus, "conflicts": conflicts,
            "instance_transfer": {"source": "BGSU correspondence/nucleotide_annotation per chain (labels transferred by BGSU)", "n_instances": len(self.idx.instances),
                                  "n_instances_evaluated": len(evaluated), "n_instances_unavailable": len(unavailable), "unavailable_instances": unavailable,
                                  "n_nucleotide_labels": n_labels, "mapping_status_counts": dict(status_counts), "n_labels_agree_with_family_track": agree,
                                  "n_labels_disagree": len(disagree), "disagreements": disagree[:200]},
            "provenance": self._provenance()}
        return self.reference

    # ------------------------------------------------------------------ 2. Motif Atlas loops
    @staticmethod
    def _strands(nts: list[dict[str, Any]]) -> list[dict[str, Any]]:
        strands: list[dict[str, Any]] = []
        prev_chain, prev_num = None, None
        for i, n in enumerate(nts):
            u = UnitId.parse(n["unit_id"])
            num = int(u.author_residue_number) if u.author_residue_number.lstrip("-").isdigit() else None
            new = prev_chain != u.auth_asym_id or num is None or prev_num is None or (num - prev_num) not in (0, 1)
            if new:
                strands.append({"auth_asym_id": u.auth_asym_id, "start_unit_id": n["unit_id"], "end_unit_id": n["unit_id"], "nucleotide_indices": [i, i], "n_nucleotides": 1})
            else:
                strands[-1]["end_unit_id"] = n["unit_id"]; strands[-1]["nucleotide_indices"][1] = i; strands[-1]["n_nucleotides"] += 1
            prev_chain, prev_num = u.auth_asym_id, num
        return strands

    def build_loops(self) -> dict[str, Any]:
        for t in LOOP_TYPES:
            rel, rmeta = fetch_motif_atlas_release(self.fetcher, t)
            self.atlas_release[t] = {"release": rel, "url": rmeta.url, "retrieved": rmeta.retrieved, "status": rmeta.status}
        motif_page_only: list[dict[str, Any]] = []
        for pdb in sorted(self.idx.by_pdb):
            loops, lmeta = fetch_entry_loops(self.fetcher, pdb)
            motifs, mmeta = fetch_entry_motifs(self.fetcher, pdb)
            evaluated = lmeta.status in ("ok", "empty") and mmeta.status in ("ok",)
            self.loop_source[pdb] = {"status": "evaluated" if evaluated else "unavailable", "loops": lmeta.to_dict(), "motifs": mmeta.to_dict(), "n_loops": len(loops)}
            by_id = {m.loop_id: m for m in motifs}
            seen = set()
            n_outside = 0
            for L in loops:
                seen.add(L.loop_id)
                m = by_id.get(L.loop_id)
                refs = [(u, self.idx.resolve_unit(u)) for u in L.units]
                # a loop with no nucleotide in any family chain (e.g. 23S, 5S or tRNA loops of a ribosome entry) is kept as a compact record,
                # without nucleotide-level data and without counting its nucleotides as unresolved; loops touching a family chain are unchanged
                if not any(ref.seq_name for _, ref in refs):
                    n_outside += 1
                    self.loops_outside.append({"loop_id": L.loop_id, "pdb_id": pdb, "loop_type": L.loop_type, "auth_asym_ids": sorted({UnitId.parse(u).auth_asym_id for u in L.units}),
                                               "n_nucleotides": len(L.units), "motif_group": m.motif_group if m else None, "atlas_status": m.atlas_status if m else "not_listed",
                                               "loop_url": BGSU_LOOP_URL.format(loop_id=L.loop_id), "reason": "no nucleotide in a family chain"})
                    continue
                nts: list[dict[str, Any]] = []
                for u, ref in refs:
                    self._note_unresolved("bgsu_motif_atlas_loop", ref, L.loop_id)
                    d = ref.to_dict(); d["unit_id"] = u
                    d["resolution"] = ("chain_not_in_family" if ref.seq_name is None else "residue_not_found" if ref.mapping_status == "unresolved_identifier" else "resolved")
                    nts.append(d)
                chains: list[dict[str, Any]] = []
                for u in L.units:
                    uu = UnitId.parse(u)
                    if not any(c["auth_asym_id"] == uu.auth_asym_id for c in chains):
                        inst = self.idx.by_auth.get((uu.pdb_id, uu.auth_asym_id))
                        chains.append({"auth_asym_id": uu.auth_asym_id, "label_asym_id": inst["label_asym_id"] if inst else None,
                                       "seq_name": inst["seq_name"] if inst else None, "entity_id": inst["entity_id"] if inst else None, "in_family": inst is not None})
                cons = sorted({n["rfam_consensus_pos"] for n in nts if n["mapping_status"] in ("mapped", "mapped_no_coordinates")})
                rec = {"loop_id": L.loop_id, "pdb_id": pdb, "loop_type": L.loop_type, "loop_url": BGSU_LOOP_URL.format(loop_id=L.loop_id),
                       "atlas_status": m.atlas_status if m else "not_listed", "motif_group": m.motif_group if m else None,
                       "motif_group_url": BGSU_GROUP_URL.format(group=m.motif_group) if m and m.motif_group else None,
                       "motif_group_assignment": m.assignment_type if m else None, "motif_group_via_loop_id": m.mapped_via_loop_id if m else None,
                       "annotation_text": m.annotation_text if m else None, "problem": m.problem if m else None, "location": m.location if m else None,
                       "n_nucleotides": len(L.units), "chains": chains, "multi_chain": len(chains) > 1,
                       "instances": [c["seq_name"] for c in chains if c["seq_name"]], "strands": self._strands(nts), "nucleotides": nts,
                       "consensus_positions": cons, "mapping_status_counts": dict(Counter(n["mapping_status"] for n in nts)),
                       "pp_category_counts": dict(Counter(n["pp_category"] for n in nts if n["mapping_status"] == "mapped")),
                       "source": {"loops_url": lmeta.url, "motifs_url": mmeta.url, "retrieved": lmeta.retrieved,
                                  "motif_atlas_release": self.atlas_release.get(L.loop_type, {}).get("release")}}
                self.loops.append(rec)
            self.loop_source[pdb]["n_loops_outside_family_chains"] = n_outside
            for m in motifs:
                if m.loop_id not in seen:
                    motif_page_only.append({"pdb_id": pdb, "loop_id": m.loop_id, "atlas_status": m.atlas_status, "motif_group": m.motif_group})
        by_type: dict[str, Any] = {}
        for t in sorted({L["loop_type"] for L in self.loops} | set(LOOP_TYPES)):
            ls = [L for L in self.loops if L["loop_type"] == t]
            groups = Counter(L["motif_group"] for L in ls if L["motif_group"])
            by_type[t] = {"n_loops": len(ls), "n_with_motif_group": sum(1 for L in ls if L["motif_group"]), "n_without_motif_group": sum(1 for L in ls if not L["motif_group"]),
                          "motif_group_assignment_counts": dict(Counter(L["motif_group_assignment"] or "none" for L in ls)),
                          "n_problematic": sum(1 for L in ls if L["atlas_status"] == "problematic"), "n_distinct_motif_groups": len(groups),
                          "motif_groups": {g: {"n_loops": n, "annotation_text": next((L["annotation_text"] for L in ls if L["motif_group"] == g and L["annotation_text"]), None),
                                               "url": BGSU_GROUP_URL.format(group=g)} for g, n in sorted(groups.items())},
                          "n_multi_chain_loops": sum(1 for L in ls if L["multi_chain"]), "n_pdb_entries": len({L["pdb_id"] for L in ls}),
                          "n_instances": len({s for L in ls for s in L["instances"]})}
        return {"schema_version": SCHEMA_VERSION, "rfam_acc": self.acc, "clen": self.clen,
                "source": {"name": "BGSU RNA 3D Motif Atlas / FR3D loops", "loops_url_template": "https://rna.bgsu.edu/rna3dhub/loops/download/{PDB}",
                           "motifs_url_template": "https://rna.bgsu.edu/rna3dhub/pdb/{PDB}/motifs", "motif_atlas_release": self.atlas_release,
                           "note": "Loop -> motif-group assignments are parsed from the HTML motifs page (no JSON endpoint). 'problematic' loops are listed by the "
                                   "Atlas with a reason and carry no motif group. Assignment types: 'member' = the loop is an instance of the group in the current "
                                   "release; 'homologous'/'geometric' = the loop is mapped to the group through a release loop (motif_group_via_loop_id)."},
                "counting_rules": {"loop_observation": "one FR3D loop id = one observation, even when its nucleotides span several chains",
                                   "junction": "J3 and J4 are separate categories; any combined junction summary is a set union of loops/instances/entries"},
                "counts": {"n_loops": len(self.loops), "n_with_motif_group": sum(1 for L in self.loops if L["motif_group"]),
                           "motif_group_assignment_counts": dict(Counter(L["motif_group_assignment"] or "none" for L in self.loops)),
                           "n_problematic": sum(1 for L in self.loops if L["atlas_status"] == "problematic"),
                           "n_distinct_motif_groups": len({L["motif_group"] for L in self.loops if L["motif_group"]}), "by_type": by_type,
                           "n_pdb_entries_evaluated": sum(1 for v in self.loop_source.values() if v["status"] == "evaluated"),
                           "n_pdb_entries_unavailable": sum(1 for v in self.loop_source.values() if v["status"] != "evaluated"),
                           "nucleotide_mapping_status_counts": dict(Counter(n["mapping_status"] for L in self.loops for n in L["nucleotides"])),
                           "n_motif_page_loops_missing_from_download": len(motif_page_only), "n_loops_outside_family_chains": len(self.loops_outside)},
                "entry_sources": self.loop_source, "motif_page_only": motif_page_only, "loops": self.loops,
                "loops_outside_family_chains": {"description": "loops of the same PDB entries with no nucleotide in a family chain (other molecules of the entry); "
                                                               "kept as compact records without nucleotide-level data", "loops": self.loops_outside},
                "provenance": self._provenance()}

    # ------------------------------------------------------------------ 3. PDBe bound molecules
    def _classify(self, ccd: str) -> tuple[str, str]:
        for cls, ids in self.config["bound_molecule_classes"].items():
            if cls.startswith("_"):
                continue
            if ccd.upper() in {x.upper() for x in ids}:
                return cls, f"listed:{cls}"
        return "ligand", "default:ligand"

    def build_ligands(self) -> dict[str, Any]:
        default_classes = set(self.config["default_track_classes"])
        rules = {r["id"]: r for r in self.config["artefact_rules"]}
        for (pdb, entity), chains in sorted(self.idx.by_entity.items()):
            ent, meta = fetch_binding_sites(self.fetcher, pdb, entity)
            key = f"{pdb}:{entity}"
            if meta.status in ("http_error", "parse_error"):
                self.ligand_source[key] = {**meta.to_dict(), "response_status": meta.status, "status": "unavailable", "n_bound_molecules": 0}
                continue
            self.ligand_source[key] = {**meta.to_dict(), "response_status": meta.status, "status": "evaluated", "n_bound_molecules": len(ent["data"]) if ent else 0}
            if not ent:
                continue
            seq_names = [c["seq_name"] for c in chains]
            if len({c["sequence_key"] for c in chains}) > 1:
                log.warning("entity %s has chains with different sequence keys", key)
            plen = int(ent.get("length") or self.idx.polymer_length(seq_names[0]))
            for site in ent["data"]:
                ccd = site["accession"].upper()
                cls, rule = self._classify(ccd)
                add = site.get("additionalData", {})
                best = add.get("bestChainId")
                best_inst = self.idx.by_label.get((pdb, best)) if best else None
                map_inst = best_inst if best_inst and best_inst["entity_id"] == entity else chains[0]
                contacts: dict[int, dict[str, Any]] = {}
                for seg in site.get("residues", []):
                    bm = seg.get("additionalData", {}).get("boundMoleculeId")
                    for i in range(int(seg["startIndex"]), int(seg["endIndex"]) + 1):
                        c = contacts.setdefault(i, {"label_seq_id": i, "bound_molecule_ids": set(), "index_type": seg.get("indexType")})
                        if bm:
                            c["bound_molecule_ids"].add(bm)
                        if i == seg["startIndex"]:
                            c["sequence_code"] = seg.get("startCode")
                        if i == seg["endIndex"]:
                            c["sequence_code"] = c.get("sequence_code") or seg.get("endCode")
                out_contacts: list[dict[str, Any]] = []
                for i in sorted(contacts):
                    ref = self.idx.resolve_seq_index(map_inst["seq_name"], i, f"{pdb}|{entity}|seq{i}")
                    self._note_unresolved("pdbe_binding_sites", ref, f"{ccd} {key}")
                    coords = [c["label_asym_id"] for c in chains if (r := self.idx.resolve_seq_index(c["seq_name"], i)).has_coordinates]
                    out_contacts.append({"label_seq_id": i, "sequence_code": contacts[i].get("sequence_code"), "index_type": contacts[i]["index_type"],
                                         "bound_molecule_ids": sorted(contacts[i]["bound_molecule_ids"]), "residue_name": ref.residue_name,
                                         "rfam_consensus_pos": ref.rfam_consensus_pos, "alignment_state": ref.alignment_state, "mapping_status": ref.mapping_status,
                                         "pp_char": ref.pp_char, "pp_category": ref.pp_category, "chains_with_coordinates": coords})
                flags: list[str] = []
                if cls in rules["ccd_is_polymer_residue"]["applies_to_classes"] and ccd in self.idx.polymer_residue_names(map_inst["seq_name"]):
                    flags.append("ccd_is_polymer_residue")
                w = int(rules["contacts_only_polymer_termini"].get("terminal_window", 2))
                if cls in rules["contacts_only_polymer_termini"]["applies_to_classes"] and out_contacts and all(c["label_seq_id"] <= w or c["label_seq_id"] > plen - w for c in out_contacts):
                    flags.append("contacts_only_polymer_termini")
                excluded_by_flag = self.config.get("exclude_flagged_from_default_track", True) and any(rules[f].get("exclude_from_default_track") for f in flags)
                rec = {"pdb_id": pdb, "entity_id": entity, "ccd_id": ccd, "name": site.get("name"), "ligand_class": cls, "class_rule": rule,
                       "scaffold_id": add.get("scaffoldId") or None, "cofactor_id": add.get("cofactorId") or None, "residue_count_reported": add.get("residueCount"),
                       "observation_level": "entity", "best_chain_label_asym_id": best, "best_chain_seq_name": best_inst["seq_name"] if best_inst else None,
                       "n_entity_chains": len(chains), "entity_chains": seq_names, "chain_attribution": "resolved" if len(chains) == 1 else "unresolved",
                       "attribution_note": None if len(chains) == 1 else "entity-level residue pooling: bestChainId is metadata, not proof of chain-specific contact",
                       "bound_molecule_ids": sorted({b for c in contacts.values() for b in c["bound_molecule_ids"]}), "n_bound_molecule_instances": len({b for c in contacts.values() for b in c["bound_molecule_ids"]}),
                       "contacts": out_contacts, "n_contact_positions": len(out_contacts),
                       "consensus_positions": sorted({c["rfam_consensus_pos"] for c in out_contacts if c["mapping_status"] in ("mapped", "mapped_no_coordinates")}),
                       "flags": flags, "in_default_track": cls in default_classes and not excluded_by_flag,
                       "exclusion_reason": None if cls in default_classes and not excluded_by_flag else (f"class:{cls}" if cls not in default_classes else "flag:" + ",".join(f for f in flags if rules[f].get("exclude_from_default_track"))),
                       "source": {"url": meta.url, "retrieved": meta.retrieved}}
                self.ligands.append(rec)
        by_ccd: dict[str, Any] = {}
        for r in self.ligands:
            d = by_ccd.setdefault(r["ccd_id"], {"name": r["name"], "ligand_class": r["ligand_class"], "entry_entities": set(), "entries": set(), "n_bound_molecule_instances": 0,
                                                 "consensus_positions": Counter(), "flags": Counter(), "in_default_track_entries": set()})
            d["entry_entities"].add(f"{r['pdb_id']}:{r['entity_id']}"); d["entries"].add(r["pdb_id"]); d["n_bound_molecule_instances"] += r["n_bound_molecule_instances"]
            for p in r["consensus_positions"]:
                d["consensus_positions"][p] += 1
            for f in r["flags"]:
                d["flags"][f] += 1
            if r["in_default_track"]:
                d["in_default_track_entries"].add(r["pdb_id"])
        summary = {ccd: {"name": d["name"], "ligand_class": d["ligand_class"], "n_entry_entities": len(d["entry_entities"]), "n_entries": len(d["entries"]), "entries": sorted(d["entries"]),
                         "n_bound_molecule_instances": d["n_bound_molecule_instances"], "n_entries_in_default_track": len(d["in_default_track_entries"]),
                         "consensus_positions_n_entries": {str(p): n for p, n in sorted(d["consensus_positions"].items())}, "flag_counts": dict(d["flags"])}
                   for ccd, d in sorted(by_ccd.items())}
        default_entries = {r["pdb_id"] for r in self.ligands if r["in_default_track"]}
        pocket = Counter(p for r in self.ligands if r["in_default_track"] for p in r["consensus_positions"])
        return {"schema_version": SCHEMA_VERSION, "rfam_acc": self.acc, "clen": self.clen,
                "source": {"name": "PDBe API v2 binding_sites", "url_template": "https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/binding_sites/{pdb_id}/{entity_id}",
                           "observation_level": "entity", "index_type": "PDB polymer sequence index == label_seq_id (validated Milestone 3A)"},
                "classification": {"config_version": self.config.get("config_version"), "default_track_classes": sorted(default_classes),
                                   "exclude_flagged_from_default_track": self.config.get("exclude_flagged_from_default_track", True),
                                   "artefact_rules": self.config["artefact_rules"], "classes": {k: v for k, v in self.config["bound_molecule_classes"].items() if not k.startswith("_")}},
                "counting_rules": {"primary_unit": "pdb_entry_entity (pdb_id, entity_id)", "contacts": "deduplicated per (entry, entity, ccd_id, label_seq_id); bound-molecule ids kept per contact",
                                   "chain_attribution": "resolved only for single-chain entities; otherwise unresolved, bestChainId retained as metadata"},
                "counts": {"n_entry_entities_evaluated": sum(1 for v in self.ligand_source.values() if v["status"] == "evaluated"),
                           "n_entry_entities_unavailable": sum(1 for v in self.ligand_source.values() if v["status"] != "evaluated"),
                           "n_bound_molecule_records": len(self.ligands), "n_distinct_ccd_ids": len(by_ccd),
                           "by_class": dict(Counter(d["ligand_class"] for d in by_ccd.values())),
                           "n_records_in_default_track": sum(1 for r in self.ligands if r["in_default_track"]), "n_entries_in_default_track": len(default_entries),
                           "n_records_flagged": sum(1 for r in self.ligands if r["flags"]), "flag_counts": dict(Counter(f for r in self.ligands for f in r["flags"])),
                           "chain_attribution": dict(Counter(r["chain_attribution"] for r in self.ligands)),
                           "contact_mapping_status_counts": dict(Counter(c["mapping_status"] for r in self.ligands for c in r["contacts"])),
                           "default_track_pocket_positions_n_entries": {str(p): n for p, n in sorted(pocket.items())}},
                "entry_entity_sources": self.ligand_source, "by_ccd": summary, "bound_molecules": self.ligands, "provenance": self._provenance()}

    # ------------------------------------------------------------------ 4. per-instance records
    def build_instance(self, inst: dict[str, Any]) -> dict[str, Any]:
        sn, pdb, entity = inst["seq_name"], inst["pdb_id"], int(inst["entity_id"])
        loops_here = [L for L in self.loops if sn in L["instances"]]
        loop_by_res: dict[int, list[dict[str, Any]]] = defaultdict(list)
        for L in loops_here:
            for n in L["nucleotides"]:
                if n["seq_name"] == sn and n["label_seq_id"] is not None:
                    loop_by_res[n["label_seq_id"]].append({"loop_id": L["loop_id"], "loop_type": L["loop_type"], "motif_group": L["motif_group"], "atlas_status": L["atlas_status"]})
        ligs = [r for r in self.ligands if r["pdb_id"] == pdb and r["entity_id"] == entity]
        lig_by_res: dict[int, list[dict[str, Any]]] = defaultdict(list)
        for r in ligs:
            attribution = "chain" if r["chain_attribution"] == "resolved" else "entity_unresolved"
            for c in r["contacts"]:
                lig_by_res[c["label_seq_id"]].append({"ccd_id": r["ccd_id"], "ligand_class": r["ligand_class"], "bound_molecule_ids": c["bound_molecule_ids"], "attribution": attribution,
                                                      "is_best_chain": r["best_chain_seq_name"] == sn, "in_default_track": r["in_default_track"], "flags": r["flags"]})
        elems = self.inst_elements.get(sn, {})
        residues = []
        for res in self.idx.residues(sn):
            i = int(res["polymer_seq_position"])
            e = elems.get(i)
            residues.append({"label_seq_id": i, "author_residue_number": res["author_residue_number"], "author_insertion_code": res["author_insertion_code"] or "",
                             "residue_name": res["residue_name"], "alignment_state": res["alignment_state"], "rfam_consensus_pos": res["rfam_consensus_pos"],
                             "insertion_anchor": res["insertion_anchor"], "pp_char": res["posterior_prob_char"] or "", "has_coordinates": bool(res["has_coordinates"]),
                             "reference_element": e["element"] if e else None, "reference_element_unit_id": e["unit_id"] if e else None,
                             "loops": loop_by_res.get(i, []), "bound_molecules": lig_by_res.get(i, [])})
        esrc = self.inst_element_source.get(sn, {"status": "unavailable"})
        lsrc = self.loop_source.get(pdb, {"status": "unavailable"})
        gsrc = self.ligand_source.get(f"{pdb}:{entity}", {"status": "unavailable"})
        return {"schema_version": SCHEMA_VERSION, "rfam_acc": self.acc, "seq_name": sn, "pdb_id": pdb, "entity_id": entity, "label_asym_id": inst["label_asym_id"],
                "auth_asym_id": inst["auth_asym_id"], "mapping_file": inst["mapping_file"],
                "sources": {"reference_structural_elements": {**esrc, "reference": f"{self.reference['reference']['pdb_id']}|1|{self.reference['reference']['auth_asym_id']}",
                                                              "note": "labels transferred by BGSU from the reference chain; for highlighting only"},
                            "motif_atlas_loops": {"status": lsrc["status"], "url": lsrc.get("loops", {}).get("url"), "motifs_url": lsrc.get("motifs", {}).get("url"),
                                                  "retrieved": lsrc.get("loops", {}).get("retrieved"), "n_loops_involving_chain": len(loops_here),
                                                  "motif_atlas_release": {t: v.get("release") for t, v in self.atlas_release.items()}},
                            "pdbe_binding_sites": {"status": gsrc["status"], "url": gsrc.get("url"), "retrieved": gsrc.get("retrieved"), "observation_level": "entity",
                                                   "n_entity_chains": len(self.idx.by_entity[(pdb, entity)]),
                                                   "chain_attribution": "resolved" if len(self.idx.by_entity[(pdb, entity)]) == 1 else "unresolved"}},
                "loops": [{"loop_id": L["loop_id"], "loop_type": L["loop_type"], "motif_group": L["motif_group"], "annotation_text": L["annotation_text"], "atlas_status": L["atlas_status"],
                           "multi_chain": L["multi_chain"], "chains": [c["auth_asym_id"] for c in L["chains"]], "consensus_positions": L["consensus_positions"]} for L in loops_here],
                "bound_molecules": [{"ccd_id": r["ccd_id"], "name": r["name"], "ligand_class": r["ligand_class"], "bound_molecule_ids": r["bound_molecule_ids"],
                                     "attribution": "chain" if r["chain_attribution"] == "resolved" else "entity_unresolved", "is_best_chain": r["best_chain_seq_name"] == sn,
                                     "in_default_track": r["in_default_track"], "flags": r["flags"], "consensus_positions": r["consensus_positions"]} for r in ligs],
                "residues": residues, "provenance": self._provenance()}

    # ------------------------------------------------------------------ 5. tracks
    def _loop_track(self, pos: int, loop_types: tuple[str, ...], observable: set[str], evaluated_inst: set[str]) -> dict[str, Any]:
        hits = [L for L in self.loops if L["loop_type"] in loop_types and pos in L["consensus_positions"]]
        nts = [(L, n) for L in hits for n in L["nucleotides"] if n["rfam_consensus_pos"] == pos and n["mapping_status"] in ("mapped", "mapped_no_coordinates")]
        instances = {n["seq_name"] for _, n in nts}
        entries = {L["pdb_id"] for L in hits}
        groups: dict[str, dict[str, set[str]]] = {}
        for L in hits:
            g = L["motif_group"] or "unassigned"
            d = groups.setdefault(g, {"loops": set(), "instances": set(), "entries": set()})
            d["loops"].add(L["loop_id"]); d["entries"].add(L["pdb_id"]); d["instances"].update(n["seq_name"] for LL, n in nts if LL["loop_id"] == L["loop_id"])
        conf = Counter(n["pp_category"] for _, n in {(n["seq_name"], n["label_seq_id"]): (L, n) for L, n in nts}.values())
        pos_obs_eval = observable & evaluated_inst
        positive_obs = instances & observable
        return {"n_loop_observations": len({L["loop_id"] for L in hits}), "n_instances": len(instances), "n_entries": len(entries),
                "n_distinct_motif_groups": len([g for g in groups if g != "unassigned"]),
                "by_motif_group": {g: {"n_loops": len(d["loops"]), "n_instances": len(d["instances"]), "n_entries": len(d["entries"])} for g, d in sorted(groups.items())},
                "n_loops_without_motif_group": len(groups.get("unassigned", {}).get("loops", ())),
                "loop_ids": sorted({L["loop_id"] for L in hits}),
                "confidence": {"high": conf.get("high", 0), "medium": conf.get("medium", 0), "low": conf.get("low", 0)},
                "denominators": {"n_observable_instances": len(observable), "n_source_evaluated_instances": len(evaluated_inst), "n_observable_and_evaluated": len(pos_obs_eval),
                                 "n_positive": len(instances), "n_positive_observable_and_evaluated": len(positive_obs & evaluated_inst),
                                 "n_not_annotated_by_source": len(pos_obs_eval - instances), "n_source_unavailable": len(observable - evaluated_inst),
                                 "n_positive_not_observable": len(instances - observable)},
                "frequency": _frac(len(positive_obs & evaluated_inst), len(pos_obs_eval))}

    def build_tracks(self) -> dict[str, Any]:
        # junction loop types: J3 and J4 always (separate tracks), plus any higher-order junction (J5, J6 ...) found in the loops
        jtypes = tuple(sorted({L["loop_type"] for L in self.loops if re.fullmatch(r"J\d+", L["loop_type"])} | {"J3", "J4"}, key=lambda t: int(t[1:])))
        n_inst = len(self.idx.instances)
        loop_eval = {i["seq_name"] for i in self.idx.instances if self.loop_source.get(i["pdb_id"], {}).get("status") == "evaluated"}
        elem_eval = {sn for sn, s in self.inst_element_source.items() if s["status"] == "evaluated"}
        lig_eval = {k for k, v in self.ligand_source.items() if v["status"] == "evaluated"}
        entity_chains = {f"{p}:{e}": [c["seq_name"] for c in cs] for (p, e), cs in self.idx.by_entity.items()}
        # per-ligand entry-level denominators: PDB entries with the bound ccd in the default track, and each entry's chains split by ligand-source status
        ccd_entries: dict[str, set[str]] = defaultdict(set)
        for r in self.ligands:
            if r["in_default_track"]:
                ccd_entries[r["ccd_id"]].add(r["pdb_id"])
        entry_chains_eval: dict[str, list[str]] = defaultdict(list)
        entry_chains_all: dict[str, list[str]] = defaultdict(list)
        for k, cs in entity_chains.items():
            entry_chains_all[k.split(":")[0]] += cs
            if k in lig_eval:
                entry_chains_eval[k.split(":")[0]] += cs
        positions: list[dict[str, Any]] = []
        for pos in range(1, self.clen + 1):
            observable = {i["seq_name"] for i in self.idx.instances if self.idx.observable(i["seq_name"], pos)}
            obs_entries = {self.idx.by_seq[s]["pdb_id"] for s in observable}
            # reference
            elem = self.reference["track"][pos - 1]
            labelled = {sn for sn, d in self.inst_elements.items() for v in d.values() if v["rfam_consensus_pos"] == pos}
            ref = {"element": elem, "reference": f"{self.reference['reference']['pdb_id']}|1|{self.reference['reference']['auth_asym_id']}",
                   "instance_transfer": {"n_instances_with_label": len(labelled), "n_instances_evaluated": len(elem_eval), "n_instances_unavailable": n_inst - len(elem_eval),
                                         "n_observable_instances": len(observable), "note": "transferred labels, not independent observations"}}
            hl = self._loop_track(pos, ("HL",), observable, loop_eval)
            il = self._loop_track(pos, ("IL",), observable, loop_eval)
            j3 = self._loop_track(pos, ("J3",), observable, loop_eval)
            j4 = self._loop_track(pos, ("J4",), observable, loop_eval)
            # every junction type present (J3, J4 and higher-order J5, J6 ... from BGSU), so no junction loop is left out of the junction track
            ju = self._loop_track(pos, jtypes, observable, loop_eval)
            by_type = {"J3": j3["n_instances"], "J4": j4["n_instances"], **{t: self._loop_track(pos, (t,), observable, loop_eval)["n_instances"] for t in jtypes if t not in ("J3", "J4")}}
            junction_union = {"n_loop_observations": ju["n_loop_observations"], "n_instances": ju["n_instances"], "n_entries": ju["n_entries"],
                              "denominators": ju["denominators"], "frequency": ju["frequency"], "loop_types": list(jtypes), "by_type_n_instances": by_type,
                              "rule": f"set union of {', '.join(jtypes)} loops/instances/entries (not a sum)"}
            # ligands (entity level)
            obs_entities = {k for k, cs in entity_chains.items() if any(c in observable for c in cs)}
            recs = [r for r in self.ligands if pos in r["consensus_positions"]]
            default = [r for r in recs if r["in_default_track"]]
            pos_ent = {f"{r['pdb_id']}:{r['entity_id']}" for r in default}
            by_ccd: dict[str, set[str]] = defaultdict(set)
            by_class: dict[str, set[str]] = defaultdict(set)
            for r in recs:
                by_class[r["ligand_class"]].add(f"{r['pdb_id']}:{r['entity_id']}")
                if r["in_default_track"]:
                    by_ccd[r["ccd_id"]].add(f"{r['pdb_id']}:{r['entity_id']}")
            conf = Counter(c["pp_category"] for r in default for c in r["contacts"] if c["rfam_consensus_pos"] == pos)
            obs_eval = obs_entities & lig_eval
            # entries (not entities: the chains of one entry count once) with each contacting ccd; an entry is eligible when a chain of a
            # ligand-source-evaluated entity has an observed nucleotide here; entries without one are not counted as negatives
            pos_entries: dict[str, set[str]] = defaultdict(set)
            for r in default:
                pos_entries[r["ccd_id"]].add(r["pdb_id"])
            by_ccd_entries = {}
            for c, pe in sorted(pos_entries.items()):
                elig = {e for e in ccd_entries[c] if any(sn in observable for sn in entry_chains_eval[e])}
                unobs = {e for e in ccd_entries[c] - elig if not any(sn in observable for sn in entry_chains_all[e])}
                by_ccd_entries[c] = {"n_positive": len(pe), "n_positive_observable_and_evaluated": len(pe & elig), "n_observable_and_evaluated": len(elig),
                                     "n_not_observable": len(unobs), "n_source_unavailable": len(ccd_entries[c]) - len(elig) - len(unobs),
                                     "frequency": _frac(len(pe & elig), len(elig))}
            lig = {"n_entry_entities": len(pos_ent), "n_entries": len({r["pdb_id"] for r in default}), "entry_entities": sorted(pos_ent),
                   "by_ccd_n_entry_entities": {c: len(s) for c, s in sorted(by_ccd.items())},
                   "by_ccd_entries": by_ccd_entries,
                   "by_class_n_entry_entities": {c: len(s) for c, s in sorted(by_class.items())},
                   "n_entry_entities_excluded_from_default": len({f"{r['pdb_id']}:{r['entity_id']}" for r in recs if not r["in_default_track"]}),
                   "confidence": {"high": conf.get("high", 0), "medium": conf.get("medium", 0), "low": conf.get("low", 0)},
                   "denominators": {"counting_unit": "pdb_entry_entity", "n_observable_entry_entities": len(obs_entities), "n_source_evaluated_entry_entities": len(lig_eval),
                                    "n_observable_and_evaluated": len(obs_eval), "n_positive": len(pos_ent), "n_positive_observable_and_evaluated": len(pos_ent & obs_eval),
                                    "n_not_annotated_by_source": len(obs_eval - pos_ent), "n_source_unavailable": len(obs_entities - lig_eval),
                                    "n_positive_not_observable": len(pos_ent - obs_entities)},
                   "frequency": _frac(len(pos_ent & obs_eval), len(obs_eval))}
            positions.append({"rfam_consensus_pos": pos, "observable": {"n_instances": len(observable), "n_entries": len(obs_entries), "n_entry_entities": len(obs_entities)},
                              "reference_structural_elements": ref, "hairpin_loops": hl, "internal_loops": il, "junction_loops_J3": j3, "junction_loops_J4": j4,
                              "junction_loops_union": junction_union, "ligand_binding_sites": lig})
        tracks_meta = {
            "reference_structural_elements": {"kind": "categorical_reference_transfer", "counting_unit": None, "source": self.reference["source"], "reference": self.reference["reference"]},
            "hairpin_loops": {"kind": "observed_3d_motif", "loop_types": ["HL"], "counting_unit": "loop observation / instance / pdb_entry", "denominator": "observable and source-evaluated instances"},
            "internal_loops": {"kind": "observed_3d_motif", "loop_types": ["IL"], "counting_unit": "loop observation / instance / pdb_entry", "denominator": "observable and source-evaluated instances"},
            "junction_loops_J3": {"kind": "observed_3d_motif", "loop_types": ["J3"], "counting_unit": "loop observation / instance / pdb_entry", "denominator": "observable and source-evaluated instances"},
            "junction_loops_J4": {"kind": "observed_3d_motif", "loop_types": ["J4"], "counting_unit": "loop observation / instance / pdb_entry", "denominator": "observable and source-evaluated instances"},
            "junction_loops_union": {"kind": "derived_union", "loop_types": list(jtypes), "rule": f"set union of all junction loop types present ({', '.join(jtypes)}); never the sum of the per-type counts"},
            "ligand_binding_sites": {"kind": "entity_level_contacts", "counting_unit": "pdb_entry_entity", "default_classes": self.config["default_track_classes"],
                                     "exclude_flagged_from_default_track": self.config.get("exclude_flagged_from_default_track", True),
                                     "denominator": "observable and source-evaluated (pdb_id, entity_id) pairs",
                                     "by_ccd_entries": "per bound ccd in contact at the position, counted per PDB entry: n_positive_observable_and_evaluated / "
                                                       "n_observable_and_evaluated = entries with a recorded contact / entries with the bound ccd that have an observed "
                                                       "nucleotide here in a chain of a source-evaluated entity; n_not_observable entries (no observed nucleotide in any "
                                                       "chain) and n_source_unavailable entries are not counted as negatives. An observed-contact frequency among "
                                                       "coordinate-observable entries, not a binding probability"},
            "motif_atlas_release": self.atlas_release}
        return {"schema_version": SCHEMA_VERSION, "rfam_acc": self.acc, "clen": self.clen, "n_instances": n_inst, "n_entries": len(self.idx.by_pdb), "n_entry_entities": len(self.idx.by_entity),
                "source_availability": {"reference_structural_elements": {"n_instances_evaluated": len(elem_eval), "n_instances_unavailable": n_inst - len(elem_eval)},
                                        "motif_atlas_loops": {"n_instances_evaluated": len(loop_eval), "n_instances_unavailable": n_inst - len(loop_eval)},
                                        "pdbe_binding_sites": {"n_entry_entities_evaluated": len(lig_eval), "n_entry_entities_unavailable": len(self.idx.by_entity) - len(lig_eval)}},
                "terminology": {"not_annotated_by_source": "observable and evaluated, but the source reports no annotation; not a confirmed biological absence",
                                "source_unavailable": "observable nucleotide but the source returned no record for this instance/entity",
                                "frequency": "n_positive_observable_and_evaluated / n_observable_and_evaluated"},
                "tracks": tracks_meta, "positions": positions, "provenance": self._provenance()}

    # ------------------------------------------------------------------ run
    def run(self, out_dir: Optional[Path] = None) -> dict[str, Any]:
        out = Path(out_dir) if out_dir else self.fam_dir
        log.info("reference structural elements")
        self.build_reference()
        log.info("motif atlas loops")
        loops = self.build_loops()
        log.info("bound molecules")
        ligands = self.build_ligands()
        log.info("instances")
        inst_sizes = 0
        for inst in self.idx.instances:
            rec = self.build_instance(inst)
            validate_json(rec, "annotation_instance")
            inst_sizes += _dump(out / "annotations" / f"{inst['seq_name']}.json", rec)
        log.info("tracks")
        tracks = self.build_tracks()
        validate_json(self.reference, "reference_structural_elements"); validate_json(loops, "loops"); validate_json(ligands, "ligands"); validate_json(tracks, "tracks")
        # per-chain loop files: what the Family Overview loads for the selected chain (the family-wide loops.json can be tens of MB for ribosomes)
        by_inst = 0
        for inst in self.idx.instances:
            sn = inst["seq_name"]
            rec = {"schema_version": SCHEMA_VERSION, "rfam_acc": self.acc, "seq_name": sn, "description": "full records of every BGSU loop with at least one nucleotide "
                   "in this chain (from loops.json; inter-chain loops included)", "loops": [L for L in self.loops if sn in L["instances"]], "provenance": self._provenance()}
            validate_json(rec, "loops_instance")
            by_inst += _dump(out / "loops_by_instance" / f"{sn}.json", rec)
        sizes = {"reference_structural_elements.json": _dump(out / "reference_structural_elements.json", self.reference), "loops.json": _dump(out / "loops.json", loops), "loops_by_instance/ (total)": by_inst,
                 "ligands.json": _dump(out / "ligands.json", ligands), "tracks.json": _dump(out / "tracks.json", tracks), "annotations/ (total)": inst_sizes}
        summary = {"rfam_acc": self.acc, "generated": self.generated, "n_instances": len(self.idx.instances), "n_entries": len(self.idx.by_pdb),
                   "reference_structural_elements": {k: self.reference[k] for k in ("n_consensus_positions_annotated",)} | {"reference": self.reference["reference"]["seq_name"],
                                                                                                                             "instance_transfer": {k: v for k, v in self.reference["instance_transfer"].items() if k != "disagreements"},
                                                                                                                             "n_conflicts": len(self.reference["conflicts"]),
                                                                                                                             "n_reference_nucleotides_not_on_consensus": len(self.reference["reference_nucleotides_not_on_consensus"])},
                   "loops": loops["counts"], "ligands": ligands["counts"], "unresolved_identifiers": self.unresolved, "n_unresolved_identifiers": len(self.unresolved),
                   "file_sizes_bytes": sizes, "provenance": self._provenance()}
        _dump(out / "annotation_summary.json", summary)
        return summary


def build_annotations(acc: str, cache_dir: Path, out_dir: Path, raw_dir: Optional[Path] = None, offline: bool = False) -> dict[str, Any]:
    from .cache import Cache
    fam_dir = out_dir / "families" / acc
    if not (fam_dir / "instances.json").exists():
        raise FileNotFoundError(f"family {acc} not built under {fam_dir}; run `build` first")
    fetcher = AnnotationFetcher(Cache(cache_dir), raw_dir, offline=offline)
    return AnnotationBuilder(fam_dir, fetcher).run()
