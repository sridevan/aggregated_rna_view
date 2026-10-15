from rna_explorer.confidence import consensus_statistics
from rna_explorer.mapping import map_aligned_sequence, build_chain_mapping
from rna_explorer.sequences import ChainInstance
from rna_explorer.stockholm import type_columns


def _inst(pdb: str, seq: str, observed: set[int]) -> ChainInstance:
    residues = [{"residue_number": i + 1, "residue_name": seq[i], "author_residue_number": i + 1, "author_insertion_code": "", "observed_ratio": 1.0 if i + 1 in observed else 0.0} for i in range(len(seq))]
    return ChainInstance(pdb, 1, "A", "A", 1, len(seq), 1, len(seq), seq, seq, residues, {}, observed)


def test_sequence_vs_structure_level() -> None:
    cm = type_columns("xxxx.x")
    a1 = map_aligned_sequence("u1", "ACGUaA", "****9*", cm)      # unique seq 1 (two chains share it)
    a2 = map_aligned_sequence("u2", "AC-U.A", "**.8.7", cm)      # unique seq 2
    c1 = build_chain_mapping(_inst("1aaa", "ACGUAA", {1, 2, 3, 4, 5, 6}), a1, cm, "RFX")
    c2 = build_chain_mapping(_inst("1aab", "ACGUAA", {1, 2}), a1, cm, "RFX")
    c3 = build_chain_mapping(_inst("1aac", "ACUA", {1, 2, 3, 4}), a2, cm, "RFX")
    s = consensus_statistics(cm, {"k1": a1, "k2": a2}, [c1, c2, c3], {"rfam_acc": "RFX"})
    sl, st = s["sequence_level"], s["structure_level"]
    assert sl["n_unique_sequences"] == 2 and st["n_instances"] == 3
    assert sl["n_match"] == [2, 2, 1, 2, 2] and sl["n_deletion"] == [0, 0, 1, 0, 0] and sl["coverage"] == [1.0, 1.0, 0.5, 1.0, 1.0]
    assert st["n_match"] == [3, 3, 2, 3, 3] and st["n_match_with_coordinates"] == [3, 3, 1, 2, 2]
    assert st["structural_coverage"][2] == round(1 / 3, 4)
    assert sl["pp_mean"][3] == round((0.975 + 0.8) / 2, 4) and sl["pp_low"][4] == 1 and sl["pp_high"][0] == 2
    assert "composite" not in json_dumps(s)


def json_dumps(o: object) -> str:
    import json
    return json.dumps(o)
