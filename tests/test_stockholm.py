"""Column typing and residue states on a synthetic alignment: insert before first match (anchor 0), '~' local end, deletion vs uncovered, trailing insertion."""
from pathlib import Path
from rna_explorer.stockholm import parse_stockholm, type_columns
from rna_explorer.mapping import map_aligned_sequence, pp_category, pp_value

STO = """# STOCKHOLM 1.0

s1              ggACG~~U-Aaa
#=GR s1 PP      **999..*.*87
s2              ..ACGcc-GA..
#=GR s2 PP      ..***99.**..
#=GC SS_cons    ..<<<~~>>>..
#=GC RF         ..xxx~~xxx..
//
"""


def test_columns_and_states(tmp_path: Path) -> None:
    p = tmp_path / "t.sto"
    p.write_text(STO)
    aln = parse_stockholm(p)
    cm = type_columns(aln.gc["RF"], aln.gc["SS_cons"])
    assert cm.clen == 6 and cm.col_type[:2] == ["insert", "insert"] and cm.col_type[5:7] == ["local_end", "local_end"]
    s1 = map_aligned_sequence("s1", aln.seqs["s1"], aln.gr["s1"]["PP"], cm)
    states = [(r.index, r.aligned_char, r.state, r.consensus_pos, r.anchor, r.offset) for r in s1.residues]
    assert states[0] == (1, "g", "insert", None, 0, 1) and states[1] == (2, "g", "insert", None, 0, 2)
    assert states[2] == (3, "A", "match", 1, None, None)
    assert states[5] == (6, "U", "match", 4, None, None)            # after the ~~ columns (gaps for s1)
    assert states[-2] == (8, "a", "insert", None, 6, 1) and states[-1] == (9, "a", "insert", None, 6, 2)
    cons = {c.consensus_pos: c.state for c in s1.consensus}
    assert cons == {1: "match", 2: "match", 3: "match", 4: "match", 5: "deletion", 6: "match"}
    s2 = map_aligned_sequence("s2", aln.seqs["s2"], aln.gr["s2"]["PP"], cm)
    st2 = {r.index: (r.state, r.consensus_pos, r.anchor, r.offset) for r in s2.residues}
    assert st2[4] == ("local_end", None, 3, 1) and st2[5] == ("local_end", None, 3, 2)
    cons2 = {c.consensus_pos: c.state for c in s2.consensus}
    assert cons2[4] == "deletion"
    assert [r.pp_char for r in s2.residues] == ["*", "*", "*", "9", "9", "*", "*"]


def test_uncovered_vs_deletion() -> None:
    cm = type_columns("xxxxxx")
    s = map_aligned_sequence("t", "--AC-G", "..**.*", cm)
    assert {c.consensus_pos: c.state for c in s.consensus} == {1: "uncovered", 2: "uncovered", 3: "match", 4: "match", 5: "deletion", 6: "match"}
    assert (s.cons_min, s.cons_max) == (3, 6)


def test_pp_bins() -> None:
    assert pp_value("*") == 0.975 and pp_value("9") == 0.9 and pp_value("0") == 0.025 and pp_value("") is None
    assert pp_category("*") == "high" and pp_category("9") == "medium" and pp_category("8") == "medium" and pp_category("7") == "low"
