"""Curated-selection builds (e.g. the 60-chain RF00177 set): selection CSV parsing and the validator scheme rebuilt from the PDBe residue listing."""
from pathlib import Path
import pytest
from rna_explorer.pdbe import scheme_from_residue_listing
from rna_explorer.pipeline import read_selection

SEL = Path(__file__).resolve().parents[1] / "data/selections/RF00177_60_chains.csv"


def test_read_selection_16s_columns_and_generic_columns(tmp_path):
    p = tmp_path / "s.csv"
    p.write_text("pdb_id,ssu_auth_chain,ssu_label_asym_id,notes\n4V9D,AA,A,x\n4v9d,BA,Y,y\n")
    assert read_selection(p) == [("4v9d", "A", "AA"), ("4v9d", "Y", "BA")]
    q = tmp_path / "g.csv"; q.write_text("pdb_id,label_asym_id\n7k00,A\n")
    assert read_selection(q) == [("7k00", "A", None)]


def test_read_selection_rejects_missing_ids_and_duplicates(tmp_path):
    p = tmp_path / "bad.csv"; p.write_text("pdb_id,label_asym_id\n7k00,\n")
    with pytest.raises(ValueError, match="needs pdb_id and label_asym_id"):
        read_selection(p)
    d = tmp_path / "dup.csv"; d.write_text("pdb_id,label_asym_id\n7k00,A\n7K00,A\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_selection(d)


@pytest.mark.skipif(not SEL.exists(), reason="RF00177 selection not present")
def test_rf00177_selection_is_the_approved_60_chain_set():
    rows = read_selection(SEL)
    assert len(rows) == 60 and len({p for p, _, _ in rows}) == 58
    assert ("7k00", "A", "A") in rows and ("4v9d", "A", "AA") in rows


def test_scheme_from_residue_listing_marks_unobserved_and_keeps_insertion_codes():
    listing = {"molecules": [{"entity_id": 1, "chains": [{"struct_asym_id": "A", "residues": [
        {"residue_number": 1, "author_residue_number": 5, "author_insertion_code": "", "residue_name": "G", "observed_ratio": 1.0},
        {"residue_number": 2, "author_residue_number": 5, "author_insertion_code": "A", "residue_name": "PSU", "observed_ratio": 1.0},
        {"residue_number": 3, "author_residue_number": 6, "author_insertion_code": None, "residue_name": "C", "observed_ratio": 0}]}]}]}
    s = scheme_from_residue_listing(listing, 1, "A", "AA")
    assert s[1] == {"pdb_strand_id": "AA", "asym_id": "A", "entity_id": "1", "pdb_seq_num": "5", "pdb_ins_code": "", "mon_id": "G", "auth_seq_num": "5"}
    assert s[2]["pdb_ins_code"] == "A" and s[2]["mon_id"] == "PSU"
    assert s[3]["auth_seq_num"] == "?" and s[3]["pdb_ins_code"] == ""
    with pytest.raises(ValueError):
        scheme_from_residue_listing(listing, 1, "B", "BB")


def test_cache_retries_transient_download_failures(tmp_path, monkeypatch):
    import requests
    from rna_explorer import cache as cache_mod
    calls = {"n": 0}
    class R:
        content = b"ok"
        def raise_for_status(self): pass
    def flaky(*a, **k):
        calls["n"] += 1
        if calls["n"] < 3:
            raise requests.exceptions.ChunkedEncodingError("Response ended prematurely")
        return R()
    monkeypatch.setattr(cache_mod.requests, "request", flaky); monkeypatch.setattr(cache_mod, "RETRY_WAIT_S", 0)
    p = cache_mod.Cache(tmp_path).fetch_bytes("https://example.org/x.cif", ".cif")
    assert p.read_bytes() == b"ok" and calls["n"] == 3
    calls["n"] = -10   # always failing: gives up after RETRIES attempts
    monkeypatch.setattr(cache_mod.requests, "request", lambda *a, **k: (_ for _ in ()).throw(requests.ConnectionError("reset")))
    with pytest.raises(requests.ConnectionError):
        cache_mod.Cache(tmp_path).fetch_bytes("https://example.org/y.cif", ".cif")


def test_display_representative_override(tmp_path):
    import json
    from rna_explorer.pipeline import apply_display_override
    rows = [{"seq_name": "a_1_A_A", "x": 1}, {"seq_name": "b_1_B_B", "x": 2}]
    auto = {**rows[0], "selection_rule": "rule"}
    p = tmp_path / "o.json"; p.write_text(json.dumps({"RFX": {"display_representative": {"seq_name": "b_1_B_B", "reason": "why", "chosen_by": "user", "date": "2026-10-17"}}}))
    d = apply_display_override("RFX", rows, auto, p)
    assert d["seq_name"] == "b_1_B_B" and d["automatic_choice"] == "a_1_A_A" and "why" in d["selection_rule"]
    assert apply_display_override("RFY", rows, auto, p) is auto                                   # no override: the rule's pick
    p.write_text(json.dumps({"RFX": {"display_representative": {"seq_name": "zzz", "reason": "r"}}}))
    with pytest.raises(ValueError):
        apply_display_override("RFX", rows, auto, p)


def test_selection_families_use_their_chain_files(tmp_path):
    import json
    from rna_explorer.pdbe import family_uses_chain_files, instance_cif
    fam = tmp_path / "F"; fam.mkdir()
    (fam / "family.json").write_text(json.dumps({"provenance": {"coordinate_source": "PDBe ModelServer per chain (first model); ..."}}))
    assert family_uses_chain_files(fam) is True
    (fam / "family.json").write_text(json.dumps({"provenance": {}})); assert family_uses_chain_files(fam) is False
    class C:
        def fetch_bytes(self, url, suffix=".bin"): return url
    inst = {"pdb_id": "7K00", "label_asym_id": "A"}
    assert instance_cif(C(), inst, True) == "https://www.ebi.ac.uk/pdbe/model-server/v1/7k00/atoms?label_asym_id=A&model_nums=1&encoding=cif"
    assert instance_cif(C(), inst, False).endswith("/7k00_updated.cif")


def test_stream_block_parses_and_keeps_nothing(tmp_path, monkeypatch):
    import hashlib, tempfile
    from rna_explorer import modifications as M
    body = b"data_x\n_entry.id X\n_chem_comp.id PSU\n"
    class R:
        content = body
        def raise_for_status(self): pass
    monkeypatch.setattr("requests.get", lambda *a, **k: R())
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    block, src = M._stream_block("https://example.org/x_updated.cif")
    assert block.find_value("_entry.id") == "X" and src["kept"] is False and src["bytes"] == len(body) and src["md5"] == hashlib.md5(body).hexdigest()
    assert list(tmp_path.iterdir()) == []                                          # the temporary file is gone
