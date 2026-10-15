from pathlib import Path
import pytest
from rna_explorer.mmcif import atom_site_observed, cif_category, entity_source, poly_seq_scheme

CIF = Path(__file__).resolve().parent.parent / "poc_rf00167" / "data" / "6p2h_updated.cif"
pytestmark = pytest.mark.skipif(not CIF.exists(), reason="POC mmCIF not present")


def test_scheme_and_observed() -> None:
    scheme = poly_seq_scheme(CIF, "A", 1)
    assert len(scheme) == 69 and scheme[1]["mon_id"] == "GTP" and scheme[1]["pdb_seq_num"] == "17"
    assert atom_site_observed(CIF, "A", 1) == set(range(1, 70))
    assert entity_source(CIF, 1)["cif_organism"] == "synthetic construct"
    assert cif_category(CIF, "_struct")[0]["title"]
