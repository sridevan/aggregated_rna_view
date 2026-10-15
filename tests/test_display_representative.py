from rna_explorer.instances import choose_display_representative


def row(**kw):
    base = dict(pdb_id="1abc", auth_asym_id="A", is_clan_best=True, experimental_method="X-ray diffraction", pdbe_2d_available=True, organism_species="Vibrio vulnificus",
                n_match=67, n_match_with_coordinates=67, fraction_clen_matched=0.65, resolution=2.0)
    base.update(kw); return base


def test_rule_prefers_xray_natural_complete_then_resolution():
    rows = [row(pdb_id="9ndd", experimental_method="Electron Microscopy", organism_species=None, fraction_clen_matched=0.96, n_match=98, n_match_with_coordinates=98, resolution=2.67),
            row(pdb_id="4tzx", resolution=2.01, fraction_clen_matched=0.677), row(pdb_id="1y26", resolution=2.1), row(pdb_id="6p2h", resolution=1.9, fraction_clen_matched=0.54),
            row(pdb_id="3rkf", pdbe_2d_available=False, resolution=1.5), row(pdb_id="4fe5", organism_species="synthetic construct", resolution=1.4),
            row(pdb_id="2b57", n_match_with_coordinates=50, resolution=1.6)]
    best = choose_display_representative(rows)
    assert best is not None and best["pdb_id"] == "4tzx" and best["eligible_count"] == 2
    assert choose_display_representative([row(experimental_method="Electron Microscopy")]) is None
