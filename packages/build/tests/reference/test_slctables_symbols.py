from pypath.inputs_v2.slctables import _protein_gene_name


def test_slc_symbol_evidence_is_exact_and_human_normalized():
    assert _protein_gene_name({"Protein name": "SV2A"}) == ["SV2A"]
    assert _protein_gene_name({"Protein name": "RhAG"}) == ["RHAG"]
    for name in ["", "pseudogene", "SLC1A1 pseudogene", "EAAC1, EAAT3", "System X-AG"]:
        assert _protein_gene_name({"Protein name": name}) == []
