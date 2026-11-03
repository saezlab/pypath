from pypath.inputs_v2.kegg import pathways_schema
from omnipath_build.silver import SilverExtractor


def test_kegg_gene_local_id_is_not_asserted_as_entrez():
    row = dict(
        pathway_id="hsa00010",
        pathway_name="test",
        taxon_id="9606",
        protein_member_kegg_ids="hsa:64689",
        protein_member_entrez_ids="64689",
        protein_member_uniprot_ids="F6S8C4",
        protein_member_names="test",
        protein_member_reaction_ids="R00001",
        reaction_ids="R00001",
    )
    ex = SilverExtractor("kegg", "pathways")
    ex.process_record(pathways_schema(row), row, "pathways:0", 0)
    protein = next(o for o in ex.entities.values() if o.entity_type == "protein")
    assert protein.namespace == "uniprot" and protein.identifier == "F6S8C4"
    assert all(x["ns"] != "entrez" for x in protein.identifiers)
    assert row["protein_member_entrez_ids"] == "64689"
