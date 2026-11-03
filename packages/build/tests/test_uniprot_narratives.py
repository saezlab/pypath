from omnipath_core.biolink import annotation_term
from pypath.inputs_v2.uniprot import _protein_description_pairs
from pypath.internals.tabular_builder import AnnotationsBuilder, CV


def test_uniprot_narratives_keep_categories_and_clean_prose():
    row = {
        "Function [CC]": "FUNCTION: Uses Ca(2+) (By similarity) (PubMed:123). {ECO:0000250|UniProtKB:P1}.",
        "Involvement in disease": "DISEASE: A [MIM:123]: Text. Note=More.",
        "Pathway": "PATHWAY: Shared text.",
        "Activity regulation": "ACTIVITY REGULATION: Shared text.",
    }
    annotations = AnnotationsBuilder(CV.from_pairs(_protein_description_pairs)).build(row)
    actual = {annotation_term(a.term): a.value for a in annotations}
    assert actual == {
        "up:Function_Annotation": "Uses Ca(2+) (By similarity).",
        "up:Disease_Annotation": "A: Text. More.",
        "up:Pathway_Annotation": "Shared text.",
        "up:Activity_Regulation_Annotation": "Shared text.",
    }
