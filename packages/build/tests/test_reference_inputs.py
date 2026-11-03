"""Reference parsers preserve source meaning through the strict Silver boundary."""

from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from omnipath_build.silver import SilverExtractor
from pypath.inputs_v2 import chebi, uniprot


def extract(record):
    extractor = SilverExtractor("test", "test")
    extractor.process_record(record, raw_payload={}, row_id="1", row_number=1)
    return extractor


def test_chebi_properties_and_relations():
    row = {
        "chebi_id": "CHEBI:1",
        "name": "chemical class",
        "formula": "H2O",
        "mass": "18",
        "charge": "0",
        "is_a": ["CHEBI:2"],
        "relationships": [
            {"type": "BFO:0000051", "target": "CHEBI:3"},
            {"type": "RO:0018033", "target": "CHEBI:4"},
            {"type": "RO:0000087", "target": "CHEBI:5"},
            {"type": "RO:0018039", "target": "CHEBI:6"},
        ],
    }
    record = chebi.molecules_schema(row)
    assert record.type == "chemical_entity"
    annotations = {format_term(a.term): a.value for a in record.annotations}
    assert annotations["has_chemical_formula"] == "H2O"
    assert annotations["chemrof:mass"] == "18"
    # A conjugate base relation is directional; no Biolink predicate keeps that.
    assert annotations["RO:0018033"] == "CHEBI:4"
    assert "RO:0000087" not in annotations and "RO:0018039" not in annotations
    relations = extract(record).relations
    assert {r.predicate: [(a["term"], a["value"]) for a in r.annotations] for r in relations} == {
        "subclass_of": [],
        "has_part": [],
        "has_chemical_role": [],
        "chemically_similar_to": [("original_predicate", "RO:0018039")],
    }
    assert {r.statement_kind for r in relations} == {"ontology"}
    assert (
        list(chebi._id_translation_rows({"chebi_id": "CHEBI:1", "inchi": "InChI=1S/test"}))[0][
            "key_type"
        ]
        == Namespace.CHEBI
    )


def test_uniprot_associations_and_narratives():
    record = uniprot.proteins_schema(
        {
            "Entry": "P1",
            "Sequence": "MA",
            "Organism (ID)": "9606",
            "Gene Ontology IDs": "GO:0005739; GO:0003674",
            "Keyword ID": "KW-0001",
            "PubMed ID": "123;456",
            "Subcellular location [CC]": "Nucleus; evidence uncertain",
            "Mass": "100",
        }
    )
    assert record.type == "protein"
    annotations = [(format_term(a.term), a.value) for a in record.annotations]
    assert ("in_taxon", "NCBITaxon:9606") in annotations
    assert ("publications", "PMID:123") in annotations
    assert ("up:Subcellular_Location_Annotation", "Nucleus; evidence uncertain") in annotations
    extracted = extract(record)
    assert len(extracted.relations) == 3
    assert {r.predicate for r in extracted.relations} == {"associated_with"}
    assert all(e.entity_type in {"protein", "ontology_class"} for e in extracted.entities.values())


def test_keywords_and_missing_ensembl_columns():
    record = uniprot.keyword_terms_schema({"id": "KW-0001", "name": "keyword", "is_a": ["KW-9993"]})
    assert record.type == "ontology_class"
    assert extract(record).relations[0].predicate == "subclass_of"
    assert uniprot.keyword_terms_schema({"id": "KW-0001", "is_obsolete": True}) is None
    assert uniprot._parse_ensembl_dr_line("DR   Ensembl; ENST1.2; -; ENSG1.3.") == [
        (Namespace.ENSG, "ENSG1"),
        (Namespace.ENSG, "ENSG1.3"),
        (Namespace.ENST, "ENST1"),
        (Namespace.ENST, "ENST1.2"),
    ]
