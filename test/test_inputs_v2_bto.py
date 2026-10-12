"""BRENDA Tissue Ontology terms: names, synonyms, definitions and hierarchy."""

from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import bto
from pypath.inputs_v2.parsers.obo import parse_obo_text

# Stanzas copied from bto.obo (data-version releases/2021-10-26); one foreign stanza added.
BTO_OBO = """format-version: 1.2
data-version: releases/2021-10-26
ontology: bto

[Term]
id: BTO:0000007
name: HEK-293 cell
def: "Established from a human primary embryonal kidney transformed by adenovirus type 5." [Deutsche_Sammlung_von_Mikroorganismen_und_Zellkulturen_GmbH:DSMZ]
synonym: "293 cell" RELATED []
synonym: "A-293 cell" RELATED []
synonym: "HEK cell" RELATED []
is_a: BTO:0000067 ! kidney cell line
is_a: BTO:0002733 ! embryonic kidney cell line

[Term]
id: BTO:0000067
name: kidney cell line
synonym: "renal cell line" RELATED []
relationship: develops_from BTO:0000671 ! kidney

[Term]
id: BTO:0000671
name: kidney
is_a: BTO:0000431 ! excretory gland
relationship: part_of BTO:0001244 ! urinary tract

[Term]
id: BTO:0001416
name: umbilical vein endothelium
alt_id: BTO:0001414
is_a: BTO:0000393 ! endothelium
relationship: part_of BTO:0001509 ! umbilical vein

[Term]
id: BTO:0000068
name: amniotic fluid
relationship: contained_in BTO:0000025 ! amniotic cavity

[Term]
id: BTO:0000311
name: culture filtrate
comment: Added as synonym to culture fluid.
is_obsolete: true

[Term]
id: UBERON:0002113
name: kidney
"""


def terms():
    return {row["id"]: row for row in parse_obo_text(BTO_OBO)}


def extract(row):
    record = bto.terms_schema(row)
    extractor = SilverExtractor("bto", "terms")
    extractor.process_record(record, raw_payload={}, row_id="1", row_number=1)
    return record, extractor


def edges(extractor):
    return sorted(
        (relation.predicate, extractor.entities[relation.object_entity_key].identifier)
        for relation in extractor.relations
    )


def test_term_has_name_synonyms_definition_and_parents():
    record, extracted = extract(terms()["BTO:0000007"])
    (entity,) = (e for e in extracted.entities.values() if e.identifier == "BTO:0000007")
    assert (entity.entity_type, entity.namespace) == ("ontology_class", "bto")
    assert entity.label == "HEK-293 cell"
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert {("name", "HEK-293 cell"), ("synonym", "293 cell"), ("synonym", "HEK cell")} <= identifiers
    assert {format_term(a.term): a.value for a in record.annotations} == {
        "description": "Established from a human primary embryonal kidney transformed by adenovirus type 5."
    }
    assert edges(extracted) == [("subclass_of", "BTO:0000067"), ("subclass_of", "BTO:0002733")]
    assert {e.namespace for e in extracted.entities.values()} == {"bto"}


def test_part_of_and_develops_from_are_edges():
    _, extracted = extract(terms()["BTO:0000067"])
    assert edges(extracted) == [("develops_from", "BTO:0000671")]
    _, extracted = extract(terms()["BTO:0000671"])
    assert edges(extracted) == [("part_of", "BTO:0001244"), ("subclass_of", "BTO:0000431")]


def test_alt_id_is_an_identifier_of_the_term():
    _, extracted = extract(terms()["BTO:0001416"])
    (entity,) = (e for e in extracted.entities.values() if e.identifier == "BTO:0001416")
    assert ("bto", "BTO:0001414") in {(i["ns"], i["id"]) for i in entity.identifiers}


def test_other_relationship_types_are_left_out():
    _, extracted = extract(terms()["BTO:0000068"])
    assert edges(extracted) == []


def test_obsolete_and_foreign_terms_give_nothing():
    assert bto.terms_schema(terms()["BTO:0000311"]) is None
    assert bto.terms_schema(terms()["UBERON:0002113"]) is None
    assert bto.terms_schema({}) is None


def test_resource_declares_the_ontology():
    assert bto.config.id == ResourceCv.BTO
    assert bto.config.license == LicenseCV.CC_BY_4_0
    assert bto.resource.terms.mapper is bto.terms_schema
    assert Namespace.BTO.value == "bto"
