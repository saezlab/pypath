"""Phenotype And Trait Ontology terms: names, synonyms, definitions, hierarchy and relationships."""

from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import pato
from pypath.inputs_v2.parsers.obo import parse_obo_text

# Stanzas copied from pato.obo (data-version releases/2026-10-01/pato.obo); the stub of a Cell
# Ontology term is the one pato-base.obo holds, which pato.obo does not.
PATO_OBO = """format-version: 1.2
data-version: releases/2026-10-01/pato.obo
ontology: pato

[Term]
id: PATO:0000169
name: viability
def: "An organismal quality inhering in a bearer or a population by virtue of the bearer's disposition to survive and develop normally or the number of surviving individuals in a given population." [PATOC:GVG]
subset: attribute_slim
subset: disposition_slim
is_a: PATO:0001995 ! organismal quality

[Term]
id: PATO:0000001
name: characteristic
alt_id: PATO:0000072
def: "The properties, attributes, or modifiers used to describe an entity or process or some aspect of an entity or process." [https://github.com/pato-ontology/pato/issues/602]
synonym: "quality" EXACT [OBO:legacy]

[Term]
id: PATO:0000300
name: immobile
def: "A mobility quality inhering in a bearer by virtue of the bearer's being incapable to move freely." [PATOC:GVG]
subset: disposition_slim
subset: value_slim
synonym: "fixed" RELATED []
is_a: PATO:0000004 ! mobility
relationship: is_opposite_of PATO:0000299 ! mobile

[Term]
id: PATO:0000299
name: mobile
def: "A mobility quality of inhering in a bearer by virtue of the bearer's disposition to move freely." [PATOC:GVG]
subset: disposition_slim
subset: value_slim
is_a: PATO:0000004 ! mobility
disjoint_from: PATO:0000300 ! immobile
relationship: is_opposite_of PATO:0000300 ! immobile

[Term]
id: PATO:0001304
name: variability of temperature
def: "A variability quality inhering in a bearer by virtue of whether the bearer exhibits temperature variation or change." [PATO:GVG]
subset: attribute_slim
is_a: PATO:0001303 ! variability
relationship: RO:0002503 PATO:0000146 ! towards temperature

[Term]
id: PATO:0000000
name: obsolete pato
is_obsolete: true

[Term]
id: CL:0000000
relationship: BFO:0000051 GO:0005634 {gci_filler="PATO:0001908", gci_relation="RO:0000053"}
"""


def terms():
    return {row["id"]: row for row in parse_obo_text(PATO_OBO)}


def extract(row):
    record = pato.terms_schema(row)
    extractor = SilverExtractor("pato", "terms")
    extractor.process_record(record, raw_payload={}, row_id="1", row_number=1)
    return record, extractor


def edges(extractor):
    return sorted(
        (relation.predicate, extractor.entities[relation.object_entity_key].identifier)
        for relation in extractor.relations
    )


def term_entity(extracted, identifier):
    (entity,) = (e for e in extracted.entities.values() if e.identifier == identifier and e.label)
    return entity


def test_term_has_name_definition_and_parent():
    record, extracted = extract(terms()["PATO:0000169"])
    entity = term_entity(extracted, "PATO:0000169")
    # the identifier keeps the prefix, exactly as BioGRID cites it
    assert (entity.entity_type, entity.namespace) == ("ontology_class", "pato")
    assert entity.label == "viability"
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert {("pato", "PATO:0000169"), ("name", "viability")} <= identifiers
    description = {format_term(a.term): a.value for a in record.annotations}["description"]
    assert description.startswith("An organismal quality inhering in a bearer or a population")
    assert edges(extracted) == [("subclass_of", "PATO:0001995")]
    assert {e.namespace for e in extracted.entities.values()} == {"pato"}


def test_synonym_and_alt_id():
    _, extracted = extract(terms()["PATO:0000001"])
    identifiers = {(i["ns"], i["id"]) for i in term_entity(extracted, "PATO:0000001").identifiers}
    assert {("synonym", "quality"), ("pato", "PATO:0000072")} <= identifiers
    assert edges(extracted) == []


def test_opposite_qualities_are_one_symmetric_edge_each_way():
    _, extracted = extract(terms()["PATO:0000300"])
    assert edges(extracted) == [("opposite_of", "PATO:0000299"), ("subclass_of", "PATO:0000004")]
    _, extracted = extract(terms()["PATO:0000299"])
    assert edges(extracted) == [("opposite_of", "PATO:0000300"), ("subclass_of", "PATO:0000004")]


def test_relationship_type_biolink_maps_broadly_keeps_the_original():
    record, extracted = extract(terms()["PATO:0001304"])  # RO:0002503, "towards"
    assert edges(extracted) == [("actively_involved_in", "PATO:0000146"), ("subclass_of", "PATO:0001303")]
    (relation,) = (r for r in record.ontology_relations if r.annotations)
    assert [(format_term(a.term), a.value) for a in relation.annotations] == [("original_predicate", "RO:0002503")]


def test_obsolete_and_foreign_terms_give_nothing():
    assert pato.terms_schema(terms()["PATO:0000000"]) is None
    assert pato.terms_schema(terms()["CL:0000000"]) is None
    assert pato.terms_schema({}) is None


def test_resource_declares_the_ontology():
    assert pato.config.id == ResourceCv.PATO
    assert pato.config.license == LicenseCV.CC_BY_3_0
    assert pato.resource.terms.mapper is pato.terms_schema
    assert Namespace.PATO.value == "pato"
