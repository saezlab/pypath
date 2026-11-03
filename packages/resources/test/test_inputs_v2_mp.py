"""Mammalian Phenotype Ontology terms: names, synonyms, definitions and hierarchy."""

from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import mp
from pypath.inputs_v2.parsers.obo import parse_obo_text

# Stanzas copied from mp.obo (data-version releases/2026-09-23/mp.obo); one foreign stanza added.
MP_OBO = """format-version: 1.2
data-version: releases/2026-09-23/mp.obo
ontology: mp

[Term]
id: MP:0001262
name: decreased body weight
def: "lower than normal average weight" [PMID:10709991]
subset: Europhenome_Terms
subset: IMPC
subset: Sanger_Terms
synonym: "decreased weight" EXACT []
synonym: "low body weight" EXACT []
synonym: "reduced body weight" EXACT []
synonym: "reduced weight" EXACT []
is_a: MP:0001259 ! abnormal body weight
is_a: MP:0001265 ! decreased body size
is_a: MP:0012322 ! decreased total tissue mass

[Term]
id: MP:0000003
name: abnormal adipose tissue morphology
alt_id: MP:0000011
def: "any structural anomaly of the connective tissue composed of fat cells enmeshed in areolar tissue" [MESH:A10.165.114, MGI:cwg]
synonym: "abnormality of adipose tissue" BROAD []
synonym: "adipose tissue abnormalities" BROAD []
synonym: "adipose tissue dysplasia" NARROW []
is_a: MP:0005375 ! adipose tissue phenotype

[Term]
id: MP:0000012
name: obsolete loss of subcutaneous adipose tissue
def: "OBSOLETE. reduction in amount or absence of adipose tissue beneath the skin" [PMID:8941642]
comment: Obsolete.  Use instead the terms decreased subcutaneous adipose tissue amount (MP:0008844) or absent subcutaneous adipose tissue (MP:0008843).
is_obsolete: true
consider: MP:0008843
consider: MP:0008844

[Term]
id: PATO:0000169
name: viability
"""


def terms():
    return {row["id"]: row for row in parse_obo_text(MP_OBO)}


def extract(row):
    record = mp.terms_schema(row)
    extractor = SilverExtractor("mp", "terms")
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


def test_term_has_name_synonyms_definition_and_parents():
    record, extracted = extract(terms()["MP:0001262"])
    entity = term_entity(extracted, "MP:0001262")
    # the identifier keeps the prefix, exactly as BioGRID cites it
    assert (entity.entity_type, entity.namespace) == ("ontology_class", "mp")
    assert entity.label == "decreased body weight"
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert {
        ("mp", "MP:0001262"),
        ("name", "decreased body weight"),
        ("synonym", "decreased weight"),
        ("synonym", "low body weight"),
    } <= identifiers
    assert {format_term(a.term): a.value for a in record.annotations}["description"] == "lower than normal average weight"
    assert edges(extracted) == [
        ("subclass_of", "MP:0001259"),
        ("subclass_of", "MP:0001265"),
        ("subclass_of", "MP:0012322"),
    ]
    assert {e.namespace for e in extracted.entities.values()} == {"mp"}


def test_alt_id_is_an_identifier_of_the_term():
    _, extracted = extract(terms()["MP:0000003"])
    entity = term_entity(extracted, "MP:0000003")
    assert entity.label == "abnormal adipose tissue morphology"
    assert ("mp", "MP:0000011") in {(i["ns"], i["id"]) for i in entity.identifiers}
    assert edges(extracted) == [("subclass_of", "MP:0005375")]


def test_obsolete_and_foreign_terms_give_nothing():
    assert mp.terms_schema(terms()["MP:0000012"]) is None
    assert mp.terms_schema(terms()["PATO:0000169"]) is None
    assert mp.terms_schema({}) is None


def test_resource_declares_the_ontology():
    assert mp.config.id == ResourceCv.MP
    assert mp.config.license == LicenseCV.CC_BY_4_0
    assert mp.resource.terms.mapper is mp.terms_schema
    assert Namespace.MP.value == "mp"
