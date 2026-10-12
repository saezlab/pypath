"""Ascomycete Phenotype Ontology terms: names, synonyms, definitions and hierarchy."""

from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import apo
from pypath.inputs_v2.parsers.obo import parse_obo_text

# Stanzas copied from apo.obo (data-version releases/2026-09-30); one foreign stanza added.
APO_OBO = """format-version: 1.2
data-version: releases/2026-09-30
ontology: apo

[Term]
id: APO:0000063
name: colony size
namespace: observable
alt_id: YPO:0000063
def: "The size of colonies of cells growing on a solid substrate." [SGD:mcc]
subset: AspGD
subset: CGD
subset: CryptoGD
subset: SGD
is_a: APO:0000060 ! colony appearance

[Term]
id: APO:0000011
name: null mutant
namespace: mutant_type
alt_id: YPO:0000011
def: "The mutation abolishes completely the function of the gene product." [SGD:RSN]
subset: SGD
synonym: "loss of function" RELATED [SGD:RSN]
synonym: "null" EXACT []
is_a: APO:0000006 ! mutant_type

[Term]
id: APO:0000245
name: decreased rate
namespace: qualifier
alt_id: curators
def: "Observed feature is decreased with respect to time compared to that in the wild type (slower)." [SGD:RSN]
subset: SGD
synonym: "slower" EXACT [SGD:RSN]
is_a: APO:0000002 ! abnormal

[Term]
id: APO:0000244
name: PCA
namespace: experiment_type
alt_id: YPO:0000244
subset: SGD
xref: BioGRID:curators
is_a: APO:0000187 ! physical interactions

[Term]
id: APO:0000017
name: observable
namespace: observable
alt_id: YPO:0000017
def: "Features of ascomycete fungal cells, cultures, or colonies that can be detected, observed, measured, or monitored." [SGD:curators]
subset: SGD

[Term]
id: APO:0000064
name: liquid culture appearance
namespace: observable
comment: This term was made obsolete because there was overlap with another term unspecified. The term unspecified will be used in cases where a more appropriate descriptor cannot be identified.
subset: AspGD
subset: SGD
is_obsolete: true

[Term]
id: GO:0008150
name: biological_process
"""


def terms():
    return {row["id"]: row for row in parse_obo_text(APO_OBO)}


def extract(row):
    record = apo.terms_schema(row)
    extractor = SilverExtractor("apo", "terms")
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
    record, extracted = extract(terms()["APO:0000063"])
    entity = term_entity(extracted, "APO:0000063")
    # the identifier keeps the prefix, exactly as BioGRID cites it
    assert (entity.entity_type, entity.namespace) == ("ontology_class", "apo")
    assert entity.label == "colony size"
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert ("apo", "APO:0000063") in identifiers
    assert ("name", "colony size") in identifiers
    assert {format_term(a.term): a.value for a in record.annotations} == {
        "description": "The size of colonies of cells growing on a solid substrate."
    }
    assert edges(extracted) == [("subclass_of", "APO:0000060")]
    assert {e.namespace for e in extracted.entities.values()} == {"apo"}


def test_synonyms_and_every_term_namespace_are_kept():
    _, extracted = extract(terms()["APO:0000011"])  # a mutant type
    entity = term_entity(extracted, "APO:0000011")
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert {("synonym", "loss of function"), ("synonym", "null")} <= identifiers
    assert edges(extracted) == [("subclass_of", "APO:0000006")]
    for accession, name in (("APO:0000245", "decreased rate"), ("APO:0000244", "PCA")):
        _, extracted = extract(terms()[accession])  # a qualifier, an experiment type
        assert term_entity(extracted, accession).label == name


def test_former_yeast_phenotype_ids_are_alt_ids_and_curator_names_are_not():
    _, extracted = extract(terms()["APO:0000063"])
    identifiers = {(i["ns"], i["id"]) for i in term_entity(extracted, "APO:0000063").identifiers}
    assert ("ypo", "YPO:0000063") in identifiers
    _, extracted = extract(terms()["APO:0000245"])
    identifiers = {(i["ns"], i["id"]) for i in term_entity(extracted, "APO:0000245").identifiers}
    assert all(value != "curators" for _, value in identifiers)
    assert ("synonym", "slower") in identifiers


def test_term_without_parent_is_a_root():
    _, extracted = extract(terms()["APO:0000017"])
    assert edges(extracted) == []
    assert term_entity(extracted, "APO:0000017").label == "observable"


def test_obsolete_and_foreign_terms_give_nothing():
    assert apo.terms_schema(terms()["APO:0000064"]) is None
    assert apo.terms_schema(terms()["GO:0008150"]) is None
    assert apo.terms_schema({}) is None


def test_resource_declares_the_ontology():
    assert apo.config.id == ResourceCv.APO
    assert apo.config.license == LicenseCV.CC_BY_4_0
    assert apo.resource.terms.mapper is apo.terms_schema
    assert Namespace.APO.value == "apo"
