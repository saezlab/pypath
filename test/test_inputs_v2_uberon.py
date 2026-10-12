"""Uberon terms: names, synonyms, definitions and hierarchy."""

from omnipath_build.silver import SilverExtractor
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import uberon
from pypath.inputs_v2.parsers.obo import parse_obo_text

# Stanzas copied from uberon-basic.obo (release 2026-10-01), with xrefs, subsets and long
# lines trimmed; the CL stanza is what a merged edition (uberon.obo) adds.
UBERON_OBO = """format-version: 1.2
data-version: uberon/releases/2026-10-01/uberon-basic.owl
ontology: uberon/uberon-basic

[Term]
id: UBERON:0002113
name: kidney
def: "A paired organ of the urinary tract that produces urine and maintains bodily fluid homeostasis." [Wikipedia:Kidney]
subset: major_organ
synonym: "reniculate kidney" NARROW SENSU [Wikipedia:Reniculate_kidney]
xref: BTO:0000671
is_a: UBERON:0000489 ! cavitated compound organ
is_a: UBERON:0005172 ! abdomen element
relationship: dc-contributor https://orcid.org/0000-0001-7655-4833
relationship: develops_from UBERON:0003918 {evidence="definitional"} ! kidney mesenchyme
relationship: fma_set_term FMA:264815
relationship: mutually_spatially_disjoint_with UBERON:0000160 ! intestine
relationship: part_of UBERON:0011143 {source="FMA"} ! upper urinary tract
relationship: present_in_taxon NCBITaxon:7742
property_value: has_relational_adjective "renal" xsd:string

[Term]
id: UBERON:0000044
name: dorsal root ganglion
alt_id: UBERON:0026602
is_a: UBERON:0002025 ! ganglion

[Term]
id: UBERON:0000000
name: obsolete processual entity
comment: "This term has been deprecated in favour of BFO:0000015 ('process' in COB)"
is_obsolete: true

[Term]
id: CL:0000066
name: epithelial cell
"""


def terms():
    return {row["id"]: row for row in parse_obo_text(UBERON_OBO)}


def extract(row):
    record = uberon.terms_schema(row)
    extractor = SilverExtractor("uberon", "terms")
    extractor.process_record(record, raw_payload={}, row_id="1", row_number=1)
    return record, extractor


def edges(extractor):
    return sorted(
        (relation.predicate, extractor.entities[relation.object_entity_key].identifier)
        for relation in extractor.relations
    )


def test_term_has_name_synonym_definition_and_hierarchy():
    record, extracted = extract(terms()["UBERON:0002113"])
    (entity,) = (e for e in extracted.entities.values() if e.identifier == "UBERON:0002113")
    assert (entity.entity_type, entity.namespace) == ("ontology_class", "uberon")
    assert entity.label == "kidney"
    assert ("synonym", "reniculate kidney") in {(i["ns"], i["id"]) for i in entity.identifiers}
    annotations = {format_term(a.term): a.value for a in record.annotations}
    assert annotations["description"].startswith("A paired organ of the urinary tract")
    assert annotations["xref"] == "BTO:0000671"
    # Taxon constraints, disjointness, contributors and FMA links are not hierarchy.
    assert edges(extracted) == [
        ("develops_from", "UBERON:0003918"),
        ("part_of", "UBERON:0011143"),
        ("subclass_of", "UBERON:0000489"),
        ("subclass_of", "UBERON:0005172"),
    ]
    assert {e.namespace for e in extracted.entities.values()} == {"uberon"}


def test_alt_id_is_an_identifier_of_the_term():
    _, extracted = extract(terms()["UBERON:0000044"])
    (entity,) = (e for e in extracted.entities.values() if e.identifier == "UBERON:0000044")
    assert ("uberon", "UBERON:0026602") in {(i["ns"], i["id"]) for i in entity.identifiers}


def test_obsolete_and_foreign_terms_give_nothing():
    assert uberon.terms_schema(terms()["UBERON:0000000"]) is None
    assert uberon.terms_schema(terms()["CL:0000066"]) is None
    assert uberon.terms_schema({}) is None


def test_resource_declares_the_ontology():
    assert uberon.config.id == ResourceCv.UBERON
    assert uberon.config.license == LicenseCV.CC_BY_3_0
    assert uberon.resource.terms.mapper is uberon.terms_schema
