"""C. elegans Phenotype Ontology terms: names, synonyms, definitions and hierarchy."""

from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import wbphenotype
from pypath.inputs_v2.parsers.obo import parse_obo_text

# Stanzas copied from wbphenotype.obo (data-version releases/2026-04-15); one foreign stanza added.
WBPHENOTYPE_OBO = """format-version: 1.2
data-version: releases/2026-04-15
ontology: wbphenotype

[Term]
id: WBPhenotype:0000050
name: embryonic lethal
def: "Animals die during embryonic development. In C. elegans, often assayed as refractile eggs that fail to hatch; when applied to large-scale RNAi screens in C. elegans, more than 10% of embryos die." [WB:cab, WB:WBPaper00000179, WB:WBPaper00004403, WB:WBPaper00004540, WB:WBPaper00004651, WB:WBPaper00004769, WB:WBPaper00005654, WB:WBPaper00024497, WB:WBPaper00024925, WB:WBPaper00025054]
subset: phenotype_slim_wb
synonym: "egg lethal" RELATED []
synonym: "Emb" BROAD []
synonym: "embryonic death" RELATED []
synonym: "inviable zygote" RELATED []
is_a: WBPhenotype:0000062 ! lethal
is_a: WBPhenotype:0000749 ! embryonic development variant

[Term]
id: WBPhenotype:0000033
name: developmental timing variant
alt_id: WBPhenotype:0000437
def: "The timing of specific developmental events in some tissues is altered relative to the timing of events in other tissues." [pmid:6494891, WB:cab]
comment: Possible XP GO:0040034.
subset: phenotype_slim_wb
synonym: "developmental timing abnormal" NARROW []
synonym: "heterochronic defect" EXACT []
is_a: WBPhenotype:0000531 ! organism development variant

[Term]
id: WBPhenotype:0000995
def: "OBSOLETE." []
comment: WBPhenotype:0000995 an alt_id for "sodium chloride chemotaxis defective" (WBPhenotype:0001084) (WS193 OBO and WS200 site) AND has term name "pos body wall contraction defect" (dump_phenotype_ace.pl)
is_obsolete: true

[Term]
id: GO:0008150
name: biological_process
"""


def terms():
    return {row["id"]: row for row in parse_obo_text(WBPHENOTYPE_OBO)}


def extract(row):
    record = wbphenotype.terms_schema(row)
    extractor = SilverExtractor("wbphenotype", "terms")
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
    record, extracted = extract(terms()["WBPhenotype:0000050"])
    entity = term_entity(extracted, "WBPhenotype:0000050")
    # the mixed-case prefix is kept, exactly as BioGRID cites it
    assert (entity.entity_type, entity.namespace) == ("ontology_class", "wbphenotype")
    assert entity.label == "embryonic lethal"
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert {
        ("wbphenotype", "WBPhenotype:0000050"),
        ("name", "embryonic lethal"),
        ("synonym", "egg lethal"),
        ("synonym", "Emb"),
    } <= identifiers
    description = {format_term(a.term): a.value for a in record.annotations}["description"]
    assert description.startswith("Animals die during embryonic development.")
    assert edges(extracted) == [
        ("subclass_of", "WBPhenotype:0000062"),
        ("subclass_of", "WBPhenotype:0000749"),
    ]
    assert {e.namespace for e in extracted.entities.values()} == {"wbphenotype"}


def test_alt_id_is_an_identifier_of_the_term():
    _, extracted = extract(terms()["WBPhenotype:0000033"])
    entity = term_entity(extracted, "WBPhenotype:0000033")
    assert entity.label == "developmental timing variant"
    assert ("wbphenotype", "WBPhenotype:0000437") in {(i["ns"], i["id"]) for i in entity.identifiers}
    assert edges(extracted) == [("subclass_of", "WBPhenotype:0000531")]


def test_obsolete_unnamed_and_foreign_terms_give_nothing():
    assert wbphenotype.terms_schema(terms()["WBPhenotype:0000995"]) is None  # obsolete, no name
    assert wbphenotype.terms_schema(terms()["GO:0008150"]) is None
    assert wbphenotype.terms_schema({}) is None


def test_prefix_is_matched_as_written():
    row = {"id": "wbphenotype:0000050", "name": "embryonic lethal"}
    assert wbphenotype.terms_schema(row) is None


def test_resource_declares_the_ontology():
    assert wbphenotype.config.id == ResourceCv.WBPHENOTYPE
    assert wbphenotype.config.license == LicenseCV.CC_BY_4_0
    assert wbphenotype.resource.terms.mapper is wbphenotype.terms_schema
    assert Namespace.WBPHENOTYPE.value == "wbphenotype"
