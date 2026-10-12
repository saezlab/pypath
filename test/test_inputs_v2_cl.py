"""Cell Ontology terms: names, synonyms, definitions and hierarchy."""

from omnipath_build.silver import SilverExtractor
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import cl
from pypath.inputs_v2.parsers.obo import parse_obo_text

# Stanzas copied from cl-basic.obo (release 2026-06-08), with xrefs trimmed; the UBERON
# stanza is what a merged edition (cl.obo) adds.
CL_OBO = """format-version: 1.2
data-version: cl/releases/2026-06-08/cl-basic.owl
ontology: cl/cl-basic

[Term]
id: CL:0002518
name: kidney epithelial cell
def: "An epithelial cell of the kidney." [GOC:tfm, KUPO:SJ]
subset: cellxgene_subset
xref: KUPO:0001019
is_a: CL:0000066 ! epithelial cell
is_a: CL:1000497 {is_inferred="true"} ! kidney cell
property_value: terms:contributor https://orcid.org/0000-0003-1980-3228
creation_date: 2011-02-08T10:46:34Z

[Term]
id: CL:0000312
name: keratinocyte
alt_id: CL:0002185
def: "An epithelial cell of stratified squamous tissues that produces keratin proteins." [GOC:tfm]
xref: BTO:0000667
is_a: CL:0000237 {is_inferred="true"} ! keratinizing barrier epithelial cell

[Term]
id: CL:0000066
name: epithelial cell
synonym: "epitheliocyte" EXACT []

[Term]
id: CL:0000005
name: neural crest derived fibroblast
synonym: "fibroblast neural crest derived" EXACT []
is_a: CL:0000057 ! fibroblast
relationship: RO:0002202 CL:0000008 ! develops from migratory cranial neural crest cell

[Term]
id: CL:0002453
name: oligodendrocyte precursor cell
synonym: "NG2 cell" RELATED [PMID:14501223]
synonym: "O2-A" RELATED OMO:0003000 [PMID:2182078]
is_a: CL:0000123 ! neuron associated cell (sensu Vertebrata)
relationship: BFO:0000050 CL:0000125 ! part of glial cell

[Term]
id: CL:0000002
name: obsolete immortal cell line cell
comment: Obsoleted in July 2013 and replaced by the CLO 'immortal cell line cell' class.
is_obsolete: true
consider: CL:0000000

[Term]
id: CP:0000000
name: obsolete CP:0000000
comment: ID was moved to the CL namespace. Use CL:0017500.
is_obsolete: true
replaced_by: CL:0017500

[Term]
id: UBERON:0002113
name: kidney
"""


def terms():
    return {row["id"]: row for row in parse_obo_text(CL_OBO)}


def extract(row):
    record = cl.terms_schema(row)
    extractor = SilverExtractor("cl", "terms")
    extractor.process_record(record, raw_payload={}, row_id="1", row_number=1)
    return record, extractor


def edges(extractor):
    return sorted(
        (relation.predicate, extractor.entities[relation.object_entity_key].identifier)
        for relation in extractor.relations
    )


def test_term_has_name_definition_and_parents():
    record, extracted = extract(terms()["CL:0002518"])
    (entity,) = (e for e in extracted.entities.values() if e.identifier == "CL:0002518")
    assert (entity.entity_type, entity.namespace) == ("ontology_class", "cl")
    assert entity.label == "kidney epithelial cell"
    annotations = {format_term(a.term): a.value for a in record.annotations}
    assert annotations["description"] == "An epithelial cell of the kidney."
    assert annotations["xref"] == "KUPO:0001019"
    assert edges(extracted) == [("subclass_of", "CL:0000066"), ("subclass_of", "CL:1000497")]
    assert {e.namespace for e in extracted.entities.values()} == {"cl"}


def test_alt_id_and_synonyms_are_identifiers_of_the_term():
    _, extracted = extract(terms()["CL:0000312"])
    (entity,) = (e for e in extracted.entities.values() if e.identifier == "CL:0000312")
    assert entity.label == "keratinocyte"
    assert ("cl", "CL:0002185") in {(i["ns"], i["id"]) for i in entity.identifiers}
    _, extracted = extract(terms()["CL:0002453"])
    (entity,) = (e for e in extracted.entities.values() if e.identifier == "CL:0002453")
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert {("synonym", "NG2 cell"), ("synonym", "O2-A")} <= identifiers


def test_develops_from_and_part_of_are_edges():
    _, extracted = extract(terms()["CL:0000005"])
    assert edges(extracted) == [("develops_from", "CL:0000008"), ("subclass_of", "CL:0000057")]
    _, extracted = extract(terms()["CL:0002453"])
    assert edges(extracted) == [("part_of", "CL:0000125"), ("subclass_of", "CL:0000123")]


def test_obsolete_and_foreign_terms_give_nothing():
    assert cl.terms_schema(terms()["CL:0000002"]) is None
    assert cl.terms_schema(terms()["CP:0000000"]) is None
    assert cl.terms_schema(terms()["UBERON:0002113"]) is None
    assert cl.terms_schema({}) is None


def test_resource_declares_the_ontology():
    assert cl.config.id == ResourceCv.CL
    assert cl.config.license == LicenseCV.CC_BY_4_0
    assert cl.resource.terms.mapper is cl.terms_schema
