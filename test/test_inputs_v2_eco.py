"""Evidence and Conclusion Ontology terms: names, synonyms, definitions and hierarchy."""

from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import eco
from pypath.inputs_v2.parsers.obo import parse_obo_text

# Stanzas copied from eco.obo (data-version releases/2026-07-10); eco.obo holds ECO's own terms
# only, so the foreign stanza at the end is added: a merged edition would carry such terms.
ECO_OBO = r"""format-version: 1.2
data-version: releases/2026-07-10
ontology: eco

[Term]
id: ECO:0000178
name: in vivo assay evidence
def: "A type of direct assay evidence resulting from a sample of microorganisms or living tissue grown in a medium." [ECO:RCT, url:https\://www.cancer.gov/publications/dictionaries/cancer-terms?cdrid=46352]
is_a: ECO:0000002 ! direct assay evidence

[Term]
id: ECO:0000181
name: in vitro assay evidence
alt_id: ECO:0005013
def: "A type of direct assay evidence resulting from a process performed outside the living organism, where the components of an organism have been isolated from their usual biological context, to permit a more detailed or more convenient analysis in an artificial setting." [ECO:RCT, url:https\://www.cancer.gov/publications/dictionaries/cancer-terms?cdrid=45733]
is_a: ECO:0000002 ! direct assay evidence

[Term]
id: ECO:0000068
name: yeast 2-hybrid evidence
def: "A type of bait-prey hybrid interaction evidence that is based on detection of protein-protein interaction by activation of a yeast reporter gene after a bait protein fused to a DNA-binding domain (which has been transfected into a yeast cell) is used to screen a cDNA library of clones fused to an activation domain." [PMID:12734586]
synonym: "yeast two-hybrid assay" RELATED []
is_a: ECO:0000025 ! bait-prey hybrid interaction evidence

[Term]
id: ECO:0005502
name: ex vivo assay evidence
def: "A type of direct assay evidence derived by studying an organ, tissues, or cells taken from an organism and studied in an external environment with the minimum alteration of natural conditions." [ECO:SN, PMID:16305312, url:http\://www.ncbi.nlm.nih.gov/pmc/articles/PMC2305722/pdf/nihms42807.pdf]
is_a: ECO:0000002 ! direct assay evidence
created_by: snadendla
creation_date: 2015-06-26T11:46:15Z

[Term]
id: ECO:0000352
name: evidence used in manual assertion
def: "A type of evidence that is used in an manual assertion." [ECO:MCC]
is_a: ECO:0000000 {is_inferred="true"} ! evidence
intersection_of: ECO:0000000 ! evidence
intersection_of: used_in ECO:0000218 ! manual assertion
relationship: used_in ECO:0000218 {is_inferred="true"} ! manual assertion
created_by: mchibucos

[Term]
id: ECO:0006163
name: nuclear magnetic resonance spectroscopy evidence
def: "A type of nuclear magnetic resonance evidence used for quantification of metabolites or for the determination of chemical structure or composition." [ECO:SN, PMID:16428685, PMID:23036848]
synonym: "NMR Spectroscopy" EXACT []
is_a: ECO:0001105 ! nuclear magnetic resonance evidence
relationship: OBI:0000312 OBI:0000623
property_value: IAO:0000234 "SIB:PG" xsd:string
created_by: snadendla

[Term]
id: ECO:0000039
name: obsolete protein assay evidence
alt_id: ECO:0005022
def: "OBSOLETE. A type of direct assay evidence resulting from determining the presence, abundance, structure, function, or activity of proteins." [PMID:18429326]
is_obsolete: true

[Term]
id: OBI:0000623
name: a term of another ontology
"""


def terms():
    return {row["id"]: row for row in parse_obo_text(ECO_OBO)}


def extract(row):
    record = eco.terms_schema(row)
    extractor = SilverExtractor("eco", "terms")
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
    record, extracted = extract(terms()["ECO:0000178"])
    entity = term_entity(extracted, "ECO:0000178")
    # the identifier keeps the prefix, exactly as the inputs cite it
    assert (entity.entity_type, entity.namespace) == ("ontology_class", "eco")
    assert entity.label == "in vivo assay evidence"
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert {("eco", "ECO:0000178"), ("name", "in vivo assay evidence")} <= identifiers
    description = {format_term(a.term): a.value for a in record.annotations}["description"]
    assert description.startswith("A type of direct assay evidence resulting from a sample of microorganisms")
    assert edges(extracted) == [("subclass_of", "ECO:0000002")]
    assert {e.namespace for e in extracted.entities.values()} == {"eco"}


def test_the_accessions_the_inputs_cite_are_named_terms():
    names = {
        "ECO:0000178": "in vivo assay evidence",
        "ECO:0000181": "in vitro assay evidence",
        "ECO:0000068": "yeast 2-hybrid evidence",
        "ECO:0005502": "ex vivo assay evidence",
    }
    for accession, name in names.items():
        _, extracted = extract(terms()[accession])
        entity = term_entity(extracted, accession)
        assert (entity.entity_type, entity.namespace, entity.label) == ("ontology_class", "eco", name)


def test_synonym_and_alt_id():
    _, extracted = extract(terms()["ECO:0000068"])
    identifiers = {(i["ns"], i["id"]) for i in term_entity(extracted, "ECO:0000068").identifiers}
    assert ("synonym", "yeast two-hybrid assay") in identifiers
    assert edges(extracted) == [("subclass_of", "ECO:0000025")]
    _, extracted = extract(terms()["ECO:0000181"])
    identifiers = {(i["ns"], i["id"]) for i in term_entity(extracted, "ECO:0000181").identifiers}
    assert ("eco", "ECO:0005013") in identifiers


def test_used_in_and_inferred_parent_qualifier():
    # `used_in` is no Biolink predicate and has no CURIE: it gives no edge; the parent keeps its id
    # although the file marks it {is_inferred="true"}
    _, extracted = extract(terms()["ECO:0000352"])
    assert edges(extracted) == [("subclass_of", "ECO:0000000")]


def test_obi_relationship_is_an_attribute_not_an_edge_to_a_foreign_term():
    record, extracted = extract(terms()["ECO:0006163"])
    assert edges(extracted) == [("subclass_of", "ECO:0001105")]
    assert ("OBI:0000312", "OBI:0000623") in {(format_term(a.term), a.value) for a in record.annotations}
    assert {e.namespace for e in extracted.entities.values()} == {"eco"}


def test_obsolete_and_foreign_terms_give_nothing():
    assert eco.terms_schema(terms()["ECO:0000039"]) is None
    assert eco.terms_schema(terms()["OBI:0000623"]) is None
    assert eco.terms_schema({}) is None


def test_resource_declares_the_ontology():
    assert eco.config.id == ResourceCv.ECO
    assert eco.config.license == LicenseCV.CC0_1_0
    assert eco.resource.terms.mapper is eco.terms_schema
    assert Namespace.ECO.value == "eco"
