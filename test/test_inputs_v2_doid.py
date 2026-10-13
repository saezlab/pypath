"""Human Disease Ontology terms: names, synonyms, definitions and hierarchy."""

from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import doid
from pypath.inputs_v2.parsers.obo import parse_obo_text

# Stanzas copied from doid.obo (data-version releases/2026-09-30/doid.obo); one foreign stanza added.
DOID_OBO = r"""format-version: 1.2
data-version: releases/2026-09-30/doid.obo
ontology: doid

[Term]
id: DOID:3702
name: cervical adenocarcinoma
def: "A cervix carcinoma that derives_from epithelial cells of glandular origin." [url:http\://en.wikipedia.org/wiki/Adenocarcinoma]
subset: DO_cancer_slim
subset: NCIthesaurus
synonym: "adenocarcinoma cervix uteri" EXACT []
synonym: "adenocarcinoma of cervix" EXACT []
synonym: "adenocarcinoma of the uterine Cervix" EXACT []
xref: NCI:C4029
xref: SNOMEDCT_US_2025_09_01:254887002
xref: UMLS_CUI:C0279672
is_a: DOID:2893 ! cervix carcinoma
is_a: DOID:299 ! adenocarcinoma

[Term]
id: DOID:8552
name: chronic myeloid leukemia
alt_id: DOID:8551
alt_id: DOID:8606
def: "A myeloid leukemia that is characterized by over production of white blood cells." [url:http\://www.cancer.gov/dictionary?CdrID=46755]
comment: OMIM mapping confirmed by DO. [SN].
subset: DO_cancer_slim
subset: DO_RAD_slim
subset: DO_rare_slim
subset: NCIthesaurus
synonym: "chronic granulocytic leukaemia" EXACT []
synonym: "chronic granulocytic leukemia" EXACT []
synonym: "chronic myelogenous leukaemia" EXACT []
synonym: "chronic myelogenous leukemia" EXACT []
synonym: "chronic myeloid leukaemia" EXACT []
synonym: "CML" EXACT OMO:0003012 []
synonym: "CML - chronic Myelogenous Leukemia" EXACT []
synonym: "Myeloid Leukemia, chronic" EXACT []
xref: EFO:0000339
xref: GARD:6105
xref: ICD9CM:205.1
xref: ICDO:9863/3
xref: KEGG:05220
xref: MESH:D015464
xref: MIM:608232
xref: NCI:C9128
xref: ORDO:521
xref: SNOMEDCT_US_2025_09_01:92818009
is_a: DOID:8692 ! myeloid leukemia

[Term]
id: DOID:0001816
name: angiosarcoma
alt_id: DOID:267
alt_id: DOID:4508
def: "A vascular cancer that derives_from the cells that line the walls of blood vessels or lymphatic vessels." [url:http\://en.wikipedia.org/wiki/Hemangiosarcoma, url:https\://en.wikipedia.org/wiki/Angiosarcoma, url:https\://ncit.nci.nih.gov/ncitbrowser/ConceptReport.jsp?dictionary=NCI_Thesaurus&ns=ncit&code=C3088, url:https\://www.ncbi.nlm.nih.gov/pubmed/23327728]
subset: DO_cancer_slim
subset: NCIthesaurus
synonym: "hemangiosarcoma" EXACT []
xref: ICDO:9120/3
xref: MESH:D006394
xref: NCI:C3088
xref: NCI:C9275
xref: SNOMEDCT_US_2025_09_01:39000009
xref: UMLS_CUI:C0018923
xref: UMLS_CUI:C0854893
is_a: DOID:175 ! vascular cancer

[Term]
id: DOID:0050001
name: obsolete Actinomadura madurae infectious disease
subset: gram-positive_bacterial_infectious_disease
is_obsolete: true

[Term]
id: UBERON:0002113
name: kidney
"""


def terms():
    return {row["id"]: row for row in parse_obo_text(DOID_OBO)}


def extract(row):
    record = doid.terms_schema(row)
    extractor = SilverExtractor("doid", "terms")
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
    record, extracted = extract(terms()["DOID:3702"])
    entity = term_entity(extracted, "DOID:3702")
    # the identifier keeps the prefix, exactly as BioGRID cites it
    assert (entity.entity_type, entity.namespace) == ("ontology_class", "doid")
    assert entity.label == "cervical adenocarcinoma"
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert {
        ("doid", "DOID:3702"),
        ("name", "cervical adenocarcinoma"),
        ("synonym", "adenocarcinoma cervix uteri"),
        ("synonym", "adenocarcinoma of cervix"),
    } <= identifiers
    annotations = {(format_term(a.term), a.value) for a in record.annotations}
    assert (
        "description",
        "A cervix carcinoma that derives_from epithelial cells of glandular origin.",
    ) in annotations
    assert ("xref", "NCI:C4029") in annotations
    assert edges(extracted) == [("subclass_of", "DOID:2893"), ("subclass_of", "DOID:299")]
    assert {e.namespace for e in extracted.entities.values()} == {"doid"}


def test_alt_ids_are_identifiers_of_the_term():
    _, extracted = extract(terms()["DOID:8552"])
    entity = term_entity(extracted, "DOID:8552")
    assert entity.label == "chronic myeloid leukemia"
    identifiers = {(i["ns"], i["id"]) for i in entity.identifiers}
    assert {("doid", "DOID:8551"), ("doid", "DOID:8606"), ("synonym", "CML")} <= identifiers
    assert edges(extracted) == [("subclass_of", "DOID:8692")]


def test_term_with_a_single_parent():
    _, extracted = extract(terms()["DOID:0001816"])
    assert term_entity(extracted, "DOID:0001816").label == "angiosarcoma"
    assert edges(extracted) == [("subclass_of", "DOID:175")]


def test_obsolete_and_foreign_terms_give_nothing():
    assert doid.terms_schema(terms()["DOID:0050001"]) is None
    assert doid.terms_schema(terms()["UBERON:0002113"]) is None
    assert doid.terms_schema({}) is None


def test_resource_declares_the_ontology():
    assert doid.config.id == ResourceCv.DOID
    assert doid.config.license == LicenseCV.CC0_1_0
    assert doid.resource.terms.mapper is doid.terms_schema
    assert Namespace.DOID.value == "doid"
