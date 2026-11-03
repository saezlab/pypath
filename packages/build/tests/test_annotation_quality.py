"""Source annotations become statements, structured forms or one vocabulary."""

import io

from omnipath_build.silver import SilverExtractor
from omnipath_core.silver_schema import format_term
from pypath.inputs_v2 import (
    chemont,
    guidetopharma,
    imm1415,
    intact,
    kegg,
    metatlas,
    mondo,
    psi_mi,
    reactome,
    rhea,
)
from pypath.inputs_v2.parsers import kegg as kegg_parser
from pypath.inputs_v2.parsers.obo import parse_obo_text


def extract(record, source="test"):
    extractor = SilverExtractor(source, "test")
    extractor.process_record(record, raw_payload={}, row_id="1", row_number=1)
    return extractor


def annotations(record):
    result = {}
    for annotation in record.annotations or []:
        result.setdefault(format_term(annotation.term), []).append(annotation.value)
    return result


def term(module, text):
    (row,) = parse_obo_text(text)
    return module.terms_schema(row)


MONDO_TERM = """[Term]
id: MONDO:0009861
name: phenylketonuria
def: "An autosomal recessive inborn error of phenylalanine metabolism." [Orphanet:716]
comment: See genetic heterogeneity of OMIM 261600.
xref: OMIM:261600 {source="MONDO:equivalentTo"}
is_a: MONDO:0017750 ! disorder of phenylalanine metabolism
relationship: BFO:0000050 MONDO:0045024 ! part of cancer or benign tumor
relationship: RO:0002007 UBERON:0000483 ! bounding layer of epithelium
relationship: RO:0002213 GO:0006366 ! positively regulates transcription by RNA polymerase II
relationship: RO:0002573 HP:0000007 ! has modifier Autosomal recessive inheritance
relationship: RO:0004001 http://identifiers.org/hgnc/7582 {source="Orphanet:544602"} ! has material basis in gain of function germline mutation in MYL1
"""


def test_mondo_relationships_between_terms_are_ontology_statements():
    record = term(mondo, MONDO_TERM)
    assert annotations(record) == {
        "description": ["An autosomal recessive inborn error of phenylalanine metabolism."],
        "rdfs:comment": ["See genetic heterogeneity of OMIM 261600."],
        "xref": ["OMIM:261600"],
        # 'has modifier' is directional; Biolink only has the symmetric related_to.
        "RO:0002573": ["HP:0000007"],
        # A gene record, not an ontology term: normalised, kept as an attribute.
        "RO:0004001": ["HGNC:7582"],
    }
    extracted = extract(record)
    assert {r.statement_kind for r in extracted.relations} == {"ontology"}
    edges = {
        (r.predicate, extracted.entities[r.object_entity_key].identifier): [
            (a["term"], a["value"]) for a in r.annotations
        ]
        for r in extracted.relations
    }
    assert edges == {
        ("subclass_of", "MONDO:0017750"): [],
        ("part_of", "MONDO:0045024"): [],
        # Biolink's part_of is broader than 'bounding layer of'.
        ("part_of", "UBERON:0000483"): [("original_predicate", "RO:0002007")],
        ("regulates", "GO:0006366"): [("original_predicate", "RO:0002213")],
    }


def test_mondo_gene_relation_keeps_its_ro_predicate():
    relation = mondo.gene_disease_association_to_entity(
        {"mondo_id": "MONDO:1", "mondo_name": "x", "relation": "RO:0004001", "hgnc_id": "7582"}
    )
    assert annotations(relation) == {"original_predicate": ["RO:0004001"]}


def test_psi_mi_xrefs_are_unescaped_and_property_values_dropped():
    record = term(
        psi_mi,
        """[Term]
id: MI:0471
name: MINT
def: "The Molecular INTeraction database." [PMID:14681455]
comment: Reference not index in medline.
xref: GO:GO\\:0044093
xref: url:http\\://mint.bio.uniroma2.it/mint/
xref: search-url: "http://mint.bio.uniroma2.it/mint/search/interaction.do?ac=${ac}"
xref: id-validation-regexp: "MINT\\-[0-9]+"
""",
    )
    assert annotations(record) == {
        "description": ["The Molecular INTeraction database."],
        "rdfs:comment": ["Reference not index in medline."],
        "xref": ["GO:0044093"],
        "url": ["http://mint.bio.uniroma2.it/mint/"],
    }


def test_chemont_book_and_article_xrefs_are_publications():
    record = term(
        chemont,
        """[Term]
id: CHEMONTID:0000002
name: Organoheterocyclic compounds
def: "Compounds containing a ring with least one carbon atom and one non-carbon atom." []
xref: ISBN:0967855098 "IUPAC. Compendium of Chemical Terminology, 2nd ed. (the 'Gold Book')."
xref: PMID:23731671
xref: Wikipedia:Heterocyclic_compound
""",
    )
    assert annotations(record)["publications"] == ["ISBN:0967855098", "PMID:23731671"]
    assert annotations(record)["xref"] == ["Wikipedia:Heterocyclic_compound"]


def participants(extracted, *terms):
    return [
        (
            relation.predicate,
            *([a["value"] for a in relation.annotations if a["term"] == t] for t in terms),
        )
        for relation in extracted.relations
        if relation.predicate in {"has_input", "has_output"}
    ]


def test_human_gem_compartments_and_directions_follow_recon3d():
    record = metatlas.reactions_schema(
        {
            "human_gem_reaction_id": "MAR04358",
            "reactants": "MAM01371:c:1||MAM02039:i:4",
            "products": "MAM01371:m:1",
            "direction": "left_to_right",
        }
    )
    assert participants(
        extract(record), "biopax:cellularLocation", "biopax:conversionDirection"
    ) == [
        ("has_input", ["cytosol"], ["LEFT-TO-RIGHT"]),
        ("has_input", ["mitochondrial intermembrane space"], ["LEFT-TO-RIGHT"]),
        ("has_output", ["mitochondrion"], ["LEFT-TO-RIGHT"]),
    ]


def test_rhea_in_and_out_are_transport_sides():
    record = rhea.reactions_schema(
        {
            "rhea_id": "29523",
            "participant_chebi": "15378||29101||15378||29101",
            "participant_role": "reactant||reactant||product||product",
            "participant_compartment": "out||out||in||in",
        }
    )
    assert participants(extract(record), "rhea:transport_side", "biopax:cellularLocation") == [
        ("has_input", ["out"], []),
        ("has_input", ["out"], []),
        ("has_output", ["in"], []),
        ("has_output", ["in"], []),
    ]


def test_reactome_coefficients_are_integers_where_integral():
    record = reactome.reactions_schema(
        {
            "entity_type": "reaction",
            "reactome_stable_id": "R-HSA-71541",
            "participant_entity_type": "chemical||chemical||chemical",
            "participant_role": "reactant||product||product",
            "participant_chebi": "15377||15378||30616",
            "participant_stoichiometry": "2.0||0.5||n",
            "direction": "LEFT-TO-RIGHT",
        }
    )
    assert [v for _, v in participants(extract(record), "stoichiometry")] == [["2"], ["0.5"], ["n"]]


def test_kegg_reaction_topics_and_uninformative_direction():
    direction, reactants, products = kegg_parser._parse_kegg_equation(
        "C00031 + C00002 <=> 2 C00008"
    )
    assert direction == "" and products == [{"kegg_id": "C00008", "stoichiometry": "2"}]
    assert kegg_parser._parse_kegg_equation("C00031 => C00008")[0] == "LEFT-TO-RIGHT"
    record = kegg.reactions_schema(
        {
            "reaction_id": "R00010",
            "ko_ids": "ko:K01194;ko:K22934",
            "rclass_ids": "rc:RC00049",
            "ec_numbers": "3.2.1.28",
        }
    )
    assert annotations(record)["has_topic"] == [
        "EC:3.2.1.28",
        "KEGG.ORTHOLOGY:K01194",
        "KEGG.ORTHOLOGY:K22934",
        "KEGG.RCLASS:RC00049",
    ]


def test_kegg_pathway_member_reactions_are_split_topics():
    record = kegg.pathways_schema(
        {
            "pathway_term_accession": "hsa00010",
            "protein_member_kegg_ids": "hsa:3101||hsa:2538",
            "protein_member_uniprot_ids": "P52789||P35575",
            "protein_member_reaction_ids": "R00299;R01786||",
        }
    )
    topics = [annotations(member).get("has_topic", []) for member in record.membership]
    assert topics == [["KEGG.REACTION:R00299", "KEGG.REACTION:R01786"], []]


def test_kegg_organism_pathway_terms_have_no_taxon_xref():
    handles = {
        "organism_list_pathways": {"hsa": io.StringIO("path:hsa00010\tGlycolysis - Homo sapiens\n")}
    }
    (row,) = kegg_parser._build_pathway_term_records(handles)
    record = kegg.pathway_ontology_terms_schema(row)
    assert annotations(record) == {"xref": ["KEGG_PATHWAY:hsa00010"]}


def test_intact_mutations_residues_and_binding_regions_are_structured():
    row = {
        "#ID(s) interactor A": "uniprotkb:P00533",
        "ID(s) interactor B": "uniprotkb:P62993",
        "Interaction type(s)": 'psi-mi:"MI:0915"(physical association)',
        "Type(s) interactor A": 'psi-mi:"MI:0326"(protein)',
        "Type(s) interactor B": 'psi-mi:"MI:0326"(protein)',
        "Feature(s) interactor A": (
            "O4'-phospho-L-tyrosine:1092-1092|mutation decreasing interaction strength:1068-1068"
            "|mutation with complex effect:700-700(376004)|variant:436-436|ha tag:c-c"
        ),
        "Feature(s) interactor B": 'sufficient binding region:55-152("MI:0117")|necessary binding region:?-?',
    }
    relation = intact.interactions_schema(row)
    subject, object_ = relation.subject.molecular_form, relation.object.molecular_form
    assert [(m["term"], m["position"]) for m in subject["modifications"]] == [
        ("O4'-phospho-L-tyrosine", 1092)
    ]
    assert [v["position"] for v in subject["variants"]] == [1068, 700, 436]
    assert [(r["type"], r["position"], r["end_position"]) for r in object_["regions"]] == [
        ("sufficient binding region", 55, 152),
        ("necessary binding region", None, None),
    ]
    # Only the tag is left as a participant feature.
    assert annotations(relation.subject)["psi_mi:participant_feature"] == ["ha tag (C-terminal)"]
    assert "psi_mi:participant_feature" not in annotations(relation.object)


def test_guidetopharma_single_precision_affinities_are_curated_decimals():
    row = {
        "Affinity Units": "pIC50",
        "Affinity Median": "6.46999979019165",
        "Original Affinity Units": "IC50",
        "Original Affinity Median nm": "340.0",
        "Original Affinity Low nm": "0.39800000190734863",
        "Original Affinity Relation": "=",
    }
    assert [m.quantity.has_numeric_value for m in guidetopharma._affinity_measurements(row)] == [
        6.47
    ]
    assert [
        m.quantity.has_numeric_value for m in guidetopharma._original_affinity_measurements(row)
    ] == [0.398, 340.0]


def test_imm1415_has_no_constant_sbo_topic():
    record = imm1415.schema(
        {"bigg.metabolite": "glc__D", "sbo": "SBO:0000247", "name": "D-Glucose"}
    )
    assert not record.annotations
