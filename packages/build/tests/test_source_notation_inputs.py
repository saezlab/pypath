"""Source notation becomes readable annotations, publications and structure."""

import io
from types import SimpleNamespace

from omnipath_build.silver import SilverExtractor
from pypath.inputs_v2 import brenda, intact, recon3d
from pypath.inputs_v2.parsers import brenda as brenda_parser
from pypath.inputs_v2.parsers import recon3d as recon3d_parser


def opener(text):
    return SimpleNamespace(result={"data": io.StringIO(text)})


def annotations(record):
    result = {}
    for annotation in record.annotations or []:
        result.setdefault(annotation.term, []).append(annotation.value)
    return result


BRENDA = """ID\t1.1.1.100
PR\t#1# Salmonella enterica <19,48>
PR\t#5# Escherichia coli {P0AEK2; source: UniProt} (#5# isoform FabG <4>) <1,4>
EN\t#1,5# E233K (#1# temperature-sensitive mutant enzyme does not allow growth
\tat 42°C <19>; #5# reductase activity is thermolabile <4>; #5# no growth defect <48>) <4,19,48>
EN\t#1,5# more (#5# 38-BamHI site <4>) <4,19>
EN\t#5# F434A (#5# the mutation leads to 53% ee (S) after 1 h))
PM\t#5# glycoprotein (#5# 3 N-glycosylation sites <19>) <19>
RF\t<4> Citation {Pubmed:111}
RF\t<19> Citation {Pubmed:222}
RF\t<48> Citation {Pubmed:}
///
"""


def test_brenda_commentary_is_split_per_protein_with_its_own_references():
    rows = list(brenda_parser.iter_molecular_forms(opener(BRENDA)))
    records = [brenda.molecular_forms_schema(row) for row in rows]
    by_row = [
        (row["protein_record_id"], row["descriptor"], annotations(record))
        for row, record in zip(rows, records)
    ]

    assert by_row[0] == (
        "1",
        "E233K",
        {
            "brenda:molecular_observation": [
                "E233K: temperature-sensitive mutant enzyme does not allow growth at 42°C"
            ],
            "in_taxon_label": ["Salmonella enterica"],
            "publications": ["PMID:222"],
        },
    )
    assert by_row[1] == (
        "5",
        "E233K",
        {
            "brenda:molecular_observation": [
                "E233K: reductase activity is thermolabile",
                "E233K: no growth defect",
            ],
            "in_taxon_label": ["Escherichia coli"],
            "brenda:protein_note": ["isoform FabG"],
            # <48> has no PubMed ID in the reference list.
            "publications": ["PMID:111"],
        },
    )
    assert records[1].molecular_form["variants"][0]["position"] == 233
    # "more" without commentary for protein 1 says nothing about it.
    assert [row["protein_record_id"] for row in rows if row["descriptor"] == "more"] == ["5"]
    assert annotations(records[2])["brenda:molecular_observation"] == ["38-BamHI site"]
    # A doubled closing parenthesis still separates the commentary.
    assert annotations(records[3])["brenda:molecular_observation"] == [
        "F434A: the mutation leads to 53% ee (S) after 1 h"
    ]
    modification = records[4].molecular_form["modifications"][0]
    assert (modification["term"], modification["description"]) == (
        "glycoprotein",
        "3 N-glycosylation sites",
    )
    values = [
        value for record in records for values in annotations(record).values() for value in values
    ]
    assert not any("#" in value or "<" in value for value in values)


def test_recon3d_reactions_use_go_compartments_biopax_directions_and_integral_coefficients():
    data = {
        "metabolites": [],
        "genes": [{"id": "1_AT1"}, {"id": "1_AT2"}, {"id": "2_AT1"}],
        "reactions": [
            {
                "id": "R1",
                "metabolites": {"atp_c": -1.0, "h_m": -2.0, "adp_c": 1.0, "x_i": 0.5},
                "lower_bound": -1000,
                "upper_bound": 0,
                "gene_reaction_rule": "(1_AT1 and 2_AT1) or (1_AT2 and 2_AT1)",
            },
            {"id": "R2", "metabolites": {}, "lower_bound": -1, "upper_bound": 1},
        ],
    }
    first, second = recon3d_parser._parse_reactions(data)
    assert first["reactants"] == "atp:c:1||h:m:2"
    assert first["products"] == "adp:c:1||x:i:0.5"
    assert (first["direction"], second["direction"]) == ("RIGHT-TO-LEFT", "REVERSIBLE")

    extracted = SilverExtractor("recon3d", "reactions")
    extracted.process_record(
        recon3d.reactions_schema(first), raw_payload={}, row_id="1", row_number=1
    )
    participants = [
        (
            relation.predicate,
            *(
                [a["value"] for a in relation.annotations if a["term"] == term]
                for term in (
                    "biopax:cellularLocation",
                    "stoichiometry",
                    "biopax:conversionDirection",
                )
            ),
        )
        for relation in extracted.relations
        if relation.predicate in {"has_input", "has_output"}
    ]
    assert participants == [
        ("has_input", ["cytosol"], ["1"], ["RIGHT-TO-LEFT"]),
        ("has_input", ["mitochondrion"], ["2"], ["RIGHT-TO-LEFT"]),
        ("has_output", ["cytosol"], ["1"], ["RIGHT-TO-LEFT"]),
        ("has_output", ["mitochondrial intermembrane space"], ["0.5"], ["RIGHT-TO-LEFT"]),
    ]
    # The gene rule is structure: one AND group of genes, no serialized rule.
    rules = [r for r in extracted.relations if r.predicate == "associated_with"]
    assert len(rules) == 1 and not rules[0].annotations
    members = [r for r in extracted.relations if r.predicate == "has_member"]
    assert sorted(extracted.entities[r.object_entity_key].identifier for r in members) == ["1", "2"]
    (group,) = recon3d_parser._parse_enzyme_complexes(data)
    assert not recon3d.enzyme_complexes_schema(group).annotations
    for gene in recon3d_parser._parse_genes(data):
        assert set(annotations(recon3d.genes_schema(gene))) == {"in_taxon"}


def mitab_row(**values):
    row = {
        "#ID(s) interactor A": "uniprotkb:P46109",
        "ID(s) interactor B": "uniprotkb:Q13191",
        "Interaction detection method(s)": 'psi-mi:"MI:0019"(coimmunoprecipitation)',
        "Publication Identifier(s)": "pubmed:10022120|mint:MINT-6731034",
        "Taxid interactor A": "taxid:9606(human)|taxid:9606(Homo sapiens)",
        "Taxid interactor B": "taxid:9606(human)",
        "Interaction type(s)": 'psi-mi:"MI:0914"(association)',
        "Source database(s)": 'psi-mi:"MI:0471"(MINT)',
        "Interaction identifier(s)": "intact:EBI-6932994",
        "Biological role(s) interactor A": 'psi-mi:"MI:0499"(unspecified role)',
        "Biological role(s) interactor B": 'psi-mi:"MI:0501"(enzyme)',
        "Experimental role(s) interactor A": 'psi-mi:"MI:0496"(bait)',
        "Experimental role(s) interactor B": 'psi-mi:"MI:0498"(prey)',
        "Type(s) interactor A": 'psi-mi:"MI:0326"(protein)',
        "Type(s) interactor B": 'psi-mi:"MI:0326"(protein)',
        "Identification method participant A": 'psi-mi:"MI:0396"(predetermined participant)',
        "Identification method participant B": "-",
        "Stoichiometry(s) interactor A": "0",
        "Stoichiometry(s) interactor B": "2",
        "Feature(s) interactor A": (
            'tandem tag:c-c|his tag:?-?(MINT-1)|flag tag:n-n("MI:0518")|mutation:429-429,430-430'
        ),
        "Feature(s) interactor B": (
            "binding-associated region:702-715(MINT-8028122)"
            "|sufficient binding region:634-913(IPR000242)"
            "|mutation with no effect:315-315"
        ),
    }
    return {**row, **values}


def test_intact_columns_have_distinct_terms_and_readable_features():
    relation = intact.interactions_schema(mitab_row())

    assert annotations(relation) == {
        "original_predicate": ["MI:0914"],
        "has_evidence_of_type": ["MI:0019"],
        "supporting_data_source": ["MI:0471"],
        "publications": ["PMID:10022120"],
    }
    subject, object_ = annotations(relation.subject), annotations(relation.object)
    # The mutation is a variant of the molecular form, not a feature annotation.
    assert len(relation.subject.molecular_form["variants"]) == 2
    assert subject == {
        "psi_mi:participant_feature": [
            "tandem tag (C-terminal)",
            "his tag",
            "flag tag (N-terminal)",
        ],
        "in_taxon": ["NCBITaxon:9606"],
        "psi_mi:experimental_role": ["MI:0496"],
        "psi_mi:participant_identification_method": ["MI:0396"],
    }
    # Binding regions and the mutation with no effect are in the molecular form.
    assert object_ == {
        "in_taxon": ["NCBITaxon:9606"],
        "psi_mi:biological_role": ["MI:0501"],
        "psi_mi:experimental_role": ["MI:0498"],
        "stoichiometry": ["2"],
    }
    form = relation.object.molecular_form
    assert [(r["type"], r["position"], r["end_position"]) for r in form["regions"]] == [
        ("binding-associated region", 702, 715),
        ("sufficient binding region", 634, 913),
    ]
    assert [(v["position"], v["description"]) for v in form["variants"]] == [
        (315, "mutation with no effect:315-315")
    ]
    assert "has_topic" not in subject and "description" not in subject
