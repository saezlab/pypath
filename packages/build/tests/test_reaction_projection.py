"""Reactions as typed entities, transport as a relation projected from them."""

import io
import types

from omnipath_core.source_attributes import (
    PROJECTED_FROM,
    SOURCE_REACTION,
    TRANSPORT_FROM,
    TRANSPORT_TO,
)
from omnipath_resolver import EntityResolver
from omnipath_build.projections import TRANSPORT_DATASET
from omnipath_build.silver import SilverExtractor
from omnipath_build.writer import ParquetWriter
from tables_fixture import read_tables

CYTOSOL, EXTRACELLULAR = "GO:0005829", "GO:0005576"


def build(path, datasets):
    """Write ``{dataset: [(row_id, record)]}`` as one resource; its tables."""
    writer = ParquetWriter(path)
    resolver = EntityResolver(library_dir=path / "no-library", defer_aliases=True)
    try:
        for dataset, rows in datasets.items():
            extractor = SilverExtractor("fixture", dataset)
            for row_id, record in rows:
                extractor.process_record(record, raw_payload={}, row_id=row_id, row_number=0)
            writer.append_observations(extractor, resolver)
        return read_tables(writer.close()["files"])
    finally:
        resolver.close()


def transports(tables):
    """(subject identifier, object identifier, from, to) of each projected relation."""
    ids = {e["entity_key"]: e["identifier"] for e in tables["entity"]}
    found = set()
    for relation in tables["relation"]:
        if relation["object_aspect_qualifier"] != ["transport"]:
            continue
        annotations = tables.children("relation_annotation", relation)
        values = {
            term: sorted(a["value"] for a in annotations if a["term"] == term)
            for term in (TRANSPORT_FROM, TRANSPORT_TO, PROJECTED_FROM)
        }
        assert relation["predicate"] == "affects"
        assert relation["causal_mechanism_qualifier"] == ["relocalization"]
        assert values[PROJECTED_FROM]
        found.add(
            (
                ids[relation["subject_entity_key"]],
                ids[relation["object_entity_key"]],
                tuple(values[TRANSPORT_FROM]),
                tuple(values[TRANSPORT_TO]),
            )
        )
    return found


def recon3d_row(reaction, reactants, products, genes, sbo):
    return {
        "bigg_reaction_id": reaction,
        "name": reaction,
        "sbo": sbo,
        "direction": "LEFT-TO-RIGHT",
        "reactants": reactants,
        "products": products,
        "reactant_chebi": "||".join("4167" for _ in reactants.split("||")),
        "product_chebi": "||".join("4167" for _ in products.split("||")),
        "gene_rule_clauses": [[gene] for gene in genes],
    }


def test_model_transport_projects_each_gene_onto_the_cargo(tmp_path):
    from pypath.inputs_v2 import recon3d

    rows = [
        ("reactions:0", recon3d.reactions_schema(
            recon3d_row("GLCt1", "glc__D:e:1", "glc__D:c:1", ["6513", "6514"], "SBO:0000185")
        )),
        # Same compartment on both sides: a conversion, nothing moved.
        ("reactions:1", recon3d.reactions_schema(
            recon3d_row("HEX1", "glc__D:c:1", "glc__D:c:1", ["3098"], "SBO:0000176")
        )),
    ]
    tables = build(tmp_path, {"reactions": rows})
    assert transports(tables) == {
        ("6513", "CHEBI:4167", (EXTRACELLULAR,), (CYTOSOL,)),
        ("6514", "CHEBI:4167", (EXTRACELLULAR,), (CYTOSOL,)),
    }
    projected = [e for e in tables["relation_evidence"] if e["dataset"] == TRANSPORT_DATASET]
    assert {e["row_id"] for e in projected} == {"reactions:0"}
    types_ = {
        a["value"] for a in tables["entity_annotation"] if a["term"] == "type"
    }
    assert types_ == {"SBO:0000185", "SBO:0000176"}


def test_compartments_keep_resolved_variants_of_a_reaction_apart(tmp_path):
    def reaction(compartments, member_places):
        return {
            "type": "molecular_activity",
            "identifiers": [{"type": "rhea", "value": "10000"}],
            "annotations": [
                {"term": "anatomical_context_qualifier", "value": c} for c in compartments
            ],
            "membership": [
                {
                    "member": {"type": "small_molecule", "identifiers": [{"type": "chebi", "value": "4167"}]},
                    "predicate": predicate,
                    "annotations": [{"term": "anatomical_context_qualifier", "value": place}],
                }
                for predicate, place in member_places
            ],
        }

    tables = build(
        tmp_path,
        {
            "reactions": [
                ("reactions:0", reaction([EXTRACELLULAR, CYTOSOL], [("has_input", EXTRACELLULAR), ("has_output", CYTOSOL)])),
                ("reactions:1", reaction(["GO:0005739", CYTOSOL], [("has_input", CYTOSOL), ("has_output", "GO:0005739")])),
                # Uptake and efflux across one membrane are one reaction.
                ("reactions:2", reaction([CYTOSOL, EXTRACELLULAR], [("has_input", CYTOSOL), ("has_output", EXTRACELLULAR)])),
            ]
        },
    )
    reactions = [e for e in tables["entity"] if e["entity_type"] == "molecular_activity"]
    assert len(reactions) == 2
    assert {e["reference_entity_key"] for e in reactions} == {"rhea:10000"}


def test_rhea_symport_moves_both_cargos_across_its_sides(tmp_path):
    from pypath.inputs_v2 import rhea
    from pypath.inputs_v2.parsers.rhea import _raw

    tsv = (
        "Reaction identifier\tEquation\tChEBI name\tChEBI identifier\tUniProt\n"
        "RHEA:28947\tL-arabinose(in) + H(+)(out) = L-arabinose(out) + H(+)(in)\t"
        "L-arabinopyranose;hydron\tCHEBI:17535;CHEBI:15378\tP0AE24\n"
        "RHEA:1\tATP + H2O = ADP + phosphate + H(+)\tATP;water;ADP;phosphate;hydron\t"
        "CHEBI:30616;CHEBI:15377;CHEBI:456216;CHEBI:43474;CHEBI:15378\tP00001\n"
    )
    rows = list(_raw(types.SimpleNamespace(result=io.StringIO(tsv))))
    assert [row["is_transport"] for row in rows] == [True, False]
    records = [(f"reactions:{i}", rhea.reactions_schema(row)) for i, row in enumerate(rows)]
    tables = build(tmp_path, {"reactions": records})
    assert transports(tables) == {
        ("P0AE24", "CHEBI:17535", ("in",), ("out",)),
        ("P0AE24", "CHEBI:15378", ("out",), ("in",)),
    }


def test_rhea_lists_each_equation_term_with_its_own_side():
    from pypath.inputs_v2.parsers.rhea import _participant_fields

    fields = _participant_fields(
        "sulfate(out) + ATP + H2O = sulfate(in) + ADP + phosphate + H(+)",
        "CHEBI:16189;CHEBI:30616;CHEBI:15377;CHEBI:456216;CHEBI:43474;CHEBI:15378".split(";"),
        "sulfate;ATP(4-);water;ADP(3-);hydrogenphosphate;hydron".split(";"),
    )
    chebi = fields["participant_chebi"].split("||")
    roles = fields["participant_role"].split("||")
    sides = fields["participant_compartment"].split("||")
    assert chebi[0] == chebi[3] == "CHEBI:16189"
    assert (roles[0], sides[0]) == ("reactant", "out")
    assert (roles[3], sides[3]) == ("product", "in")
    assert sides[4:] == ["", "", ""]
    assert _participant_fields("2 H2O2 = O2 + 2 H2O", ["a", "b", "c"], ["", "", ""])[
        "participant_stoichiometry"
    ] == "2||1||2"


def test_catalysis_stated_apart_pairs_with_its_own_source_reaction(tmp_path):
    """Two source reactions resolve to one entity; each catalyst moves only its own cargo."""

    def reaction(stable_id, cargo):
        return {
            "type": "molecular_activity",
            "identifiers": [{"type": "rhea", "value": "21292"}],
            "annotations": [
                {"term": "anatomical_context_qualifier", "value": EXTRACELLULAR},
                {"term": "anatomical_context_qualifier", "value": CYTOSOL},
                {"term": SOURCE_REACTION, "value": stable_id},
            ],
            "membership": [
                {
                    "member": {"type": "small_molecule", "identifiers": [{"type": "chebi", "value": cargo}]},
                    "predicate": predicate,
                    "annotations": [{"term": "anatomical_context_qualifier", "value": place}],
                }
                for predicate, place in (("has_input", EXTRACELLULAR), ("has_output", CYTOSOL))
            ],
        }

    def catalysis(stable_id, catalyst):
        return {
            "subject": {"type": "protein", "identifiers": [{"type": "uniprot", "value": catalyst}]},
            "predicate": "catalyzes",
            "object": {
                "type": "molecular_activity",
                "identifiers": [{"type": "rhea", "value": "21292"}],
                "annotations": [
                    {"term": "anatomical_context_qualifier", "value": CYTOSOL},
                    {"term": "anatomical_context_qualifier", "value": EXTRACELLULAR},
                ],
            },
            "annotations": [{"term": SOURCE_REACTION, "value": stable_id}],
        }

    tables = build(
        tmp_path,
        {
            "reactions": [("reactions:0", reaction("R-HSA-1", "4167")), ("reactions:1", reaction("R-HSA-2", "28757"))],
            "controls": [("controls:0", catalysis("R-HSA-1", "P11166")), ("controls:1", catalysis("R-HSA-2", "P22732"))],
        },
    )
    assert len([e for e in tables["entity"] if e["entity_type"] == "molecular_activity"]) == 1
    assert {(s, o) for s, o, *_ in transports(tables)} == {("P11166", "CHEBI:4167"), ("P22732", "CHEBI:28757")}


def test_pathwaycommons_transport_needs_a_catalysis_mediator():
    from pypath.inputs_v2 import pathwaycommons as pc

    def participant(identifier, kind):
        return {"id": identifier, "types": [kind], "names": [identifier], "unification_xrefs": []}

    row = {
        "interaction_type": "controls-transport-of-chemical",
        "participant_a": participant("SLC2A1", "ProteinReference"),
        "participant_b": participant("chebi:4167", "SmallMoleculeReference"),
        "data_sources": ["Reactome"],
        "pubmed": [],
        "mediator_ids": "http://bioregistry.io/reactome:R-HSA-70153;reactome:Catalysis5850",
    }

    def qualifiers(row):
        return {
            str(a.term.name if hasattr(a.term, "name") else a.term): str(getattr(a.value, "text", a.value))
            for a in pc.map_interaction(row).annotations
            if "qualifier" in str(getattr(a.term, "name", a.term))
        }

    assert qualifiers(row) == {
        "object_aspect_qualifier": "transport",
        "causal_mechanism_qualifier": "relocalization",
    }
    regulated = {**row, "mediator_ids": "reactome:Control12"}
    assert qualifiers(regulated) == {"object aspect qualifier": "localization"}


def test_error_page_in_place_of_data_is_recognised(tmp_path):
    from pypath.share.downloads import _is_error_page

    page = tmp_path / "substrates.tsv"
    page.write_text("<!DOCTYPE HTML PUBLIC>\n<html><title>500 Internal Server Error</title>")
    data = tmp_path / "data.tsv"
    data.write_text("2.A.1.1.51\tCHEBI:5584;hydron\n")
    html = tmp_path / "families.html"
    html.write_text("<html></html>")
    assert _is_error_page(page, "tsv")
    assert not _is_error_page(data, "tsv")
    assert not _is_error_page(html, "html")
