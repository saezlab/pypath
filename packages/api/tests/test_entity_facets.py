"""Scalar facets retain usable counts when entities are searched by stored CURIEs."""

import pytest
from fastapi.testclient import TestClient

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_core.keys import entity_key
from omnipath_core.fixtures import write_resource


@pytest.fixture
def reference_facets_engine(tmp_path):
    for source, records in [
        (
            "a",
            [
                ("gene", "entrez", "11052", "9606", "entrez:11052"),
                ("protein", "entrez", "11052", "9606", "entrez:11052"),
                ("protein", "uniprot", "Q16630", "9606", "entrez:11052"),
            ],
        ),
        (
            "b",
            [
                ("rna_product", "entrez", "11052", "9606", "entrez:11052"),
                ("transcript", "enst", "ENST12345", "10090", "entrez:11052"),
                ("protein", "uniprot", "P99999", "9606", None),
            ],
        ),
    ]:
        entities = [
            dict(
                entity_key=entity_key(kind, ns, identifier),
                entity_type=kind,
                namespace=ns,
                identifier=identifier,
                taxon=taxon,
                label="CPSF6" if ref else "Unrelated protein",
                reference_entity_key=ref,
                gene_reference_keys=[ref] if ref else [],
                has_hierarchy=False,
                parent_count=0,
                child_count=0,
                identifiers=[],
                annotations=[],
                evidence=[],
            )
            for kind, ns, identifier, taxon, ref in records
        ]
        write_resource(tmp_path / "resources" / source / "1", entities)
    return ParquetQueries(tmp_path)


def counts(rows):
    return {(row["facetName"], row["facetValue"]): row["scopedCount"] for row in rows}


def test_gene_reference_facets_include_all_matching_source_types(reference_facets_engine):
    engine = reference_facets_engine
    results = engine.search_entities_api(query="entrez:11052")["entities"]
    assert len(results) == 5
    facets = counts(engine.get_scoped_entity_facets({"query": "entrez:11052"}))
    assert facets == {
        ("entity_type", "gene"): 1,
        ("entity_type", "protein"): 2,
        ("entity_type", "rna_product"): 1,
        ("entity_type", "transcript"): 1,
        ("taxonomy_id", "9606"): 4,
        ("taxonomy_id", "10090"): 1,
        ("source", "a"): 3,
        ("source", "b"): 2,
    }
    with TestClient(create_app(queries=engine)) as client:
        response = client.post("/entities/scoped-facets", json={"query": "entrez:11052"})
        assert response.status_code == 200
        assert counts(response.json()) == facets


@pytest.mark.parametrize("query", ["uniprot:Q16630", "uniprot|Q16630", "uniprot:P99999"])
def test_native_curie_facets_match_entity_search(reference_facets_engine, query):
    engine = reference_facets_engine
    result = engine.search_entities_api(query=query)["entities"]
    assert len(result) == 1
    assert result[0]["canonicalIdentifier"] == query.replace("|", ":").split(":")[-1]
    facets = counts(engine.get_scoped_entity_facets({"query": query}))
    source = "b" if "P99999" in query else "a"
    assert facets == {
        ("entity_type", "protein"): 1,
        ("taxonomy_id", "9606"): 1,
        ("source", "a"): int(source == "a"),
        ("source", "b"): int(source == "b"),
    }


def test_reference_facets_keep_source_cross_filtering_and_other_filters(reference_facets_engine):
    engine = reference_facets_engine
    selected = counts(engine.get_scoped_entity_facets({"query": "entrez:11052", "sources": ["a"]}))
    assert selected == {
        ("entity_type", "gene"): 1,
        ("entity_type", "protein"): 2,
        ("taxonomy_id", "9606"): 3,
        ("source", "a"): 3,
        ("source", "b"): 2,
    }
    proteins = counts(
        engine.get_scoped_entity_facets(
            {
                "query": "entrez:11052",
                "sources": ["a"],
                "entityTypes": ["protein"],
            }
        )
    )
    assert proteins == {
        ("entity_type", "protein"): 2,
        ("taxonomy_id", "9606"): 2,
        ("source", "a"): 2,
        ("source", "b"): 0,
    }
    mouse = counts(
        engine.get_scoped_entity_facets(
            {
                "query": "entrez:11052",
                "ncbi_tax_id": ["10090"],
            }
        )
    )
    assert mouse == {
        ("entity_type", "transcript"): 1,
        ("taxonomy_id", "10090"): 1,
        ("source", "a"): 0,
        ("source", "b"): 1,
    }
