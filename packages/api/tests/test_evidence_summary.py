"""Evidence summaries of many relations at once: sites and PubMed ids per source."""

import duckdb
import pytest
from fastapi.testclient import TestClient

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_api.settings import Settings
from omnipath_core.fixtures import write_resource

KINASE, PLN, PRODUCT, OTHER = [str(i) * 64 for i in range(1, 5)]


def publications(*values):
    return [dict(term="publications", value=value, scope="relation") for value in values]


def form(*modifications, protein=None):
    return dict(protein_entity_key=protein, modifications=list(modifications))


def phospho(residue, position):
    return dict(term="MOD:00696", residue=residue, position=position)


def site(side="object", **fields):
    return dict(
        dict(
            side=side, term=None, residue=None, position=None, end_position=None, description=None
        ),
        **fields,
    )


S16, T17 = (site(term="MOD:00696", residue=r, position=p) for r, p in [("S", 16), ("T", 17)])


def entities():
    return [
        dict(
            entity_key=key,
            namespace="entrez",
            identifier=identifier,
            entity_type="protein",
            label=label,
            taxon="9606",
            reference_entity_key=f"entrez:{identifier}",
            gene_reference_keys=[f"entrez:{identifier}"],
        )
        for key, identifier, label in [
            (KINASE, "5566", "PRKACA"),
            (PLN, "5350", "PLN"),
            (PRODUCT, "5350", "PLN product"),
            (OTHER, "7", "OTHER"),
        ]
    ]


def relation(key, subject, evidence, annotations=(), source="a"):
    return dict(
        relation_key=key,
        subject_entity_key=subject,
        object_entity_key=PLN,
        subject_reference_entity_key="entrez:5566" if subject == KINASE else "entrez:7",
        object_reference_entity_key="entrez:5350",
        subject_type="protein",
        object_type="protein",
        subject_label="PRKACA" if subject == KINASE else "OTHER",
        object_label="PLN",
        predicate="biolink:affects",
        category="interaction",
        taxon="9606",
        sources=[source],
        annotations=list(annotations),
        evidence=evidence,
    )


@pytest.fixture
def engine(tmp_path):
    # Resource a states the phosphorylation with sites from two sources; b states the
    # same relation (same key) without sites; a also has a relation without evidence.
    write_resource(
        tmp_path / "resources" / "a" / "1",
        entities(),
        [
            relation(
                "phospho",
                KINASE,
                [
                    dict(
                        source="psp",
                        row_id="1",
                        object_molecular_form=form(phospho("S", 16), protein=PRODUCT),
                        annotations=publications("PMID:1", "PMID:2"),
                    ),
                    dict(
                        source="psp",
                        row_id="2",
                        object_molecular_form=form(phospho("S", 16), phospho("T", 17)),
                        annotations=publications("PMID:2", "DOI:10.1/x", "PMID:3; PMID:6"),
                    ),
                    dict(
                        source="signor",
                        row_id="3",
                        object_molecular_form=form(dict(description="PLN_S16")),
                        annotations=publications("4"),
                    ),
                ],
                annotations=publications("PMID:9"),
            ),
            relation("bare", OTHER, []),
        ],
    )
    write_resource(
        tmp_path / "resources" / "b" / "1",
        entities(),
        [
            relation(
                "phospho",
                KINASE,
                [
                    dict(
                        source="b",
                        row_id="1",
                        subject_molecular_form=form(phospho("S", 338)),
                        annotations=publications("PMID:5"),
                    )
                ],
                source="b",
            )
        ],
    )
    queries = ParquetQueries(tmp_path)
    yield queries
    queries.close()


def test_summary_per_relation_and_source(engine):
    result = engine.get_relation_evidence_summaries(["phospho", "bare", "unknown", "phospho"])
    assert result["evidence"] == {
        "phospho": [
            dict(source="psp", sites=[S16, T17], pmids=["1", "2", "3", "6"]),
            dict(source="signor", sites=[site(description="PLN_S16")], pmids=["4"]),
            dict(
                source="b",
                sites=[site("subject", term="MOD:00696", residue="S", position=338)],
                pmids=["5"],
            ),
        ],
        "bare": [],
        "unknown": [],
    }
    assert engine.get_relation_evidence_summaries([])["evidence"] == {}


def test_summary_matches_the_evidence_endpoint(engine):
    """The summary holds what a client read from GET /relations/{pk}/evidence."""
    summaries = engine.get_relation_evidence_summaries(["phospho", "bare"])["evidence"]
    for pk, summary in summaries.items():
        expected = {}
        for item in engine.get_relation_evidence(pk)["evidence"]:
            sites, pmids = expected.setdefault(item["source"], (set(), set()))
            for side in ("subject", "object"):
                for mod in (item[f"{side}MolecularForm"] or {}).get("modifications") or []:
                    sites.add(
                        (side, mod.get("residue"), mod.get("position"), mod.get("description"))
                    )
            for record in item["recordAttributes"]:
                if record["term"] == "publications":
                    for part in record["value"].split(";"):
                        if part.strip().removeprefix("PMID:").isdigit():
                            pmids.add(part.strip().removeprefix("PMID:"))
        found = {
            item["source"]: (
                {(s["side"], s["residue"], s["position"], s["description"]) for s in item["sites"]},
                set(item["pmids"]),
            )
            for item in summary
        }
        assert found == expected


def test_records_without_annotations_have_the_relations(engine, tmp_path):
    path = tmp_path / "resources" / "a" / "1" / "relation_evidence.parquet"
    with duckdb.connect() as db:
        db.execute(
            f"COPY (SELECT * REPLACE (CASE WHEN source = 'signor' THEN NULL ELSE annotations END"
            f" AS annotations) FROM read_parquet('{path}')) TO '{path}.new' (FORMAT PARQUET)"
        )
    path.with_suffix(".parquet.new").replace(path)
    engine.reload_resources()
    signor = engine.get_relation_evidence_summaries(["phospho"])["evidence"]["phospho"][1]
    assert signor == dict(source="signor", sites=[site(description="PLN_S16")], pmids=["9"])


def test_molecular_form_filters_keep_matching_occurrences(engine):
    result = engine.get_relation_evidence_summaries(
        ["phospho"], filters={"protein_entity_keys": [PRODUCT]}
    )
    assert result["evidence"] == {"phospho": [dict(source="psp", sites=[S16], pmids=["1", "2"])]}


def test_resources_narrow_the_summary(engine):
    result = engine.get_relation_evidence_summaries(["phospho"], resources=["b"])
    assert [item["source"] for item in result["evidence"]["phospho"]] == ["b"]


def test_too_many_relations(engine):
    with pytest.raises(ValueError, match="At most 1000"):
        engine.get_relation_evidence_summaries([str(i) for i in range(1001)])


def test_http(engine):
    client = TestClient(create_app(queries=engine, settings=Settings(admin_secret="test")))
    response = client.post("/relations/evidence", json={"relation_pks": ["phospho", "bare"]})
    assert response.status_code == 200
    body = response.json()
    assert body["evidence"]["bare"] == []
    assert body["evidence"]["phospho"][0] == dict(
        source="psp", sites=[S16, T17], pmids=["1", "2", "3", "6"]
    )
    assert "X-Query-Duration-Ms" in response.headers
    # The single-relation endpoints stay as they were.
    assert client.get("/relations/phospho/evidence").status_code == 200
    assert client.get("/relations/phospho").json()["relation"]["relationPk"] == "phospho"
    too_many = client.post(
        "/relations/evidence", json={"relation_pks": [str(i) for i in range(1001)]}
    )
    assert too_many.status_code == 422
