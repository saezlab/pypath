"""TEMPORARY: the old Python export against the relations step file (removed after)."""
import pyarrow as pa
import pyarrow.parquet as pq

from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.release_database import serve_newest
from test_client import canonical, service  # noqa: F401


def test_old_and_new_export(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    cases = [
        {"entity_ids": ["TP53"]}, {"entity_ids": ["TP53"], "sources": ["alpha"]},
        {"entity_ids": ["TP53"], "gene_mode": "source"}, {"entity_ids": ["TP53"], "gene_mode": "both"},
        {"entity_ids": ["caffeine", "O'Brien"], "gene_mode": "both"}, {"entity_ids": ["gene:entrez:7157"]},
        {"entity_ids": ["uniprot:P04637"]}, {"reference_entity_keys": ["entrez:7157"]},
        {"entity_ids": ["TP53"], "reference_entity_keys": ["entrez:7157"]},
        {"entity_ids": ["TP53"], "scope_entity_ids": ["caffeine"]},
        {"predicates": ["biolink:interacts_with"]},
        {"interaction_types": ["protein"], "partner_types": ["small_molecule"]},
        {"taxonomy_ids": ["9606"], "relation_categories": ["interaction"]}, {"sources": ["beta"]},
        {"signs": [1]}, {"sign": -1}, {"is_directed": True}, {"min_evidence_count": 2}, {},
        {"protein_entity_keys": ["x"]}, {"entity_ids": ["TP53"], "protein_entity_keys": ["x"]},
    ]
    with api.query_scope():
        rows = api._db.execute(
            f"SELECT DISTINCT subject_molecular_form.protein_entity_key FROM {api._served().alias}.main.relation_evidence "
            "WHERE subject_molecular_form.protein_entity_key IS NOT NULL"
        ).fetchall()
    cases += [{"protein_entity_keys": [r[0]]} for r in rows] + [{"entity_ids": ["TP53"], "protein_entity_keys": [r[0]]} for r in rows]
    for filters in cases:
        with api.query_scope():
            old = pq.read_table(pa.BufferReader(api._old_export_slice(filters, None)[0]))
            new = pq.read_table(pa.BufferReader(api.export_slice(filters, None)[0]))
        assert new.schema.names == old.schema.names, filters
        assert canonical(new) == canonical(old), filters
        assert new.num_rows == old.num_rows, filters
    print(len(cases), "cases")
