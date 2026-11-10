"""Relation exports (the step file relations) and the relation queries' parameters."""

from __future__ import annotations

from typing import Any


# Relations one export returns; whole resources are downloads of their published tables.
RELATION_EXPORT_LIMIT = 100_000
PARQUET = "application/vnd.apache.parquet"


def relation_params(filters):
    """A request's relation filters (normalized) as the parameters of the relation
    queries."""
    from omnipath_api.queries.filters import normalize_filters

    f = normalize_filters(filters)
    lists = (
        "reference_entity_keys", "sources", "predicates", "relation_categories", "taxonomy_ids",
        "interaction_types", "entity_types", "partner_types", "signs", "object_aspect_qualifier",
        "object_direction_qualifier", "causal_mechanism_qualifier", "protein_entity_keys",
        "transcript_entity_keys", "isoform_identifiers", "scope_entity_ids",
    )
    params = {name: (list(f.get(name) or []) or None) for name in lists}
    tokens = list(f.get("entity_ids") or []) + list(f.get("entity_pks") or [])
    params.update(
        entities=tokens or None,
        gene_mode=f.get("gene_mode") or "any",
        scope_endpoint_mode=f.get("scope_endpoint_mode") or "any",
        molecular_endpoint_mode=f.get("molecular_endpoint_mode") or "any",
        is_directed=f.get("is_directed"),
        sign=f.get("sign"),
        min_evidence_count=f.get("min_evidence_count"),
    )
    # Qualifier values as stored (Biolink enumerations); an unknown one is an error.
    from omnipath_core.biolink import annotation_value

    for term in ("object_aspect_qualifier", "object_direction_qualifier", "causal_mechanism_qualifier"):
        if params[term] is not None:
            params[term] = [annotation_value(term, value) for value in params[term]]
    if params["signs"] is not None:
        params["signs"] = [int(s) for s in params["signs"]]
    if params["taxonomy_ids"] is not None:
        params["taxonomy_ids"] = [str(t) for t in params["taxonomy_ids"]]
    return params


class ExportQueries:
    """Relation exports as Parquet."""

    def export_slice(
        self,
        filters: dict[str, Any] | None = None,
        resources: list[str] | None = None,
        limit: int = RELATION_EXPORT_LIMIT,
    ) -> tuple[bytes, str]:
        """A filtered relation slice as Parquet, with its annotations and evidence: the
        query ``relations``, which the clients run too."""
        import io

        import pyarrow.parquet as pq

        params = {**relation_params(filters), "max_relations": int(limit)}
        table = self.run_query("relations", params, resources, output="table")
        sink = io.BytesIO()
        pq.write_table(table, sink, compression="zstd")
        return sink.getvalue(), PARQUET
