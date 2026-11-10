"""Export query service, including domain SQL and result shaping."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any


from omnipath_api.queries.store.connection import sql_literal


logger = logging.getLogger(__name__)


def copy_parquet_command(select_sql: str, target_path: Path) -> str:
    return (
        f"COPY ({select_sql}) TO {sql_literal(str(target_path))} (FORMAT PARQUET, COMPRESSION ZSTD)"
    )


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
    """Export queries over the release database and the shaping contract."""

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

    def _old_export_slice(
        self,
        filters: dict[str, Any] | None = None,
        resources: list[str] | None = None,
        limit: int = RELATION_EXPORT_LIMIT,
    ) -> tuple[bytes, str]:
        """A filtered relation slice as Parquet, with its annotations and evidence."""
        from omnipath_api.queries.molecular import form_match_sql, has_form_filters
        from omnipath_api.queries.filters import normalize_filters
        from omnipath_core.measurements import QUANTITY_STRUCT

        filters = normalize_filters(filters)
        scope = resources or filters["sources"] or None
        # A plain entity selection from its range of the sorted endpoint table.
        fast = self._entity_relation_sql(filters, resources)
        selection, params = fast if fast is not None else self._relation_selection(filters, resources)
        # Two phases, like the relation pages: the matching relations first, then only
        # their annotations and evidence, read per resource by sorted relation id (a join
        # against every resource's evidence table read millions of nested rows).
        pairs = self._db.execute(
            f"SELECT resource, relation_id FROM ({selection}) LIMIT {int(limit)}", params
        ).fetchall()
        names = [f"quantity_{field.name}" for field in QUANTITY_STRUCT]
        quantity = (
            f"CASE WHEN {' AND '.join(f'{n} IS NULL' for n in names)} THEN NULL ELSE struct_pack("
            + ", ".join(f"{field.name} := quantity_{field.name}" for field in QUANTITY_STRUCT)
            + ") END"
        )
        form, form_params = "TRUE", []
        if has_form_filters(filters):
            form, form_params = form_match_sql(filters)
        # Each part is read on its own into a temporary table (read inside one statement,
        # DuckDB scanned whole tables for them), and the evidence once.
        staged = {
            "__export_relations": self._lookup_sql("relation", "relation_id", pairs),
            "__export_annotations": self._lookup_sql("relation_annotation", "relation_id", pairs),
            "__export_evidence": self._lookup_sql(
                "relation_evidence", "relation_id", pairs, where=form, params=form_params
            ),
        }
        try:
            for name, (part, part_params) in staged.items():
                self._db.execute(f"CREATE OR REPLACE TEMP TABLE {name} AS {part}", part_params)
            # The products the evidence names, then only those entities, by key: a join
            # of every entity row read the whole entity table.
            product_keys = [
                row[0]
                for row in self._db.execute(
                    """SELECT DISTINCT unnest([subject_molecular_form.protein_entity_key,
                        subject_molecular_form.transcript_entity_key,
                        object_molecular_form.protein_entity_key,
                        object_molecular_form.transcript_entity_key]) AS k FROM __export_evidence"""
                ).fetchall()
                if row[0]
            ]
            self._db.execute(
                f"CREATE OR REPLACE TEMP TABLE __export_products AS SELECT * FROM {self._table('entity', scope)} "
                f"WHERE entity_key IN ({', '.join(sql_literal(k) for k in product_keys) or 'NULL'})"
            )
            # Rows as the published tables nest them: each relation with its annotations
            # and evidence (only the matching occurrences under a molecular form filter),
            # and the product entities its evidence names.
            sql = f"""WITH selected AS (SELECT * FROM __export_relations),
            annotations AS (
                SELECT resource, relation_id, list(struct_pack(term := term, value := value,
                    quantity := {quantity}, source := source, dataset := dataset, scope := scope)
                    ORDER BY ordinal) AS annotations
                FROM __export_annotations GROUP BY ALL
            ), evidence AS (
                SELECT resource, relation_id, list(struct_pack(source := source, dataset := dataset,
                    row_id := row_id, upstream_id := upstream_id, annotations := annotations,
                    subject_molecular_form := subject_molecular_form,
                    object_molecular_form := object_molecular_form) ORDER BY ordinal) AS evidence
                FROM __export_evidence GROUP BY ALL
            ), rows AS (
                SELECT s.* EXCLUDE (evidence_count),
                    coalesce(a.annotations, []) AS annotations, coalesce(e.evidence, []) AS evidence,
                    {"len(coalesce(e.evidence, []))" if form_params else "s.evidence_count"} AS evidence_count
                FROM selected s LEFT JOIN annotations a USING (resource, relation_id)
                LEFT JOIN evidence e USING (resource, relation_id)
            ), product_keys AS (
                SELECT DISTINCT resource, relation_id, unnest([f.subject_molecular_form.protein_entity_key,
                    f.subject_molecular_form.transcript_entity_key,
                    f.object_molecular_form.protein_entity_key,
                    f.object_molecular_form.transcript_entity_key]) AS entity_key
                FROM rows, unnest(rows.evidence) AS occurrences(f)
            ), products AS (
                -- each row's own products: one list shared by every row repeated it per row
                SELECT k.resource, k.relation_id,
                    list(e ORDER BY e.entity_key, e.resource) AS referenced_product_records
                FROM product_keys k JOIN __export_products e USING (entity_key)
                GROUP BY ALL
            ) SELECT rows.* EXCLUDE (resource, relation_id),
                coalesce(products.referenced_product_records, []) AS referenced_product_records
            FROM rows LEFT JOIN products USING (resource, relation_id)"""
            return self._write_parquet(sql, []), PARQUET
        finally:
            for name in (*staged, "__export_products"):
                self._db.execute(f"DROP TABLE IF EXISTS {name}")

    def _write_parquet(self, sql: str, params: list[Any]) -> bytes:
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            self._db.execute(copy_parquet_command(sql, tmp_path), params)
            return tmp_path.read_bytes()
        finally:
            tmp_path.unlink(missing_ok=True)
