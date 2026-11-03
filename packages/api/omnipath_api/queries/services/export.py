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


class ExportQueries:
    """Export queries over ParquetQueries' storage and shaping contract."""

    def export_slice(
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
        selection, params = self._relation_selection(filters, resources)
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
        relations, relation_params = self._lookup_sql("relation", "relation_id", pairs)
        annotations, annotation_params = self._lookup_sql(
            "relation_annotation", "relation_id", pairs
        )
        evidence, evidence_params = self._lookup_sql(
            "relation_evidence", "relation_id", pairs, where=form, params=form_params
        )
        # Rows as the published tables nest them: each relation with its annotations and
        # evidence (only the matching occurrences under a molecular form filter), and
        # the product entities its evidence names.
        sql = f"""WITH selected AS MATERIALIZED ({relations}),
        annotations AS (
            SELECT resource, relation_id, list(struct_pack(term := term, value := value,
                quantity := {quantity}, source := source, dataset := dataset, scope := scope)
                ORDER BY ordinal) AS annotations
            FROM ({annotations}) GROUP BY ALL
        ), evidence AS (
            SELECT resource, relation_id, list(struct_pack(source := source, dataset := dataset,
                row_id := row_id, upstream_id := upstream_id, annotations := annotations,
                subject_molecular_form := subject_molecular_form,
                object_molecular_form := object_molecular_form) ORDER BY ordinal) AS evidence
            FROM ({evidence}) GROUP BY ALL
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
            FROM product_keys k JOIN {self._table("entity", scope)} e USING (entity_key)
            GROUP BY ALL
        ) SELECT rows.* EXCLUDE (resource, relation_id),
            coalesce(products.referenced_product_records, []) AS referenced_product_records
        FROM rows LEFT JOIN products USING (resource, relation_id)"""
        return self._write_parquet(
            sql, [*relation_params, *annotation_params, *evidence_params]
        ), PARQUET

    def _write_parquet(self, sql: str, params: list[Any]) -> bytes:
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            self._db.execute(copy_parquet_command(sql, tmp_path), params)
            return tmp_path.read_bytes()
        finally:
            tmp_path.unlink(missing_ok=True)
