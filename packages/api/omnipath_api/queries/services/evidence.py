"""Evidence query service, including domain SQL and result shaping."""

from __future__ import annotations

import logging
import json
from typing import Any


from omnipath_api.queries.annotations import (
    split_compiled_annotations,
)


logger = logging.getLogger(__name__)


class EvidenceQueries:
    """Evidence queries over ParquetQueries' storage and shaping contract."""

    def get_relation_evidence(self, relation_pk, resources=None, filters=None):
        from omnipath_api.queries.molecular import matching_evidence, has_form_filters
        from omnipath_api.queries.filters import normalize_filters
        import hashlib

        rows = self._fetch_dicts(
            f"SELECT resource, relation_id FROM {self._table('relation', resources)} "
            "WHERE relation_key = ? ORDER BY resource",
            [relation_pk],
        )
        self._children("relation_annotation", rows, "annotations")
        self._children("relation_evidence", rows, "evidence")
        filters = normalize_filters(filters)
        out, annotations = [], []
        for row in rows:
            compiled = split_compiled_annotations(row.get("annotations") or [])
            annotations.extend(compiled["relation"])
            items = row.get("evidence") or []
            for index, record in enumerate(items):
                if has_form_filters(filters) and not matching_evidence([record], filters):
                    continue
                item_annotations = record.get("annotations")
                parts = (
                    split_compiled_annotations(item_annotations)
                    if item_annotations is not None
                    else compiled
                )
                identity = json.dumps(
                    [
                        row["resource"],
                        row["relation_id"],
                        record.get("source"),
                        record.get("dataset"),
                        record.get("row_id"),
                        index,
                    ]
                )
                out.append(
                    {
                        "relationPk": relation_pk,
                        "relationEvidencePk": relation_pk
                        + ":"
                        + hashlib.sha256(identity.encode()).hexdigest(),
                        "source": record.get("source") or "unknown",
                        "dataset": record.get("dataset"),
                        "rowId": record.get("row_id"),
                        "upstreamId": record.get("upstream_id"),
                        "subjectMolecularForm": record.get("subject_molecular_form"),
                        "objectMolecularForm": record.get("object_molecular_form"),
                        "subjectAttributes": parts["subject"],
                        "objectAttributes": parts["object"],
                        "recordAttributes": parts["relation"],
                        "evidence": parts["relation"],
                    }
                )
        return {"evidence": out, "annotations": annotations}

    def get_relation_payloads(
        self, relation_pk: str, resources: list[str] | None = None
    ) -> dict[str, Any]:
        """Fetch raw input JSON/string records from evidence_payloads.parquet for inspection."""
        rows = self._fetch_dicts(
            f"SELECT source, row_id, payload_json FROM {self._table('evidence_payloads', resources)} "
            "WHERE relation_key = ?",
            [relation_pk],
        )
        out = []
        for r in rows:
            raw = r.get("payload_json")
            parsed = raw
            if isinstance(raw, str):
                try:
                    parsed = json.loads(raw)
                except json.JSONDecodeError:
                    logger.warning(
                        "Non-JSON evidence payload for relation %s, row %s",
                        relation_pk,
                        r.get("row_id"),
                    )
                    parsed = {"raw": raw}
            out.append(
                {
                    "source": str(r.get("source") or "unknown"),
                    "rowId": str(r.get("row_id") or ""),
                    "payload": parsed,
                }
            )
        return {"relationPk": relation_pk, "payloads": out}
