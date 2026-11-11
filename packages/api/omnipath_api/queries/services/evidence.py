"""Relation evidence from the step files; raw source payloads from their files."""

from __future__ import annotations

import json
import logging
import time
from typing import Any

logger = logging.getLogger(__name__)

# Relations one evidence summary request may name (the query checks it too).
EVIDENCE_SUMMARY_LIMIT = 1000


def _form_params(filters):
    """The molecular form filters of a request, as the evidence queries take them."""
    from omnipath_api.queries.filters import normalize_filters

    filters = normalize_filters(filters)
    return {
        "protein_entity_keys": filters.get("protein_entity_keys") or None,
        "transcript_entity_keys": filters.get("transcript_entity_keys") or None,
        "isoform_identifiers": filters.get("isoform_identifiers") or None,
        "molecular_endpoint_mode": filters.get("molecular_endpoint_mode") or "any",
    }


class EvidenceQueries:
    """A relation's evidence, evidence summaries of many, raw payloads."""

    def get_relation_evidence(self, relation_pk, resources=None, filters=None):
        text = self.run_query(
            "relation_evidence", {"relation_key": str(relation_pk), **_form_params(filters)}, resources
        )
        return json.loads(text)

    def get_relation_evidence_summaries(self, relation_pks, resources=None, filters=None):
        started = time.perf_counter()
        keys = [str(k) for k in dict.fromkeys(relation_pks or [])]
        text = self.run_query(
            "relation_evidence_batch", {"relation_keys": keys, **_form_params(filters)}, resources
        )
        return {**json.loads(text), "elapsed_ms": round((time.perf_counter() - started) * 1000, 2)}

    def get_relation_payloads(
        self, relation_pk: str, resources: list[str] | None = None
    ) -> dict[str, Any]:
        """Fetch raw input JSON/string records from evidence_payloads.parquet for inspection."""
        rows = self._fetch_dicts(
            f"SELECT source, row_id, payload_json FROM {self._file_table('evidence_payloads', resources)} "
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
