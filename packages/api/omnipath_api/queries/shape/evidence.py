"""Shape raw DuckDB evidence and annotations into API models."""

from __future__ import annotations

from typing import Any


def shape_evidence_item(item: dict[str, Any]) -> dict[str, Any]:
    """Format single evidence item with source, dataset, row_id, upstream_id, and annotations."""
    raw_anns = item.get("annotations") or []
    annotations: list[dict[str, Any]] = []
    if isinstance(raw_anns, list):
        for a in raw_anns:
            if isinstance(a, dict) and a.get("term"):
                annotations.append(
                    {
                        "term": str(a.get("term")),
                        "value": a.get("value"),
                        **({"quantity": a["quantity"]} if a.get("quantity") is not None else {}),
                        "source": a.get("source"),
                        "dataset": a.get("dataset"),
                        "scope": a.get("scope"),
                    }
                )

    return {
        "source": str(item.get("source") or ""),
        "dataset": str(item.get("dataset") or ""),
        "rowId": str(item.get("row_id") or ""),
        "upstreamId": str(item.get("upstream_id") or ""),
        "annotations": annotations,
        "subjectMolecularForm": item.get("subject_molecular_form"),
        "objectMolecularForm": item.get("object_molecular_form"),
        "molecularForm": item.get("molecular_form"),
    }
