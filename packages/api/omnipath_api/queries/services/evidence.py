"""Evidence query service, including domain SQL and result shaping."""

from __future__ import annotations

import logging
import json
import re
import time
from collections import defaultdict
from typing import Any


from omnipath_api.queries.annotations import (
    split_compiled_annotations,
)


logger = logging.getLogger(__name__)

# Relations one evidence summary request may name.
EVIDENCE_SUMMARY_LIMIT = 1000
SITE_FIELDS = ("term", "residue", "position", "end_position", "description")


def _pmids(value):
    """PubMed ids of a publications value: ``PMID:1``, ``1`` or a list of them split by
    ``;``, ``,`` or ``|``; other references (DOIs) are left out."""
    ids = (
        re.sub(r"^PMID:\s*", "", part.strip(), flags=re.I)
        for part in re.split(r"[;,|]", str(value or ""))
    )
    return [i for i in ids if i.isdigit()]


def _summarize(records):
    """Per evidence source, in first-seen order: the modification sites of either
    endpoint's molecular form and the PubMed ids its records state."""
    by_source: dict[str, dict[str, dict]] = {}
    for record in records:
        item = by_source.setdefault(record["source"] or "unknown", {"sites": {}, "pmids": {}})
        for side in ("subject", "object"):
            for mod in record[f"{side}_modifications"] or []:
                site = {"side": side, **{field: mod.get(field) for field in SITE_FIELDS}}
                item["sites"][tuple(site.values())] = site
        for value in record["publications"]:
            item["pmids"].update(dict.fromkeys(_pmids(value)))
    return [
        {"source": source, "sites": list(item["sites"].values()), "pmids": list(item["pmids"])}
        for source, item in by_source.items()
    ]


class EvidenceQueries:
    """Evidence queries over the release database and the shaping contract."""

    def get_relation_evidence(self, relation_pk, resources=None, filters=None):
        from omnipath_api.queries.molecular import matching_evidence, has_form_filters
        from omnipath_api.queries.filters import normalize_filters
        import hashlib

        rows = [
            {"resource": r["resource"], "relation_id": r["relation_id"]}
            for r in self._relations_by_keys([relation_pk], resources)
        ]
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

    def get_relation_evidence_summaries(self, relation_pks, resources=None, filters=None):
        """{relation pk: [{source, sites, pmids}]} for up to 1,000 relations: the sites and
        PubMed ids that ``get_relation_evidence`` returns for each, from one read.

        Every resource stating a relation is read, as for one relation's evidence; only
        the occurrences matching molecular form filters count. An unknown pk has none.
        """
        from omnipath_api.queries.molecular import form_match_sql, has_form_filters
        from omnipath_api.queries.filters import normalize_filters

        started = time.perf_counter()
        keys = list(dict.fromkeys(relation_pks or []))
        if len(keys) > EVIDENCE_SUMMARY_LIMIT:
            raise ValueError(f"At most {EVIDENCE_SUMMARY_LIMIT} relation_pks per request")
        filters = normalize_filters(filters)
        found = []
        if keys:
            found = self._relations_by_keys(keys, resources)
        key_of = {(r["resource"], r["relation_id"]): r["relation_key"] for r in found}
        if not key_of:
            elapsed = round((time.perf_counter() - started) * 1000, 2)
            return {"evidence": {key: [] for key in keys}, "elapsed_ms": elapsed}
        form = form_match_sql(filters) if has_form_filters(filters) else ("TRUE", [])
        sql, params = self._lookup_sql("relation_evidence", "relation_id", list(key_of), *form)
        # Only the nested columns the summary needs: decoding them is most of the cost.
        records = self._fetch_dicts(
            f"""SELECT resource, relation_id, source, annotations IS NULL AS inherits,
                [a.value FOR a IN coalesce(annotations, []) IF a.term = 'publications'
                    AND coalesce(lower(a.scope), 'relation') = 'relation'] AS publications,
                subject_molecular_form.modifications AS subject_modifications,
                object_molecular_form.modifications AS object_modifications
            FROM ({sql}) ORDER BY resource, relation_id, ordinal""",
            params,
        )
        # A record without its own annotations has its relation's, as in the evidence.
        inheriting = sorted({(r["resource"], r["relation_id"]) for r in records if r["inherits"]})
        inherited = defaultdict(list)
        for row in self._lookup(
            "relation_annotation",
            "relation_id",
            inheriting,
            fields="relation_id, value",
            extra="AND term = 'publications' AND coalesce(lower(scope), 'relation') = 'relation' "
            "ORDER BY relation_id, ordinal",
        ):
            inherited[(row["resource"], row["relation_id"])].append(row["value"])
        by_key = defaultdict(list)
        for record in records:
            pair = (record["resource"], record["relation_id"])
            if record["inherits"]:
                record["publications"] = inherited[pair]
            by_key[key_of[pair]].append(record)
        return {
            "evidence": {key: _summarize(by_key.get(key, [])) for key in keys},
            "elapsed_ms": round((time.perf_counter() - started) * 1000, 2),
        }

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
