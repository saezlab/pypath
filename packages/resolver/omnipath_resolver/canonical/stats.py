"""Resolution statistics per entity type and per matching rule."""

from __future__ import annotations

from collections import defaultdict
from typing import Any
import sqlite3
import tempfile
from pathlib import Path

from .match import Match
from .policy import CHEMICAL, EntityPolicy

# Entity types whose resolution is not built yet; unmatched ones are out of scope rather than
# failures. A microRNA that does match a gene still counts as resolved.
RESOLUTION_PENDING_TYPES = frozenset({"microrna"})


def _normalized_type(entity_type: str | None) -> str:
    return str(entity_type or "").strip().lower().replace(" ", "_").replace("-", "_")


class ResolutionTracker:
    """Counts unique entity keys by outcome.

    * ``resolved``        – matched a library node
    * ``structure``       – a chemical identified by one valid InChIKey (supplied or derived
      from SMILES) that the reference does not hold yet
    * ``unresolved``      – class is matched against a library, nothing matched
    * ``not_applicable``  – class has no library (cv_term, complex, generic …), or its
      resolution does not exist yet (microRNAs)
    """

    def __init__(self) -> None:
        self.input_entities = 0
        self.resolved_entities = 0
        self.structure_entities = 0
        self.unresolved_entities = 0
        self.not_applicable_entities = 0
        self.by_entity_type: dict[str, dict[str, int]] = defaultdict(
            lambda: {
                "input": 0,
                "resolved": 0,
                "structure": 0,
                "unresolved": 0,
                "not_applicable": 0,
            }
        )
        self.by_rule: dict[str, int] = defaultdict(int)
        self._storage = tempfile.TemporaryDirectory(prefix="omnipath-resolution-")
        self._keys = sqlite3.connect(str(Path(self._storage.name) / "keys.sqlite"))
        self._keys.execute("PRAGMA cache_size=-4096")
        self._keys.execute(
            "CREATE TABLE outcomes (key TEXT PRIMARY KEY, entity_type TEXT, status TEXT, rule TEXT) WITHOUT ROWID"
        )

    def record(self, entity_key: str, entity_type: str, policy: EntityPolicy, match: Match) -> None:
        if not policy.matched:
            status = "not_applicable"
        elif match.matched:
            status = "resolved"
        elif _normalized_type(entity_type) in RESOLUTION_PENDING_TYPES:
            status = "not_applicable"
        elif policy.library == CHEMICAL and match.canonical_namespace == "inchikey":
            status = "structure"
        else:
            status = "unresolved"
        rule = f"{policy.entity_class}/{match.resolved_by}"

        previous = self._keys.execute(
            "SELECT entity_type,status,rule FROM outcomes WHERE key=?", [entity_key]
        ).fetchone()
        if previous is not None:
            prev_type, prev_status, prev_rule = previous
            if (prev_type, prev_status) == (entity_type, status):
                return
            self._apply(prev_type, prev_status, -1)
            self.by_rule[prev_rule] -= 1
        self._keys.execute(
            "INSERT OR REPLACE INTO outcomes VALUES (?,?,?,?)",
            [entity_key, entity_type, status, rule],
        )
        self._apply(entity_type, status, 1)
        self.by_rule[rule] += 1

    def add_counts(self, entity_type: str, status: str, rule: str, count: int) -> None:
        """Import already deduplicated aggregate outcomes."""
        self._apply(entity_type, status, count)
        self.by_rule[rule] += count

    def export_keys(self, path: str | Path) -> Path:
        """Stream distinct outcomes to a portable shard without retaining keys in RAM."""
        import pyarrow as pa
        import pyarrow.parquet as pq

        fields = ("key", "entity_type", "status", "rule")
        schema = pa.schema([(name, pa.string()) for name in fields])
        cursor = self._keys.execute("SELECT key,entity_type,status,rule FROM outcomes")
        with pq.ParquetWriter(path, schema, compression="zstd") as writer:
            while batch := cursor.fetchmany(8192):
                writer.write_table(
                    pa.Table.from_pylist([dict(zip(fields, row)) for row in batch], schema=schema)
                )
        return Path(path)

    def close(self) -> None:
        self._keys.close()
        self._storage.cleanup()

    def _apply(self, entity_type: str, status: str, delta: int) -> None:
        self.input_entities += delta
        if status == "resolved":
            self.resolved_entities += delta
        elif status == "structure":
            self.structure_entities += delta
        elif status == "unresolved":
            self.unresolved_entities += delta
        else:
            self.not_applicable_entities += delta
        counts = self.by_entity_type[entity_type]
        counts["input"] += delta
        counts[status] += delta

    def summary(self) -> dict[str, Any]:
        return {
            "scope": "unique_entity_keys",
            "input_entities": self.input_entities,
            "resolved_entities": self.resolved_entities,
            "structure_entities": self.structure_entities,
            "unresolved_entities": self.unresolved_entities,
            "not_applicable_entities": self.not_applicable_entities,
            "lookup_entities": self.resolved_entities
            + self.structure_entities
            + self.unresolved_entities,
            "by_entity_type": {k: dict(v) for k, v in sorted(self.by_entity_type.items())},
            "by_rule": {k: v for k, v in sorted(self.by_rule.items()) if v},
        }
