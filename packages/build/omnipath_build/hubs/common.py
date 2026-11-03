"""Row helpers shared by every hub source."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Any, Callable

from .schema import CHEMICAL_TAXON, as_values
from .writer import HubParquetWriter


def explode_record(
    writer: HubParquetWriter,
    *,
    hub_id: str,
    hub_type: str,
    backend: str,
    taxonomy_id: str,
    fields: dict[str, Any],
) -> bool:
    if not writer.add_identity(hub_type, hub_id, taxonomy_id, backend):
        return False
    for source_type, raw in fields.items():
        for value in as_values(raw):
            if value == hub_id and source_type == hub_type:
                continue
            if not writer.add(source_type, value, hub_id, taxonomy_id, backend):
                return False
    return True


def emit_from_records(
    writer: HubParquetWriter,
    records: Iterable[dict[str, Any]],
    *,
    hub_column: str,
    hub_type: str,
    backend: str,
    columns: dict[str, str],
    taxonomy_id: str = CHEMICAL_TAXON,
    extra: Callable[[dict[str, Any]], dict[str, Any]] | None = None,
) -> None:
    for record in records:
        if writer.full:
            return
        hub_ids = as_values(record.get(hub_column))
        if not hub_ids:
            continue
        fields = {source_type: record.get(column) for source_type, column in columns.items()}
        if extra:
            fields.update(extra(record))
        for hub_id in hub_ids:
            if not explode_record(
                writer,
                hub_id=hub_id,
                hub_type=hub_type,
                backend=backend,
                taxonomy_id=taxonomy_id,
                fields=fields,
            ):
                return


def iter_tsv_dicts(lines: Iterable[str]) -> Iterator[dict[str, str]]:
    header: list[str] | None = None
    for line in lines:
        stripped = line.rstrip("\n")
        if not stripped:
            continue
        if header is None:
            if stripped.startswith("#"):
                continue
            header = stripped.split("\t")
            continue
        fields = stripped.split("\t")
        yield {
            key: fields[index] if index < len(fields) else "" for index, key in enumerate(header)
        }
