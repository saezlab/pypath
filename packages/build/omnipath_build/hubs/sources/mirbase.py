"""miRBase hub from HTML <br>-separated TSV tables."""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from ..stream import iter_url_lines
from ..writer import HubParquetWriter

MIRNA_URL = "https://www.mirbase.org/download/CURRENT/database_files/mirna.txt"
MATURE_URL = "https://www.mirbase.org/download/CURRENT/database_files/mirna_mature.txt"

SPECIES_PREFIX = {"9606": "hsa"}


def _iter_html_rows(lines: Iterable[str]) -> Iterator[list[str]]:
    buf = ""
    for chunk in lines:
        buf += chunk.replace("<BR>", "<br>").replace("<Br>", "<br>")
        while "<br>" in buf:
            row, buf = buf.split("<br>", 1)
            row = row.strip().strip("<>/p")
            if "\t" in row:
                yield row.split("\t")
    if "\t" in buf:
        yield buf.strip().strip("<>/p").split("\t")


def emit(writer: HubParquetWriter, taxonomy_id: str = "9606") -> None:
    prefix = SPECIES_PREFIX.get(str(taxonomy_id), "hsa")
    seen: set[str] = set()
    for row in _iter_html_rows(iter_url_lines(MIRNA_URL)):
        if writer.full:
            return
        name = row[2].strip() if len(row) > 2 else ""
        if prefix and not name.startswith(prefix):
            continue
        hub_id = row[1].strip() if len(row) > 1 else ""
        if not hub_id:
            continue
        if hub_id not in seen:
            seen.add(hub_id)
            if not writer.add_identity("mirbase", hub_id, taxonomy_id, "mirbase"):
                return
        for index in (2, 3):
            value = row[index].strip() if len(row) > index else ""
            if value and not writer.add("mir-pre", value, hub_id, taxonomy_id, "mirbase"):
                return
    for row in _iter_html_rows(iter_url_lines(MATURE_URL)):
        if writer.full:
            return
        name = row[1].strip() if len(row) > 1 else ""
        if prefix and not name.startswith(prefix):
            continue
        hub_id = row[3].strip() if len(row) > 3 else ""
        if not hub_id:
            continue
        if hub_id not in seen:
            seen.add(hub_id)
            if not writer.add_identity("mirbase", hub_id, taxonomy_id, "mirbase"):
                return
        for index in (1, 2):
            value = row[index].strip() if len(row) > index else ""
            if value and not writer.add("mir-mat-name", value, hub_id, taxonomy_id, "mirbase"):
                return
