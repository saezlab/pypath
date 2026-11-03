"""HMDB hub from the metabolites XML inside the zip download."""

from __future__ import annotations

import io
import xml.etree.ElementTree as ET
from collections.abc import Iterable, Iterator

from ..common import emit_from_records
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

URL = "https://rescued.omnipathdb.org/hmdb_metabolites.zip"
_NS = "{http://www.hmdb.ca}"


class _TextBytes(io.RawIOBase):
    def __init__(self, source: Iterable[str]) -> None:
        self._it = iter(source)
        self._buf = b""

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:  # type: ignore[no-untyped-def]
        wanted = len(buffer)
        while len(self._buf) < wanted:
            try:
                self._buf += next(self._it).encode("utf-8")
            except StopIteration:
                break
        out = self._buf[:wanted]
        self._buf = self._buf[wanted:]
        buffer[: len(out)] = out
        return len(out)


def _iter_metabolites(lines: Iterable[str]) -> Iterator[dict[str, str]]:
    context = ET.iterparse(_TextBytes(lines), events=("end",))
    for _, elem in context:
        if elem.tag != f"{_NS}metabolite":
            continue

        def text(tag: str) -> str:
            child = elem.find(f"{_NS}{tag}")
            return (child.text or "").strip() if child is not None else ""

        yield {
            "accession": text("accession"),
            "name": text("name"),
            "inchikey": text("inchikey"),
            "inchi": text("inchi"),
            "smiles": text("smiles"),
            "chebi_id": text("chebi_id"),
            "pubchem_compound_id": text("pubchem_compound_id"),
            "kegg_id": text("kegg_id"),
            "drugbank_id": text("drugbank_id"),
            "cas_registry_number": text("cas_registry_number"),
        }
        elem.clear()


def emit(writer: HubParquetWriter) -> None:
    emit_from_records(
        writer,
        _iter_metabolites(iter_url_lines(URL)),
        hub_column="accession",
        hub_type="hmdb",
        backend="hmdb",
        columns={
            "name": "name",
            "inchikey": "inchikey",
            "inchi": "inchi",
            "smiles": "smiles",
            "chebi": "chebi_id",
            "pubchem": "pubchem_compound_id",
            "kegg": "kegg_id",
            "drugbank": "drugbank_id",
            "cas": "cas_registry_number",
        },
    )
