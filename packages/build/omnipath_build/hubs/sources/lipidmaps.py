"""LIPID MAPS hub from the SDF inside the zip download."""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from ..common import emit_from_records
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

URL = "https://lipidmaps.org/files/?file=LMSD&ext=sdf.zip"


def _iter_sdf_records(lines: Iterable[str]) -> Iterator[dict[str, str]]:
    record: dict[str, str] = {}
    field: str | None = None
    buf: list[str] = []
    for raw in lines:
        line = raw.rstrip("\n")
        if line.startswith("$$$$"):
            if field is not None:
                record[field] = "\n".join(buf).strip()
            if record:
                yield record
            record, field, buf = {}, None, []
            continue
        if line.startswith("> <") and line.endswith(">"):
            if field is not None:
                record[field] = "\n".join(buf).strip()
            field = line[3:-1]
            buf = []
            continue
        if field is not None:
            if line == "":
                record[field] = "\n".join(buf).strip()
                field, buf = None, []
            else:
                buf.append(line)
    if field is not None:
        record[field] = "\n".join(buf).strip()
    if record:
        yield record


def emit(writer: HubParquetWriter) -> None:
    emit_from_records(
        writer,
        _iter_sdf_records(iter_url_lines(URL)),
        hub_column="LM_ID",
        hub_type="lipidmaps",
        backend="lipidmaps",
        columns={
            "name": "NAME",
            "lipid_shorthand": "ABBREVIATION",
            "systematic_name": "SYSTEMATIC_NAME",
            "inchikey": "INCHI_KEY",
            "inchi": "INCHI",
            "smiles": "SMILES",
            "chebi": "CHEBI_ID",
            "pubchem": "PUBCHEM_CID",
            "hmdb": "HMDB_ID",
            "swisslipids": "SWISSLIPIDS_ID",
            "synonym": "SYNONYMS",
        },
        extra=lambda row: {"name": row.get("NAME") or row.get("COMMON_NAME")},
    )
