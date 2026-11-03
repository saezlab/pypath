"""ChEBI hub from chebi.obo.gz, parsed incrementally."""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from typing import Any

from ..common import emit_from_records
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz"

_RE_SYNONYM = re.compile(r'^"((?:\\.|[^"\\])*)"')
_KEGG_RE = re.compile(r"(?:kegg(?:\s+compound)?|kegg\.compound)[:\s]+(C\d+)", re.I)
_PUBCHEM_RE = re.compile(r"(?:pubchem(?:\s+compound)?|cid)[:\s]+(\d+)", re.I)
_HMDB_RE = re.compile(r"(HMDB\d+)", re.I)
_LIPIDMAPS_RE = re.compile(r"(?:lipid\s*maps|lipidmaps)[:\s]+([A-Z0-9]+)", re.I)
_CAS_RE = re.compile(r"\b(\d{2,7}-\d{2}-\d)\b")
_PROPERTY_RE = re.compile(r"^(\S+)\s+\"(.*?)\"")


def _iter_obo_terms(lines: Iterable[str]) -> Iterator[dict[str, Any]]:
    term: dict[str, Any] | None = None
    for raw in lines:
        line = raw.strip()
        if line == "[Term]":
            if term and term.get("chebi_id"):
                yield term
            term = {
                "chebi_id": "",
                "name": "",
                "alt_ids": [],
                "synonyms": [],
                "inchikey": "",
                "inchi": "",
                "smiles": "",
                "formula": "",
                "kegg_compound": [],
                "pubchem_compound": [],
                "hmdb": [],
                "lipidmaps": [],
                "cas": [],
            }
            continue
        if term is None or not line or line.startswith("["):
            continue
        if line.startswith("id: "):
            term["chebi_id"] = line[4:].strip()
        elif line.startswith("name: "):
            term["name"] = line[6:].strip()
        elif line.startswith("alt_id: "):
            term["alt_ids"].append(line[8:].strip())
        elif line.startswith("synonym: "):
            match = _RE_SYNONYM.match(line[9:].strip())
            if match:
                term["synonyms"].append(match.group(1).replace('\\"', '"'))
        elif line.startswith("xref: "):
            xref = line[6:].split("!")[0].strip()
            if match := _KEGG_RE.search(xref):
                term["kegg_compound"].append(match.group(1).upper())
            if match := _PUBCHEM_RE.search(xref):
                term["pubchem_compound"].append(match.group(1))
            if match := _HMDB_RE.search(xref):
                term["hmdb"].append(match.group(1).upper())
            if match := _LIPIDMAPS_RE.search(xref):
                term["lipidmaps"].append(match.group(1).upper())
            if match := _CAS_RE.search(xref):
                term["cas"].append(match.group(1))
        elif line.startswith("property_value: "):
            match = _PROPERTY_RE.match(line[16:].strip())
            if not match:
                continue
            key, value = match.group(1), match.group(2)
            if key.endswith("inchi_key_string") or key.endswith("inchikey"):
                term["inchikey"] = value
            elif key.endswith("inchi_string") or key.endswith("inchi"):
                term["inchi"] = value
            elif key.endswith("smiles_string") or key.endswith("smiles"):
                term["smiles"] = value
            elif key.endswith("formula") or "empirical_formula" in key:
                term["formula"] = value
    if term and term.get("chebi_id"):
        yield term


def emit(writer: HubParquetWriter) -> None:
    emit_from_records(
        writer,
        _iter_obo_terms(iter_url_lines(URL)),
        hub_column="chebi_id",
        hub_type="chebi",
        backend="chebi",
        columns={
            "inchikey": "inchikey",
            "inchi": "inchi",
            "smiles": "smiles",
            "name": "name",
            "kegg": "kegg_compound",
            "pubchem": "pubchem_compound",
            "hmdb": "hmdb",
            "lipidmaps": "lipidmaps",
            "cas": "cas",
            "formula": "formula",
        },
        extra=lambda row: {
            "chebi": row.get("alt_ids") or [],
            "synonym": row.get("synonyms") or [],
        },
    )
