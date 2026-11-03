"""RefMet nomenclature hub, retaining source IDs and exact structure claims."""

from __future__ import annotations

import csv
from collections.abc import Iterable

from ..common import emit_from_records
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

URL = "https://www.metabolomicsworkbench.org/databases/refmet/refmet_download.php"


def emit(writer: HubParquetWriter, records: Iterable[dict] | None = None) -> None:
    # Supplying records permits construction from the cached native export.
    if records is None:
        records = csv.DictReader(iter_url_lines(URL))
    records = ({k.strip(): v for k, v in row.items() if k is not None} for row in records)
    emit_from_records(
        writer,
        records,
        hub_column="refmet_id",
        hub_type="refmet",
        backend="refmet",
        columns={
            "name": "refmet_name",
            "inchikey": "inchi_key",
            "pubchem": "pubchem_cid",
            "chebi": "chebi_id",
            "hmdb": "hmdb_id",
            "lipidmaps": "lipidmaps_id",
            "kegg": "kegg_id",
        },
    )
