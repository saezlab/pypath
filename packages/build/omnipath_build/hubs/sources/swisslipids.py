"""SwissLipids hub from the gzip TSV export."""

from __future__ import annotations

from ..common import emit_from_records, iter_tsv_dicts
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

URL = "https://swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv"


def emit(writer: HubParquetWriter) -> None:
    emit_from_records(
        writer,
        iter_tsv_dicts(iter_url_lines(URL, encoding="latin-1")),
        hub_column="Lipid ID",
        hub_type="swisslipids",
        backend="swisslipids",
        columns={
            "name": "Name",
            "lipid_shorthand": "Abbreviation*",
            "inchikey": "InChI key (pH7.3)",
            "inchi": "InChI (pH7.3)",
            "smiles": "SMILES (pH7.3)",
            "chebi": "CHEBI",
            "lipidmaps": "LIPID MAPS",
            "hmdb": "HMDB",
            "metanetx": "MetaNetX",
            "synonym": "Synonyms*",
        },
    )
