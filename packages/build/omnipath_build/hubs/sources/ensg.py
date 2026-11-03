"""Ensembl gene hub from NCBI gene2ensembl.gz (gzip-streamed)."""

from __future__ import annotations

from ..common import explode_record
from ..schema import as_values
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"


def emit(writer: HubParquetWriter) -> None:
    seen: set[tuple[str, str]] = set()
    for line in iter_url_lines(URL):
        if writer.full:
            return
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue
        taxon = fields[0] if fields[0].isdigit() else "0"
        entrez = "" if fields[1] == "-" else fields[1].strip()
        ensg = "" if fields[2] in {"", "-"} else fields[2].split(".", 1)[0]
        ensp = "" if fields[6] in {"", "-"} else fields[6].split(".", 1)[0]
        if not ensg or not ensg.upper().startswith("ENS"):
            continue
        key = (taxon, ensg)
        field_map = {"entrez": entrez, "ensp": ensp}
        if key in seen:
            for source_type, raw in field_map.items():
                for value in as_values(raw):
                    if not writer.add(source_type, value, ensg, taxon, "gene2ensembl"):
                        return
            continue
        seen.add(key)
        if not explode_record(
            writer,
            hub_id=ensg,
            hub_type="ensg",
            backend="gene2ensembl",
            taxonomy_id=taxon,
            fields=field_map,
        ):
            return
