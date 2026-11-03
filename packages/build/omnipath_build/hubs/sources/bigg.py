"""BiGG metabolite hub from universal metabolites TSV."""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Iterable

from ..common import explode_record
from ..schema import CHEMICAL_TAXON, as_values
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

URL = "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt"

_DB = {
    "CHEBI": "chebi",
    "Human Metabolome Database": "hmdb",
    "KEGG Compound": "kegg",
    "KEGG Drug": "kegg",
    "KEGG Glycan": "kegg_glycan",
    "MetaNetX (MNX) Chemical": "metanetx",
    "SEED Compound": "seed",
    "LipidMaps": "lipidmaps",
    "BioCyc": "biocyc",
    "Reactome Compound": "reactome",
    "InChI Key": "inchikey",
}
_RE_URL_ID = re.compile(r".*/([^/]+)$")


def _parse_links(db_links: str) -> dict[str, set[str]]:
    result: dict[str, set[str]] = defaultdict(set)
    if not db_links:
        return result
    for entry in db_links.split(";"):
        entry = entry.strip()
        if ": " not in entry:
            continue
        db_name, url = entry.split(": ", 1)
        id_type = _DB.get(db_name.strip())
        if id_type is None:
            continue
        match = _RE_URL_ID.match(url.strip())
        if match:
            result[id_type].add(match.group(1))
    return result


def emit(writer: HubParquetWriter, lines: Iterable[str] | None = None) -> None:
    header_skipped = False
    for line in lines if lines is not None else iter_url_lines(URL):
        if writer.full:
            return
        if not header_skipped:
            header_skipped = True
            if "universal_bigg_id" in line.lower() or line.lower().startswith("bigg_id"):
                continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        bigg_id = fields[1].strip()
        if not bigg_id:
            continue
        xrefs: dict = dict(_parse_links(fields[4].strip() if len(fields) > 4 else ""))
        if len(fields) > 2 and fields[2].strip():
            xrefs["name"] = fields[2].strip()
        model_id = fields[0].strip()
        if model_id and model_id != bigg_id:
            existing = {str(item) for item in as_values(xrefs.get("bigg"))}
            existing.add(model_id)
            xrefs["bigg"] = existing
        if not explode_record(
            writer,
            hub_id=bigg_id,
            hub_type="bigg",
            backend="bigg",
            taxonomy_id=CHEMICAL_TAXON,
            fields=xrefs,
        ):
            return
