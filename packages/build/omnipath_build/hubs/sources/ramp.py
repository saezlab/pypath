"""Normalize cached RaMP source assertions without consulting resolved outputs."""

from ..common import explode_record

NAMESPACES = {
    "hmdb": "hmdb",
    "lipidmaps": "lipidmaps",
    "cas": "cas",
    "chebi": "chebi",
    "chemspider": "chemspider",
    "ensembl": "ensg",
    "entrez": "entrez",
    "gene_symbol": "genesymbol",
    "kegg": "kegg",
    "reactome": "reactome",
    "wiki": "wikipathways",
    "pfocr": "pfocr",
    "pubchem": "pubchem",
    "chembl": "chembl",
    "uniprot": "uniprot",
    "en": "ec",
    "brenda": "ec",
}


def assertions(dataset, row):
    """Yield domain, native ID, taxon and source fields; names never imply identity."""
    ident = row.get("rampId") or row.get("ramp_id")
    if not ident or not ident.startswith(("RAMP_C_", "RAMP_G_")):
        return
    hub = "ramp" if ident.startswith("RAMP_C_") else "ramp_gene"
    fields = {}
    if dataset == "analyte":
        fields["name"] = row.get("common_name")
    elif dataset in ("source", "metabolite_class", "chem_props"):
        ns = NAMESPACES.get(
            str(row.get("IDtype") or row.get("source") or row.get("chem_data_source") or "").lower()
        )
        value = row.get("sourceId") or row.get("class_source_id") or row.get("chem_source_id")
        if ns and value:
            fields[ns] = (
                str(value).split(":", 1)[-1] if ns != "kegg" or hub == "ramp" else str(value)
            )
        fields["name"] = row.get("commonName") or row.get("common_name")
        if dataset == "chem_props":
            fields.update(
                inchikey=row.get("inchi_key"), inchi=row.get("inchi"), smiles=row.get("iso_smiles")
            )
    elif dataset != "analytehaspathway":
        return
    yield hub, ident, "0", fields


def emit(writer, records, *, domain="chemical"):
    selected = "ramp" if domain == "chemical" else "ramp_gene"
    for dataset, row in records:
        if writer.full:
            return
        for hub, ident, taxon, fields in assertions(dataset, row):
            if hub == selected:
                explode_record(
                    writer,
                    hub_id=ident,
                    hub_type=hub,
                    backend="ramp",
                    taxonomy_id=taxon,
                    fields=fields,
                )


def cached_records():
    """Read the native RaMP download, never a resolved resource projection."""
    import re
    import sqlite3
    from contextlib import closing
    from ...discovery import setup_pypath_cache

    cache = setup_pypath_cache() / "ramp"
    candidates = list(cache.glob("RaMP_SQLite_v*.sqlite"))
    if candidates:
        path = max(candidates, key=lambda p: tuple(int(n) for n in re.findall(r"\d+", p.stem)))
        with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)) as con:
            con.row_factory = sqlite3.Row
            for table in (
                "analyte",
                "source",
                "metabolite_class",
                "chem_props",
                "analytehaspathway",
            ):
                for row in con.execute(f'SELECT * FROM "{table}"'):
                    yield table, dict(row)
        return
    from pypath.inputs_v2 import rampdb

    opener = rampdb.download.open()
    for table in ("analyte", "source", "metabolite_class", "chem_props", "analytehaspathway"):
        for row in rampdb.parser(opener, table=table):
            yield table, row


def emit_chemical(writer):
    emit(writer, cached_records(), domain="chemical")


def emit_gene(writer):
    emit(writer, cached_records(), domain="gene_protein")
