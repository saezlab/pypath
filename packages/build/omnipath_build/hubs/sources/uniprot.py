"""UniProt AC hub from the full FTP idmapping file plus secondary ACs.

Streams ``idmapping.dat.gz`` (all organisms) and keeps only the protein
reference identifier types from ``pypath.inputs_v2.uniprot``:
UniProt, entry name, gene name / synonym, Entrez, Ensembl, HGNC.
PDB, EMBL, GI, DrugBank, GO, UniRef, and similar xrefs are dropped.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Iterator

from ..stream import iter_url_lines
from ..writer import HubParquetWriter

IDMAPPING_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/idmapping/idmapping.dat.gz"
)
SEC_AC_URL = (
    "https://ftp.expasy.org/databases/uniprot/current_release/"
    "knowledgebase/complete/docs/sec_ac.txt"
)

# Keep in sync with pypath.inputs_v2.uniprot.PROTEIN_REFERENCE_KEY_TYPES
# and omnipath-utils resolver_protein (ensg/ensp/enst split).
IDTYPE_MAP = {
    "UniProtKB-ID": "uniprot_entry",
    "Gene_Name": "genesymbol",
    "Gene_Synonym": "genesymbol-syn",
    "GeneID": "entrez",
    "Ensembl": "ensg",
    "Ensembl_TRS": "enst",
    "Ensembl_PRO": "ensp",
    "HGNC": "hgnc",
    "RefSeq": "refseq_protein",
    "RefSeq_NT": "refseq",
    "EMBL-CDS": "genbank",
}

_ENSEMBL_FTP_TYPES = frozenset({"Ensembl", "Ensembl_TRS", "Ensembl_PRO"})

PROTEIN_SOURCE_TYPES: frozenset[str] = frozenset(
    {
        "uniprot",
        "uniprot-sec",
        "uniprot_entry",
        "genesymbol",
        "genesymbol-syn",
        "entrez",
        "ensembl",
        "ensg",
        "enst",
        "ensp",
        "hgnc",
        "refseq",
        "refseq_protein",
        "genbank",
    }
)


def _strip_ensembl_version(value: str) -> str:
    if "." not in value:
        return value
    head, tail = value.rsplit(".", 1)
    return head if tail.isdigit() else value


def emit_uniprot_secondary(writer: HubParquetWriter) -> None:
    for index, line in enumerate(iter_url_lines(SEC_AC_URL)):
        if writer.full:
            return
        if index < 30:
            continue
        fields = line.split()
        if len(fields) != 2:
            continue
        if not writer.add("uniprot-sec", fields[0].strip(), fields[1].strip(), "0", "uniprot"):
            return


def _iter_idmapping_lines(lines: Iterable[str] | None) -> Iterator[str]:
    if lines is not None:
        yield from lines
        return
    local = os.environ.get("OMNIPATH_BUILD_FTP_FILE")
    if local:
        import gzip

        with gzip.open(local, mode="rt", encoding="utf-8", errors="replace") as handle:
            yield from handle
        return
    yield from iter_url_lines(IDMAPPING_URL)


def emit(
    writer: HubParquetWriter,
    lines: Iterable[str] | None = None,
    *,
    include_secondary: bool | None = None,
) -> None:
    """Stream the full UniProt idmapping file onto primary accessions.

    Rows are grouped by AC so ``NCBI_TaxID`` can be applied before flush.
    Only ``IDTYPE_MAP`` types are written. ``include_secondary`` defaults to
    true when reading the live FTP stream.
    """
    if include_secondary is None:
        include_secondary = lines is None

    source = _iter_idmapping_lines(lines)
    closer = getattr(source, "close", None)
    current_ac = ""
    taxid = "0"
    pending: list[tuple[str, str]] = []

    def flush() -> bool:
        if not current_ac:
            return True
        if not writer.add_identity("uniprot", current_ac, taxid, "uniprot_ftp"):
            return False
        for source_type, value in pending:
            if not writer.add(source_type, value, current_ac, taxid, "uniprot_ftp"):
                return False
        return True

    try:
        for line in source:
            if writer.full:
                return
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            accession, ftp_type, value = (
                fields[0].strip(),
                fields[1].strip(),
                fields[2].strip(),
            )
            if not accession or not value:
                continue
            if accession != current_ac:
                if not flush():
                    return
                current_ac = accession
                taxid = "0"
                pending = []
            if ftp_type == "NCBI_TaxID":
                if value.isdigit():
                    taxid = value
                continue
            source_type = IDTYPE_MAP.get(ftp_type)
            if source_type is None:
                continue
            if ftp_type in _ENSEMBL_FTP_TYPES:
                value = _strip_ensembl_version(value)
            pending.append((source_type, value))
        if not writer.full:
            flush()
    finally:
        if closer is not None:
            closer()

    if include_secondary and not writer.full:
        emit_uniprot_secondary(writer)
