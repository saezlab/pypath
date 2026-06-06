#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""UniProt FTP ID mapping file access.

Downloads and parses the UniProt idmapping.dat.gz files -- the most
comprehensive source for UniProt cross-references. Available as a
complete dump (all organisms, ~18GB) or per-organism files (~100MB each).
"""

from __future__ import annotations

import gzip
import logging
import os
from collections import defaultdict

from pypath.share.downloads import dm

_log = logging.getLogger(__name__)

FTP_BASES = (
    'https://ftp.uniprot.org/pub/databases/uniprot',
    'https://ftp.ebi.ac.uk/pub/databases/uniprot',
    'https://ftp.expasy.org/databases/uniprot',
)

IDMAPPING_PATH = 'current_release/knowledgebase/idmapping'

# ID type names in the idmapping.dat file -> our canonical names
IDTYPE_MAP = {
    'UniProtKB-ID': 'uniprot_entry',
    'Gene_Name': 'genesymbol',
    'GeneID': 'entrez',
    'RefSeq': 'refseqp',
    'RefSeq_NT': 'refseqn',
    'GI': 'gi',
    'PDB': 'pdb',
    'GO': 'go',
    'UniRef100': 'uniref100',
    'UniRef90': 'uniref90',
    'UniRef50': 'uniref50',
    'UniParc': 'uniparc',
    'EMBL-CDS': 'embl',
    'EMBL': 'embl_id',
    'Ensembl': 'ensg',
    'Ensembl_TRS': 'enst',
    'Ensembl_PRO': 'ensp',
    'HGNC': 'hgnc',
    'KEGG': 'kegg',
    'NCBI_TaxID': '_taxid',  # special: not an ID mapping, it's the organism
    'ChEMBL': 'chembl',
    'DrugBank': 'drugbank',
    'Reactome': 'reactome',
    'STRING': 'string',
}


_CODES = {
    9606: 'HUMAN', 10090: 'MOUSE', 10116: 'RAT',
    559292: 'YEAST', 83333: 'ECOLI', 7227: 'DROME',
    7955: 'DANRE', 6239: 'CAEEL', 9031: 'CHICK',
    3702: 'ARATH', 44689: 'DICDI', 284812: 'SCHPO',
}


def organism_urls(ncbi_tax_id: int) -> list[str]:
    """Get URLs for per-organism idmapping files (primary + mirrors)."""

    code = _CODES.get(ncbi_tax_id)

    if not code:
        return []

    path = f'{IDMAPPING_PATH}/by_organism/{code}_{ncbi_tax_id}_idmapping.dat.gz'

    return [f'{base}/{path}' for base in FTP_BASES]


def idmapping_stream(
    ncbi_tax_id: int | None = None,
    id_types: set[str] | None = None,
):
    """Stream ID mapping records from a UniProt FTP file.

    Args:
        ncbi_tax_id: If given, download the per-organism file. If None,
            download the full idmapping.dat.gz (very large!).
        id_types: If given, only yield rows for these ID type names
            (as they appear in the file, e.g. 'Gene_Name', 'GeneID').
            Pass None to get all types.

    Yields:
        Tuples of (uniprot_ac, id_type_name, id_value).
    """

    if ncbi_tax_id:
        urls = organism_urls(ncbi_tax_id)
    else:
        urls = [f'{base}/{IDMAPPING_PATH}/idmapping.dat.gz' for base in FTP_BASES]

    if not urls:
        _log.warning('No FTP URL for taxid %s', ncbi_tax_id)
        return

    try:
        path = dm.download(urls[0], fallback_urls=urls[1:], connecttimeout=10)
    except Exception as e:
        _log.error('UniProt FTP download failed: %s', e)
        return

    if not path or os.path.getsize(path) == 0:
        _log.error('All UniProt FTP mirrors failed')
        return

    _log.info('Parsing %s', path)

    with gzip.open(path, 'rt') as f:

        for line in f:
            parts = line.rstrip('\n').split('\t')

            if len(parts) != 3:
                continue

            uniprot_ac, id_type_name, id_value = parts

            if id_types and id_type_name not in id_types:
                continue

            yield uniprot_ac, id_type_name, id_value


def idmapping(
    id_type: str,
    ncbi_tax_id: int = 9606,
) -> dict[str, set[str]]:
    """Load a mapping table from the FTP idmapping file.

    Args:
        id_type: The target ID type name as it appears in the file
            (e.g. 'Gene_Name', 'GeneID', 'Ensembl').
        ncbi_tax_id: Organism.

    Returns:
        Dict mapping UniProt AC -> set of target IDs.
    """

    data = defaultdict(set)

    for uniprot_ac, id_type_name, id_value in idmapping_stream(
        ncbi_tax_id=ncbi_tax_id,
        id_types={id_type},
    ):
        data[uniprot_ac].add(id_value)

    _log.info(
        'FTP idmapping: %d UniProt ACs with %s for organism %d',
        len(data), id_type, ncbi_tax_id,
    )

    return dict(data)


def all_id_types(ncbi_tax_id: int = 9606) -> set[str]:
    """Get all ID type names present in the file for an organism."""

    types = set()

    for _, id_type_name, _ in idmapping_stream(ncbi_tax_id=ncbi_tax_id):
        types.add(id_type_name)

    return types


def full_idmapping_urls() -> list[str]:
    """URLs for the complete idmapping.dat.gz (all organisms)."""
    return [
        f"{base}/{IDMAPPING_PATH}/idmapping.dat.gz"
        for base in FTP_BASES
    ]


def idmapping_full_stream(
    id_types: set[str] | None = None,
    include_taxids: bool = True,
):
    """Stream ALL records from the complete idmapping.dat.gz.

    This file is ~18GB compressed. Lines are streamed from gzip.

    Args:
        id_types: If given, only yield rows for these ID type names.
            NCBI_TaxID rows are always included if include_taxids is True.
        include_taxids: Always yield NCBI_TaxID rows (for organism assignment).

    Yields:
        Tuples of (uniprot_ac, id_type_name, id_value).
    """
    urls = full_idmapping_urls()

    try:
        path = dm.download(urls[0], fallback_urls=urls[1:], connecttimeout=10)
    except Exception as e:
        _log.error("Full idmapping download failed: %s", e)
        return

    if not path or os.path.getsize(path) <= 1000:
        _log.error("All mirrors failed for full idmapping.dat.gz")
        return

    _log.info("Streaming full idmapping from %s", path)

    filter_types = set(id_types) if id_types else None
    if filter_types and include_taxids:
        filter_types.add("NCBI_TaxID")

    count = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue

            uniprot_ac, id_type_name, id_value = parts

            if filter_types and id_type_name not in filter_types:
                continue

            yield uniprot_ac, id_type_name, id_value
            count += 1

            if count % 10_000_000 == 0:
                _log.info("Streamed %dM records", count // 1_000_000)

    _log.info("Finished streaming: %d total records", count)


def idmapping_full_chunks(
    id_types: set[str] | None = None,
    n_chunks: int = 16,
    workdir: str | None = None,
) -> list[str]:
    """Download the full idmapping.dat.gz and split it into N parallel chunks.

    The complete ``idmapping.dat.gz`` (~18 GB, all organisms -- including the
    long tail that has no per-organism file) is downloaded once via the
    dlmachine download manager (cached, so re-runs do not re-fetch). It is then
    decompressed, filtered to the requested ID types (``NCBI_TaxID`` is always
    kept, for organism resolution) and split into ``n_chunks`` line-aligned,
    AC-contiguous chunk files with a single C pipeline
    (``pigz``/``zcat`` | ``awk`` | ``split -n l/N``).

    ``split -n l/N`` makes each chunk a contiguous byte range, so every UniProt
    AC's block stays inside a single chunk (apart from a handful straddling
    chunk boundaries). That lets a consumer resolve the organism inline per AC
    (see `parse_idmapping_chunk`) instead of a global join.

    Args:
        id_types: ID type names to keep, as they appear in the file (e.g.
            'Gene_Name', 'GeneID', 'Ensembl'). ``None`` keeps every type.
        n_chunks: Number of chunk files to produce (= the consumer's degree of
            parallelism).
        workdir: Directory for the large intermediate + chunk files. A
            temporary directory is created when ``None``. The caller owns
            cleanup of ``workdir``.

    Returns:
        Sorted list of chunk file paths (decompressed, filtered, TSV). Empty if
        the download failed.
    """

    import glob
    import shlex
    import shutil
    import subprocess
    import tempfile

    urls = full_idmapping_urls()
    try:
        path = dm.download(urls[0], fallback_urls=urls[1:], connecttimeout=10)
    except Exception as e:
        _log.error("Full idmapping download failed: %s", e)
        return []

    if not path or os.path.getsize(path) <= 1000:
        _log.error("All mirrors failed for full idmapping.dat.gz")
        return []

    workdir = workdir or tempfile.mkdtemp(prefix="idmapping_chunks_")
    os.makedirs(workdir, exist_ok=True)
    filtered = os.path.join(workdir, "filtered.dat")
    chunk_prefix = os.path.join(workdir, "chunk_")

    decomp = "pigz -dc" if shutil.which("pigz") else "zcat"
    if id_types:
        wanted = " ".join(sorted(set(id_types) | {"NCBI_TaxID"}))
        awk_prog = (
            'BEGIN{n=split(w,a," ");for(i=1;i<=n;i++)W[a[i]]=1} ($2 in W)'
        )
        cmd = (
            f"{decomp} {shlex.quote(path)} | "
            f"awk -F'\\t' -v w={shlex.quote(wanted)} {shlex.quote(awk_prog)} "
            f"> {shlex.quote(filtered)}"
        )
    else:
        cmd = f"{decomp} {shlex.quote(path)} > {shlex.quote(filtered)}"

    _log.info("Decompressing + filtering idmapping (%s)...", decomp)
    subprocess.run(["bash", "-c", cmd], check=True)

    _log.info("Splitting filtered idmapping into %d chunks...", n_chunks)
    subprocess.run(
        [
            "split", "-n", f"l/{n_chunks}", "-d", "-a", "4",
            filtered, chunk_prefix,
        ],
        check=True,
    )
    os.unlink(filtered)  # free space before the consumer reads the chunks

    chunks = sorted(glob.glob(chunk_prefix + "*"))
    _log.info("Produced %d chunk files in %s", len(chunks), workdir)

    return chunks


def parse_idmapping_chunk(
    path: str,
    id_types: set[str] | None = None,
):
    """Yield ``(uniprot_ac, ncbi_tax_id, id_type_name, id_value)`` from a chunk.

    The chunk must be AC-contiguous (as produced by `idmapping_full_chunks`):
    every line for one UniProt AC appears consecutively. The organism is
    resolved inline from each AC's ``NCBI_TaxID`` line, so no global join is
    needed. ``NCBI_TaxID`` rows are not yielded themselves; when ``id_types`` is
    given, only those types are yielded.
    """

    prev_ac = None
    group = []

    def rows_for_group():
        if not group:
            return []
        taxid = 0
        for tname, tval in group:
            if tname == "NCBI_TaxID":
                try:
                    taxid = int(tval)
                except ValueError:
                    taxid = 0
                break
        out = []
        for tname, tval in group:
            if tname == "NCBI_TaxID":
                continue
            if id_types and tname not in id_types:
                continue
            out.append((prev_ac, taxid, tname, tval))
        return out

    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            ac, tname, tval = parts
            if ac != prev_ac:
                for row in rows_for_group():
                    yield row
                group = []
                prev_ac = ac
            group.append((tname, tval))
        for row in rows_for_group():
            yield row


def stream_full_idmapping(block_size: int = 1 << 20, path: str | None = None):
    """Yield decompressed bytes blocks of the full idmapping.dat.

    Downloads the complete ``idmapping.dat.gz`` once via the dlmachine download
    manager (cached, all mirrors tried as fallbacks under one cache entry) and
    streams its **decompressed** bytes in ``block_size`` chunks — suitable for
    piping straight into a Postgres ``COPY ... FROM STDIN`` (the file is already
    tab-separated ``ac<TAB>id_type_label<TAB>id_value``, the default COPY text
    format). Decompression uses ``pigz``/``zcat`` when present (fast), falling
    back to Python ``gzip`` (portable). Only the caller's host needs a
    decompressor — nothing is required on the database server.

    Args:
        block_size: Bytes per yielded block.
        path: Stream this already-downloaded ``.gz`` directly, skipping the
            download manager (e.g. a pre-staged file).

    Yields:
        ``bytes`` blocks of the decompressed file.
    """

    import shlex
    import shutil
    import subprocess

    if not path:
        urls = full_idmapping_urls()
        try:
            path = dm.download(
                urls[0], fallback_urls=urls[1:], connecttimeout=10,
            )
        except Exception as e:
            _log.error("Full idmapping download failed: %s", e)
            return

    if not path or not os.path.exists(path) or os.path.getsize(path) <= 1000:
        _log.error("All mirrors failed for full idmapping.dat.gz")
        return

    decomp = (
        "pigz -dc" if shutil.which("pigz")
        else "zcat" if shutil.which("zcat")
        else None
    )

    if decomp:
        _log.info("Streaming idmapping via %s", decomp)
        proc = subprocess.Popen(
            ["bash", "-c", f"{decomp} {shlex.quote(path)}"],
            stdout=subprocess.PIPE,
            bufsize=block_size,
        )
        try:
            while True:
                block = proc.stdout.read(block_size)
                if not block:
                    break
                yield block
        finally:
            proc.stdout.close()
            if proc.wait() != 0:
                _log.warning("Decompressor exited with %s", proc.returncode)
    else:
        _log.info("Streaming idmapping via Python gzip (no zcat/pigz)")
        with gzip.open(path, "rb") as f:
            while True:
                block = f.read(block_size)
                if not block:
                    break
                yield block
