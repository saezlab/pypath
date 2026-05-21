"""Parser for PFOCR GMT exports from WikiPathways."""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.share.downloads import download_and_open


_CURRENT_PFOCR_URL = 'https://data.wikipathways.org/pfocr/current/'
_PFOCR_FILE_RE = re.compile(r'href="(pfocr-\d+-(chemical-)?gmt-([A-Za-z_]+)\.gmt)"')
_CURRENT_PFOCR_INDEX_FILENAME = 'current-index.html'
_CURRENT_PFOCR_HTML: str | None = None
_TAXONOMY_IDS = {
    'Caenorhabditis_elegans': '6239',
    'Danio_rerio': '7955',
    'Drosophila_melanogaster': '7227',
    'Homo_sapiens': '9606',
    'Mus_musculus': '10090',
    'Saccharomyces_cerevisiae': '4932',
}


def _current_pfocr_html() -> str:
    """Download the current PFOCR index through the shared download manager."""

    global _CURRENT_PFOCR_HTML

    if _CURRENT_PFOCR_HTML is not None:
        return _CURRENT_PFOCR_HTML

    opener = download_and_open(
        url=_CURRENT_PFOCR_URL,
        filename=_CURRENT_PFOCR_INDEX_FILENAME,
        subfolder='pfocr',
        large=False,
        encoding='utf-8',
        default_mode='r',
        force_download=True,
    )

    try:
        _CURRENT_PFOCR_HTML = opener.result
    finally:
        opener.close()

    return _CURRENT_PFOCR_HTML


def current_pfocr_filename(*, species: str = 'Homo_sapiens', data_type: str = 'gene', **_kwargs: Any) -> str:
    """Resolve the current PFOCR GMT filename for a species and data type."""
    chemical = data_type == 'chemical'
    html = _current_pfocr_html()

    matches = [
        filename
        for filename, chemical_part, file_species in _PFOCR_FILE_RE.findall(html)
        if file_species == species and bool(chemical_part) == chemical
    ]

    if not matches:
        raise RuntimeError(f'Could not resolve current PFOCR {data_type} GMT for {species}.')

    return sorted(matches)[-1]


def current_pfocr_url(**kwargs: Any) -> str:
    """Resolve the current PFOCR GMT URL."""
    return f'{_CURRENT_PFOCR_URL}{current_pfocr_filename(**kwargs)}'


def iter_gmt(
    opener,
    *,
    species: str = 'Homo_sapiens',
    data_type: str = 'gene',
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one PFOCR figure-to-entity association record per GMT member.

    GMT fields are ``PFOCR figure id``, ``figure/pathway caption`` and then one
    or more entity identifiers. Gene GMT identifiers are NCBI Gene IDs; chemical
    GMT identifiers are ChEBI accessions.
    """
    handle = _first_handle(opener)
    if handle is None:
        return

    for line in handle:
        if isinstance(line, bytes):
            line = line.decode('utf-8', 'ignore')
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 3:
            continue

        pfocr_id, caption, *identifiers = parts
        pfocr_id = pfocr_id.strip()
        if not pfocr_id:
            continue

        for identifier in identifiers:
            identifier = identifier.strip()
            if not identifier:
                continue
            yield {
                'pfocr_id': pfocr_id,
                'caption': caption.strip(),
                'identifier': identifier,
                'data_type': data_type,
                'species': species,
                'taxonomy_id': _TAXONOMY_IDS.get(species, ''),
            }


def _first_handle(opener):
    if not opener or not getattr(opener, 'result', None):
        return None
    return next(iter(opener.result.values()), None) if isinstance(opener.result, dict) else opener.result
