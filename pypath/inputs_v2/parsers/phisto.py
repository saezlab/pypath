"""
Parser for PHISTO pathogen–host protein interactions.

One row is an experimentally reported interaction between a pathogen protein
and a human protein, with the pathogen taxon, detection method and PubMed ID.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv

COLUMNS = (
    'Pathogen',
    'Taxonomy ID',
    'Uniprot ID [pathogen]',
    'Pathogen Protein',
    'Uniprot ID [human]',
    'Human Protein',
    'Experimental Method',
    'Pubmed ID',
)

def iter_interactions(
    opener, **kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield validated PHISTO CSV records."""
    yield from iter_csv(
        opener,
        strict=True,
        required=COLUMNS,
        **kwargs,
    )
# this will give a list of dictionaries where each dict corresponds to a row
