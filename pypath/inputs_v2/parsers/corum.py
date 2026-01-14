"""
CORUM (Comprehensive Resource of Mammalian protein complexes) parser.

Parses CORUM's TSV format and yields flat dictionaries.
"""

from __future__ import annotations

from collections.abc import Generator
import csv
from typing import Any


def _raw(opener, organism: int = 9606, **_kwargs: object) -> Generator[dict[str, Any], None, None]:
    """
    Parse CORUM complex data.

    Args:
        opener: File opener from download_and_open
        organism: NCBI taxonomy ID to filter by (default: 9606 for human)

    Yields:
        Dictionary for each complex record
    """
    if not opener or not opener.result:
        return

    # CORUM comes as a zip containing allComplexes.txt
    file_handle = None
    if isinstance(opener.result, dict):
        for name, handle in opener.result.items():
            if 'allComplexes' in name:
                file_handle = handle
                break
    else:
        file_handle = opener.result

    if not file_handle:
        return

    reader = csv.DictReader(file_handle, delimiter='\t')

    for rec in reader:
        # Filter by organism if specified
        cplex_organism = rec.get('Organism', '')
        # CORUM uses organism names, need to check the name or ID
        if organism == 9606 and 'Human' not in cplex_organism and 'sapiens' not in cplex_organism.lower():
            continue
        elif organism == 10090 and 'Mouse' not in cplex_organism and 'musculus' not in cplex_organism.lower():
            continue
        elif organism == 10116 and 'Rat' not in cplex_organism and 'norvegicus' not in cplex_organism.lower():
            continue

        yield rec
