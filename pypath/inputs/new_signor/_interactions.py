#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
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

from __future__ import annotations

from collections.abc import Generator

import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy
from pypath.share.downloads import dm

__all__ = [
    'signor_interactions_raw',
    'signor_interactions',
]


def signor_interactions_raw(organism: int = 9606) -> Generator[str]:
    """
    Download SIGNOR interaction data in causalTab (MITAB) format.

    Args:
        organism: NCBI taxonomy ID (9606 for human, 10090 for mouse, 10116 for rat)

    Yields:
        Lines from the causalTab file
    """
    if isinstance(organism, int):
        if organism in taxonomy.taxids:
            _organism = taxonomy.taxids[organism]
        else:
            raise ValueError(f'Unknown organism: {organism}')
    else:
        _organism = organism

    if _organism not in {'human', 'rat', 'mouse'}:
        raise ValueError(f'Organism {_organism} not supported by SIGNOR')

    url = urls.urls['signor']['all_url_new']

    # Download file with POST form data
    file_path = dm.download(
        url,
        filename=f'signor_{_organism}_causalTab.txt',
        subfolder='signor',
        query={
            'organism': _organism,
            'format': 'causalTab',
            'submit': 'Download',
        },
        post=True,
    )

    # Read and yield lines from the downloaded file
    if file_path:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:  # Skip empty lines
                    yield line


def signor_interactions(
    organism: int = 9606,
) -> Generator:
    """
    Download and parse SIGNOR interactions in causalTab format.

    Args:
        organism: NCBI taxonomy ID (9606 for human, 10090 for mouse, 10116 for rat)
        raw_records: If True, yield raw lines; if False, yield parsed MITAB records

    Yields:
        Interaction records (raw or parsed)
    """
    from .. import mitab

    data = signor_interactions_raw(organism=organism)

    yield from mitab.mitab_interactions(
        data=data,
        organism=organism,
        skip_header=True,
    )
