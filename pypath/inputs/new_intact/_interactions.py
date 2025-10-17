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

"""
Download and process interaction data from IntAct database.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.share.downloads import download_and_open
from pypath.inputs import mitab
from pypath.inputs.mitab import MitabInteraction

__all__ = [
    'intact_interactions_raw',
    'intact_interactions',
]


def intact_interactions_raw(organism: int = 9606) -> Generator[str, None, None]:
    """
    Download IntAct interaction data in PSI-MITAB 2.7 format.

    Downloads and extracts zip file, yields lines for low memory usage.
    The human dataset is ~1 GB compressed, ~8 GB uncompressed.

    Args:
        organism: NCBI taxonomy ID (9606 for human)

    Yields:
        Lines from the MITAB file
    """
    if organism != 9606:
        raise ValueError(f'Currently only human (9606) is supported for IntAct')

    url = 'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip'

    # Download and open zip file
    opener = download_and_open(
        url=url,
        filename='human.zip',
        subfolder='intact',
        large=True,
        ext='zip',
    )

    # opener.result is a dict with filename -> file handle
    # Get the first file from the zip
    for file_handle in opener.result.values():
        yield from file_handle
        break  # Only process the first file


def intact_interactions(
    organism: int = 9606,
) -> Generator[MitabInteraction, None, None]:
    """
    Download and parse IntAct interactions in MITAB format.

    Args:
        organism: NCBI taxonomy ID (9606 for human)
        raw_records: If True, yield raw lines; if False, yield parsed MITAB records

    Yields:
        Interaction records (raw or parsed)
    """
    data = intact_interactions_raw(organism=organism)
    # Parse using mitab utilities
    yield from mitab.mitab_interactions(
        data=data,
        organism=organism,
        skip_header=True,
    )
