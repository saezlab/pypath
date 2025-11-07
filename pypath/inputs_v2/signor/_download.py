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
Raw data download functions for SIGNOR database.

This module handles downloading data files from SIGNOR without parsing.
"""

from __future__ import annotations

import csv
from collections.abc import Generator

import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy
from pypath.share.downloads import dm

__all__ = [
    'download_complexes',
    'download_protein_families',
    'download_phenotypes',
    'download_stimuli',
    'download_interactions',
]


def download_complexes() -> Generator[dict]:
    """
    Download SIGNOR complex data.

    Yields:
        Dictionaries with keys: complex_id, name, components (list)
    """
    url = urls.urls['signor']['complexes']

    file_path = dm.download(
        url,
        filename='signor_complexes.txt',
        subfolder='signor',
        query={'submit': 'Download complex data'},
        post=True,
    )

    if file_path:
        with open(file_path, 'r') as f:
            reader = csv.reader(f, delimiter=';')
            # Skip header
            next(reader, None)

            for row in reader:
                if len(row) >= 3:
                    yield {
                        'complex_id': row[0].strip(),
                        'name': row[1].strip(' "\n\r'),
                        'components': [c.strip(' "\n\r') for c in row[2].split(',')],
                    }


def download_protein_families() -> Generator[dict]:
    """
    Download SIGNOR protein family data.

    Yields:
        Dictionaries with keys: family_id, name, members (list)
    """
    url = urls.urls['signor']['complexes']

    file_path = dm.download(
        url,
        filename='signor_protein_families.txt',
        subfolder='signor',
        query={'submit': 'Download protein family data'},
        post=True,
    )

    if file_path:
        with open(file_path, 'r') as f:
            reader = csv.reader(f, delimiter=';')
            # Skip header
            next(reader, None)

            for row in reader:
                if len(row) >= 3:
                    yield {
                        'family_id': row[0].strip(),
                        'name': row[1].strip(' "\n\r') if len(row) > 1 else '',
                        'members': [m.strip(' "\n\r') for m in row[2].split(',')],
                    }


def download_phenotypes() -> Generator[dict]:
    """
    Download SIGNOR phenotype data.

    Yields:
        Dictionaries with keys: phenotype_id, name, description
    """
    url = urls.urls['signor']['complexes']

    file_path = dm.download(
        url,
        filename='signor_phenotypes.txt',
        subfolder='signor',
        query={'submit': 'Download phenotype data'},
        post=True,
    )

    if file_path:
        with open(file_path, 'r') as f:
            reader = csv.reader(f, delimiter=';')
            # Skip header
            next(reader, None)

            for row in reader:
                if len(row) >= 3:
                    yield {
                        'phenotype_id': row[0].strip(),
                        'name': row[1].strip(' "\n\r'),
                        'description': row[2].strip(' "\n\r'),
                    }


def download_stimuli() -> Generator[dict]:
    """
    Download SIGNOR stimulus data.

    Yields:
        Dictionaries with keys: stimulus_id, name, description
    """
    url = urls.urls['signor']['complexes']

    file_path = dm.download(
        url,
        filename='signor_stimuli.txt',
        subfolder='signor',
        query={'submit': 'Download stimulus data'},
        post=True,
    )

    if file_path:
        with open(file_path, 'r') as f:
            reader = csv.reader(f, delimiter=';')
            # Skip header
            next(reader, None)

            for row in reader:
                if len(row) >= 3:
                    yield {
                        'stimulus_id': row[0].strip(),
                        'name': row[1].strip(' "\n\r'),
                        'description': row[2].strip(' "\n\r'),
                    }


def download_interactions(organism: int = 9606) -> Generator[str]:
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
