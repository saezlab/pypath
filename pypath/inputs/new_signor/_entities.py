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

import csv
from collections.abc import Generator

import pypath.resources.urls as urls
from pypath.share.downloads import dm
from ._records import SignorComplex, SignorProteinFamily, SignorPhenotype, SignorStimulus

__all__ = [
    'signor_complexes',
    'signor_protein_families',
    'signor_phenotypes',
    'signor_stimuli',
]


def signor_complexes(organism: int = 9606) -> Generator[SignorComplex]:
    """
    Download SIGNOR complex data.

    Args:
        organism: NCBI taxonomy ID (9606 for human, 10090 for mouse, 10116 for rat)

    Yields:
        SignorComplex named tuples
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
                    complex_id = row[0].strip()
                    name = row[1].strip(' "\n\r')
                    components = [c.strip(' "\n\r') for c in row[2].split(',')]

                    yield SignorComplex(
                        complex_id=complex_id,
                        name=name,
                        components=components
                    )


def signor_protein_families(organism: int = 9606) -> Generator[SignorProteinFamily]:
    """
    Download SIGNOR protein family data.

    Args:
        organism: NCBI taxonomy ID (9606 for human, 10090 for mouse, 10116 for rat)

    Yields:
        SignorProteinFamily named tuples
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
                    family_id = row[0].strip()
                    name = row[1].strip(' "\n\r') if len(row) > 1 else ''
                    # Members are in the 3rd column, comma-separated
                    members = [c.strip(' "\n\r') for c in row[2].split(',')]

                    yield SignorProteinFamily(
                        family_id=family_id,
                        name=name,
                        members=members
                    )


def signor_phenotypes(organism: int = 9606) -> Generator[SignorPhenotype]:
    """
    Download SIGNOR phenotype data.

    Args:
        organism: NCBI taxonomy ID (9606 for human, 10090 for mouse, 10116 for rat)

    Yields:
        SignorPhenotype named tuples
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
                    phenotype_id = row[0].strip()
                    name = row[1].strip(' "\n\r')
                    description = row[2].strip(' "\n\r')

                    yield SignorPhenotype(
                        phenotype_id=phenotype_id,
                        name=name,
                        description=description
                    )


def signor_stimuli(organism: int = 9606) -> Generator[SignorStimulus]:
    """
    Download SIGNOR stimulus data.

    Args:
        organism: NCBI taxonomy ID (9606 for human, 10090 for mouse, 10116 for rat)

    Yields:
        SignorStimulus named tuples
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
                    stimulus_id = row[0].strip()
                    name = row[1].strip(' "\n\r')
                    description = row[2].strip(' "\n\r')

                    yield SignorStimulus(
                        stimulus_id=stimulus_id,
                        name=name,
                        description=description
                    )
