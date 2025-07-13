#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
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

from typing import Literal, NamedTuple

import collections
import csv
import io

import pypath.share.curl as curl
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping
import pypath.resources.urls as urls
import pypath.share.session as session

_log = session.Logger(name = 'gutmgene_input')._log


class GutmgeneRaw(NamedTuple):
    index: str
    pmid: str
    gut_microbiota: str
    gut_microbiota_ncbi_id: str
    rank: str
    strain: str
    gene: str
    gene_id: str
    alteration: str
    throughput: str
    associative_mode: str
    organism: str
    sample: str
    experimental_method: str
    measurement_technique: str
    description: str
    condition: str
    doid: str


class GutmgeneAnnotation(NamedTuple):
    gut_microbiota: str
    gut_microbiota_ncbi_id: str
    rank: str
    strain: str
    alteration: str
    throughput: str
    associative_mode: str
    sample: str
    experimental_method: str
    measurement_technique: str
    description: str
    condition: str
    doid: str
    pmid: str


def gutmgene_raw(organism: Literal['human', 'mouse'] = 'human') -> list[tuple]:
    """
    Gut microbiota genes from the gut microbiota gene database (gutMGene).

    Args:
        organism: Organism ID or name; human and mouse are available.

    Returns:
        A list of named tuples containing information about gut microbiota
        genes in human.
    """

    organism_ = taxonomy.ensure_common_name(organism, lower = True)

    if organism not in ('human', 'mouse'):

        err = f'`organism` must be either `human` or `mouse`, not `{organism}`.'
        _log(err)
        raise ValueError(err)

    # Use the rescued file since the original URLs are no longer accessible
    url = urls.urls['gutmgene']['url_rescued']
    c = curl.Curl(url, silent = False, large = True, encoding = 'iso-8859-1')

    result = set()

    # Parse CSV content
    csv_content = '\n'.join(c.result)
    csv_reader = csv.reader(io.StringIO(csv_content))

    # Skip header
    next(csv_reader, None)

    for row in csv_reader:

        # Skip if not enough columns
        if len(row) < 18:
            continue

        # Filter by organism (column 11 is human/mouse)
        if row[11].lower() != organism_:
            continue

        # Handle empty fields
        row = [None if not field.strip() else field.strip() for field in row]

        try:
            result.add(GutmgeneRaw(*row))
        except (TypeError, ValueError):
            # Skip malformed lines
            continue

    return list(result)


def gutmgene_annotations(
        organism: Literal['human', 'mouse'] = 'human',
    ) -> dict[str, set[GutmgeneAnnotation]]:
    """
    Microbial effectors of human or mouse genes from the gutMGene database.

    Args:
        organism:
            Organism ID or name; human and mouse are available.

    Return:
        A dict of sets of named tuples representing microbial relationships;
        top level keys are UniProt IDs.
    """

    raw = gutmgene_raw(organism)
    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)
    result = collections.defaultdict(set)

    for rec in raw:

        # Use gene symbol for mapping
        uniprots = mapping.map_name(
            rec.gene,
            'genesymbol',
            'uniprot',
            ncbi_tax_id = ncbi_tax_id,
        )

        for uniprot in uniprots:

            result[uniprot].add(
                GutmgeneAnnotation(
                    gut_microbiota = rec.gut_microbiota,
                    gut_microbiota_ncbi_id = rec.gut_microbiota_ncbi_id,
                    rank = rec.rank,
                    strain = rec.strain,
                    alteration = rec.alteration,
                    throughput = rec.throughput,
                    associative_mode = rec.associative_mode,
                    sample = rec.sample,
                    experimental_method = rec.experimental_method,
                    measurement_technique = rec.measurement_technique,
                    description = rec.description,
                    condition = rec.condition,
                    doid = rec.doid,
                    pmid = rec.pmid,
                )
            )

    return dict(result)
