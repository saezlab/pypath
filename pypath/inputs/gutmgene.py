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

import pypath.share.curl as curl
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping
import pypath.resources.urls as urls
import pypath.share.session as session

_log = session.Logger(name = 'gutmgene_input')._log


class GutmgeneRaw(NamedTuple):
    microbe_taxon: str
    microbe_ncbi_tax_id: str
    gut_microbiota_id: str
    classification: str
    genesymbol: str
    entrez: str
    effect: str
    throughput: str


class GutmgeneAnnotation(NamedTuple):
    microbe_taxon: str
    microbe_ncbi_tax_id: str
    gut_microbiota_id: str
    classification: str
    effect: str
    throughput: str


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

        err = '`organism` must be either `human` or `mouse`, not `{organism}`.'
        _log(err)
        raise ValueError(err)

    url = urls.urls['gutmgene'][f'url_{organism_}']
    c = curl.Curl(url, silent = False, large = True)

    result = set()

    for l in c.result:

        if l.startswith('"'):
            continue

        l = l.replace('"', '')
        l = l.strip().split('\t')
        l = (None if not i else i for i in l)

        if l:

            result.add(GutmgeneRaw(*l))

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

        uniprots = mapping.map_name(
            rec.genesymbol,
            'genesymbol',
            'uniprot',
            ncbi_tax_id = ncbi_tax_id,
        )

        for uniprot in uniprots:

            result[uniprot].add(
                GutmgeneAnnotation(
                    microbe_taxon = rec.microbe_taxon,
                    microbe_ncbi_tax_id = rec.microbe_ncbi_tax_id,
                    gut_microbiota_id = rec.gut_microbiota_id,
                    classification = rec.classification,
                    effect = rec.effect,
                    throughput = rec.throughput,
                )
            )

    return dict(result)
