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

from collections.abc import Generator

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.inputs.common as inputs_common
import pypath.utils.mapping as mapping


def _slc_raw(table = 1, sheet = 1) -> list[list]:

    url = urls.urls['slc'][f'table_s{table}']
    c = curl.Curl(url, large = True, silent = False)
    content = inputs_common.read_xls(c.fileobj.name, sheet = sheet)

    return content


def _to_namedtuples(table: list[list], name: str) -> list[tuple]:

    record = collections.namedtuple(f'Slc{name}', table[0])

    return [record(*line) for line in table[1:]]


def slc_annotation() -> list[list]:

    return _to_namedtuples(
        _slc_raw(table = 1, sheet = 1),
        'Annotation',
    )


def slc_substrate_ontology() -> list[list]:

    return _to_namedtuples(
        _slc_raw(table = 1, sheet = 2),
        'SubstrateOntology',
    )


def slc_chebi_mapping() -> list[tuple]:

    return _to_namedtuples(
        [
            r[:-7]
            for r in _slc_raw(table = 2, sheet = 1)
        ],
        'ChebiMapping',
    )


def _pmids(pmids: str) -> str | None:

    return ';'.join(sorted(
        pp
        for p in str(pmids).split(',')
        if (pp := p.strip())
    )) or None


def slc_localization_annotations() -> dict[str, set[tuple]]:

    SlcLocalization = collections.namedtuple(
        'SlcLocalization',
        [
            'localization',
            'pmids',
        ],
    )

    result = collections.defaultdict(set)

    for raw in slc_annotation():

        uniprots = mapping.map_name(
            raw.HGNC_symbol,
            'genesymbol',
            'uniprot',
        )

        for uniprot in uniprots:

            result[uniprot].add(
                SlcLocalization(
                    localization = raw.Subcellular_localization,
                    pmids = _pmids(raw.Subcellular_localization_PMID),
                )
            )

    return dict(result)


def slc_interactions() -> Generator[tuple, None, None]:

    SlcInteraction = collections.namedtuple(
        'SlcInteraction',
        [
            'transporter',
            'substrate',
            'role',
            'pmids',
            'localization',
            'transport_mechanism',
            'substrate_class',
        ],
    )

    SlcSubstrate = collections.namedtuple(
        'SlcSubstrate',
        [
            'slc_name',
            'chebi',
            'label',
            'synonyms',
        ],
    )

    SlcTransporter = collections.namedtuple(
        'SlcTransporter',
        [
            'uniprot',
            'genesymbol',
            'entrez',
            'ensg',
        ]
    )

    substrates = {}

    for sub in slc_chebi_mapping():

        substrates[sub.Substrate_name_annotation] = SlcSubstrate(
            slc_name = sub.Substrate_name_annotation,
            chebi = sub.chebi_id,
            label = sub.chebi_term,
            synonyms = sorted(sub.synonym.split(', ')),
        )

    transporters = {}
    NOTHING = {'', 'Unknown', 'None'}

    for raw in slc_annotation():

        if (genesymbol := raw.HGNC_symbol) not in transporters:

            uniprots = mapping.map_name(
                genesymbol,
                'genesymbol',
                'uniprot',
            )

            transporters[genesymbol] = {
                SlcTransporter(
                    uniprot = uniprot,
                    genesymbol = genesymbol,
                    entrez = raw.Entrez_ID,
                    ensg = raw.Ensembl_ID,
                )
                for uniprot in uniprots
            }

        for role, i in (('substrate', 4), ('coupled_ion', 6)):

            for partner in raw[i].split('; '):

                if (
                    partner in NOTHING or
                    not (substrate := substrates.get(partner))
                ):

                    continue

                for transporter in transporters[genesymbol]:

                    yield SlcInteraction(
                        transporter = transporter,
                        substrate = substrate,
                        role = role,
                        pmids = _pmids(raw[i + 1]),
                        localization = raw.Subcellular_localization,
                        transport_mechanism = raw.Transport_mechanism,
                        substrate_class = raw.Substrate_class,
                    )
