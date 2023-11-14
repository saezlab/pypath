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

import re
import csv
import collections
import itertools

import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy
import pypath.internals.intera as intera
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.share.session as session
import pypath.core.entity as entity

_logger = session.Logger(name = 'scconnect_input')
_log = _logger._log


def scconnect_annotations(organism = 9606):
    """
    Ligand and receptor annotations from scConnect
    (https://github.com/JonETJakobsson/scConnect).

    Args
        organism (int,str): Name or identifier of the organism. Human, mouse,
            fruitfly, zebrafish, chicken, C. elegans, Xenopus tropicalis,
            yeast and Anolis carolinensis are available.
    """

    _organism = taxonomy.ensure_ensembl_name(organism)
    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

    if not _organism:

        msg = 'Could not recognize organism: `%s`.' % str(organism)
        _log(msg)
        raise ValueError(msg)

    record = collections.namedtuple(
        'ScconnectAnnotation',
        (
            'role',
            'family',
            'type',
            'inferred_from',
        ),
    )

    result = collections.defaultdict(set)
    reinf = re.compile('inferred from (\w+)')

    for role in ('ligand', 'receptor'):

        url = urls.urls['scconnect']['annot'] % (_organism, role)
        c = curl.Curl(url, silent = False, large = True)
        tab = csv.DictReader(c.result)

        for rec in tab:

            typ = rec.get('type', None)
            family = rec.get('family', None)
            inferred_from = rec.get('comment', None)

            if inferred_from:

                inferred_from = reinf.match(inferred_from)

                if inferred_from:

                    inferred_from = inferred_from.groups()[0].lower()

            annot = record(
                role = role,
                family = family,
                type = typ,
                inferred_from = inferred_from,
            )

            genesymbols = (
                rec.get('gene', rec.get('preprogene')).
                strip('[\']').
                split(',')
            )

            for gs in genesymbols:

                components = gs.split('|')

                uniprots = [
                    mapping.map_name(
                        gs_comp.strip(),
                        'genesymbol',
                        'uniprot',
                        ncbi_tax_id = ncbi_tax_id,
                    )
                    for gs_comp in components
                ]

                for ups in itertools.product(*uniprots):

                    if len(ups) > 1:

                        cplex = intera.Complex(
                            components = ups,
                            sources = 'scConnect',
                        )

                        result[cplex].add(annot)

                    else:

                        result[ups[0]].add(annot)

    return dict(result)


def scconnect_complexes(organism = 9606):
    """
    Protein complexes from scConnect
    (https://github.com/JonETJakobsson/scConnect).

    Args
        organism (int,str): Name or identifier of the organism. Human, mouse,
            fruitfly, zebrafish, chicken, C. elegans, Xenopus tropicalis,
            yeast and Anolis carolinensis are available.
    """

    annot = scconnect_annotations(organism = organism)

    return {
        cplex
        for cplex in annot.keys()
        if entity.Entity._is_complex(cplex)
    }


def scconnect_interactions():
    """
    Ligand-receptor interactions from scConnect
    (https://github.com/JonETJakobsson/scConnect).

    Returns
        (list): List of interactions, each represented as a named tuple.
            Proteins and protein complexes are translated to UniProt IDs,
            small molecule IDs are left intact.
    """


    def process_partner(rec, partner):

        organisms = [
            taxonomy.ensure_ncbi_tax_id(org)
                if org not in {'', 'None', 'Unknown'} else
            _const.NOT_ORGANISM_SPECIFIC
            for org in rec['%s_species' % partner].split('|')
        ]
        id_field = (
            'target_uniprot'
                if partner == 'target' else
            'ligand_gene_symbol'
        )
        id_type = 'uniprot' if partner == 'target' else 'genesymbol'

        for organism in set(organisms):

            ids_raw = [
                _id
                for _id, _org in
                zip(
                    rec[id_field].split('|'),
                    organisms
                )
                if _org == organism or organism is None
            ]

            ids = (
                [
                    mapping.map_name(
                        _id,
                        id_type,
                        'uniprot',
                        ncbi_tax_id = organism,
                    )
                    for _id in ids_raw
                ]
                    if organism else
                [(rec[partner],)]
            )

            ids = [
                (
                    (
                        intera.Complex(
                            components = _ids,
                            sources = 'scConnect',
                        ),
                        'complex',
                    )
                        if len(_ids) > 1 else
                    (
                        _ids[0],
                        'protein' if organism else 'small_molecule',
                    )
                )
                for _ids in itertools.product(*ids)
            ]

            for _id, entity_type in ids:

                yield _id, entity_type, organism


    url = urls.urls['scconnect']['intera']
    c = curl.Curl(url, silent = False, large = True)
    tab = csv.DictReader(c.result)

    record = collections.namedtuple(
        'ScconnectInteraction',
        (
            'ligand_id',
            'target_id',
            'ligand_organism',
            'target_organism',
            'ligand_type',
            'target_type',
            'effect',
            'references',
        ),
    )

    result = []

    for rec in tab:

        targets = process_partner(rec, 'target')
        ligands = process_partner(rec, 'ligand')

        for ligand_target in itertools.product(ligands, targets):

            (
                (ligand, ligand_type, ligand_organism),
                (target, target_type, target_organism)
            ) = ligand_target

            result.append(
                record(
                    ligand_id = ligand,
                    target_id = target,
                    ligand_organism = ligand_organism,
                    target_organism = target_organism,
                    ligand_type = ligand_type,
                    target_type = target_type,
                    effect = rec['action'],
                    references = rec['pubmed_id'],
                )
            )

    return result
