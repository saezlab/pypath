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

import csv
import itertools
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera
import pypath.share.session as session
import pypath.core.entity as entity

_logger = session.Logger(name = 'cellinker_input')
_log = _logger._log


CellinkerInteraction = collections.namedtuple(
    'CellinkerInteraction',
    (
        'ligand',
        'receptor',
        'ligand_location',
        'receptor_location',
        'resources',
        'pmids',
        'type',
    ),
)


def cellinker_complexes_raw(organism = 9606):
    """
    Downloads protein complex data from the Cellinker database
    (http://www.rna-society.org/cellinker/).

    Args
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.

    Returns
        (list): List of tuples each describing a protein complex with its
            role (ligand or receptor), components, localization and
            Cellinker ID.
    """

    CellinkerComplex = collections.namedtuple(
        'CellinkerComplex',
        (
            'role',
            'cellinker_id',
            'components',
            'location',
        ),
    )

    CellinkerComplexComponent = collections.namedtuple(
        'CellinkerComplexComponent',
        (
            'genesymbol',
            'entrez',
        ),
    )

    organism_common = taxonomy.ensure_common_name(organism).lower()

    if organism_common not in {'human', 'mouse'}:

        msg = (
            'Unknown organism: %s (%s). Only human and mouse '
            'are available.' % (str(organism_common), str(organism))
        )
        _log(msg)
        raise ValueError(msg)

    url = urls.urls['cellinker_rescued']['complex'] % organism_common
    c = curl.Curl(url, large = True, silent = False)

    result = []

    _ = next(c.result)

    for r in c.result:

        r = r.split(',')

        components = tuple(
            CellinkerComplexComponent(
                genesymbol = r[i],
                entrez = r[i + 1],
            )
            for i in range(3, 13, 2)
            if r[i]
        )

        result.append(
            CellinkerComplex(
                role = 'ligand' if r[0] else 'receptor',
                cellinker_id = r[0] or r[1],
                components = components,
                location = r[13],
            )
        )

    return result


def components_to_complex(components, organism = None):
    """
    Converts a set of components to `pypath.internals.intera.Complex`
    objects.

    Args
        components (tuple): Components of a complex, as returned by
            `cellinker_complexes_raw`.
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available. Optional, because organism can be
            guessed from the identifiers.

    Returns
        (set): A set of `pypath.internals.intera.Complex` objects.
    """

    if not organism:

        organism = (
            9606 if (
                all(c.genesymbol.upper() == c.genesymbol for c in components)
            ) else 10090
        )

    _organism = taxonomy.ensure_ncbi_tax_id(organism)

    if _organism not in {9606, 10090}:

        msg = (
            'Unknown organism: %s (%s). Only human and mouse '
            'are available.' % (str(_organism), str(organism))
        )
        _log(msg)
        raise ValueError(msg)

    result = set()

    for uniprots in itertools.product(*(
        _cellinker_uniprots(c.genesymbol, c.entrez, _organism)
        for c in components
    )):

        result.add(
            intera.Complex(
                components = uniprots,
                ncbi_tax_id = _organism,
                sources = 'Cellinker',
            )
        )

    return result


def cellinker_complexes(organism = 9606):
    """
    Protein complex information from the Cellinker database
    (http://www.rna-society.org/cellinker/).

    Args
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.

    Returns
        (dict): A dict of complexes, with string representations as keys and
            `pypath.internals.intera.Complex` objects as values.
    """

    result = {}

    for c in cellinker_complexes_raw(organism = organism):

        for cplex in components_to_complex(c.components, organism = organism):

            result[cplex.__str__()] = cplex

    return result


def cellinker_lr_interactions_raw(organism = 9606):
    """
    Ligand-receptor interactions from the Cellinker database
    (http://www.rna-society.org/cellinker/).

    Args
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.

    Returns
        (list): A list of dicts, each representing an interaction records
            as it is provided by the database.
    """

    return _cellinker_interactions_raw(organism = organism)


def cellinker_smol_interactions_raw(organism = 9606):
    """
    Small molecule ligand-protein receptor interactions from the Cellinker
    database (http://www.rna-society.org/cellinker/).

    Args
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.

    Returns
        (list): A list of dicts, each representing an interaction records
            as it is provided by the database.
    """

    return _cellinker_interactions_raw(dataset = 'smol', organism = organism)


def _cellinker_interactions_raw(dataset = 'lr', organism = 9606):
    """
    Downloads either the ligand-receptor or the small molecule ligand-receptor
    dataset from the Cellinker database.

    Args
        dataset (str): Either `lr` or `smol`, meainng protein-protein or
            small molecule-protein ligand-receptor interactions.
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.

    Returns
        (list): A list of dicts, each representing an interaction records
            as it is provided by the database.
    """

    if dataset not in {'lr', 'smol'}:

        msg = 'Unknown Cellinker interaction dataset: `%s`.' % str(dataset)
        _log(msg)
        raise ValueError(msg)

    org_name_type = 'latin' if dataset == 'lr' else 'common'

    organisms_allowed = (
        {'Homo sapiens', 'Mus musculus'}
            if org_name_type == 'latin' else
        {'Human', 'Mouse'}
    )

    _organism = getattr(taxonomy, 'ensure_%s_name' % org_name_type)(organism)

    if _organism not in organisms_allowed:

        msg = (
            'Unknown organism: %s (%s). Only human and mouse '
            'are available.' % (str(_organism), str(organism))
        )
        _log(msg)
        raise ValueError(msg)

    if org_name_type == 'common': _organism = _organism.lower()

    url = urls.urls['cellinker_rescued'][dataset] % _organism
    c = curl.Curl(url, large = True, silent = False)

    result = list(csv.DictReader(c.result, delimiter = '\t'))

    return result


def cellinker_lr_interactions(organism = 9606):
    """
    Ligand-receptor interactions from the Cellinker database
    (http://www.rna-society.org/cellinker/).

    Args
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.

    Returns
        (set): A set of tuples, each representing a preprocessed Cellinker
            interaction. The proteins are represented by their UniProt IDs,
            while the protein complexes by `Complex` objects.
    """

    db_names = {
        'IUPHAR': 'Guide2Pharma',
        'CellphoneDB': 'CellPhoneDB',
    }

    result = set()

    raw = cellinker_lr_interactions_raw(organism = organism)
    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)
    complexes = dict(
        (
            c.cellinker_id,
            components_to_complex(c.components, organism = ncbi_tax_id)
        )
        for c in cellinker_complexes_raw(organism = organism)
    )

    for r in raw:

        ligands = _cellinker_uniprots(
            r['Ligand_symbol'],
            r['Ligand_id'],
            ncbi_tax_id,
            complexes = complexes,
        )
        receptors = _cellinker_uniprots(
            r['Receptor_symbol'],
            r['Receptor_id'],
            ncbi_tax_id,
            complexes = complexes,
        )

        resources = ';'.join(
            db_names.get(db, db)
            for db in r['Other.DB'].split(';')
        ) or None

        for ligand, receptor in itertools.product(ligands, receptors):

            result.add(
                CellinkerInteraction(
                    ligand = ligand,
                    receptor = receptor,
                    ligand_location = r['Ligand_location'],
                    # yes, labels are not consistent
                    receptor_location = r['Receptor.location'],
                    resources = resources,
                    pmids = r['Pmubmed.ID'] or None, # typo
                    type = r['Type'],
                )
            )

    return result


def cellinker_smol_interactions(organism = 9606):
    """
    Small molecule ligand-protein receptor interactions from the Cellinker
    database (http://www.rna-society.org/cellinker/).

    Args
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.

    Returns
        (set): A set of tuples, each representing a preprocessed Cellinker
            interaction. The proteins are represented by their UniProt IDs,
            the small molecules by PubChem CIDs, while the protein complexes
            by `Complex` objects.
    """

    db_names = {
        'IUPHAR': 'Guide2Pharma',
        'CellphoneDB': 'CellPhoneDB',
    }

    result = set()

    raw = cellinker_smol_interactions_raw(organism = organism)
    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)
    complexes = dict(
        (
            c.cellinker_id,
            components_to_complex(c.components, organism = ncbi_tax_id)
        )
        for c in cellinker_complexes_raw(organism = organism)
    )

    for r in raw:

        if not r['ligand_pubchem_cid']:

            continue

        ligands = (r['ligand_pubchem_cid'],)

        receptors = _cellinker_uniprots(
            r['Receptor_symbol'],
            r['Receptor_id'],
            ncbi_tax_id,
            complexes = complexes,
        )

        resources = ';'.join(
            db_names.get(db, db)
            for db in r['Other.DB'].split(';')
        ) or None

        for ligand, receptor in itertools.product(ligands, receptors):

            result.add(
                CellinkerInteraction(
                    ligand = ligand,
                    receptor = receptor,
                    ligand_location = None,
                    receptor_location = r['Receptor_location'],
                    resources = resources,
                    pmids = r['pubmed_id'] or None,
                    type = r['Type'],
                )
            )

    return result


def cellinker_annotations(organism = 9606, entity_type = None):
    """
    Ligand and receptor annotations from the Cellinker database
    (http://www.rna-society.org/cellinker/).

    Args
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.
        entity_type (str): Either `protein` or `complex`. If `None`, both
            proteins and protein complexes will be included.

    Returns
        (dict): A dict of sets of tuples, keys are UniProt IDs for proteins
            and `Complex` objects for protein complexes. The tuples are
            annotations with ligand or receptor role, localization and
            type.
    """

    CellinkerAnnotation = collections.namedtuple(
        'CellinkerAnnotation',
        (
            'role',
            'location',
            'type',
        ),
    )

    ia = cellinker_lr_interactions(organism = organism)

    result = collections.defaultdict(set)

    for i in ia:

        for role in ('ligand', 'receptor'):

            this_entity = getattr(i, role)
            this_entity_type = entity.Entity._get_entity_type(this_entity)

            if not entity_type or entity_type == this_entity_type:

                result[this_entity].add(
                    CellinkerAnnotation(
                        role = role,
                        location = getattr(i, '%s_location' % role),
                        type = i.type,
                    )
                )

    return dict(result)


def cellinker_protein_annotations(organism = 9606):
    """
    Ligand and receptor annotations from the Cellinker database
    (http://www.rna-society.org/cellinker/).

    Args
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.

    Returns
        (dict): A dict of sets of tuples, keys are UniProt IDs. The tuples
            are annotations with ligand or receptor role, localization and
            type.
    """

    return cellinker_annotations(organism = organism, entity_type = 'protein')


def cellinker_complex_annotations(organism = 9606):
    """
    Ligand and receptor annotations from the Cellinker database
    (http://www.rna-society.org/cellinker/).

    Args
        organism (int,str): Name or identifier of the organism. Only mouse
            and human are available.

    Returns
        (dict): A dict of sets of tuples, keys are `Complex` objects.
            The tuples are annotations with ligand or receptor role,
            localization and type.
    """

    return cellinker_annotations(organism = organism, entity_type = 'complex')


def _cellinker_uniprots(gsymbol, entrez, ncbi_tax_id, complexes = None):
    """
    Translates the Gene Symbols and Entrez Gene IDs to UniProt IDs.

    Returns
        (set): Set of UniProt IDs.
    """

    return (
        complexes[entrez]
            if complexes and entrez in complexes else
        (
            mapping.map_name(
                gsymbol,
                'genesymbol',
                'uniprot',
                ncbi_tax_id = ncbi_tax_id,
            ) |
            mapping.map_name(
                entrez,
                'entrez',
                'uniprot',
                ncbi_tax_id = ncbi_tax_id,
            )
        )
    )
