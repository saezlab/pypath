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

from future.utils import iteritems
from past.builtins import xrange, range

import re
import csv
import collections
import itertools

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera
import pypath.share.common as common
import pypath.share.session as session
import pypath.inputs.uniprot_db as uniprot_db

_logger = session.Logger(name = 'cellphonedb_input')
_log = _logger._log

CellPhoneDBAnnotation = collections.namedtuple(
    'CellPhoneDBAnnotation',
    (
        'receptor',
        'receptor_class',
        'peripheral',
        'secreted',
        'secreted_class',
        'transmembrane',
        'integrin',
    )
)


def cellphonedb_ligands_receptors():
    """
    Retrieves the set of ligands and receptors from CellPhoneDB.
    Returns tuple of sets.
    """

    receptors = set()
    ligands   = set()

    proteins = cellphonedb_protein_annotations()
    complexes = cellphonedb_complex_annotations()

    for _id, annot in itertools.chain(
        iteritems(proteins),
        iteritems(complexes)
    ):
        if annot.receptor:
            receptors.add(_id)

        if annot.secreted or (
            not annot.receptor and (
                annot.transmembrane or
                annot.peripheral
            )
        ):
            ligands.add(_id)

    return ligands, receptors


def _cellphonedb_annotations(url, name_method):

    replacements = {
        '_add': '',
        ' | ': ',',
        ' ': '_',
    }

    def get_bool(rec, attr):

        return attr in rec and rec[attr].upper() == 'TRUE'

    def get_desc(rec, attr):

        desc = '%s_desc' % attr

        value = (
            ''
                if (
                    attr in rec and rec[attr].upper() == 'FALSE' or
                    attr not in rec and not rec[desc]
                ) else
            rec[desc]
                if rec[desc] else
            attr
        )

        for pattern, repl in iteritems(replacements):

            value = value.replace(pattern, repl)

        value = value.lower().split(',') if value else None

        return tuple(sorted(common.to_set(value)))

    record = CellPhoneDBAnnotation

    annot = {}

    c = curl.Curl(url, large = True, silent = False)
    tab = list(csv.DictReader(c.result))

    for rec in tab:

        names = name_method(rec)

        if isinstance(names, (str, intera.Complex)):

            names = (names,)

        for name in names:

            annot[name] = record(
                receptor = get_bool(rec, 'receptor'),
                receptor_class = get_desc(rec, 'receptor'),
                peripheral = get_bool(rec, 'peripheral'),
                secreted = get_bool(rec, 'secreted'),
                secreted_class = get_desc(rec, 'secreted'),
                transmembrane = get_bool(rec, 'transmembrane'),
                integrin = get_bool(rec, 'integrin'),
            )

    return annot


def cellphonedb_protein_annotations(add_complex_annotations = False):
    """
    :arg bool add_complex_annotations:
        Deprecated because results wrong annotations.
        Copy the annotations of complexes to each of their member proteins.
    """

    def name_method(rec):
        uniprot = rec['uniprot']
        uniprot = _cellphonedb_hla(uniprot)
        uniprot = mapping.map_names(uniprot, 'uniprot', 'uniprot')

        return uniprot

    protein_annotations = _cellphonedb_annotations(
        url = urls.urls['cellphonedb_git']['proteins'],
        name_method = name_method,
    )

    return protein_annotations


def _cellphonedb_hla(uniprot):
    """
    Returns *set*.
    """

    uniprots = None

    # for HLA genes in the uniprot column we have "HLA..." gene symbols
    # but not always in the standard format (with dash)
    if uniprot.startswith('HLA') and '-' not in uniprot:
        genesymbol = 'HLA-%s' % uniprot[3:]
        uniprots = mapping.map_name(genesymbol, 'genesymbol', 'uniprot')

    return uniprots or {uniprot}


def cellphonedb_complex_annotations():

    def get_uniprots(rec):

        return tuple(
            uniprot
            for uniprot in
            (rec['uniprot_%u' % i] for i in xrange(1, 5))
            if uniprot
        )

    def get_stoichiometry(rec):

        if not rec['stoichiometry']:
            return get_uniprots(rec)

        return tuple(
            mapping.map_name0(genesymbol, 'genesymbol', 'uniprot')
            for genesymbol in rec['stoichiometry'].split(';')
        )

    def name_method(rec):
        comp = get_stoichiometry(rec)

        cplex = intera.Complex(
            name = rec['complex_name'],
            components = comp,
            sources = 'CellPhoneDB',
            ids = rec['complex_name'],
        )

        return cplex

    return _cellphonedb_annotations(
        url = urls.urls['cellphonedb_git']['complexes'],
        name_method = name_method,
    )


def _cellphonedb_get_entity(name, complexes):

    if name in complexes:
        return (complexes[name],)

    if '_by' in name:
        _log(f'Ignoring entity: `{name}`.')
        return ()

    if ':' in name:
        name = name.split(':')[1]

    if '_' in name:
        name = mapping.map_name0(name, 'name-entry', 'name')

    if not uniprot_db.is_uniprot(name):
        uniprot = mapping.map_name0(name, 'genesymbol', 'uniprot')
        name = uniprot or name

    name = _cellphonedb_hla(name)

    return (name,) if isinstance(name, str) else name


def cellphonedb_interactions():
    """
    Interactions between ligands and receptors from CellPhoneDB.

    Yields:
        Named tuples representing interactions.
    """

    def get_type(entity):

        return (
            'ligand'
                if entity in ligands else
            'receptor'
                if entity in receptors else
            'unknown'
        )

    def get_bool(rec, attr):

        return attr in rec and rec[attr].upper() == 'TRUE'

    CellphonedbInteraction = collections.namedtuple(
        'CellphonedbInteraction',
        [
            'id_a',
            'id_b',
            'sources',
            'references',
            'interaction_type',
            'type_a',
            'type_b',
            'is_ppi',
        ]
    )

    repmid = re.compile(r'PMID: ([0-9]+)')
    recomma = re.compile(r'[,;]')

    ligands, receptors = cellphonedb_ligands_receptors()
    complexes = dict(
        (
            _id,
            cplex
        )
        for cplex in cellphonedb_complexes().values()
        for _id in cplex.ids['CellPhoneDB']
    )

    url = urls.urls['cellphonedb_git']['interactions']
    c = curl.Curl(url, silent = False, large = True)
    reader = csv.DictReader(c.result)

    for rec in reader:

        _partner_a = _cellphonedb_get_entity(
            rec['partner_a'],
            complexes = complexes,
        )
        _partner_b = _cellphonedb_get_entity(
            rec['partner_b'],
            complexes = complexes,
        )
        _is_ppi = get_bool(rec, 'is_ppi')

        for partner_a, partner_b in itertools.product(_partner_a, _partner_b):

            type_a = get_type(partner_a)
            type_b = get_type(partner_b)
            rev = type_b == 'ligand' and type_a == 'receptor'
            _type_a = type_b if rev else type_a
            _type_b = type_a if rev else type_b

            sources = (
                'CellPhoneDB'
                    if rec['annotation_strategy'] == 'curated' else
                '%s;CellPhoneDB' % (
                    ';'.join(
                        recomma.split(
                            rec['annotation_strategy'].replace(
                                'guidetopharmacology.org',
                                'Guide2Pharma'
                            )
                        )
                    )
                )
            )
            refs   = ';'.join(repmid.findall(rec['source']))

            yield (
                CellphonedbInteraction(
                    id_a = partner_b if rev else partner_a,
                    id_b = partner_a if rev else partner_b,
                    sources = sources,
                    references = refs,
                    interaction_type = '%s-%s' % (_type_a, _type_b),
                    type_a = _type_a,
                    type_b = _type_b,
                    is_ppi = _is_ppi,
                )
            )


def cellphonedb_complexes():

    annot = cellphonedb_complex_annotations()

    complexes = {}

    for cplex in annot.keys():
        key = cplex.__str__()

        if key in annot:
            cplex.add_attr('CellPhoneDB', annot[key])

        complexes[key] = cplex

    return complexes
