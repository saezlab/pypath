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
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy
import pypath.share.session as session

_logger = session.Logger(name = 'cellcall_input')
_log = _logger._log


def cellcall_download(extended = False, mouse = False):
    """
    Downloads a ligand-receptor-TF pathway dataset from CellCall
    (https://github.com/ShellyCoder/cellcall). This function downloads a
    single dataset, to download multiple datasets, see
    ``cellcall_download_all``. The pathway identifiers refer to KEGG
    pathways.

    Args
        extended (bool): CellCall has core and extended datasets, if this
            argument is True, the extended dataset will be retrieved.
        mouse (bool): CellCall has human and homology inferred mouse
            datasets. If this argument is True, the homology inferred
            dataset will be retrieved with mouse identifiers.

    Returns
        A list of dicts, each is a record as it provided by the CellCall
        database.
    """

    dataset = '%s%s' % (
        '_homology' if mouse else '',
        '_extended' if extended else '',
    )

    url = urls.urls['cellcall']['url'] % dataset
    c = curl.Curl(url, large = True, silent = False)

    return list(csv.DictReader(c.result, delimiter = '\t'))


def cellcall_download_all(extended = True, human = True, mouse = True):
    """
    Downloads ligand-receptor-TF pathway data from CellCall
    (https://github.com/ShellyCoder/cellcall). CellCall has core (high
    confidence) and extended datasets, human and homology inferred mouse
    datasets, 4 datasets in total. By default all these are downloaded
    here, with the parameters you can exclude the extended part and select
    the organism. The pathway identifiers refer to KEGG pathways.

    Args
        extended (bool): Use also the extended datasets.
        human (bool): Include human interactions.
        mouse (bool): Include mouse interactions.

    Returns
        A list of dicts, each is a record as it provided by the CellCall
        database.
    """

    result = []

    for ext, homo in itertools.product(*((True, False),) * 2):

        if (
            (extended or not ext) and
            (
                (human and not homo) or
                (mouse and homo)
            )
        ):

            dataset = cellcall_download(extended = ext, mouse = homo)
            _ = [
                (
                    rec.update(extended = ext),
                    rec.update(organism = 10090 if homo else 9606),
                )
                for rec in dataset
            ]
            result.extend(dataset)

    return result


def cellcall_interactions(extended = False, organism = 9606):
    """
    Ligand-receptor interactions from the CellCall database
    (https://github.com/ShellyCoder/cellcall).

    Args
        extended (bool): Include not only the core but also the extended
            set of interactions.
        organism (int,str): The organism to use, human (9606) and mouse
            (10090) are supported.

    Returns
        List of named tuples, each describing a ligand-receptor interaction.
    """

    record = collections.namedtuple(
        'CellcallInteraction',
        (
            'ligand_uniprot',
            'receptor_uniprot',
            'core',
        ),
    )


    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

    if ncbi_tax_id not in {9606, 10090}:

        msg = 'Unknown organism: `%s`.' % str(organism)
        _log(msg)
        raise ValueError(msg)

    human = ncbi_tax_id == 9606
    mouse = ncbi_tax_id == 10090
    raw = cellcall_download_all(
        extended = extended,
        human = human,
        mouse = mouse,
    )

    result = set()
    unmapped = set()

    for r in raw:

        ligands = mapping.map_name(
            r['Ligand_ID'],
            'entrez',
            'uniprot',
            ncbi_tax_id = ncbi_tax_id,
        )

        receptors = mapping.map_name(
            r['Receptor_ID'],
            'entrez',
            'uniprot',
            ncbi_tax_id = ncbi_tax_id,
        )

        if not ligands:

            unmapped.add(r['Ligand_ID'])

        if not receptors:

            unmapped.add(r['Receptor_ID'])

        for lig_up, rec_up in itertools.product(ligands, receptors):

            result.add(
                record(
                    ligand_uniprot = lig_up,
                    receptor_uniprot = rec_up,
                    core = not r['extended'],
                )
            )

    _log(
        'Could not find UniProt IDs for %u '
        'CellCall proteins.' % len(unmapped)
    )

    return list(result)


def cellcall_annotations(extended = False, organism = 9606):
    """
    Ligand and receptor annotations from the CellCall database
    (https://github.com/ShellyCoder/cellcall).

    Args
        extended (bool): Include not only the core but also the extended
            set of interactions.
        organism (int,str): The organism to use, human (9606) and mouse
            (10090) are supported.

    Returns
        Dict of annotations, keys are UniProt IDs, values are sets of
        annotations.
    """

    record = collections.namedtuple(
        'CellcallAnnotation',
        (
            'role',
        ),
    )


    interactions = cellcall_interactions(
        extended = extended,
        organism = organism,
    )
    result = collections.defaultdict(set)

    for i in interactions:

        result[i.ligand_uniprot].add(
            record(role = 'ligand')
        )
        result[i.receptor_uniprot].add(
            record(role = 'receptor')
        )

    return dict(result)
