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

import sys
import re
import collections
import itertools
import bs4
import csv

import pypath.inputs.common as inputs_common
import pypath.share.progress as progress
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common


def signor_interactions(
    organism = 9606,
    raw_records = False,
    expand_families = 0
):
    """
    Downloads the full dataset from SIGNOR (https://signor.uniroma2.it/).
    Returns the records with the most important fields.
    If ``raw_records`` is `True` it returns the table split to list of
    lists but unchanged content.

    Args
        organism (int, str): The NCBI Taxonomy ID or name of the organism.
            Human (9606), mouse (10090) and rat (10116) are available.
        raw_records (bool): Process the records or return them raw,
            as they are.
        expand_families (int): Expand protein families up to this size.
            Zero or one means no expansion.

    Return
        list: A list with processed records as named tuples or dicts of
            raw records if ``raw_records`` is True.
    """


    def process_name(name):

        isoform = None

        if name in families:

            main = (
                families[name]
                    if len(families[name]) <= expand_families else
                ()
            )

        elif name in complexes_by_id:

            main = complexes_by_id[name]

        else:

            main, isoform = inputs_common._try_isoform(name)
            main = (main,)

        return main, isoform


    SignorInteraction = collections.namedtuple(
        'SignorInteraction',
        (
            'source',
            'target',
            'source_isoform',
            'target_isoform',
            'source_type',
            'target_type',
            'effect',
            'mechanism',
            'ncbi_tax_id',
            'pubmeds',
            'direct',
            'ptm_type',
            'ptm_residue',
            'ptm_motif',
        )
    )

    families = signor_protein_families(organism = organism)
    complexes = signor_complexes(organism = organism)

    complexes_by_id = collections.defaultdict(set)

    for cplex in complexes.values():

        for cplex_id in cplex.ids['SIGNOR']:

            complexes_by_id[cplex_id].add(cplex)

    if isinstance(organism, int):

        if organism in taxonomy.taxids:

            _organism = taxonomy.taxids[organism]

        else:

            sys.stdout.write('\t:: Unknown organism: `%u`.\n' % organism)
            return []

    else:

        _organism = organism

    if _organism not in {'human', 'rat', 'mouse'}:

        return []

    url = urls.urls['signor']['all_url_new']
    binary_data = [
        (b'organism', _organism.encode('utf-8')),
        (b'format', b'csv'),
        (b'submit', b'Download'),
    ]

    c = curl.Curl(
        url,
        silent = False,
        large = True,
        follow = True,
        timeout = 180,
        binary_data = binary_data,
        return_headers = True,
    )

    reader = csv.DictReader(c.result, delimiter = '\t')

    if raw_records:

        return list(reader)

    result = []

    for line in reader:

        sources, source_isoform = process_name(line['IDA'])
        targets, target_isoform = process_name(line['IDB'])

        for source, target in itertools.product(sources, targets):

            this_record = SignorInteraction(
                source = source,
                target = target,
                source_isoform = source_isoform,
                target_isoform = target_isoform,
                source_type = line['TYPEA'],
                target_type = line['TYPEB'],
                effect = line['EFFECT'],
                mechanism = line['MECHANISM'],
                ncbi_tax_id = line['TAX_ID'],
                pubmeds = line['PMID'],
                direct = line['DIRECT'] == 'YES',
                ptm_type = line['MECHANISM'],
                ptm_residue = line['RESIDUE'],
                ptm_motif = line['SEQUENCE'],
            )

            result.append(this_record)

    return result


def signor_enzyme_substrate(organism = 9606):
    """
    Loads and processes Signor PTMs.
    Returns dict of dicts.
    """
    reres = re.compile(r'([A-Za-z]{3})([0-9]+)')
    result = []
    aalet = dict((k.lower().capitalize(), v)
                 for k, v in iteritems(common.aaletters))

    data = signor_interactions(organism = organism)

    for d in data:

        resm = reres.match(d.ptm_residue)

        if resm is not None:
            aa = aalet[resm.groups()[0].capitalize()]
            aanum = int(resm.groups()[1])
            typ = d.ptm_type,
            inst = d.ptm_motif.upper()
            result.append({
                'typ': d.ptm_type,
                'resnum': aanum,
                'instance': inst,
                'substrate': d.target,
                'start': aanum - 7,
                'end': aanum + 7,
                'kinase': d.source,
                'resaa': aa,
                'motif': inst,
                'enzyme_isoform': d.source_isoform,
                'substrate_isoform': d.target_isoform,
                'references': {d.pubmeds} if d.pubmeds != 'Other' else set()
            })

    return result


def signor_pathways(**kwargs):
    """
    Obtains pathway annotations from Signor.
    """

    url = urls.urls['signor']['list_url']
    baseurl = urls.urls['signor']['all_url_new']

    proteins_pathways = {}
    interactions_pathways = {}

    c = curl.Curl(url, silent = True)

    soup = bs4.BeautifulSoup(c.result, 'html.parser')

    pathway_names = [
        (opt['value'], opt.text)
        for opt in soup.find(
            'select', {'name': 'pathway_list'}
        ).findAll('option')
    ]

    prg = progress.Progress(
        len(pathway_names),
        'Downloading data from Signor',
        1,
        percent = False
    )

    for short, full in pathway_names:

        prg.step()

        if not short:

            continue

        binary_data = [
            (b'pathway_list', short.encode('ascii')),
            (b'submit', b'Download')
        ]

        c_pw = curl.Curl(
            baseurl,
            silent = True,
            binary_data = binary_data,
            encoding = 'utf-8',
        )

        #csv.DictReader(c_pw.result)

        sep = '@#@#@'
        lines = inputs_common.csv_sep_change(
            c_pw.result,
            '\t',
            sep
        ).split('\n')[1:]

        data = list(
            filter(
                lambda l:
                    len(l) > 6,
                map(
                    lambda l:
                        l.strip().split(sep),
                    lines
                )
            )
        )

        proteins_pathways[full] = set()
        interactions_pathways[full] = set()

        for row in data:

            for uniprot1, uniprot2 in itertools.product(
                mapping.map_name(row[4], 'uniprot', 'uniprot'),
                mapping.map_name(row[8], 'uniprot', 'uniprot'),
            ):

                proteins_pathways[full].add(uniprot1)
                proteins_pathways[full].add(uniprot2)

                interactions_pathways[full].add((uniprot1, uniprot2))

    prg.terminate()

    return proteins_pathways, interactions_pathways


def signor_pathway_annotations():

    SignorPathway = collections.namedtuple(
        'SignorPathway', ['pathway']
    )


    result = collections.defaultdict(set)

    proteins, interactions = signor_pathways()

    for pathway, uniprots in iteritems(proteins):

        record = SignorPathway(pathway = pathway)

        for uniprot in uniprots:

            result[uniprot].add(record)

    return dict(result)


def signor_protein_families(organism = 9606):
    #TODO: implement organism

    families = {}

    url = urls.urls['signor']['complexes']
    c = curl.Curl(
        url,
        binary_data = [(b'submit', b'Download protein family data')],
        large = True,
    )
    _ = next(c.result)

    for rec in c.result:

        rec = rec.split(';')
        components = [u.strip('\n\r" ') for u in rec[2].split(',')]
        families[rec[0]] = components

    return families


def signor_complexes(organism = 9606):
    #TODO: implement organism


    def process_on_hold(on_hold, complexes_by_id, complexes):

        on_hold_next = []

        for name, components, id_ in on_hold:

            components = [
                [comp.components for comp in complexes_by_id[comp_id]]
                    if comp_id in complexes_by_id else
                ((comp_id,),)
                for comp_id in components
            ]

            for components0 in itertools.product(*components):

                this_components = list(itertools.chain(*components0))

                if any(
                    comp.startswith('SIGNOR-C') for comp in this_components
                ):

                    on_hold_next.append((name, this_components, id_))

                else:

                    cplex = intera.Complex(
                        name = name.replace('"', '').strip(),
                        components = this_components,
                        sources = 'SIGNOR',
                        ids = id_,
                    )

                    complexes[cplex.__str__()] = cplex
                    complexes_by_id[id_].add(cplex)

        return on_hold_next, complexes_by_id, complexes


    complexes = {}
    on_hold = []

    families = signor_protein_families(organism = organism)

    url = urls.urls['signor']['complexes']
    c = curl.Curl(
        url,
        binary_data = [(b'submit', b'Download complex data')],
        large = True,
    )
    _ = next(c.result)

    complexes_by_id = collections.defaultdict(set)

    for rec in c.result:

        rec = rec.split(';')
        components = [u.strip('\n\r" ') for u in rec[2].split(',')]

        components = [
            families[comp] if comp in families else [comp]
            for comp in components
        ]

        for this_components in itertools.product(*components):

            # some complex contains other complexes
            if any(comp.startswith('SIGNOR-C') for comp in this_components):

                on_hold.append((rec[1], this_components, rec[0]))

            else:

                cplex = intera.Complex(
                    name = rec[1].replace('"', '').strip(),
                    components = this_components,
                    sources = 'SIGNOR',
                    ids = rec[0],
                )

                complexes[cplex.__str__()] = cplex
                complexes_by_id[rec[0]].add(cplex)

    while True:

        # complexes are defined recursively
        count_on_hold = len(on_hold)
        on_hold, complexes_by_id, complexes = (
            process_on_hold(on_hold, complexes_by_id, complexes)
        )

        if len(on_hold) == count_on_hold:

            break

    return complexes
