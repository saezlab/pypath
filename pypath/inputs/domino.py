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

from past.builtins import xrange, range

import re
import itertools
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.progress as progress
import pypath.inputs.ontology as ontology
import pypath.internals.intera as intera
import pypath.internals.resource as resource_internals
import pypath.core.evidence as evidence

_logger = session.Logger(name = 'domino_input')
_log = _logger._log


def get_domino(none_values = False, outfile = None):
    """
    Returns
        A list of records with the following fields:
        header = ['uniprot_A', 'uniprot_B', 'isoform_A', 'isoform_B', #3
        'exp_method', 'references', 'taxon_A', 'taxon_B', #7
        'role_A', 'role_B', 'binding_site_range_A', 'binding_site_range_B', #11
        'domains_A', 'domains_B', 'ptm_residue_A', 'ptm_residue_B', #15
        'ptm_type_mi_A', 'ptm_type_mi_B', 'ptm_type_A', 'ptm_type_B', #19
        'ptm_res_name_A', 'ptm_res_name_B', 'mutations_A', 'mutations_B', #23
        'mutation_effects_A', 'mutation_effects_B', 'domains_interpro_A', #26
        'domains_interpro_B', 'negative'] #28
    """

    DominoRecord = collections.namedtuple(
        'DominoRecord',
        (
            'uniprot_A',
            'uniprot_B',
            'isoform_A',
            'isoform_B',
            'exp_method',
            'references',
            'taxon_A',
            'taxon_B',
            'role_A',
            'role_B',
            'binding_site_range_A',
            'binding_site_range_B',
            'domains_A',
            'domains_B',
            'ptm_residue_A',
            'ptm_residue_B',
            'ptm_type_mi_A',
            'ptm_type_mi_B',
            'ptm_type_A',
            'ptm_type_B',
            'ptm_res_name_A',
            'ptm_res_name_B',
            'mutations_A',
            'mutations_B',
            'mutation_effects_A',
            'mutation_effects_B',
            'domains_interpro_A',
            'domains_interpro_B',
            'negative',
        ),
    )

    result = []
    taxid = re.compile(r'taxid:(.*)\([a-zA-Z ]*\)')
    miont = re.compile(r'MI:[0-9]{4}\((.*)\)')
    binds = re.compile(r'([-0-9]*);.*')
    domai = re.compile(r'.*;.*;.*\((.*)\)')
    dipro = re.compile(r'.*;.*;.+:(IPR[0-9]*).*')
    ptmrs = re.compile(r'([-0-9]*);.*')
    ptmmi = re.compile(r'[0-9]*;(MI:[0-9]*)\(.*\);.*;.*')
    ptmrn = re.compile(
        r'.*sequence:[\s]*[0-9]+-[0-9]+[\s]*:[\s]*([A-Z]{10,}).*')
    ptmty = re.compile(r'[0-9]*;MI:[0-9]*\((.*)\);.*;.*')
    refrs = re.compile(r'(pubmed|doi):["]*([-0-9a-zA-Z\.\(\)/]*)["]*')
    url = urls.urls['domino']['rescued']
    c = curl.Curl(url, silent = False, large = True)
    data = c.result
    _ = next(data)

    for r in data:

        r = r.strip().split('\t')

        if len(r) < 39:

            continue

        this_row = [
            None if ':' not in r[0] else r[0].split(':')[1].split('-')[0],
            None if ':' not in r[1] else r[1].split(':')[1].split('-')[0],
            '1'  if '-' not in r[0] else r[0].split('-')[1],
            '1'  if '-' not in r[1] else r[1].split('-')[1],
            miont.match(r[6]).groups(1)[0] if miont.match(r[6]) else None,
            refrs.match(r[8]).groups(1)[1] if refrs.match(r[8]) else None,
            taxid.match(r[9]).groups(1)[0] if taxid.match(r[9]) else None,
            taxid.match(r[10]).groups(1)[0] if taxid.match(r[10]) else None,
            miont.match(r[11]).groups(1)[0] if miont.match(r[11]) else None,
            miont.match(r[17]).groups(1)[0] if miont.match(r[16]) else None,
            ';'.join(
                binds.match(x).groups(1)[0] if binds.match(x) else ''
                for x in r[32].split(',')
            ),
            ';'.join(
                binds.match(x).groups(1)[0] if binds.match(x) else ''
                for x in r[33].split(',')
            ),
            ';'.join(
                domai.match(x).groups(1)[0] if domai.match(x) else ''
                for x in r[32].split(',')
            ),
            ';'.join(
                domai.match(x).groups(1)[0] if domai.match(x) else ''
                for x in r[33].split(',')
            ),
            ';'.join(
                ptmrs.match(x).groups(1)[0] if ptmrs.match(x) else ''
                for x in r[34].split('|')
            ),
            ';'.join(
                ptmrs.match(x).groups(1)[0] if ptmrs.match(x) else ''
                for x in r[35].split('|')
            ),
            ';'.join(
                ptmmi.match(x).groups(1)[0] if ptmmi.match(x) else ''
                for x in r[34].split('|')
            ),
            ';'.join(
                ptmmi.match(x).groups(1)[0] if ptmmi.match(x) else ''
                for x in r[35].split('|')
            ),
            ';'.join(
                ptmty.match(x).groups(1)[0] if ptmty.match(x) else ''
                for x in r[34].split('|')
            ),
            ';'.join(
                ptmty.match(x).groups(1)[0] if ptmty.match(x) else ''
                for x in r[35].split('|')
            ),
            ';'.join(
                ptmrn.match(x).groups(1)[0] if ptmrn.match(x) else ''
                for x in r[34].split('|')
            ),
            ';'.join(
                ptmrn.match(x).groups(1)[0] if ptmrn.match(x) else ''
                for x in r[35].split('|')
            ),
            ';'.join(
                ptmrs.match(x).groups(1)[0] if ptmrs.match(x) else ''
                for x in r[36].split('|')
            ), ';'.join(
                ptmrs.match(x).groups(1)[0] if ptmrs.match(x) else ''
                for x in r[37].split('|')
            ),
            ';'.join(
                ptmty.match(x).groups(1)[0] if ptmty.match(x) else ''
                for x in r[36].split('|')
            ),
            ';'.join(
                ptmty.match(x).groups(1)[0] if ptmty.match(x) else ''
                for x in r[37].split('|')
            ),
            dipro.match(r[32]).groups(1)[0] if dipro.match(r[32]) else '',
            dipro.match(r[33]).groups(1)[0] if dipro.match(r[33]) else '',
            '0' if r[38].strip() == '-' else '1',
        ]

        if not none_values:

            this_row = ['' if x is None else x for x in this_row]
            this_row = DominoRecord(*this_row)

        result.append(this_row)

    if outfile:

        _log('Saving data into `%s`.' % outfile)

        with open(outfile, 'w') as outf:

            outf.write('\t'.join(header) + '\n')

            for r in result:

                outf.write(
                    '\t'.join('' if x is None else x for x in r) + '\n'
                )

    return result


def domino_interactions():

    domino = get_domino()

    interactions = [
        l for l in domino
        if (
            l[0] and
            l[1] and
            ''.join(l[5]) and
            ''.join([
                l[i]
                for i in range(10, 12) + range(14, 22) + range(24, 26)
            ]) and
            l[28] != '1'
        )
    ]

    return interactions


def domino_ddi():

    domi = domino_enzsub()

    return domi['ddi']


def domino_enzsub():
    """
    Returns
        A dict of two elements: `ddi` contains domain-domain, while `dmi`
        domain-motif interactions. The latter includes protein-PTM
        interactions.
    """

    domino_resource = resource_internals.EnzymeSubstrateResource(
        name = 'DOMINO',
        input_method = 'domino.domino_enzsub',
    )

    domino = get_domino()

    try:

        miont = ontology.ontology('MI')

    except:

        miont = {}

    dmi = []
    ddi = []
    prg = progress.Progress(len(domino), 'Processing DOMINO', 11)

    ptm_types = {
        "o4'-phospho-tyrosine": 'phosphorylation',
        'phosphorylated residue': 'phosphorylation',
        'o-phospho-threonine': 'phosphorylation',
        'o-phospho-serine': 'phosphorylation',
        'n6-methyl-lysine': 'methylation',
        'n6,n6,n6-trimethyl-lysine': 'trimethylation',
        'n6,n6-dimethyl-lysine': 'dimethylation',
        'acetylated residue': 'acetylation',
    }

    for l in domino:

        prg.step()

        if (
            (
                l[14].strip() != '' or
                l[15].strip() != '' or
                (
                    l[10] != '' and
                    l[11] != ''
                )
            ) and
            len(l[0]) > 0 and
            len(l[1]) > 0
        ):

            uniprot1 = l[0]
            uniprot2 = l[1]

            # ptms
            if (
                '-' not in l[14] and
                '-' not in l[15]
            ):

                ptmre12 = [int(x) for x in l[14].split(';')] if l[14] else []
                ptmre21 = [int(x) for x in l[15].split(';')] if l[15] else []
                ptmty12 = l[16].split(';') if l[16] else [None] * len(ptmre12)
                ptmty12 = [
                    ptm_types[miont[x]] if x in miont else None
                    for x in ptmty12
                ]
                ptmrn12 = l[20].split(';') if l[20] else [None] * len(ptmre12)

                ptmrn12 = [
                    None
                        if (
                            x is None or
                            x == '' or
                            len(x) < min(ptmre12[i] - 1, 11)
                        ) else
                    x[10]
                        if ptmre12[i] > 10 else
                    x[ptmre12[i] - 1]
                    for i, x in enumerate(ptmrn12)
                ]
                ptmty21 = l[17].split(';') if l[17] else [None] * len(ptmre12)
                ptmty21 = [
                    ptm_types[miont[x]] if x in miont else None
                    for x in ptmty21
                ]
                ptmrn21 = l[21].split(';') if l[21] else [None] * len(ptmre21)
                ptmrn21 = [
                    None
                        if (
                            x is None or
                            x == '' or
                            len(x) < min(ptmre21[i] - 1, 11)
                        ) else
                    x[10]
                        if ptmre21[i] > 10 else
                    x[ptmre21[i] - 1]
                    for i, x in enumerate(ptmrn21)
                ]

                for i, resnum in enumerate(ptmre12):

                    res = intera.Residue(resnum, ptmrn12[i], uniprot2)
                    ptm = intera.Ptm(
                        uniprot2,
                        typ = ptmty12[i] or 'unknown',
                        residue = res,
                        evidences = evidence.Evidence(
                            resource = domino_resource,
                        ),
                    )
                    dom = intera.Domain(uniprot1)
                    dm = intera.DomainMotif(
                        domain = dom,
                        ptm = ptm,
                        evidences = evidence.Evidence(
                            resource = domino_resource,
                            references = l[5].split(';'),
                        ),
                    )
                    dmi.append(dm)

            # binding sites
            if l[10] and l[11]:

                try:

                    bssrt1 = [
                        int(x.split('-')[0])
                        for x in l[10].split(';')
                        if x != '' and x != '0'
                    ]
                    bsend1 = [
                        int(x.split('-')[1])
                        for x in l[10].split(';')
                        if x != '' and x != '0'
                    ]
                    bssrt2 = [
                        int(x.split('-')[0])
                        for x in l[11].split(';')
                        if x != '' and x != '0'
                    ]
                    bsend2 = [
                        int(x.split('-')[1])
                        for x in l[11].split(';')
                        if x != '' and x != '0'
                    ]

                except:

                    sys.stdout.write('Error processing line:\n')
                    sys.stdout.write(l)
                    sys.stdout.write('\n')
                    sys.stdout.flush()

                    return None

                bs1 = []
                bs2 = []

                if l[26]:

                    for i, n in enumerate(bssrt1):

                        bs1.append(
                            intera.Domain(
                                protein = uniprot1,
                                domain = l[26],
                                start = bssrt1[i],
                                end = bsend1[i],
                                domain_id_type = 'interpro',
                                isoform = l[2],
                            )
                        )

                else:

                    for i, n in enumerate(bssrt1):

                        mot = intera.Motif(
                            protein = uniprot1,
                            start = bssrt1[i],
                            end = bsend1[i],
                            isoform = l[2],
                        )
                        bs1.append(
                            intera.Ptm(
                                protein = uniprot1,
                                motif = mot,
                                evidences = evidence.Evidence(
                                    resource = domino_resource,
                                ),
                                isoform = l[2],
                            )
                        )

                if l[27]:

                    for i, n in enumerate(bssrt2):

                        bs2.append(
                            intera.Domain(
                                protein = uniprot2,
                                domain = l[27],
                                start = bssrt2[i],
                                end = bsend2[i],
                                domain_id_type = 'interpro',
                                isoform = l[3],
                            )
                        )

                else:

                    for i, n in enumerate(bssrt2):

                        mot = intera.Motif(
                            protein = uniprot2,
                            start = bssrt2[i],
                            end = bsend2[i],
                            isoform = l[3],
                        )
                        bs2.append(
                            intera.Ptm(
                                protein = uniprot2,
                                motif = mot,
                                evidences = evidence.Evidence(
                                    resource = domino_resource,
                                ),
                            )
                        )

                for one, two in itertools.product(bs1, bs2):


                    if (
                        one.__class__.__name__ == 'Domain' and
                        two.__class__.__name__ == 'Domain'
                    ):

                        dd = intera.DomainDomain(
                            one,
                            two,
                            sources = 'DOMINO',
                        )
                        ddi.append(dd)

                    if (
                        one.__class__.__name__ == 'Domain' and
                        two.__class__.__name__ == 'Ptm'
                    ):

                        dm = intera.DomainMotif(
                            domain = one,
                            ptm = two,
                            evidences = evidence.Evidence(
                                resource = domino_resource,
                                references = l[6].split(';')
                            ),
                        )
                        dmi.append(dm)

                    if (
                        two.__class__.__name__ == 'Domain' and
                        one.__class__.__name__ == 'Ptm'
                    ):

                        dm = intera.DomainMotif(
                            domain = two,
                            ptm = one,
                            evidences = evidence.Evidence(
                                resource = domino_resource,
                                references = l[6].split(';')
                            ),
                        )
                        dmi.append(dm)

    prg.terminate()

    return {'ddi': ddi, 'dmi': dmi}
