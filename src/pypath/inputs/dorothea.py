#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Helps to translate from the mouse data to human data
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from future.utils import iteritems

import re
import csv
import collections
import itertools
import functools
import pyreadr

import pypath.share.curl as curl
import pypath.resources.urls as urls


def get_dorothea_old(
        levels = {'A', 'B'},
        only_curated = False
    ):
    """
    Retrieves TF-target interactions from TF regulons.

    :param set levels:
        Confidence levels to be used.
    :param bool only_curated:
            Retrieve only literature curated interactions.

    Details
    -------
    TF regulons is a comprehensive resource of TF-target interactions
    combining multiple lines of evidences: literature curated databases,
    ChIP-Seq data, PWM based prediction using HOCOMOCO and JASPAR matrices
    and prediction from GTEx expression data by ARACNe.

    For details see https://github.com/saezlab/DoRothEA.
    """

    url = urls.urls['dorothea']['url'] % (
        'all' if 'E' in levels else
        'ABCD' if 'D' in levels else
        'ABC' if 'C' in levels else
        'AB' if 'B' in levels else
        'A'
    )

    c = curl.Curl(url, silent = False, large = True)
    _ = next(c.result)

    return (
        list(
            itertools.chain(
                ll[:4],
                (s == 'TRUE' for s in ll[4:8]),
                ll[-4:],
                [','.join(s for s in ll[-4:] if s)]
                if not only_curated else ll[8]
            )
        )
        for ll in (
            l.strip('\n\r').split('\t') for l in c.result
        ) if (
            ll[3] in levels and
            not only_curated or ll[4] == 'TRUE'
        )
    )


def get_dorothea(
        levels = {'A', 'B'},
        only_curated = False
    ):
    """
    Retrieves TF-target interactions from TF regulons.

    :param set levels:
        Confidence levels to be used.
    :param bool only_curated:
            Retrieve only literature curated interactions.

    Details
    -------
    TF regulons is a comprehensive resource of TF-target interactions
    combining multiple lines of evidences: literature curated databases,
    ChIP-Seq data, PWM based prediction using HOCOMOCO and JASPAR matrices
    and prediction from GTEx expression data by ARACNe.

    For details see https://github.com/saezlab/DoRothEA.
    """

    def process_sources(sources):
        
        upper = (
            'jaspar',
            'trred',
            'kegg',
            'trrust',
            'tred',
            'trrd',
            'hocomoco',
            'fantom4',
        )

        special_case = {
            'tfact': 'TFactS',
            'oreganno': 'ORegAnno',
            'reviews': 'DoRothEA-reviews',
            'HOCOMOCO_v11': 'HOCOMOCO-v11',
            'JASPAR_v2018': 'JASPAR-v2018',
        }

        revia = re.compile(r',|_via_')

        sources = functools.reduce(
            lambda s, r: s.replace(r, r.upper()),
            upper,
            sources,
        )

        sources = functools.reduce(
            lambda s, r: s.replace(*r),
            iteritems(special_case),
            sources,
        )

        return ','.join(
            '%s_DoRothEA' % s
            for s in revia.split(sources)
        )


    evidence_types = (
        'chipSeq',
        'TFbindingMotif',
        'coexpression',
        'curateddatabase'
    )

    DorotheaInteraction = collections.namedtuple(
        'DorotheaInteraction',
        [
            'tf',
            'target',
            'effect',
            'level',
            'curated',
            'chipseq',
            'predicted',
            'coexp',
            'curated_sources',
            'chipseq_sources',
            'predicted_sources',
            'coexp_sources',
            'all_sources',
            'pubmed',
            'kegg_pathways',
        ]
    )

    url = urls.urls['dorothea_git']['url']

    c = curl.Curl(
        url,
        silent = False,
        large = True,
        files_needed = ['database.csv'],
    )

    reader = csv.DictReader(c.result['database.csv'])

    for rec in reader:
        # process only the ones of the requested levels or if curated
        if (
            rec['score'] not in levels and
            not (
                only_curated and
                rec['is_evidence_curateddatabase'] == 'TRUE'
            )
        ):

            continue

        rec = dict(
            (k, v if v not in {'-', 'none'} else '')
            for k, v in iteritems(rec)
        )

        yield DorotheaInteraction(
            **dict(zip(
                DorotheaInteraction._fields,
                itertools.chain(
                    # TF, target, effect, score
                    (
                        rec[key] for key in
                        ('TF', 'target', 'effect', 'score')
                    ),
                    # boolean values for curated, chipseq, motif pred.
                    # and coexp
                    (
                        rec['is_evidence_%s' % key] == 'TRUE'
                        for key in evidence_types
                    ),
                    # databases & datasets
                    (
                        rec['which_%s' % key]
                        for key in evidence_types
                    ),
                    # all data sources (or only the curated ones)
                    (
                        process_sources(
                            ','.join(
                            rec[key]
                                for key in
                                ('which_%s' % evt for evt in evidence_types)
                                if rec[key]
                            )
                                if not only_curated else
                            rec['which_curateddatabase']
                        ),
                    ),
                    # PubMed and KEGG pw
                    (
                        rec['pubmedID_from_curated_resources'],
                        rec['kegg_pathway'],
                    )
                )
            ))
        )


def dorothea_rda_raw():

    url = urls.urls['dorothea_git']['rda']

    c = curl.Curl(url, silent = False, large = True)
    rdata_path = c.fileobj.name
    c.fileobj.close()

    rdata = pyreadr.read_r(rdata_path)['entire_database']

    return rdata


# synonyms
tfregulons_interactions = get_dorothea
get_tfregulons = get_dorothea
dorothea_interactions = get_dorothea
