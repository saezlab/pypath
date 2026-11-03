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

import re
import csv
import collections
import itertools
import functools
import pyreadr

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.session as session
import pypath.utils.taxonomy as taxonomy


_logger = session.Logger(name = 'dorothea_input')


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


_resources_upper = (
    'jaspar',
    'trred',
    'kegg',
    'trrust',
    'tred',
    'trrd',
    'hocomoco',
    'fantom4',
    'pazar',
)

_resources_special_case = {
    'tfact': 'TFactS',
    'tf_act': 'TFactS',
    'htri_db': 'HTRIdb',
    'int_act': 'IntAct',
    'fantom_4': 'FANTOM4',
    'oreganno': 'ORegAnno',
    'reviews': 'DoRothEA-reviews',
    'HOCOMOCO_v11': 'HOCOMOCO-v11',
    'hocomoco_v11': 'HOCOMOCO-v11',
    'JASPAR_v2018': 'JASPAR-v2018',
    'remap': 'ReMap',
    'gtex': 'ARACNe-GTEx',
    'nfi_regulome_db': 'NFIRegulomeDB',
    'tf_e': 'TFe',
    'reg_network': 'RegNetwork',
}


def _process_resources(sources):

    if sources == 'none':

        return ''

    revia = re.compile(r',|_via_')

    sources = functools.reduce(
        lambda s, r: s.replace(r, r.upper()),
        _resources_upper,
        sources,
    )

    sources = functools.reduce(
        lambda s, r: s.replace(*r),
        iteritems(_resources_special_case),
        sources,
    )

    return ','.join(revia.split(sources))



def get_dorothea_old(
        levels = {'A', 'B'},
        only_curated = False
    ):
    """
    Retrieves TF-target interactions from DoRothEA.

    :param set levels:
        Confidence levels to be used.
    :param bool only_curated:
            Retrieve only literature curated interactions.

    Details
    -------
    DoRothEA is a comprehensive resource of TF-target interactions
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


def dorothea_old_csv(
        levels = {'A', 'B'},
        only_curated = False
    ):
    """
    Retrieves TF-target interactions from DoRothEA.

    :param set levels:
        Confidence levels to be used.
    :param bool only_curated:
            Retrieve only literature curated interactions.

    Details
    -------
    Note: this method processes DoRothEA from an old CSV file generated in
    2018. For an up to date version of DoRothEA please use the
    ``dorothea_interactions`` method.
    DoRothEA is a comprehensive resource of TF-target interactions
    combining multiple lines of evidences: literature curated databases,
    ChIP-Seq data, PWM based prediction using HOCOMOCO and JASPAR matrices
    and prediction from GTEx expression data by ARACNe.

    For details see https://github.com/saezlab/DoRothEA.
    """

    evidence_types = (
        'chipSeq',
        'TFbindingMotif',
        'coexpression',
        'curateddatabase'
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
                        _process_resources(
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


def dorothea_rda_raw(organism = 9606):
    """


    :param int,str organism:
        Name or NCBI Taxonomy ID of the organism. Human and mouse are
        supported. If `None`, the human interactions will be returned
        with additional details included.
    """

    _organism = taxonomy.ensure_ncbi_tax_id(organism)

    if _organism not in (9606, 10090, None):

        msg = (
            'DoRothEA: invalid organism: `%s`. Only human and mouse '
            'are supported.' % str(organism)
        )
        _logger._log(msg)
        raise ValueError(msg)

    fname = (
        'entire_database'
            if _organism is None else
        'dorothea_%s' % ('hs' if _organism == 9606 else 'mm')
    )

    url = urls.urls['dorothea_git']['rda'] % fname

    c = curl.Curl(url, silent = False, large = True)
    rdata_path = c.fileobj.name
    c.fileobj.close()

    rdata = None

    try:
        rdata = pyreadr.read_r(rdata_path)[fname]
    except pyreadr.custom_errors.LibrdataError as e:
        _logger._log(
            'Could not parse DoRothEA data from Rdata file: '
            '`%s`. '
            'Make sure your `pyreadr` installation supports the xz '
            'compression.' % e.args[0]
        )

    return rdata


def dorothea_full_raw(organism = 9606):
    """
    DoRothEA data as it is provided in the R package.

    Args
        organism (int,str): Name or NCBI Taxonomy ID of the organism. The
            complete DoRothEA database (with all the details about the
            original sources) is available only for human.

    Returns
        (pandas.DataFrame): A data frame of TF-target interactions from
            DoRothEA.
    """

    _organism = taxonomy.ensure_ncbi_tax_id(organism)

    if _organism != 9606 and organism:

        msg = (
            'DoRothEA: invalid organism: `%s`. The full database is '
            'available only for human. To have them for other organisms, '
            'you can load DoRothEA in a Network object and use the homology '
            'translation function.' % str(organism)
        )
        _logger._log(msg)
        raise ValueError(msg)

    dorothea_full = dorothea_rda_raw(organism = None)

    return dorothea_full


def dorothea_interactions(
        organism = 9606,
        levels = {'A', 'B', 'C', 'D'},
        only_curated = False,
        confidence_pairwise = True,
    ):
    """
    Retrieves TF-target interactions from TF regulons.

    :param int,str organism:
        Name or NCBI Taxonomy ID of the organism. Only human is available.
    :param set levels:
        Confidence levels to be used.
    :param bool only_curated:
        Retrieve only literature curated interactions.

    Details
    -------
    TF regulons is a comprehensive resource of TF-target interactions
    combining multiple lines of evidences: literature curated databases,
    ChIP-Seq data, PWM based prediction using HOCOMOCO and JASPAR matrices
    and prediction from GTEx expression data by ARACNe. As KEGG is not longer
    part of the public version of DoRothEA the `kegg_pathways` field is
    always empty.

    For details see https://github.com/saezlab/DoRothEA.
    """

    evidence_types = (
        'curated',
        'chip_seq',
        'tfbs', # inferred
        'inferred', # coexp
    )

    df = dorothea_full_raw(organism = organism)
    df = df[df.confidence.isin(levels)]

    if only_curated:

        df = df[df.is_evidence_curated]

    for rec in df.itertuples():

        yield(
            DorotheaInteraction(
                **dict(zip(
                    DorotheaInteraction._fields,
                    itertools.chain(
                        # TF, target, effect, score
                        (
                            rec.tf,
                            rec.target,
                            int(rec.mor),
                            rec.confidence,
                        ),

                        # boolean values for curated, chipseq, motif pred.
                        # and coexp
                        (
                            getattr(rec, 'is_evidence_%s' % evt)
                            for evt in evidence_types
                        ),
                        # databases & datasets
                        (
                            _process_resources(
                                getattr(rec, 'which_%s' % evt)
                            )
                            for evt in evidence_types
                        ),
                        # all data sources (or only the curated ones)
                        (
                            _process_resources(
                                ','.join(
                                    getattr(rec, key)
                                    for key in
                                    (
                                        'which_%s' % evt
                                        for evt in evidence_types
                                    )
                                    if getattr(rec, key) != 'none'
                                )
                                    if not only_curated else
                                rec.which_curated
                            ),
                        ),
                        # PubMed and KEGG pw
                        (
                            rec.pubmed_id if rec.pubmed_id.isdigit() else '',
                            '',
                        ),
                    )
                ))
            )
        )


# synonyms
dorothea_interactions_old = dorothea_old_csv
tfregulons_interactions_old = dorothea_old_csv
get_tfregulons = dorothea_rda_raw
tfregulons_interactions = dorothea_interactions
