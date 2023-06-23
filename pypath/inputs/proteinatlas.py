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

import csv
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.inputs.common as inputs_common
import pypath.inputs.science as science


def get_proteinatlas(normal = True, pathology = True, cancer = True):

    result = {
        'normal':    collections.defaultdict(dict),
        'pathology': collections.defaultdict(dict),
    }


    def line(l):

        return l.strip('\n\r').split('\t')


    if normal:

        c = curl.Curl(urls.urls['proteinatlas']['normal'],
                    silent = False, large = True)
        fp = list(c.result.values())[0]
        hdr = line(fp.readline().strip())

        for l in fp:
            l = line(l)

            uniprots = mapping.map_name(l[1], 'genesymbol', 'uniprot')
            tissue = '%s:%s' % (l[2], l[3])

            for u in uniprots:
                result['normal'][tissue][u] = (l[4], l[5].strip())


    if cancer or pathology:

        c = curl.Curl(urls.urls['proteinatlas']['pathology'],
                    silent = False, large = True)
        fp = list(c.result.values())[0]
        hdr = line(fp.readline())

        for l in fp:

            l = line(l)
            uniprots = mapping.map_name(l[1], 'genesymbol', 'uniprot')
            tissue   = l[2]

            values = dict(
                (h, float(l[i + 3]) if '.' in l[i + 3] else int(l[i + 3]))
                for i, h in enumerate(hdr[3:])
                if len(l) and len(l[i + 3].strip())
            )

            for u in uniprots:

                result['pathology'][tissue][u] = values

    return dict((k, dict(v)) for k, v in iteritems(result))


def proteinatlas_annotations(normal = True, pathology = True, cancer = True):

    LEVELS = ('Not detected', 'Low', 'Medium', 'High')


    ProteinatlasAnnotation = collections.namedtuple(
        'ProtainatlasAnnotation',
        [
            'organ',
            'tissue',
            'level',
            'status',
            'n_not_detected',
            'n_low',
            'n_medium',
            'n_high',
            'prognostic',
            'favourable',
            'score',
            'pathology',
        ],
    )
    ProteinatlasAnnotation.__new__.__defaults__ = (
        (None,) * 4 + (False, False, None, False)
    )


    def n_or_none(ex, key):

        return ex[key] if key in ex else None


    data = get_proteinatlas(
        normal = normal,
        pathology = pathology,
        cancer = cancer,
    )

    result = collections.defaultdict(set)

    if normal:

        for tissue, gex in iteritems(data['normal']):

            organ = tissue

            if ':' in tissue:
                organ, tissue = tissue.split(':')

            organ = organ.strip()
            tissue = tissue.strip()

            for uniprot, ex in iteritems(gex):
                uniprots = mapping.map_name(uniprot, 'uniprot', 'uniprot')

                for _uniprot in uniprots:
                    result[_uniprot].add(
                        ProteinatlasAnnotation(
                            organ = organ,
                            tissue = tissue,
                            level = ex[0],
                            status = ex[1],
                        )
                    )

    if pathology or cancer:
        for condition, gex in iteritems(data['pathology']):
            for uniprot, ex in iteritems(gex):
                try:
                    effect, score = next(
                        i for i in iteritems(ex) if i[0] not in LEVELS
                    )
                    prognostic = not effect.startswith('unprognostic')
                    favourable = not effect.endswith('unfavourable')

                except StopIteration:
                    prognostic, favourable, score = None, None, None

                uniprots = mapping.map_name(uniprot, 'uniprot', 'uniprot')

                for _uniprot in uniprots:
                    result[_uniprot].add(
                        ProteinatlasAnnotation(
                            organ = condition,
                            tissue = condition,
                            level = max(
                                (i for i in iteritems(ex) if i[0] in LEVELS),
                                key = lambda i: i[1],
                                default = (None,),
                            )[0],
                            status = None,
                            n_not_detected = n_or_none(ex, 'Not detected'),
                            n_low = n_or_none(ex, 'Low'),
                            n_medium = n_or_none(ex, 'Medium'),
                            n_high = n_or_none(ex, 'High'),
                            prognostic = prognostic,
                            favourable = favourable,
                            score = score,
                            pathology = True,
                        )
                    )

    return dict(result)


def proteinatlas_subcellular_annotations():

    ProteinatlasSubcellularAnnotation = collections.namedtuple(
        'ProteinatlasSubcellularAnnotation',
        [
            'location',
            'status',
        ],
    )

    url = urls.urls['proteinatlas']['subcell']

    c = curl.Curl(
        url,
        large = True,
        silent = False,
        default_mode = 'r',
    )
    reader = csv.DictReader(
        c.files_multipart['subcellular_location.tsv'],
        delimiter = '\t',
    )

    result = collections.defaultdict(set)

    for rec in reader:
        uniprots = mapping.map_name(rec['Gene name'], 'genesymbol', 'uniprot')

        for uniprot in uniprots:
            for status in ('Enhanced', 'Supported', 'Uncertain'):
                if not rec[status]:
                    continue

                for location in rec[status].split(';'):
                    result[uniprot].add(ProteinatlasSubcellularAnnotation(
                        location = location,
                        status = status,
                    ))

    return dict(result)


def proteinatlas_secretome_annotations():

    ProteinatlasSecretomeAnnotation = collections.namedtuple(
        'ProteinatlasSecretomeAnnotation',
        [
            'mainclass',
            'secreted',
        ],
    )


    url = urls.urls['proteinatlas']['secretome']
    path = science.science_download(url)
    reader = inputs_common.read_xls(path)[1:]
    result = collections.defaultdict(set)

    for rec in reader:

        for uniprot_original in rec[2].split(','):

            uniprots = mapping.map_name(
                uniprot_original,
                'uniprot',
                'uniprot',
            )

            for uniprot in uniprots:
                result[uniprot].add(ProteinatlasSecretomeAnnotation(
                    mainclass = rec[3],
                    secreted = 'secreted' in rec[3].lower(),
                ))

    return dict(result)
