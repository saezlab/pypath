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
import collections

import pandas as pd

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping


BiogpsDataset = collections.namedtuple(
    'BiogpsDataset',
    (
        'organism',
        'label',
        'url',
    ),
)


def biogps_datasets():
    """
    List the expression profile datasets available from BioGPS.

    Returns
        (list): Named tuples with the label, organism and URL of the datasets.
    """

    biogps_urls = urls.urls['biogps']
    url = biogps_urls['url']
    datasets = biogps_urls['datasets']


    return [
        BiogpsDataset(
            organism = label.split('_', maxsplit = 1)[0],
            label = label,
            url = url % fname,
        )
        for label, fname in iteritems(datasets)
    ]


def biogps_download(dataset):
    """
    Retrieve one BioGPS expression profile dataset.

    Args
        dataset (str,BiogpsDataset): Either the label of a dataset or a
            `BiogpsDataset` object. For a list of available datasets, see
            `biogps_datasets`.

    Returns
        (pandas.DataFrame): A data frame of expression values, columns are
            tissues, cell types or cell lines, rows are microarray probes
            representing transcripts.
    """

    biogps_urls = urls.urls['biogps']
    label = dataset.label if hasattr(dataset, 'label') else dataset
    url = biogps_urls['url'] % biogps_urls['datasets'][label]
    c = curl.Curl(url, silent = False, large = True)

    the_file = (
        common.first(c.result.values())
            if isinstance(c.result, dict) else
        c.result
    )

    fileobj = the_file if hasattr(the_file, 'name') else c.fileobj
    sep = ',' if fileobj.name.endswith('csv') else '\t'
    header = next(the_file).split(sep)
    header = [
        re.sub(
            '\W+',
            '_',
            re.sub('^\d+', '', h).strip('')
        ).strip('_')
        for h in header
    ]
    header[0] = 'probe'

    hdr_dup = {}

    for i in xrange(len(header)):

        field = header[i]
        header[i] = (
            field
                if field not in hdr_dup else
            '%s_%u' % (field, hdr_dup[field])
        )
        hdr_dup[field] = hdr_dup.get(field, 0)
        hdr_dup[field] += 1

    content = [
        [
            value if j == 0 else float(value)
            for j, value in enumerate(row.strip().split(sep))
        ]
        for row in the_file
    ]

    df = pd.DataFrame(content, columns = header)

    return df


def biogps_download_all(organism = None, exclude = None, only = None):
    """
    Downloads all expression data from BioGPS.

    Args
        organism (str,int): Restrict the download to this organism. Human,
            mouse and rat datasets are available.
        exclude (str,set): One or more datasets to exclude. Should be
            dataset labels as shown by `biogps_datasets`.
        only (str,set): Restrict the download only to these datasets.

    Returns
        (dict): Dict of pandas data frames, keys are dataset labels, values
            are data frames of expression values, columns are tissues, cell
            types or cell lines, rows are microarray probes representing
            transcripts.
    """

    if organism:

        organism = taxonomy.ensure_common_name(organism).lower()

    exclude = common.to_set(exclude)
    only = common.to_set(only)
    datasets = biogps_datasets()

    result = {}

    for dataset in datasets:

        if (
            (
                not organism or
                dataset.organism == organism
            ) and (
                not only or
                dataset.label in only
            ) and
            dataset not in exclude
        ):

            result[dataset.label] = biogps_download(dataset)

    return result


def biogps_annotations(organism = 9606, exclude = None, only = None):
    """
    Expression data from BioGPS compiled in the annotation format commonly
    used in this module.

    Args
        organism (int,str): Name or ID of the organism. This function
            is able to process data for one organism at once.
        exclude (str,set): One or more datasets to exclude. Should be
            dataset labels as shown by `biogps_datasets`.
        only (str,set): Restrict the download only to these datasets.

    Returns
        (dict): Dict of annotations, keys are UniProt IDs, values are
            sets of named tuples of annotations.
    """

    BiogpsAnnotation = collections.namedtuple(
        'BiogpsAnnotation',
        (
            'dataset',
            'sample',
            'probe',
            'expression',
        ),
    )

    data = biogps_download_all(
        organism = organism,
        exclude = exclude,
        only = only,
    )

    result = collections.defaultdict(set)

    for label, df in iteritems(data):

        for row in df.itertuples(index = False):

            uniprots = mapping.map_name(
                row[0],
                'affy',
                'uniprot',
                ncbi_tax_id = organism,
            )

            for uniprot in uniprots:

                result[uniprot].update(
                    {
                        BiogpsAnnotation(
                            dataset = label,
                            sample = sample,
                            probe = row[0],
                            expression = expr,
                        )
                        for sample, expr in zip(row._fields[1:], row[1:])
                    }
                )

    return result
