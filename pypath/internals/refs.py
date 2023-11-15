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
import os
import sys
import webbrowser
import pandas as pd
try:
    import cPickle as pickle
except:
    import pickle

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls
from pypath.inputs import pubmed as pubmed_input
import pypath.share.cache as cache
import pypath.inputs.pubmed as pubmed


class Reference(object):

    __slots__ = ['pmid']

    def __init__(self, pmid):
        self.pmid = str(pmid).strip()

    def __eq__(self, other):
        return self.pmid == other.pmid

    def __hash__(self):
        return hash(self.pmid)

    def open(self):
        pubmed_input.open_pubmed(self.pmid)

    def __str__(self):
        return self.pmid

    def info(self):
        return pubmed_input.get_pubmeds([self.pmid])

    def __repr__(self):

        return '<Reference: %s>' % self.pmid


def get_pubmed_data(
        pp,
        cachefile = None,
        htp_threshold = 20
    ):
    """
    For one PyPath object, obtains metadata for all PubMed IDs
    through NCBI E-utils.

    :param pp:
        ``pypath.PyPath`` object
    :param htp_threshold:
        The number of interactions for one reference
        above the study considered to be high-throughput.
    """


    if cachefile is None:

        cachefile = cache.cache_item('pubmed_cache')

    if htp_threshold is not None:
        pp.htp_stats()

    pubmeds = common.unique_list(
        common.flat_list([[r.pmid for r in e['references']]
                         for e in pp.graph.es]))

    if htp_threshold is not None:
        pubmeds = set(pubmeds) - pp.htp[htp_threshold]['htrefs']

    notpmid = [i for i in pubmeds if not i.isdigit()]

    sys.stdout.write('\t:: Number of non PubMed ID references: %u\n' %
                     len(notpmid))

    pmdata = {}
    if os.path.exists(cachefile):
        sys.stdout.write('\t:: Loading data previously downloaded '
                         'from PubMed, from file `%s`\n' % cachefile)
        pmdata = pickle.load(open(cachefile, 'rb'))

    missing = list(set(pubmeds) - set(pmdata.keys()))
    sys.stdout.write('\t:: Downloading data from PubMed about %s papers\n' %
                     len(missing))
    cached_pubmeds_len = len(pmdata)
    pmdata_new = pubmed_input.get_pubmeds(missing)
    pmdata.update(pmdata_new)

    sys.stdout.write('\t:: Saving PubMed data to file `%s`\n' % cachefile)

    if len(pmdata) > cached_pubmeds_len:
        pickle.dump(pmdata, open(cachefile, 'wb'))

    pmdata = dict(i for i in pmdata.items() if i[0] in pubmeds)

    points = []
    earliest = []

    for e in pp.graph.es:

        for s, rs in iteritems(e['refs_by_source']):

            pms = [
                r.pmid for r in rs
                if (htp_threshold is None or r.pmid not in pp.htp[
                    htp_threshold]['htrefs']
                    ) and r.pmid in pmdata and 'pubdate' in pmdata[r.pmid]
            ]
            if len(pms) > 0:
                yrs = [int(pmdata[pm]['pubdate'][:4]) for pm in pms]
                earliest.append((s, 0, min(yrs), '', e.index))
                for pm in pms:
                    points.append((s, pm, int(pmdata[pm]['pubdate'][:4]),
                                   pmdata[pm]['source'], e.index))

    points = common.unique_list(points)
    earliest = common.unique_list(earliest)

    points = pd.DataFrame.from_records(points)
    earliest = pd.DataFrame.from_records(earliest)
    points.columns = ['database', 'pmid', 'year', 'journal', 'eid']
    earliest.columns = ['database', 'none', 'year', 'none', 'eid']

    return points, earliest
