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
import collections
import bs4

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping
import pypath.utils.reflists as reflists


def tcdb_families():

    retag = re.compile(r'<.*>')
    rethe = re.compile(r'^[tT]he (.*) (?:Super)?[Ff]amily')

    url = urls.urls['tcdb']['url_families']

    c = curl.Curl(url, large = False, silent = False)

    lines = bs4.BeautifulSoup(c.result, features = 'lxml').find('p').text

    return dict(
        (
            tcid,
            rethe.sub(r'\g<1>', family.replace('\t', ' '))
        )
        for tcid, family in
        (
            retag.sub('', line.strip()).split('\t', maxsplit = 1)
            for line in lines.strip().split('\n')
        )
    )


def tcdb_classes():

    refam = re.compile(r'(\d\.[A-Z]\.\d+)')
    retab = re.compile(r'\t+')

    url = urls.urls['tcdb']['url_acc2tc']

    c = curl.Curl(url, large = True, silent = False)

    result = {}

    for line in c.result:

        if not line:

            continue

        ac, tc = retab.split(line.rstrip())
        family = refam.search(tc).groups()[0]

        result[ac] = (tc, family)

    return result


def tcdb_annotations(organism = 9606):


    TcdbAnnotation = collections.namedtuple(
        'TcdbAnnotation',
        [
            'family',
            'tcid',
        ]
    )


    families = tcdb_families()
    classes = tcdb_classes()
    result = collections.defaultdict(set)

    for ac, (tc, family) in iteritems(classes):

        uniprots = mapping.map_name(
            ac,
            'uniprot',
            'uniprot',
            ncbi_tax_id = organism,
        )

        for uniprot in uniprots:

            if reflists.check(uniprot, 'uniprot', ncbi_tax_id = organism):

                result[uniprot].add(
                    TcdbAnnotation(
                        family = families[family],
                        tcid = tc,
                    )
                )

    return dict(result)


from collections.abc import Generator
import re
import pathlib as pl
import collections
import pandas as pd
import pyparsing
from pypath_common import _misc as _common
from pypath.share import curl
from pypath.utils import mapping, taxonomy

url = "https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py"
c = curl.Curl(url)
TC2Sub = c.result
url = "https://www.tcdb.org/cgi-bin/projectv/public/acc2tcid.py"
c = curl.Curl(url)
TC2Uni = c.result

#ID conversion
dic_TCDB_Uni = {}
for line in TC2Uni.strip().split("\n"):
    uniprot,tcdb = line.split("\t")
    if tcdb not in dic_TCDB_Uni:
        dic_TCDB_Uni[tcdb] = []
    dic_TCDB_Uni[tcdb].append(uniprot)

Uni2Che = []
for line in TC2Sub.strip().split("\n"):
    tcdb, substrate = line.split("\t")
    if tcdb in dic_TCDB_Uni:
        uniprot = ";".join(dic_TCDB_Uni[tcdb])
        Uni2Che.append((uniprot, substrate))


#protein location
import requests
url = "https://rest.uniprot.org/uniprotkb/search?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cft_intramem%2Ccc_subcellular_location%2Cft_transmem%2Cft_topo_dom&format=tsv&query=%28%28taxonomy_id%3A10090%29%29&size=500"
from requests.adapters import HTTPAdapter, Retry
import re
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

interactions = []
for batch, total in get_batch(url):
    for line in batch.text.splitlines()[1:]:
       interactions.append(line)
    print(f'{len(interactions)} / {total}')

fileObject = open('data/uniprot_subLocation.txt', 'w')
for i in interactions:
        fileObject.write(i)
        fileObject.write("\n")
fileObject.close()

#SLC table
import requests
url = "https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209652&file=msb209652-sup-0003-TableEV2.xlsx"
c = curl.Curl(url)
SLC_table = c.result

url = "https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209652&file=msb209652-sup-0002-TableEV1.xlsx"
