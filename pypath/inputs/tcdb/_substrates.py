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

from collections.abc import Generator
import collections
import itertools
import csv
from pypath.share import curl
from pypath.resources import urls
from pypath.inputs import uniprot
import pandas as pd
__all__ = ['tcdb_substrate']
TcdbSubstrate = collections.namedtuple(
    'TcdbSubstrate',
    ['transporter_uniprot', 'substrate_id', 'substrate_name', 'location'],
)


def tcdb_substrate() -> Generator[TcdbSubstrate, None, None]:

    url_substrates = urls.urls['tcdb']['url_substrates']
    c = curl.Curl(url_substrates, large = True)
    tc_to_substrate = c.result
    url_acc2tc = urls.urls['tcdb']['url_acc2tc']
    c = curl.Curl(url_acc2tc, large = True)
    tc_to_uniprot = c.result

    #ID conversion
    tcdb_uniprot = collections.defaultdict(list)
    for line in tc_to_uniprot:
        uniprot_symbol, tcdb = line.strip('\n').split('\t')
        if uniprot_symbol:
            tcdb_uniprot[tcdb].append(uniprot_symbol.upper())

    all_uniprots = set().union(*tcdb_uniprot.values())
    all_uniprots = [
        u.split('.', maxsplit = 1)[0]
        for u in all_uniprots
        if uniprot.valid_uniprot(u)
    ]

    uniprot_locations = {}
    chunk_size = 98

    for i in range(0, len(all_uniprots), chunk_size):

        uniprot_locations.update(
            uniprot.uniprot_locations(
                accession=all_uniprots[i : i + chunk_size],
                organism=None,
                reviewed=None,
            )
        )

    for line in tc_to_substrate:

        tcdb, substrates = line.strip('\n').split('\t')

        for substrate in substrates.split('|'):

            substrate_id, substrate_name = substrate.split(';')

            for transporter_uniprot in tcdb_uniprot[tcdb]:

                yield TcdbSubstrate(
                    transporter_uniprot,
                    substrate_id,
                    substrate_name,
                    uniprot_locations.get(transporter_uniprot),
                )


if __name__ == "__main__":
    count = 0
    with open('tcdb_substrates.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['transporter_uniprot', 'substrate_id', 'substrate_name', 'location'])

        for record in tcdb_substrate():
            writer.writerow(record)
            count += 1

            if count % 1000 == 0:
                print(f"finished {count} lines...")
