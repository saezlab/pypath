#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
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

import urllib.request
import time
import itertools
import collections
import json
import importlib as imp


_urls = (
    'https://omnipathdb.org/',
    'http://localhost:33333/',
)

_queries = {
    'interactions': [
        {
            'types': (
                'post_transcriptional',
                'transcriptional',
                'post_translational',
                'mirna_transcriptional',
                'lncrna_post_transcriptional',
            ),
        },
        {
            'datasets': (
                'omnipath',
                'kinaseextra',
                'pathwayextra',
                'ligrecextra',
                'tf_target',
                'dorothea',
            ),
        },
        {
            'types': 'transcriptional',
            'datasets': 'dorothea',
            'dorothea_levels': 'A,B,C,D',
        },
    ]
}


class ServerTest(object):

    def __init__(self, outfile = None, urls = None, queries = None):

        self.outfile = outfile or (
            'omnipath-server-test-%s.tsv' % time.strftime('%Y-%m-%d %H:%M:%S')
        )
        self.urls = urls or _urls
        self.queries = queries or _queries


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        import importlib as imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def main(self):

        self.add_databases()
        self.generate_targets()
        self.retrieve()
        self.export()


    def add_databases(self):

        databases = collections.defaultdict(set)

        for url in self.urls:

            databases_url = '%sdatabases?format=json' % url

            con = urllib.request.urlopen(databases_url)

            for typ, dbs in json.loads(con.read()).items():

                databases[typ].update(set(dbs))

        for typ, dbs in databases.items():

            self.queries['interactions'].append({
                'resources': tuple(dbs),
                'types': typ,
            })


    def generate_targets(self):

        self.targets = []

        for query_type, params in self.queries.items():

            for param in params:

                keys, values = zip(*param.items())

                values = [
                    (val,) if isinstance(val, str) else val
                    for val in values
                ]

                for this_values in itertools.product(*values):

                    this_url_param = '%s?%s' % (
                        query_type,
                        '&'.join(
                            '%s=%s' % (key, value)
                            for key, value in zip(keys, this_values)
                        ),
                    )

                    self.targets.append(this_url_param)


    def retrieve(self):

        self.result = []

        for target in self.targets:

            self.result.append(
                [target] +
                [
                    self.retrieve_one('%s%s' % (url, target))
                    for url in self.urls
                ]
            )


    def retrieve_one(self, url):

        print('Retrieving %s' % url)

        con = urllib.request.urlopen(url)

        if con.getcode() == 200:

            content = con.read()

            if not content.startswith(b'Something is not'):

                return len(content.split(b'\n')) - 2


    def export(self, outfile = None):

        outfile = outfile or self.outfile

        with open(outfile, 'w') as fp:

            fp.write('\t%s\n' % '\t'.join(self.urls))

            fp.write('\n'.join(
                '\t'.join(str(i) for i in line)
                for line in self.result
            ))