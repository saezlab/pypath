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

from pypath.omnipath.server import run

run.Rest(33333, serverclass = run.TableServer)

# without tables:
# run.Rest(33333, serverclass = run.TableServer, only_tables = [])

# run.Rest(33333, serverclass = run.TableServer, only_tables = {'annotations', 'intercell', 'interactions'})

#srv = server.Rest(
    #33333,
    #serverclass = server.TableServer,
    #only_tables = 'annotations',
#)

# Examples:
#
# curl 'http://omnipathdb.org/interactions?datasets=tfregulons&sources=FOXA1,FOXA2,FOXA3,FOXB1,FOXB2,FOXC1,FOXH1&genesymbols=1&fields=sources,tfregulons_level' > omnipath_webservice_test_get.tsv

# curl 'http://localhost:33333/interactions/FOXA1,FOXA2,FOXA3,FOXB1,FOXB2,FOXC1,FOXH1/SSBP2,ARAP2/AND/?datasets=tfregulons&genesymbols=1&fields=sources,tfregulons_level' > omnipath_webservice_test_get_path.tsv

# curl 'http://localhost:33333/interactions' -d 'datasets=tfregulons&sources=FOXA1,FOXA2,FOXA3,FOXB1,FOXB2,FOXC1,FOXH1&genesymbols=1&fields=sources,tfregulons_level' -H 'Content-Type: application/x-www-form-urlencoded' -X 'POST' > omnipath_webservice_test_post_form.tsv

# curl 'http://localhost:33333/interactions' -d '{"datasets": "tfregulons", "sources": "FOXA1,FOXA2,FOXA3,FOXB1,FOXB2,FOXC1,FOXH1", "genesymbols": "1", "fields": "sources,tfregulons_level"}' -H 'Content-Type: application/json' -X 'POST'  > omnipath_webservice_test_post_json.tsv
