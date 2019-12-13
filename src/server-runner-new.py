#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2014-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

from pypath import server

server.Rest(33333, serverclass = server.TableServer)

#server.Rest(33333, serverclass = server.TableServer, only_tables = 'intercell')

srv = server.Rest(
    33333,
    serverclass = server.TableServer,
    only_tables = 'annotations',
)

# Examples:
#
# curl 'http://omnipathdb.org/interactions?datasets=tfregulons&sources=FOXA1,FOXA2,FOXA3,FOXB1,FOXB2,FOXC1,FOXH1&genesymbols=1&fields=sources,tfregulons_level' > omnipath_webservice_test_get.tsv

# curl 'http://localhost:33333/interactions/FOXA1,FOXA2,FOXA3,FOXB1,FOXB2,FOXC1,FOXH1/SSBP2,ARAP2/AND/?datasets=tfregulons&genesymbols=1&fields=sources,tfregulons_level' > omnipath_webservice_test_get_path.tsv

# curl 'http://localhost:33333/interactions' -d 'datasets=tfregulons&sources=FOXA1,FOXA2,FOXA3,FOXB1,FOXB2,FOXC1,FOXH1&genesymbols=1&fields=sources,tfregulons_level' -H 'Content-Type: application/x-www-form-urlencoded' -X 'POST' > omnipath_webservice_test_post_form.tsv

# curl 'http://localhost:33333/interactions' -d '{"datasets": "tfregulons", "sources": "FOXA1,FOXA2,FOXA3,FOXB1,FOXB2,FOXC1,FOXH1", "genesymbols": "1", "fields": "sources,tfregulons_level"}' -H 'Content-Type: application/json' -X 'POST'  > omnipath_webservice_test_post_json.tsv
