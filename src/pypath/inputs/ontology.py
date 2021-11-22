#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import pypath.share.session as session
import pypath.resources.urls as urls
import pypath.inputs.ebi as ebi

_logger = session.Logger(name = 'ontology_input')


def ontology(ontology, fields = None):
    """
    Downloads an ontology using the EBI Ontology Lookup Service.

    Args:
        ontology (str): The ID of an ontology available in EBI OLS. For a
            full list, call ``listof_ontologies``.
        fields (list): Additional fields to include, apart from the default
            ones.
    """

    _logger._log('Retrieving ontology `%s` from EBI OLS.' % ontology)

    url = urls.urls['ols']['url'] + '/%s/terms' % ontology.lower()
    _fields = sorted(['label', 'obo_id'] + (fields or []))
    result = ebi.ebi_rest(url = url, fields = _fields)

    if not fields:

        result = dict(
            (i.obo_id, i.label)
            for i in result
        )
        result.pop(None, None)

    return result


def listof_ontologies(fields = None, full_config = False):
    """
    Returns a list of available ontologies in the EBI Ontology Lookup Service.

    Args:
        fields (list): Additional fields to include, apart from the default
            ones.
        full_config (bool): Keep the full config dict.
    """

    url = urls.urls['ols']['url']
    _fields = sorted(['ontologyId', 'config'] + (fields or []))
    result = ebi.ebi_rest(url, fields = _fields)

    if not fields and not full_config:

        result = dict(
            (i.ontologyId, i.config['title'])
            for i in result
        )
        result.pop(None, None)

    return result
