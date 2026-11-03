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

from typing import Dict, List, Optional, Union

import glom

import pypath.share.session as session
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.inputs.ebi as ebi

_logger = session.Logger(name = 'ontology_input')


def ontology(
        ontology: str,
        fields: Optional[Union[List, str, dict]] = None,
    ) -> Union[List[tuple], Dict[str, str]]:
    """
    Downloads an ontology using the EBI Ontology Lookup Service.

    Args
        ontology: The ID of an ontology available in EBI OLS. For a
            full list, call ``listof_ontologies``.
        fields: Additional fields to include, apart from the default
            ones. Can be a dict of glom specs or simply keys in the
            terms section of the OLS Terms query response.
    """

    _logger._log('Retrieving ontology `%s` from EBI OLS.' % ontology)

    url = urls.urls['ols']['url'] + '/%s/terms' % ontology.lower()

    _fields = {
        'label': ('_embedded.terms', ['label']),
        'obo_id': ('_embedded.terms', ['obo_id']),
    }

    if not isinstance(fields, dict):

        fields = dict(
            (
                f.rsplit('.', maxsplit = 1)[-1],
                ('_embedded.terms', [glom.Coalesce(f, default = None)])
            )
            for f in common.to_list(fields)
        )

    _fields.update(fields)

    result = ebi.ebi_rest(
        url = url,
        fields = _fields,
        page_param = 'page',
        page_field = 'page.number',
        paginate = True,
    )

    if not fields:

        result = dict(
            (i.obo_id, i.label)
            for i in result
        )
        result.pop(None, None)

    return result


def listof_ontologies(
        fields: Optional[Union[List, str]] = None,
        full_config: bool = False,
    ) -> Union[List[tuple], Dict[str, str]]:
    """
    Returns a list of available ontologies in the EBI Ontology Lookup Service.

    Args
        fields: Additional field(s) to include, apart from the default
            ones.
        full_config: Keep the full config dict.

    Return
        If no extra fields and no full config are requested: a dict with
        ontology abbreviations as keys and ontology names as values.
        Otherwise a list of named tuples, each representing an ontology
        described by the requested fields.
    """

    _fields = {
        'config': ('_embedded.ontologies', ['config']),
        'ontologyId': ('_embedded.ontologies', ['ontologyId']),
    }

    _fields.update(
        dict(
            (f, ('_embedded.ontologies', [f]))
            for f in common.to_list(fields)
        )
    )

    result = ebi.ebi_rest(
        url = urls.urls['ols']['url'],
        fields = _fields,
        page_param = 'page',
        page_field = 'page.number',
        paginate = True,
    )

    if not fields and not full_config:

        result = dict(
            (i.ontologyId, i.config['title'])
            for i in result
        )
        result.pop(None, None)

    return result
