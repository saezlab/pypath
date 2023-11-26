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

from __future__ import annotations

import copy

import pypath.share.settings as settings
import pypath.share.session as session
import pypath_common._constants as _const
import pypath.inputs.uniprot as uniprot_input
import pypath.inputs.unichem as unichem_input

_logger = session.Logger(name = 'input_formats')

__all__ = [
    'NetworkInput',
    'ReadList',
]


class NetworkInput:


    def __init__(
            self,
            name = "unknown",
            separator = None,
            id_col_a = 0,
            id_col_b = 1,
            id_type_a = "uniprot",
            id_type_b = "uniprot",
            entity_type_a = "protein",
            entity_type_b = "protein",
            is_directed = False,
            sign = False,
            input = None,
            references = None,
            extra_edge_attrs = None,
            extra_node_attrs_a = None,
            extra_node_attrs_b = None,
            header = False,
            taxon_a = 9606,
            taxon_b = 9606,
            ncbi_tax_id = 9606,
            interaction_type = 'post_translational',
            positive_filters = None,
            negative_filters = None,
            mark_source  =  None,
            mark_target  =  None,
            input_args = None,
            curl_args = None,
            must_have_references = True,
            huge = False,
            resource = None,
            unique_fields = None,
            expand_complexes = None,
            data_model = None,
            allow_loops = None,
            only_default_organism = False,
            dataset = None,
        ):
        """
        :param str mark_source:
            Creates a boolean vertex attribute and sets it True for the
            source vertex of directed interactions from this particular
            resource.
        :param str mark_target:
            Same as ``mark_source`` but for target vertices.
        """

        self.entity_type_a = entity_type_a
        self.entity_type_b = entity_type_b
        self.id_col_a = id_col_a
        self.id_col_b = id_col_b
        self.id_type_a = id_type_a
        self.id_type_b = id_type_b
        self.is_directed = is_directed
        self.input = input
        self.extra_edge_attrs = extra_edge_attrs or {}
        self.extra_node_attrs_a = extra_node_attrs_a or {}
        self.extra_node_attrs_b = extra_node_attrs_b or {}
        self.name = name
        self.separator = separator
        self.header = header
        self.refs = references or None
        self.sign = sign
        self.taxon_a = taxon_a
        self.taxon_b = taxon_b
        self.ncbi_tax_id = ncbi_tax_id
        self.interaction_type = interaction_type
        self.positive_filters = positive_filters or []
        self.negative_filters = negative_filters or []
        self.input_args = input_args or {}
        self.curl_args = curl_args or {}
        self.must_have_references = must_have_references and bool(references)
        self.huge = huge
        self.resource = self.name if resource is None else resource
        self.mark_source = mark_source
        self.mark_target = mark_target
        self.unique_fields = unique_fields or set()
        self.expand_complexes = expand_complexes
        self.data_model = data_model
        self.allow_loops = allow_loops
        self.only_default_organism = only_default_organism
        self.dataset = dataset


    def _field(self, value, cls):

        return value if isinstance(value, cls) else cls(compact = value)


class ReadList:


    def __init__(
            self,
            name = 'unknown',
            separator = None,
            id_col = 0,
            id_type = 'uniprot',
            entity_type = 'protein',
            input = None,
            extra_attrs = None,
            header = False,
        ):

        self.entity_type = entity_type
        self.id_col = id_col
        self.id_type = id_type
        self.input = input
        self.extra_attrs = extra_attrs or {}
        self.name = name
        self.separator = separator
        self.header = header
