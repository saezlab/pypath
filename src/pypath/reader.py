#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
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

import os
import itertools
import pandas as pd
import pypath.common as common
import pypath.dataio as dataio
import pypath.curl as curl

FIELDS = {
    'nameColA',
    'nameColB',
    'nameTypeA',
    'nameTypeB',
    'typeA',
    'typeB',
    'isDirected',
    'sign',
    'references',
    'taxonA',
    'taxonB',
    'ncbiTaxId',
    'interactionType',
    'resource'
    }

UNIQUE_FIELDS = {
    'nameColA',
    'nameColB'
}


class ReaderBase(object):


    def __init__(self, settings):
        self.settings = settings


    def __iter__(self):
        for i, row in enumerate(self.resource):
            if i < self.settings.header or not row:
                continue

            if isinstance(row, common.basestring):
                if row[0] == "#":
                    continue

                row = row.strip('\n\r ').split(self.settings.separator)

            yield row


    def setup_resource(self):
        self.input = self.settings.inFile

        if callable(self.input):
            self.resource = self.input(**self.settings.inputArgs)

        elif isinstance(self.input,common.basestring):
            if hasattr(dataio, self.input):
                self.resource = getattr(dataio, self.input)(
                    **self.settings.inputArgs
                )
            elif (
                os.path.exists(self.input) or
                curl.is_url(self.input)
            ):
                c = curl.Curl(self.input, **self.settings.curlArgs)
                self.resource = c.result

        elif hasattr(self.input, '__iter__'):
            self.resource = self.input

        else:
            self.resource = []




class Reader(ReaderBase):
    """
    Reads network data from an iterable according to
    an :py:class:``pypath.input_formats.NetworkInput`` object.
    """
    def __init__(self, settings):
        """
        :arg input_formats.NetworkInput settings:
            NetworkInput instance
        """
        ReaderBase.__init__(self, settings)


    def iter_rows(self):
        for row in ReaderBase.__iter__(self):
            yield Row(row)


    def __iter__(self):
        for row in self.iter_rows():
            self.new_row(row)
            for values in itertools.product(*self.field_processors):
                yield values


    def iter_dicts(self):
        for values in self:
            yield dict(zip(self.fields, values))


    def get_dataframe(self):
        return pd.DataFrame(list(self), columns=self.fields)

    
    def new_row(self, row):
        for fp in self.field_processors:
            fp.new_row(row)
    

    def setup_fields(self):
        self.fields = FIELDS
        self.fields = self.fields & set(self.settings.extraEdgeAttrs.keys())
        self.fields = self.fields & set(self.settings.extraNodeAttrsA.keys())
        self.fields = self.fields & set(self.settings.extraNodeAttrsB.keys())
        self.fields = sorted(self.fields)
        self.unique_fields = UNIQUE_FIELDS
        self.unique_fields = (
            self.unique_fields & set(self.settings.unique_fields)
        )


    def setup_field_processors(self):
        self.field_processors = [

            FieldProcessor(getattr(self.settings, field),
                           field in self.unique_fields)

            for field in self.fields
        ]


class FieldProcessor(object):
    def __init__(self, field, is_unique = False):
        """
        :arg field:
            Field processing definition
        :arg callable method:
            Custome method which accepts a row as an argument and returns
            the processed field value
        """
        self.field = field
        self.single_value = is_unique
        self.failed = False
        self.setup_method()


    def __iter__(self):
        fields = self.process()
        if isinstance(fields, common.simple_types):
            fields = (fields,)
        if self.single_value:
            for field in self.process():
                yield field
        
        else:
            yield set(fields)



    def setup_method(self):

        if isinstance(self._method, common.basestring):
            self._method = self.str_method

        elif callable(self.field):
            self._method = self.field

        elif isinstance(self.field, int):
            self.i = self.field
            self._method = self.index_method

        elif isinstance(self.field, (tuple,list)):
            self._method = self.tuple_method

        elif isinstance(self.field, dict):
            self._method = self.dict_method


    def str_method(self, row = None):
        return self.field


    def index_method(self, row = None):
        return self.row[self.i]


    def tuple_method(self, row = None):
        self.i = self.field[0]
        value = self.index_method()

        if isinstance(self.field[1], common.basestring) and value:
            value = value.split(self.field[1])

        else:
            value = ()

        if len(self.field) > 2:
            value = bool(set(value) & set(self.field[2]))

        return value


    def dict_method(self):
        self.i = self.field['col']
        value = self.index_method()
        mapping = self.field['dict']

        return mapping[value] if value in mapping else None



    def new_row(self, row):
        self.row = row


    def process(self):
        self.process_common()
        return self._method(self.row)


    def process_common(self):
        pass




class Row(object):
    __slots__ = ["row", "failed"]

    def __init__(self, row):
        self.row = row
        self.failed = False


    def __getitem__(self, i):
        """
        :arg int i:
            Index of the field
        """
        try:
            return self.row[i]

        except IndexError:
            self.failed = True
            return ""





class EdgeAttribute(FieldProcessor):
    def __init__(self,
                 separator=None,
                 isDirected=False,
                 sign=False,
                 references=False,
                 must_have_references=True,
                 extraEdgeAttrs={},
                 interactionType='PPI',
                 positiveFilters=[],
                 negativeFilters=[],
                 inFile=None,
                 header=False
                 ):

        Reader.__init__(self,
                 separator=None,
                 isDirected=False,
                 sign=False,
                 references=False,
                 must_have_references=True,
                 extraEdgeAttrs={},
                 interactionType='PPI',
                 positiveFilters=[],
                 negativeFilters=[],
                 inFile = None,
                 header = False
                        )


        self.isDirected = isDirected
        self.extraEdgeAttrs = extraEdgeAttrs
        self.separator = separator
        self.refs = references
        self.must_have_references = must_have_references and bool(references)
        self.sign = sign
        self.intType = interactionType
        self.positiveFilters = positiveFilters
        self.negativeFilters = negativeFilters
        self.inFile = inFile
        self.header = header



class NodeAttribute(FieldProcessor):
    def __init__(self,
                 name="unknown",
                 separator=None,
                 nameColA=0,
                 nameColB=1,
                 nameTypeA="uniprot",
                 nameTypeB="uniprot",
                 typeA="protein",
                 typeB="protein",
                 inFile=None,
                 extraNodeAttrsA={},
                 extraNodeAttrsB={},
                 header=False,
                 taxonA=9606,
                 taxonB=9606,
                 ncbiTaxId=False,
                 mark_source=None,
                 mark_target=None,
                 inputArgs={},
                 huge=False,
                 resource=None,
                 positiveFilters=[],
                 negativeFilters=[]
                 ):

        self.typeA = typeA
        self.typeB = typeB
        self.nameColA = nameColA
        self.nameColB = nameColB
        self.nameTypeA = nameTypeA
        self.nameTypeB = nameTypeB
        self.extraNodeAttrsA = extraNodeAttrsA
        self.extraNodeAttrsB = extraNodeAttrsB
        self.name = name
        self.separator = separator
        self.taxonA = taxonA
        self.taxonB = taxonB
        self.ncbiTaxId = ncbiTaxId
        self.positiveFilters = positiveFilters
        self.negativeFilters = negativeFilters
        self.inputArgs = inputArgs
        self.huge = huge
        self.resource = self.name if resource is None else resource
        self.mark_source = mark_source
        self.mark_target = mark_target
        self.inFile = inFile
        self.header = header


class MergeTable():
    def __init__(self):
        pass


class CreateIgraphObject():
    def __init__(self):
        pass
