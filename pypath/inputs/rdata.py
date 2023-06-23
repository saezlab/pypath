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

import rdata
import numpy as np


def _rdata_data_frame_get_rownames(robj):

    for i, attr in enumerate(robj.attributes.value):

        if (
            attr.value and
            attr.value[1] and (
                (
                    hasattr(attr.value[1], 'tag') and
                    attr.value[1].tag and (
                        (
                            attr.value[1].tag.referenced_object and
                            attr.value[1].tag.referenced_object.value and
                            attr.value[1].tag.referenced_object.value.\
                                value == b'row.names'
                        ) or (
                            attr.value[1].tag.value and
                            attr.value[1].tag.value.value == b'row.names'
                        )
                    )
                ) or (
                    attr.value[1].tag and
                    attr.value[1].tag.referenced_object and
                    attr.value[1].tag.referenced_object.value and (
                        attr.value[1].tag.referenced_object.value.value ==
                        b'row.names'
                    )
                ) or (
                    attr.value[1].tag and
                    attr.value[1].tag.value and
                    attr.value[1].tag.value.value == b'row.names'
                )
            )
        ):

            break

    rownames = (
        attr.value[1].value[0].value
            if attr.value[0].value[0].value == b'data.frame' else
        attr.value[0].value
            if (
                attr.value[1].value[0].value[0].value ==
                b'data.frame'
            ) else
        []
    )

    return [rn.value.decode('utf-8') for rn in rownames]


def _rdata_list_get_names(robj):

    return [
        item.value.decode('utf-8')
        for item in robj.attributes.value[0].value
    ]


def _patch_rdata():

    def parse_R_object(self, reference_list=None):
        """
        Parse a R object.
        """

        if reference_list is None:
            # Index is 1-based, so we insert a dummy object
            reference_list = [None]

        info_int = self.parse_int()

        info = rdata.parser._parser.parse_r_object_info(info_int)

        tag = None
        attributes = None
        referenced_object = None

        tag_read = False
        attributes_read = False
        add_reference = False

        if info.type == rdata.parser._parser.RObjectType.SYM:
            # Read Char
            value = self.parse_R_object(reference_list)
            # Symbols can be referenced
            add_reference = True

        elif info.type in [
            rdata.parser._parser.RObjectType.LIST,
            rdata.parser._parser.RObjectType.LANG
        ]:

            tag = None
            if info.attributes:
                raise NotImplementedError('Attributes not suported for LIST')
            elif info.tag:
                tag = self.parse_R_object(reference_list)
                tag_read = True

            # Read CAR and CDR
            car = self.parse_R_object(reference_list)
            cdr = self.parse_R_object(reference_list)
            value = (car, cdr)

        elif info.type == rdata.parser._parser.RObjectType.CHAR:

            length = self.parse_int()
            if length > 0:
                value = self.parse_string(length=length)
            else:
                value = b''

        elif info.type == rdata.parser._parser.RObjectType.LGL:

            length = self.parse_int()

            value = np.empty(length, dtype=rdata.parser._parser.np.bool_)

            for i in range(length):
                value[i] = self.parse_bool()

        elif info.type == rdata.parser._parser.RObjectType.INT:
            length = self.parse_int()

            value = rdata.parser._parser.np.empty(
                length,
                dtype=rdata.parser._parser.np.int64
            )

            for i in range(length):
                value[i] = self.parse_int()

        elif info.type == rdata.parser._parser.RObjectType.REAL:
            length = self.parse_int()

            value = np.empty(length, dtype=rdata.parser._parser.np.double)

            for i in range(length):
                value[i] = self.parse_double()

        elif info.type == rdata.parser._parser.RObjectType.CPLX:
            length = self.parse_int()

            value = np.empty(length, dtype=rdata.parser._parser.np.complex_)

            for i in range(length):
                value[i] = self.parse_complex()

        elif info.type in [
                rdata.parser._parser.RObjectType.STR,
                rdata.parser._parser.RObjectType.VEC,
                rdata.parser._parser.RObjectType.EXPR
            ]:
            length = self.parse_int()

            value = [None] * length

            for i in range(length):
                value[i] = self.parse_R_object(reference_list)

        elif info.type == rdata.parser._parser.RObjectType.NILVALUE:
            value = None

        elif info.type == rdata.parser._parser.RObjectType.REF:
            value = None
            referenced_object = reference_list[info.reference]

        else:
            raise NotImplementedError(f'Type {info.type} not implemented')

        if info.tag and not tag_read:
            rdata.parser._parser.warnings.warn(
                f'Tag not implemented for type {info.type} and ignored'
            )
        if info.attributes and not attributes_read:
            attributes = self.parse_R_object(reference_list)

        result = rdata.parser._parser.RObject(
            info=info,
            tag=tag,
            attributes=attributes,
            value=value,
            referenced_object=referenced_object,
        )

        if add_reference:
            reference_list.append(result)

        return result

    rdata.parser._parser.Parser.parse_R_object = parse_R_object


_patch_rdata()
