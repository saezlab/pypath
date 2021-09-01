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
#                  Sebastian Lobentanzer
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import pypath.inputs.cancerdrugs_db as cancerdrugs_db

def test_cancerdrugs_db():
    # test annotations
    a = cancerdrugs_db.cancerdrugs_db_annotations()
    b = a.get('CHEMBL3301610')
    c = next(b.__iter__()).label

    assert 100 < len(a) < 1000
    assert c == 'Abemaciclib'

    # test interactions
    d = cancerdrugs_db.cancerdrugs_db_interactions()
    e = next(d.__iter__())

    assert 500 < len(d) < 20000
    # how to test for specific interactions when they are not 1 to 1?