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

import pypath.omnipath.databases.define as define


def class_and_param(label):
    """
    Retrieves the database definition of a built in database.
    """
    
    manager = get_manager()
    
    return manager.class_and_param(label)


def build(label):
    """
    Builds a built in database and returns the instance.
    This is not the preferred method to get a database instance.
    Unless there is a strong reason, both built in and user defined databases
    should be managed by the ``pypath.omnipath.app`` module.
    """
    
    manager = get_manager()
    
    return manager.build(label)


def get_manager():
    
    if 'manager' not in globals():
        
        init_manager()
    
    return globals()['manager']


def init_manager():
    
    globals()['manager'] = define.DatabaseDefinitionManager()
