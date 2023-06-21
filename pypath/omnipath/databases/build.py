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

def build(dbclass, dbdef):
    """
    Builds a database following the instructions in a ``DatabaseDefinition``
    object, using the database class or method ``dbclass``.
    
    This is not the preferred method to get a database instance.
    Unless there is a strong reason, both built in and user defined databases
    should be managed by the ``pypath.omnipath.app`` module.
    """
    
    dbclass = dbclass if callable(dbclass) else dbclass.get_class()
    
    build_method = (
        dbclass
            if not dbdef.get('init') else
        getattr(dbclass, dbdef.get('init'))
    )
    
    build_args = dbdef.get('args') or {}
    
    db = build_method(**build_args)
    
    prep = dbdef.get('prepare') or {}
    
    for var, method in prep.items():
        
        locals()[var] = getattr(db, method)()
    
    workflow = dbdef.get('workflow') or {}
    
    for step in workflow:
        
        method = step['method']
        args = dict(
            (
                argname,
                locals()[val] if val in locals() else val
            )
            for argname, val in step['args']
        )
        
        getattr(db, method)(**args)
    
    return db
