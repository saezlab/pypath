#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
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

import pypath.resources._network as _netres
import pypath.resources as _resources
import pypath.internals.resource as _resfmt
from pypath.resources._network import choose_dataset, dorothea_expand_levels

_co = _resources.get_controller()

for _dataset_label in dir(_netres):

    _dataset = getattr(_netres, _dataset_label)

    if not isinstance(_dataset, _resfmt.NetworkDataset):

        continue

    for _label, _resource in _dataset.items():

        _resource.resource_attrs['license'] = _co.license(_resource.name)

    globals()[_dataset_label] = _dataset
