#!/bin/bash

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
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

rm *.pyc
rm pypath/*.pyc
rm -r pypath.egg-info

sed -i 's/\([0-9]*\.[0-9]*\.\)\([0-9]*\)/echo \1$\(\(\2+1\)\)/ge' pypath/__version__

python setup.py sdist bdist_wheel
