#!/bin/bash

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#
#
#  Tested on OS X 10.11.3
#

CONDA3URL="http://repo.continuum.io/archive/Anaconda3-4.1.1-MacOSX-x86_64.sh"
CONDA3INS="Anaconda3-4.1.1-MacOSX-x86_64.sh"
CONDA2URL="http://repo.continuum.io/archive/Anaconda2-4.1.1-MacOSX-x86_64.sh"
CONDA2INS="Anaconda2-4.1.1-MacOSX-x86_64.sh"
CONDA="conda"
CONDAPIP="pip"
PYPATHURL="http://pypath.omnipathdb.org/releases/latest/pypath-latest.tar.gz"

if [[ $PYMAINVER == "3" ]];
then
    curl -LO $CONDA3URL
    chmod +x $CONDA3INS
    bash ./$CONDA3INS
else
    curl -LO $CONDA2URL
    chmod +x $CONDA2INS
    bash ./$CONDA2INS
fi

$CONDA install -c vgauthier cairo=1.12.18
$CONDA install -c richlewis pycairo=1.10.0
$CONDA install pymysql
$CONDA install graphviz

$CONDAPIP install fabric3
$CONDAPIP install pygraphviz
$CONDAPIP install pysftp
$CONDAPIP install future
$CONDAPIP install bioservices
$CONDAPIP install -i https://pypi.anaconda.org/pypi/simple python-igraph
$CONDAPIP install $PYPATHURL
