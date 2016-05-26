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

# installing HomeBrew first:

USER=`whoami`
HOME="/Users/$USER"
LOCAL="$HOME/local"
LOCALBIN="$LOCAL/bin"
LOCALPIP="$LOCALBIN/pip2.7"
PYPATHURL="http://pypath.omnipathdb.org/releases/latest/pypath-latest.tar.gz"

if [ ! -d $LOCAL ];
    then mkdir $LOCAL;
fi

export PATH="$LOCALBIN:$PATH"

cd $LOCAL

# downloading and extracting homebrew:
curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1

brew update
brew install python py2cairo homebrew/science/igraph graphviz

$LOCALPIP install --upgrade pip
$LOCALPIP install python-igraph
$LOCALPIP install pysftp
$LOCALPIP install fabric
$LOCALPIP install pandas
$LOCALPIP install suds-jurko
$LOCALPIP install bioservices
$LOCALPIP install pymysql
$LOCALPIP install pygraphviz

$LOCALPIP install $PYPATHURL
