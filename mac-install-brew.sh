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
BREWDIR="$LOCAL/brew"

if ! -d $LOCAL;
    then mkdir $LOCAL;
fi

mkdir $BREWDIR
cd $BREWDIR

# downloading and extracting homebrew:
curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
