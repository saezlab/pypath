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

$USAGE="\tUsage: $0 [-h] [-p <2|3> (Python version)]\n"
PYMAINVER="2"

while getopts ":hp:" opt;
do
    case $opt in
        h)
            echo -en $USAGE;
            exit 0
            ;;
        p)
            PYMAINVER=${OPTARG}
            ;;
        ?)
            echo -en $USAGE;
            exit 2
            ;;
    esac
done

if [[ $PYMAINVER == "3" ]];
    then
        PYVER="3.5";
        PYCAIRONAME="py3cairo";
        PYTHONNAME="python3";
    else
        PYVER="2.7";
        PYMAINVER="2";
        PYCAIRONAME="py2cairo";
        PYTHONNAME="python";
fi

USER=`whoami`
HOME="/Users/$USER"
LOCAL="$HOME/local"
LOCALBIN="$LOCAL/bin"
LOCALPIP="$LOCALBIN/pip$PYVER"
LOCALPYPATH="$LOCAL/lib/python$PYVER/site-packages"
PYPATHURL="http://pypath.omnipathdb.org/releases/latest/pypath-latest.tar.gz"

if [ ! -d $LOCAL ];
    then mkdir $LOCAL;
fi

export PATH="$LOCALBIN:$PATH"

# autocompletion is highly useful
# we add it to the .pyhthonrc
# feel free to remove later if u don't like
cat << EOF >> .pythonrc
import readline
import rlcompleter
if 'libedit' in readline.__doc__:
    readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")
EOF

cat << EOF >> .bash_profile
export PYTHONSTARTUP=~/.pythonrc
EOF

# downloading and extracting homebrew:
cd $LOCAL
curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1
cd ~

brew update
brew install $PYTHONNAME $PYCAIRONAME homebrew/science/igraph graphviz

$LOCALPIP install --upgrade pip
$LOCALPIP install python-igraph
$LOCALPIP install pysftp
$LOCALPIP install fabric
$LOCALPIP install pandas
$LOCALPIP install scipy
$LOCALPIP install suds-jurko
$LOCALPIP install bioservices
$LOCALPIP install pymysql
$LOCALPIP install pygraphviz

$LOCALPIP install $PYPATHURL

# adding local paths and python paths permantently
cat << EOF >> .bash_profile
export PYTHONPATH="$LOCALPYPATH:\$PYTHONPATH"
export PATH="$LOCALBIN:\$PATH"
EOF

# with Py3 this will be necessary:
if [[ $PYMAINVER == "3" ]];
    then $LOCALPIP install git+https://github.com/brentp/fishers_exact_test.git;
fi
