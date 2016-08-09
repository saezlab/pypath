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

USAGE="Usage:\n\t$0\n\t\t[-h (show help and exit)]\n\t\t[-p <2|3> (Python version)]\n\t\t"\
"[-t (run tests only)]\n\t\t-c <Anaconda path, e.g. ~/anaconda3>\n\n"
PYMAINVER="3"

INSTALL=true
TESTS=true


while getopts ":htp:c:" opt;
do
    case $opt in
        h)
            echo -en "$USAGE"
            exit 0
            ;;
        t)
            INSTALL=false
            ;;
        c)
            CONDAROOT="${OPTARG}"
            ;;
        p)
            PYMAINVER="${OPTARG}"
            ;;
        ?)
            echo -en "$USAGE"
            exit 2
            ;;
    esac
done

if [ -z "${CONDAROOT+x}" ];
then
    if [[ $PYMAINVER == "3" ]];
    then
        CONDAINS="Anaconda3-4.1.1-MacOSX-x86_64.sh"
        CONDABIN="$HOME/anaconda3/bin"
        CONDAURL="http://repo.continuum.io/archive/Anaconda3-4.1.1-MacOSX-x86_64.sh"
        PYFABRIC="fabric3"
    else
        CONDAINS="Anaconda2-4.1.1-MacOSX-x86_64.sh"
        CONDABIN="$HOME/anaconda2/bin"
        CONDAURL="http://repo.continuum.io/archive/Anaconda2-4.1.1-MacOSX-x86_64.sh"
        PYFABRIC="fabric"
    fi
else
    CONDABIN="$CONDAROOT/bin"
    CONDA="$CONDABIN/conda"
    CONDAPIP="$CONDABIN/pip"
fi 

CONDA="$CONDABIN/conda"
CONDAPIP="$CONDABIN/pip"
CONDAPY="$CONDABIN/python"
PYPATHURL="http://pypath.omnipathdb.org/releases/latest/pypath-latest.tar.gz"

if [[ "$INSTALL" == "true" ]];
then
    echo -en "\n\n===[ Attempting to install pypath and all its dependencies with the help of Anaconda. ]===\n\n"
    echo -en "\t Note: this method works on most of the Mac computers. Watch out for errors, and the test results post installation.\n\t"\
" This will take a couple of minutes. Now relax, and hope the best.\n\n"
    if [ ! -d $CONDABIN ];
    then
        if [ ! -f $CONDAINS ];
        then
            curl -LO $CONDAURL
        fi
        chmod +x $CONDAINS
        bash ./$CONDAINS -b
    fi

    $CONDA install -y -c vgauthier cairo=1.12.18
    $CONDA install -y -c richlewis pycairo=1.10.0
    $CONDA install -y pymysql
    $CONDA install -y graphviz
    $CONDA install -y -c bioconda python-igraph
    $CONDA install -y -c omnia pygraphviz

    if [[ $PYMAINVER == "3" ]];
    then
        $CONDAPIP install git+https://github.com/brentp/fishers_exact_test.git
    fi
    $CONDAPIP install $PYFABRIC
    $CONDAPIP install pysftp
    $CONDAPIP install future
    $CONDAPIP install bioservices
    $CONDAPIP install $PYPATHURL
fi

# beginning part `TESTS`

if [[ "$TESTS" == "true" ]];
then
    # check and report:

    echo -en "\n\n\t===[ Testing installation ]===\n\n"

    declare -a modules=(cairo igraph future numpy scipy pandas suds bioservices pymysql pygraphviz fisher pysftp fabric pypath)

    for mod in "${modules[@]}"
    do
        $CONDAPY -c "import $mod" >/dev/null 2>&1
        if [[ $? == 0 ]];
        then
            echo -en "\t [ OK ] Module \`$mod\` for Python $PYMAINVER is available.\n"
        else
            echo -en "\t [ !! ] Module \`$mod\` for Python $PYMAINVER failed to install.\n"
        fi
    done

    $CONDAPY -c "import pypath; pa = pypath.PyPath()" >/dev/null 2>&1
    if [[ $? == 0 ]];
    then
        echo -en "\t [ OK ] Congratulations! You have pypath installed! :)\n"\
            "\t        The authors wish you interesting findings in your analysis.\n"\
            "\t        If you experience any issue, don't hesitate to contact omnipath@googlegroups.com.\n\n"
    else
        echo -en "\t [ !! ] You have no pypath installed or some issue avoids the module to load.\n"\
            "\t        Please check the list above for failed items. Contact omnipath@googlegroups.com if you need help.\n\n"
    fi

# end of part `TESTS`
fi

echo -en "\n===[ Tasks complete. Please report any issue to omnipath@googlegroups.com. Bye. ]===\n\n"
