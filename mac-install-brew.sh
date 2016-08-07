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

$USAGE="\tUsage: $0 [-h (show help and exit)] [-p <2|3> (Python version)] "\
    "[-t (run tests only)] [-u (uninstall everything)] [-f (do not ask confirmation at uninstall)] "\
    "[-m (uninstall Python modules)] [-b (uninstall HomeBrew and formulas)] [-e (remove environment changes)]\n"
PYMAINVER="2"
INSTALL=true
TESTS=true
UNINSTM=false
UNINSTB=false
UNINSTE=false
UCONFIRM=true

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
        t)
            INSTALL=false
            ;;
        u)
            INSTALL=false
            UNINSTM=true
            UNINSTB=true
            UNINSTE=true
            TESTS=false
            ;;
        m)
            INSTALL=false
            UNINSTM=true
            TESTS=false
            ;;
        b)
            INSTALL=false
            UNINSTB=true
            TESTS=false
            ;;
        e)
            INSTALL=false
            UNINSTE=true
            TESTS=false
            ;;
        f)
            UCONFIRM=false
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
        PYFABRIC="fabric3";
    else
        PYVER="2.7";
        PYMAINVER="2";
        PYCAIRONAME="py2cairo";
        PYTHONNAME="python";
        PYFABRIC="fabric"
fi

USER=`whoami`
LOCAL="$HOME/local"
LOCALBIN="$LOCAL/bin"
LOCALPIP="$LOCALBIN/pip$PYVER"
#LOCALPYPATH="$LOCAL/lib/python$PYVER/site-packages"
PYPATHURL="http://pypath.omnipathdb.org/releases/latest/pypath-latest.tar.gz"

if [ ! -d $LOCAL ];
    then mkdir $LOCAL;
fi

export PATH="$LOCALBIN:$PATH"

# beginning part `INSTALL`
if [[ "$INSTALL" = "true" ]];
then

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
    if [ ! `type brew >/dev/null 2>&1` ];
    then
        cd $LOCAL
        curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1
        cd ~
    fi

    brew update
    brew install $PYTHONNAME $PYCAIRONAME homebrew/science/igraph graphviz

    $LOCALPIP install --upgrade pip
    $LOCALPIP install python-igraph
    $LOCALPIP install pysftp
    $LOCALPIP install $PYFABRIC
    $LOCALPIP install future
    $LOCALPIP install pandas
    $LOCALPIP install scipy
    $LOCALPIP install suds-jurko
    $LOCALPIP install bioservices
    $LOCALPIP install pymysql
    $LOCALPIP install pygraphviz

    $LOCALPIP install $PYPATHURL

    # adding local paths and python paths permantently
    cat << EOF >> .bash_profile
    #export PYTHONPATH="$LOCALPYPATH:\$PYTHONPATH"
    export PATH="$LOCALBIN:\$PATH"
EOF

    # with Py3 this will be necessary:
    if [[ $PYMAINVER == "3" ]];
        then $LOCALPIP install git+https://github.com/brentp/fishers_exact_test.git;
    fi

# end of part `INSTALL`
fi

# beginning part `TESTS`

if [[ "$TESTS" == "true" ]];
then
    # check and report:

    echo -en "\n\n====[ Testing installation ]====\n\n"

    if [ ! `type brew >/dev/null 2>&1` ];
    then
        echo -en "\t [!!] Failed to install HomeBrew. This prevents the installation of essential dependencies.\n"
    else
        echo -en "\t [OK] HomeBrew is available.\n"
    fi

    declare -a formulas=($PYTHONNAME $PYCAIRONAME homebrew/science/igraph graphviz)

    for frm in "${modules[@]}"
    do
    hasfrm=`brew ls --versions $frm`
    if [[ $hasfrm"o" == "o" ]];
    then
        echo -en "\t [!!] Brew formula `$frm` is missing.\n"
    else
        echo -en "\t [OK] Brew formula `$frm` is available.\n"
    fi
    done

    declare -a modules=(cairo igraph future numpy scipy pandas suds bioservices pymysql pygraphviz fisher pysftp $PYFABRIC pypath)

    for mod in "${modules[@]}"
    do
        $PYTHONNAME -c "import $mod" >/dev/null 2>&1
        if [[ $? == 0 ]];
        then
            echo -en "\t [OK] Module `$mod` for Python $PYVER is available.\n"
        else
            echo -en "\t [!!] Module `$mod` for Python $PYVER failed to install.\n"
        fi
    done

    $PYTHONNAME -c "import pypath; pa = pypath.PyPath()" >/dev/null 2>&1
    if [[ $? == 0 ]];
    then
        echo -en "\t [OK] Congratulations! You have pypath installed! :)\n"
    else
        echo -en "\t [!!] You have no pypath installed or some issue avoids the module to load.\n"\
            "\t      Please check the list above for failed items. Contact omnipath@googlegroups.com if you need help.\n"
    fi

# end of part `TESTS`
fi

# beginning part `UNINSTALL`

declare -a modules=(pypath $PYFABRIC pysftp fisher pygraphviz pymysql bioservices suds pandas scipy numpy future igraph cairo)
mcols=`for $(seq 0 $((${#modules[@]} - 1))); do echo "($i) ${modules[$i]}"; done; echo "(a) all"; echo "() none"`
mcols=`echo $mcols | columns`
mmsg=`echo -en "===[ Remove the following Python $PYVER modules? ]===\n$mcols\n[Answer e.g. \"1,2,3\" or \"a\"] "`

declare -a formulas=(graphviz homebrew/science/igraph $PYCAIRONAME $PYTHONNAME)
fcols=`for $(seq 0 $((${#modules[@]} - 1))); do echo "($i) ${formulas[$i]}"; done; echo "(a) all"; echo "() none"`
fcols=`echo $fcols | columns`
fmsg=`echo -en "===[ Remove the following HomeBrew formulas? ]===\n$fcols\n[Answer e.g. \"1,2,3\" or \"a\"] "`

# uninstalling python modules:
if [[ "$UNINSTM" == "true" ]];
then
    declare -a msel=(a)
    if [[ "$UCONFIRM" == "true" ]];
    then
        (IFS=","; read -p "$mmsg" -r -a msel)
    fi
    if [[ ${#msel[@]} -gt 0 ]];
    then
        if [[ ${msel[0]} == "a" ]];
        then
            msel=$(seq 0 $((${#modules[@] - 1})))
        else
            for i in ${msel[@]};
            do
                if [[ $i -lt ${#modules[@]} ]];
                then
                    $LOCALPIP --uninstall ${modules[$i]}
                fi
            done
        fi
    fi
fi

# uninstalling homebrew formulas:
if [[ "$UNINSTF" == "true" ]];
then
    declare -a fsel=(a)
    if [[ "$UCONFIRM" == "true" ]];
    then
        (IFS=","; read -p "$fmsg" -r -a fsel)
    fi
    if [[ ${#fsel[@]} -gt 0 ]];
    then
        if [[ ${fsel[0]} == "a" ]];
        then
            fsel=$(seq 0 $((${#formulas[@] - 1})))
        else
            for i in ${fsel[@]};
            do
                if [[ $i -lt ${#formulas[@]} ]];
                then
                    brew uninstall ${formulas[$i]}
                fi
            done
        fi
    fi
fi

# end of part `UNINSTALL`
