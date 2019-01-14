#!/bin/bash

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2018
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolas Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#
#
#  Tested on OS X 10.11.3
#

# installing HomeBrew first:

USAGE="Usage:\n\t$0\n\t\t[-h (show help and exit)]\n\t\t[-p <2|3> (Python version)]\n\t\t"\
"[-c (do not install cairo)]\n\t\t[-g (do not install graphviz)]\n\t\t"\
"[-t (run tests only)]\n\t\t[-u (uninstall everything)]\n\t\t[-f (do not ask confirmation at uninstall)]\n\t\t"\
"[-m (uninstall Python modules)]\n\t\t[-b (uninstall HomeBrew and formulas)]\n\t\t[-e (remove environment changes)]\n"
PYMAINVER="2"
INSTALL=true
TESTS=true
ICAIRO=true
IGRAPHVIZ=true
UNINSTM=false
UNINSTB=false
UNINSTE=false
UCONFIRM=true

while getopts ":htumbefp:" opt;
do
    case $opt in
        h)
            echo -en "$USAGE";
            exit 0
            ;;
        p)
            PYMAINVER="${OPTARG}"
            ;;
        t)
            INSTALL=false
            ;;
        u)
            INSTALL=false
            UNINSTM=true
            UNINSTF=true
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
            UNINSTF=true
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
        c)
            ICAIRO=false
            ;;
        g)
            IGRAPHVIZ=false
            ;;
        ?)
            echo -en "$USAGE";
            exit 2
            ;;
    esac
done

if [[ $PYMAINVER == "3" ]];
    then
        PYVER="3.7";
        PYCAIRONAME="py3cairo";
        PYTHONNAME="python3";
        PYFABRIC="fabric3";
    else
        PYVER="2.7";
        PYMAINVER="2";
        PYCAIRONAME="py2cairo";
        PYTHONNAME="python";
        PYFABRIC="fabric";
fi

USER=`whoami`
LOCAL="$HOME/local"
LOCALBIN="$LOCAL/bin"
LOCALPIP="$LOCALBIN/pip$PYVER"
PYPATHURL="http://pypath.omnipathdb.org/releases/latest/pypath-latest.tar.gz"

if [ ! -d $LOCAL ];
    then mkdir $LOCAL;
fi

export PATH="$LOCALBIN:$PATH"

# beginning part `INSTALL`

PYTHONRCCONTENT=''\
'import readline# pypath added\n'\
'import rlcompleter# pypath added\n'\
'if "libedit" in readline.__doc__:# pypath added\n'\
'    readline.parse_and_bind("bind ^I rl_complete")# pypath added\n'\
'else:# pypath added\n'\
'    readline.parse_and_bind("tab: complete")# pypath added'
BASHPROFPYRC='export PYTHONSTARTUP=~/.pythonrc # pypath added'
BASHPROFLOCP='export PATH="'$LOCALBIN':$PATH" # pypath added'

if [[ "$INSTALL" = "true" ]];
then
    echo -en "\n\n===[ Attempting to install pypath and all its dependencies with the help of HomeBrew. ]===\n\n"
    echo -en "\t Note: this method works on most of the Mac computers.\n"\
"\t Watch out for errors, and the test results post installation.\n"\
"\t This will last at least 10 mins. Now relax, and hope the best.\n\n"
    if [[ ! -f .pythonrc || "$(grep 'tab:[[:space:]]\?complete' .pythonrc)" == "" ]];
    then
        # autocompletion is highly useful
        # we add it to the .pyhthonrc
        # feel free to remove later if u don't like
        echo -en "$PYTHONRCCONTENT" >> .pythonrc
    fi

    if [[ ! -f .bash_profile || "$(grep '.pythonrc' .bash_profile)" == "" ]];
    then
        # adding this to .bash_profile to use .pythonrc
        echo -en "$BASHPROFPYRC" >> .bash_profile
    fi
    if [[ ! -f .bash_profile || "$(grep $LOCALBIN .bash_profile)" == "" ]];
    then
        # adding this to .bash_profile to have local in path
        echo -en "$BASHPROFLOCP" >> .bash_profile
    fi

    # downloading and extracting homebrew:
    type brew >/dev/null 2>&1
    if [[ $? != 0 ]];
    then
        cd $LOCAL
        curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1
        cd ~
    fi
    
    # update repositories
    brew update
    # installing curl and openssl
    brew install curl
    brew install openssl
    brew install --force libxml2
    brew install --force libxslt
    brew link libxml2 --force
    brew link libxslt --force
    # obtaining a recent python distribution from brew:
    brew install $PYTHONNAME
    # optionally install (py)cairo:
    if [[ "$ICAIRO" == "true" ]];
        then brew install $PYCAIRONAME;
    fi
    # optionally install graphviz:
    if [[ "$IGRAPHVIZ" == "true" ]];
        then brew install homebrew/science/igraph graphviz;
    fi
    
    # maybe we just installed new python, update the python/pip version
    PYVER=$(
        [[ `$PYTHONNAME --version` =~ (3\.[0-9])\.[0-9] ]] &&
        echo ${BASH_REMATCH[1]}
    )
    LOCALPIP="$LOCALBIN/pip$PYVER"
    
    # installing another python modules by pip
    $LOCALPIP install --upgrade pip
    if [[ `$LOCALPIP show pycurl` ]];
        then $LOCALPIP uninstall pycurl
    fi
    if [[ `$LOCALPIP show lxml` ]];
        then $LOCALPIP uninstall lxml
    fi
    export LD_LIBRARY_PATH=/usr/local/opt/curl/lib
    export CPATH=/usr/local/opt/openssl/include
    export LIBRARY_PATH=/usr/local/opt/openssl/lib
    $LOCALPIP --no-cache-dir install pycurl
    STATIC_DEPS=true $LOCALPIP --no-cache-dir install lxml
    $LOCALPIP install python-igraph
    $LOCALPIP install pysftp
    $LOCALPIP install $PYFABRIC
    $LOCALPIP install future
    $LOCALPIP install pandas
    $LOCALPIP install scipy
    $LOCALPIP install suds-jurko
    $LOCALPIP install bioservices
    $LOCALPIP install pymysql
    $LOCALPIP install tqdm
    # optionally install pygraphviz
    if [[ "$IGRAPHVIZ" == "true" ]];
        then $LOCALPIP install pygraphviz;
    fi

    $LOCALPIP install $PYPATHURL

    # adding local paths and python paths permantently
    cat << EOF >> .bash_profile
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

    echo -en "\n\n\t===[ Testing installation ]===\n\n"

    type brew >/dev/null 2>&1

    if [[ $? -eq 1 ]];
    then
        echo -en "\t [ !! ] Failed to install HomeBrew. This prevents the installation of essential dependencies.\n"
    else
        echo -en "\t [ OK ] HomeBrew is available.\n"
    fi

    declare -a formulas=(curl openssl $PYTHONNAME $PYCAIRONAME homebrew/science/igraph graphviz)

    for frm in "${formulas[@]}"
    do
        hasfrm="$(brew ls --versions $(echo $frm | awk 'BEGIN{FS="/"}{print $NF}'))"
    if [[ $hasfrm"o" == "o" ]];
    then
        echo -en "\t [ !! ] Brew formula \`$frm\` is missing.\n"
    else
        echo -en "\t [ OK ] Brew formula \`$frm\` is available.\n"
    fi
    done

    declare -a modules=(cairo igraph future numpy scipy pandas suds bioservices pymysql pygraphviz fisher pysftp fabric tqdm pypath)

    for mod in "${modules[@]}"
    do
        $PYTHONNAME -c "import $mod" >/dev/null 2>&1
        if [[ $? == 0 ]];
        then
            echo -en "\t [ OK ] Module \`$mod\` for Python $PYVER is available.\n"
        else
            echo -en "\t [ !! ] Module \`$mod\` for Python $PYVER failed to install.\n"
        fi
    done

    $PYTHONNAME -c "import pypath; pa = pypath.PyPath()" >/dev/null 2>&1
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

# beginning part `UNINSTALL`

declare -a modules=(pypath $PYFABRIC pysftp fisher pygraphviz pymysql bioservices suds pandas scipy numpy future igraph cairo tqdm)
mcols="$(for i in $(seq 0 $((${#modules[@]} - 1))); do echo "($i) ${modules[$i]}"; done; echo "(a) all"; echo "() none")"
mcols="$(echo $mcols | column)"
mmsg="$(echo -en "\n===>[ Remove the following Python $PYVER modules? ]<===\n\n$mcols\n[Answer e.g. \"1,2,3\" or \"a\"] ")"

declare -a formulas=(graphviz homebrew/science/igraph $PYCAIRONAME $PYTHONNAME)
fcols="$(for i in $(seq 0 $((${#formulas[@]} - 1))); do echo "($i) ${formulas[$i]}"; done; echo "(a) all"; echo "() none")"
fcols="$(echo $fcols | column)"
fmsg="$(echo -en "\n===>[ Remove the following HomeBrew formulas? ]<===\n\n$fcols\n[Answer e.g. \"1,2,3\" or \"a\"] ")"

bmsg="$(echo -en "\n===>[ Do you want to remove HomeBrew? ]<===\n\n [ y/N ] ")"

# uninstalling python modules:
if [[ "$UNINSTM" == "true" ]];
then

# make sure python/pip version is correct
PYVER=$(
    [[ `$PYTHONNAME --version` =~ (3\.[0-9])\.[0-9] ]] &&
    echo ${BASH_REMATCH[1]}
)
LOCALPIP="$LOCALBIN/pip$PYVER"

    if [[ "$UCONFIRM" == "true" ]];
    then
        ifs_default=$IFS
        IFS=","
        read -p "$mmsg" -r -a msel
        IFS=$ifs_default
    else
        declare -a msel=(a)
    fi
    if [[ ${#msel[@]} -gt 0 ]];
    then
        echo -en "\n\n\t===[ Uninstalling Python modules ]===\n\n"
        if [[ ${msel[0]} == "a" ]];
        then
            msel=$(seq 0 $((${#modules[@] - 1})))
        fi
        for i in ${msel[@]};
        do
            if [[ $i -lt ${#modules[@]} ]];
            then
                $LOCALPIP uninstall -y ${modules[$i]}
            fi
        done
    fi
fi

# uninstalling homebrew formulas:
if [[ "$UNINSTF" == "true" ]];
then
    if [[ "$UCONFIRM" == "true" ]];
    then
        ifs_default=$IFS
        IFS=","
        read -p "$fmsg" -r -a fsel
        IFS=$ifs_default
    else
        declare -a fsel=(a)
    fi
    if [[ ${#fsel[@]} -gt 0 ]];
    then
        if [[ ${fsel[0]} == "a" ]];
        then
            fsel=$(seq 0 $((${#formulas[@] - 1})))
        fi
        echo -en "\n\n===[ Uninstalling HomeBrew formulas ]===\n\n"
        for i in ${fsel[@]};
        do
            if [[ $i -lt ${#formulas[@]} ]];
            then
                brew remove --force $(echo ${formulas[$i]} | awk 'BEGIN{FS="/"}{print $NF}')
                brew rm $(join <(brew leaves) <(brew deps ${formulas[$i]}))
            fi
        done
    fi
    if [[ "$UCONFIRM" == "true" ]];
    then
        read -p "$bmsg" -n 1 UNINSTB
        if [[ "$UNINSTB" =~ ^(y|Y)$ ]];
        then
            UNINSTB=true
        else
            UNINSTB=false
        fi
    fi
    if [[ "$UNINSTB" == "true" ]];
    then
        echo -en "\n===[ Uninstalling HomeBrew ]===\n\n"
        curl -L 'https://raw.githubusercontent.com/Homebrew/install/master/uninstall' > brew-uninstall.rb
        chmod +x brew-uninstall.rb
        ./brew-uninstall.rb -qf
    fi
fi

if [[ "$UNINSTE" == "true" ]];
then
    echo -en "\n\n===[ Restoring environment (~/.pythonrc, ~/.bash_profile) ]===\n\n"
    if [ -f .pythonrc ];
    then
        echo "$(tr '\n' '\a' < .pythonrc)" | sed 's/'"$(echo -en "$PYTHONRCCONTENT" | tr '\n' '\a')"'//g' | tr '\a' '\n' > .pythonrc_tmp
        mv .pythonrc_tmp .pythonrc
    fi
    if [ -f .bash_profile ];
    then
        BASHPROFPYRCE="$(echo $BASHPROFPYRC | sed -e 's/[]\/$*.^|[]/\\&/g')"
        echo "$(tr '\n' '\a' < .bash_profile)" | sed 's/'"$(echo -en "$BASHPROFPYRCE" | tr '\n' '\a')"'//g' | tr '\a' '\n' > .bash_profile_tmp
        mv .bash_profile_tmp .bash_profile
        BASHPROFLOCPE="$(echo $BASHPROFLOCP | sed -e 's/[]\/$*.^|[]/\\&/g')"
        echo "$(tr '\n' '\a' < .bash_profile)" | sed 's/'"$(echo -en "$BASHPROFLOCPE" | tr '\n' '\a')"'//g' | tr '\a' '\n' > .bash_profile_tmp
        mv .bash_profile_tmp .bash_profile
    fi
fi

# end of part `UNINSTALL`

unset LOCALBIN
unset LOCALPIP
unset PYPATHURL

echo -en "\n===[ Tasks complete. Please report any issue to omnipath@googlegroups.com. Bye. ]===\n\n"
