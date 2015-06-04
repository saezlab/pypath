#!/bin/bash

#
#  This file is part of the `bioigraph` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
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
# thanks for http://user.astro.columbia.edu/~williams/mac/python.html
# and https://truongtx.me/2013/11/08/macports-from-home-directory/
#

#

#     This script attempts to install the bioigraph module with all its
#     dependencies on Mac OS X systems. If everything goes well it does
#     not need user intervention. This is how it works:
# 
#     1: Creates local installation root under $HOME/local
#     2: Adds $HOME/local/bin to $PATH, and the local python dirs 
#        to $PYTHONPATH
#     3: Downloads and installs MacPorts from macports.org
#     4: Installs Python 2.7 (py27) and pycairo (py27-cairo) by port
#     5: Installs igraph and other dependencies by easy_install into the
#        local python package dir
#     4: Finally does the same with bioigraph
# 
#     If you already have Python 2.7 installed, you can load all the 
#     newly installed packages until the local dir is in your $PYTHONPATH.
#     The script attempts to add the path permanently by appending it to
#     ~/.bashrc. 
#     If you have Python <= 2.6 installed, bioigraph most probably won't
#     work under that version. Although as Python 2.7 is installed by
#     port, that will be available on your system.
# 
#     This method is designed to work on most of the various OS X systems,
#     and intends not to have any prerequisite.
#     Depending on your specific environment, simpler, faster and more
#     disk space saving methods might be possible.

PYVER=`python --version | sed 's/.*\([0-9]\.[0-9]\).*/\1/p'`
BIOIGRAPHSRC="http://www.ebi.ac.uk/~denes/54b510889336eb2591d8beff/bioigraph-0.1.25.tar.gz"
BUILDDIR="$HOME/build"
LOCAL="$HOME/local"
LOCALBIN="$LOCAL/bin"
LOCPYHOME="$LOCAL/lib/python2.7"
LOCPYPATH="$LOCPYHOME/site-packages"
MACPORTSVER="2.3.3"
export PYTHONPATH="$LOCPYPATH:$PYTHONPATH"
export PYTHONHOME=$LOCPYHOME
export PATH="$LOCALBIN:$PATH"
mkdir -p $LOCALBIN
mkdir -p $LOCPYPATH
mkdir -p $BUILDDIR

cd $BUILDDIR

echo -en "\n\n## paths modified by bioigraph installer ##\n\n" >> ~/.bash_profile

cat >> ~/.bash_profile <<- EOF
if [ -d $LOCALBIN ]; then
    export PATH=$LOCALBIN:"\${PATH}"
fi
EOF

curl -L $BIOIGRAPHSRC > bioigraph.tar.gz
curl -L https://distfiles.macports.org/MacPorts/MacPorts-$MACPORTSVER.tar.bz2 > macports.tar.bz2
tar -vxjf macports.tar.bz2
mv MacPorts-* macports
cd macports
./configure --disable-readline --prefix=$LOCAL --with-install-user=`id -un` --with-install-group="everyone"
make
make install
cd ..
port -v selfupdate
port install curl-ca-bundle
port install gcc5
export CC="$LOCALBIN/gcc-mp-5"
port install curl
port install python27
PY27="$LOCAL/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7"
PY27DIR="$PY27/site-packages"
PY27DYN="$PY27/lib-dynload/"
DYLDLIB="$LOCAL/lib"
export PYTHONPATH="$PY27:$PY27DIR:$PY27DYN:$PYTHONPATH"
export DYLD_FALLBACK_LIBRARY_PATH="$DYLDLIB"
port install py27-cairo
port install py27-fabric

easy_install-2.7 --install-dir=$PY27DIR ipython
easy_install-2.7 --install-dir=$PY27DIR igraph
easy_install-2.7 --install-dir=$PY27DIR pandas
easy_install-2.7 --install-dir=$PY27DIR requests-cache
easy_install-2.7 --install-dir=$PY27DIR bioigraph.tar.gz

cat >> ~/.bash_profile <<- EOF
if [ -d $PY27DIR ]; then
    export PYTHONPATH="$PY27DIR:\$PYTHONPATH"
fi
EOF

cat >> ~/.bash_profile <<- EOF
if [ -d $DYLDLIB ]; then
    export DYLD_LIBRARY_PATH="$DYLDLIB"
fi
EOF

echo -en "\n\n## end: paths modified by bioigraph installer ##\n\n" >> ~/.bash_profile

cd $HOME

exit 0