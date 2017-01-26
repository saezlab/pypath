#!/bin/bash

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2017 - EMBL-EBI
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

#
# The main difficulty of installing pypath on Mac is that
# it depends on igraph, which depends on pycairo, and
# on cairo consequently. The latter 2 can not be installed
# just by pip, so we need to install them beforehand.
#
# This method relies on the built in Python distributed with
# recent versions of OS X.
#
# Requirements:
#
#     # Python 2.7 installed on your system and available
#       in your path
#     # XCode installed on your system, or being available
#       to give permissions for installation when it
#       asks for it.
#
# We do not use here any package management system, but
# we compile some libraries from source.
#
# The versions used here are the recent ones as of 20/05/2016.
#

#
# For troubleshooting pixman and cairo compilation,
# please see this guide:
# https://github.com/Automattic/node-canvas/wiki/Installation---OSX
#

# setting some variables

echo -en "\n\n\t [!!] Please note: this is not the suggested way to install pypath.\n"\
    "\t      Have you considered \`mac-install-brew.sh\` or \`mac-install-conda.sh\`?\n"\
    "\t      If you want to use this script anyways, please open in your editor,\n"\
    "\t      so you know what are you doing.\n"\
    "\t      Now exiting.\n\n"
exit 0

USER=`whoami`
HOME="/Users/$USER"
LOCAL="$HOME/local"
LOCALBIN="$LOCAL/bin"
BUILDDIR="$HOME/pypath_build"
PYVER=`python --version 2>&1 | sed 's/.*\([0-9]\.[0-9]\.[0-9]\{1,2\}\).*/\1/g'`
PYVERSHORT=`python --version 2>&1 | sed 's/.*\([0-9]\.[0-9]\)\.[0-9].*/\1/g'`

# modify these versions if you find more recent ones:
PYCAIROVER="1.10.0"
CAIROVER="1.14.6"
PXMANVER="0.34.0"
PKGCFGVER="0.28"
LIBPNGVER="1.6.21"
GVIZVER="2.38.0"
AMVER="1.15"
LTOOLVER="2.4.6"

# URLs
PKGCFGURL="https://pkg-config.freedesktop.org/releases/pkg-config-$PKGCFGVER.tar.gz"
PYCAIROURL="http://cairographics.org/releases/py2cairo-$PYCAIROVER.tar.bz2"
PXMANURL="http://www.cairographics.org/releases/pixman-$PXMANVER.tar.gz"
CAIROURL="http://cairographics.org/releases/cairo-$CAIROVER.tar.xz"
LIBPNGURL="http://downloads.sourceforge.net/project/libpng/libpng16/$LIBPNGVER/libpng-$LIBPNGVER.tar.xz"
PYPATHURL="http://pypath.omnipathdb.org/releases/latest/pypath-latest.tar.gz"
# GVIZURL="http://www.graphviz.org/pub/graphviz/stable/SOURCES/graphviz-$GVIZVER.tar.gz"
GVIZURL="https://github.com/ellson/graphviz/archive/master.tar.gz"

# local Python installation paths
LOCALPYHOME="/Users/$USER/Library/Python/$PYVERSHORT"
LOCALPYHOME2="$LOCAL/lib/python2.7"
LOCALPYPATH="$LOCALPYHOME/site-packages"
LOCALPYBIN="$LOCALPYHOME/bin"
LOCALPYPATH2="$LOCALPYHOME2/site-packages"
LOCALPYBIN2="$LOCALPYHOME2/bin"

export PYTHONPATH="$LOCALPYPATH:$LOCALPYPATH2:$PYTHONPATH"
export PATH="$LOCALPYBIN:$LOCALPYBIN2:$LOCALBIN:$PATH"

cd ~

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

# creating local Python install dir if does not exist
if ! [ -d $LOCALPYDIR ];
  then mkdir -p $LOCALPYDIR;
fi

# doing all downloads and builds in this dir
# later will be easy to remove
mkdir $BUILDDIR
cd $BUILDDIR

# downloading pip install script
curl -L -O https://bootstrap.pypa.io/get-pip.py

# installing pip
python2.7 ./get-pip.py --user

# adding the local Python dir to the path
# otherwise Python won`t find it and gives
# `no module named ...` errors
export PYTHONPATH="$LOCALPYPATH:$PYTHONPATH"

# updating setuptools
pip install --upgrade --user setuptools

# installing pkg-config
curl -L $PKGCFGURL -o pkgconfig.tar.gz
tar -xzf pkgconfig.tar.gz
cd "pkg-config-$PKGCFGVER"
./configure --prefix=$LOCAL --with-internal-glib
make install
cd ..
# pkg-config installed

# installing pixman
curl -L $PXMANURL -o pixman.tar.gz
tar -xzf pixman.tar.gz
cd "pixman-$PXMANVER"
./configure --prefix=$LOCAL --disable-dependency-tracking
make install
cd ..
# pixman installed

# installing libpng
curl -L $LIBPNGURL -o libpng16.tar.xz
tar -xf libpng16.tar.xz
cd libpng-$LIBPNGVER
./configure --prefix=$LOCAL --disable-dependency-tracking
make install
cd ..
# libpng installed

# installing cairo
curl -L $CAIROURL -o cairo.tar.xz
tar -xf cairo.tar.xz
cd cairo-$CAIROVER
./configure --prefix=$LOCAL --disable-dependency-tracking
make install
cd ..
# cairo installed

# installing pycairo
export PKG_CONFIG_PATH="$LOCAL/lib/pkgconfig"
export ARCHFLAGS="-arch x86_64"
curl -L $PYCAIROURL -o pycairo.tar.bz2
tar -xjf pycairo.tar.bz2
cd py2cairo-$PYCAIROVER
./waf configure --prefix=$LOCAL
./waf build
./waf install
cd ..
export PYTHONPATH="$LOCAL/lib/python2.7/site-packages:$PYTHONPATH"
# pycairo installed

# installing autoconf, automake and libtool
# this is needed to build graphviz from git
curl -L http://ftpmirror.gnu.org/autoconf/autoconf-latest.tar.gz -o autoconf.tar.gz
tar xzf autoconf.tar.gz
cd autoconf-*
./configure --prefix=$LOCAL
make install
cd ..

curl -L http://ftpmirror.gnu.org/automake/automake-$AMVER.tar.xz -o automake.tar.xz
tar -xf automake.tar.xz
cd automake-$AMVER
./configure --prefix=$LOCAL
make install
cd ..

curl -L http://ftpmirror.gnu.org/libtool/libtool-$LTOOLVER.tar.xz -o libtool.tar.xz
tar -xf libtool.tar.xz
cd libtool-$LTOOLVER
./configure --prefix=$LOCAL
make install
cd ..
# automake, autoconf and libtool are installed

# installing graphviz -- optional
curl -L $GVIZURL --retry 5 -o graphviz.tar.gz
tar -xzf graphviz.tar.gz
cd "graphviz-master"
sh ./autogen.sh
./configure --prefix=$LOCAL
make install
cd ..
# graphviz installed

# installing igraph
pip install --user python-igraph

# optional dependencies
pip install --user pysftp
pip install --user fabric
pip install --user scipy
pip install --user pandas
pip install --user https://bitbucket.org/jurko/suds/get/tip.zip
pip install --user bioservices
pip install --user pymysql
pip install --global-option=build_ext --global-option="-I$LOCAL/include/" --global-option="-L$LOCAL/lib" --user pygraphviz

# installing pypath
pip install --user $PYPATHURL

#removing build dir
cd ~
rm -Rf $BUILDDIR

# adding local paths and python paths permantently
cat << EOF >> .bash_profile
export PYTHONPATH="$LOCALPYPATH2:$LOCALPYPATH:\$PYTHONPATH"
export PATH="$LOCALPYBIN:$LOCALPYBIN2:$LOCALBIN:\$PATH"
EOF


exit 0
