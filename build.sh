#!/bin/bash

rm src/*.pyc
rm src/pypath/*.pyc
rm -r src/pypath.egg-info

sed -i 's/\([0-9]*\.[0-9]*\.\)\([0-9]*\)/echo \1$\(\(\2+1\)\)/ge' src/pypath/__version__

python2 setup.py sdist
