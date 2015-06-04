#!/bin/bash

rm src/*.pyc
rm src/bioigraph/*.pyc
rm -r src/bioigraph.egg-info

sed -i 's/\([0-9]*\.[0-9]*\.\)\([0-9]*\)/echo \1$\(\(\2+1\)\)/ge' src/bioigraph/__version__

python2 setup.py sdist
