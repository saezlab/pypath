#!/bin/bash

rm src/*.pyc
rm src/bioigraph/*.pyc
rm -r src/bioigraph.egg-info

sed -i 's/\(_MICRO = \)\([0-9]*\)/echo \1$\(\(\2+1\)\)/ge' setup.py

python2 setup.py sdist
