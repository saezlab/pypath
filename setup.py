#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

__revision__ = "$Id$"
import os
from setuptools import setup, find_packages
import imp

from setuptools.command.install import install as _install


_version = imp.load_source(
    '_version',
    os.path.join('src', 'pypath', '_version.py')
)
__version__ = _version.__version__


with open('README.rst') as f:
    readme = f.read()
with open('HISTORY.rst') as f:
    history = f.read()


setup(
    cmdclass = {'install': __install},
    name = 'pypath',
    version = __version__,
    maintainer = 'Dénes Türei, Nicolas Palacio',
    maintainer_email = 'turei.denes@gmail.com',
    author = 'Dénes Türei, Nicolas Palacio',
    author_email = 'turei.denes@gmail.com',
    long_description = readme + '\n\n' + history,
    keywords = [
        'protein', 'mRNA', 'DNA', 'signaling',
        'SignaLink', 'Signor', 'InnateDB', 'IntAct', 'Reactome',
        'MPPI', 'NCI-PID', 'DIP', 'MatrixDB', 'PANTHER',
        'PhosphoSite', 'PhosphoPoint', 'DEPOD', 'SPIKE', 'KEGG',
        'Autophagy', 'ARN', 'NRF2', 'NRF2ome', 'Guide to Pharmacology',
        'regulation', 'phosphorylation', 'kinase', 'phosphatase',
        'dephosphorylation', 'directed graph',
    ],
    description = 'Molecular networks in Python',
    license = 'GPLv3',
    platforms = ['Linux', 'Unix', 'MacOSX', 'Windows'],
    url = ['http://pypath.omnipathdb.org/'],
    download_url = ['http://pypath.omnipathdb.org/releases/'],
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    # package installation
    package_dir = {'':'src'},
    packages = list(set(find_packages() + ['pypath', 'pypath.data'])),
    include_package_data = True,
    # dependency_links = deplinks
    install_requires = [
        'python-igraph',
        'beautifulsoup4',
        'configparser',
        'pyopenssl',
        'ndg-httpsclient',
        'chembl_webresource_client',
        'pyasn1',
        'numpy',
        'scipy',
        'matplotlib',
        'statsmodels',
        'pycurl',
        'lxml',
        'xlrd',
        'httplib2',
        'future',
        'tqdm',
    ],
    extras_require = {
        'tests': [
            'pytest',
        ],
        'docs': [
            'sphinx',
        ]
    }
)
