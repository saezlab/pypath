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
import sys
import os
from setuptools import setup, find_packages
import glob
import imp

from setuptools.command.install import install as _install
import pip
import pip._internal

_version = imp.load_source('_version', os.path.join('src', 'pypath', '_version.py'))
__version__ = _version.__version__

class __install(_install):
    
    def run(self):
        
        pip._internal.main(['install', 'git+git://github.com/fabioticconi/fishers_exact_test@fix-setup-deps'])
        
        _install.run(self)

def which(exe):
    '''
    Checks if executable is available.
    Source:
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    '''
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    
    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext
    
    fpath, fname = os.path.split(exe)
    if fpath:
        if is_exe(exe):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, exe)
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return True
    return False

def query_yes_no(question, default="yes"):
    """
    From http://stackoverflow.com/questions/3041986/python-command-line-yes-no-input
    
    Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)
    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

metainfo = {
    'authors': {
        'Türei':('Dénes Türei','turei.denes@gmail.com'),
    },
    'version': __version__,
    'license': 'GPLv3',
    'download_url': ['http://pypath.omnipathdb.org'],
    'url': ['http://pypath.omnipathdb.org'],
    'description': 'Work with molecular networks in Python igraph',
    'platforms': ['Linux', 'Unix', 'MacOSX', 'Windows'],
    'keywords': ['protein', 'mRNA', 'DNA', 'signaling',
                 'SignaLink', 'Signor', 'InnateDB', 'IntAct', 'Reactome',
                 'MPPI', 'NCI-PID', 'DIP', 'MatrixDB', 'PANTHER',
                 'PhosphoSite', 'PhosphoPoint', 'DEPOD', 'SPIKE', 'KEGG',
                 'Autophagy', 'ARN', 'NRF2', 'NRF2ome', 'Guide to Pharmacology', 
                 'regulation',
                 'phosphorylation', 'kinase', 'phosphatase',
                 'dephosphorylation', 'directed graph'],
    'classifiers': [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Natural Language :: English',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Information Analysis',
    'Topic :: Scientific/Engineering :: Mathematics']
}

with open('README.rst') as f:
    readme = f.read()
with open('HISTORY.rst') as f:
    history = f.read()

deps = [
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
    'tqdm'
]

#deplinks = [
    #'https://github.com/fabioticconi/fishers_exact_test/tarball/fix-setup-deps#egg=fisher-0.1.5'
#]

setup(
    cmdclass = {'install': __install},
    name = 'pypath',
    version = metainfo['version'],
    maintainer = metainfo['authors']['Türei'][0],
    maintainer_email = metainfo['authors']['Türei'][1],
    author = metainfo['authors']['Türei'][0],
    author_email = metainfo['authors']['Türei'][1],
    long_description = readme + '\n\n' + history,
    keywords = metainfo['keywords'],
    description = metainfo['description'],
    license = metainfo['license'],
    platforms = metainfo['platforms'],
    url = metainfo['url'],
    download_url = metainfo['download_url'],
    classifiers = metainfo['classifiers'],
    # package installation
    package_dir = {'':'src'},
    packages = list(set(find_packages() + ['pypath', 'pypath.data'])),
    include_package_data = True,
    install_requires = deps
    # dependency_links = deplinks
)
