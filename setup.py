# -*- coding: utf-8 -*-
#from ez_setup import use_setuptools
#use_setuptools()
__revision__ = "$Id$"
import sys
import os
from setuptools import setup, find_packages
import glob
import imp

common = imp.load_source('common', os.path.join('src', 'pypath', 'common.py'))
__version__ = common.__version__

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
    'Türei':('Dénes Türei','denes@ebi.ac.uk'),
    },
    'version': __version__,
    'license': 'LGPL',
    'download_url': ['http://157.181.231.40/~denes/pypath'],
    'url': ['http://157.181.231.40/~denes/pypath'],
    'description': 'Work with molecular networks in Python igraph',
    'platforms': ['Linux', 'Unix', 'MacOSX', 'Windows'],
    'keywords': ['graph', 'network', 'protein', 'mRNA', 'DNA', 'signaling',
                 'SignaLink', 'Signor', 'InnateDB', 'IntAct', 'Reactome',
                 'MPPI', 'NCI-PID', 'DIP', 'MatrixDB', 'PANTHER',
                 'PhosphoSite', 'PhosphoPoint', 'DEPOD', 'SPIKE', 'KEGG',
                 'Autophagy', 'ARN', 'NRF2', 'NRF2ome', 'Guide to Pharmacology', 
                 'regulation',
                 'phosphorylation', 'kinase', 'phosphatase',
                 'dephosphorylation', 'directed graph'],
    'classifiers': [
    'Development Status :: 2 - Pre-Alpha',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: GNU Library or Lesser General Public License (LGPL)',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 2.7',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Information Analysis',
    'Topic :: Scientific/Engineering :: Mathematics']
}

#package_data = [
    #'pypath/data/alzpw-ppi.csv', 
    #'pypath/data/arn_014.csv', 
    #'pypath/data/arn.csv', 
    #'pypath/data/arn_lr01.csv', 
    #'pypath/data/ca1.csv', 
    #'pypath/data/cancer_gene_census.csv', 
    #'pypath/data/cell-map-edge-attributes.txt', 
    #'pypath/data/dd_refs.csv', 
    #'pypath/data/depod-refs.csv', 
    #'pypath/data/dip_human_core_processed.csv', 
    #'pypath/data/entrez_uniprot.csv', 
    #'pypath/data/gdsc.sif', 
    #'pypath/data/gold_standard.csv', 
    #'pypath/data/gold_standard.xlsx', 
    #'pypath/data/innatedb.csv', 
    #'pypath/data/intact_filtered.csv', 
    #'pypath/data/krishnadev_atg_1.tab', 
    #'pypath/data/krishnadev_atg.tab', 
    #'pypath/data/krishnadev_vegeredmeny.csv', 
    #'pypath/data/kshirsagar_atg_1.tab', 
    #'pypath/data/kshirsagar_atg.tab', 
    #'pypath/data/kshirsagar_vegeredmeny.csv', 
    #'pypath/data/macrophage-strict.csv', 
    #'pypath/data/matrixdb_core.csv', 
    #'pypath/data/mppi_human_rep.csv', 
    #'pypath/data/nci-pid-strict.csv', 
    #'pypath/data/netpath_refs.csv', 
    #'pypath/data/nrf2ome.csv', 
    #'pypath/data/phosphopoint.csv', 
    #'pypath/data/phosphosite_human_hc.csv', 
    #'pypath/data/phosphosite_human_noref.csv', 
    #'pypath/data/salmonella_atg.tar.gz', 
    #'pypath/data/sec_ac.txt', 
    #'pypath/data/shlecker_atg_1.tab', 
    #'pypath/data/shlecker_vegeredmeny.csv', 
    #'pypath/data/signor_ppi.tsv', 
    #'pypath/data/slk01human.csv', 
    #'pypath/data/spike_hc.csv', 
    #'pypath/data/swissprot2.csv', 
    #'pypath/data/swissprot-gsymbol-name.csv', 
    #'pypath/data/trembl2.csv', 
    #'pypath/data/uniprot-all-human.tab',
    #'pypath/data/intogene_cancerdrivers.tsv'
#]

with open('README.rst') as f:
    readme = f.read()
with open('HISTORY.rst') as f:
    history = f.read()

# choosing module for mysql access:
deps = [
    'python-igraph',
    'beautifulsoup4',
    'configparser',
    'pyopenssl',
    'ndg-httpsclient',
    'chembl_webresource_client',
    'pyasn1',
    'fisher',
    'statsmodels',
    'pycurl',
    'lxml',
    'xlrd',
    'httplib2',
    'future'
]

#mysql = 'pymysql'
#if which('mysql') and which('mysql_config'):
    #mysql_alt = query_yes_no('Looks like MySQL is installed on your system. \n'\
        #'Do you want to use MySQL-python instead of pymysql?')
    #if mysql_alt:
        #mysql = 'MySQL-python'

#deps.append(mysql)

setup(
    name = 'pypath',
    version = __version__,
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
)
