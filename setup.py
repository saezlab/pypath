#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

__revision__ = "$Id$"
import os
from setuptools import setup, find_packages
import importlib as imp

from setuptools.command.install import install as _install

version_mod_path = os.path.join('src', 'pypath', '_version.py')
spec = imp.util.spec_from_file_location('_version', version_mod_path)
_version = imp.util.module_from_spec(spec)
spec.loader.exec_module(_version)

__version__ = _version.__version__


with open('README.rst', 'r') as f:
    readme = f.read()
with open('HISTORY.rst', 'r') as f:
    history = f.read()

with open('desc.rst', 'w') as f:
    f.write(readme + '\n' + history)


ENTRY_POINTS = {
    'console_scripts': [
        'bio2bel_omnipath = pypath.omnipath.bel:main',
    ],
    'bio2bel': [
        'omnipath = pypath.omnipath.bel',
    ],
}


def read_requirements():

    with open('requirements.txt', 'r') as fp:

        requirements = [
            name.strip()
            for name in fp
            if name and not name.startswith('-')
        ]

    return requirements


setup(
    name = 'pypath-omnipath',
    version = __version__,
    maintainer = 'Dénes Türei, Sebastian Lobentanzer, Ahmet Rifaioglu',
    maintainer_email = 'turei.denes@gmail.com',
    author = (
        'Dénes Türei, Nicolàs Palacio, Olga Ivanova, '
        'Sebastian Lobentanzer, Ahmet Rifaioglu'
    ),
    author_email = 'turei.denes@gmail.com',
    long_description = readme + '\n' + history,
    long_description_content_type = 'text/x-rst; charset=UTF-8',
    keywords = sorted({
        'protein', 'mRNA', 'miRNA', 'DNA', 'signaling',
        'SignaLink', 'SIGNOR', 'InnateDB', 'IntAct', 'Reactome',
        'MPPI', 'NCI-PID', 'DIP', 'MatrixDB', 'PANTHER',
        'PhosphoSite', 'PhosphoPoint', 'DEPOD', 'SPIKE', 'KEGG',
        'Autophagy', 'ARN', 'NRF2ome', 'Guide to Pharmacology',
        'UniProt', 'BioPAX', 'Ensembl', 'Surfaceome',
        'Exocarta', 'Vesiclepedia', 'Matrisome', 'Human Protein Atlas',
        'Compleat', 'CORUM', 'ComplexPortal', 'BioGRID', 'STRING',
        'ICELLNET', 'Cell Surface Protein Atlas', 'COSMIC',
        'Cancer Gene Census', 'IntOGen', 'TopDB', 'iTALK',
        'Human Plasma Membrane Receptome', 'EMBRACE', 'ELM', 'phospho.ELM',
        'CancerSEA', 'ComPPI', 'CellPhoneDB', 'DGIdb', 'DisGeNet',
        'PAZAR', 'ORegAnno', 'TRED', 'DoRothEA', 'TRRD', 'CPAD',
        'regulation', 'phosphorylation', 'kinase', 'phosphatase',
        'dephosphorylation', 'directed graph', 'annotations', 'cancer',
        'complexes', 'intercellular communication', 'HGNC', 'GPCRdb',
        'MSigDB', 'GSEA', 'Phobius', 'Phosphatome', 'NetPath',
        'gene', 'gene symbol', 'mouse', 'rat', 'HomoloGene',
        'integrin', 'adhesion', 'receptor', 'ligand', 'transporter',
        'ion channel', 'disease', 'activity flow', 'transcription', 'PPI',
        'subcellular localization', 'pathway', 'signaling pathway',
        'cell-cell communication', 'scRNA-Seq', 'Boolean modeling',
        'causal reasoning', 'systems biomedicine', 'molecular biology',
        'prior knowledge', 'literature knowledge', 'systems biology',
        'SBGN', 'Cellinker', 'scConnect', 'CellChatDB', 'CellTalkDB',
        'talklr',
    }),
    description = 'Molecular signaling prior knowledge in Python',
    license = 'GPLv3',
    platforms = ['Linux', 'Unix', 'MacOSX', 'Windows'],
    url = 'https://omnipathdb.org/',
    download_url = 'https://pypath.omnipathdb.org/releases/',
    project_url = ('Git repo', 'https://github.com/saezlab/pypath'),
    classifiers = [
        'Development Status :: 4 - Beta',
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
    install_requires = read_requirements(),
    extras_require = {
        'tests': [
            'pytest',
        ],
        'docs': [
            'sphinx',
        ],
        'bel': [
            'pybel',
            'bio2bel',
            'click',
        ],
        'graph': [
            'python-igraph',
        ],
    },
    entry_points = ENTRY_POINTS,
)
