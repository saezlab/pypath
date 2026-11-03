#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

import sys
from future.utils import iteritems

import codecs
import bs4
import textwrap

import pypath.omnipath.server._html as _html
import pypath.resources.data_formats as data_formats
import pypath.resources.urls as urls
import pypath.share.session as session_mod

__all__ = ['descriptions', 'gen_html', 'write_html']

if 'long' not in __builtins__:
    long = int

if 'unicode' not in __builtins__:
    unicode = str


_logger = session_mod.Logger(name = 'descriptions')
_log = _logger._log


descriptions = {
    'HuRI': {
        'year': 2016,
        'releases': [2012, 2014, 2016],
        'recommend':
        'very large, quality controlled, unbiased yeast-2-hybrid screening',
        'label': 'HuRI HI-III',
        'full_name': 'Human Reference Interactome',
        'urls': {
            'articles':
            ['http://www.cell.com/cell/abstract/S0092-8674(14)01422-6'],
            'webpages': [
                'http://interactome.dfci.harvard.edu/H_sapiens/',
                'http://www.interactome-atlas.org/',
            ],
        },
        'pubmeds': [25416956],
        'emails':
        [('Michael_Calderwood@dfci.harvard.edu', 'Michael Calderwood')],
        'type': 'high-throughput',
        'subtype': 'yeast 2 hybrid',
        'omnipath': False,
        'license': {
            'name':
            'No license. "This dataset is freely available to the research community through the search engine or via download."',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'HuRI Lit-BM': {
        'year': 2017,
        'releases': [2013, 2017],
        'label': 'HuRI Lit-BM-17',
        'full_name': 'Human Reference Interactome Literature Benchmark',
        'urls': {
            'articles':
            ['http://www.cell.com/cell/abstract/S0092-8674(14)01422-6'],
            'webpages': [
                'http://interactome.dfci.harvard.edu/H_sapiens/',
                'http://www.interactome-atlas.org/',
            ],
        },
        'authors': ['CCSB'],
        'pubmeds': [25416956],
        'descriptions': [
            u'''
            High-quality non-systematic Literature dataset. In 2013, we extracted interaction data from BIND, BioGRID, DIP, HPRD, MINT, IntAct, and PDB to generate a high-quality binary literature dataset comprising ~11,000 protein-protein interactions that are binary and supported by at least two traceable pieces of evidence (publications and/or methods) (Rolland et al Cell 2014). Although this dataset does not result from a systematic investigation of the interactome search space and should thus be used with caution for any network topology analyses, it represents valuable interactions for targeted studies and is freely available to the research community through the search engine or via download.
            '''
        ],
        'emails':
        [('Michael_Calderwood@dfci.harvard.edu', 'Michael Calderwood')],
        'type': 'high-throughput',
        'subtype': 'yeast 2 hybrid',
        'omnipath': False,
        'pypath': {
            'get': ['pypath.dataio.get_lit_bm_13()'],
            'data': ['pypath.urls.urls[\'hid\'][\'lit-bm-13\']'],
            'input': ['pypath.data_formats.interaction_misc[\'lit13\']']
        },
        'license': {
            'name':
            'No license. "This dataset is freely available to the research community through the search engine or via download."',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'ELM': {
        'year': 2014,
        'releases': [2003, 2008, 2009, 2012, 2013, 2014, 2016],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['ELM Consortium'],
        'label': 'ELM',
        'color': '',
        'urls': {
            'webpages': ['http://elm.eu.org/'],
            'articles': [
                'http://nar.oxfordjournals.org/content/40/D1/D242.long',
                'http://nar.oxfordjournals.org/content/42/D1/D259.long',
                'http://nar.oxfordjournals.org/content/44/D1/D294.long'
            ],
            'omictools':
            ['http://omictools.com/eukaryotic-linear-motif-resource-tool']
        },
        'annot': ['domain', 'residue'],
        'recommend':
        'structural details: domain-motif relationships; very high confidence',
        'pubmeds': [22110040, 24214962, 26615199],
        'emails': [('feedback@elm.eu.org', 'ELM Team'), ('gibson@embl.de',
                                                         'Toby Gibson')],
        'type': 'literature curated',
        'subtype': 'post-translational modifications',
        'data_integration': 'dynamic',
        'descritpions': [
            u'''
            Ideally, each motif class has multiple example instances of this motif annotated, whereby an instance is described as a match to the regular expression pattern of the ELM motif class in a protein sequence. For each instance entry, ideally, multiple sources of experimental evidence are recorded (identifying participant, detecting motif presence and detecting interaction), and, following annotation best practices, a reliability score is given by the annotator.
            '''
        ],
        'omnipath': True,
        'license': {
            'name': 'ELM Software License Agreement, non-free',
            'url': 'http://elm.eu.org/media/Elm_academic_license.pdf',
            'commercial_use': False
        },
        'pypath': {
            'data': [
                'pypath.urls.urls[\'ielm_domains\'][\'url\']',
                'pypath.urls.urls[\'elm_class\'][\'url\']',
                'pypath.urls.urls[\'elm_inst\'][\'url\']',
                'urls.urls[\'elm_int\'][\'url\']'
            ],
            'format': [
                'pypath.data_formats.ptm[\'elm\']',
                'pypath.data_formats.omnipath[\'elm\']'
            ],
            'input': [
                'pypath.dataio.get_elm_domains()',
                'pypath.dataio.get_elm_classes()',
                'pypath.dataio.get_elm_instances()'
            ],
            'intr': ['pypath.dataio.get_elm_interactions()'],
            'dmi': ['pypath.pypath.PyPath.load_elm()']
        }
    },
    'LMPID': {
        'year': 2015,
        'releases': [2015],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['Bose Institute'],
        'label': 'LMPID',
        'color': '',
        'urls': {
            'webpages':
            ['http://bicresources.jcbose.ac.in/ssaha4/lmpid/index.php'],
            'articles':
            ['http://database.oxfordjournals.org/content/2015/bav014.long'],
            'omictools': [
                'http://omictools.com/linear-motif-mediated-protein-interaction-database-tool'
            ]
        },
        'pubmeds': [25776024],
        'emails': [('ssaha4@jcbose.ac.in', 'Sudipto Saha')],
        'type': 'literature curated',
        'subtype': 'post-translational modifications',
        'recommend':
        'structural details: domain-motif relationships; similar to ELM, but more recent and larger',
        'annot': ['domain', 'mechanism'],
        'descriptions': [
            u'''
            LMPID (Linear Motif mediated Protein Interaction Database) is a manually curated database which provides comprehensive experimentally validated information about the LMs mediating PPIs from all organisms on a single platform. About 2200 entries have been compiled by detailed manual curation of PubMed abstracts, of which about 1000 LM entries were being annotated for the first time, as compared with the Eukaryotic LM resource.
            '''
        ],
        'omnipath': True,
        'pypath': {
            'input': ['pypath.dataio.load_lmpid()'],
            'intr': ['pypath.dataio.lmpid_interactions()'],
            'data': ['pypath.data/LMPID_DATA_pubmed_ref.xml'],
            'format': [
                'pypath.data_formats.ptm[\'lmpid\']',
                'pypath.data_formats.omnipath[\'lmpid\']'
            ],
            'dmi': [
                'pypath.dataio.lmpid_dmi()',
                'pypath.pypath.PyPath().process_dmi(source = \'LMPID\')'
            ]
        },
        'license': {
            'name':
            'No license. If you are using this database please cite Sarkar 2015.',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'PDZBase': {
        'year': 2004,
        'releases': [2004],
        'authors': ['Weinstein Group'],
        'urls': {
            'webpages': ['http://abc.med.cornell.edu/pdzbase'],
            'articles': [
                'http://bioinformatics.oxfordjournals.org/content/21/6/827.long'
            ],
            'omictools': ['http://omictools.com/pdzbase-tool']
        },
        'pubmeds': [15513994],
        'taxons': ['human'],
        'color': None,
        'label': 'PDZBase',
        'annot': ['domain'],
        'recommend':
        'a handful specific interactions for proteins with PDZ domain',
        'descriptions': [
            u'''
            PDZBase is a database that aims to contain all known PDZ-domain-mediated protein-protein interactions. Currently, PDZBase contains approximately 300 such interactions, which have been manually extracted from &gt;200 articles.
            PDZBase currently contains ∼300 interactions, all of which have been manually extracted from the literature, and have been independently verified by two curators. The extracted information comes from in vivo (co-immunoprecipitation) or in vitro experiments (GST-fusion or related pull-down experiments). Interactions identified solely from high throughput methods (e.g. yeast two-hybrid or mass spectrometry) were not included in PDZBase. Other prerequisites for inclusion in the database are: that knowledge of the binding sites on both interacting proteins must be available (for instance through a truncation or mutagenesis experiment); that interactions must be mediated directly by the PDZ-domain, and not by any other possible domain within the protein.
            '''
        ],
        'emails': [('haw2002@med.cornell.edu', 'Harel Weinstein'),
                   ('pdzbase@med.cornell.edu', 'PDZBase Team')],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'pypath': {
            'intr': ['pypath.dataio.get_pdzbase()'],
            'data': [
                'pypath.urls.urls[\'pdzbase\']',
                'pypath.urls.urls[\'pdz_details\']'
            ],
            'format': [
                'pypath.data_formats.pathway[\'pdz\']',
                'pypath.data_formats.omnipath[\'pdz\']'
            ]
        },
        'license': {
            'name': 'No license.',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        },
        'pathguide': 160
    },
    'Guide2Pharma': {
        'year': 2015,
        'releases': [2007, 2008, 2009, 2011, 2013, 2014, 2015, 2016],
        'size': None,
        'authors': None,
        'label': 'Guide to Pharmacology',
        'full_name': 'Guide to Pharmacology',
        'color': None,
        'pubmeds': [24234439],
        'urls': {
            'webpages': ['http://www.guidetopharmacology.org/'],
            'articles': [
                'http://nar.oxfordjournals.org/content/42/D1/D1098.long',
                'http://onlinelibrary.wiley.com/doi/10.1111/j.1476-5381.2011.01649_1.x/full'
            ],
            'omictools': [
                'http://omictools.com/international-union-of-basic-and-clinical-pharmacology-british-pharmacological-society-guide-to-pharmacology-tool'
            ]
        },
        'recommend':
        'one of the strongest in ligand-receptor interactions; still does not contain everything, therefore worth to combine with larger activity flow resources like Signor',
        'descriptions': [
            u'''
            Presently, the resource describes the interactions between target proteins and 6064 distinct ligand entities (Table 1). Ligands are listed against targets by their action (e.g. activator, inhibitor), and also classified according to substance types and their status as approved drugs. Classes include metabolites (a general category for all biogenic, non-peptide, organic molecules including lipids, hormones and neurotransmitters), synthetic organic chemicals (e.g. small molecule drugs), natural products, mammalian endogenous peptides, synthetic and other peptides including toxins from non-mammalian organisms, antibodies, inorganic substances and other, not readily classifiable compounds.
            The new database was constructed by integrating data from IUPHAR-DB and the published GRAC compendium. An overview of the curation process is depicted as an organizational flow chart in Figure 2. New information was added to the existing relational database behind IUPHAR-DB and new webpages were created to display the integrated information. For each new target, information on human, mouse and rat genes and proteins, including gene symbol, full name, location, gene ID, UniProt and Ensembl IDs was manually curated from HGNC, the Mouse Genome Database (MGD) at Mouse Genome Informatics (MGI), the Rat Genome Database (RGD), UniProt and Ensembl, respectively. In addition, ‘Other names’, target-specific fields such as ‘Principal transduction’, text from the ‘Overview’ and ‘Comments’ sections and reference citations (downloaded from PubMed; http://www.ncbi.nlm.nih.gov/pubmed) were captured from GRAC and uploaded into the database against a unique Object ID.
            '''
        ],
        'emails':
        [('enquiries@guidetopharmacology.org', 'Guide to Pharmacology Team'),
         ('cdsouthan@hotmail.com', 'Cristopher Southan')],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'license': {
            'name': 'CC-Attribution-ShareAlike-3.0',
            'url': 'http://creativecommons.org/licenses/by-sa/3.0/',
            'commercial_use': True
        },
        'pypath': {
            'data': ['pypath.data_formats'],
            'format': [
                'pypath.data_formats.pathway[\'guide2pharma\']',
                'pypath.data_formats.omnipath[\'guide2pharma\']'
            ]
        },
        'pathguide': 345
    },
    'phosphoELM': {
        'year': 2010,
        'releases': [2004, 2007, 2010],
        'urls': {
            'webpages': ['http://phospho.elm.eu.org/'],
            'articles': [
                'http://www.biomedcentral.com/1471-2105/5/79',
                'http://nar.oxfordjournals.org/content/36/suppl_1/D240.full',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D261'
            ],
            'omictools': ['http://omictools.com/phospho-elm-tool']
        },
        'pubmeds': [15212693, 17962309, 21062810],
        'annot': ['mechanism', 'residue'],
        'recommend':
        'one of the largest kinase-substrate databases; substantial number of specific proteins and interactions, with more receptors than PhosphoSite',
        'descriptions': [
            u'''
            Phospho.ELM http://phospho.elm.eu.org is a new resource containing experimentally verified phosphorylation sites manually curated from the literature and is developed as part of the ELM (Eukaryotic Linear Motif) resource. Phospho.ELM constitutes the largest searchable collection of phosphorylation sites available to the research community. The Phospho.ELM entries store information about substrate proteins with the exact positions of residues known to be phosphorylated by cellular kinases. Additional annotation includes literature references, subcellular compartment, tissue distribution, and information about the signaling pathways involved as well as links to the molecular interaction database MINT. Phospho.ELM version 2.0 contains 1,703 phosphorylation site instances for 556 phosphorylated proteins. (Diella 2004)
            ''', u'''
            Phospho.ELM is a manually curated database of eukaryotic phosphorylation sites. The resource includes data collected from published literature as well as high-throughput data sets. The current release of the Phospho.ELM data set (version 7.0, July 2007) contains 4078 phospho-protein sequences covering 12,025 phospho-serine, 2,362 phospho-threonine and 2,083 phospho-tyrosine sites with a total of 16,470 sites.
            For each phospho-site we report if the phosphorylation evidence has been identified by small-scale analysis (low throughput; LTP) that typically focus on one or a few proteins at a time or by large-scale experiments (high throughput; HTP), which mainly apply MS techniques. It is noteworthy that in our data set there is a small overlap between instances identified by LTP and HTP experiments. (Diella 2007)
            ''', u'''
            The current release of the Phospho.ELM data set (version 9.0) contains more than 42,500 non-redundant instances of phosphorylated residues in more than 11,000 different protein sequences (3370 tyrosine, 31 754 serine and 7449 threonine residues). For each phosphosite we report whether the phosphorylation evidence has been identified by small-scale analyses (low-throughput, LTP) and/or by large-scale experiments (high-throughput, HTP), which mainly apply MS techniques. The majority of the protein instances from Phospho. ELM are vertebrate (mostly Homo sapiens (62%) and Mus musculus (16%)) though 22% are from other species, mainly Drosophila melanogaster (13%) and Caenorhabditis elegans (7%). In total, more than 300 different kinases have been annotated and a document providing additional information about all kinases annotated in Phospho.ELM can be found at http://phospho.elm.eu.org/kinases.html. (Dinkel 2010)
            '''
        ],
        'emails': [('toby.gibson@embl.de', 'Toby Gibson')],
        'type': 'Literature curated',
        'data_integration': 'dynamic',
        'subtype': 'PTM',
        'label': 'phospho.ELM',
        'omnipath': True,
        'license': {
            'name': 'phospho.ELM Academic License, non-free',
            'url':
            'http://phospho.elm.eu.org/dumps/Phospho.Elm_AcademicLicense.pdf',
            'commercial_use': False
        },
        'pypath': {
            'format': [
                'pypath.data_formats.ptm[\'phelm\']',
                'pypath.data_formats.omnipath[\'phelm\']'
            ],
            'data': [
                'pypath.urls.urls[\'p_elm\'][\'psites\']',
                'urls.urls[\'p_elm_kin\'][\'url\']'
            ],
            'intr': ['pypath.dataio.phelm_interactions()'],
            'input': [
                'pypath.dataio.get_phosphoelm()',
                'pypath.dataio.get_phelm_kinases()',
                'pypath.dataio.phelm_psites()'
            ],
            'ptm': ['pypath.pypath.PyPath().load_phosphoelm()']
        },
        'pathguide': 138
    },
    'DOMINO': {
        'year': 2006,
        'releases': [2006],
        'authors': ['Cesareni Group'],
        'urls': {
            'webpages':
            ['http://mint.bio.uniroma2.it/domino/search/searchWelcome.do'],
            'articles':
            ['http://nar.oxfordjournals.org/content/35/suppl_1/D557.long'],
            'omictools': ['http://omictools.com/domino-tool']
        },
        'pubmeds': [17135199],
        'taxons': [
            'human', 'yeast', 'C. elegans', 'mouse', 'rat', 'HIV',
            'D. melanogaster', 'A. thaliana', 'X. laevis', 'B. taurus',
            'G. gallus', 'O. cuniculus', 'Plasmodium falciparum'
        ],
        'annot': ['experiment'],
        'recommend':
        'rich details and many specific information; discontinued, Signor from the same lab is larger and newer, and contains most of its data',
        'descriptions': [
            u'''
            DOMINO aims at annotating all the available information about domain-peptide and domain–domain interactions. The core of DOMINO, of July 24, 2006 consists of more than 3900 interactions extracted from peer-reviewed articles and annotated by expert biologists. A total of 717 manuscripts have been processed, thus covering a large fraction of the published information about domain–peptide interactions. The curation effort has focused on the following domains: SH3, SH2, 14-3-3, PDZ, PTB, WW, EVH, VHS, FHA, EH, FF, BRCT, Bromo, Chromo and GYF. However, interactions mediated by as many as 150 different domain families are stored in DOMINO.
            ''', u'''
            DOMINO is an open-access database comprising more than 3900 annotated experiments describing interactions mediated by protein-interaction domains. The curation effort aims at covering the interactions mediated by the following domains (SH3, SH2, 14-3-3, PDZ, PTB, WW, EVH, VHS, FHA, EH, FF, BRCT, Bromo, Chromo, GYF). However, interactions mediated by as many as 150 different domain families are stored in DOMINO.
            ''', u'''
            The curation process follows the PSI-MI 2.5 standard but with special emphasis on the mapping of the interaction to specific protein domains of both participating proteins. This is achieved by paying special attention to the shortest protein fragment that was experimentally verified as sufficient for the interaction. Whenever the authors report only the name of the domain mediating the interaction (i.e. SH3, SH2 ...), without stating the coordinates of the experimental binding range, the curator may choose to enter the coordinates of the Pfam domain match in the protein sequence. Finally whenever the information is available, any mutation or posttranslational modification affecting the interaction affinity is noted in the database.
            '''
        ],
        'emails': [('giovanni.cesareni@uniroma2.it', 'Gianni Cesareni')],
        'type': 'Literature curated',
        'subtype': 'PTM',
        'omnipath': True,
        'license': {
            'name': 'CC-Attribution-2.5',
            'url': 'http://creativecommons.org/licenses/by/2.5',
            'commercial_use': True
        },
        'pypath': {
            'data': ['pypath.urls.urls[\'domino\'][\'url\']'],
            'input':
            ['pypath.dataio.get_domino()', 'pypath.dataio.get_domino_ddi()'],
            'format': [
                'pypath.data_formats.ptm[\'domino\']',
                'pypath.data_formats.omnipath[\'domino\']'
            ],
            'intr': ['pypath.dataio.domino_interactions()'],
            'ptm': ['pypath.dataio.get_domino_ptms()'],
            'dmi': [
                'pypath.dataio.get_domino_dmi()',
                'pypath.pypath.PyPath().load_domino_dmi()'
            ],
            'ddi': ['pypath.dataio.get_domino_ddi()']
        },
        'pathguide': 239
    },
    'dbPTM': {
        'year': 2015,
        'releases': [2005, 2009, 2012, 2015],
        'authors': ['ISBLab'],
        'urls': {
            'webpages': ['http://dbptm.mbc.nctu.edu.tw/'],
            'articles': [
                'http://nar.oxfordjournals.org/content/41/D1/D295.long',
                'http://www.biomedcentral.com/1756-0500/2/111',
                'http://nar.oxfordjournals.org/content/34/suppl_1/D622.long',
                'http://nar.oxfordjournals.org/content/44/D1/D435.long'
            ],
            'omictools': ['http://omictools.com/dbptm-tool']
        },
        'pubmeds': [16381945, 19549291, 23193290, 26578568],
        'taxons': ['human', 'Metazoa', 'Bacteria', 'plants', 'yeast'],
        'annot': ['mechanism', 'residue'],
        'recommend':
        'integrates many small efforts; beside phosphorylations provides all types of PTMs and enzyme-substrate relationships',
        'descriptions': [
            u'''
            Due to the inaccessibility of database contents in several online PTM resources, a total 11 biological databases related to PTMs are integrated in dbPTM, including UniProtKB/SwissProt, version 9.0 of Phospho.ELM, PhosphoSitePlus, PHOSIDA, version 6.0 of O-GLYCBASE, dbOGAP, dbSNO, version 1.0 of UbiProt, PupDB, version 1.1 of SysPTM and release 9.0 of HPRD.
            With the high throughput of MS-based methods in post-translational proteomics, this update also includes manually curated MS/MS-identified peptides associated with PTMs from research articles through a literature survey. First, a table list of PTM-related keywords is constructed by referring to the UniProtKB/SwissProt PTM list (http://www.uniprot.org/docs/ptmlist.txt) and the annotations of RESID (28). Then, all fields in the PubMed database are searched based on the keywords of the constructed table list. This is then followed by downloading the full text of the research articles. For the various experiments of proteomic identification, a text-mining system is developed to survey full-text literature that potentially describes the site-specific identification of modified sites. Approximately 800 original and review articles associated with MS/MS proteomics and protein modifications are retrieved from PubMed (July 2012). Next, the full-length articles are manually reviewed for precisely extracting the MS/MS peptides along with the modified sites. Furthermore, in order to determine the locations of PTMs on a full-length protein sequence, the experimentally verified MS/MS peptides are then mapped to UniProtKB protein entries based on its database identifier (ID) and sequence identity. In the process of data mapping, MS/MS peptides that cannot align exactly to a protein sequence are discarded. Finally, each mapped PTM site is attributed with a corresponding literature (PubMed ID).
            ''', u'''
            dbPTM was developed as a comprehensive database of experimentally verified PTMs from several databases with annotations of potential PTMs for all UniProtKB protein entries. For this tenth anniversary of dbPTM, the updated resource includes not only a comprehensive dataset of experimentally verified PTMs, supported by the literature, but also an integrative interface for accessing all available databases and tools that are associated with PTM analysis. As well as collecting experimental PTM data from 14 public databases, this update manually curates over 12,000 modified peptides, including the emerging S-nitrosylation, S-glutathionylation and succinylation, from approximately 500 research articles, which were retrieved by text mining.
            '''
        ],
        'emails': [('francis@saturn.yzu.edu.tw', 'Hsien-Da Huang'),
                   ('bryan@mail.nctu.edu.tw', 'Hsien-Da Huang')],
        'type': 'Literature curated',
        'subtype': 'PTM',
        'omnipath': True,
        'data_integration': 'dynamic',
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        },
        'pypath': {
            'format': [
                'pypath.data_formats.ptm[\'dbptm\']',
                'pypath.data_formats.omnipath[\'dbptm\']'
            ],
            'data': [
                'pypath.urls.urls[\'dbptm_benchmark\'][\'urls\']',
                'pypath.urls.urls[\'dbptm\'][\'urls\']'
            ],
            'input': ['pypath.dataio.get_dbptm()'],
            'intr': ['pypath.dataio.dbptm_interactions()'],
            'ptm': ['pypath.pypath.PyPath().load_dbptm()']
        }
    },
    'SIGNOR': {
        'year': 2015,
        'releases': [2015,2019],
        'urls': {
            'webpages': ['http://signor.uniroma2.it/'],
            'articles': [
                'http://nar.oxfordjournals.org/content/44/D1/D548',
                'http://f1000research.com/posters/1098359'
            ],
            'omictools':
            ['http://omictools.com/signaling-network-open-resource-tool']
        },
        'full_name': 'Signaling Network Open Resource',
        'pubmeds': [26467481],
        'annot': ['mechanism', 'pathway'],
        'recommend':
        'provides effect sign for an unprecedented number of interactions; large and recent curation effort; many specific entities; PTMs with enzymes',
        'descriptions': [
            u'''
            SIGNOR, the SIGnaling Network Open Resource, organizes and stores in a structured format signaling information published in the scientific literature. The captured information is stored as binary causative relationships between biological entities and can be represented graphically as activity flow. The entire network can be freely downloaded and used to support logic modeling or to interpret high content datasets. The core of this project is a collection of more than 11000 manually-annotated causal relationships between proteins that participate in signal transduction. Each relationship is linked to the literature reporting the experimental evidence. In addition each node is annotated with the chemical inhibitors that modulate its activity. The signaling information is mapped to the human proteome even if the experimental evidence is based on experiments on mammalian model organisms.
            '''
        ],
        'authors': ['Cesareni Group'],
        'label': 'SIGNOR',
        'color': '',
        'data_import': ['SignaLink3', 'PhosphoSite'],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'emails': [('perfetto@live.it', 'Livia Perfetto')],
        'license': {
            'name': 'CC-Attribution-ShareAlike 4.0',
            'url': 'https://creativecommons.org/licenses/by-sa/4.0/',
            'commercial_use': True
        },
        'pypath': {
            'data': ['pypath.urls.urls[\'signor\'][\'all_url\']'],
            'format': ['pypath.data_formats.pathway[\'signor\']'],
            'intr': ['pypath.dataio.signor_interactions()'],
            'ptm': [
                'pypath.dataio.load_signor_ptms()',
                'pypath.pypath.PyPath().load_signor_ptms()'
            ]
        }
    },
    'HuPho': {
        'year': 2015,
        'releases': [2012, 2015],
        'urls': {
            'webpages': ['http://hupho.uniroma2.it/'],
            'articles': [
                'http://onlinelibrary.wiley.com/doi/10.1111/j.1742-4658.2012.08712.x/full'
            ],
            'omictools':
            ['http://omictools.com/human-phosphatase-portal-tool']
        },
        'pubmeds': [22804825],
        'descriptions': [
            u'''
            In order to offer a proteome-wide perspective of the phosphatase interactome, we have embarked on an extensive text-mining-assisted literature curation effort to extend phosphatase interaction information that was not yet covered by protein–protein interaction (PPI) databases. Interaction evidence captured by expert curators was annotated in the protein interaction database MINT according to the rapid curation standard. This data set was next integrated with protein interaction information from three additional major PPI databases, IntAct, BioGRID and DIP. These databases are part of the PSIMEx consortium and adopt a common data model and common controlled vocabularies, thus facilitating data integration. Duplicated entries were merged and redundant interactions have been removed.
            As a result, from the HuPho website it is possible to explore experimental evidence from 718 scientific articles reporting 4600 experiments supporting protein interactions where at least one of the partners is a phosphatase. Since some interactions are supported by more than one piece of evidence, the actual number of non-redundant interactions is smaller, 2500 at the time of writing this paper. Moreover, 199 phosphatases have at least one reported ligand, while 53 have none. Interaction evidence is fairly evenly distributed in the four PSIMEx resources suggesting a substantial lack of overlap among the data curated by each database.
            '''
        ],
        'notes': [
            u'''
            The database is dynamically updated, so is up to date at any given time. That's why it is marked as up to date in 2015, despite it has no new release after 2012.
            '''
        ],
        'authors': ['Cesareni Group'],
        'label': 'HuPho',
        'full_name': 'Human Phosphatase Portal',
        'color': '',
        'type': 'high throughput and literature curated',
        'subtype': 'post-translational modification',
        'omnipath': False,
        'emails': [('perfetto@live.it', 'Livia Perfetto')],
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'SignaLink3': {
        'year': 2015,
        'releases': [2010, 2012, 2016],
        'size': 0,
        'authors': ['NetBiol Group'],
        'label': 'SignaLink',
        'color': '',
        'pubmeds': [20542890, 23331499],
        'urls': {
            'webpages': ['http://signalink.org/'],
            'articles': [
                'http://bioinformatics.oxfordjournals.org/content/26/16/2042.long',
                'http://www.biomedcentral.com/1752-0509/7/7'
            ],
            'omictools': ['http://omictools.com/signalink-tool']
        },
        'taxons': ['human', 'D. melanogaster', 'C. elegans'],
        'annot': ['pathway'],
        'recommend':
        'one of the largest resources with effect sign; due to its specific, biochemically defined pathways suitable for cross-talk analysis',
        'descriptions': [
            u'''
            In each of the three organisms, we first listed signaling proteins and interactions from reviews (and from WormBook in C.elegans) and then added further signaling interactions of the listed proteins. To identify additional interactions in C.elegans, we examined all interactions (except for transcription regulation) of the signaling proteins listed in WormBase and added only those to SignaLink that we could manually identify in the literature as an experimentally verified signaling interaction. For D.melanogaster, we added to SignaLink those genetic interactions from FlyBase that were also reported in at least one yeast-2-hybrid experiment. For humans, we manually checked the reliability and directions for the PPIs found with the search engines iHop and Chilibot.
            SignaLink assigns proteins to signaling pathways using the full texts of pathway reviews (written by pathway experts). While most signaling resources consider 5–15 reviews per pathway, SignaLink uses a total of 170 review papers, i.e. more than 20 per pathway on average. Interactions were curated from a total of 941 articles (PubMed IDs are available at the website). We added a small number of proteins based on InParanoid ortholog clusters. For curation, we used a self-developed graphical tool and Perl/Python scripts. The current version of SignaLink was completed in May 2008 based on WormBase (version 191), FlyBase (2008.6), Ensembl, UniProt and the publications listed on the website.
            The curation protocol of SignaLink (Fig. 1A) contains several steps aimed specifically at reducing data and curation errors. We used reviews as a starting point, manually looked up interactions three times, and manually searched for interactions of known signaling proteins with no signaling interactions so far in the database.
            '''
        ],
        'notes': [
            u'''
            For OmniPath we used the literature curated part of version 3 of SignaLink, which is unpublished yet. Version 2 is publicly available, and format definitions in pypath exist to load the version 2 alternatively.
            '''
        ],
        'emails': [('korcsmaros@gmail.com', 'Tamas Korcsmaros'),
                   ('tamas.korcsmaros@tgac.ac.uk', 'Tamas Korcsmaros')],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'license': {
            'name': 'CC-Attribution-NonCommercial-ShareAlike-3.0',
            'url': 'http://creativecommons.org/licenses/by-nc-sa/3.0/',
            'commercial_use': False
        },
        'pathguide': 320
    },
    'NRF2ome': {
        'year': 2013,
        'releases': [2013],
        'size': {
            'nodes': None,
            'edges': None,
        },
        'authors': ['NetBiol Group'],
        'label': 'NRF2ome',
        'color': '',
        'urls': {
            'webpages': ['http://nrf2.elte.hu/'],
            'articles': [
                'http://www.hindawi.com/journals/omcl/2013/737591/',
                'http://www.sciencedirect.com/science/article/pii/S0014579312003912'
            ]
        },
        'pubmeds': [22641035, 23710289],
        'taxons': ['human'],
        'recommend':
        'specific details about NRF2 related oxidative stress signaling; connections to transcription factors',
        'descriptions': [
            u'''
            From Korcsmaros 2010: ... we first listed signaling proteins and interactions from reviews and then added further signaling interactions of the listed proteins. We used reviews as a starting point, manually looked up interactions three times, and manually searched for interactions of known signaling proteins with no signaling interactions so far in the database.
            '''
        ],
        'emails': [('korcsmaros@gmail.com', 'Tamas Korcsmaros'),
                   ('tamas.korcsmaros@tgac.ac.uk', 'Tamas Korcsmaros')],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'license': {
            'name': 'CC-Attribution-NonCommercial-ShareAlike-3.0',
            'url': 'http://creativecommons.org/licenses/by-nc-sa/3.0/',
            'commercial_use': False
        }
    },
    'ARN': {
        'year': 2014,
        'releases': [2014],
        'size': 0,
        'authors': ['NetBiol Group'],
        'label': 'ARN',
        'color': '',
        'pubmeds': [25635527],
        'urls': {
            'webpages': ['http://autophagyregulation.org/'],
            'articles': [
                'http://www.tandfonline.com/doi/full/10.4161/15548627.2014.994346'
            ]
        },
        'taxons': ['human'],
        'annot': ['pathway'],
        'recommend':
        'well curated essential interactions in autophagy regulation; connections to transcription factors',
        'descriptions': [
            u'''
            From Korcsmaros 2010: ... we first listed signaling proteins and interactions from reviews and then added further signaling interactions of the listed proteins. We used reviews as a starting point, manually looked up interactions three times, and manually searched for interactions of known signaling proteins with no signaling interactions so far in the database.
            '''
        ],
        'emails': [('korcsmaros@gmail.com', 'Tamas Korcsmaros'),
                   ('tamas.korcsmaros@tgac.ac.uk', 'Tamas Korcsmaros')],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'license': {
            'name': 'CC-Attribution-NonCommercial-ShareAlike-3.0',
            'url': 'http://creativecommons.org/licenses/by-nc-sa/3.0/',
            'commercial_use': False
        }
    },
    'HPRD': {
        'year': 2010,
        'releases': [2002, 2005, 2009, 2010],
        'urls': {
            'webpages': ['http://www.hprd.org/'],
            'articles': [
                'http://genome.cshlp.org/content/13/10/2363.long',
                'http://nar.oxfordjournals.org/content/34/suppl_1/D411.long',
                'http://nar.oxfordjournals.org/content/37/suppl_1/D767.long'
            ],
            'omictools':
            ['http://omictools.com/human-protein-reference-database-tool']
        },
        'annot': ['mechanism'],
        'recommend':
        'one of the largest kinase-substrate resources; provides large amount of specific information; discontinued',
        'pubmeds': [14525934, 16381900, 18988627],
        'descriptions': [
            u'''
            The information about protein-protein interactions was cataloged after a critical reading of the published literature. Exhaustive searches were done based on keywords and medical subject headings (MeSH) by using Entrez. The type of experiments that served as the basis for establishing protein-protein interactions was also annotated. Experiments such as coimmunoprecipitation were designated in vivo, GST fusion and similar “pull-down” type of experiments were designated in vitro, and those identified by yeast two-hybrid were annotated as yeast two-hybrid.
            Posttranslational modifications were annotated based on the type of modification, site of modification, and the modified residue. In addition, the upstream enzymes that are responsible for modifications of these proteins were reported if described in the articles. The most commonly known and the alternative subcellular localization of the protein were based on the literature. The sites of expression of protein and/or mRNA were annotated based on published studies.
            '''
        ],
        'full_name': 'Human Protein Reference Database',
        'emails': [('pandey@jhmi.edu', 'Akhilesh Pandey')],
        'type': 'literature curated',
        'subtype': 'post-translational modification',
        'omnipath': True,
        'license': {
            'name':
            'No license. Everything in HPRD is free as long as it is not used for commercial purposes. Commercial entitites will have to pay a fee under a licensing arrangement which will be used to make this database even better. Commercial users should send an e-mail for details. This model of HPRD is similar to the SWISS-PROT licensing arrangement. We do not have any intentions to profit from HPRD. Our goal is to promote science by creating the infrastructure of HPRD. We hope to keep it updated with the assistance of the entire biomedical community. Any licensing fee, if generated, will be used to annotate HPRD better and to add more entries and features.',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        },
        'pypath': {
            'data': [
                'urls.urls[\'hprd_all\'][\'url\']',
                'urls.urls[\'hprd_all\'][\'ptm_file\']'
            ],
            'input': ['pypath.dataio.get_hprd()'],
            'intr': ['pypath.dataio.hprd_interactions()'],
            'format': [
                'pypath.data_formats.ptm[\'hprd\']',
                'pypath.data_formats.omnipath[\'hprd\']'
            ],
            'ptm': ['pypath.dataio.get_hprd_ptms()']
        },
        'pathguide': 14
    },
    'ACSN': {
        'year': 2015,
        'releases': [2008, 2014, 2015, 2016],
        'authors': ['Curie'],
        'urls': {
            'webpages': ['https://acsn.curie.fr'],
            'articles': [
                'http://www.nature.com/oncsis/journal/v4/n7/full/oncsis201519a.html',
                'http://msb.embopress.org/content/4/1/0174.long'
            ]
        },
        'pubmeds': [26192618, 18319725],
        'taxons': ['human'],
        'recommend':
        'the third largest process description resource; focused on signaling pathways; relationships mapped into a topological space',
        'descriptions': [
            u'''
            The map curator studies the body of literature dedicated to the biological process or molecular mechanism of interest. The initial sources of information are the major review articles from high-impact journals that represent the consensus view on the studied topic and also provide a list of original references. The map curator extracts information from review papers and represents it in the form of biochemical reactions in CellDesigner. This level of details reflects the ‘canonical’ mechanisms. Afterwards, the curator extends the search and analyses original papers from the list provided in the review articles and beyond. This information is used to enrich the map with details from the recent discoveries in the field. The rule for confident acceptance and inclusion of a biochemical reaction or a process is the presence of sufficient evidences from more than two studies, preferably from different scientific groups. The content of ACSN is also verified and compared with publicly available databases such as REACTOME, KEGG, WikiPathways, BioCarta, Cell Signalling and others to ensure comprehensive representation of consensus pathways and links on PMIDs of original articles confirmed annotated molecular interactions.
            ''', u'''
            CellDesigner 3.5 version was used to enter biological facts from a carefully studied selection of papers (see the whole bibliography on the web site with Supplementary information). Whenever the details of a biological fact could not be naturally expressed with CellDesigner standard notations, it was fixed and some solution was proposed. For example, we added a notation (co‐factor) to describe all the components intervening in the transcription of genes mediated by the E2F family proteins.
            '''
        ],
        'emails': [('andrei.zinovyev@curie.fr', 'Andrei Zinovyev')],
        'type': 'literature curated',
        'subtype': 'reaction',
        'full_name': 'Atlas of Cancer Signalling Networks',
        'omnipath': False,
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True

        },
        'pypath': {
            'data': [
                'pypath.urls.urls[\'acsn\'][\'url\']',
                'pypath.data_formats.files[\'acsn\'][\'ppi\']',
                'pypath.data_formats.files[\'acsn\'][\'names\']',
                'pypath.urls.urls[\'acsn\'][\'biopax_l3\']'
            ],
            'input': [
                'pypath.dataio.get_acsn()', 'pypath.dataio.get_acsn_effects()',
                'pypath.dataio.acsn_biopax()',
                'pypath.pypath.PyPath().acsn_effects()'
            ],
            'intr': ['pypath.dataio.acsn_ppi()'],
            'format': [
                'pypath.data_formats.reaction[\'acsn\']',
                'pypath.data_formats.reaction_misc[\'acsn\']'
            ]
        }
    },
    'DeathDomain': {
        'year': 2012,
        'releases': [2011, 2012],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['Myoungji University'],
        'label': 'DeathDomain',
        'color': '',
        'taxons': ['human'],
        'pubmeds': [22135292],
        'urls': {
            'articles': ['http://nar.oxfordjournals.org/content/40/D1/D331'],
            'webpages': ['http://deathdomain.org/']
        },
        'license': {
            'name':
            'No license. Please cite the following paper when you use Death Domain database in your publications, which is very important to sustain our service: Kwon et al. 2012',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        },
        'emails': [('hyunho@ynu.ac.kr', 'Hyun Ho Park')],
        'files': {
            'articles': ['DeathDomain_Kwon2011.pdf'],
            'data': {
                'raw': ['deathdomain.tsv'],
                'processed': ['deathdomain.sif']
            }
        },
        'taxons': ['human'],
        'size': {
            'nodes': 99,
            'edges': 175
        },
        'identifiers': ['GeneSymbol'],
        'annot': ['experiment'],
        'recommend':
        'focused deep curation effort on death domain superfamily proteins; many specific relationships',
        'descriptions': [
            u'''
            The PubMed database was used as the primary source for collecting information and constructing the DD database. After finding synonyms for each of the 99 DD superfamily proteins using UniProtKB and Entrez Gene, we obtained a list of articles using each name of the proteins and its synonyms on a PubMed search, and we selected the articles that contained evidence for physical binding among the proteins denoted. We also manually screened information that was in other databases, such as DIP, IntAct, MINT, STRING and Entrez Gene. All of the 295 articles used for database construction are listed on our database website.
            '''
        ],
        'notes': [
            u'''
            Detailful dataset with many references. Sadly the data can be extracted only by parsing HTML. It doesn't mean more difficulty than parsing XML formats, just these are not intended to use for this purpose.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'data_integration': 'static',
        'pypath': {
            'data': ['pypath.data/dd_refs.csv'],
            'format': [
                'pypath.data_formats.pathway[\'death\']',
                'pypath.data_formats.omnipath[\'death\']'
            ]
        },
        'pathguide': 442
    },
    'TRIP': {
        'year': 2014,
        'releases': [2010, 2012],
        'urls': {
            'articles': [
                'http://www.plosone.org/article/info:doi/10.1371/journal.pone.0047165',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D356.full',
                r'http://link.springer.com/article/10.1007%2Fs00424-013-1292-2'
            ],
            'webpages': ['http://www.trpchannel.org'],
            'omictools': [
                'http://omictools.com/transient-receptor-potential-channel-interacting-protein-database-tool'
            ]
        },
        'label': 'TRIP',
        'full_name':
        'Mammalian Transient Receptor Potential Channel-Interacting Protein Database',
        'emails': [('jhjeon2@snu.ac.kr', 'Ju-Hong Jeon')],
        'size': {
            'nodes': 468,
            'edges': 744
        },
        'pubmeds': [20851834, 23071747, 23677537],
        'files': {
            'articles': ['TRIP_Shin2012.pdf'],
            'data': {
                'raw': [],
                'processed': ['trip.sif']
            }
        },
        'taxons': ['human', 'mouse', 'rat'],
        'identifiers': ['GeneSymbol'],
        'recommend':
        'high number of specific interactions; focused on TRP channels',
        'descriptions': [
            u'''
            The literature on TRP channel PPIs found in the PubMed database serve as the primary information source for constructing the TRIP Database. First, a list of synonyms for the term ‘TRP channels’ was constructed from UniprotKB, Entrez Gene, membrane protein databases (Supplementary Table S2) and published review papers for nomenclature. Second, using these synonyms, a list of articles was obtained through a PubMed search. Third, salient articles were collected through a survey of PubMed abstracts and subsequently by search of full-text papers. Finally, we selected articles that contain evidence for physical binding among the proteins denoted. To prevent omission of relevant papers, we manually screened information in other databases, such as DIP, IntAct, MINT, STRING, BioGRID, Entrez Gene, IUPHAR-DB and ISI Web of Knowledge (from Thomson Reuters). All 277 articles used for database construction are listed in our database website.
            '''
        ],
        'notes': [
            u'''
            Good manually curated dataset focusing on TRP channel proteins, with ~800 binary interactions. The provided formats are not well suitable for bioinformatics use because of the non standard protein names, with greek letters and only human understandable formulas. Using HTML processing from 5-6 different tables, with couple hundreds lines of code, one have a chance to compile a usable table.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'pypath': {
            'input': [
                'pypath.dataio.trip_process()', 'pypath.dataio.take_a_trip()',
                'pypath.dataio.trip_process_table()',
                'pypath.dataio.trip_get_uniprot()',
                'pypath.dataio.trip_find_uniprot()'
            ],
            'intr': ['pypath.dataio.trip_interactions()'],
            'data': [
                'pypath.urls.urls[\'trip\'][\'intr\']',
                'pypath.urls.urls[\'trip\'][\'show\']',
                'pypath.urls.urls[\'trip\'][\'json\']',
            ],
            'format': [
                'pypath.data_formats.pathway[\'trip\']',
                'pypath.data_formats.omnipath[\'trip\']'
            ]
        },
        'license': {
            'name': 'CC-Attribution-ShareAlike-3.0',
            'url': 'http://creativecommons.org/licenses/by-nc-sa/3.0/',
            'commercial_use': True
        },
        'pathguide': 409
    },
    'Awan2007': {
        'year': 2007,
        'size': 0,
        'authors': ['Wang Group'],
        'label': 'Awan 2007',
        'color': '',
        'data_import': ['BioCarta', 'CA1'],
        'contains': ['BioCarta', 'CA1'],
        'urls': {
            'articles':
            ['http://www.cancer-systemsbiology.org/Papers/iet-sb2007.pdf']
        },
        'emails': [('Edwin.Wang@cnrc-nrc.gc.ca', 'Edwin Wang')],
        'pubmeds': [17907678],
        'descriptions': [
            u'''
            To construct the human cellular signalling network, we manually curated signalling pathways from literature. The signalling data source for our pathways is the BioCarta database (http://www.biocarta.com/genes/allpathways.asp), which, so far, is the most comprehensive database for human cellular signalling pathways. Our curated pathway database recorded gene names and functions, cellular locations of each gene and relationships between genes such as activation, inhibition, translocation, enzyme digestion, gene transcription and translation, signal stimulation and so on. To ensure the accuracy and the consistency of the database, each referenced pathway was cross-checked by different researchers and finally all the documented pathways were checked by one researcher. In total, 164 signalling pathways were documented (supplementary Table 2). Furthermore, we merged the curated data with another literature-mined human cellular signalling network. As a result, the merged network contains nearly 1100 proteins (SupplementaryNetworkFile). To construct a signalling network, we considered relationships of proteins as links (activation or inactivation as directed links and physical interactions in protein complexes as neutral links) and proteins as nodes.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': False,
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'Cui2007': {
        'year': 2007,
        'authors': ['Wang Group'],
        'label': 'Cui 2007',
        'color': '',
        'data_import': ['Awan2007', 'CancerCellMap'],
        'pubmeds': [18091723],
        'contains': ['Awan2007', 'CancerCellMap', 'CA1', 'BioCarta'],
        'urls': {
            'articles': ['http://msb.embopress.org/content/3/1/152'],
            'webpages': []
        },
        'emails': [('Edwin.Wang@cnrc-nrc.gc.ca', 'Edwin Wang')],
        'files': {
            'articles': ['Cui2007.pdf'],
            'data': {
                'raw': ['cui-network.xls', ],
                'processed': ['cui.sif']
            }
        },
        'identifiers': ['EntrezGene'],
        'size': {
            'edges': 4249,
            'nodes': 1528
        },
        'taxons': ['human'],
        'descriptions': [
            u'''
            To build up the human signaling network, we manually curated the signaling molecules (most of them are proteins) and the interactions between these molecules from the most comprehensive signaling pathway database, BioCarta (http://www.biocarta.com/). The pathways in the database are illustrated as diagrams. We manually recorded the names, functions, cellular locations, biochemical classifications and the regulatory (including activating and inhibitory) and interaction relations of the signaling molecules for each signaling pathway. To ensure the accuracy of the curation, all the data have been crosschecked four times by different researchers. After combining the curated information with another literature‐mined signaling network that contains ∼500 signaling molecules (Ma'ayan et al, 2005)[this is the CA1], we obtained a signaling network containing ∼1100 proteins (Awan et al, 2007). We further extended this network by extracting and adding the signaling molecules and their relations from the Cancer Cell Map (http://cancer.cellmap.org/cellmap/), a database that contains 10 manually curated signaling pathways for cancer. As a result, the network contains 1634 nodes and 5089 links that include 2403 activation links (positive links), 741 inhibitory links (negative links), 1915 physical links (neutral links) and 30 links whose types are unknown (Supplementary Table 9). To our knowledge, this network is the biggest cellular signaling network at present.
            ''', u'''
            From Awan 2007: To construct the human cellular signalling network, we manually curated signalling pathways from literature. The signalling data source for our pathways is the BioCarta database (http://www.biocarta.com/genes/allpathways.asp), which, so far, is the most comprehensive database for human cellular signalling pathways. Our curated pathway database recorded gene names and functions, cellular locations of each gene and relationships between genes such as activation, inhibition, translocation, enzyme digestion, gene transcription and translation, signal stimulation and so on. To ensure the accuracy and the consistency of the database, each referenced pathway was cross-checked by different researchers and finally all the documented pathways were checked by one researcher.
            '''
        ],
        'notes': [
            u'''
            Excellent signaling network with good topology for all those who doesn't mind to use data of unknown origin. Supposedly a manually curated network, but data files doesn't include article references. Merging CA1 network with CancerCellMap and BioCarta (also without references) makes the origin of the data untraceable.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': False,
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'BioCarta': {
        'year': 2006,
        'releases': [2006],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['Community'],
        'label': 'BioCarta',
        'color': '',
        'urls': {
            'webpages': [
                'http://www.biocarta.com/',
                'http://cgap.nci.nih.gov/Pathways/BioCarta_Pathways'
            ],
            'articles': []
        },
        'emails':
        [('info@biocarta.com', 'BioCarta Scientific Advisory Board')],
        'taxons': ['human'],
        'descriptions': [
            u'''
            Community built pathway database based on expert curation.
            '''
        ],
        'notes': [
            u'''
            This resource includes a huge number of pathways, each curated by experts from a few reviews. The data is not available for download from the original webpage, only from second hand, for example from NCI-PID, in NCI-XML format. However, these files doesn't contain any references, which makes problematic the use of the BioCarta dataset. Also, some pathways are reviewed long time ago, possibly outdated.
            ''', '''
            The Company and the website looks like was abandoned around 2003-2006.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': False,
        'license': {
            'name':
            'BioCarta webpage Terms and Conditions of Use (pathways are not owned by BioCarta and are free to use)',
            'url':
            'http://web.archive.org/web/20150207091158/http://biocarta.com/legal/terms.asp',
            'commercial_use': False
        }
    },
    'TLR': {
        'urls': {
            'articles':
            ['https://www.embopress.org/doi/10.1038/msb4100057']
        },
        'type': 'literature curated',
        'subtype': 'model',
        'omnipath': False,
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'CA1': {
        'year': 2005,
        'releases': [2005],
        'size': {
            'nodes': 545,
            'edges': 1259
        },
        'authors': ['Iyengar Lab'],
        'label': 'Ma\'ayan 2005',
        'full_name': 'Human Hippocampal CA1 Region Neurons Signaling Network',
        'color': '',
        'pubmeds': [16099987],
        'urls': {
            'articles':
            ['http://www.sciencemag.org/content/309/5737/1078.full'],
            'webpages': []
        },
        'emails': [('ravi.iyengar@mssm.edu', 'Ravi Iyengar')],
        'taxons': ['human', 'mouse'],
        'annot': ['localization', 'mechanism'],
        'recommend':
        'among the largest resources providing effect signs; more than a decade old dataset, used in many studies',
        'descriptions': [
            u'''
            We used published research literature to identify the key components of signaling pathways and cellular machines, and their binary interactions. Most components (~80%) have been described in hippocampal neurons or related neuronal cells. Other components are from other cells, but are included because they are key components in processes known to occur in hippocampal neurons, such as translation. We then established that these interactions were both direct and functionally relevant. All of the connections were individually verified by at least one of the authors of this paper by reading the relevant primary paper(s). We developed a system made of 545 components (nodes) and 1259 links (connections). We used arbitrary but consistent rules to sort components into various groups. For instance, transcription factors are considered a as part of the transcriptional machinery, although it may also be equally valid to consider them as the most downstream component of the central signaling network. Similarly the AMPA receptor-channel (AMPAR) is considered part of the ion channels in the electrical response system since its activity is essential to defining the postsynaptic response, although it binds to and is activated by glutamate, and hence can be also considered a ligand gated receptor-channel in the plasma membrane. The links were specified by two criteria: function and biochemical mechanism. Three types of functional links were specified. This follows the rules used for representation of pathways in Science’s STKE (S1). Links may be activating, inhibitory or neutral. Neutral links do not specify directionality between components, and are mostly used to represent scaffolding and anchoring undirected or bidirectional interactions. The biochemical specification includes defining the reactions as non-covalent binding interactions or enzymatic reactions. Within the enzymatic category, reactions were further specified as phosphorylation, dephosphorylation, hydrolysis, etc. These two criteria for specification are independent and were defined for all interactions. For the analyses in this study we only used the functional criteria: activating, inhibitory or neutral specifications. We chose papers that demonstrated direct interactions that were supported by either biochemical or physiological effects of the interactions. From these papers we identified the components and interactions that make up the system we analyzed. During this specification process we did not consider whether these interactions would come together to form higher order organizational units. Each component and interaction was validated by a reference from the primary literature (1202 papers were used). A list of authors who read the papers to validate the components and interactions is provided under authors contributions.
            '''
        ],
        'notes': [
            u'''
            One of the earliest manually curated networks, available in easily accessible tabular format, including UniProt IDs and PubMed references.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'data_integration': 'dynamic',
        'pypath': {
            'intr': ['pypath.urls.urls[\'ca1\'][\'url\']'],
            'data': ['pypath.dataio.get_ca1()'],
            'format': [
                'pypath.data_formats.pathway[\'ca1\']',
                'pypath.data_formats.omnipath[\'ca1\']'
            ]
        },
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'CancerCellMap': {
        'year': 2006,
        'relases': [2004, 2006],
        'urls': {
            'articles': [],
            'webpages': [
                'http://www.pathwaycommons.org/archives/PC1/last_release-2011/tab_delim_network/by_source/',
                'http://web.archive.org/web/20130729025707/http://cancer.cellmap.org/cellmap/home.do'
            ]
        },
        'authors': ['Bader Lab'],
        'emails': [('gary.bader@utoronto.ca', 'Gary Bader')],
        'descriptions': [
            u'''
            Manually curated data, unpublished. A team of M.Sc. and Ph.D. biologists at the Institute of Bioinformatics in Bangalore, India read original research papers and hand-entered the pathway data into our database. The quality of the Cancer Cell Map pathways is very high. Half of the pathways were reviewed by experts at Memorial Sloan-Kettering Cancer Center and were found to contain only a few errors, which were subsequently fixed. A pathway is a collection of all genes/proteins that have been described as pathway members in any publication and all the interactions between them that can be found described in the literature.
            '''
        ],
        'notes': [
            u'''
            One of the earliest manually curated datasets, now only available from second hand, e.g. from PathwayCommons. Included in many other resources. Contains binary interactions with PubMed references.
            '''
        ],
        'taxons': ['human', 'mouse', 'rat'],
        'annot': ['localization'],
        'recommend':
        'old literature curated resource; effect signs for some interactions',
        'type': 'literature curated',
        'subtype': 'interaction',
        'omnipath': True,
        'pypath': {
            'intr': ['pypath.dataio.get_ccmap()'],
            'data': [
                'pypath.urls.urls[\'ccmap\'][\'nodes\']',
                'pypath.urls.urls[\'ccmap\'][\'edges\']'
            ],
            'format': [
                'pypath.data_formats.interaction[\'ccmap\']',
                'pypath.data_formats.omnipath[\'ccmap\']'
            ]
        },
        'data_integration': 'dynamic',
        'license': {
            'name': 'CC-Attribution-2.5',
            'url': 'http://creativecommons.org/licenses/by/2.5/',
            'commercial_use': True
        },
        'pathguide': 223
    },
    'CARFMAP': {
        'year': 2015,
        'urls': {
            'articles': [
                'http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0143274'
            ],
            'webpages': ['http://visionet.erc.monash.edu.au/CARFMAP/']
        },
        'pubmeds': [26673252],
        'emails': [('hieu.nim@monash.edu', 'Hieu T Nim'),
                   ('sarah.boyd@monash.edu', 'Sarah E Boyd')],
        'license': {
            'name': 'CC-Attribution-4.0',
            'url': 'http://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        },
        'omnipath': False,
        'type': 'Literature curated',
        'subtype': 'Pathway'
    },
    'HSN': {
        'year': 2014,
        'releases': [2009, 2010, 2011, 2012, 2013, 2014],
        'nodes': 6300,
        'edges': 63000,
        'authors': ['Wang Group'],
        'full_name': 'Human Signaling Network version 6',
        'label': 'HumanSignalingNetwork',
        'color': '',
        'data_import': ['Cui2007', 'BioCarta', 'CST', 'NCI-PID', 'iHOP'],
        'urls': {
            'webpages': [
                'http://www.cancer-systemsbiology.org/dataandsoftware.htm',
                'http://www.bri.nrc.ca/wang/'
            ],
            'articles': []
        },
        'emails': [('Edwin.Wang@cnrc-nrc.gc.ca', 'Edwin Wang')],
        'taxons': ['human', 'mouse', 'rat'],
        'contains': [
            'Cui2007', 'CancerCellMap', 'Awan2007', 'NetPath', 'CA1',
            'NCI-PID', 'BioCarta', 'CST', 'iHOP'
        ],
        'descriptions': [
            u'''
            Composed from multiple manually curated datasets, and contains own manual cuartion effort. Methods are unclear, and the dataset has not been published in reviewed paper. Based on the Cui et al 2007.
            Wang Lab has manually curated human signaling data from literature since 2005. The data sources include BioCarta, CST Signaling pathways, NCI Pathway Interaction Database, iHOP, and many review papers. The contents are updated every year.
            iHOP is not literature curated, but is a literature mining platform.
            '''
        ],
        'notes': [
            u'''
            This network aims to merge multiple manually curated networks. Unfortunately a precise description of the sources and methods is missing. Also, the dataset does not include the references. Moreover, the data file misses header and key, so users can only guess about the meaning of columns and values.
            '''
        ],
        'type': 'literature curated',
        'data_integration': 'dynamic',
        'subtype': 'activity flow',
        'omnipath': False,
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        },
        'pypath': {
            'data': ['pypath.urls.urls[\'hsn\']'],
            'format': ['pypath.data_formats.interaction_misc[\'hsn\']'],
            'intr': ['pypath.dataio.get_hsn()']
        }
    },
    'Ataxia': {
        'year': 2010,
        'urls': {
            'webpages': ['http://franklin.imgen.bcm.tmc.edu/ppi/tables.html'],
            'articles':
            ['http://hmg.oxfordjournals.org/content/20/3/510.long']
        },
        'taxons': ['human'],
        'authors': ['Shaw Lab'],
        'descriptions': [
            u'''
            In order to expand the interaction dataset, we added relevant direct protein–protein interactions from currently available human protein–protein interaction networks (Rual et al., 2005; Stelzl et al., 2005). We also searched public databases, including BIND (Bader et al., 2003), DIP (Xenarios et al., 2002), HPRD (Peri et al., 2003), MINT (Zanzoni et al., 2002), and MIPS (Pagel et al., 2005), to identify literature-based binary interactions involving the 54 ataxia-associated baits and the 561 interacting prey proteins. We identified 4796 binary protein–protein interactions for our Y2H baits and prey proteins (Table S4) and incorporated them in the Y2H protein–protein interaction map (Figures 4A–4C).
            '''
        ],
        'notes': [
            u'''
            The Ataxia network doesn't contain original manual curation effort. The integrated data are very old.
            '''
        ],
        'type': 'high-throughput',
        'subtype': 'interaction',
        'omnipath': False,
        'emails': [('Tong_Hao@dfci.harvard.edu', 'Tong Hao'),
                   ('barabasi@gmail.com', 'Albert-Laszlo Barabasi')],
        'license': {
            'url': 'http://creativecommons.org/licenses/by-nc/2.5',
            'name': 'CC-Attribution-2.5',
            'commercial_use': True
        }
    },
    'Reactome': {
        'year': 2016,
        'releases': [2004, 2008, 2010, 2012, 2014, 2015, 2016],
        'urls': {
            'webpages': ['http://reactome.org/'],
            'articles': [
                'http://nar.oxfordjournals.org/content/33/suppl_1/D428.long',
                'http://nar.oxfordjournals.org/content/37/suppl_1/D619.long',
                'http://onlinelibrary.wiley.com/doi/10.1002/pmic.201100066/abstract',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D691.long',
                'http://genomebiology.com/content/8/3/R39',
                'http://nar.oxfordjournals.org/content/42/D1/D472.long',
                'http://www.mdpi.com/2072-6694/4/4/1180/htm'
            ],
            'omictools': ['http://omictools.com/reactome-tool']
        },
        'pubmeds':
        [15608231, 18981052, 21751369, 21067998, 24213504, 24243840],
        'annot': ['localization', 'mechanism', 'pathway'],
        'recommend':
        'the largest process description resource; peerless coverage; both signaling and metabolism',
        'descriptions': [
            u'''
            Once the content of the module is approved by the author and curation staff, it is peer-reviewed on the development web-site, by one or more bench biologists selected by the curator in consultation with the author. The peer review is open and the reviewers are acknowledged in the database by name. Any issues raised in the review are resolved, and the new module is scheduled for release.
            '''
        ],
        'notes': [
            u'''
            No binary interactions can be exported programmatically from any format of the Reactome dataset. Reactome's curation method doesn't cover binary interactions, the inferred lists on the webpage are based on automatic expansion of complexes and reactions, and thus are unreliable. In lack of information, references cannot be assigned to interactions.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'process description',
        'data_integration': 'dynamic',
        'omnipath': False,
        'emails': [('help@reactome.org', 'Reactome Team'),
                   ('hhe@ebi.ac.uk', 'Henning Hermjakob')],
        'license': {
            'name': 'CC-Attribution-4.0',
            'url': 'http://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        },
        'pypath': {
            'data': [
                'pypath.urls.urls[\'reactome\'][\'sbml\']',
                'pypath.urls.urls[\'reactome\'][\'biopax_l2\']',
                'pypath.urls.urls[\'reactome\'][\'biopax_l3\']'
            ],
            'format': ['pypath.data_formats.reaction_misc[\'reactome\']'],
            'intr': ['pypath.dataio.reactome_interactions()'],
            'input': [
                'pypath.dataio.reactome_sbml()',
                'pypath.dataio.reactome_biopax()'
            ]
        },
        'pathguide': 103
    },
    'Li2012': {
        'year': 2012,
        'urls': {
            'articles': ['http://genome.cshlp.org/content/22/7/1222']
        },
        'pubmeds': [22194470],
        'label': 'Li 2012',
        'authors': ['Wang Lab'],
        'emails': [('edwin.wang@cnrc-nrc.gc.ca', 'Edwin Wang')],
        'taxons': ['human'],
        'descriptions': [
            u'''
            Human phosphotyrosine signaling network.
            We manually collected the experimentally determined human TK–substrate interactions and substrate–SH2/PTB domain interactions from the literature (see Supplemental Materials), as well as the Phospho.ELM and PhosphoSitePlus databases. [71 references, 585 circuits]
            '''
        ],
        'type': 'high-throughput',
        'data_integration': 'static',
        'subtype': 'yeast 2 hybrid',
        'omnipath': False,
        'pypath': {
            'input': ['pypath.dataio.get_li2012()'],
            'data': [
                'pypath.urls.urls[\'li2012\'][\'file\']'
                'pypath.data/li2012.csv'
            ],
            'intr': ['pypath.dataio.li2012_interactions()'],
            'dmi': ['pypath.dataio.li2012_dmi()'],
            'ptm': ['pypath.dataio.li2012_phospho()'],
            'format': ['pypath.data_formats.ptm_misc[\'li2012\']']
        },
        'license': {
            'name': 'No license.',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'Zaman2013': {
        'year': 2013,
        'urls': {
            'articles': [
                'http://www.sciencedirect.com/science/article/pii/S2211124713004695'
            ],
            'webpages': []
        },
        'pubmeds': [24075989],
        'label': 'Zaman 2013',
        'authors': ['Wang Lab'],
        'emails': [('edwin.wang@cnrc-nrc.gc.ca', 'Edwin Wang')],
        'contains': [
            'HumanSignalingNetwork', 'Cui2007', 'CA1', 'BioCarta', 'Awan2007',
            'Li2012', 'NCI-PID'
        ],
        'descriptions': [
            u'''
            The human signaling network (Version 4, containing more than 6,000 genes and more than 50,000 relations) includes our previous data obtained from manually curated signaling networks (Awan et al., 2007; Cui et al., 2007; Li et al., 2012) and by PID (http://pid.nci.nih.gov/) and our recent manual curations using the iHOP database (http://www.ihop-net.org/UniPub/iHOP/).
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': False,
        'license': {
            'name': 'No license.',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        }
    },
    'AlzPathway': {
        'year': 2015,
        'releases': [2012, 2015],
        'size': 0,
        'authors': ['Tokyo Bioinf'],
        'label': 'AlzPathway',
        'color': '',
        'urls': {
            'articles': [
                'http://www.biomedcentral.com/1752-0509/6/52',
                'http://onlinelibrary.wiley.com/doi/10.1038/clpt.2013.37/epdf',
                r'http://link.springer.com/protocol/10.1007%2F978-1-4939-2627-5_25'
            ],
            'webpages': ['http://alzpathway.org/AlzPathway.html']
        },
        'pubmeds': [23511713, 22647208],
        'emails': [('ogishima@sysmedbio.org', 'Soichi Ogishima'),
                   ('info@alzpathway.org', 'Soichi Ogishima')],
        'license':
        u'''Licensed under a Creative Commons Attribution 3.0 Unported License. Cite Mizuno et al''',
        'descriptions': [
            u'''
            We collected 123 review articles related to AD accessible from PubMed. We then manually curated these review articles, and have built an AD pathway map by using CellDesigner. Molecules are distinguished by the following types: proteins, complexes, simple molecules, genes, RNAs, ions, degraded products, and phenotypes. Gene symbols are pursuant to the HGNC symbols. Reactions are also distinguished by the following categories: state transition, transcription, translation, heterodimer association, dissociation, transport, unknown transition, and omitted transition. All the reactions have evidences to the references in PubMed ID using the MIRIAM scheme. All the references used for constructing the AlzPathway are listed in the ‘References for AlzPathway’. Cellular types are distinguished by the followings: neuron, astrocyte, and microglial cells. Cellular compartments are also distinguished by the followings: brain blood barrier, presynaptic, postsynaptic, and their inner cellular localizations.
            '''
        ],
        'notes': [
            u'''
            References can be fetched only from XML formats, not from the SIF file. Among approx. 150 protein-protein interactions, also contains interactions of many small molecules, denoted by pubchem IDs.
            '''
        ],
        'type': 'literature curated',
        'data_integration': 'static',
        'subtype': 'activity flow',
        'omnipath': False,
        'pypath': {
            'data': ['pypath.data/alzpw-ppi.csv'],
            'format': ['pypath.data_formats.interaction[\'alz\']']
        },
        'license': {
            'name': 'CC-Attribution-3.0',
            'url': 'http://creativecommons.org/licenses/by/3.0/',
            'commercial_use': True
        }
    },
    'MPPI': {
        'year': 2005,
        'releases': [2000, 2005],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['MIPS Munich'],
        'license': {
            'name':
            'No license. "You are free to use the database as you please including full download of'
            ' the dataset for your own analyses as long as you cite the source properly (Pagel et al. 2005)."',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False # To be clarified
        },
        'label': 'MPPI',
        'full_name': 'The MIPS Mammalian Protein-Protein Interaction Database',
        'urls': {
            'articles':
            ['http://bioinformatics.oxfordjournals.org/content/21/6/832'],
            'webpages': ['http://mips.helmholtz-muenchen.de/proj/ppi/'],
            'omictools': [
                'http://omictools.com/mips-mammalian-protein-protein-interaction-database-tool'
            ]
        },
        'emails': [('p.pagel@wzw.tum.de', 'Philipp Pagel'),
                   ('d.frishman@helmholtz-muenchen.de', 'Dmitrij Frishman')],
        'pubmeds': [15531608],
        'taxons': ['human', 'mammalia'],
        'data_integration': 'static',
        'recommend':
        'small, literature curated interaction dataset; complements well activity flow resources',
        'descriptions': [
            u'''
            The first and foremost principle of our MPPI database is to favor quality over completeness. Therefore, we decided to include only published experimental evidence derived from individual experiments as opposed to large-scale surveys. High-throughput data may be integrated later, but will be marked to distinguish it from evidence derived from individual experiments.
            '''
        ],
        'notes': [
            u'''
            This database contains hundreds of interactions curated manually from original papers. The format is perfect, with UniProt IDs, and PubMed references.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'interaction',
        'omnipath': True,
        'pypath': {
            'data': ['pypath.data/mppi_human_rep.csv'],
            'format': [
                'pypath.data_formats.interaction[\'mppi\']',
                'pypath.data_formats.omnipath[\'mppi\']'
            ]
        },
        'pathguide': 319
    },
    'Negatome': {
        'year': 2013,
        'relases': [2009],
        'urls': {
            'articles': [
                'http://nar.oxfordjournals.org/content/38/suppl_1/D540.long',
                'http://nar.oxfordjournals.org/content/42/D1/D396.long'
            ],
            'webpages':
            ['http://mips.helmholtz-muenchen.de/proj/ppi/negatome/'],
            'omictools': ['http://omictools.com/negatome-tool']
        },
        'pubmeds': [24214996, 19920129],
        'emails': [('d.frishman@helmholtz-muenchen.de', 'Dmitrij Frishman')],
        'descriptions': [
            u'''
            Annotation of the manual dataset was performed analogous to the annotation of protein–protein interactions and protein complexes in previous projects published by our group. Information about NIPs was extracted from scientific literature using only data from individual experiments but not from high-throughput experiments. Only mammalian proteins were considered. Data from high-throughput experiments were omitted in order to maintain the highest possible standard of reliability.
            '''
        ],
        'type': 'literature curated',
        'data_integration': 'static',
        'subtype': 'negative',
        'omnipath': True,
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        },
        'pypath': {
            'data': ['pypath.data/negatome_manual.csv'],
            'format': ['pypath.data_formats.negative[\'negatome\']'],
            'misc': [
                'pypath.pypath.PyPath().negative_report()',
                'pypath.pypath.PyPath().apply_negative()'
            ]
        },
        'pathguide': 444
    },
    'Macrophage': {
        'year': 2010,
        'urls': {
            'articles': ['http://www.biomedcentral.com/1752-0509/4/63'],
            'webpages': ['http://www.macrophages.com/macrophage-pathways'],
            'omictools':
            ['http://omictools.com/macrophage-pathway-knowledgebase-tool']
        },
        'emails': [('tom.freeman@roslin.ed.ac.uk', 'Tom Freeman')],
        'pubmeds': [20470404],
        'annot': ['localization', 'mechanism'],
        'recommend':
        'medium size resource with effect signs; high ratio of specific interactions; interactions confirmed in macrophages',
        'descriptions': [
            u'''
            Ongoing analysis of macrophage-related datasets and an interest in consolidating our knowledge of a number of signalling pathways directed our choice of pathways to be mapped (see Figure 1). Public and propriety databases were initially used as resources for data mining, but ultimately all molecular interaction data was sourced from published literature. Manual curation of the literature was performed to firstly evaluate the quality of the evidence supporting an interaction and secondly, to extract the necessary and additional pieces of information required to 'understand' the pathway and construct an interaction diagram. We have drawn pathways based on our desire to model pathways active in a human macrophage and therefore all components have been depicted using standard human gene nomenclature (HGNC). However, our understanding of the pathway components and the interactions between them, have been drawn largely from a consensus view of literature knowledge. As such the pathways presented here are based on data derived from a range of different cellular systems and mammalian species (human and mouse).
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': True,
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/legalcode',
            'commercial_use': True
        },
        'data_integration': 'static',
        'pypath': {
            'data': ['pypath.data/macrophage_strict.csv'],
            'format': ['pypath.data_formats.pathway[\'macrophage\']']
        }
    },
    'NetPath': {
        'year': 2015,
        'releases': [2010, 2011, 2012, 2013, 2014, 2015],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['Pandey Lab', 'IOB Bangalore'],
        'emails': [('pandey@jhmi.edu', 'Akhilesh Pandey')],
        'license': {
            'name': 'CC-Attribution 2.5',
            'url': 'http://creativecommons.org/licenses/by/2.5/',
            'commercial_use': True
            },
        'label': 'NetPath',
        'color': '',
        'data_import': ['CancerCellMap'],
        'includes': ['CancerCellMap'],
        'annot': ['experiment', 'mechanism', 'pathway'],
        'urls': {
            'articles': [
                'http://genomebiology.com/content/11/1/R3',
                'http://database.oxfordjournals.org/content/2011/bar032.long',
                'http://database.oxfordjournals.org/content/2011/bar021.long',
                'http://www.omicsonline.com/0974-276X/JPB-04-184.php',
                'http://www.biomedcentral.com/1756-0500/4/408',
                'http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3830942/',
                'http://link.springer.com/article/10.1007/s12079-012-0168-0',
                'http://www.hindawi.com/journals/jst/2012/376470/',
                'http://link.springer.com/article/10.1007/s12079-012-0186-y',
                'http://link.springer.com/article/10.1007/s12079-012-0181-3',
                'http://link.springer.com/article/10.1007/s12079-013-0197-3',
                'http://www.tandfonline.com/doi/full/10.3109/15419061.2013.791683',
                'http://link.springer.com/article/10.1007/s12079-013-0200-z',
                'http://database.oxfordjournals.org/content/2014/bau007.long',
                'http://link.springer.com/article/10.1007/s12079-014-0224-z',
                'http://www.hindawi.com/journals/jst/2014/962962/'
            ],
            'webpages': ['http://netpath.org/']
        },
        'pubmeds': [
            20067622, 21959865, 21742767, 21996254, 24255551, 22684822,
            22649723, 23255051, 23161412, 23504413, 23631681, 23606317,
            24573880, 24584707, 24829797
        ],
        'recommend':
        'the smallest process description resource, but represents a high quality literature curation effort; provides pathway annotations',
        'descriptions': [
            u'''
            The initial annotation process of any signaling pathway involves gathering and reading of review articles to achieve a brief overview of the pathway. This process is followed by listing all the molecules that arereported to be involved in the pathway under annotation. Information regarding potential pathway authorities are also gathered at this initial stage. Pathway experts are involved in initial screening of the molecules listed to check for any obvious omissions. In the second phase, annotators manually perform extensive literature searches using search keys, which include all the alter native names of the molecules involved, the name of the pathway, the names of reactions, and so on. In addition, the iHOP resource is also used to perform advanced PubMed-based literature searches to collect the reactions that were reported to be implicated in a given pathway. The collected reactions are manually entered using the PathBuilder annotation interface, which is subjected to an internal review process involving PhD level scientists with expertise in the areas of molecular biology, immunology and biochemistry. However, there are instances where a molecule has been implicated in a pathway in a published report but the associated experimental evidence is either weak or differs from experiments carried out by other groups. For this purpose, we recruit several investigators as pathway authorities based on their expertise in individual signaling pathways. The review by pathway authorities occasionally leads to correction of errors or, more commonly, to inclusion of additional information. Finally, the pathway authorities help in assessing whether the work of all major laboratories has been incorporated for the given signaling pathway.
            '''
        ],
        'notes': [
            u'''
            Formats are unclear. The tab delimited format contains the pathway memberships of genes, PubMed references, but not the interaction partners! The Excel file is very weird, in fact it is not an excel table, and contains only a few rows from the tab file. The PSI-MI XML is much better. By writing a simple parser, a lot of details can be extracted.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'process description',
        'data_integration': 'dynamic',
        'omnipath': False,
        'license': {
            'name': 'CC-Attribution-2.5',
            'url': 'http://creativecommons.org/licenses/by/2.5/',
            'commercial_use': True
        },
        'pypath': {
            'input': ['pypath.dataio.netpath_names()'],
            'intr': ['pypath.dataio.netpath()'],
            'data': [
                'pypath.urls.urls[\'netpath_psimi\']',
                'pypath.urls.urls[\'netpath_names\']'
            ],
            'format': ['pypath.data_formats.interaction[\'netpath\']']
        },
        'pathguide': 315
    },
    'InnateDB': {
        'year': 2015,
        'urls': {
            'articles': [
                'http://msb.embopress.org/content/4/1/218.long',
                'http://www.biomedcentral.com/1752-0509/4/117',
                'http://nar.oxfordjournals.org/content/41/D1/D1228.long'
            ],
            'webpages': ['http://www.innatedb.com/'],
            'omictools': ['http://omictools.com/innatedb-tool']
        },
        'authors': ['Brinkman Lab', 'Hancock Lab', 'Lynn Group'],
        'emails': [('innatedb-mail@sfu.ca', 'InnateDB Team'),
                   ('david.lynn@sahmri.com', 'David Lynn')],
        'pubmeds': [20727158, 23180781, 18766178],
        'releases': [2008, 2010, 2013, 2014, 2015],
        'annot': ['experiment', 'mechanism'],
        'recommend':
        'middle size interaction resource; special focus on immune pathways',
        'descriptions': [
            u'''
            InnateDB (www.innatedb.com) is a database and integrated analysis platform specifically designed to facilitate systems-level analyses of the mammalian innate immune response (Lynn et al. 2008; 2010, 2013). To enrich our knowledge of innate immunity networks and pathways, the InnateDB curation team has contextually annotated >25,000 human and mouse innate immunity-relevant molecular interactions through the review of >5,000 biomedical articles. Curation adheres to the MIMIx guidelines and new interactions are added weekly. Importantly, interactions are curated between molecules with a documented role in an innate immunity relevant biological process or pathway and all other interactors regardless of whether the interacting molecule has any known role in innate immunity. This approach captures interactions between the innate immune system and other systems.

            InnateDB is not limited to data on the innate immune system. It is a comprehensive database of human, mouse and bovine molecular interactions and pathways, consisting of more than 300,000 molecular interactions and 3,000+ pathways, integrated from major public molecular interaction and pathway databases. InnateDB is also an analysis platform offering user-friendly bioinformatics tools, including pathway and ontology analysis, network visualization and analysis and the ability to upload and analyze user-supplied gene expression or other quantitative data in a network and/or pathway context. The platform has a global profile and is utilised by >10,000 users per annum and is widely cited. A mirror of the site hosted in Australia is also available at innatedb.sahmri.com.

            Note that new interactions and gene annotations are added to InnateDB on an almost weekly database so the data is being continuously updated.
            '''
        ],
        'notes': [
            u'''
            Probably the largest manually curated binary protein interaction dataset, developed by a dedicated full time team of curators. Formats are clear and accessible, comprising UniProt IDs, PubMed references, experimental evidences and mechanisms.
            '''
        ],
        'type': 'literature curated',
        'data_integration': 'static',
        'subtype': 'interaction',
        'omnipath': True,
        'license': {
            'name': 'Design Science License',
            'url': 'http://www.innatedb.com/license.jsp',
            'commercial_use': True
        },
        'pypath': {
            'format': [
                'pypath.data_formats.interaction[\'innatedb\']',
                'pypath.data_formats.omnipath[\'innatedb\']'
            ],
            'data': ['pypath.data/innatedb.csv']
        },
        'pathguide': 264
    },
    'CORUM': {
        'year': 2012,
        'releases': [2007, 2009],
        'urls': {
            'articles': [
                'http://nar.oxfordjournals.org/content/36/suppl_1/D646.long',
                'http://nar.oxfordjournals.org/content/38/suppl_1/D497.long'
            ],
            'webpages':
            ['http://mips.helmholtz-muenchen.de/genre/proj/corum'],
            'omictools': [
                'http://omictools.com/comprehensive-resource-of-mammalian-protein-complexes-tool'
            ]
        },
        'full_name': 'Comprehensive Resource of Mammalian protein complexes',
        'emails': [('andreas.ruepp@helmholtz-muenchen.de', 'Andreas Ruepp')],
        'pubmeds': [19884131, 17965090],
        'taxons': ['human', 'mouse', 'rat'],
        'descriptions': [
            u'''
            The CORUM database is a collection of experimentally verified mammalian protein complexes. Information is manually derived by critical reading of the scientific literature from expert annotators. Information about protein complexes includes protein complex names, subunits, literature references as well as the function of the complexes.
            In order to provide a high-quality dataset of mammalian protein complexes, all entries are manually created. Only protein complexes which have been isolated and characterized by reliable experimental evidence are included in CORUM. To be considered for CORUM, a protein complex has to be isolated as one molecule and must not be a construct derived from several experiments. Also, artificial constructs of subcomplexes are not taken into account. Since information from high-throughput experi ments contains a significant fraction of false-positive results, this type of data is excluded. References for relevant articles were mainly found in general review articles, cross-references to related protein complexes within analysed literature and comments on referenced articles in UniProt.
            '''
        ],
        'notes': [
            u'''
            CORUM is not part of the OmniPath pathways network, because we did not applied any complex expansion. But it has an interface built in the pypath module.
            '''
        ],
        'type': 'literature curated',
        'data_integration': 'dynamic',
        'subtype': 'complexes',
        'omnipath': False,
        'license': {
            'name': 'No license.',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        },
        'pypath': {
            'data': ['pypath.urls.urls[\'corum\'][\'url\']'],
            'input': [
                'pypath.dataio.get_corum()',
                'pypath.pypath.PyPath().load_corum()'
            ]
        },
        'pathguide': 322
    },
    'CST': {
        'year': 2015,
        'releases': [2005, 2015],
        'nodes': None,
        'edges': None,
        'authors': ['CST'],
        'label': 'CST Pathways',
        'color': '',
        'full_name': 'Cell Signaling Technology Pathways',
        'urls': {
            'articles': [],
            'webpages': [
                'http://www.cellsignal.com/common/content/content.jsp?id=science-pathways'
            ]
        },
        'emails': [('info@cellsignal.com', 'Cell Signaling Technology')],
        'descriptions': [
            u'''
            On these resource pages you can find signaling pathway diagrams, research overviews, relevant antibody products, publications, and other research resources organized by topic. The pathway diagrams associated with these topics have been assembled by CST scientists and outside experts to provide succinct and current overviews of selected signaling pathways.
            '''
        ],
        'notes': [
            u'''
            The pathway diagrams are based on good quality, manually curated data, probably from review articles. However, those are available only in graphical (PDF and InDesign) formats. There is no programmatic way to obtain the interactions and references, as it was confirmed by the authors, who I contacted by mail. Wang's HumanSignalingNetwork includes the data from this resource, which probably has been entered manually, but Wang's data doesn't have source annotations, despite it's compiled from multiple sources. The date of the beginning of this project is estimated using the Internet wayback machine.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'omnipath': False,
        'license': {
            'name': 'No license.',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        },
    },
    'DIP': {
        'year': 2016,
        'releases': [
            2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
            2011, 2012, 2013, 2014, 2015, 2016
        ],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['UCLA', 'Eisenberg Group'],
        'label': 'DIP',
        'full_name': 'Database of Interacting Proteins',
        'color': '',
        'urls': {
            'articles': [
                'http://nar.oxfordjournals.org/content/28/1/289.long',
                'http://nar.oxfordjournals.org/content/29/1/239.long',
                'http://nar.oxfordjournals.org/content/30/1/303.long',
                'http://nar.oxfordjournals.org/content/32/suppl_1/D449.full'
            ],
            'webpages': ['http://dip.doe-mbi.ucla.edu/dip/Main.cgi'],
            'omictools':
            ['http://omictools.com/database-of-interacting-proteins-tool']
        },
        'pubmeds': [10592249, 11125102, 11752321, 14681454],
        'annot': ['mechanism', 'experiment'],
        'recommend':
        'one of the earliest interaction databases; remarkable literature curation effort',
        'descriptions': [
            u'''
            In the beginning (near 2000), it was a entirely manually curated database:
            Currently protein–protein interactions are entered into the DIP only following publication in peer-reviewed journals. Entry is done manually by the curator, followed by automated tests that show the proteins and citations exist. Interactions are double-checked by a second curator and flagged accordingly in the database.
            From 2001, it contains high-throughput interactions:
            Because the reliability of experimental evidence varies widely, methods of quality assessment have been developed and utilized to identify the most reliable subset of the interactions. This CORE set can be used as a reference when evaluating the reliability of high-throughput protein-protein interaction data sets, for development of prediction methods, as well as in the studies of the properties of protein interaction networks.
            '''
        ],
        'notes': [
            u'''
            The 'core' dataset contains manually curated interactions from small-scale studies. Interactions are well annotated with PubMed IDs, evidences, and mechanism (binding, chemical reaction, etc). The format is esily accessible (MITAB).
            '''
        ],
        'type': 'literature curated',
        'subtype': 'interaction',
        'data_integration': 'static',
        'omnipath': True,
        'emails': [('david@mbi.ucla.edu', 'David Eisenberg')],
        'license': {
            'name': 'CC-Attribution-NoDerivs-3.0',
            'url': 'http://creativecommons.org/licenses/by-nd/3.0/',
            'commercial_use': True
        },
        'pypath': {
            'format': [
                'pypath.data_formats.interaction[\'dip\']',
                'pypath.data_formats.interaction[\'omnipath\']'
            ],
            'data': ['pypath.data/dip_human_core_processed.csv']
        },
        'pathguide': 3
    },
    'DEPOD': {
        'year': 2016,
        'releases': [2013, 2014, 2016],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['EMBL & EMBL-EBI'],
        'label': 'DEPOD',
        'full_name': 'Human Dephosphorylation Database',
        'urls': {
            'articles': [
                'http://stke.sciencemag.org/content/6/275/rs10.long',
                'http://nar.oxfordjournals.org/content/43/D1/D531.long'
            ],
            'webpages': ['http://www.koehn.embl.de/depod/index.php'],
            'omictools':
            ['http://omictools.com/dephosphorylation-database-tool']
        },
        'pubmeds': [23674824, 25332398],
        'taxons': ['human'],
        'annot': ['residue', 'experiment'],
        'recommend':
        'the only resource focusing on dephosphorylation interactions; most of these can not be found anywhere else',
        'descriptions': [
            u'''
            DEPOD the human DEPhOsphorylation Database (version 1.0) is a manually curated database collecting human active phosphatases, their experimentally verified protein and non-protein substrates and dephosphorylation site information, and pathways in which they are involved. It also provides links to popular kinase databases and protein-protein interaction databases for these phosphatases and substrates. DEPOD aims to be a valuable resource for studying human phosphatases and their substrate specificities and molecular mechanisms; phosphatase-targeted drug discovery and development; connecting phosphatases with kinases through their common substrates; completing the human phosphorylation/dephosphorylation network.
            '''
        ],
        'notes': [
            u'''
            Nice manually curated dataset with PubMed references, in easily accessible MITAB format with UniProt IDs, comprises 832 dephosphorylation reactions on protein substrates, and few hundreds on small molecules.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'post-translational modification',
        'omnipath': True,
        'emails': [('koehn@embl.de', 'Maja Kohn')],
        'license': {
            'name': 'CC-Attribution-NonCommercial 4.0',
            'url': 'https://creativecommons.org/licenses/by-nc/4.0/',
            'commercial_use': False
        },
        'pypath': {
            'format': [
                'pypath.data_formats.ptm[\'depod\']',
                'pypath.data_formats.omnipath[\'depod\']'
            ],
            'data': [
                'pypath.data/depod-refs.csv',
                'pypath.urls.urls[\'depod\'][\'url\']'
            ],
            'input': ['pypath.dataio.get_depod()'],
            'ptm': ['pypath.pypath.PyPath().load_depod_dmi()']
        }
    },
    'PhosphoPoint': {
        'year': 2008,
        'urls': {
            'articles': [
                'http://bioinformatics.oxfordjournals.org/content/24/16/i14.long'
            ],
            'webpages': [
                'http://kinase.bioinformatics.tw/',
                'http://web.archive.org/web/20140530070613/http://kinase.bioinformatics.tw/'
            ],
            'omictools': ['http://omictools.com/phosphopoint-tool']
        },
        'taxons': ['human'],
        'pubmeds': [18689816],
        'descriptions': [
            u'''
            We have integrated three existing databases, including Phospho.ELM (release 6.0, total 9236 phosphorylation sites), HPRD (release 6, total 8992 phosphorylation sites), SwissProt (release 51.5, total 6529 phosphorylation sites), and our manually curated 400 kinase–substrate pairs, which are primarily from review articles.
            Among these phosphorylation sites, 7843 (6152+995+696) are from high-throughput (HTP) screening, 6329 (3828+1152+1349) are from low-throughput (LTP) analysis, and only 679 (420+97+162) are both from HTP and LTP screening. One special note is that there are 887 phosphorylation sites, which do not have annotation from literature in the SwissProt database and it is not possible distinguish whether these are from HTP or LTP.
            '''
        ],
        'notes': [
            u'''
            It contains 400 manually curated interactions and much more from HTP methods. The manually curated set can not be distinguished in the data formats offered.
            '''
        ],
        'type': 'literature curated and prediction',
        'subtype': 'post-translational modification',
        'data_integration': 'static',
        'omnipath': False,
        'emails': [('cyhuang5@ym.edu.tw', 'Chi-Ying F. Huang'),
                   ('kmchao@csie.ntu.edu.tw', 'Kun-Mao Chao')],
        'license': {
            'name': 'CC-Attribution-ShareAlike 4.0',
            'url': 'https://creativecommons.org/licenses/by-sa/4.0/',
            'commercial_use': True
        },
        'pypath': {
            'format': ['pypath.data_formats.ptm_misc[\'ppoint\']'],
            'data': ['pypath.data/phosphopoint.csv'],
            'misc': ['pypath.pypath.PyPath().phosphopoint_directions()'],
            'input': ['pypath.dataio.phosphopoint_directions()']
        },
        'pathguide': 252
    },
    'PANTHER': {
        'year': 2016,
        'releases':
        [2000, 2001, 2002, 2003, 2005, 2006, 2010, 2011, 2012, 2014, 2016],
        'urls': {
            'articles': [
                'http://link.springer.com/protocol/10.1007%2F978-1-60761-175-2_7#section=82252&page=1',
                'http://nar.oxfordjournals.org/content/35/suppl_1/D247.long'
            ],
            'webpages': ['http://www.pantherdb.org/'],
            'omictools': [
                'http://omictools.com/protein-analysis-through-evolutionary-relationships-tool'
            ]
        },
        'full_name': 'Pathway Analysis Through Evolutionary Relationships',
        'pubmeds': [17130144, 19597783],
        'recommend': '',
        'descriptions': [
            u'''
            References are captured at three levels. First, each pathway as a whole requires a reference. For signaling pathways, at least three references, usually review papers, are required in order to provide a more objective view of the scope of the pathway. For metabolic pathways, a textbook reference is usually sufficient. Second, references are often associated to each molecule class in the pathway. Most of these references are OMIM records or review papers. Third, references are provided to support association of specific protein sequences with a particular molecule class, e.g., the SWISS-PROT sequence P53_HUMAN annotated as an instance of the molecule class ‘‘P53’’ appearing in the pathway class ‘‘P53 pathway’’. These are usually research papers that report the experimental evidence that a particular protein or gene participates in the reactions represented in the pathway diagram.
            There are three major properties that make this infrastructure differ from other pathway curation systems, such as from Reactome and EcoCyc. First, the pathway diagrams are drawn with CellDesigner software. There are two advantages to using CellDesigner. First, controlled graphical notations are used to draw the pathway diagram, and the software automatically creates a computational representation that is compatible with the SBML standard. Second, a pathway diagram can be viewed with an exact, one-to-one correspondence with the ontological representation of the pathways stored in the back-end. The second property is that the scope of the pathway is defined first based on literature, and pathway components (proteins, genes, RNAs) are treated as ontology terms, or molecule classes, rather than specific instances. This means that multiple proteins from the same organism or different organisms can potentially play the same given role in a pathway. The advantage is that the work flow is more similar to the thinking process of the biologists who are the users of our curation software module. The third major property is that the curation software is designed to be simple enough to be used directly by bench biologists after a brief training course. All other pathway databases we are aware of employ highly trained curators, who of course cannot be experts in all areas of biology. The current set of PANTHER pathways has been curated by more than 40 different external experts from the scientific community; they must only have demonstrated their expertise with publications in the relevant field.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'process description',
        'omnipath': False,
        'emails': [('feedback@pantherdb.org', 'Panther Team'),
                   ('paul.thomas@sri.com', 'Paul Thomas')],
        'pathguide': 164,
        'license': {
            'name': 'GNU-GPL v2',
            'url': 'https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html',
            'commercial_use': True
        },
        'pathguide': 164
    },
    'PhosphoSite': {
        'year': 2016,
        'releases': [2011, 2015, 2016],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['CST'],
        'label': 'PhosphoSite',
        'full_name': 'PhosphoSitePlus',
        'color': '',
        'urls': {
            'articles': [
                'http://onlinelibrary.wiley.com/doi/10.1002/pmic.200300772/abstract',
                'http://nar.oxfordjournals.org/content/40/D1/D261.long',
                'http://nar.oxfordjournals.org/content/43/D1/D512.long'
            ],
            'webpages': ['http://www.phosphosite.org/homeAction.do'],
            'omictools': ['http://omictools.com/phosphositeplus-tool']
        },
        'pubmeds': [15174125, 22135298, 25514926],
        'taxons': ['human', 'mouse', 'eubacteria', 'eukarya'],
        'annot': ['experiment', 'mechanism', 'resudue'],
        'recommend':
        'the largest kinase-substrate database; huge amount of literature referenced',
        'descriptions': [
            u'''
            PSP integrates both low- and high-throughput (LTP and HTP) data sources into a single reliable and comprehensive resource. Nearly 10,000 journal articles , including both LTP and HTP reports, have been manually curated by expert scientists from over 480 different journals since 2001.
            Information from nearly 13 000 papers and 600 different journals characterizing modification sites with LTP methods has been curated into PSP.
            Information is gathered from published literature and other sources. Published literature is searched semi-automatically with multiple intelligent search algorithms to identify reports that potentially identify phosphorylation sites in human, mouse or other species. Each identified report is then scanned by our highly trained curatorial staff (all with PhDs and extensive research experience in cell biology or related disciplines) to select only those papers that either identify new physiological phosphorylation sites or those that illuminate the biological function of the phosphorylation event. Records that are selected for inclusion into PhosphoSite are placed in the curatorial queue for processing. Note: while we gather records that describe both in vitro and in vivo phosphorylation events, we only finally submit records about in vitro sites when we have additional hard evidence that the site is also phosphorylated in vivo.
            '''
        ],
        'type': 'literature curated and high throughput',
        'subtype': 'post-translational modification',
        'data_integration': 'dynamic',
        'omnipath': True,
        'emails': [('phornbeck@cellsignal.com', 'Peter V Hornbeck'),
                   ('EditorPhosphoSite@cellsignal.com', 'PhosphoSite Team')],
        'license': {
            'name': 'CC-NonCommercial-ShareAlike',
            'url': 'http://creativecommons.org/licenses/by-nc-sa/3.0/',
            'commercial_use': False
        },
        'pypath': {
            'format': [
                'pypath.data_formats.ptm[\'psite\']',
                'pypath.data_formats.omnipath[\'psite\']',
                'pypath.data_formats.ptm_misc[\'psite_noref\']'
            ],
            'data': [
                'pypath.urls.urls[\'psite_bp\'][\'url\']',
                'pypath.urls.urls[\'psite_kin\'][\'url\']',
                'pypath.urls.urls[\'psite_p\'][\'url\']',
                'pypath.urls.urls[\'psite_reg\'][\'url\']',
                'pypath.data_formats.files[\'phosphosite\'][\'curated\']',
                'pypath.data_formats.files[\'phosphosite\'][\'noref\']'
            ],
            'input': [
                'pypath.dataio.get_phosphosite()',
                'pypath.dataio.phosphosite_directions()',
                'pypath.dataio.get_psite_phos()',
                'pypath.dataio.get_psite_p()', 'pypath.dataio.get_psite_reg()'
            ],
            'intr': [
                'pypath.dataio.get_phosphosite_curated()',
                'pypath.dataio.get_phosphosite_noref()'
            ],
            'misc': ['pypath.pypath.PyPath().phosphosite_directions()'],
            'ptm': ['pypath.pypath.PyPath().load_psite_phos()']
        },
        'pathguide': 82
    },
    'SPIKE': {
        'year': 2012,
        'releases': [2008, 2011, 2012],
        'urls': {
            'articles': [
                'http://www.biomedcentral.com/1471-2105/9/110',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D793.full.html'
            ],
            'webpages': ['http://www.cs.tau.ac.il/~spike/'],
            'omictools': [
                'http://omictools.com/signaling-pathways-integrated-knowledge-engine-tool'
            ]
        },
        'authors': ['Shamir Group', 'Shiloh Group'],
        'pubmeds': [18289391, 21097778],
        'full_name': 'Signaling Pathway Integrated Knowledge Engine',
        'annot': ['experiment', 'mechanism'],
        'recommend':
        'one of the 3 large resources providing effect signs; combined with Signor and SignaLink gives the major part of the literature curated activity flow network',
        'descriptions': [
            u'''
            SPIKE’s data on relationships between entities come from three sources: (i) Highly curated data submitted directly to SPIKE database by SPIKE curators and experts in various biomedical domains. (ii) Data imported from external signaling pathway databaes. At present, SPIKE database imports such data from Reactome, KEGG, NetPath and The Transcription Factor Encyclopedia (http://www.cisreg.ca/cgi-bin/tfe/home.pl). (iii) Data on protein–protein interactions (PPIs) imported either directly from wide-scale studies that recorded such interactions [to date,PPI data were imported from Stelzl et al., Rual et al. and Lim et al.] or from external PPI databases [IntAct and MINT]. Relationship data coming from these different sources vary greatly in their quality and this is reflected by a quality level attribute, which is attached to each relationship in SPIKE database (Supplementary Data). Each relationship in SPIKE is linked to at least one PubMed reference that supports it.
            As of August 2010, the SPIKE database contains 20 412 genes/proteins, 542 complexes (327 of high quality), 320 protein families (167 of high quality) and 39 small molecules. These entities are linked by 34 338 interactions (of which 2400 are of high quality) and 6074 regulations (4420 of high quality). These are associated with 5873 journal references in total.
            Each of the maps is constructed by a domain expert; typically the same expert will also be responsible later for keeping it up-to-date. The expert reads the relevant literature and identifies those interactions and regulations that are pertinent to the pathway.
            The regulations and interactions in the database are assigned quality values ranging from 1 to 4. In general, relationships (regulations and interactions) derived from highly focused biochemical studies are assigned high quality (2 or 1) while those derived from high-throughput experiments are assigned lower quality (4 or 3). The curator uses best judgment to assign a quality level. For example, relationships mentioned in two independent research reports, or cited repeatedly in reviews written by leading authorities will get quality 1. Relationships with cited concrete references and those imported en masse from external curated signaling DBs are initially assigned quality 2 but later can be changed to the highest quality after the curator has read and was convinced by the cited papers. Data imported from protein-protein interaction DBs and datasets are assigned quality 3 or 4, depending on the experimental technique.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'data_integration': 'static',
        'omnipath': True,
        'emails': [('rshamir@tau.ac.il', 'Ron Shamir')],
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        },
        'pypath': {
            'format': [
                'pypath.data_formats.pathway[\'spike\']',
                'pypath.data_formats.omnipath[\'spike\']'
            ],
            'data': ['pypath.data/spike_hc.csv']
        },
        'pathguide': 317
    },
    'NCI-PID': {
        'year': 2012,
        'releases': [2008, 2012],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['NCI'],
        'label': 'NCI-PID',
        'full_name': 'NCI-Nature Pathway Interaction Database',
        'color': '',
        'urls': {
            'webpages': [
                'http://pid.nci.nih.gov/index.shtml',
                'http://www.ndexbio.org/',
                'http://web.archive.org/web/20150908030438/http://pid.nci.nih.gov/index.shtml'
            ],
            'articles':
            ['http://nar.oxfordjournals.org/content/37/suppl_1/D674.long']
        },
        'pubmeds': [18832364],
        'data_import': ['BioCarta', 'Reactome'],
        'contains': ['BioCarta', 'Reactome'],
        'taxons': ['human'],
        'annot': ['experiment', 'localization', 'mechanism', 'pathway'],
        'recommend':
        'one of the largest process description resources; represents remarkable curation effort; discontinued',
        'descriptions': [
            u'''
            In curating, editors synthesize meaningful networks of events into defined pathways and adhere to the PID data model for consistency in data representation: molecules and biological processes are annotated with standardized names and unambiguous identifiers; and signaling and regulatory events are annotated with evidence codes and references. To ensure accurate data representation, editors assemble pathways from data that is principally derived from primary research publications. The majority of data in PID is human; however, if a finding discovered in another mammal is also deemed to occur in humans, editors may decide to include this finding, but will also record that the evidence was inferred from another species. Prior to publication, all pathways are reviewed by one or more experts in a field for accuracy and completeness.
            '''
        ],
        'notes': [
            u'''
            From the NCI-XML interactions with references, directions and signs can be extracted. Complexes are ommited.
            ''', u'''
            From the end of 2015, the original NCI-PID webpage is not accessible anymore, and the data is available through the NDEx webserver and API.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'process description',
        'omnipath': False,
        'emails': [('yanch@mail.nih.gov', 'Chunhua Yan')],
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False
        },
        'pypath': {
            'data': ['pypath.urls.urls[\'nci-pid\'][\'biopax_l3\']'],
            'format': ['pypath.data_formats.reaction_misc[\'nci_pid\']'],
            'input': ['pypath.dataio.pid_biopax()'],
            'intr': ['pypath.dataio.pid_interactions()']
        },
        'pathguide': 119
    },
    'WikiPathways': {
        'year': 2016,
        'releases': [2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016],
        'urls': {
            'webpages':
            ['http://www.wikipathways.org/index.php/WikiPathways'],
            'articles': ['http://nar.oxfordjournals.org/content/40/D1/D1301'],
            'omictools': ['http://omictools.com/wikipathways-tool']
        },
        'recommend':
        'process description resource developed in a community effort; high proportion of specific proteins',
        'descriptions': [
            u'''
            The goal of WikiPathways is to capture knowledge about biological pathways (the elements, their interactions and layout) in a form that is both human readable and amenable to computational analysis.
            '''
        ],
        'notes': [
            u'''
            The data is not accessible. Interactions are available in BioPAX format, but without references.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'process description',
        'omnipath': False,
        'emails': [('thomaskelder@gmail.com', 'Thomas Kelder'),
                   ('apico@gladstone.ucsf.edu', 'Alex Pico')],
        'license': {
            'name': 'CC-Attribution-3.0',
            'url': 'http://creativecommons.org/licenses/by/3.0/',
            'commercial_use': True
        },
        'pathguide': 237
    },
    'ConsensusPathDB': {
        'year': 2015,
        'releases': [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015],
        'urls': {
            'webpages': ['http://cpdb.molgen.mpg.de/CPDB'],
            'articles': [
                'http://nar.oxfordjournals.org/content/37/suppl_1/D623.long',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D712.long',
                'http://nar.oxfordjournals.org/content/41/D1/D793.long'
            ],
            'omictools': ['http://omictools.com/consensuspathdb-tool']
        },
        'taxons': ['human', 'mouse', 'yeast'],
        'pubmeds': [18940869, 23143270, 21071422],
        'descriptions': [
            '''
            Interaction data in ConsensusPathDB currently originates from 12 interaction databases and comprises physical interactions, biochemical reactions and gene regulations. Importantly, the source of physical entities and interactions is always recorded, which allows linking to the original data in the source database.
            In order to assess the content overlap of the source databases and to reduce redundancy, we have applied a method to merge identical physical entities and identify similar interactions. The method is straightforward and efficient for the integration of networks from any single species. Simple physical entities of the same type (genes, proteins, transcripts, metabolites) are compared on the basis of common database identifiers like UniProt, Ensembl, Entrez, ChEBI, etc. Since different databases tend to annotate physical entities with different identifier types (e.g. some databases annotate proteins with UniProt identifiers, others with Ensembl identifiers), we first translated the annotations to a uniform identifier type, which is a UniProt entry name in case of proteins, Ensembl gene ID in case of genes and transcripts, and KEGG/ChEBI ID in case of metabolites. Protein complexes are compared according to their individual protein composition. Simple physical entities with the same identifier, and complexes with the same composition, are merged in ConsensusPathDB. Information provided by the according source databases for the merged entities is stored in a complementary manner.
            Functional interactions of physical entities are also compared with each other. Here, we distinguish between primary and secondary interaction participants. Primary participants are substrates and products in case of biochemical reactions, interactors in case of physical interactions and target genes in case of gene regulation. All other participants, e.g. enzymes and interaction modifiers, are secondary participants. If the primary participants of two or more interactions match, these interactions are considered similar. Two similar interactions may have different stoichiometry, modification and/or localization of the participants. To allow for flexibility, similar interactions are marked as such in the database, but the decision whether they should be considered identical despite mismatching details is left to the user and depends on his specific problem. Moreover, ConsensusPathDB does not provide any additional quality control filters. All interactions provided by the different database sources are treated in the same way.
            '''
        ],
        'notes': [
            u'''
            ConsensusPathDB comprises data from 32 resources. The format is easy to use, tab delimited text file, with UniProtKB names and PubMed IDs. However, the dataset is extremely huge, and several databases containing HTP data is included.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'activity flow',
        'emails': [('kamburov@molgen.mpg.de', 'Atanas Kamburov')],
        'license': {
            'name':
            'Constituting resources carry their own licenses. "Due to several licensing issues, we are not allowed to release the complete integrated network (including signaling, metabolism and gene regulation)."',
            'url': 'http://cpdb.molgen.mpg.de/CPDB/tutorial#moreinfo.lic',
            'commercial_use': False
        },
        'pathguide': 275
    },
    'KEGG': {
        'year': 2016,
        'releases': [
            2000, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015,
            2016
        ],
        'urls': {
            'webpages': ['http://www.genome.jp/kegg/'],
            'articles':
            ['http://nar.oxfordjournals.org/content/28/1/27.long'],
            'omictools': ['http://omictools.com/kegg-tool']
        },
        'full_name': 'Kyoto Encyclopedia of Genes and Genomes',
        'notes': [
            u'''
            From 2011, KEGG data is not freely available. The downloadable KGML files contain binary interactions, most of them between large complexes. No references available.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'process description',
        'omnipath': False,
        'emails': [('kanehisa@kuicr.kyoto-u.ac.jp', 'Minoru Kaneshia')],
        'license': {
            'name': 'KEGG License, non-free',
            'url': 'http://www.genome.jp/kegg/legal.html',
            'commercial_use': False
        },
        'pathguide': 16,
        'pypath': {
            'data': [
                'pypath.urls.urls[\'kegg_pws\'][\'kgml_url\']',
                'pypath.urls.urls[\'kegg_pws\'][\'list_url\']'
            ],
            'input': ['pypath.dataio.kegg_pathways()'],
            'misc': ['pypath.pypath.kegg_directions()']
        }
    },
    'BioGRID': {
        'year': 2016,
        'releases': [
            2003, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015,
            2016
        ],
        'label': 'BioGRID',
        'authors': ['Tyers Lab'],
        'urls': {
            'webpages': ['http://thebiogrid.org/'],
            'articles': [
                'http://genomebiology.biomedcentral.com/articles/10.1186/gb-2003-4-3-r23',
                'http://nar.oxfordjournals.org/content/34/suppl_1/D535.long',
                'http://nar.oxfordjournals.org/content/36/suppl_1/D637.long',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D698.long',
                'http://nar.oxfordjournals.org/content/41/D1/D816.long',
                'http://nar.oxfordjournals.org/content/43/D1/D470.long'
            ],
            'omictools': [
                'http://omictools.com/biological-general-repository-for-interaction-datasets-tool'
            ]
        },
        'pubmeds':
        [25428363, 23203989, 21071413, 18000002, 16381927, 12620108],
        'type': 'high throughput',
        'subtype': 'interaction',
        'recommend':
        'one of the largest interaction resources; contains both HTP and literature curated interactions',
        'full_name': 'Biological General Repository for Interaction Datasets',
        'omnipath': True,
        'emails': [('biogridadmin@gmail.com', 'BioGRID Team'),
                   ('md.tyers@umontreal.ca', 'Michael Tyers')],
        'license': {
            'name': 'MIT License',
            'url': 'https://biogrid-downloads.nyc3.digitaloceanspaces.com/LICENSE.txt',
            'commercial_use': True
        },
        'pypath': {
            'format': [
                'pypath.data_formats.interaction[\'biogrid\']',
                'pypath.data_formats.omnipath[\'biogrid\']'
            ],
            'data': ['pypath.urls.urls[\'biogrid\']'],
            'intr': ['pypath.dataio.biogrid_interactions()']
        },
        'pathguide': 7
    },
    'STRING': {
        'year': 2016,
        'releases': [2016, 2015, 2013, 2011, 2009, 2007, 2005, 2003, 2000],
        'urls': {
            'webpages': ['http://string-db.org/'],
            'articles': [
                'http://nar.oxfordjournals.org/content/43/D1/D447.long',
                'http://nar.oxfordjournals.org/content/41/D1/D808.long',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D561.long',
                'http://nar.oxfordjournals.org/content/37/suppl_1/D412.long',
                'http://nar.oxfordjournals.org/content/35/suppl_1/D358.long',
                'http://nar.oxfordjournals.org/content/33/suppl_1/D433.long',
                'http://nar.oxfordjournals.org/content/31/1/258.long',
                'http://nar.oxfordjournals.org/content/28/18/3442.long'
            ],
            'omictools': ['http://omictools.com/string-tool']
        },
        'pubmeds': [
            25352553, 23203871, 21045058, 18940858, 17098935, 15608232,
            12519996, 10982861
        ],
        'authors': ['Bork Lab'],
        'label': 'STRING',
        'type': 'high-throughput and prediction',
        'subtype': 'interaction',
        'omnipath': False,
        'emails': [('bork@embl.de', 'Peer Bork'),
                   ('lars.juhl.jensen@cpr.ku.dk', 'Lars Juhl Jensen'),
                   ('mering@imls.uzh.ch', 'Christian von Mering')],
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'http://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        },
        'pathguide': 93
    },
    'MINT': {
        'year': 2015,
        'releases': [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015],
        'label': 'MINT',
        'urls': {
            'webpages': ['http://mint.bio.uniroma2.it/mint/Welcome.do'],
            'articles':
            ['http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1751541/']
        },
        'type': 'literature curated and high-throughput',
        'subtype': 'interaction',
        'full_name': 'Molecular Interaction Database',
        'omnipath': False,
        'emails': [('livia.perfetto@live.it', 'Livia Perfetto')],
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'http://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        },
        'pathguide': 17
    },
    'IntAct': {
        'year': 2016,
        'releases':
        [2003, 2006, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['EBI'],
        'label': 'IntAct',
        'full_name': 'IntAct Molecular Interaction Database',
        'data_import': ['InnateDB', 'MINT'],
        'urls': {
            'articles': [
                'http://nar.oxfordjournals.org/content/42/D1/D358.long',
                'http://nar.oxfordjournals.org/content/40/D1/D841.long',
                'http://nar.oxfordjournals.org/content/38/suppl_1/D525.long',
                'http://nar.oxfordjournals.org/content/35/suppl_1/D561.long',
                'http://nar.oxfordjournals.org/content/32/suppl_1/D452.long'
            ],
            'webpages': ['http://www.ebi.ac.uk/intact/'],
            'omictools': ['http://omictools.com/intact-tool']
        },
        'pubmeds': [14681455, 17145710, 19850723, 22121220, 24234451],
        'annot': ['experiment', 'mechanism'],
        'recommend':
        'the largest interaction resource; highest coverage on the human proteome; miscore confidence score',
        'descriptions': [
            u'''
            The information within the IntAct database primarily consists of protein–protein interaction (PPI) data. The majority of the PPI data within the database is annotated to IMEx standards, as agreed by the IMEx consortium. All such records contain a full description of the experimental conditions in which the interaction was observed. This includes full details of the constructs used in each experiment, such as the presence and position of tags, the minimal binding region defined by deletion mutants and the effect of any point mutations, referenced to UniProtKB, the underlying protein sequence database. Protein interactions can be described down to the isoform level, or indeed to the post-translationally cleaved mature peptide level if such information is available in the publication, using the appropriate UniProtKB identifiers.
            Each entry in IntAct is peer reviewed by a senior curator, and not released until accepted by that curator. Additional rule-based checks are run at the database level, and manually fixed when necessary. Finally, on release of the data, the original author of each publication is contacted and asked to comment on the representation of their data; again manual updates are made to the entry should the author highlight any errors.
            All binary interactions evidences in the IntAct database, including those generated by Spoke expansion of co-complex data, are clustered to produce a non-redundant set of protein pairs (R. C. Jimenez et al., manuscript in preparation). Each binary pair is then scored, using a simple addition of the cumulated value of a weighted score for the interaction detection method and the interaction type for each interaction evidence associated with that binary pair, as described using the PSI-MI CV terms. The scores are given in Table 1, all children of each given parent receives that score. Only experimental data is scored, inferred interactions, for example, would be excluded. Any low confidence data or data manually tagged by a curator for exclusion from the process, would not be scored. Isoforms and post-processed protein chains are regarded as distinct proteins for scoring purposes.
            '''
        ],
        'notes': [
            u'''
            We can not draw a sharp distinction between low and high throughput methods, and I can agree, that this is not the only and best measure of quality considering experimental data. I see that IntAct came up with a good solution to estimate the confidence of interactions. The mi-score system gives a comprehensive way to synthetize information from multiple experiments, and weight interactions according to experimental methods, interaction type, and number of evidences.
            '''
        ],
        'type': 'literature curated and high-throughput',
        'subtype': 'interaction',
        'omnipath': True,
        'emails': [('orchard@ebi.ac.uk', 'Sandra Orchard'),
                   ('hhe@ebi.ac.uk', 'Henning Hermjakob')],
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        },
        'pathguide': 111,
        'pypath': {
            'data': ['pypath.data/intact_filtered.csv'],
            'format': ['pypath.data_formats.interaction_misc[\'intact\']']
        }
    },
    'MatrixDB': {
        'year': 2015,
        'releases': [2009, 2011, 2015],
        'urls': {
            'articles': [
                'http://bioinformatics.oxfordjournals.org/content/25/5/690.long',
                'http://nar.oxfordjournals.org/content/43/D1/D321.long',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D235.long'
            ],
            'webpages': ['http://matrixdb.univ-lyon1.fr/'],
            'omictools': ['http://omictools.com/matrixdb-tool']
        },
        'pubmeds': [19147664, 20852260, 25378329],
        'taxons': ['mammalia'],
        'annot': ['experiment'],
        'recommend':
        'small, literature curated interaction resource; many interactions for receptors and extracellular proteins',
        'descriptions': [
            u'''
            Protein data were imported from the UniProtKB/Swiss-Prot database (Bairoch et al., 2005) and identified by UniProtKB/SwissProt accession numbers. In order to list all the partners of a protein, interactions are associated by default to the accession number of the human protein. The actual source species used in experiments is indicated in the page reporting interaction data. Intracellular and membrane proteins were included to obtain a comprehensive network of the partners of extracellular molecules. Indeed, ECM proteins and GAGs bind to a number of membrane proteins or cell-associated proteoglycans and some of them interact with intracellular partners upon internalization (Dixelius et al., 2000). ECM proteins were identified by the UniProtKB/Swiss-Prot keyword ‘extracellular matrix’ and by the GO terms ‘extracellular matrix’, ‘proteinaceous extracellular matrix’ and their child terms. The proteins annotated with the GO terms ‘extracellular region’ and ‘extracellular space’, which are used for proteins found in biological fluids, were not included because circulating molecules do not directly contribute to the extracellular scaffold. Additionally, 96 proteins were manually (re-)annotated through literature curation. MatrixDB integrates 1378 interactions from the Human Protein Reference Database (HPRD, Prasad et al., 2009), 211 interactions from the Molecular INTeraction database (MINT, Chatr-Aryamontri et al., 2007), 46 interactions from the Database of Interacting Proteins (DIP, Salwinski et al., 2004), 232 interactions from IntAct (Kerrien et al., 2007a) and 839 from BioGRID (Breitkreutz et al., 2008) involving at least one extracellular biomolecule of mammalian origin. We added 283 interactions from manual literature curation and 65 interactions from protein and GAG array experiments.
            ''', u'''
            Interaction data stored in MatrixDB are (i) experimentally determined in the laboratory using surface plasmon resonance (SPR) binding assays, including protein and glycosaminoglycan arrays probed by SPR imaging, (ii) extracted from the literature by manual curation and (iii) imported from other interaction databases belonging to the IMEx consortium [IntAct, DIP, MINT, BioGRID], as well as from the Human Protein Reference Database. Imported data are restricted to interactions involving at least one extracellular protein.
            ''', u'''
            The content of MatrixDB has been updated with new interaction data manually curated by the MatrixDB team, and by importing interaction data from four interaction databases of the IMEx consortium via The Proteomics Standard Initiative Common QUery InterfaCe (PSICQUIC), a community standard for computational access to molecular-interaction data resources.  In the current release MatrixDB contains 904 interactions supported by 1244 experiments, which have been manually curated from 237 publications, compared to 490 interactions supported by 847 experiments in the previous version of the database. This is the MatrixDB ‘core’ data set.
            '''
        ],
        'notes': [
            u'''
            Very nice! Note: The interactions imported from IMEX databases or any other database, are collected separately, in the PSICQUIC-extended dataset. The MatrixDB-core dataset is curated manually by the MatrixDB team.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'interaction',
        'data_integration': 'static',
        'omnipath': True,
        'emails': [('matrixdb@ibcp.fr', 'MatrixDB Team'),
                   ('sylvie.ricard-blum@ibcp.fr', 'Sylvie Ricard-Blum')],
        'license': {
            'name': 'CC BY 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        },
        'pypath': {
            'format': [
                'pypath.data_formats.interaction[\'matrixdb\']',
                'pypath.data_formats.omnipath[\'matrixdb\']'
            ],
            'data': ['pypath.data/matrixdb_core.csv']
        },
        'pathguide': 298
    },
    'PathwayCommons': {
        'year': 2016,
        'releases': [2010, 2011, 2012, 2013, 2014, 2015, 2016],
        'urls': {
            'webpages': [
                'http://www.pathwaycommons.org/pc2/',
                'http://www.pathwaycommons.org/about/'
            ],
            'articles':
            ['http://nar.oxfordjournals.org/content/39/suppl_1/D685.long'],
            'omictools': ['http://omictools.com/pathway-commons-tool']
        },
        'label': 'PathwayCommons',
        'authors': ['Bader Lab', 'MSKCC cBio'],
        'notes': [
            u'''
            Pathway Commons is a collection of publicly available pathway information from multiple organisms. It provides researchers with convenient access to a comprehensive collection of biological pathways from multiple sources represented in a common language for gene and metabolic pathway analysis.
            ''', u'''
            Pathway Commons integrates a number of pathway and molecular interaction databases supporting BioPAX and PSI-MI formats into one large BioPAX model, which can be queried using our web API (documented below). This API can be used by computational biologists to download custom subsets of Pathway Commons for analysis, or can be used to incorporate powerful biological pathway and network information retrieval and query functionality into websites, scripts and software. For computational biologists looking for comprehensive biological pathway data for analysis, we also make available batch downloads of the data in several formats.
            ''', u'''
            Warehouse data (canonical molecules, ontologies) are converted to BioPAX utility classes, such as EntityReference, ControlledVocabulary, EntityFeature sub-classes, and saved as the initial BioPAX model, which forms the foundation for integrating pathway data and for id-mapping.
            Pathway and binary interaction data (interactions, participants) are normalized next and merged into the database. Original reference molecules are replaced with the corresponding BioPAX warehouse objects.
            '''
        ],
        'data_import': [
            'Reactome',
            'NCI-PID',
            'CancerCellMap',
            'BioCarta',
            'HPRD',
            'PhosphoSite',
            'PANTHER',
            'DIP',
            'IntAct',
            'BioGRID',
            'BIND',
            'CORUM',
        ],
        'license': {
            'name': 'Constituting databases carry their own licenses.',
            'url': 'https://www.pathwaycommons.org/pc2/datasources',
            'commercial_use': False # Partially
        },
        'type': 'combined',
        'subtype': 'interaction',
        'pubmeds': [21071392],
        'omnipath': False,
        'emails': [('gary.bader@utoronto.ca', 'Gary Bader')],
        'pathguide': 293
    },
    'Laudanna': {
        'year': 2014,
        'releases': [2014],
        'urls': {
            'webpages': ['http://dp.univr.it/~laudanna/LCTST/downloads/']
        },
        'full_name':
        'Compiled Datasets for Network Analysis from Laudanna Lab',
        'authors': ['Laudanna Lab'],
        'type': 'combined',
        'subtype': 'mixed',
        'notes': [
            '''
            Data sets are compiled from public data-bases and from literature and manually curated for accuracy. They are intended for network reconstruction, topological and multidimensional analysis in cell biology.
            '''
        ],
        'omnipath': False,
        'data_import': [
            'BioGRID', 'ConsensusPathDB', 'dbPTM', 'DIP',
            'HumanSignalingNetwork', 'IntAct', 'MINT', 'MPPI',
            'PathwayCommons', 'phospho.ELM', 'PhosphoPoint', 'PhosphoSite',
            'SignaLink'
        ],
        'emails': [('giovanni.scardoni@gmail.com', 'Giovanni Scardoni')],
        'license': {
            'name': 'CC-Attribution-NonCommercial 4.0',
            'url': 'https://creativecommons.org/licenses/by-nc/4.0/',
            'commercial_use': False
        },
        'pypath': {
            'data': [
                'pypath.urls.urls[\'sigflow\']',
                'pypath.urls.urls[\'sigdir\']'
            ],
            'input': [
                'pypath.dataio.get_laudanna_directions()',
                'pypath.dataio.get_laudanna_effects()'
            ],
            'misc': [
                'pypath.pypath.PyPath().laudanna_directions()',
                'pypath.pypath.PyPath().laudanna_effects()'
            ]
        }
    },
    'ORegAnno': {
        'year': 2016,
        'urls': {
            'articles': [
                'https://academic.oup.com/nar/article/44/D1/D126/2502683',
            ],
            'webpages': [
                'http://www.oreganno.org/',
            ]
        },
        'emails': [('ogriffit@genome.wustl.edu', 'Obi L. Griffith')],
        'type': 'literature curated & high throughput',
        'subtype': 'transcription regulation',
        'full_name': 'Open Regulatory Annotation',
        'omnipath': False,
        'dorothea': True,
        'pubmeds': [26578589],
        'license': {
            'name': 'GNU LGPLv3',
            'url': 'http://www.gnu.org/licenses/license-list.html#LGPLv3',
            'commercial_use': True
        },
        'notes': [
            '''
            One of the largest TF-target databases. Covers at least 18
            organisms and contains data from literature curation and many
            screening technologies and in silico prediction.
            '''
        ],
        'pypath': {
            'data': [
                'pypath.urls.urls[\'oreganno\']',
            ],
            'input': [
                'pypath.dataio.get_oreganno()',
            ],
        },
        'omictools': ['https://omictools.com/oreganno-tool'],
    },
    'PAZAR': {
        'year': 2009,
        'releases': [2007],
        'urls': {
            'articles': [
                (
                    'https://genomebiology.biomedcentral.com/'
                    'articles/10.1186/gb-2007-8-10-r207'
                ),
                'https://academic.oup.com/nar/article/37/suppl_1/D54/1010368',
            ],
            'webpages': [
                'http://www.pazar.info/',
            ],
        },
        'emails': [('wyeth@cmmt.ubc.ca', 'Wyeth Wasserman')],
        'authors': ['Wasserman Lab'],
        'type': 'literature curated & high throughput',
        'subtype': 'transcription regulation',
        'full_name': (
            'A Public Database of Transcription Factor '
            'and Regulatory Sequence Annotation'
        ),
        'omnipath': False,
        'dorothea': True,
        'pubmeds': [18971253, 17916232],
        'license': {
            'name': 'GNU LGPLv3',
            'url': 'http://www.gnu.org/licenses/license-list.html#LGPLv3',
            'commercial_use': True
        },
        'notes': [
            '''
            One of the oldest and largest TF-target databases. From the
            Wasserman Lab, who also developed JASPAR and many other tools.
            Unfortunately the website is down at the moment (April 2019).
            Which was, by the way, in that time (2007) a super nice and
            innovative design for a molecular database webpage.
            '''
        ],
        'pypath': {
            'data': [
                'pypath.urls.urls[\'pazar\']',
            ],
            'input': [
                'pypath.dataio.get_pazar()',
            ],
        },
        'omictools': ['http://omictools.com/pazar-tool'],
    },
# From here on, only license info
    'TRRUST': {
        'license': {
            'name': 'CC-Attribution-ShareAlike 4.0',
            'url': 'http://creativecommons.org/licenses/by-nc-sa/4.0/',
            'commercial_use': True
        }
    },
    'HTRIdb': {
        'license': {
            'name': 'GNU LGPLv3',
            'url': 'http://www.gnu.org/licenses/license-list.html#LGPLv3',
            'commercial_use': True
        }
    },
    'NFIRegulomeDB': {
        'license': {
            'name': 'GNU LGPLv3',
            'url': 'http://www.gnu.org/licenses/license-list.html#LGPLv3',
            'commercial_use': True
        }
    },
    'TFe': {
        'license': {
            'name': 'CC-Attribution-ShareAlike 3.0',
            'url': 'https://creativecommons.org/licenses/by-sa/3.0/',
            'commercial_use': True
        }
    },
    'TRRD': {
        'license': {# Via TfactS
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True # Again, via TfactS
        }
    },
    'TRED': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'Fantom4': {
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'http://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        }
    },
    'RegNetwork': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'TfactS': {
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        }
    },
    'GeneOntology': {
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        }
    },
    'Membranome': {
        'license': {
            'name': 'Apache 2.0',
            'url': 'http://www.apache.org/licenses/LICENSE-2.0.html',
            'commercial_use': True
        }
    },
    'Exocarta': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'Vesiclepedia': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'Matrisome': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'CSPA': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'HPMR': {
        'license': {
            'name': 'CC-Attribution-NonCommercial 4.0',
            'url': 'https://creativecommons.org/licenses/by-nc/4.0/',
            'commercial_use': False
        }
    },
    'LOCATE': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'CellPhoneDB': {
        'license': {
            'name': 'MIT License',
            'url': 'https://github.com/Teichlab/cellphonedb/blob/master/LICENSE',
            'commercial_use': True
        }
    },
    'ComPPI': {
        'license': {
            'name': 'CC-Attribution-ShareAlike 4.0',
            'url': 'https://creativecommons.org/licenses/by-sa/4.0/',
            'commercial_use': True
        }
    },
    'Kirouac2010': {
        'license': {
            'name': 'CC-Attribution-NonCommercial-NoDerivs 3.0',
            'url': 'https://creativecommons.org/licenses/by-nc-nd/3.0/',
            'commercial_use': False
        }
    },
    'Ramilowski2015': {
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'http://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        }
    },
    'Adhesome': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'Integrins': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'Zhong2015': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'HGNC': {
        'license': {
            'name': 'Custom license',
            'url': 'https://www.genenames.org/about/',
            'commercial_use': True
        }
    },
    'TopDB': {
        'license': {
            'name': 'CC-Attribution-NonCommercial 4.0',
            'url': 'https://creativecommons.org/licenses/by-nc/4.0/',
            'commercial_use': False
        }
    },
    'OPM': {
        'license': {
            'name': 'Apache 2.0',
            'url': 'http://www.apache.org/licenses/LICENSE-2.0.html',
            'commercial_use': True
        }
    },
    'Havugimana': {
        'license': {
            'name': 'CC-Attribution-ShareAlike 4.0',
            'url': 'https://creativecommons.org/licenses/by-sa/4.0/',
            'commercial_use': True
        }
    },
    'Compleat': {
        'license': {
            'name': 'CC-Attribution-NonCommercial 4.0',
            'url': 'https://creativecommons.org/licenses/by-nc/4.0/',
            'commercial_use': False
        }
    },
    'ComplexPortal': {
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        }
    },
    'HPA': {
        'license': {
            'name': 'CC-Attribution-ShareAlike 3.0',
            'url': 'https://creativecommons.org/licenses/by-sa/3.0/',
            'commercial_use': True
        }
    },
    'PDB': {
        'license': {
            'name': 'Custom license',
            'url': 'https://www.wwpdb.org/about/privacy',
            'commercial_use': True
        }
    },
    'Humap': {
        'license': {
            'name': 'CC0-Attribution',
            'url': 'https://creativecommons.org/share-your-work/public-domain/cc0/',
            'commercial_use': True
        }
    },
    'PhosphoNetworks': {
        'license': {
            'name': 'Custom license',
            'url': 'https://phosphonetworks.org/about.html',
            'commercial_use': False
        }
    },
    'MIMP': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'Li2012': {
        'license': {
            'name': 'CC-Attribution-NonCommercial 3.0',
            'url': 'https://creativecommons.org/licenses/by-nc/3.0/',
            'commercial_use': False
        }
    },
    'UniProt': {
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'http://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        }
    },
    'miRBase': {
        'license': {
            'name': 'CC0',
            'url': 'https://creativecommons.org/share-your-work/public-domain/cc0/',
            'commercial_use': True
        }
    },
    'mir2Disease': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'miRDeatdhDB': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'TransmiR': {
        'license': {
            'name': 'CC-Attribution-NonCommercial 4.0',
            'url': 'https://creativecommons.org/licenses/by-nc/4.0/',
            'commercial_use': False
        }
    },
    'ENCODE': {
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        }
    },
    'lncRNADisease': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'lncrnadb': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'miRecords': {
        'license': {
            'name': 'Custom license',
            'url': 'http://c1.accurascience.com/miRecords/copyright.php',
            'commercial_use': False
        }
    },
    'miRTarBase': {
        'license': {
            'name': 'Custom license',
            'url': 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/LICENSE',
            'commercial_use': False
        }
    },
    'HIPPIE': {
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'https://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True
        }
    },
    'CPAD': {
        'license': {
            'name': 'Custom license',
            'url': 'https://www.iitm.ac.in/bioinfo/CPAD/',
            'commercial_use': False
        }
    },
    'IntOGen': {
        'license': {
            'name': 'CC-Attribution-NonCommercial 4.0',
            'url': 'http://creativecommons.org/licenses/by-nc/4.0/',
            'commercial_use': False
        }
    },
    'COSMIC': {
        'license': {
            'name': 'Custom license',
            'url': 'https://cancer.sanger.ac.uk/cosmic/license',
            'commercial_use': False
        }
    },
    'DGIdb': {
        'license': {
            'name': 'MIT license',
            'url': 'https://github.com/griffithlab/dgi-db/blob/master/LICENSE',
            'commercial_use': True
        }
    },
    'DisGeNet': {
        'license': {
            'name': 'CC-Attribution-ShareAlike-NonCommercial 4.0',
            'url': 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
            'commercial_use': False
        }
    },
    'kinase.com': {
        'license': {
            'name': 'Custom license',
            'url': 'http://kinase.com/about/Disclaimer.html',
            'commercial_use': False
        }
    },
    'phosphatome': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'Vaquerizas2009': {# PENDING
        #'license': {
        #    'name': '',
        #    'url': '',
        #    'commercial_use': False
        #}
    },
    'ProtMapper': {
        'license': {
            'name': 'CC-Attribution-ShareAlike-NonCommercial 4.0',
            'url': 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
            'commercial_use': False
        }
    },
    'ABS': {
        'license': {
            'name': 'GNU-GPLv2',
            'url': 'http://genome.crg.es/main/GNU-GPL.html',
            'commercial_use': True
        }
    },
    'LRdb': {
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False,
        },
        'emails': [('jacques.colinge@inserm.fr', 'Jacques Colinge'),],
    },
    'MSigDB': {
        'license': {
            'name': 'CC-Attribution 4.0',
            'url': 'http://creativecommons.org/licenses/by/4.0/',
            'commercial_use': True,
        },
        'emails': [('gsea-team@broadinstitute.org', 'GSEA Team'),],
    },
    'Baccin2019': {
        'license': {
            'name': 'No license',
            'url': 'http://www.gnu.org/licenses/license-list.html#NoLicense',
            'commercial_use': False,
        },
        'emails': [
            ('lars.velten@embl.de', 'Lars Velten'),
            ('a.trumpp@dkfz-heidelberg.de', 'Andreas Trumpp'),
            ('s.haas@dkfz-heidelberg.de', 'Simon Haas'),
        ],
    },
}

pypath_methods = {
    'data': 'Data source (URLs and files)',
    'format': 'Data format definition',
    'intr': 'Interactions',
    'input': 'Data input methods',
    'ptm': 'Enzyme-substrate relationships and PTMs',
    'dmi': 'Domain-motif interactions',
    'ddi': 'Domain-domain interactions',
    'misc': 'Miscellaneous'
}


def gen_html():
    '''
    Generates a HTML page from the `descriptions` array.
    This HTML is provided by the webservice under `/info`,
    or can be saved locally with `write_html()`.
    '''
    # Header
    title = 'Metadata about signaling pathway resources'
    doc = (
        '<div class="yellowbox box">\n'
        '<p>\n'
            '<em>\n'
            'Information on this page has been last revised in Nov 2016.\n'
            'As of Oct 2019 we are working on updating and extending this\n'
            'page and will publish the new version soon.\n'
            'About updates of the OmniPath database content please refer to\n'
            '<a href="http://archive.omnipathdb.org/README.txt">\n'
                'our archive.\n'
            '</a>\n'
            '</em>\n'
        '</p>\n'
        '</div>\n'
        '<p>This collection was created during the construction '
        'of OmniPath when we considered more than 50 resources and '
        'selected the ones containing literature curation effort. '
        'OmniPath is a network '
        'of signaling pathways intending to '
        'combine all high quality, manually curated efforts. The '
        'descriptions here cite the relevant sentences '
        'about the curation protocols from the original articles and webpages. '
        'URLs pointing to the articles and the webpages, and some '
        'additional metadata are provided where available. '
        'The resources with green title are included by default in '
        'OmniPath. <span class="code">pypath</span> methods are listed '
        ' where available, to know more please look at <a '
        'target="_blank" href="http://pypath.omnipathdb.org/">'
        'pypath documentation.</a> This list is only about network '
        'resources. <span class="code">pypath</span> is able to '
        'process and integrate many other resources, please see '
        'the paper and the documentation to know more.</p>'
        '<p class="small">We searched for license information '
        'in the main, About, Download and FAQ sections of the webpages, '
        'and run Google searches for the database name and license. '
        'Where we could not find anything about licensing, we assumed '
        'no license. Unfortunately due to todays restrictive copyright '
        'legislations, users don\'t have the freedom to use, modify and '
        'redistribute the data without a license explicitely granting '
        'these to them. Despite the clear intention from the authors to '
        'make their data public, and statements on the webpage like '
        '"free to use" or "available for download".</p>\n'
    )
    doc += '\t<h2>Contents</h2>\n'
    doc += '\t<ul>\n'
    # Table of Content
    for k, v in sorted(descriptions.items(), key=lambda x: x[0].lower()):
        doc += '\t\t\t<li><a href="#%s" class="%s">%s</a></li>\n' % \
            (k, 'omnipath' if 'omnipath' in v and v['omnipath'] else 'base',
             v['label'] if 'label' in v else k)
    doc += '\t</ul>\n'
    # Sections
    for k, v in sorted(descriptions.items(), key=lambda x: x[0].lower()):
        doc += u'\t\t<br>\n\t\t<h2 id="%s" class="%s">%s%s</h2>\n' % \
            (k, 'omnipath' if 'omnipath' in v and v['omnipath'] else 'base',
             v['label'] if 'label' in v else k,
             (u' – %s' % (v['full_name'],)) if 'full_name' in v else '')
        doc += '\t\t\t<p><b>Category || Subcategory &gt;&gt;&gt;</b> %s || %s</p>\n' % \
            (v['type'].capitalize() if 'type' in v else 'Undefined',
                v['subtype'].capitalize() if 'subtype' in v else 'Undefined')
        if 'year' in v:
            doc += '\t\t\t<h3>Last released: %u<\h3>\n' % v['year']
        if 'releases' in v:
            doc += '\t\t\t<p><b>Released in years: </b>%s</p>\n' % \
                ', '.join(['%u' % y for y in v['releases']])
        if 'authors' in v and v['authors'] is not None:
            doc += '\t\t\t<p><b>Created by </b>%s</p>\n' % ', '.join(v[
                'authors'])
        if 'emails' in v and v['emails'] is not None:
            doc += '\t\t\t<p><b>Contact: </b></p>\n\n\t\t\t\t<ul>\n%s\n' % \
                ''.join(['\t\t\t\t<li><a href="mailto:%s">%s &lt;%s&gt;</li>\n' %
                         (em[0], em[1], em[0]) for em in v['emails']])
            doc += '\t\t\t\t</ul>\n'
        if 'license' in v:
            try:
                doc += '\t\t\t<p><b>License:</b> %s%s%s</p>\n' % (
                    ('<a href="%s" target="_blank">' % v['license']['url']) if
                    'url' in v['license'] else '', v['license']['name'], '</a>'
                    if 'url' in v['license'] else '')
            except KeyError:
                _log('Wrong license format for `%s`.' % k)
        if 'urls' in v:
            for uk, uv in iteritems(v['urls']):
                if len(uv) > 0 and uk != 'omictools':
                    try:
                        doc += '\t\t\t<h3>%s</h3>\n' % (uk.capitalize())
                        doc += '\t\t\t<ul>\n'
                        for a in uv:
                            doc += (
                                '\t\t\t\t<li><a href="%s" '
                                'target="_blank">%s</a></li>\n' % (
                                    a, a
                                )
                            )
                        doc += '\t\t\t</ul>\n'
                    except UnicodeDecodeError:
                        sys.stdout.write('UnicdeDecodeError at %s\n' % k)
                        sys.stdout.flush()
        if 'pubmeds' in v:
            doc += '\t\t\t<h3>PubMed</h3>\n'
            doc += '\t\t\t<ul>\n'
            for pmid in v['pubmeds']:
                doc += '\t\t\t\t<li><a href="%s" '\
                    'target="_blank">%s</a></li>\n' % (
                        'http://www.ncbi.nlm.nih.gov/pubmed/%u' % pmid,
                        'http://www.ncbi.nlm.nih.gov/pubmed/%u' % pmid
                    )
            doc += '\t\t\t</ul>\n'
        if ('urls' in v and 'omictools' in v['urls']) or 'pathguide' in v:
            doc += '\t\t\t<h3>Collections</h3>\n\t\t\t<ul>'
            if 'omictools' in v['urls']:
                doc += '\t\t\t<li><a href="%s" target="_blank">OmicTools</a></li>\n' % \
                    v['urls']['omictools'][0]
            if 'pathguide' in v:
                doc += '\t\t\t<li><a href="%s" target="_blank">PathGuide</a></li>\n' % \
                    (urls.urls['pathguide']['url'] % v['pathguide'])
            doc += '\t\t\t</ul>\n'
        if 'taxons' in v:
            doc += '<p><b>Taxons: </b><em>%s</em></p>' % \
                ', '.join(['%s%s' % (t[0].upper(), t[1:])
                           for t in v['taxons']])
        if 'size' in v and type(v['size']) is dict and \
                v['size']['nodes'] is not None and v['size']['edges'] is not None:
            doc += '<p><b>Nodes: </b>%s, <b>Edges:</b>%s</p>' % (
                v['size']['nodes'], v['size']['edges'])
        if 'data_import' in v:
            doc += '\t\t\t<p><b>Direct data import from: </b>%s</p>\n' % \
                ', '.join(v['data_import'])
        if 'includes' in v:
            doc += '\t\t\t<p><b>Includes data from: </b>%s</p>\n' % \
                ', '.join(v['includes'])
        if 'descriptions' in v or 'notes' in v:
            doc += '\t\t\t<h3>Quotes</h3>\n'
            if 'descriptions' in v:
                doc += '\t\t\t\t<div class="quotebox box">\n'
                pars = v['descriptions'][0].split('\n')
                for p in pars:
                    p = p.strip()
                    if len(p) > 0:
                        doc += '\t\t\t\t<p>%s</p>\n' % p
                doc += '\t\t\t\t</div>\n'
            if 'notes' in v:
                doc += '\t\t\t\t<div class="quotebox box">\n'
                pars = v['notes'][0].split('\n')
                for p in pars:
                    p = p.strip()
                    if len(p) > 0:
                        doc += '\t\t\t\t<p>%s</p>\n' % p
                doc += '\t\t\t\t</div>\n'
        if 'data_integration' in v:
            doc += '\t\t\t<p><b>Data integration in '\
                '<span class="code">pypath:</span></b> %s</p>' % \
                v['data_integration']
        if 'pypath' in v:
            doc += '\t\t\t<h3>Methods in <span class="code">pypath'\
                '</span></h3>\n'
            doc += '\t\t\t\t<div class="codebox box">\n'
            for cat in sorted(pypath_methods.keys()):
                name = pypath_methods[cat]
                if cat in v['pypath']:
                    doc += '\t\t\t\t\t<p>%s</p>\n\t\t\t\t\t<ul>\n' % name
                    for met in v['pypath'][cat]:
                        doc += '\t\t\t\t\t\t<li><span class="code">%s'\
                            '</span></li>\n' % met
                    doc += '\t\t\t\t\t</ul>\n'
            doc += '\t\t\t\t</div>\n'

    return _html.default_template(doc, title, title)


def write_html(filename='resources.html'):
    '''
    Saves the HTML descriptions to custom local file.
    '''
    html = gen_html()
    # with codecs.open(filename, encoding = 'utf-8', mode = 'w') as f:
    with open(filename, 'w') as f:
        f.write(html)


def resource_list_latex(filename='resource-list.tex',
                        latex_hdr=True,
                        fontsize=8,
                        font='HelveticaNeueLTStd-LtCn'):
    '''
    Generates Supplementary Table 3 (The list of the 52 resources considered) for the article.
    '''
    _latex_hdr = r'''\documentclass[a4paper,%upt]{extarticle}
        \usepackage{fontspec}
        \usepackage{xunicode}
        \usepackage{polyglossia}
        \setdefaultlanguage{english}
        \usepackage{xltxtra}
        \usepackage{microtype}
        \usepackage[margin=5pt,portrait,paperwidth=15cm,paperheight=18cm]{geometry}
        \usepackage{amsmath}
        \usepackage{amssymb}
        \usepackage{textcomp}
        \usepackage[usenames,dvipsnames,svgnames,table]{xcolor}
        \usepackage{color}
        \usepackage{booktabs}
        \usepackage{tabularx}
        \setmainfont{%s}
        \definecolor{grey875}{gray}{0.125}
        \begin{document}
        \color{grey875}
        \thispagestyle{empty}
        \vfill
    ''' % (fontsize, font) if latex_hdr else ''
    _latex_end = r'''
            \end{document}
        ''' if latex_hdr else ''
    tex = r'''\begin{tabularx}{0.94\textwidth}{>{\raggedright\scriptsize\arraybackslash\hsize=.15\hsize}X>{\raggedright\scriptsize\arraybackslash\hsize=.35\hsize}X>{\raggedright\scriptsize\arraybackslash\hsize=.15\hsize}X>{\raggedright\scriptsize\arraybackslash\hsize=.15\hsize}X}
    \toprule
    Resource name & Class, subclass & Resource name & Class, subclass \\
    \midrule
    '''
    res = sorted(
        [(v['label'] if 'label' in v else k,
          '%s, %s' % (v['type'].capitalize(), v['subtype'].capitalize())
          if 'type' in v and 'subtype' in v else '')
         for k, v in iteritems(descriptions)],
        key=lambda x: x[0].lower())
    if len(res) % 2 != 0:
        res.append('')
    res2 = zip(res[:int(len(res) / 2)], res[int(len(res) / 2):])
    for r in res2:
        tex += r'%s & %s & %s & %s \\' % (
            r[0][0].replace('&', '\&'),
            r[0][1].replace('&', '\&'),
            (r[1][0] if len(r[1]) else '').replace('&', '\&'),
            (r[1][1] if len(r[1]) else '').replace('&', '\&'),
        ) + '\n'
    tex += r'\bottomrule' + '\n'
    tex += r'\end{tabularx}' + '\n'
    with open(filename, 'w') as f:
        f.write('%s%s%s' % (_latex_hdr if latex_hdr else '', tex, _latex_end
                            if latex_hdr else ''))


def export_licenses(outfile = 'licenses.tsv'):

    hdr = [
        'Name',
        'License',
        'License URL',
        'Contact',
    ]
    rows = []

    for k, v in iteritems(descriptions):

        name = v['label'] if 'label' in v else k
        license_name = v['license']['name'] if 'license' in v else ''
        license_url = (
            v['license']['url']
                if 'license' in v and 'url' in v['license'] else
            ''
        )
        emails = (
            ','.join('%s <%s>' % tuple(reversed(e)) for e in v['emails'])
            if 'emails' in v else ''
        )

        rows.append([
            name,
            license_name,
            license_url,
            emails,
        ])

    with open(outfile, 'w') as fp:

        _ = fp.write('\t'.join(hdr) + '\n')
        _ = fp.write('\n'.join(
            '\t'.join(row) for row in rows
        ))
