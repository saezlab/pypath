#!/usr/bin/env python2
# -*coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
#
#  Copyright (c) 2014-2015 EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

#http://www.ijbs.com/v06p0051.htm
#http://www.nature.com/cddis/journal/v4/n8/full/cddis2013292a.html

import codecs
import bs4

__all__ = ['descriptions', 'gen_html', 'write_html']

descriptions = {
    'Lit13': {
        'year': 2013,
        'label': 'Lit-BM-13: High-quality non-systematic Literature dataset',
        'urls': {
            'articles': [
                'http://www.cell.com/cell/abstract/S0092-8674(14)01422-6'
            ],
            'webpages': [
                'http://interactome.dfci.harvard.edu/H_sapiens/'
            ]
        },
        'authors': ['CCSB'],
        'pubmeds': [25416956],
        'descriptions': [
            u'''
            In 2013, we extracted interaction data from BIND, BioGRID, DIP, HPRD, MINT, IntAct, and PDB to generate a high-quality binary literature dataset comprising ~11,000 protein-protein interactions that are binary and supported by at least two traceable pieces of evidence (publications and/or methods) (Rolland et al Cell 2014). Although this dataset does not result from a systematic investigation of the interactome search space and should thus be used with caution for any network topology analyses, it represents valuable interactions for targeted studies and is freely available to the research community through the search engine or via download. 
            '''
        ],
        'emails': [('Michael_Calderwood@dfci.harvard.edu', 'Michael Calderwood')],
        'type': 'high-throughput',
        'subtype': 'yeast 2 hybrid',
        'omnipath': False
    },
    'ELM': {
        'year': 2014,
        'releases': [2003, 2008, 2009, 2012, 2013, 2014],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['ELM Consortium'],
        'label': 'ELM',
        'color': '',
        'urls': {
            'webpages': [
                'http://elm.eu.org/'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/40/D1/D242.long',
                'http://nar.oxfordjournals.org/content/42/D1/D259.long'
            ]
        },
        'pubmeds': [22110040, 24214962],
        'emails': [('feedback@elm.eu.org', 'ELM Team'), ('gibson@embl.de', 'Toby Gibson')],
        'type': 'literature curated',
        'subtype': 'post-translational modifications',
        'omnipath': True
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
            'webpages': [
                'http://bicresources.jcbose.ac.in/ssaha4/lmpid/index.php'
            ],
            'articles': [
                'http://database.oxfordjournals.org/content/2015/bav014.long'
            ]
        },
        'pubmeds': [25776024],
        'emails': [('ssaha4@jcbose.ac.in', 'Sudipto Saha')],
        'type': 'literature curated',
        'subtype': 'post-translational modifications',
        'omnipath': True
    },
    'PDZBase': {
        'year': 2004,
        'releases': [2004],
        'authors': ['Weinstein Group'],
        'urls': {
            'webpages': [
                'http://abc.med.cornell.edu/pdzbase'
            ],
            'articles': [
                'http://bioinformatics.oxfordjournals.org/content/21/6/827.long'
            ]
        },
        'pubmeds': [15513994],
        'taxons': ['human'],
        'color': None,
        'label': 'PDZBase',
        'descriptions': [
            u'''
            PDZBase is a database that aims to contain all known PDZ-domain-mediated protein-protein interactions. Currently, PDZBase contains approximately 300 such interactions, which have been manually extracted from &gt;200 articles.
            PDZBase currently contains ∼300 interactions, all of which have been manually extracted from the literature, and have been independently verified by two curators. The extracted information comes from in vivo (co-immunoprecipitation) or in vitro experiments (GST-fusion or related pull-down experiments). Interactions identified solely from high throughput methods (e.g. yeast two-hybrid or mass spectrometry) were not included in PDZBase. Other prerequisites for inclusion in the database are: that knowledge of the binding sites on both interacting proteins must be available (for instance through a truncation or mutagenesis experiment); that interactions must be mediated directly by the PDZ-domain, and not by any other possible domain within the protein.
            '''
        ],
        'emails': [('haw2002@med.cornell.edu', 'Harel Weinstein'), ('pdzbase@med.cornell.edu' , 'PDZBase Team')],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True
    },
    'Guide2Pharmacology': {
        'year': 2015,
        'releases': [2007, 2008, 2009, 2011, 2015],
        'size': None,
        'authors': None,
        'label': 'Guide to Pharmacology',
        'full_name': 'Guide to Pharmacology',
        'color': None,
        'pubmeds': [24234439],
        'urls': {
            'webpages': [
                'http://www.guidetopharmacology.org/'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/42/D1/D1098.long',
                'http://onlinelibrary.wiley.com/doi/10.1111/j.1476-5381.2011.01649_1.x/full'
            ]
        },
        'descriptions': [
            u'''
            Presently, the resource describes the interactions between target proteins and 6064 distinct ligand entities (Table 1). Ligands are listed against targets by their action (e.g. activator, inhibitor), and also classified according to substance types and their status as approved drugs. Classes include metabolites (a general category for all biogenic, non-peptide, organic molecules including lipids, hormones and neurotransmitters), synthetic organic chemicals (e.g. small molecule drugs), natural products, mammalian endogenous peptides, synthetic and other peptides including toxins from non-mammalian organisms, antibodies, inorganic substances and other, not readily classifiable compounds. 
            The new database was constructed by integrating data from IUPHAR-DB and the published GRAC compendium. An overview of the curation process is depicted as an organizational flow chart in Figure 2. New information was added to the existing relational database behind IUPHAR-DB and new webpages were created to display the integrated information. For each new target, information on human, mouse and rat genes and proteins, including gene symbol, full name, location, gene ID, UniProt and Ensembl IDs was manually curated from HGNC, the Mouse Genome Database (MGD) at Mouse Genome Informatics (MGI), the Rat Genome Database (RGD), UniProt and Ensembl, respectively. In addition, ‘Other names’, target-specific fields such as ‘Principal transduction’, text from the ‘Overview’ and ‘Comments’ sections and reference citations (downloaded from PubMed; http://www.ncbi.nlm.nih.gov/pubmed) were captured from GRAC and uploaded into the database against a unique Object ID.
            '''
        ],
        'emails': [('enquiries@guidetopharmacology.org', 'Guide to Pharmacology Team'), ('tony.harmar@ed.ac.uk', 'Tony Harmar')],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True
    },
    'phosphoELM': {
        'year': 2010,
        'releases': [2004, 2007, 2010],
        'urls': {
            'webpages': [
                'http://phospho.elm.eu.org/'
            ],
            'articles': [
                'http://www.biomedcentral.com/1471-2105/5/79',
                'http://nar.oxfordjournals.org/content/36/suppl_1/D240.full',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D261'
            ]
        },
        'pubmeds': [15212693, 17962309, 21062810],
        'descriptions': [
            u'''
            Phospho.ELM http://phospho.elm.eu.org is a new resource containing experimentally verified phosphorylation sites manually curated from the literature and is developed as part of the ELM (Eukaryotic Linear Motif) resource. Phospho.ELM constitutes the largest searchable collection of phosphorylation sites available to the research community. The Phospho.ELM entries store information about substrate proteins with the exact positions of residues known to be phosphorylated by cellular kinases. Additional annotation includes literature references, subcellular compartment, tissue distribution, and information about the signaling pathways involved as well as links to the molecular interaction database MINT. Phospho.ELM version 2.0 contains 1703 phosphorylation site instances for 556 phosphorylated proteins. (Diella 2004)
            ''',
            u'''
            Phospho.ELM is a manually curated database of eukaryotic phosphorylation sites. The resource includes data collected from published literature as well as high-throughput data sets. The current release of the Phospho.ELM data set (version 7.0, July 2007) contains 4078 phospho-protein sequences covering 12 025 phospho-serine, 2362 phospho-threonine and 2083 phospho-tyrosine sites with a total of 16 470 sites. 
            For each phospho-site we report if the phosphorylation evidence has been identified by small-scale analysis (low throughput; LTP) that typically focus on one or a few proteins at a time or by large-scale experiments (high throughput; HTP), which mainly apply MS techniques. It is noteworthy that in our data set there is a small overlap between instances identified by LTP and HTP experiments. (Diella 2007)
            ''',
            u'''
            The current release of the Phospho.ELM data set (version 9.0) contains more than 42 500 non-redundant instances of phosphorylated residues in more than 11 000 different protein sequences (3370 tyrosine, 31 754 serine and 7449 threonine residues). For each phosphosite we report whether the phosphorylation evidence has been identified by small-scale analyses (low-throughput, LTP) and/or by large-scale experiments (high-throughput, HTP), which mainly apply MS techniques. The majority of the protein instances from Phospho. ELM are vertebrate (mostly Homo sapiens (62%) and Mus musculus (16%)) though 22% are from other species, mainly Drosophila melanogaster (13%) and Caenorhabditis elegans (7%). In total, more than 300 different kinases have been annotated and a document providing additional information about all kinases annotated in Phospho.ELM can be found at http://phospho.elm.eu.org/kinases.html. (Dinkel 2010)
            '''
        ],
        'emails': [('toby.gibson@embl.de', 'Toby Gibson')],
        'type': 'Literature curated',
        'subtype': 'PTM',
        'omnipath': True
    },
    'DOMINO': {
        'year': 2006,
        'releases': [2006],
        'authors': ['Cesareni Group'],
        'urls': {
            'webpages': [
                'http://mint.bio.uniroma2.it/domino/search/searchWelcome.do'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/35/suppl_1/D557.long'
            ]
        },
        'pubmeds': [17135199],
        'taxons': ['human', 'Metazoa', 'yeast'],
        'descriptions': [
            u'''
            DOMINO aims at annotating all the available information about domain-peptide and domain–domain interactions. The core of DOMINO, of July 24, 2006 consists of more than 3900 interactions extracted from peer-reviewed articles and annotated by expert biologists. A total of 717 manuscripts have been processed, thus covering a large fraction of the published information about domain–peptide interactions. The curation effort has focused on the following domains: SH3, SH2, 14-3-3, PDZ, PTB, WW, EVH, VHS, FHA, EH, FF, BRCT, Bromo, Chromo and GYF. However, interactions mediated by as many as 150 different domain families are stored in DOMINO.
            The curation process follows the PSI-MI 2.5 standard but with special emphasis on the mapping of the interaction to specific protein domains of both participating proteins. This is achieved by paying special attention to the shortest protein fragment that was experimentally verified as sufficient for the interaction. Whenever the authors report only the name of the domain mediating the interaction (i.e. SH3, SH2...), without stating the coordinates of the experimental binding range, the curator may choose to enter the coordinates of the Pfam domain match in the protein sequence. Finally whenever the information is available, any mutation or post-translational modification affecting the interaction affinity is noted in the database.
            '''
        ],
        'emails': [('giovanni.cesareni@uniroma2.it', 'Gianni Cesareni')],
        'type': 'Literature curated',
        'subtype': 'PTM',
        'omnipath': True
    },
    'dbPTM': {
        'year': 2012,
        'releases': [2005, 2009, 2012],
        'authors': ['ISBLab'],
        'urls': {
            'webpages': [
                'http://dbptm.mbc.nctu.edu.tw/'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/41/D1/D295.long',
                'http://www.biomedcentral.com/1756-0500/2/111',
                'http://nar.oxfordjournals.org/content/34/suppl_1/D622.long'
            ]
        },
        'pubmeds': [16381945, 19549291, 23193290],
        'taxons': ['human'],
        'descriptions': [
            u'''
            Due to the inaccessibility of database contents in several online PTM resources, a total 11 biological databases related to PTMs are integrated in dbPTM, including UniProtKB/SwissProt, version 9.0 of Phospho.ELM, PhosphoSitePlus, PHOSIDA, version 6.0 of O-GLYCBASE, dbOGAP, dbSNO, version 1.0 of UbiProt, PupDB, version 1.1 of SysPTM and release 9.0 of HPRD.
            With the high throughput of MS-based methods in post-translational proteomics, this update also includes manually curated MS/MS-identified peptides associated with PTMs from research articles through a literature survey. First, a table list of PTM-related keywords is constructed by referring to the UniProtKB/SwissProt PTM list (http://www.uniprot.org/docs/ptmlist.txt) and the annotations of RESID (28). Then, all fields in the PubMed database are searched based on the keywords of the constructed table list. This is then followed by downloading the full text of the research articles. For the various experiments of proteomic identification, a text-mining system is developed to survey full-text literature that potentially describes the site-specific identification of modified sites. Approximately 800 original and review articles associated with MS/MS proteomics and protein modifications are retrieved from PubMed (July 2012). Next, the full-length articles are manually reviewed for precisely extracting the MS/MS peptides along with the modified sites. Furthermore, in order to determine the locations of PTMs on a full-length protein sequence, the experimentally verified MS/MS peptides are then mapped to UniProtKB protein entries based on its database identifier (ID) and sequence identity. In the process of data mapping, MS/MS peptides that cannot align exactly to a protein sequence are discarded. Finally, each mapped PTM site is attributed with a corresponding literature (PubMed ID).
            ''',
            u'''
            Webpage was inaccessile on the 23rd of October 2015.
            '''
        ],
        'emails': [('francis@saturn.yzu.edu.tw', 'Hsien-Da Huang'), 
            ('bryan@mail.nctu.edu.tw', 'Hsien-Da Huang')],
        'type': 'Literature curated',
        'subtype': 'PTM',
        'omnipath': True
    },
    'Signor': {
        'year': 2015,
        'releases': [2015],
        'urls': {
            'webpages': ['http://signor.uniroma2.it/'],
            'articles': ['http://f1000research.com/posters/1098359']
        },
        'descriptions': [
            u'''
            SIGNOR, the SIGnaling Network Open Resource, organizes and stores in a structured format signaling information published in the scientific literature. The captured information is stored as binary causative relationships between biological entities and can be represented graphically as activity flow. The entire network can be freely downloaded and used to support logic modeling or to interpret high content datasets. The core of this project is a collection of more than 11000 manually-annotated causal relationships between proteins that participate in signal transduction. Each relationship is linked to the literature reporting the experimental evidence. In addition each node is annotated with the chemical inhibitors that modulate its activity. The signaling information is mapped to the human proteome even if the experimental evidence is based on experiments on mammalian model organisms.
            '''
        ],
        'authors': ['Cesareni Group'],
        'label': 'Signor',
        'color': '',
        'data_import': ['SignaLink2', 'PhosphoSite'],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True,
        'emails': [('perfetto@live.it', 'Livia Perfetto')],
    },
    'HuPho': {
        'year': 2015,
        'releases': [2012, 2015],
        'urls': {
            'webpages': ['http://hupho.uniroma2.it/'],
            'articles': ['http://onlinelibrary.wiley.com/doi/10.1111/j.1742-4658.2012.08712.x/full']
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
        'color': '',
        'type': 'high throughput and literature curated',
        'subtype': 'post-translational modification',
        'omnipath': False,
        'emails': [('perfetto@live.it', 'Livia Perfetto')]
    },
    'SignaLink2': {
        'year': 2015,
        'releases': [2010, 2012, 2015],
        'size': 0,
        'authors': ['NetBiol Group'],
        'label': 'SignaLink',
        'color': '',
        'pubmeds': [20542890, 23331499],
        'urls': {
            'webpages': [
                'http://signalink.org/'
            ],
            'articles': [
                'http://bioinformatics.oxfordjournals.org/content/26/16/2042.long',
                'http://www.biomedcentral.com/1752-0509/7/7'
            ]
        },
        'taxons': [
            'human',
            'D. melanogaster',
            'C. elegans'
        ],
        'descriptions': [
            u'''
            In each of the three organisms, we first listed signaling proteins and interactions from reviews (and from WormBook in C.elegans) and then added further signaling interactions of the listed proteins. To identify additional interactions in C.elegans, we examined all interactions (except for transcription regulation) of the signaling proteins listed in WormBase and added only those to SignaLink that we could manually identify in the literature as an experimentally verified signaling interaction. For D.melanogaster, we added to SignaLink those genetic interactions from FlyBase that were also reported in at least one yeast-2-hybrid experiment. For humans, we manually checked the reliability and directions for the PPIs found with the search engines iHop and Chilibot. 
            SignaLink assigns proteins to signaling pathways using the full texts of pathway reviews (written by pathway experts). While most signaling resources consider 5–15 reviews per pathway, SignaLink uses a total of 170 review papers, i.e. more than 20 per pathway on average. Interactions were curated from a total of 941 articles (PubMed IDs are available at the website). We added a small number of proteins based on InParanoid ortholog clusters. For curation, we used a self-developed graphical tool and Perl/Python scripts. The current version of SignaLink was completed in May 2008 based on WormBase (version 191), FlyBase (2008.6), Ensembl, UniProt and the publications listed on the website. 
            The curation protocol of SignaLink (Fig. 1A) contains several steps aimed specifically at reducing data and curation errors. We used reviews as a starting point, manually looked up interactions three times, and manually searched for interactions of known signaling proteins with no signaling interactions so far in the database. 
            '''
        ],
        'emails': [('korcsmaros@gmail.com', 'Tamas Korcsmaros'), ('tamas.korcsmaros@tgac.ac.uk', 'Tamas Korcsmaros')],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True
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
            'webpages': [
                'http://nrf2.elte.hu/'
            ],
            'articles': [
                'http://www.hindawi.com/journals/omcl/2013/737591/',
                'http://www.sciencedirect.com/science/article/pii/S0014579312003912'
            ]
        },
        'pubmeds': [22641035, 23710289],
        'taxons': [
            'human'
        ],
        'descriptions': [
            u'''
            From Korcsmaros 2010: ... we first listed signaling proteins and interactions from reviews and then added further signaling interactions of the listed proteins. We used reviews as a starting point, manually looked up interactions three times, and manually searched for interactions of known signaling proteins with no signaling interactions so far in the database.
            '''
        ],
        'emails': [('korcsmaros@gmail.com', 'Tamas Korcsmaros'), ('tamas.korcsmaros@tgac.ac.uk', 'Tamas Korcsmaros')],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True
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
            'webpages': [
                'http://autophagy-regulation.org/'
            ],
            'articles': [
                'http://www.tandfonline.com/doi/full/10.4161/15548627.2014.994346'
            ]
        },
        'taxons': [
            'human'
        ],
        'descriptions': [
            u'''
            From Korcsmaros 2010: ... we first listed signaling proteins and interactions from reviews and then added further signaling interactions of the listed proteins. We used reviews as a starting point, manually looked up interactions three times, and manually searched for interactions of known signaling proteins with no signaling interactions so far in the database.
            '''
        ],
        'emails': [('korcsmaros@gmail.com', 'Tamas Korcsmaros'), ('tamas.korcsmaros@tgac.ac.uk', 'Tamas Korcsmaros')],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True
    },
    'HPRD': {
        'year': 2010,
        'releases': [2002, 2005, 2009, 2010],
        'urls': {
            'webpages': [
                'http://www.hprd.org/'
            ],
            'articles': [
                'http://genome.cshlp.org/content/13/10/2363.long',
                'http://nar.oxfordjournals.org/content/34/suppl_1/D411.long',
                'http://nar.oxfordjournals.org/content/37/suppl_1/D767.long'
            ]
        },
        'pubmeds': [14525934, 16381900, 18988627],
        'descriptions': [
            u'''
            The information about protein-protein interactions was cataloged after a critical reading of the published literature. Exhaustive searches were done based on keywords and medical subject headings (MeSH) by using Entrez. The type of experiments that served as the basis for establishing protein-protein interactions was also annotated. Experiments such as coimmunoprecipitation were designated in vivo, GST fusion and similar “pull-down” type of experiments were designated in vitro, and those identified by yeast two-hybrid were annotated as yeast two-hybrid.
            Posttranslational modifications were annotated based on the type of modification, site of modification, and the modified residue. In addition, the upstream enzymes that are responsible for modifications of these proteins were reported if described in the articles. The most commonly known and the alternative subcellular localization of the protein were based on the literature. The sites of expression of protein and/or mRNA were annotated based on published studies.
            '''
        ],
        'emails': [('pandey@jhmi.edu', 'Akhilesh Pandey')],
        'type': 'literature curated',
        'subtype': 'post-translational modification',
        'omnipath': True
    },
    'ACSN': {
        'year': 2015,
        'releases': [2008, 2015],
        'authors': ['Curie'],
        'urls': {
            'webpages': [
                'https://acsn.curie.fr'
            ],
            'articles': [
                'http://www.nature.com/oncsis/journal/v4/n7/full/oncsis201519a.html',
                'http://msb.embopress.org/content/4/1/0174.long'
            ]
        },
        'pubmeds': [
            26192618,
            18319725
        ],
        'taxons': [
            'human'
        ],
        'descriptions': [
            u'''
            The map curator studies the body of literature dedicated to the biological process or molecular mechanism of interest. The initial sources of information are the major review articles from high-impact journals that represent the consensus view on the studied topic and also provide a list of original references. The map curator extracts information from review papers and represents it in the form of biochemical reactions in CellDesigner. This level of details reflects the ‘canonical’ mechanisms. Afterwards, the curator extends the search and analyses original papers from the list provided in the review articles and beyond. This information is used to enrich the map with details from the recent discoveries in the field. The rule for confident acceptance and inclusion of a biochemical reaction or a process is the presence of sufficient evidences from more than two studies, preferably from different scientific groups. The content of ACSN is also verified and compared with publicly available databases such as REACTOME, KEGG, WikiPathways, BioCarta, Cell Signalling and others to ensure comprehensive representation of consensus pathways and links on PMIDs of original articles confirmed annotated molecular interactions.
            ''',
            u'''
            CellDesigner 3.5 version was used to enter biological facts from a carefully studied selection of papers (see the whole bibliography on the web site with Supplementary information). Whenever the details of a biological fact could not be naturally expressed with CellDesigner standard notations, it was fixed and some solution was proposed. For example, we added a notation (co‐factor) to describe all the components intervening in the transcription of genes mediated by the E2F family proteins. 
            '''
        ],
        'emails': [('andrei.zinovyev@curie.fr', 'Andrei Zinovyev')],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True
    },
    'DOMINO': {
        'urls': {
            'articles': [
                'http://nar.oxfordjournals.org/content/35/suppl_1/D557.full'
            ],
            'webpages': [
                'http://mint.bio.uniroma2.it/domino/'
            ]
        },
        'pubmeds': [17135199],
        'descriptions': [
            u'''
            DOMINO is an open-access database comprising more than 3900 annotated experiments describing interactions mediated by protein-interaction domains. The curation effort aims at covering the interactions mediated by the following domains (SH3, SH2, 14-3-3, PDZ, PTB, WW, EVH, VHS, FHA, EH, FF, BRCT, Bromo, Chromo, GYF). However, interactions mediated by as many as 150 different domain families are stored in DOMINO.
            ''',
            u'''
            The curation process follows the PSI-MI 2.5 standard but with special emphasis on the mapping of the interaction to specific protein domains of both participating proteins. This is achieved by paying special attention to the shortest protein fragment that was experimentally verified as sufficient for the interaction. Whenever the authors report only the name of the domain mediating the interaction (i.e. SH3, SH2 ...), without stating the coordinates of the experimental binding range, the curator may choose to enter the coordinates of the Pfam domain match in the protein sequence. Finally whenever the information is available, any mutation or posttranslational modification affecting the interaction affinity is noted in the database.
            '''
        ],
        'taxons': [
            'human',
            'yeast',
            'C. elegans',
            'mouse',
            'rat',
            'HIV',
            'D. melanogaster',
            'A. thaliana',
            'X. laevis',
            'B. taurus',
            'G. gallus',
            'O. cuniculus',
            'Plasmodium falciparum'
        ],
        'emails': [('Cesareni@uniroma2.it', 'Gianni Cesareni')],
        'type': 'literature curated',
        'subtype': 'post-translational modification',
        'omnipath': True
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
        'taxons': [
            'human'
        ],
        'pubmeds': [22135292],
        'urls': {
            'articles': [
                'http://nar.oxfordjournals.org/content/40/D1/D331'
            ],
            'webpages': [
                'http://deathdomain.org/'
            ]
        },
        'emails': [('hyunho@ynu.ac.kr', 'Hyun Ho O')],
        'files': {
            'articles': [
                'DeathDomain_Kwon2011.pdf'
            ],
            'data': {
                'raw': [
                    'deathdomain.tsv'
                ],
                'processed': [
                    'deathdomain.sif'
                ]
            }
        },
        'taxons': [
            'human'
        ],
        'size': {
            'nodes': 99,
            'edges': 175
        },
        'identifiers': [
            'GeneSymbol'
        ],
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
        'subtype': 'pathway',
        'omnipath': True
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
            'webpages': [
                'http://www.trpchannel.org'
            ]
        },
        'emails': [('jhjeon2@snu.ac.kr', 'Ju-Hong Jeon')],
        'size': {
            'nodes': 468,
            'edges': 744
        },
        'pubmeds': [20851834, 23071747, 23677537],
        'files': {
            'articles': [
                'TRIP_Shin2012.pdf'
            ],
            'data': {
                'raw': [],
                'processed': [
                    'trip.sif'
                ]
            }
        },
        'taxons': [
            'human',
            'mouse',
            'rat'
        ],
        'identifiers': [
            'GeneSymbol'
        ],
        'descriptions': [
            u'''
            The literature on TRP channel PPIs found in the PubMed database serve as the primary information source for constructing the TRIP Database. First, a list of synonyms for the term ‘TRP channels’ was constructed from UniprotKB, Entrez Gene, membrane protein databases (Supplementary Table S2) and published review papers for nomenclature. Second, using these synonyms, a list of articles was obtained through a PubMed search. Third, salient articles were collected through a survey of PubMed abstracts and subsequently by search of full-text papers. Finally, we selected articles that contain evidence for physical binding among the proteins denoted. To prevent omission of relevant papers, we manually screened information in other databases, such as DIP, IntAct, MINT, STRING, BioGRID, Entrez Gene, IUPHAR-DB and ISI Web of Knowledge (from Thomson Reuters). All 277 articles used for database construction are listed in our database website.
            '''
        ],
        'notes': [
            u'''
            Good manually curated dataset focusing on TRP channel proteins, with ~800 binary interactions. The provided formats are not well suitable for bioinformatics use because of the non standard protein names, with greek letters and only human understandable formulas. Using HTML processing, and processing the data from 5-6 different tables, with couple hundreds of lines of code, one have a chance to compile a usable data table. 
            '''
        ],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True
    },
    'Awan2007': {
        'year': 2007,
        'size': 0,
        'authors': ['Wang Group'],
        'label': 'Awan 2007',
        'color': '',
        'data_import': ['BioCarta', 'CA1'],
        'contains': [
            'BioCarta',
            'CA1'
        ],
        'urls': {
            'articles': [
                'http://www.cancer-systemsbiology.org/Papers/iet-sb2007.pdf'
            ]
        },
        'emails': [('Edwin.Wang@cnrc-nrc.gc.ca', 'Edwin Wang')],
        'pubmeds': [17907678],
        'descriptions': [
            u'''
            To construct the human cellular signalling network, we manually curated signalling pathways from literature. The signalling data source for our pathways is the BioCarta database (http://www.biocarta.com/genes/allpathways.asp), which, so far, is the most comprehensive database for human cellular signalling pathways. Our curated pathway database recorded gene names and functions, cellular locations of each gene and relationships between genes such as activation, inhibition, translocation, enzyme digestion, gene transcription and translation, signal stimulation and so on. To ensure the accuracy and the consistency of the database, each referenced pathway was cross-checked by different researchers and finally all the documented pathways were checked by one researcher. In total, 164 signalling pathways were documented (supplementary Table 2). Furthermore, we merged the curated data with another literature-mined human cellular signalling network. As a result, the merged network contains nearly 1100 proteins (SupplementaryNetworkFile). To construct a signalling network, we considered relationships of proteins as links (activation or inactivation as directed links and physical interactions in protein complexes as neutral links) and proteins as nodes.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': False
    },
    'Cui2007': {
        'year': 2007,
        'authors': ['Wang Group'],
        'label': 'Cui 2007',
        'color': '',
        'data_import': ['Awan2007', 'CancerCellMap'],
        'pubmeds': [18091723],
        'contains': [
            'Awan2007',
            'CancerCellMap',
            'CA1',
            'BioCarta'
        ],
        'urls': {
            'articles': [
                'http://msb.embopress.org/content/3/1/152'
            ],
            'webpages': []
        },
        'emails': [('Edwin.Wang@cnrc-nrc.gc.ca', 'Edwin Wang')],
        'files': {
            'articles': [
                'Cui2007.pdf'
            ],
            'data': {
                'raw': [
                    'cui-network.xls',
                ],
                'processed': [
                    'cui.sif'
                ]
            }
        },
        'identifiers': [
            'EntrezGene'
        ],
        'size': {
            'edges': 4249,
            'nodes': 1528
        },
        'taxons': [
            'human'
        ],
        'descriptions': [
            u'''
            To build up the human signaling network, we manually curated the signaling molecules (most of them are proteins) and the interactions between these molecules from the most comprehensive signaling pathway database, BioCarta (http://www.biocarta.com/). The pathways in the database are illustrated as diagrams. We manually recorded the names, functions, cellular locations, biochemical classifications and the regulatory (including activating and inhibitory) and interaction relations of the signaling molecules for each signaling pathway. To ensure the accuracy of the curation, all the data have been crosschecked four times by different researchers. After combining the curated information with another literature‐mined signaling network that contains ∼500 signaling molecules (Ma'ayan et al, 2005)[this is the CA1], we obtained a signaling network containing ∼1100 proteins (Awan et al, 2007). We further extended this network by extracting and adding the signaling molecules and their relations from the Cancer Cell Map (http://cancer.cellmap.org/cellmap/), a database that contains 10 manually curated signaling pathways for cancer. As a result, the network contains 1634 nodes and 5089 links that include 2403 activation links (positive links), 741 inhibitory links (negative links), 1915 physical links (neutral links) and 30 links whose types are unknown (Supplementary Table 9). To our knowledge, this network is the biggest cellular signaling network at present.
            ''',
            u'''
            From Awan 2007: To construct the human cellular signalling network, we manually curated signalling pathways from literature. The signalling data source for our pathways is the BioCarta database (http://www.biocarta.com/genes/allpathways.asp), which, so far, is the most comprehensive database for human cellular signalling pathways. Our curated pathway database recorded gene names and functions, cellular locations of each gene and relationships between genes such as activation, inhibition, translocation, enzyme digestion, gene transcription and translation, signal stimulation and so on. To ensure the accuracy and the consistency of the database, each referenced pathway was cross-checked by different researchers and finally all the documented pathways were checked by one researcher.
            '''
        ],
        'notes': [
            u'''
            Excellent signaling network with good topology for all those who doesn't mind to use data of unknown origin. Supposedly a manually curated network, but data files doesn't include article references. Merging CA1 network with CancerCellMap and BioCarta (also without references) makes the origin of the data untraceable.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': False
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
                'http://www.biocarta.com/'
            ],
            'articles': []
        },
        'emails': [('info@biocarta.com', 'BioCarta Scientific Advisory Board')],
        'taxons': [
            'human'
        ],
        'descriptions': [
            u'''
            Community built pathway database based on expert curation.
            '''
        ],
        'notes': [
            u'''
            This resource includes a huge number of pathways, each curated by experts from a few reviews. The data is not available for download from the original webpage, only from second hand, for example from NCI-PID, in NCI-XML format. However, these files doesn't contain any references, which makes problematic the use of the BioCarta dataset. Also, some pathways are reviewed long time ago, possibly outdated.
            ''',
            '''
            The Company and the website looks like was abandoned around 2003-2006.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': False
    },
    'TLR': {
        'urls': {
            'articles': [
                'http://msb.embopress.org/content/2/1/2006.0015.long'
            ]
        },
        'type': 'literature curated',
        'subtype': 'model',
        'omnipath': False
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
        'color': '',
        'pubmeds': [16099987],
        'urls': {
            'articles': [
                'http://www.sciencemag.org/content/309/5737/1078.full'
            ],
            'webpages': []
        },
        'emails': [('ravi.iyengar@mssm.edu', 'Ravi Iyengar')],
        'taxons': [
            'human',
            'mouse'
        ],
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
        'subtype': 'pathway',
        'omnipath': True
    },
    'CancerCellMap': {
        'urls': {
            'articles': [],
            'webpages': [
                'http://www.pathwaycommons.org/pc-snapshot/current-release/tab_delim_network/by_source/'
            ]
        },
        'authors': ['Bader Lab'],
        'emails': [('gary.bader@utoronto.ca', 'Gary Bader')],
        'descriptions': [
            u'''
            Manually curated data, unpublished. A team of M.Sc. and Ph.D. biologists at the Institute of Bioinformatics in Bangalore, India read original research papers and hand-entered the pathway data into our database. The quality of the Cancer Cell Map pathways is very high. Half of the pathways were reviewed by experts at Memorial Sloan-Kettering Cancer Center and were found to contain only a few errors, which were subsequently fixed.
            '''
        ],
        'notes': [
            u'''
            One of the earliest manually curated datasets, now only available from second hand, e.g. from PathwayCommons. Included in many other resources. Contains binary interactions with PubMed references.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True
    },
    'HSN': {
        'year': 2014,
        'releases': [2009, 2010, 2011, 2012, 2013, 2014],
        'nodes': 6300,
        'edges': 63000,
        'authors': ['Wang Group'],
        'label': 'HumanSignalingNetwork',
        'color': '',
        'data_import': ['Cui2007', 'BioCarta', 'CST', 'NCI-PID', 'iHOP'],
        'urls':{
            'webpages': [
                'http://www.cancer-systemsbiology.org/dataandsoftware.htm',
                'http://www.bri.nrc.ca/wang/'
            ],
            'articles': []
        },
        'emails': [('Edwin.Wang@cnrc-nrc.gc.ca', 'Edwin Wang')],
        'taxons': [
            'human',
            'mouse',
            'rat'
        ],
        'contains': [
            'Cui2007', 
            'CancerCellMap',
            'Awan2007',
            'NetPath',
            'CA1',
            'NCI-PID',
            'BioCarta',
            'CST',
            'iHOP'
        ],
        'descriptions': [
            u'''
            Composed from multiple manually curated datasets, and contains own manual cuartion effort. Methods are unclear, and the dataset has not been published in reviewed paper. Based on the Cui et al 2007. 
            Wang Lab has manually curated human signaling data from literature since 2005. The data sources include BioCarta, CST Signaling pathways, NCI Pathway Interaction Database, IHOP, and many review papers. The contents are updated every year. 
            iHOP is not literature curated, but is a literature mining platform. 
            '''
        ],
        'notes': [
            u'''
            This network aims to merge multiple manually curated networks. Unfortunately a precise description of the sources and methods is missing. Also, the dataset doesn't include the references. Moreover, the data file misses header and key, so users can only guess about the meaning of columns and values.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': False
    },
    'Ataxia': {
        'year': 2010,
        'urls':{
            'webpages': [
                'http://franklin.imgen.bcm.tmc.edu/ppi/tables.html'
            ],
            'articles': [
                'http://hmg.oxfordjournals.org/content/20/3/510.long'
            ]
        },
        'taxons': [
            'human'
        ],
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
        'emails': [('Tong_Hao@dfci.harvard.edu', 'Tong Hao'), ('barabasi@gmail.com', 'Albert-Laszlo Barabasi')]
    },
    'Reactome': {
        'urls': {
            'webpages': [
                'http://reactome.org/'
            ],
            'articles': [
                'http://genomebiology.com/content/8/3/R39',
                'http://nar.oxfordjournals.org/content/42/D1/D472.long'
            ]
        },
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
        'subtype': 'reaction network',
        'omnipath': False,
        'emails': [('help@reactome.org', 'Reactome Team'), ('hhe@ebi.ac.uk', 'Henning Hermjakob')]
    },
    'Li2012': {
        'year': 2012,
        'urls': {
            'articles': [
                'http://genome.cshlp.org/content/22/7/1222'
            ]
        },
        'pubmeds': [22194470],
        'label': 'Li 2012',
        'authors': ['Wang Lab'],
        'emails': [('edwin.wang@cnrc-nrc.gc.ca', 'Edwin Wang')],
        'taxons': [
            'human'
        ],
        'descriptions': [
            u'''
            Human phosphotyrosine signaling network. 
            We manually collected the experimentally determined human TK–substrate interactions and substrate–SH2/PTB domain interactions from the literature (see Supplemental Materials), as well as the Phospho.ELM and PhosphoSitePlus databases. [71 references, 585 circuits]
            '''
        ],
        'type': 'high-throughput',
        'subtype': 'yeast 2 hybrid',
        'omnipath': False
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
            'HumanSignalingNetwork',
            'Cui2007',
            'CA1',
            'BioCarta',
            'Awan2007',
            'Li2012',
            'NCI-PID'
        ],
        'descriptions': [
            u'''
            The human signaling network (Version 4, containing more than 6,000 genes and more than 50,000 relations) includes our previous data obtained from manually curated signaling networks (Awan et al., 2007; Cui et al., 2007; Li et al., 2012) and by PID (http://pid.nci.nih.gov/) and our recent manual curations using the iHOP database (http://www.ihop-net.org/UniPub/iHOP/).
            '''
        ],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': False
    },
    'AlzPathway': {
        'year': 2012,
        'releases': [2012],
        'size': 0,
        'authors': ['Tokyo Bioinf'],
        'label': 'AlzPathway',
        'color': '',
        'urls': {
            'articles': [
                'http://www.biomedcentral.com/1752-0509/6/52',
                'http://onlinelibrary.wiley.com/doi/10.1038/clpt.2013.37/epdf'
            ],
            'webpages': [
                'http://alzpathway.org/AlzPathway.html'
            ]
        },
        'pubmeds': [23511713, 22647208],
        'emails': [('ogishima@sysmedbio.org', 'Soichi Ogishima'), ('info@alzpathway.org', 'Soichi Ogishima')],
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
        'subtype': 'pathway',
        'omnipath': True
    },
    'MPPI': {
        'year': 2005,
        'releases': [2000, 2005],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['MIPS Munich'],
        'label': 'MPPI',
        'color': '',
        'urls': {
            'articles': [
                'http://bioinformatics.oxfordjournals.org/content/21/6/832'
            ],
            'webpages': [
                'http://mips.helmholtz-muenchen.de/proj/ppi/'
            ]
        },
        'emails': [('p.pagel@wzw.tum.de', 'Philipp Pagel'), ('d.frishman@helmholtz-muenchen.de', 'Dmitrij Frishman')],
        'pubmeds': [15531608],
        'taxons': [
            'human',
            'mammalia'
        ],
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
        'omnipath': True
    },
    'Negatome': {
        'year': 2013,
        'relases': [2009],
        'urls': {
            'articles': [
                'http://nar.oxfordjournals.org/content/38/suppl_1/D540.long',
                'http://nar.oxfordjournals.org/content/42/D1/D396.long'
            ],
            'webpages': [
                'http://mips.helmholtz-muenchen.de/proj/ppi/negatome/'
            ]
        },
        'pubmeds': [24214996, 19920129],
        'emails': [('d.frishman@helmholtz-muenchen.de', 'Dmitrij Frishman')],
        'descriptions': [
            u'''
            Annotation of the manual dataset was performed analogous to the annotation of protein–protein interactions and protein complexes in previous projects published by our group. Information about NIPs was extracted from scientific literature using only data from individual experiments but not from high-throughput experiments. Only mammalian proteins were considered. Data from high-throughput experiments were omitted in order to maintain the highest possible standard of reliability.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'negative',
        'omnipath': True
    },
    'Macrophage': {
        'year': 2010,
        'urls': {
            'articles': [
                'http://www.biomedcentral.com/1752-0509/4/63'
            ],
            'webpages': []
        },
        'emails': [('tom.freeman@roslin.ed.ac.uk', 'Tom Freeman')],
        'pubmeds': [20470404],
        'descriptions': [
            u'''
            Ongoing analysis of macrophage-related datasets and an interest in consolidating our knowledge of a number of signalling pathways directed our choice of pathways to be mapped (see Figure 1). Public and propriety databases were initially used as resources for data mining, but ultimately all molecular interaction data was sourced from published literature. Manual curation of the literature was performed to firstly evaluate the quality of the evidence supporting an interaction and secondly, to extract the necessary and additional pieces of information required to 'understand' the pathway and construct an interaction diagram. We have drawn pathways based on our desire to model pathways active in a human macrophage and therefore all components have been depicted using standard human gene nomenclature (HGNC). However, our understanding of the pathway components and the interactions between them, have been drawn largely from a consensus view of literature knowledge. As such the pathways presented here are based on data derived from a range of different cellular systems and mammalian species (human and mouse).
            '''
        ],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True
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
        'label': 'NetPath',
        'color': '',
        'data_import': [
            'CancerCellMap'
        ],
        'includes': ['CancerCellMap'],
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
            'webpages': [
                'http://netpath.org/'
            ]
        },
        'pubmeds': [20067622, 21959865, 21742767, 21996254, 24255551, 22684822,
            22649723, 23255051, 23161412, 23504413, 23631681, 23606317, 24573880,
            24584707, 24829797],
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
        'subtype': 'reaction network',
        'omnipath': True
    },
    'InnateDB': {
        'year': 2015,
        'urls': {
            'articles': [
                'http://msb.embopress.org/content/4/1/218.long',
                'http://www.biomedcentral.com/1752-0509/4/117',
                'http://nar.oxfordjournals.org/content/41/D1/D1228.long'
            ],
            'webpages': [
                'http://www.innatedb.com/'
            ]
        },
        'authors': ['Brinkman Lab', 'Hancock Lab', 'Lynn Group'],
        'emails': [('innatedb-mail@sfu.ca', 'InnateDB Team'), ('david.lynn@teagasc.ie', 'David Lynn')],
        'license': 'http://innatedb.com/license.jsp',
        'pubmeds': [
            20727158,
            23180781,
            18766178
        ],
        'releases': [2008, 2010, 2013, 2014, 2015],
        'descriptions': [
            u'''
            To date, the InnateDB curation team has reviewed more than 1000 publications and curated more than 3500 innate immunity-relevant interactions, richly annotating them in terms of the experimental evidence and the context in which they occur.
            To date, InnateDB manual curation has prioritized molecules that are well-described members of key innate immunity signaling pathways, including the TLR pathways, the NF-kB pathway, MAPK signaling pathway, JNK signaling pathway, NOD-like receptor pathway and the RIG-I antiviral pathway (Kanneganti et al, 2007; Lee and Kim, 2007; Thompson and Locarnini, 2007). We have then curated experimentally verified interactions between these molecules and any other molecule, regardless of whether the interacting molecule has any known role in innate immunity.
            InnateDB project has had a full-time curation team employed for more than three years. As of February 15th 2010, there were 11,786 InnateDB-curated molecular interactions in InnateDB (>3,000 published articles reviewed). Currently, InnateDB only curates interactions involving human and mouse molecules, with the majority of curated interactions (72% or 8,569 interactions) involving human molecules (although there has been no specific focus on human as opposed to mouse). Additionally, there are 1,005 hybrid interactions involving both human and mouse participants. Curated interactions are primarily protein-protein interactions (9,244 interactions), however, there are also almost 2,500 protein-DNA interactions and a small, but important, number of RNA interactions (mainly microRNAs).
            As of September 2012, our curation team has reviewed >4300 publications, and >18 000 interactions of relevance to innate immunity have been annotated.
            '''
        ],
        'notes': [
            u'''
            Probably the largest manually curated binary protein interaction dataset, developed by a dedicated full time team of curators. Formats are clear and accessible, comprising UniProt IDs, PubMed references, experimental evidences and mechanisms.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'interaction',
        'omnipath': True
    },
    'CORUM': {
        'year': 2009,
        'releases': [2007, 2009],
        'urls': {
            'articles': [
                'http://nar.oxfordjournals.org/content/36/suppl_1/D646.long',
                'http://nar.oxfordjournals.org/content/38/suppl_1/D497.long'
            ],
            'webpages': [
                'http://mips.helmholtz-muenchen.de/genre/proj/corum'
            ]
        },
        'emails': [('andreas.ruepp@helmholtz-muenchen.de', 'Andreas Ruepp')],
        'pubmeds': [19884131, 17965090],
        'taxons': [
            'human',
            'mouse',
            'rat'
        ],
        'descriptions': [
            u'''
            The CORUM database is a collection of experimentally verified mammalian protein complexes. Information is manually derived by critical reading of the scientific literature from expert annotators. Information about protein complexes includes protein complex names, subunits, literature references as well as the function of the complexes.
            In order to provide a high-quality dataset of mammalian protein complexes, all entries are manually created. Only protein complexes which have been isolated and characterized by reliable experimental evidence are included in CORUM. To be considered for CORUM, a protein complex has to be isolated as one molecule and must not be a construct derived from several experiments. Also, artificial constructs of subcomplexes are not taken into account. Since information from high-throughput experi ments contains a significant fraction of false-positive results, this type of data is excluded. References for relevant articles were mainly found in general review articles, cross-references to related protein complexes within analysed literature and comments on referenced articles in UniProt.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'complexes',
        'omnipath': True
    },
    'CST': {
        'year': 2015,
        'releases': [2005, 2015],
        'nodes': None,
        'edges': None,
        'authors': ['CST'],
        'label': 'CST pathways',
        'color': '',
        'full_name': 'Cell Signaling Technology pathways',
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
        'subtype': 'pathway',
        'omnipath': False
    },
    'DIP': {
        'year': 2014,
        'releases': [2000, 2001, 2002, 2004, 2011, 2014],
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
            'webpages': [
                'http://dip.doe-mbi.ucla.edu/dip/Main.cgi'
            ]
        },
        'pubmeds': [10592249, 11125102, 11752321, 14681454],
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
        'omnipath': True,
        'emails': [('david@mbi.ucla.edu', 'David Eisenberg')]
    },
    'DEPOD': {
        'year': 2014,
        'releases': [2013, 2014],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['EMBL & EMBL-EBI'],
        'label': 'DEPOD',
        'color': '',
        'urls': {
            'articles': [
                'http://stke.sciencemag.org/content/6/275/rs10.long',
                'http://nar.oxfordjournals.org/content/43/D1/D531.long'
            ],
            'webpages': [
                'http://www.koehn.embl.de/depod/index.php'
            ]
        },
        'pubmeds': [23674824, 25332398],
        'taxons': [
            'human'
        ],
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
        'emails': [('koehn@embl.de', 'Maja Kohn')]
    },
    'PhosphoPoint': {
        'year': 2008,
        'urls': {
            'articles': [
                'http://bioinformatics.oxfordjournals.org/content/24/16/i14.long'
            ],
            'webpages': [
                'http://kinase.bioinformatics.tw/'
            ]
        },
        'taxons': [
            'human'
        ],
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
        'omnipath': False,
        'emails': [('cyhuang5@ym.edu.tw' , 'Chi-Ying F. Huang'), ('kmchao@csie.ntu.edu.tw', 'Kun-Mao Chao')]
    },
    'PANTHER': {
        'year': 2014,
        'releases': [2000, 2001, 2002, 2003, 2005, 2006, 2010, 2011, 2012, 2014],
        'urls': {
            'articles': [
                'http://link.springer.com/protocol/10.1007%2F978-1-60761-175-2_7#section=82252&page=1',
                'http://nar.oxfordjournals.org/content/35/suppl_1/D247.long'
            ],
            'webpages': [
                'http://www.pantherdb.org/'
            ]
        },
        'pubmeds': [17130144, 19597783],
        'descriptions': [
            u'''
            References are captured at three levels. First, each pathway as a whole requires a reference. For signaling pathways, at least three references, usually review papers, are required in order to provide a more objective view of the scope of the pathway. For metabolic pathways, a textbook reference is usually sufficient. Second, references are often associated to each molecule class in the pathway. Most of these references are OMIM records or review papers. Third, references are provided to support association of specific protein sequences with a particular molecule class, e.g., the SWISS-PROT sequence P53_HUMAN annotated as an instance of the molecule class ‘‘P53’’ appearing in the pathway class ‘‘P53 pathway’’. These are usually research papers that report the experimental evidence that a particular protein or gene participates in the reactions represented in the pathway diagram.
            There are three major properties that make this infrastructure differ from other pathway curation systems, such as from Reactome and EcoCyc. First, the pathway diagrams are drawn with CellDesigner software. There are two advantages to using CellDesigner. First, controlled graphical notations are used to draw the pathway diagram, and the software automatically creates a computational representation that is compatible with the SBML standard. Second, a pathway diagram can be viewed with an exact, one-to-one correspondence with the ontological representation of the pathways stored in the back-end. The second property is that the scope of the pathway is defined first based on literature, and pathway components (proteins, genes, RNAs) are treated as ontology terms, or molecule classes, rather than specific instances. This means that multiple proteins from the same organism or different organisms can potentially play the same given role in a pathway. The advantage is that the work flow is more similar to the thinking process of the biologists who are the users of our curation software module. The third major property is that the curation software is designed to be simple enough to be used directly by bench biologists after a brief training course. All other pathway databases we are aware of employ highly trained curators, who of course cannot be experts in all areas of biology. The current set of PANTHER pathways has been curated by more than 40 different external experts from the scientific community; they must only have demonstrated their expertise with publications in the relevant field.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'reaction network',
        'omnipath': False,
        'emails': [('feedback@pantherdb.org', 'Panther Team'), ('paul.thomas@sri.com', 'Paul Thomas')]
    },
    'PhosphoSite': {
        'year': 2015,
        'releases': [2012, 2015],
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
                'http://nar.oxfordjournals.org/content/40/D1/D261.long',
                'http://nar.oxfordjournals.org/content/43/D1/D512.long'
            ],
            'webpages': [
                'http://www.phosphosite.org/homeAction.do'
            ]
        },
        'pubmeds': [22135298, 25514926],
        'taxons': [
            'human',
            'mouse',
            'eubacteria',
            'eukarya'
        ],
        'descriptions': [
            u'''
            PSP integrates both low- and high-throughput (LTP and HTP) data sources into a single reliable and comprehensive resource. Nearly 10,000 journal articles , including both LTP and HTP reports, have been manually curated by expert scientists from over 480 different journals since 2001.
            Information from nearly 13 000 papers and 600 different journals characterizing modification sites with LTP methods has been curated into PSP.
            Information is gathered from published literature and other sources. Published literature is searched semi-automatically with multiple intelligent search algorithms to identify reports that potentially identify phosphorylation sites in human, mouse or other species. Each identified report is then scanned by our highly trained curatorial staff (all with PhDs and extensive research experience in cell biology or related disciplines) to select only those papers that either identify new physiological phosphorylation sites or those that illuminate the biological function of the phosphorylation event. Records that are selected for inclusion into PhosphoSite are placed in the curatorial queue for processing. Note: while we gather records that describe both in vitro and in vivo phosphorylation events, we only finally submit records about in vitro sites when we have additional hard evidence that the site is also phosphorylated in vivo.
            '''
        ],
        'type': 'literature curated and high throughput',
        'subtype': 'post-translational modification',
        'omnipath': True,
        'emails': [('phornbeck@cellsignal.com', 'Paul Hornbeck'), ('EditorPhosphoSite@cellsignal.com', 'PhosphoSite Team')]
    },
    'SPIKE': {
        'year': 2012,
        'releases': [2008, 2011, 2012],
        'urls': {
            'articles': [
                'http://www.biomedcentral.com/1471-2105/9/110',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D793.full.html'
            ],
            'webpages': [
                'http://www.cs.tau.ac.il/~spike/'
            ]
        },
        'authors': ['Shamir Group', 'Shiloh Group'],
        'pubmeds': [18289391, 21097778],
        'descriptions': [
            u'''
            SPIKE’s data on relationships between entities come from three sources: (i) Highly curated data submitted directly to SPIKE database by SPIKE curators and experts in various biomedical domains. (ii) Data imported from external signaling pathway databaes. At present, SPIKE database imports such data from Reactome, KEGG, NetPath and The Transcription Factor Encyclopedia (http://www.cisreg.ca/cgi-bin/tfe/home.pl). (iii) Data on protein–protein interactions (PPIs) imported either directly from wide-scale studies that recorded such interactions [to date,PPI data were imported from Stelzl et al., Rual et al. and Lim et al.] or from external PPI databases [IntAct and MINT]. Relationship data coming from these different sources vary greatly in their quality and this is reflected by a quality level attribute, which is attached to each relationship in SPIKE database (Supplementary Data). Each relationship in SPIKE is linked to at least one PubMed reference that supports it.
            As of August 2010, the SPIKE database contains 20 412 genes/proteins, 542 complexes (327 of high quality), 320 protein families (167 of high quality) and 39 small molecules. These entities are linked by 34 338 interactions (of which 2400 are of high quality) and 6074 regulations (4420 of high quality). These are associated with 5873 journal references in total.
            Each of the maps is constructed by a domain expert; typically the same expert will also be responsible later for keeping it up-to-date. The expert reads the relevant literature and identifies those interactions and regulations that are pertinent to the pathway.
            The regulations and interactions in the database are assigned quality values ranging from 1 to 4. In general, relationships (regulations and interactions) derived from highly focused biochemical studies are assigned high quality (2 or 1) while those derived from high-throughput experiments are assigned lower quality (4 or 3). The curator uses best judgment to assign a quality level. For example, relationships mentioned in two independent research reports, or cited repeatedly in reviews written by leading authorities will get quality 1. Relationships with cited concrete references and those imported en masse from external curated signaling DBs are initially assigned quality 2 but later can be changed to the highest quality after the curator has read and was convinced by the cited papers. Data imported from protein-protein interaction DBs and datasets are assigned quality 3 or 4, depending on the experimental technique.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'pathway',
        'omnipath': True,
        'emails': [('rshamir@tau.ac.il', 'Ron Shamir')]
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
        'full_name': 'Pathway Interaction Database',
        'color': '',
        'urls': {
            'webpages': [
                'http://pid.nci.nih.gov/index.shtml'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/37/suppl_1/D674.long'
            ]
        },
        'pubmeds': [18832364],
        'data_import': ['BioCarta', 'Reactome'],
        'contains': [
            'BioCarta',
            'Reactome'
        ],
        'taxons': [
            'human'
        ],
        'descriptions': [
            u'''
            In curating, editors synthesize meaningful networks of events into defined pathways and adhere to the PID data model for consistency in data representation: molecules and biological processes are annotated with standardized names and unambiguous identifiers; and signaling and regulatory events are annotated with evidence codes and references. To ensure accurate data representation, editors assemble pathways from data that is principally derived from primary research publications. The majority of data in PID is human; however, if a finding discovered in another mammal is also deemed to occur in humans, editors may decide to include this finding, but will also record that the evidence was inferred from another species. Prior to publication, all pathways are reviewed by one or more experts in a field for accuracy and completeness.
            '''
        ],
        'notes': [
            u'''
            From the NCI-XML interactions with references, directions and signs can be extracted. Complexes are ommited.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'reaction network',
        'omnipath': True,
        'emails': [('yanch@mail.nih.gov', 'Chunhua Yan')]
    },
    'WikiPathways': {
        'urls': {
            'webpages': [
                'http://www.wikipathways.org/index.php/WikiPathways'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/40/D1/D1301'
            ]
        },
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
        'subtype': 'reaction network',
        'omnipath': False,
        'emails': [('thomaskelder@gmail.com', 'Thomas Kelder'), ('apico@gladstone.ucsf.edu', 'Alex Pico')]
    },
    'ConsensusPathDB': {
        'year': 2015,
        'releases': [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015],
        'urls': {
            'webpages': [
                'http://cpdb.molgen.mpg.de/CPDB'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/37/suppl_1/D623.long',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D712.long',
                'http://nar.oxfordjournals.org/content/41/D1/D793.long'
            ]
        },
        'taxons': [
            'human',
            'mouse',
            'yeast'
        ],
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
        'subtype': 'pathway',
        'emails': [('kamburov@molgen.mpg.de', 'Atanas Kamburov')]
    },
    'KEGG': {
        'year': 2015,
        'releases': [2000, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015],
        'urls': {
            'webpages': [
                'http://www.genome.jp/kegg/'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/28/1/27.long'
            ]
        },
        'notes': [
            u'''
            From 2011, KEGG data is not freely available. The downloadable KGML files contain binary interactions, most of them between large complexes. No references available.
            '''
        ],
        'type': 'literature curated',
        'subtype': 'reaction network',
        'omnipath': False,
        'emails': [('kanehisa@kuicr.kyoto-u.ac.jp', 'Minoru Kaneshia')]
    },
    'BioGRID': {
        'year': 2015,
        'releases': [2003, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015],
        'label': 'BioGRID',
        'authors': ['Tyers Lab'],
        'urls': {
            'webpages': [
                'http://thebiogrid.org/'
            ],
            'articles': [
                'http://genomebiology.biomedcentral.com/articles/10.1186/gb-2003-4-3-r23',
                'http://nar.oxfordjournals.org/content/34/suppl_1/D535.long',
                'http://nar.oxfordjournals.org/content/36/suppl_1/D637.long',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D698.long',
                'http://nar.oxfordjournals.org/content/41/D1/D816.long',
                'http://nar.oxfordjournals.org/content/43/D1/D470.long'
            ]
        },
        'pubmeds': [25428363, 23203989, 21071413, 18000002, 16381927, 12620108],
        'type': 'high throughput',
        'subtype': 'interaction',
        'omnipath': False,
        'emails': [('biogridadmin@gmail.com', 'BioGRID Team'), ('md.tyers@umontreal.ca', 'Michael Tyers')]
    },
    'STRING': {
        'year': 2015,
        'releases': [2015, 2013, 2011, 2009, 2007, 2005, 2003, 2000],
        'urls': {
            'webpages': [
                'http://string-db.org/'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/43/D1/D447.long',
                'http://nar.oxfordjournals.org/content/41/D1/D808.long',
                'http://nar.oxfordjournals.org/content/39/suppl_1/D561.long',
                'http://nar.oxfordjournals.org/content/37/suppl_1/D412.long',
                'http://nar.oxfordjournals.org/content/35/suppl_1/D358.long',
                'http://nar.oxfordjournals.org/content/33/suppl_1/D433.long',
                'http://nar.oxfordjournals.org/content/31/1/258.long',
                'http://nar.oxfordjournals.org/content/28/18/3442.long'
            ]
        },
        'pubmeds': [25352553, 23203871, 21045058, 18940858, 17098935,
            15608232, 12519996, 10982861],
        'authors': ['Bork Lab'],
        'label': 'STRING',
        'type': 'high-throughput and prediction',
        'subtype': 'interaction',
        'omnipath': False,
        'emails': [('bork@embl.de', 'Peer Bork'), ('lars.juhl.jensen@cpr.ku.dk', 'Lars Juhl Jensen'), 
            ('mering@imls.uzh.ch', 'Christian von Mering')]
    },
    'MINT': {
        'year': 2015,
        'releases': [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015],
        'label': 'MINT',
        'urls': {
            'webpages': [
                'http://mint.bio.uniroma2.it/mint/Welcome.do'
            ],
            'articles': [
                'http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1751541/'
            ]
        },
        'type': 'literature curated and high-throughput',
        'subtype': 'interaction',
        'omnipath': False,
        'emails': [('livia.perfetto@live.it', 'Livia Perfetto')]
    },
    'IntAct': {
        'year': 2015,
        'releases': [2003, 2006, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['EBI'],
        'label': 'IntAct',
        'color': '',
        'data_import': ['InnateDB', 'MINT'],
        'urls': {
            'articles': [
                'http://nar.oxfordjournals.org/content/42/D1/D358.long',
                'http://nar.oxfordjournals.org/content/40/D1/D841.long',
                'http://nar.oxfordjournals.org/content/38/suppl_1/D525.long',
                'http://nar.oxfordjournals.org/content/35/suppl_1/D561.long',
                'http://nar.oxfordjournals.org/content/32/suppl_1/D452.long'
            ],
            'webpages': [
                'http://www.ebi.ac.uk/intact/'
            ]
        },
        'pubmeds': [14681455, 17145710, 19850723, 22121220, 24234451],
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
        'omnipath': False,
        'emails': [('orchard@ebi.ac.uk', 'Sandra Orchard'), ('hhe@ebi.ac.uk', 'Henning Hermjakob')]
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
            'webpages': [
                'http://matrixdb.ibcp.fr/'
            ]
        },
        'pubmeds': [19147664, 20852260, 25378329],
        'taxons': [
            'mammalia'
        ],
        'descriptions': [
            u'''
            Protein data were imported from the UniProtKB/Swiss-Prot database (Bairoch et al., 2005) and identified by UniProtKB/SwissProt accession numbers. In order to list all the partners of a protein, interactions are associated by default to the accession number of the human protein. The actual source species used in experiments is indicated in the page reporting interaction data. Intracellular and membrane proteins were included to obtain a comprehensive network of the partners of extracellular molecules. Indeed, ECM proteins and GAGs bind to a number of membrane proteins or cell-associated proteoglycans and some of them interact with intracellular partners upon internalization (Dixelius et al., 2000). ECM proteins were identified by the UniProtKB/Swiss-Prot keyword ‘extracellular matrix’ and by the GO terms ‘extracellular matrix’, ‘proteinaceous extracellular matrix’ and their child terms. The proteins annotated with the GO terms ‘extracellular region’ and ‘extracellular space’, which are used for proteins found in biological fluids, were not included because circulating molecules do not directly contribute to the extracellular scaffold. Additionally, 96 proteins were manually (re-)annotated through literature curation. MatrixDB integrates 1378 interactions from the Human Protein Reference Database (HPRD, Prasad et al., 2009), 211 interactions from the Molecular INTeraction database (MINT, Chatr-Aryamontri et al., 2007), 46 interactions from the Database of Interacting Proteins (DIP, Salwinski et al., 2004), 232 interactions from IntAct (Kerrien et al., 2007a) and 839 from BioGRID (Breitkreutz et al., 2008) involving at least one extracellular biomolecule of mammalian origin. We added 283 interactions from manual literature curation and 65 interactions from protein and GAG array experiments.
            ''',
            u'''
            Interaction data stored in MatrixDB are (i) experimentally determined in the laboratory using surface plasmon resonance (SPR) binding assays, including protein and glycosaminoglycan arrays probed by SPR imaging, (ii) extracted from the literature by manual curation and (iii) imported from other interaction databases belonging to the IMEx consortium [IntAct, DIP, MINT, BioGRID], as well as from the Human Protein Reference Database. Imported data are restricted to interactions involving at least one extracellular protein.
            ''',
            u'''
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
        'omnipath': True,
        'emails': [('matrixdb@ibcp.fr', 'MatrixDB Team'), ('sylvie.ricard-blum@ibcp.fr', 'Sylvie Ricard-Blum')]
    }
}


def gen_html():
    head = u'''<!DOCTYPE html>
    <html lang="en">
        <head>
            <meta http-equiv="X-UA-Compatible" content="IE=edge">
            <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
            <meta name="dc.language" content="en">
            <meta name="viewport" content="width=device-width, 
            initial-scale=1.0">
            <style type="text/css">
                @font-face{
                    font-family:'HelveticaNeueLT Pro';
                    src:url('http://www.ebi.ac.uk/web_guidelines/fonts/helveticaneueltprolt.eot');
                    src:url('http://www.ebi.ac.uk/web_guidelines/fonts/helveticaneueltprolt.eot?#iefix') format('embedded-opentype'),
                        url('http://www.ebi.ac.uk/web_guidelines/fonts/helveticaneueltprolt.otf') format('opentype'),
                        local('☺'),
                        url('http://www.ebi.ac.uk/web_guidelines/fonts/helveticaneueltprolt.woff') format('woff'),
                        url('http://www.ebi.ac.uk/web_guidelines/fonts/helveticaneueltprolt.ttf') format('truetype'),
                        url('http://www.ebi.ac.uk/web_guidelines/fonts/helveticaneueltprolt.svg#HelveticaNeueLTPro-Lt') format('svg');
                    font-weight:normal;
                    font-style:normal
                }
                .quotebox {
                    border-color: #6FA6A9;
                    border-style: solid;
                    border-width: 1px 1px 1px 3px;
                    margin: 0.75em 1em;
                    padding: 0px 0.75em;
                    background-color: #E4ECEC;
                    color: #444444;
                }
                body {
                    margin: 6ex 23.61% 18ex 14.59%;
                    width: 61.8%;
                    color: #646567;
                    font-family: "HelveticaNeueLT Pro", "Helvetica Neue", Helvetica, Arial, sans-serif;
                }
                a:link, a:active, a:visited {
                    color: #646567;
                    text-decoration: dotted;
                }
                h1 {
                    color: #6FA6A9;
                }
                h2.omnipath, a.omnipath:link, a.omnipath:active, a.omnipath:visited {
                    color: #6EA945;
                    text-decoration: none;
                }
                h2.base, a.base:link, a.base:active, a.base:visited {
                    color: #646567;
                    text-decoration: none;
                }
            </style>
            <title>Signaling Pathway Resources</title>
        </head>
        <body>
        '''
    foot = u'''
        <br>
        <br>
        <p> <a href="http://www.ebi.ac.uk/~denes">Denes Turei, 2015.</a> 
            Comments: denes@ebi.ac.uk </p>
        <p>
            <a href="https://validator.w3.org/check/referer">
                <img 
                    src="data:image/png;base64,
                        iVBORw0KGgoAAAANSUhEUgAAAFAAAAAPCAIAAAD8q9/YAAAA6ElEQVRIieVXUQqEIBAdl74D57cD
                        BB6mY3iojuFhgg7Qrx+eYD/eIjJF2rYVrI9B5KXjPGdUUlprqgkNEXnvnw7jJjBzg94y9OfddW4O
                        1pz3cxHacSKi19Nh3I3qBXduXvc7N8PESMEUAqUlOjDBbw4WUwSfxXaGhZJl6GHpFgjmDHDy0/MP
                        GftigjWwQ2s1heMgbH23fXfbZdMSrGnHCe2+k6sEQxgSm/JrpgSIMpvAkq/Yl/KlpWBUaWxBxrqN
                        TCzmn7xnhUhTGjfraIaV1tp7X8k7zMzVP0t/j09JPx3GTWBmVdvf0hvVzH8tNHhSIAAAAABJRU5E
                        rkJggg=="
                     alt="Valid HTML5" />
            </a>
            <a href="https://jigsaw.w3.org/css-validator/check/referer">
                <img 
                    src="data:image/png;base64,
                        iVBORw0KGgoAAAANSUhEUgAAAFAAAAAPCAIAAAD8q9/YAAAA7ElEQVRIieWX0Q2EIAyGC2EBuoAD
                        MIcrOIejOIcrOIcDuAAjmHvopWlaLycS9YHvwdQ/hPSnrUQXY4SWCACQc347jYdAxECRH5f67fap
                        T/NWv89NrEMHAP7tNJ6mOcNBve9Tz+3N8T71pMjOJ7FmFqjH5BSQIkWl2AWlaMOEtA3ClTwCFZSy
                        Dh17UIFdw8iDuOb52LDlVz3rv3acd5o3VXOrwFFfFHF2hv24+HHh3masUgp3KQCkeWOT55UidIXJ
                        FT9JtDPM5i9X2Fbv7wAfKqW4GGPOuZF7GBGbu5aaM/xt6bfTeAhEdK39LX0ALc+PofdKE1EAAAAA
                        SUVORK5CYII="
                     alt="Valid CSS3" />
            </a>
        </p>
        \t</body>\n
    </html>\n'''
    doc = u'''<img src="data:image/png;base64,
        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" alt="EMBL-EBI"/>'''
    doc += '\t<h1>Metadata about signaling pathway resources</h1>\n'
    doc += '\t<p>This collection was created during the construction '\
        'of OmniPath, a network of signaling pathways intending to '\
        'merge all high quality, manually curated efforts. The '\
        'descriptions here aim to cite the relevant sentences '\
        'about the curation protocols from the original articles. '\
        'URLs pointing to the articles and the webpages, and some '\
        'additional metadata is provided where available. Note: '\
        'the green resources are included by default in '\
        'OmniPath.</p>\n'
    doc += '\t<h2>Contents</h2>\n'
    doc += '\t<ul>\n'
    for k, v in sorted(descriptions.items(), key = lambda x: x[0].lower()):
        doc += '\t\t\t<li><a href="#%s" class="%s">%s</a></li>\n' % \
            (k, 'omnipath' if 'omnipath' in v and v['omnipath'] else 'base', 
            v['full_name'] if 'full_name' in v else v['label'] if 'label' in v else k)
    doc += '\t</ul>\n'
    for k, v in sorted(descriptions.items(), key = lambda x: x[0].lower()):
        doc += '\t\t<br>\n\t\t<h2 id="%s" class="%s">%s</h2>\n' % \
            (k, 'omnipath' if 'omnipath' in v and v['omnipath'] else 'base',
            v['full_name'] if 'full_name' in v else v['label'] if 'label' in v else k)
        doc += '\t\t\t<p><b>Category || Subcategory &gt;&gt;&gt;</b> %s || %s</p>\n' % \
            (v['type'].capitalize() if 'type' in v else 'Undefined', 
                v['subtype'].capitalize() if 'subtype' in v else 'Undefined')
        if 'year' in v:
            doc += '\t\t\t<h3>Last updated: %u<\h3>\n'%v['year']
        if 'releases' in v:
            doc += '\t\t\t<p><b>Updated in years: </b>%s</p>\n' % \
                ', '.join(['%u'%y for y in v['releases']])
        if 'authors' in v and v['authors'] is not None:
            doc += '\t\t\t<p><b>Created by </b>%s</p>\n' % ', '.join(v['authors'])
        if 'emails' in v and v['emails'] is not None:
            doc += '\t\t\t<p><b>Contact: </b></p>\n\n\t\t\t\t<ul>\n%s\n' % \
                ''.join(['\t\t\t\t<li><a href="mailto:%s">%s &lt;%s&gt;</li>\n' % (em[0], em[1], em[0]) for em in v['emails']])
            doc += '\t\t\t\t</ul>\n'
        for uk, uv in v['urls'].iteritems():
            if len(uv) > 0:
                doc += '\t\t\t<h3>%s</h3>\n' % (uk.capitalize())
                doc += '\t\t\t<ul>\n'
                for a in uv:
                    doc += '\t\t\t\t<li><a href="%s" target="_blank">%s</a></li>\n' % (
                        a,a)
                doc += '\t\t\t</ul>\n'
        if 'pubmeds' in v:
            doc += '\t\t\t<h3>PubMed</h3>\n'
            doc += '\t\t\t<ul>\n'
            for pmid in v['pubmeds']:
                doc += '\t\t\t\t<li><a href="%s" target="_blank">%s</a></li>\n' % (
                    'http://www.ncbi.nlm.nih.gov/pubmed/%u' % pmid, 
                    'http://www.ncbi.nlm.nih.gov/pubmed/%u' % pmid
                    )
            doc += '\t\t\t</ul>\n'
        if 'taxons' in v:
            doc += '<p><b>Taxons: </b><em>%s</em></p>' % \
                ', '.join(['%s%s'%(t[0].upper(), t[1:]) for t in v['taxons']])
        if 'size' in v and type(v['size']) is dict:
            doc += '<p><b>Nodes: </b>%s, <b>Edges:</b>%s</p>' % (
                v['size']['nodes'],v['size']['edges'])
        if 'data_import' in v:
            doc += '\t\t\t<p><b>Direct data import from: </b>%s</p>\n' % \
                ', '.join(v['data_import'])
        if 'includes' in v:
            doc += '\t\t\t<p><b>Includes data from: </b>%s</p>\n' % \
                ', '.join(v['includes'])
        doc += '\t\t\t<h3>Quotes</h3>\n'
        if 'descriptions' in v:
            doc += '\t\t\t\t<div class="quotebox">\n'
            pars = v['descriptions'][0].split('\n')
            for p in pars:
                p = p.strip()
                if len(p) > 0:
                    doc += '\t\t\t\t<p>%s</p>\n' % p
            doc += '\t\t\t\t</div>\n'
        if 'notes' in v:
            doc += '\t\t\t\t<div class="quotebox">\n'
            pars = v['notes'][0].split('\n')
            for p in pars:
                p = p.strip()
                if len(p) > 0:
                    doc += '\t\t\t\t<p>%s</p>\n' % p
            doc += '\t\t\t\t</div>\n'
    html = head + doc + foot
    soup = bs4.BeautifulSoup(html, from_encoding = 'utf-8')
    return soup.prettify()

def write_html(filename = 'resources.html'):
    html = gen_html()
    with codecs.open(filename, encoding = 'utf-8', mode = 'w') as f:
        f.write(html)