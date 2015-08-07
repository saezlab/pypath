#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
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

data = {
    'CancerCellMap': {
        'year': 2006,
        'size': 0,
        'authors': ['Bader Lab', 'Sander Group'],
        'label': 'Cancer Cell Map',
        'color': ''
    },
    'MatrixDB': {
        'year': None,
        'size': None,
        'authors': None,
        'label': None,
        'color': None
    },
    'SPIKE': {
        'year': 2012,
        'releases': [2008, 2011],
        'size': 0,
        'authors': [''],
        'label': 'SPIKE',
        'color': '',
        'data_import': ['KEGG', 'Reactome', 'HPRD', 'IntAct', 'NetPath', 'MINT', 'TFE', 'Rual2005', 'Ataxia']
    },
    'DIP': {
        'year': 2014,
        'size': 0,
        'authors': ['UCLA', 'Eisenberg Group'],
        'label': 'DIP',
        'color': ''
    },
    'HPRD': {
        'year': None,
        'size': None,
        'authors': None,
        'label': None,
        'color': None
    },
    'PDZBase': {
        'year': None,
        'size': None,
        'authors': None,
        'label': None,
        'color': None
    },
    'dbPTM': {
        'year': None,
        'size': None,
        'authors': None,
        'label': None,
        'color': None
    },
    'InnateDB': {
        'year': 2014,
        'size': 0,
        'authors': ['Brinkman Lab', 'Hancock Lab'],
        'label': 'InnateDB',
        'color': ''
    },
    'ACSN': {
        'year': 2015,
        'size': 0,
        'authors': ['Curie'],
        'label': 'ACSN',
        'color': ''
    },
    'DOMINO': {
        'year': 2006,
        'releases': [2006],
        'size': 0,
        'authors': ['Cesareni Group'],
        'label': 'DOMINO',
        'color': ''
    },
    'Signor': {
        'year': 2015,
        'size': 0,
        'authors': ['Cesareni Group'],
        'label': 'Signor',
        'color': '',
        'data_import': ['SignaLink2', 'PhosphoSite']
    },
    'Macrophage': {
        'year': 2010,
        'releases': [2006, 2010],
        'size': 0,
        'authors': ['Raza', 'Oda'],
        'label': 'Macrophage',
        'color': ''
    },
    'NetPath': {
        'year': 2010,
        'size': 0,
        'authors': ['Pandey Lab', 'IOB Bangalore'],
        'label': 'NetPath',
        'color': '',
        'data_import': ['CancerCellMap']
    },
    'ELM': {
        'year': 2014,
        'releases': [2003, 2008, 2009, 2012, 2013, 2014],
        'size': 0,
        'authors': 'ELM Consortium',
        'label': 'ELM',
        'color': ''
    },
    'SignaLink2': {
        'year': 2012,
        'size': 0,
        'authors': ['NetBiol Group'],
        'label': 'SignaLink',
        'color': ''
    },
    'NRF2ome': {
        'year': 2012,
        'size': 0,
        'authors': ['NetBiol Group'],
        'label': 'NRF2ome',
        'color': ''
    },
    'DEPOD': {
        'year': 2014,
        'size': 0,
        'authors': ['EMBL & EBI'],
        'label': 'DEPOD',
        'color': ''
    },
    'MPPI': {
        'year': 2005,
        'size': 0,
        'authors': ['MIPS Munich'],
        'label': 'MPPI',
        'color': ''
    },
    'Guide2Pharmacology': {
        'year': None,
        'size': None,
        'authors': None,
        'label': None,
        'color': None
    },
    'IntAct': {
        'year': 2015,
        'size': 0,
        'authors': ['EBI'],
        'label': None,
        'color': None,
        'data_import': ['InnateDB']
    },
    'AlzPathway': {
        'year': 2012,
        'size': 0,
        'authors': ['Tokyo Bioinf'],
        'label': 'AlzPathway',
        'color': ''
    },
    'PhosphoSite': {
        'year': 2015,
        'size': 0,
        'authors': ['CST'],
        'label': 'PhosphoSite',
        'color': ''
    },
    'NCI-PID': {
        'year': 2014,
        'size': 0,
        'authors': ['NCI'],
        'label': 'NCI-PID',
        'color': ''
    },
    'DeathDomain': {
        'year': 2012,
        'releases': [2011, 2012],
        'size': 0,
        'authors': [''],
        'label': 'DeathDomain',
        'color': ''
    },
    'ARN': {
        'year': 2014,
        'releases': [2014],
        'size': 0,
        'authors': ['NetBiol Group'],
        'label': 'ARN',
        'color': ''
    },
    'BioCarta': {
        'year': 2006,
        'releases': [2006],
        'size': 0,
        'authors': ['Community'],
        'label': 'BioCarta',
        'color': ''
    },
    'CA1': {
        'year': 2005,
        'releases': [2005],
        'size': 0,
        'authors': ['Ma\'ayan Lab'],
        'label': 'Ma\'ayan 2005',
        'color': ''
    },
    'Cui2007': {
        'year': 2007,
        'size': 0,
        'authors': ['Wang Group'],
        'label': 'Cui 2007',
        'color': '',
        'data_import': ['Awan2007', 'CancerCellMap']
    },
    'Awan2007': {
        'year': 2007,
        'size': 0,
        'authors': ['Wang Group'],
        'label': 'Awan 2007',
        'color': '',
        'data_import': ['BioCarta', 'CA1']
    },
    'CST': {
        'year': None,
        'size': 0,
        'authors': ['CST'],
        'label': 'CST pathways',
        'color': ''
    },
    'HSN': {
        'year': 2014,
        'size': 0,
        'authors': ['Wang Group'],
        'label': 'HumanSignalingNetwork',
        'color': '',
        'data_import': ['Cui2007', 'BioCarta', 'CST', 'NCI-PID']
    }
}


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
        ]
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
                'http://nar.oxfordjournals.org/content/40/D1/D242.long'
                'http://nar.oxfordjournals.org/content/42/D1/D259.long'
            ]
        },
        'pubmeds': [22110040, 24214962]
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
        ]
    },
    'SignaLink2': {
        'year': 2012,
        'releases': [2010, 2012],
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
            SignaLink assigns proteins to signaling pathways using the full texts of pathway reviews (written by pathway experts). While most signaling resources consider 5–15 reviews per pathway, SignaLink uses a total of 170 review papers, i.e. more than 20 per pathway on average. Interactions were curated from a total of 941 articles (PubMed IDs are available at the website). We added a small number of proteins based on InParanoid ortholog clusters. For curation, we used a self-developed graphical tool and Perl/Python scripts. The current version of SignaLink was completed in May 2008 based on WormBase (version 191), FlyBase (2008.6), Ensembl (49), UniProt (87) and the publications listed on the website. 
            The curation protocol of SignaLink (Fig. 1A) contains several steps aimed specifically at reducing data and curation errors. We used reviews as a starting point, manually looked up interactions three times, and manually searched for interactions of known signaling proteins with no signaling interactions so far in the database. 
            '''
        ]
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
        ]
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
        ]
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
        ]
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
        ]
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
            The curation process follows the PSI-MI 2.5 standard but with special emphasis on the mapping of the interaction to specific protein domains of both participating proteins. This is achieved by paying special attention to the shortest protein fragment that was experimentally verified as sufficient for the interaction. Whenever the authors report only the name of the domain mediating the interaction (i.e. SH3, SH2 ...), without stating the coordinates of the experimental binding range, the curator may choose to enter the coordinates of the Pfam domain match in the protein sequence. Finally whenever the information is available, any mutation or post- translational modification affecting the interaction affinity is noted in the database.
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
        ]
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
        ]
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
        ]
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
        'pubmeds': [17907678],
        'descriptions': [
            u'''
            To construct the human cellular signalling network, we manually curated signalling pathways from literature. The signalling data source for our pathways is the BioCarta database (http://www.biocarta.com/genes/allpathways.asp), which, so far, is the most comprehensive database for human cellular signalling pathways. Our curated pathway database recorded gene names and functions, cellular locations of each gene and relationships between genes such as activation, inhibition, translocation, enzyme digestion, gene transcription and translation, signal stimulation and so on. To ensure the accuracy and the consistency of the database, each referenced pathway was cross-checked by different researchers and finally all the documented pathways were checked by one researcher. In total, 164 signalling pathways were documented (supplementary Table 2). Furthermore, we merged the curated data with another literature-mined human cellular signalling network. As a result, the merged network contains nearly 1100 proteins (SupplementaryNetworkFile). To construct a signalling network, we considered relationships of proteins as links (activation or inactivation as directed links and physical interactions in protein complexes as neutral links) and proteins as nodes.
            '''
        ]
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
        ]
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
            This resource includes a huge number of pathways, each curated by experts form a few reviews. The data is not available for download from the original webpage, only from second hand, for example from NCI-PID, in NCI-XML format. However, these files doesn't contain any references, which makes problematic the use of the BioCarta dataset. Also, some pathways are reviewed long time ago, possibly outdated.
            '''
        ]
    },
    'TLR': {
        'urls': {
            'articles': [
                'http://msb.embopress.org/content/2/1/2006.0015.long'
            ]
        }
    },
    'CA1': {
        'year': 2005,
        'releases': [2005],
        'size': {
            'nodes': 545,
            'edges': 1259
        },
        'authors': ['Ma\'ayan Lab'],
        'label': 'Ma\'ayan 2005',
        'color': '',
        'pubmeds': [16099987],
        'urls': {
            'articles': [
                'http://www.sciencemag.org/content/309/5737/1078.full'
            ],
            'webpages': []
        },
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
        ]
    },
    'CancerCellMap': {
        'urls': {
            'articles': [],
            'webpages': [
                'http://www.pathwaycommons.org/pc-snapshot/current-release/tab_delim_network/by_source/'
            ]
        },
        'descriptions': [
            u'''
            Manually curated data, unpublished. A team of M.Sc. and Ph.D. biologists at the Institute of Bioinformatics in Bangalore, India read original research papers and hand-entered the pathway data into our database. The quality of the Cancer Cell Map pathways is very high. Half of the pathways were reviewed by experts at Memorial Sloan-Kettering Cancer Center and were found to contain only a few errors, which were subsequently fixed.
            '''
        ],
        'notes': [
            u'''
            One of the first manually curated datasets, now only available from second hand, e.g. from PathwayCommons. Included in many other resources. Contains binary interactions with PubMed references.
            '''
        ]
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
        ]
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
        'descriptions': [
            u'''
            In order to expand the interaction dataset, we added relevant direct protein–protein interactions from currently available human protein–protein interaction networks (Rual et al., 2005; Stelzl et al., 2005). We also searched public databases, including BIND (Bader et al., 2003), DIP (Xenarios et al., 2002), HPRD (Peri et al., 2003), MINT (Zanzoni et al., 2002), and MIPS (Pagel et al., 2005), to identify literature-based binary interactions involving the 54 ataxia-associated baits and the 561 interacting prey proteins. We identified 4796 binary protein–protein interactions for our Y2H baits and prey proteins (Table S4) and incorporated them in the Y2H protein–protein interaction map (Figures 4A–4C).
            '''
        ],
        'notes': [
            u'''
            The Ataxia network doesn't contain original manual curation effort. The integrated data are very old. 
            '''
        ]
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
        ]
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
        'taxons': [
            'human'
        ],
        'descriptions': [
            u'''
            Human phosphotyrosine signaling network. 
            We manually collected the experimentally determined human TK–substrate interactions and substrate–SH2/PTB domain interactions from the literature (see Supplemental Materials), as well as the Phospho.ELM and PhosphoSitePlus databases. [71 references, 585 circuits]
            '''
        ]
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
        ]
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
        'descriptions': [
            u'''
            We collected 123 review articles related to AD accessible from PubMed. We then manually curated these review articles, and have built an AD pathway map by using CellDesigner. Molecules are distinguished by the following types: proteins, complexes, simple molecules, genes, RNAs, ions, degraded products, and phenotypes. Gene symbols are pursuant to the HGNC symbols. Reactions are also distinguished by the following categories: state transition, transcription, translation, heterodimer association, dissociation, transport, unknown transition, and omitted transition. All the reactions have evidences to the references in PubMed ID using the MIRIAM scheme. All the references used for constructing the AlzPathway are listed in the ‘References for AlzPathway’. Cellular types are distinguished by the followings: neuron, astrocyte, and microglial cells. Cellular compartments are also distinguished by the followings: brain blood barrier, presynaptic, postsynaptic, and their inner cellular localizations.
            '''
        ],
        'notes': [
            u'''
            References can be fetched only from XML formats, not from the SIF file. Among approx. 150 protein-protein interactions, also contains interactions of many small molecules, denoted by pubchem IDs.
            '''
        ]
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
        ]
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
        'descriptions': [
            u'''
            Annotation of the manual dataset was performed analogous to the annotation of protein–protein interactions and protein complexes in previous projects published by our group. Information about NIPs was extracted from scientific literature using only data from individual experiments but not from high-throughput experiments. Only mammalian proteins were considered. Data from high-throughput experiments were omitted in order to maintain the highest possible standard of reliability.
            '''
        ]
    },
    'Macrophage': {
        'year': 2010,
        'urls': {
            'articles': [
                'http://www.biomedcentral.com/1752-0509/4/63'
            ],
            'webpages': []
        },
        'pubmeds': [20470404],
        'descriptions': [
            u'''
            Ongoing analysis of macrophage-related datasets and an interest in consolidating our knowledge of a number of signalling pathways directed our choice of pathways to be mapped (see Figure 1). Public and propriety databases were initially used as resources for data mining, but ultimately all molecular interaction data was sourced from published literature. Manual curation of the literature was performed to firstly evaluate the quality of the evidence supporting an interaction and secondly, to extract the necessary and additional pieces of information required to 'understand' the pathway and construct an interaction diagram. We have drawn pathways based on our desire to model pathways active in a human macrophage and therefore all components have been depicted using standard human gene nomenclature (HGNC). However, our understanding of the pathway components and the interactions between them, have been drawn largely from a consensus view of literature knowledge. As such the pathways presented here are based on data derived from a range of different cellular systems and mammalian species (human and mouse).
            '''
        ]
    },
    'NetPath': {
        'year': 2011,
        'releases': [2010, 2011],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['Pandey Lab', 'IOB Bangalore'],
        'label': 'NetPath',
        'color': '',
        'data_import': [
            'CancerCellMap'
        ],
        'includes': ['CancerCellMap'],
        'urls': {
            'articles': [
                'http://genomebiology.com/content/11/1/R3',
                'http://database.oxfordjournals.org/content/2011/bar032.long'
            ],
            'webpages': [
                'http://netpath.org/'
            ]
        },
        'pubmeds': [20067622, 21959865],
        'descriptions': [
            u'''
            The initial annotation process of any signaling pathway involves gathering and reading of review articles to achieve a brief overview of the pathway. This process is followed by listing all the molecules that arereported to be involved in the pathway under annotation. Information regarding potential pathway authorities are also gathered at this initial stage. Pathway experts are involved in initial screening of the molecules listed to check for any obvious omissions. In the second phase, annotators manually perform extensive literature searches using search keys, which include all the alter native names of the molecules involved, the name of the pathway, the names of reactions, and so on. In addition, the iHOP resource is also used to perform advanced PubMed-based literature searches to collect the reactions that were reported to be implicated in a given pathway. The collected reactions are manually entered using the PathBuilder annotation interface, which is subjected to an internal review process involving PhD level scientists with expertise in the areas of molecular biology, immunology and biochemistry. However, there are instances where a molecule has been implicated in a pathway in a published report but the associated experimental evidence is either weak or differs from experiments carried out by other groups. For this purpose, we recruit several investigators as pathway authorities based on their expertise in individual signaling pathways. The review by pathway authorities occasionally leads to correction of errors or, more commonly, to inclusion of additional information. Finally, the pathway authorities help in assessing whether the work of all major laboratories has been incorporated for the given signaling pathway.
            '''
        ],
        'notes': [
            u'''
            Formats are unclear. The tab delimited format contains the pathway memberships of genes, PubMed references, but not the interaction partners! The Excel file is very weird, in fact it is not an excel table, and contains only a few rows from the tab file. The PSI-MI XML is much better. By writing a simple parser, a lot of details can be extracted.
            '''
        ]
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
        'pubmeds': [
            20727158,
            23180781,
            18766178
        ],
        'releases': [2008, 2010, 2013],
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
        ]
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
        ]
    },
    'CST': {
        'year': 2015,
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
        'descriptions': [
            u'''
            On these resource pages you can find signaling pathway diagrams, research overviews, relevant antibody products, publications, and other research resources organized by topic. The pathway diagrams associated with these topics have been assembled by CST scientists and outside experts to provide succinct and current overviews of selected signaling pathways.
            '''
        ],
        'notes': [
            u'''
            The pathway diagrams are based on good quality, manually curated data, probably from review articles. However, those are available only in graphical (PDF and InDesign) formats. There is no programmatic way to obtain the interactions and references, as it was confirmed by the authors, who I contacted by mail. Wang's HumanSignalingNetwork includes the data from this resource, which probably has been entered manually, but Wang's data doesn't have source annotations, despite it's compiled from multiple sources.
            '''
        ]
    },
    'DIP': {
        'year': 2014,
        'releases': [2000, 2011],
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
        ]
    },
    'DEPOD': {
        'year': 2014,
        'releases': [2013, 2014],
        'size': {
            'nodes': None,
            'edges': None
        },
        'authors': ['EMBL & EBI'],
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
            DEPOD - the human DEPhOsphorylation Database (version 1.0) is a manually curated database collecting human active phosphatases, their experimentally verified protein and non-protein substrates and dephosphorylation site information, and pathways in which they are involved. It also provides links to popular kinase databases and protein-protein interaction databases for these phosphatases and substrates. DEPOD aims to be a valuable resource for studying human phosphatases and their substrate specificities and molecular mechanisms; phosphatase-targeted drug discovery and development; connecting phosphatases with kinases through their common substrates; completing the human phosphorylation/dephosphorylation network.
            '''
        ],
        'notes': [
            u'''
            Nice manually curated dataset with PubMed references, in easily accessible MITAB format with UniProt IDs, comprises 832 dephosphorylation reactions on protein substrates, and few hundreds on small molecules.
            '''
        ]
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
        ]
    },
    'PANTHER': {
        'year': 2009,
        'releases': [2006, 2009],
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
        ]
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
        ]
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
        'pubmeds': [18289391, 21097778],
        'descriptions': [
            u'''
            SPIKE’s data on relationships between entities come from three sources: (i) Highly curated data submitted directly to SPIKE database by SPIKE curators and experts in various biomedical domains. (ii) Data imported from external signaling pathway databaes. At present, SPIKE database imports such data from Reactome, KEGG, NetPath and The Transcription Factor Encyclopedia (http://www.cisreg.ca/cgi-bin/tfe/home.pl). (iii) Data on protein–protein interactions (PPIs) imported either directly from wide-scale studies that recorded such interactions [to date,PPI data were imported from Stelzl et al., Rual et al. and Lim et al.] or from external PPI databases [IntAct and MINT (19)]. Relationship data coming from these different sources vary greatly in their quality and this is reflected by a quality level attribute, which is attached to each relationship in SPIKE database (Supplementary Data). Each relationship in SPIKE is linked to at least one PubMed reference that supports it.
            As of August 2010, the SPIKE database contains 20 412 genes/proteins, 542 complexes (327 of high quality), 320 protein families (167 of high quality) and 39 small molecules. These entities are linked by 34 338 interactions (of which 2400 are of high quality) and 6074 regulations (4420 of high quality). These are associated with 5873 journal references in total.
            Each of the maps is constructed by a domain expert; typically the same expert will also be responsible later for keeping it up-to-date. The expert reads the relevant literature and identifies those interactions and regulations that are pertinent to the pathway.
            The regulations and interactions in the database are assigned quality values ranging from 1 to 4. In general, relationships (regulations and interactions) derived from highly focused biochemical studies are assigned high quality (2 or 1) while those derived from high-throughput experiments are assigned lower quality (4 or 3). The curator uses best judgment to assign a quality level. For example, relationships mentioned in two independent research reports, or cited repeatedly in reviews written by leading authorities will get quality 1. Relationships with cited concrete references and those imported en masse from external curated signaling DBs are initially assigned quality 2 but later can be changed to the highest quality after the curator has read and was convinced by the cited papers. Data imported from protein-protein interaction DBs and datasets are assigned quality 3 or 4, depending on the experimental technique.
            '''
        ]
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
        ]
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
        ]
    },
    'ConsensusPathDB': {
        'year': 2014,
        'releases': [2007, 2014],
        'urls': {
            'webpages': [
                'http://cpdb.molgen.mpg.de/CPDB'
            ],
            'articles': [
                'http://nar.oxfordjournals.org/content/41/D1/D793',
                
            ]
        },
        'taxons': [
            'human',
            'mouse',
            'yeast'
        ],
        'pubmeds': [18940869, 23143270],
        'descriptions': [
            '''
            
            '''
        ],
        'notes': [
            u'''
            ConsensusPathDB comprises data from 32 resources. The format is easy to use, tab delimited text file, with UniProtKB names and PubMed IDs. However, the dataset is extremely huge, and several databases containing HTP data is included.
            '''
        ]
    },
    'KEGG': {
        'urls': {
            'webpages': [
                'http://www.genome.jp/kegg/'
            ],
            'articles': [
                ''
            ]
        },
        'notes': [
            u'''
            From 2011, KEGG data is not freely available. The downloadable KGML files contain binary interactions, most of them between large complexes. No references available.
            '''
        ]
    },
    'IntAct': {
        'year': 2015,
        'releases': [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015],
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
            The information within the IntAct database primarily consists of protein–protein interaction (PPI) data. The majority of the PPI data within the database is annotated to IMEx standards, as agreed by the IMEx consortium. All such records contain a full description of the experimental conditions in which the interaction was observed. This includes full details of the constructs used in each experiment, such as the presence and position of tags, the minimal binding region defined by deletion mutants and the effect of any point mutations, referenced to UniProtKB (2), the underlying protein sequence database. Protein interactions can be described down to the isoform level, or indeed to the post-translationally cleaved mature peptide level if such information is available in the publication, using the appropriate UniProtKB identifiers.
            Each entry in IntAct is peer reviewed by a senior curator, and not released until accepted by that curator. Additional rule-based checks are run at the database level, and manually fixed when necessary. Finally, on release of the data, the original author of each publication is contacted and asked to comment on the representation of their data; again manual updates are made to the entry should the author highlight any errors.
            All binary interactions evidences in the IntAct database, including those generated by Spoke expansion of co-complex data, are clustered to produce a non-redundant set of protein pairs (R. C. Jimenez et al., manuscript in preparation). Each binary pair is then scored, using a simple addition of the cumulated value of a weighted score for the interaction detection method and the interaction type for each interaction evidence associated with that binary pair, as described using the PSI-MI CV terms. The scores are given in Table 1, all children of each given parent receives that score. Only experimental data is scored, inferred interactions, for example, would be excluded. Any low confidence data or data manually tagged by a curator for exclusion from the process, would not be scored. Isoforms and post-processed protein chains are regarded as distinct proteins for scoring purposes.
            '''
        ],
        'notes': [
            u'''
            We can not draw a sharp distinction between low and high throughput methods, and I can agree, that this is not the only and best measure of quality considering experimental data. I see that IntAct came up with a good solution to estimate the confidence of interactions. The mi-score system gives a comprehensive way to synthetize information from multiple experiments, and weight interactions according to experimental methods, interaction type, and number of evidences.
            '''
        ]
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
            Interaction data stored in MatrixDB are (i) experimentally determined in the laboratory using surface plasmon resonance (SPR) binding assays, including protein and glycosaminoglycan arrays probed by SPR imaging (11), (ii) extracted from the literature by manual curation and (iii) imported from other interaction databases belonging to the IMEx consortium [IntAct (12), DIP (13), MINT (14), BioGRID (15)], as well as from the Human Protein Reference Database (16). Imported data are restricted to interactions involving at least one extracellular protein.
            ''',
            u'''
            The content of MatrixDB has been updated with new interaction data manually curated by the MatrixDB team, and by importing interaction data from four interaction databases of the IMEx consortium via The Proteomics Standard Initiative Common QUery InterfaCe (PSICQUIC), a community standard for computational access to molecular-interaction data resources (25).  In the current release MatrixDB contains 904 interactions supported by 1244 experiments, which have been manually curated from 237 publications, compared to 490 interactions supported by 847 experiments in the previous version of the database (2). This is the MatrixDB ‘core’ data set.
            '''
        ],
        'notes': [
            u'''
            Very nice! Note: The interactions imported from IMEX databases or any other database, are collected separately, in the PSICQUIC-extended dataset. The MatrixDB-core dataset is curated manually by the MatrixDB team.
            '''
        ]
    }
}


def gen_html():
    head = u'''<!DOCTYPE html>
    <html>
        <head>
            <meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1">
            <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
            <meta http-equiv="Content-Language" content="en">
            <meta name="dc.language" content="en">
            <meta name="viewport" content="width=device-width, 
            initial-scale=1.0">
            <style type="text/css">
                .quotebox {
                    border-color: #7AADBD;
                    border-style: solid;
                    border-width: 1px 1px 1px 3px;
                    margin: 0.75em 1em;
                    padding: 0px 0.75em;
                    background-color: #F6F9FC;
                    color: #444444;
                }
                body {
                    margin: 6ex 23.61% 18ex 14.59%;
                    width: 61.8%;
                }
            </style>
        </head>
        <body>
        '''
    foot = u'''\t</body>\n
    </html>\n'''
    doc = u''
    doc += '\t<h1>Metadata about signaling pathway resources</h1>\n'
    doc += '\t<p>This collection was created during the construction '\
        'of OmniPath, a network of signaling pathways intending to '\
        'merge all high quality, manually curated efforts. The '\
        'descriptions here aim to cite the relevant sentences '\
        'about the curation protocols from the original articles. '\
        'URLs pointing to the articles and the webpages, and some '\
        'additional metadata is provided where available. Note: '\
        'not all resources listed here were finally included in '\
        'OmniPath.</p>\n'
    doc += '\t<h2>Contents</h2>\n'
    doc += '\t<ul>\n'
    for k, v in sorted(descriptions.items(), key = lambda x: x[0].lower()):
        doc += '\t\t\t<li><a href="#%s">%s</a></li>\n' % \
            (k, v['full_name'] if 'full_name' in v else v['label'] if 'label' in v else k)
    doc += '\t</ul>\n'
    for k, v in sorted(descriptions.items(), key = lambda x: x[0].lower()):
        doc += '\t\t<h2 id="%s">%s</h2>\n' % \
            (k, v['full_name'] if 'full_name' in v else v['label'] if 'label' in v else k)
        if 'year' in v:
            doc += '\t\t\t<h3>Last updated: %u<\h3>\n'%v['year']
        if 'releases' in v:
            doc += '\t\t\t<p><b>Updated in years: </b>%s</p>\n' % \
                ', '.join(['%u'%y for y in v['releases']])
        if 'authors' in v and v['authors'] is not None:
            doc += '\t\t\t<p><b>Created by </b>%s</p>\n' % ', '.join(v['authors'])
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