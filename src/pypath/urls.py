#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import os

import pypath.common as common

urls = {
    'uniprot_pdb': {
        'label': 'Getting PDB IDs of 3D structures for UniProtIDs',
        'url': 'http://www.uniprot.org/docs/pdbtosp.txt'
    },
    'uniprot_basic': {
        'label': 'URL for UniProt queries',
        'url': 'http://www.uniprot.org/uniprot/',
        'lists': 'http://www.uniprot.org/uploadlists/'
    },
    'corum': {
        'label':
        'CORUM is a database of protein complexes, downloadable in csv format',
        'url_old':
        'http://mips.helmholtz-muenchen.de/genre/proj/corum/allComplexes.csv',
        'url':
        'http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt'
    },
    'pfam_pdb': {
        'label': 'PDB-Pfam mapping and names of Pfam domains',
        'url':
        'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt'
    },
    'pfam_up': {
        'label': 'Mapping Pfam regions to UniProt',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/'
        'Pfam-A.regions.tsv.gz'
    },
    '3dcomplexes_contact': {
        'label': 'This file contains the topology definition of the complexes',
        'url': 'http://shmoo.weizmann.ac.il/elevy/3dcomplexV4/dataV4/'
        'contactDefinition.txt'
    },
    '3dcomplexes_correspondancy': {
        'label': 'This is the dictionary of chain names',
        'url': 'http://shmoo.weizmann.ac.il/elevy/3dcomplexV4/dataV4/'
        'pdb_chain_corresV2.txt'
    },
    'pdb_chains': {
        'label': 'Corresponding UniProt IDs and residue numbers for each chain'
        'in PDB structures',
        'url':
        'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_'
        'uniprot.tsv.gz'
    },
    'complex_portal': {
        'label': 'Complexes curated by IntAct',
        'url':
        'ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/psi25/'
    },
    'pisa_interfaces': {
        'label': 'Base URL for download interface data from PDBe PISA',
        'url': 'http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?'
    },
    'catalytic_sites': {
        'label': 'Catalytic Site Atlas',
        'url': 'http://www.ebi.ac.uk/thornton-'
        'srv/databases/CSA/downloads/CSA_2_0_121113.txt'
    },
    '3did_ddi': {
        'label': 'Domain-domain interactions derived from 3D structures',
        'url': 'http://3did.irbbarcelona.org/download/current/3did_flat.gz'
    },
    '3did_dmi': {
        'label': 'Domain-motif interactions from 3DID',
        'url': 'http://3did.irbbarcelona.org/download/current/3did_dmi_flat.gz'
    },
    'swissprot_full': {
        'label': 'Full UniProt/Swissprot database',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/'
        'current_release/uniprot_sprot.dat.gz'
    },
    'instruct_human': {
        'label':
        'Protein interactome networks annotated to 3D structural resolution',
        'url': 'http://instruct.yulab.org/download/sapiens.sin'
    },
    'i3d_human': {
        'label': 'Interactome3D representative dataset for human proteins',
        'url':
        'http://interactome3d.irbbarcelona.org/user_data/human/download/'
        'representative/interactions.dat'
    },
    'instruct_offsets': {
        'label': 'Offsets between PDB chains and UniProt sequences',
        'url': 'http://instruct.yulab.org/download/indexing_uniprot2pdb.txt'
    },
    'compleat': {
        'label': 'Curated and inferred complexes from multiple databases',
        'url': 'http://www.flyrnai.org/compleat/ComplexDownload'
        '?requestType=complexDownload&org=Human'
    },
    'switches.elm': {
        'label': 'Curated data on molecular switches in cellular regulation',
        'url': 'http://switches.elm.eu.org/downloads/switches.ELM-v1.txt'
    },
    'ols': {
        'label': 'WSDL interface for the Ontology Lookup Service',
        'url': 'http://www.ebi.ac.uk/ontology-lookup/OntologyQuery.wsdl'
    },
    'comppi': {
        'label': 'Compartmentalized PPI database',
        'url': 'http://comppi.linkgroup.hu/downloads'
    },
    'elm_int': {
        'label': 'Interactions between ELM instances and globular domains',
        'url': 'http://elm.eu.org/interactions/as_tsv'
    },
    'elm_depr': {
        'label': 'Deprecated ELM class names',
        'url': 'http://elm.eu.org/infos/browse_renamed.tsv'
    },
    'pdbsws': {
        'label': 'PDB-UniProt residue level mapping',
        'url': 'http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl'
    },
    'pdb_align': {
        'label': 'PDB-UniProt residue level mapping',
        'url':
        'http://pdb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query='
    },
    'pepcyber': {
        'label': 'MySQL injection to PEPCyber website :)',
        'url':
        'http://www.pepcyber.org/PPEP/search_result.php?domain=Any&ppbd_symbol'
        '=Any&search_field=symbol&query_value=%27+OR+1&binding_sequence=&go_id='
        'Any&Submit=Search'
    },
    'pepcyber_details': {
        'label': 'Interaction details from pepcyber',
        'url': 'http://www.pepcyber.org/PPEP/idetail.php?iid=%u'
    },
    'pdzbase': {
        'label': 'Manually curated interactions of PDZ domain proteins',
        'url': 'http://abc.med.cornell.edu/pdzbase/allinteractions'
    },
    'pdz_details': {
        'label': 'Details of interactions in PDZbase',
        'url': 'http://abc.med.cornell.edu/pdzbase/interaction_detail/%u'
    },
    'psite_reg': {
        'label': 'PhosphoSite annotated regulatory sites',
        'url': 'http://www.phosphosite.org/downloads/Regulatory_sites.gz'
    },
    'psite_bp': {
        'label': 'PhosphoSite kinase substrates in BioPAX format',
        'url': 'http://www.phosphosite.org/downloads/Kinase_substrates.owl.gz'
    },
    'psite_ac': {
        'label': 'PhosphoSite acetylation sites',
        'url':
        'http://www.phosphosite.org/downloads/Acetylation_site_dataset.gz'
    },
    'psite_kin': {
        'label': 'PhosphoSite kinase-substrate interactions',
        'url':
        'http://www.phosphosite.org/downloads/Kinase_Substrate_Dataset.gz'
    },
    'psite_me': {
        'label': 'PhosphoSite methylation sites',
        'url':
        'http://www.phosphosite.org/downloads/Methylation_site_dataset.gz'
    },
    'psite_ga': {
        'label': 'PhosphoSite O-GalNAc sites',
        'url': 'http://www.phosphosite.org/downloads/O-GalNAc_site_dataset.gz'
    },
    'psite_gl': {
        'label': 'PhosphoSite O-GlcNAc sites',
        'url': 'http://www.phosphosite.org/downloads/O-GlcNAc_site_dataset.gz'
    },
    'psite_p': {
        'label': 'PhosphoSite phosphorylation sites',
        'url':
        'http://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz'
    },
    'psite_sm': {
        'label': 'Sumoylation sites',
        'url':
        'http://www.phosphosite.org/downloads/Sumoylation_site_dataset.gz'
    },
    'psite_ub': {
        'label': 'Ubiquitination sites',
        'url':
        'http://www.phosphosite.org/downloads/Ubiquitination_site_dataset.gz'
    },
    'proteomic_ielm': {
        'label': 'Proteomic iELM',
        'url': 'http://i.elm.eu.org/test_submit/'
    },
    'ielm_domains': {
        'label': 'List of domains form iELM',
        'url': 'http://i.elm.eu.org/domains/'
    },
    'domino': {
        'label': 'Domino PPI and domain-motif database in MI-TAB format',
        'url': 'ftp://mint.bio.uniroma2.it/pub/domino/release/mitab/'
        '2009-10-22/2009-10-22-domino-full-binary.mitab26'
    },
    'hprd_all': {
        'label': 'HPRD all data in flat files',
        'url': 'http://www.hprd.org/RELEASE9/HPRD_FLAT_FILES_041310.tar.gz',
        'ptm_file': 'FLAT_FILES_072010/POST_TRANSLATIONAL_MODIFICATIONS.txt',
        'int_file': 'FLAT_FILES_072010/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt'
    },
    'p_elm': {
        'label': 'phosphoELM',
        'url':
        'http://phospho.elm.eu.org/dumps/phosphoELM_vertebrate_latest.dump.tgz',
        'psites': 'phosphoELM_vertebrate_'
    },
    'p_elm_kin': {
        'label': 'List of kinases from phosphoELM',
        'url': 'http://phospho.elm.eu.org/kinases.html'
    },
    'elm_inst': {
        'label': 'List of ELM instances',
        'url': 'http://elm.eu.org/elms/browse_instances.tsv?q=*'
    },
    'elm_class': {
        'label': 'List of ELM classes',
        'url': 'http://elm.eu.org/elms/browse_elms.tsv'
    },
    'dbptm_benchmark': {
        'label': 'dbPTM is a PTM database compiled from multiple other dbs',
        'urls': [
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/N-linked_Glycosylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/N6-succinyllysine.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/O-linked_Glycosylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Phosphorylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Phosphorylation_CDK1.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Phosphorylation_CK2.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Phosphorylation_MAPK1.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Phosphorylation_PKA.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Phosphorylation_PKC.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/S-nitrosylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Amidation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Hydroxylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Acetylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Methylation.tgz'
        ]
    },
    'dbptm': {
        'label': 'dbPTM is a PTM database compiled from multiple other dbs',
        'urls': [
            'http://dbptm.mbc.nctu.edu.tw/download/N-linked.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/O-linked.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/C-linked.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Phosphorylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Acetylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Methylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Myristoylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Palmitoylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Prenylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Carboxylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Sulfation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Ubiquitylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Sumoylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/download/Nitrosylation.tgz'
        ]
    },
    'phosnw': {
        'label': 'Human kinase-substrate relationships and phosphosites',
        'url': 'http://phosphonetworks.org/download/highResolutionNetwork.csv'
    },
    'kegg_pws': {
        'label': 'KEGG Pathways',
        'list_url': 'http://www.genome.jp/kegg/pathway.html',
        'kgml_url':
        'http://www.kegg.jp/kegg-bin/download?entry=%s&format=kgml',
        'biopax_l3': 'http://www.pathwaycommons.org/archives/PC2/'
        'v8/PathwayCommons.8.kegg.BIOPAX.owl.gz'
    },
    'depod': {
        'label': 'Dephosphorylation substrates and sites',
        'urls': [
            'http://depod.bioss.uni-freiburg.de/download/DEPOD_201408'
            '_human_phosphatase-substrate.txt',
            'http://depod.bioss.uni-freiburg.de/download/DEPOD_201405'
            '_human_phosphatase-substrate.mitab'
        ]
    },
    'mimp': {
        'label': 'Kinase-substrate relationships',
        'url': 'http://mimp.baderlab.org/fetch_data/phosphorylation_data'
        '.tab/phosphorylation_data.tab'
    },
    'unip_iso': {
        'label': 'Isoform sequences from UniProt',
        'url': 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/'
        'knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz'
    },
    'kinclass': {
        'label': 'Kinase families and groups',
        'url':
        'http://kinase.com/static/colt/data/human/kinome/tables/Table%20S2.txt'
    },
    'protdb_exp': {
        'label': 'Expression data from ProteomicsDB',
        'url': 'https://www.proteomicsdb.org/proteomicsdb/logic/api/'
        'proteinexpression.xsodata/InputParams(PROTEINFILTER='
        '%%27%s%%27,MS_LEVEL=%u,TISSUE_ID_SELECTION=%%27%%27,'
        'TISSUE_CATEGORY_SELECTION=%%27tissue;fluid%%27,SCOPE_SELECTION=%u,'
        'GROUP_BY_TISSUE=1,CALCULATION_METHOD=0,EXP_ID=-1)/Results?'
        '$select=UNIQUE_IDENTIFIER,TISSUE_ID,TISSUE_NAME,TISSUE_SAP_SYNONYM,'
        'SAMPLE_ID,SAMPLE_NAME,AFFINITY_PURIFICATION,EXPERIMENT_ID,'
        'EXPERIMENT_NAME,EXPERIMENT_SCOPE,EXPERIMENT_SCOPE_NAME,PROJECT_ID,'
        'PROJECT_NAME,PROJECT_STATUS,UNNORMALIZED_INTENSITY,'
        'NORMALIZED_INTENSITY,MIN_NORMALIZED_INTENSITY,'
        'MAX_NORMALIZED_INTENSITY,SAMPLES&$format=json',
        'subs': ['uniprot', 'ms_level', 'scope']
    },
    'protdb_tis': {
        'label': 'Get all tissues where a given protein is expressed',
        'url': 'https://www.proteomicsdb.org/proteomicsdb/logic/api/'
        'proteinspertissue.xsodata/InputParams(TISSUE_ID='
        '%%27%s%%27,CALCULATION_METHOD=0,SWISSPROT_ONLY=%u,'
        'NO_ISOFORM=%u)/Results?$select=ENTRY_NAME,UNIQUE_IDENTIFIER,DATABASE,'
        'PROTEIN_DESCRIPTION,PEPTIDES,SAMPLE_NAME,SAMPLE_DESCRIPTION,'
        'UNNORMALIZED_EXPRESSION,NORMALIZED_EXPRESSION&$format=xml',
        'subs': ['bto', 'swissprot_only', 'isoform']
    },
    'abs': {
        'label': '',
        'url': 'http://genome.crg.es/datasets/abs2005/data/abs.gff'
    },
    'uniprot_sec': {
        'label': 'Secondary UniProt ACs',
        'url': 'ftp://ftp.uniprot.org/pub/databases/uniprot/'
        'knowledgebase/docs/sec_ac.txt'
    },
    'uniprot_idmap_ftp': {
        'label': 'Human ID mapping from UniProt',
        'url': 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/'
        'knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz'
    },
    'uniprot_idmap': {
        'label': 'ID mapping from UniProt',
        'url': 'http://uniprot.org/mapping/'
    },
    'pazar': {
        'label': 'TF-target gene lists from PAZAR',
        'url': 'http://www.pazar.info/tftargets/tftargets.zip'
    },
    'htri': {
        'label': 'TF-target gene lists from HTRI',
        'url':
        'http://www.lbbc.ibb.unesp.br/htri/consulta?type=1&all=true&down=3'
        '&iconss1.x=57&iconss1.y=48',
        'init_url': 'http://www.lbbc.ibb.unesp.br/htri/pagdown.jsp'
    },
    'oreganno_old': {
        'label': 'TF-target gene lists from ORegAnno, previous version',
        'url': 'http://www.oreganno.org/oregano/htdocs'
        '/data/oreganno_UCSC_08Nov10.txt.gz'
    },
    'oreganno': {
        'label': 'TF-target gene lists from ORegAnno',
        'url': 'http://py-gi1.stanford.edu:8080/oregano/htdocs/'
        'downloads/ORegAnno_Combined_2015.09.16.tsv'
    },
    'cpdb': {
        'label': 'All human interactions from ConsensusPathDB',
        'url':
        'http://cpdb.molgen.mpg.de/download/ConsensusPathDB_human_PPI.gz'
    },
    'goa': {
        'label': 'UniProt GO annotations from GOA',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/%s/'
        'goa_%s.gaf.gz'
    },
    'quickgo': {
        'label': 'UniProt GO annotations from QuickGO',
        'url': 'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&'
        'limit=-1%s&termUse=%s&tax=%u&col=proteinID,goID,goName,aspect'
    },
    'goslim_gen': {
        'label': 'Generic GOSlim from GO Consortium',
        'url':
        'http://www.geneontology.org/ontology/subsets/goslim_generic.obo'
    },
    'go': {
        'label': 'The whole Gene Ontology',
        'url': 'http://purl.obolibrary.org/obo/go.obo'
    },
    'netpath_names': {
        'label': 'NetPath numeric pathway IDs can be translated to '
        'pathway names only by extracting them from HTML....',
        'url': 'http://www.netpath.org/browse'
    },
    'netpath_psimi': {
        'label': 'Batch download of NetPath pathways in PSI-MI format',
        'url': 'http://www.netpath.org/download/zipped/PSI-MI.zip'
    },
    'netpath_bp': {
        'label': 'Netpath pathways one by one in BioPAX level 3 format',
        'biopax_l3': 'http://www.netpath.org/data/biopax/NetPath_%u.owl'
    },
    'proteomemap': {
        'label':
        'Human Proteome Map: Mass-spec expression data in healthy human tissues',
        'url':
        'http://www.humanproteomemap.org/Download_HPM/HPM_protein_level_'
        'expression_matrix_Kim_et_al_052914.csv'
    },
    'proteinatlas': {
        'label': 'Human Protein Atlas: Immuncytochemistry expression data in '
        'healthy human cells or cancer cells or cell lines',
        'normal': 'http://www.proteinatlas.org/download/normal_tissue.csv.zip',
        'cancer': 'http://www.proteinatlas.org/download/cancer.csv.zip'
    },
    'lincs-compounds': {
        'label':
        'List of small molecules in LINCS with synonyms and PubChem IDs',
        'url': 'http://lincs.hms.harvard.edu/db/datasets/20000/smallmolecules'
        '?search=&output_type=.csv'
    },
    'pubmed-eutils': {
        'label': 'Retrieving summary of PubMed records',
        'url': 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',
        'conv':
        'http://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=%s&format=json'
    },
    'pubmed': {
        'label': 'PubMed baseurl',
        'url': 'http://www.ncbi.nlm.nih.gov/pubmed/%s'
    },
    'hpmr': {
        'label': 'Human Plasma Membrane Receptome',
        'url':
        'http://receptome.stanford.edu/hpmr/SearchDB/findGenes.asp?textName=*'
    },
    'tfcensus': {
        'label': 'A census of human transcription factors (Vaquerizas 2009)',
        'url': 'http://www.nature.com/nrg/journal/v10/n4/extref/nrg2538-s3.txt'
    },
    'tfcheckp': {
        'label': 'A comprehensive list of transcription factors',
        'url':
        'http://www.tfcheckpoint.org/data/TFCheckpoint_download_180515.txt'
    },
    'gtp': {
        'label': 'Guide to Pharmacology: ligand-receptor interactions',
        'url': 'http://www.guidetopharmacology.org/DATA/interactions.csv'
    },
    'giant': {
        'label': 'JSON network from GIANT by gene based query',
        'init_url': 'http://giant.princeton.edu/',
        'url':
        'http://giant.princeton.edu/contexts/integration/query/%u/?%s&prior=0.1'
    },
    'msigdb': {
        'label': 'Molecular Signatures Database',
        'login1': 'http://software.broadinstitute.org/gsea/login.jsp',
        'login2':
        'http://software.broadinstitute.org/gsea/j_spring_security_check',
        'url': 'http://www.broadinstitute.org/gsea/msigdb/download_file.jsp?'
        'filePath=/resources/msigdb/5.0/%s.all.v5.0.%s.gmt',
        'coll': 'http://www.broadinstitute.org/gsea/msigdb/collections.jsp',
        'url_stem': 'http://www.broadinstitute.org/gsea/%s',
        'one_set':
        'http://www.broadinstitute.org/gsea/msigdb/download_geneset.jsp?'
        'geneSetName=%s&fileType=txt'
    },
    'disgenet': {
        'label': 'Disease-gene associations',
        'url':
        'http://www.disgenet.org/ds/DisGeNET/results/%s_gene_disease_associations.tar.gz',
        'datasets': ['curated', 'literature', 'befree', 'all']
    },
    'hsn': {
        'label': 'The Wang Human Signaling Network',
        'url': 'http://www.cancer-systemsbiology.org/HuamnSignalingNet_v6.csv'
    },
    'li2012': {
        'label': 'Human Phosphotyrosine Signaling Network',
        'url': 'http://genome.cshlp.org/content/suppl/2011/12/27/'
        'gr.128819.111.DC1/Supplementary_files_S1-S5.xls'
    },
    'trip': {
        'label': 'The TRP channel database',
        'base': 'http://www.trpchannel.org/',
        'show': 'http://www.trpchannel.org/proteins/show?id=%s',
        'intr': 'http://www.trpchannel.org/interactions/show?trp=%s&'
        'interactor=%s',
        'url': 'http://www.trpchannel.org/20141116.csv',
        'json': 'http://www.trpchannel.org/proteins/getjson'
    },
    'signor': {
        'label': 'SIGNOR pathways',
        'list_url': 'http://signor.uniroma2.it/downloads.php',
        'base_url': 'http://signor.uniroma2.it/%s',
        'all_url': 'http://signor.uniroma2.it/DownloadServlet?action=all',
        'all_url_new': 'http://signor.uniroma2.it/download_entity.php'
    },
    'hiii14': {
        'label': 'Rolland et al 2014, Human Interactome II',
        'url':
        'http://www.cell.com/cms/attachment/2021150198/2041343422/mmc3.xlsx'
    },
    'kinome': {
        'label': 'List of human kinases',
        'url':
        'http://kinase.com/static/colt/data/human/kinome/tables/Kincat_Hsap.08.02.xls'
    },
    'dgidb': {
        'label': 'Druggable genes compiled from multiple resources',
        'main_url': 'http://dgidb.genome.wustl.edu/search_categories',
        'url': 'http://dgidb.genome.wustl.edu/druggable_gene_categories/%s?'
        'sources=BaderLabGenes,CarisMolecularIntelligence,FoundationOneGenes,GO,'
        'GuideToPharmacologyGenes,HopkinsGroom,MskImpact,RussLampel,dGene'
    },
    'reactome': {
        'label': 'The Reactome reaction network database',
        'sbml':
        'http://www.reactome.org/download/current/homo_sapiens.2.sbml.gz',
        'biopax_l3': 'http://www.reactome.org/download/current/biopax.zip',
        'biopax_l2': 'http://www.reactome.org/download/current/biopax2.zip'
    },
    'acsn': {
        'label': 'Atlas of Cancer Signaling Networks',
        'biopax_l3': 'https://acsn.curie.fr/files/acsn_v1.1.owl',
        'sif': 'https://acsn.curie.fr/files/acsn_ppi.sif'
    },
    'nci-pid': {
        'label': 'National Cancer Institute -- Pathway Interaction Database',
        'biopax_l3':
        'ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_3/NCI-Nature_Curated.bp3.owl.gz'
    },
    'panther': {
        'label': 'PANTHER pathways',
        'biopax_l3':
        'ftp://ftp.pantherdb.org//pathway/current_release/BioPAX.tar.gz'
    },
    'laudanna': {
        'label':
        'Directionality and effect information inferred from multiple sources by Laudanna Lab',
        'sigflow':
        'http://dp.univr.it/~laudanna/LCTST/downloads/files/SignalingFlow.EA',
        'sigdir':
        'http://dp.univr.it/~laudanna/LCTST/downloads/files/SignalingDirection.EA'
    },
    'biogrid': {
        'label': 'BioGRID version 2 tab format',
        'url':
        'http://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-MV-Physical-LATEST.tab2.zip'
    },
    'wang': {
        'label':
        'Human Signaling Network compiled by Wang Group, version 6 (2014)',
        'url':
        'http://www.bri.nrc.ca/wang/cancerMap/HumanSignalingNetwork_v6.csv'
    },
    'graphviz': {
        'label': 'List of Graphviz attributes',
        'url': 'http://www.graphviz.org/doc/info/attrs.html'
    },
    'hid': {
        'label': 'CCSB Human Interactome Project',
        'lit-bm-13':
        'http://interactome.dfci.harvard.edu/H_sapiens/download/Lit-BM-13.tsv',
        'hi-ii-14':
        'http://interactome.dfci.harvard.edu/H_sapiens/download/HI-II-14.tsv',
        'hi-i-05':
        'http://interactome.dfci.harvard.edu/H_sapiens/download/HI-I-05.tsv'
    },
    'ca1': {
        'label': 'Supplementary Online Materials for Ma\'ayan 2005',
        'url': 'http://science.sciencemag.org/highwire/filestream/586741/'
        'field_highwire_adjunct_files/1/Maayan_SOM_External_Files.zip'
    },
    'ccmap': {
        'label': 'Cancer Cell Map from PathwayCommons 2011 snapshot',
        'edges':
        'http://www.pathwaycommons.org/archives/PC1/last_release-2011/tab_delim_network/by_source/cell-map-edge-attributes.txt.zip',
        'nodes':
        'http://www.pathwaycommons.org/archives/PC1/last_release-2011/tab_delim_network/by_source/cell-map-node-attributes.txt.zip'
    },
    'pathguide': {
        'label':
        'Collection of metabolic and signaling pathway and molecular interaction resources',
        'url':
        'http://pathguide.org/fullrecord.php?organisms=all&availability=all&standards=all&order=alphabetic&DBID=%u'
    },
    'cgc': {
        'label': 'Cancer Gene Census: list of cancer related (driver) genes',
        'host': 'sftp-cancer.sanger.ac.uk',
        'file': '/files/grch38/cosmic/v76/cancer_gene_census.csv'
    },
    'havugimana': {
        'label': 'Census of human soluble protein complexes',
        'url':
        'http://www.cell.com/cms/attachment/2021768736/2041631145/mmc3.xls'
    },
    'matrixdb': {
        'label': 'MatrixDB in house curated interactions, PSI-MI tab format',
        'old_url': 'http://matrixdb.ibcp.fr/download/matrixdb_CORE.tab.gz',
        'url': 'http://matrixdb.univ-lyon1.fr/download/matrixdb_CORE.tab.gz'
    },
    'innatedb': {
        'label': 'InnateDB PSI-MI tab',
        'url':
        'http://innatedb.com/download/interactions/innatedb_ppi.mitab.gz'
    },
    'dip': {
        'label': 'DIP PSI-MI tab',
        'login': 'http://dip.doe-mbi.ucla.edu/dip/Login.cgi',
        'url': 'http://dip.mbi.ucla.edu/dip/file?'
        'ds=current&fn=Hsapi20160430%s&ff=txt',
        'ik': 'http://dip.doe-mbi.ucla.edu/dip/DIPview.cgi?IK=%u'
    },
    'vaquerizas2009': {
        'label': 'A census of human transcription factors: function, '
        'expression and evolution; Supplementary Table S3',
        'url': 'http://www.nature.com/nrg/journal/v10/n4/extref/nrg2538-s3.txt'
    },
    'spike': {
        'label': 'SPIKE database XML',
        'url': 'http://www.cs.tau.ac.il/~spike/download/LatestSpikeDB.xml.zip'
    },
    'mppi': {
        'label': 'MIPS PPI full download',
        'url': 'http://mips.helmholtz-muenchen.de/proj/ppi/data/mppi.gz'
    },
    'negatome': {
        'label': 'Negatome manually curated non-interacting protein pairs',
        'manual':
        'http://mips.helmholtz-muenchen.de/proj/ppi/negatome/manual.txt'
    },
    'macrophage': {
        'label': 'Macrophage Pathways; Raza 2010, Supplementary Materials 1',
        'url': r'https://static-content.springer.com/esm/art%3A10.1186%2F1752-'
        r'0509-4-63/MediaObjects/12918_2010_452_MOESM2_ESM.XLS'
    },
    'intact': {
        'label': 'IntAct entire database PSI-MI tab',
        'mitab': 'ftp://ftp.ebi.ac.uk/pub/databases/intact/'
        'current/psimitab/intact.zip'
    },
    'death': {
        'label': 'DeathDomain webpage',
        'url': 'http://www.deathdomain.org/proteins/show?family=%s'
    },
    'string': {
        'label': 'STRING',
        'actions': 'http://string-db.org/download/'
        'protein.actions.v10.5/%u.protein.actions.v10.5.txt.gz',
        'links': 'http://string-db.org/download/protein.links.detailed.v10/%u'
        '.protein.links.detailed.v10.txt.gz'
    },
    'wikipw': {
        'label': 'WikiPathways human biopax',
        'biopax_l3': 'http://wikipathways.org//wpi/cache/wikipathways_'
        'Homo_sapiens_Curation-AnalysisCollection__owl.zip'
    },
    'pwcommons': {
        'label': 'PathwayCommons binary SIF files',
        'url': 'http://www.pathwaycommons.org/archives/PC2/v8/'
        'PathwayCommons.8.%s.BINARY_SIF.hgnc.txt.sif.gz'
    },
    'homologene': {
        'label': 'NCBI HomoloGene data, recent release',
        'url': 'ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data'
    },
    'mirbase' : {
        'label': 'miRBase: miRNA main reference database',
        'aliases': 'ftp://mirbase.org/pub/mirbase/CURRENT/aliases.txt.gz'
    },
    'mir2dis': {
        'label': 'miR2Disease experimentally validated'\
            'miRNA-target interactions',
        'url': 'http://watson.compbio.iupui.edu:8080'\
            '/miR2Disease/download/miRtar.txt'
    },
    'mirdeathdb': {
        'label': 'miRDeathDB experimentally verified'\
            'miRNA-target interactions',
        'url': 'http://www.rna-world.org/mirdeathdb/data'\
            '/miRDeathDB_all_data.txt'
    },
    'mirecords': {
        'label': 'miRecords experimentally validated'\
            'miRNA-target interactions',
        'url': 'http://c1.accurascience.com/miRecords/download_data.php?v=4'
    },
    'mirtarbase': {
        'label': 'miRTarBase experimentally validated'\
            'miRNA-target interactions',
        'strong': 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/'\
            '6.1/miRTarBase_SE_WR.xls',
        'all': 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/'\
            '6.1/miRTarBase_MTI.xlsx'
    },
    'lncdisease': {
        'label': 'lncRNA and disease database: experimentally '\
            'verified lncRNA interactions',
        'url': 'http://210.73.221.6/files/images/ldd/data2.txt'
    },
    'lncrnadb': {
        'label': 'Literature curated lncRNA interactions',
        'url': 'http://lncrna.com/rest/all/nomenclature/species/association'
    },
    'transmir': {
        'label': 'Literature curated TF-miRNA interactions',
        'url': 'http://www.cuilab.cn/files/images/transmir/'\
            'transmir_v1.2.txt'
    },
    'encode': {
        'label': 'Interaction data from the Nature ENCODE project',
        'tf-mirna': 'http://encodenets.gersteinlab.org/enets10.TF-miRNA.txt'
    },
    'imweb': {
        'label': 'Imweb interaction score/prediction',
        'url': 'https://www.intomics.com/inbio/map/api/'\
            'get_data?file=InBio_Map_core_2016_09_12.tar.gz',
        'login': 'https://www.intomics.com/inbio/api/login_guest?ref=&_=%u',
        'refresh': 'https://www.intomics.com/inbio/api/refresh?_=%u'
    }
}

files = {
    'signalink': {
        'edges': 'signalink3_edges.tsv',
        'nodes': 'signalink3_nodes.tsv'
    },
    'acsn': {
        'names': os.path.join(common.ROOT, 'data', 'acsn_names.gmt'),
        'ppi': os.path.join(common.ROOT, 'data', 'acsn_ppi.txt')
    },
    'phosphopoint': {
        'data': os.path.join(common.ROOT, 'data', 'phosphopoint.csv')
    },
    'phosphosite': {
        'curated': os.path.join('cache', 'phosphosite_curated.pickle'),
        'noref': os.path.join('cache', 'phosphosite_noref.pickle')
    }
}
