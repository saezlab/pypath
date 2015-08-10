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

# external modules:
import os

# from bioigraph:
import input_formats
import common

__all__ = ['urls', 'mapList', 'otherMappings', 'refLists', 
           'best', 'good', 'ugly', 'transcription', 'negative', 'gdsc_comp_target', 'cgc', 
           'mapListUniprot', 'mapListBasic']

ROOT = common.ROOT

urls = {
    'uniprot_pdb': {
        'label': 'Getting PDB IDs of 3D structures for UniProtIDs',
        'url': 'http://www.uniprot.org/docs/pdbtosp.txt'
    },
    'uniprot_basic': {
        'label': 'URL for UniProt queries',
        'url': 'http://www.uniprot.org/uniprot/'
    },
    'corum': {
        'label': 'CORUM is a database of protein complexes, downloadable in csv format',
        'url': 'http://mips.helmholtz-muenchen.de/genre/proj/corum/allComplexes.csv'
    },
    'pfam_pdb': {
        'label': 'PDB-Pfam mapping and names of Pfam domains',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt'
    },
    'pfam_up': {
        'label': 'Mapping Pfam regions to UniProt',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/'\
            'Pfam-A.regions.tsv.gz'
    },
    '3dcomplexes_contact': {
        'label': 'This file contains the topology definition of the complexes',
        'url': 'http://shmoo.weizmann.ac.il/elevy/3dcomplexV4/dataV4/'\
            'contactDefinition.txt'
    },
    '3dcomplexes_correspondancy': {
        'label': 'This is the dictionary of chain names',
        'url': 'http://shmoo.weizmann.ac.il/elevy/3dcomplexV4/dataV4/'\
            'pdb_chain_corresV2.txt'
    },
    'pdb_chains': {
        'label': 'Corresponding UniProt IDs and residue numbers for each chain'\
            'in PDB structures',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_'\
            'uniprot.tsv.gz'
    },
    'complex_portal': {
        'label': 'Complexes curated by IntAct',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/psi25/'
    },
    'pisa_interfaces': {
        'label': 'Base URL for download interface data from PDBe PISA',
        'url': 'http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?'
    },
    'catalytic_sites': {
        'label': 'Catalytic Site Atlas',
        'url': 'http://www.ebi.ac.uk/thornton-'\
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
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/'\
            'current_release/uniprot_sprot.dat.gz'
    },
    'instruct_human': {
        'label': 'Protein interactome networks annotated to 3D structural resolution',
        'url': 'http://instruct.yulab.org/download/sapiens.sin'
    },
    'i3d_human': {
        'label': 'Interactome3D representative dataset for human proteins',
        'url': 'http://interactome3d.irbbarcelona.org/user_data/human/download/'\
            'representative/interactions.dat'
    },
    'instruct_offsets': {
        'label': 'Offsets between PDB chains and UniProt sequences',
        'url': 'http://instruct.yulab.org/download/indexing_uniprot2pdb.txt'
    },
    'compleat': {
        'label': 'Curated and inferred complexes from multiple databases',
        'url': 'http://www.flyrnai.org/compleat/ComplexDownload'\
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
        'url': 'http://pdb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query='
    },
    'pepcyber': {
        'label': 'MySQL injection to PEPCyber website :)',
        'url': 'http://www.pepcyber.org/PPEP/search_result.php?domain=Any&ppbd_symbol'\
            '=Any&search_field=symbol&query_value=%27+OR+1&binding_sequence=&go_id='\
            'Any&Submit=Search'
    },
    'pepcyber_details': {
        'label': 'interaction details from pepcyber',
        'url': 'http://www.pepcyber.org/PPEP/idetail.php?iid=%u'
    },
    'pdzbase': {
        'label': 'manually curated interactions of PDZ domain proteins',
        'url': 'http://abc.med.cornell.edu/pdzbase/allinteractions'
    },
    'pdz_details': {
        'label': 'details of interactions in PDZbase',
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
        'url': 'http://www.phosphosite.org/downloads/Acetylation_site_dataset.gz'
    },
    'psite_kin': {
        'label': 'PhosphoSite kinase-substrate interactions',
        'url': 'http://www.phosphosite.org/downloads/Kinase_Substrate_Dataset.gz'
    },
    'psite_met': {
        'label': 'PhosphoSite methylation sites',
        'url': 'http://www.phosphosite.org/downloads/Methylation_site_dataset.gz'
    },
    'psite_gal': {
        'label': 'PhosphoSite O-GalNAc sites',
        'url': 'http://www.phosphosite.org/downloads/O-GalNAc_site_dataset.gz'
    },
    'psite_glc': {
        'label': 'PhosphoSite O-GlcNAc sites',
        'url': 'http://www.phosphosite.org/downloads/O-GlcNAc_site_dataset.gz'
    },
    'psite_p': {
        'label': 'PhosphoSite phosphorylation sites',
        'url': 'http://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz'
    },
    'psite_sum': {
        'label': 'Sumoylation sites',
        'url': 'http://www.phosphosite.org/downloads/Sumoylation_site_dataset.gz'
    },
    'psite_ub': {
        'label': 'Ubiquitination sites',
        'url': 'http://www.phosphosite.org/downloads/Ubiquitination_site_dataset.gz'
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
        'url': 'ftp://mint.bio.uniroma2.it/pub/domino/release/mitab/'\
            '2009-10-22/2009-10-22-domino-full-binary.mitab26'
    },
    'hprd_all': {
        'label': 'HPRD all data in flat files',
        'url': 'http://www.hprd.org/RELEASE9/HPRD_FLAT_FILES_041310.tar.gz',
        'ptm_file': 'FLAT_FILES_072010/POST_TRANSLATIONAL_MODIFICATIONS.txt'
    },
    'p_elm': {
        'label': 'phosphoELM',
        'url': 'http://phospho.elm.eu.org/dumps/phosphoELM_vertebrate_latest.dump.tgz',
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
    'dbptm': {
        'label': 'dbPTM is a PTM database compiled from multiple other dbs',
        'urls': [
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/N-linked.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/O-linked.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/C-linked.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Phosphorylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Acetylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Methylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Myristoylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Palmitoylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Prenylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Carboxylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Sulfation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Ubiquitylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Sumoylation.tgz',
            'http://dbptm.mbc.nctu.edu.tw/Benchmark/Nitrosylation.tgz'
        ]
    },
    'phosnw': {
        'label': 'Human kinase-substrate relationships and phosphosites',
        'url': 'http://phosphonetworks.org/download/highResolutionNetwork.csv'
    },
    'depod': {
        'label': 'Dephosphorylation substrates and sites',
        'urls': [
            'http://www.koehnlab.de/depod/download/DEPOD_201408'\
                '_human_phosphatase-substrate.txt',
            'http://www.koehnlab.de/depod/download/DEPOD_201405'\
                '_human_phosphatase-substrate.mitab'
        ]
    },
    'mimp': {
        'label': 'Kinase-substrate relationships',
        'url': 'http://mimp.baderlab.org/fetch_data/phosphorylation_data'\
            '.tab/phosphorylation_data.tab'
    },
    'unip_iso': {
        'label': 'Isoform sequences from UniProt',
        'url': 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/'\
            'knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz'
    },
    'kinclass': {
        'label': 'Kinase families and groups',
        'url': 'http://kinase.com/static/colt/data/human/kinome/tables/Table%20S2.txt'
    },
    'acsn': {
        'label': 'Atlas of Cancer Signaling Networks',
        'url': 'https://acsn.curie.fr/files/acsn_ppi.sif'
    },
    'protdb_exp': {
        'label': 'Expression data from ProteomicsDB',
        'url': 'https://www.proteomicsdb.org/proteomicsdb/logic/api/'\
            'proteinexpression.xsodata/InputParams(PROTEINFILTER='\
            '%%27%s%%27,MS_LEVEL=%u,TISSUE_ID_SELECTION=%%27%%27,'\
            'TISSUE_CATEGORY_SELECTION=%%27tissue;fluid%%27,SCOPE_SELECTION=%u,'\
            'GROUP_BY_TISSUE=1,CALCULATION_METHOD=0,EXP_ID=-1)/Results?'\
            '$select=UNIQUE_IDENTIFIER,TISSUE_ID,TISSUE_NAME,TISSUE_SAP_SYNONYM,'\
            'SAMPLE_ID,SAMPLE_NAME,AFFINITY_PURIFICATION,EXPERIMENT_ID,'\
            'EXPERIMENT_NAME,EXPERIMENT_SCOPE,EXPERIMENT_SCOPE_NAME,PROJECT_ID,'\
            'PROJECT_NAME,PROJECT_STATUS,UNNORMALIZED_INTENSITY,'\
            'NORMALIZED_INTENSITY,MIN_NORMALIZED_INTENSITY,'\
            'MAX_NORMALIZED_INTENSITY,SAMPLES&$format=json',
        'subs': ['uniprot', 'ms_level', 'scope']
    },
    'protdb_tis': {
        'label': 'Get all tissues where a given protein is expressed',
        'url': 'https://www.proteomicsdb.org/proteomicsdb/logic/api/'\
            'proteinspertissue.xsodata/InputParams(TISSUE_ID='\
            '%%27%s%%27,CALCULATION_METHOD=0,SWISSPROT_ONLY=%u,'\
            'NO_ISOFORM=%u)/Results?$select=ENTRY_NAME,UNIQUE_IDENTIFIER,DATABASE,'\
            'PROTEIN_DESCRIPTION,PEPTIDES,SAMPLE_NAME,SAMPLE_DESCRIPTION,'\
            'UNNORMALIZED_EXPRESSION,NORMALIZED_EXPRESSION&$format=xml',
        'subs': ['bto', 'swissprot_only', 'isoform']
    },
    'abs': {
        'label': '',
        'url': 'http://genome.crg.es/datasets/abs2005/data/abs.gff'
    },
    'uniprot_sec': {
        'label': 'Secondary UniProt ACs',
        'url': 'ftp://ftp.uniprot.org/pub/databases/uniprot/'\
            'knowledgebase/docs/sec_ac.txt'
    },
    'uniprot_idmap_ftp': {
        'label': 'Human ID mapping from UniProt',
        'url': 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/'\
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
        'url': 'http://www.lbbc.ibb.unesp.br/htri/consulta?type=1&all=true&down=3'\
            '&iconss1.x=57&iconss1.y=48',
        'init_url': 'http://www.lbbc.ibb.unesp.br/htri/pagdown.jsp'
    },
    'oreganno': {
        'label': 'TF-target gene lists from ORegAnno',
        'url': 'http://www.oreganno.org/oregano/htdocs'\
            '/data/oreganno_UCSC_08Nov10.txt.gz'
    },
    'cpdb': {
        'label': 'All human interactions from ConsensusPathDB',
        'url': 'http://cpdb.molgen.mpg.de/download/ConsensusPathDB_human_PPI.gz'
    },
    'goa': {
        'label': 'UniProt GO annotations from GOA',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/%s/'\
            'gene_association.goa_%s.gz'
    },
    'quickgo': {
        'label': 'UniProt GO annotations from QuickGO',
        'url': 'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&'\
            'limit=-1%s&termUse=%s&tax=%u&col=proteinID,goID,goName,aspect'
    },
    'goslim_gen': {
        'label': 'Generic GOSlim from GO Consortium',
        'url': 'http://www.geneontology.org/ontology/subsets/goslim_generic.obo'
    },
    'netpath_names': {
        'label': 'NetPath numeric pathway IDs can be translated to '\
            'pathway names only by extracting them from HTML....',
        'url': 'http://www.netpath.org/browse'
    },
    'netpath_psimi': {
        'label': 'Batch download of NetPath pathways in PSI-MI format',
        'url': 'http://www.netpath.org/download/zipped/PSI-MI.zip'
    },
    'proteomemap': {
        'label': 'Human Proteome Map: Mass-spec expression data in healthy human tissues',
        'url': 'http://www.humanproteomemap.org/Download_HPM/HPM_protein_level_'\
            'expression_matrix_Kim_et_al_052914.csv'
    },
    'proteinatlas_normal': {
        'label': 'Human Protein Atlas: Immuncytochemistry expression data in '\
            'healthy human cells',
        'url': 'http://www.proteinatlas.org/download/normal_tissue.csv.zip'
    },
    'proteinatlas_cancer': {
        'label': 'Human Protein Atlas: Immuncytochemistry expression data in '\
            'human tumour cells',
        'url': 'http://www.proteinatlas.org/download/cancer.csv.zip'
    },
    'lincs-compounds': {
        'label': 'List of small molecules in LINCS with synonyms and PubChem IDs',
        'url': 'http://lincs.hms.harvard.edu/db/datasets/20000/smallmolecules'\
            '?search=&output_type=.csv'
    },
    'pubmed-eutils': {
        'label': 'Retrieving summary of PubMed records',
        'url': 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    },
    'pubmed': {
        'label': 'PubMed baseurl',
        'url': 'http://www.ncbi.nlm.nih.gov/pubmed/%s'
    },
    'hpmr': {
        'label': 'Human Plasma Membrane Receptome',
        'url': 'http://receptome.stanford.edu/hpmr/SearchDB/findGenes.asp?textName=*'
    },
    'tfcensus': {
        'label': 'A census of human transcription factors (Vaquerizas 2009)',
        'url': 'http://www.nature.com/nrg/journal/v10/n4/extref/nrg2538-s3.txt'
    },
    'gtp': {
        'label': 'Guide to Pharmacology: ligand-receptor interactions',
        'url': 'http://www.guidetopharmacology.org/DATA/interactions.csv'
    },
    'giant': {
        'label': 'JSON network from GIANT by gene based query',
        'init_url': 'http://giant.princeton.edu/',
        'url': 'http://giant.princeton.edu/contexts/integration/query/%u/?%s&prior=0.1'
    },
    'msigdb': {
        'label': 'Molecular Signatures Database', 
        'login1': 'http://www.broadinstitute.org/gsea/login.jsp',
        'login2': 'http://www.broadinstitute.org/gsea/j_spring_security_check',
        'url': 'http://www.broadinstitute.org/gsea/msigdb/download_file.jsp?'\
            'filePath=/resources/msigdb/5.0/%s.all.v5.0.%s.gmt',
        'coll': 'http://www.broadinstitute.org/gsea/msigdb/collections.jsp',
        'url_stem': 'http://www.broadinstitute.org/gsea/%s',
        'one_set': 'http://www.broadinstitute.org/gsea/msigdb/download_geneset.jsp?'\
            'geneSetName=%s&fileType=txt'
    },
    'disgenet': {
        'label': 'Disease-gene associations',
        'url': 'http://www.disgenet.org/ds/DisGeNET/results/%s_gene_disease_associations.tar.gz',
        'datasets': ['curated', 'literature', 'befree', 'all']
    },
    'hsn': {
        'label': 'The Wang Human Signaling Network',
        'url': 'http://www.cancer-systemsbiology.org/HuamnSignalingNet_v6.csv'
    },
    'li2012': {
        'label': 'Human Phosphotyrosine Signaling Network',
        'url': 'http://genome.cshlp.org/content/22/7/1222/suppl/DC1',
        'file': 'li2012.csv'
    },
    'trip': {
        'label': 'The TRP channel database',
        'base': 'http://www.trpchannel.org/',
        'show': 'http://www.trpchannel.org/proteins/show?id=%s',
        'intr': 'http://www.trpchannel.org/interactions/show?trp=%s&'\
            'interactor=%s',
        'url': 'http://www.trpchannel.org/20141116.csv',
        'json': 'http://www.trpchannel.org/proteins/getjson'
    }
}

# this is all what is needed to load the resources 
# included in the bioigraph package
mapList = [
    {
        "one": "uniprot-sec",
        "two": "uniprot-pri",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data',
                                                'sec_ac.txt'),0,1,None,header=0)
    },
    {
        "one": "trembl",
        "two": "genesymbol",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data',
                    'trembl3.tab'),0,1,"\t",header=0)
    },
    {
        "one": "genesymbol",
        "two": "swissprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data',
                    'swissprot3.tab'),1,0,"\t",header=0)
    },
    {
        "one": "genesymbol",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data',
                'uniprot3.tab'),1,0,"\t",header=0,bi=True)
    },
    {
        "one": "genesymbol-fallback",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data', 
                    'human-genesymbol-all.tab'),1,0,"\t",header=0)
    },
    {
        "one": "refseq",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data', 
                    'uniprot-refseq-human-1.tab'),1,0,"\t",header=0)
    },
    {
        "one": "entrez",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data', 'entrez_uniprot.csv'),1,0,";",header=0)
    },
    {
        "one": "hgnc",
        "two": "uniprot",
        "typ": "protein",
        "src": "mysql",
        "par": input_formats.MysqlMapping("hgnc_new","gsy","u","mapping",None,bi=True)
    },
    {
        "one": "uniprot",
        "two": "hgncapprov",
        "typ": "protein",
        "src": "mysql",
        "par": input_formats.MysqlMapping("hgnc_prim","u","gsy","mapping",None,bi=True)
    },
    {
        "one": "uniprot",
        "two": "hgnc",
        "typ": "protein",
        "src": "mysql",
        "par": input_formats.MysqlMapping("hgnc_names","u","gsy","mapping",None)
    }
    #,
    #{
        #"one": "uniprot",
        #"two": "hgnc",
        #"typ": "protein",
        #"src": "mysql",
        #"par": input_formats.UniprotMapping("hgnc_names","u","gsy","mapping",None)
    #}
]

mapListUniprot = {
    ('embl', 'uniprot'): input_formats.UniprotMapping('embl'), 
    ('genesymbol', 'uniprot'): input_formats.UniprotMapping('genesymbol', 
        bi = True), 
    ('genesymbol-syn', 'uniprot'): input_formats.UniprotMapping('genesymbol-syn'), 
    ('entrez', 'uniprot'): input_formats.UniprotMapping('entrez'), 
    ('hgnc', 'uniprot'): input_formats.UniprotMapping('hgnc'), 
    ('enst', 'uniprot'): input_formats.UniprotMapping('enst'), 
    ('refseqp', 'uniprot'): input_formats.UniprotMapping('refseqp'),
    ('uniprot-entry', 'uniprot'): input_formats.UniprotMapping('uniprot-entry'),
    ('protein-name', 'uniprot'): input_formats.UniprotMapping('protein-name'),
    ('protein-name-all', 'uniprot'): input_formats.UniprotMapping('protein-name', 
        swissprot = None)
}

mapListBasic = {
    ('uniprot-sec', 'uniprot-pri'): input_formats.FileMapping(
        'get_uniprot_sec', 0, 1, None, header = 0),
    ('genesymbol', 'trembl'): input_formats.UniprotMapping('genesymbol', 
        swissprot = 'no', bi = True),
    ('genesymbol', 'swissprot'): input_formats.UniprotMapping('genesymbol'),
    ('genesymbol-syn', 'swissprot'): input_formats.UniprotMapping('genesymbol-syn'),
    ('genesymbol', 'uniprot'): 
        input_formats.UniprotMapping('genesymbol', bi = True, swissprot = None)
}

# this is all what is needed for corrections of unirpot ids 
# i.e. to get primary swissprot id for all proteins
mapListUniprotOld = [
    {
        "one": "uniprot-sec",
        "two": "uniprot-pri",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data',
                                                'sec_ac.txt'),0,1,None,header=0)
    },
    {
        "one": "trembl",
        "two": "genesymbol",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data',
                    'trembl3.tab'),0,1,"\t",header=0)
    },
    {
        "one": "genesymbol",
        "two": "swissprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data',
                    'swissprot3.tab'),1,0,"\t",header=0)
    },
    {
        "one": "genesymbol-fallback",
        "two": "uniprot",
        "typ": "protein",
        "src": "file",
        "par": input_formats.FileMapping(os.path.join(ROOT, 'data', 
                    'human-genesymbol-all.tab'),1,0,"\t",header=0)
    }
]

otherMappings = [
    {
        "one": "entrez",
        "two": "uniprot",
        "typ": "protein",
        "src": "mysql",
        "par": input_formats.MysqlMapping("geneid","geneid","u","mapping","ncbi")
    },
    {
        "one": "uniprot",
        "two": "genesymbol",
        "typ": "protein",
        "src": "mysql",
        "par": input_formats.MysqlMapping("uniprot_gs","u","gs","mapping","ncbi")
    }
    ]

refLists = [
        input_formats.ReferenceList('uniprot','protein',9606,
                              os.path.join(ROOT, 'data','uniprot-all-human.tab'))
    ]

best = {
    'slk': input_formats.ReadSettings(name="SignaLink2", separator=",", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein",isDirected=(7,['1','2']), 
                sign=(6,'1','-1'),
                inFile=os.path.join(ROOT, 'data', 'slk01human.csv'),references=(9,':'),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "netbiol_effect": 8, 
                    "is_direct": 6,
                    "is_directed": 7},
                extraNodeAttrsA={
                    "slk_pathways": (4, ":"),
                    "gene_name": 2},
                extraNodeAttrsB={
                    "slk_pathways": (5, ":"),
                    "gene_name": 3}),
    'ccmap': input_formats.ReadSettings(name="CancerCellMap", 
                separator="\t", nameColA=3, nameColB=4,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein",
                isDirected = False, sign = False, ncbiTaxId=9606,
                inFile=os.path.join(ROOT, 'data', 'cell-map-edge-attributes.txt'),references=(6, ";"),
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={},
                header=True),
    'spike': input_formats.ReadSettings(name="SPIKE", separator="\t", nameColA=1, nameColB=3,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected=(4,['1']), 
                sign = (7, '1', '2'),
                inFile=os.path.join(ROOT, 'data', 'spike_hc.csv'),references=(5, ";"),ncbiTaxId=9606,
                extraEdgeAttrs={
                    'spike_effect': 7, 
                    'spike_mechanism': 11},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'mppi': input_formats.ReadSettings(name="MPPI", separator="|", nameColA=2, nameColB=6,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=False,sign=False,
                inFile=os.path.join(ROOT, 'data', 'mppi_human_rep.csv'),references=(0, ";"),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "mppi_evidences": (1, ";")},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'psite': input_formats.ReadSettings(name="PhosphoSite", separator="\t", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=True, sign=False,
                inFile=os.path.join(ROOT, 'data', 'phosphosite_human_hc.csv'),references=(5, ";"),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "psite_evidences": (4, ";")},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'dip': input_formats.ReadSettings(name="DIP", separator="\t", 
                nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein",isDirected=False,sign=False,
                inFile=os.path.join(ROOT, 'data', 'dip_human_core_processed.csv'),
                references=(2, ";"),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "dip_methods": (4, ";"),
                    "dip_type": (3, ";")},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'netpath': input_formats.ReadSettings(name = "NetPath", separator = None,
                nameColA = 1, nameColB = 3,
                nameTypeA="entrez", nameTypeB="entrez",
                typeA="protein", typeB="protein", isDirected=False, sign=False,
                inFile='netpath', references=(4, ";"), ncbiTaxId=9606,
                extraEdgeAttrs={
                    "netpath_methods": (5, ";"),
                    "netpath_type": (6, ";"),
                    "netpath_pathways": (7, ';')},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'guide2pharma': input_formats.ReadSettings(name = "Guide2Pharmacology", 
                separator = None, nameColA = 0, nameColB = 1,
                nameTypeA = "uniprot", nameTypeB = "uniprot",
                typeA = "protein", typeB = "protein", isDirected = True, 
                sign = (2, 1, -1), inFile = 'guide2pharma', 
                references=(3, "|"), ncbiTaxId = 9606,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {'g2p_ligand': 4},
                extraNodeAttrsB = {'g2p_receptor': 4}),
    'ca1': input_formats.ReadSettings(name="CA1", separator=";", 
                nameColA=1, nameColB=6,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", 
                isDirected=(10,['_','+']), sign=(10,'+','_'),
                header=True, 
                inFile=os.path.join(ROOT, 'data', 'ca1.csv'),
                references=(12, ";"), ncbiTaxId=9606,
                extraEdgeAttrs={
                    "ca1_effect": 10,
                    "ca1_type": 11},
                extraNodeAttrsA={
                    "ca1_location": 4,
                    "ca1_function": 3},
                extraNodeAttrsB={
                    "ca1_location": 9,
                    "ca1_function": 8}),
    'nci_pid': input_formats.ReadSettings(name="NCI-PID", 
                separator="\t", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=(2,['1','-1']), 
                inFile=os.path.join(ROOT, 'data', 'nci-pid-strict.csv'),
                references=(4, ";"),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "pid_effect": 2,
                    "pid_evidence": (5, ";"),
                    "pid_pathways": (6, ";")},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'arn': input_formats.ReadSettings(name="ARN", separator=",", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=(3,['1','2']), sign=(4,'1','-1'),
                inFile=os.path.join(ROOT, 'data', 'arn_curated.csv'),
                references=(7, ":"),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "netbiol_effect": 4,
                    "is_direct": 2,
                    "is_directed": 3
                    },
                extraNodeAttrsA={
                    "atg": 5
                    },
                extraNodeAttrsB={
                    "atg": 6
                    }),
    'nrf2': input_formats.ReadSettings(name="NRF2ome", 
                separator=",", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", 
                isDirected=(3,['1','2']), sign=(4,'1','-1'),
                inFile=os.path.join(ROOT, 'data', 'nrf2ome.csv'),
                references=(5, ":"),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "netbiol_effect": 4,
                    "is_direct": 2,
                    "is_directed": 3
                    },
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'macrophage': input_formats.ReadSettings(name="Macrophage", 
                separator=";", nameColA=0, nameColB=1,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected=(3,['1']),
                sign=(2,'Activation','Inhibition'),
                inFile=os.path.join(ROOT, 'data', 'macrophage-strict.csv'),
                references=(5, ","),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "macrophage_type": (2, ","),
                    "macrophage_location": (4, ",")},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'innatedb': input_formats.ReadSettings(name="InnateDB", 
                separator=";", nameColA=0, nameColB=2,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=False, sign=False,
                inFile=os.path.join(ROOT, 'data', 'innatedb.csv'),
                references=(4, ":"),ncbiTaxId=9606,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'alz': input_formats.ReadSettings(name="AlzPathway", 
                separator="\t", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=False, sign=False,
                inFile=os.path.join(ROOT, 'data', 'alzpw-ppi.csv'),references=(8, ";"),ncbiTaxId=9606,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'depod': input_formats.ReadSettings(name="DEPOD", separator=";", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=True, sign=False,
                inFile=os.path.join(ROOT, 'data', 'depod-refs.csv'),references=(2, "|"),ncbiTaxId=9606,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'matrixdb': input_formats.ReadSettings(name="MatrixDB", separator=";", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=False, sign=False,
                inFile=os.path.join(ROOT, 'data', 'matrixdb_core.csv'),references=(2, "|"),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "matrixdb_methods": (3, '|')
                    },
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'intact': input_formats.ReadSettings(name="IntAct", separator=",", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=False, sign=False,
                inFile=os.path.join(ROOT, 'data', 'intact_filtered.csv'),references=(2, ";"),ncbiTaxId=9606,
                extraEdgeAttrs={
                    "intact_methods": (3, ';')
                    },
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'death': input_formats.ReadSettings(name="DeathDomain", 
                separator="\t", nameColA=0, nameColB=1,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected=False, sign=False,
                inFile=os.path.join(ROOT, 'data', 'dd_refs.csv'),
                references=(3, ";"), ncbiTaxId = 9606,
                extraEdgeAttrs={
                    "dd_methods": (2, ';')
                    },
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'signor': input_formats.ReadSettings(name="Signor", separator="\t", 
                nameColA = 2, nameColB = 6,
                nameTypeA="uniprot", nameTypeB="uniprot",
                positiveFilters = [(22, 'YES')], # only direct interactions
                negativeFilters = [(9, 'transcriptional regulation')], # exclude TF-target interactions
                typeA="protein", typeB="protein", 
                isDirected = (8, ['up-regulates', 'up-regulates activity', 
                    'up-regulates quantity by stabilization',
                    'down-regulates', 'down-regulates activity', 
                    'down-regulates quantity by destabilization']), 
                sign = (8, ['up-regulates', 'up-regulates activity', 
                    'up-regulates quantity by stabilization'],
                    ['down-regulates', 'down-regulates activity', 
                    'down-regulates quantity by destabilization']), 
                ncbiTaxId=9606,
                inFile=os.path.join(ROOT, 'data', 'signor_22052015.tab'),references=(21, ";"),header=True,
                extraEdgeAttrs={
                    "signor_mechanism": (9, ';')
                    },
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'pdz': input_formats.ReadSettings(name="PDZBase", separator = None, nameColA=1,
                nameColB=4, nameTypeA="uniprot", nameTypeB="uniprot",
                typeA = "protein", typeB = "protein", isDirected = 1, sign = False,
                ncbiTaxId = {'col': 5, 'dict': {'human': 9606}}, 
                inFile = 'get_pdzbase', references = 6, header = False,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'lmpid': input_formats.ReadSettings(name="LMPID", separator = None, nameColA = 0,
                nameColB = 1, nameTypeA="uniprot", nameTypeB="uniprot",
                typeA = "protein", typeB = "protein", isDirected = 0, sign = False,
                ncbiTaxId = 9606, 
                inFile = 'lmpid_interactions', references = (2, ';'), header = False,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'phelm': input_formats.ReadSettings(name="phosphoELM", separator = None, nameColA=0,
                nameColB=1, nameTypeA="uniprot", nameTypeB="uniprot",
                typeA = "protein", typeB = "protein", isDirected = 1, sign = False,
                ncbiTaxId = {'col': 3, 'dict': {'Homo sapiens': 9606}}, 
                inFile = 'phelm_interactions', references = (2, ';'), header = False,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'elm': input_formats.ReadSettings(name="ELM", separator = None, nameColA=2,
                nameColB=3, nameTypeA="uniprot", nameTypeB="uniprot",
                typeA = "protein", typeB = "protein", isDirected = 0, sign = False,
                ncbiTaxId = {'A': {'col': 11, 'dict': {'"9606"(Homo sapiens)': 9606}}, 
                             'B': {'col': 12, 'dict': {'"9606"(Homo sapiens)': 9606}}},
                inFile = 'get_elm_interactions', references = (10, ','), header = False,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'domino': input_formats.ReadSettings(name="DOMINO", separator = None, nameColA=0,
                nameColB=1, nameTypeA="uniprot", nameTypeB="uniprot",
                typeA = "protein", typeB = "protein", isDirected = 0, sign = False,
                ncbiTaxId = {'A': {'col': 6, 'dict': {'9606': 9606}}, 
                             'B': {'col': 7, 'dict': {'9606': 9606}}},
                inFile = 'get_domino_interactions', references = (5, ';'), 
                header = False,
                extraEdgeAttrs={'domino_methods': (4, ';')},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'dbptm': input_formats.ReadSettings(name="dbPTM", separator = None, nameColA=0,
                nameColB=1, nameTypeA="genesymbol", nameTypeB="uniprot",
                typeA = "protein", typeB = "protein", isDirected = 1, sign = False,
                ncbiTaxId = 9606,
                inFile = 'dbptm_interactions', references = (2, ';'), header = False,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'hprd': input_formats.ReadSettings(name="HPRD", separator = None, nameColA = 6,
                nameColB = 3, nameTypeA = "genesymbol", nameTypeB = "refseqp",
                typeA = "protein", typeB = "protein", isDirected = 1, sign = False,
                ncbiTaxId = 9606,
                inFile = 'hprd_interactions', references = (10, ','), header = False,
                extraEdgeAttrs={'hprd_mechanism': 8},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'trip': input_formats.ReadSettings(name="TRIP", separator = None, nameColA = 1,
                nameColB = 0, nameTypeA = "uniprot", nameTypeB = "uniprot",
                typeA = "protein", typeB = "protein", 
                isDirected = (4, ['stimulation', 'inhibition']), 
                sign = (4, 'stimulation', 'inhibition'),
                ncbiTaxId = 9606,
                inFile = 'trip_interactions', references = (2, ';'), header = False,
                extraEdgeAttrs={'trip_methods': (3, ';')},
                extraNodeAttrsA={},
                extraNodeAttrsB={})
}

good = {
    'psite_noref': input_formats.ReadSettings(name="PhosphoSite_noref", separator="\t", 
                nameColA=0, nameColB=1, nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=True, 
                sign=False, ncbiTaxId=9606,
                inFile=os.path.join(ROOT, 'data', 'phosphosite_human_noref.csv'),
                references=False,
                extraEdgeAttrs={
                    "psite_evidences": (4, ";")},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'ppoint': input_formats.ReadSettings(name="PhosphoPoint", separator=";", 
                nameColA=1, nameColB=3,
                nameTypeA="entrez", nameTypeB="entrez",
                typeA="protein", typeB="protein", isDirected=0, header=True, ncbiTaxId=9606, 
                inFile=os.path.join(ROOT, 'data', 'phosphopoint.csv'), 
                references = False, sign=False,
                extraEdgeAttrs={
                    "phosphopoint_category": 4
                    },
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'pnetworks': input_formats.ReadSettings(name="PhosphoNetworks", 
                separator = None, nameColA=0,
                nameColB=1, nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA = "protein", typeB = "protein", isDirected = 1, sign = False,
                ncbiTaxId = 9606,
                inFile = 'pnetworks_interactions', references = False, header = False,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'mimp': input_formats.ReadSettings(name = "MIMP", 
                separator = None, nameColA = 0,
                nameColB = 1, nameTypeA = "genesymbol", nameTypeB = "genesymbol",
                typeA = "protein", typeB = "protein", isDirected = 1, sign = False,
                ncbiTaxId = 9606,
                inFile = 'mimp_interactions', references = False, header = False,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {})
}

ugly = {
    'hsn': input_formats.ReadSettings(name="Wang", separator=",", 
                nameColA = 0, nameColB = 2,
                nameTypeA="entrez", nameTypeB="entrez",
                typeA="protein", typeB="protein", isDirected = (4, ['Pos', 'Neg']), 
                sign = (4, 'Pos', 'Neg'),
                inFile = 'get_hsn',references = False, ncbiTaxId = 9606,
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'li2012': input_formats.ReadSettings(name = "Li2012", separator = False, 
                nameColA = 0, nameColB = 1,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected = 1, 
                sign = False,
                inFile = 'li2012_interactions', references = False, ncbiTaxId = 9606,
                extraEdgeAttrs = {
                    'li2012_mechanism': 3,
                    'li2012_route': 2
                },
                extraNodeAttrsA = {},
                extraNodeAttrsB = {}),
    'acsn': input_formats.ReadSettings(name="ACSN", 
                separator = None, nameColA=0,
                nameColB=2, nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA = "protein", typeB = "protein", 
                isDirected = (1, ['CATALYSIS', 'UNKNOWN_CATALYSIS', 'INHIBITION',
                    'PHYSICAL_STIMULATION', 'TRIGGER', 'activates', 
                    'UNKNOWN_POSITIVE_INFLUENCE', 'inhibits', 'MODULATION']), 
                sign = (1, ['PHYSICAL_STIMULATION', 'TRIGGER', 'activates',
                    'UNKNOWN_POSITIVE_INFLUENCE'], ['INHIBITION', 'inhibits']),
                ncbiTaxId = 9606,
                inFile = 'get_acsn', references = False, header = False,
                extraEdgeAttrs={
                    'acsn_effect': 1
                },
                extraNodeAttrsA={},
                extraNodeAttrsB={}),
    'hi3': input_formats.ReadSettings(name="HI-III", 
                separator = None, nameColA=1,
                nameColB=3, nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA = "protein", typeB = "protein", 
                isDirected = False, 
                sign = False,
                ncbiTaxId = 9606,
                inFile = '/home/denes/Dokumentumok/pw/data/hi3.tsv', 
                references = False, header = True,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {}),
    'lit13': input_formats.ReadSettings(name="Lit-BM-13", 
                separator = None, nameColA=1,
                nameColB=3, nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA = "protein", typeB = "protein", 
                isDirected = False, 
                sign = False,
                ncbiTaxId = 9606,
                inFile = '/home/denes/Dokumentumok/pw/data/Lit-BM-13.tsv', 
                references = False, header = True,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {}),
    'cpdb': input_formats.ReadSettings(name="CPDB", 
                separator = None, nameColA = 0,
                nameColB = 1, nameTypeA="uniprot-entry", nameTypeB="uniprot-entry",
                typeA = "protein", typeB = "protein", 
                isDirected = False, 
                sign = False,
                ncbiTaxId = 9606,
                inFile = 'get_cpdb', 
                references = (3, ','), header = False,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {})
}

transcription = {
    'abs': input_formats.ReadSettings(name="ABS", 
                separator = None, nameColA = 0,
                nameColB = 1, nameTypeA = "genesymbol", nameTypeB="embl_id",
                typeA = "protein", typeB = "protein", 
                isDirected = True, 
                sign = False,
                ncbiTaxId = 9606,
                inFile = 'get_abs', 
                interactionType = 'TF',
                references = False, header = False,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {}),
    'encode_dist': input_formats.ReadSettings(name = "ENCODE_distal", 
                separator = None, nameColA = 0,
                nameColB = 2, nameTypeA = "genesymbol", nameTypeB = "genesymbol",
                typeA = "protein", typeB = "protein", 
                isDirected = True, 
                sign = False,
                ncbiTaxId = 9606,
                inFile = 'http://encodenets.gersteinlab.org/enets3.Distal.txt', 
                interactionType = 'TF',
                references = False, header = False,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {}),
    'encode_prox': input_formats.ReadSettings(name = "ENCODE_proximal", 
                separator = None, nameColA = 0,
                nameColB = 2, nameTypeA = "genesymbol", nameTypeB = "genesymbol",
                typeA = "protein", typeB = "protein", 
                isDirected = True, 
                sign = False,
                ncbiTaxId = 9606,
                inFile = \
                    'http://encodenets.gersteinlab.org/enets2.Proximal_filtered.txt', 
                interactionType = 'TF',
                references = False, header = False,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {}),
    'pazar': input_formats.ReadSettings(name = "PAZAR", 
                separator = None, nameColA = 0,
                nameColB = 1, nameTypeA = "enst", nameTypeB = "ensg",
                typeA = "protein", typeB = "protein", 
                isDirected = True, 
                sign = False,
                ncbiTaxId = 9606,
                inFile = 'get_pazar', 
                interactionType = 'TF',
                references = 2, header = False,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {}),
    'htri': input_formats.ReadSettings(name = "HTRI", 
                separator = None, nameColA = 0,
                nameColB = 1, nameTypeA = "entrez", nameTypeB = "entrez",
                typeA = "protein", typeB = "protein", 
                isDirected = True, 
                sign = False,
                ncbiTaxId = 9606,
                inFile = 'get_htri', 
                interactionType = 'TF',
                references = 2, header = False,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {}),
    'oreganno': input_formats.ReadSettings(name = "ORegAnno", 
                separator = None, nameColA = 0, 
                nameColB = 1, nameTypeA = "genesymbol", nameTypeB = "genesymbol",
                typeA = "protein", typeB = "protein", 
                isDirected = True, 
                sign = False,
                ncbiTaxId = 9606,
                inFile = 'get_oreganno', 
                interactionType = 'TF',
                references = 2, header = False,
                extraEdgeAttrs = {},
                extraNodeAttrsA = {},
                extraNodeAttrsB = {})
}

negative = {
    'negatome': input_formats.ReadSettings(name="Negatome", separator="\t", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=0, 
                inFile=os.path.join(ROOT, 'data', 'negatome_manual.csv'),
                extraEdgeAttrs={
                    "references": (2, ';'),
                    "negatome_methods": (3, ';')
                    },
                extraNodeAttrsA={},
                extraNodeAttrsB={})
}

slk = input_formats.ReadSettings(name="SignaLink2", separator=",", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'slk01human.csv'),
                extraEdgeAttrs={
                    "is_stimulation": 8, 
                    "is_direct": 6,
                    "is_directed": 7,
                    "references": (9, ":")},
                extraNodeAttrsA={
                    "slk_pathways": (4, ":"),
                    "gene_name": 2},
                extraNodeAttrsB={
                    "slk_pathways": (5, ":"),
                    "gene_name": 3})

cui = input_formats.ReadSettings(name="Cui2007", separator=";", nameColA=1, nameColB=4,
                nameTypeA="entrez", nameTypeB="entrez",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'cui.sif'),
                extraEdgeAttrs={
                    "effect":6},
                extraNodeAttrsA={
                    "location": 2},
                extraNodeAttrsB={
                    "location": 5})

ca1 = input_formats.ReadSettings(name="CA1", separator=";", nameColA=1, nameColB=6,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=1, header=True,
                inFile=os.path.join(ROOT, 'data', 'ca1.csv'),
                extraEdgeAttrs={
                    "effect": 10,
                    "references": 12,
                    "type": 11},
                extraNodeAttrsA={
                    "location": 4,
                    "function": 3},
                extraNodeAttrsB={
                    "location": 9,
                    "function": 8})

ataxia = input_formats.ReadSettings(name="Ataxia", separator=";", nameColA=1, nameColB=3,
                nameTypeA="entrez", nameTypeB="entrez",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'ataxia.csv'),
                extraEdgeAttrs={
                    "sub_source":4},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

macrophage = input_formats.ReadSettings(name="Macrophage", separator=";", nameColA=0, nameColB=1,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'macrophage.sif'),
                extraEdgeAttrs={
                    "macrophage_type": 2,
                    "macrophage_location": 4},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

# this was Wang's network
hsn = input_formats.ReadSettings(name="Wang", separator=";", 
                nameColA=1, nameColB=3,
                nameTypeA="entrez", nameTypeB="entrez",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'ccmap.sif'),
                extraEdgeAttrs={
                    "ccmap_effect": 4},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

ccmap2 = input_formats.ReadSettings(name="CancerCellMap", separator="\t", 
                nameColA=3, nameColB=4,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'cell-map-edge-attributes.txt'),
                extraEdgeAttrs={
                    "references": (6, ";")},
                extraNodeAttrsA={},
                extraNodeAttrsB={},
                header=True)

spike = input_formats.ReadSettings(name="SPIKE", separator="\t", nameColA=1, nameColB=3,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'spike_hc.csv'),
                extraEdgeAttrs={
                    "references": (5, ";")},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

mppi = input_formats.ReadSettings(name="MPPI", separator="\t", nameColA=2, nameColB=6,
                nameTypeA="uniprot", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'mppi_human.csv'),
                extraEdgeAttrs={
                    "references": (0, ";"),
                    "mppi_evidences": (1, ";")},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

psite = input_formats.ReadSettings(name="PhosphoSite", separator="\t", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'phosphosite_human_hc.csv'),
                extraEdgeAttrs={
                    "references": (5, ";"),
                    "psite_evidences": (4, ";")},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

panther = input_formats.ReadSettings(name="Panther", separator=";", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'panther1.csv'),
                extraEdgeAttrs={
                    "panther_type": 3},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

tlr = input_formats.ReadSettings(name="TLR", separator=";", nameColA=0, nameColB=1,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected=0, 
                inFile=os.path.join(ROOT, 'data', 'tlr.csv'),
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

trip = input_formats.ReadSettings(name="TRIP", separator=";", nameColA=0, nameColB=1,
                nameTypeA="genesymbol", nameTypeB="genesymbol",
                typeA="protein", typeB="protein", isDirected=0, 
                inFile=os.path.join(ROOT, 'data', 'trip.sif'),
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

alz = input_formats.ReadSettings(name="AlzPathway", separator="\t", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=0, 
                inFile=os.path.join(ROOT, 'data', 'alzpw-ppi.csv'),
                extraEdgeAttrs={
                    "refrences": (8, ";")
                    },
                extraNodeAttrsA={},
                extraNodeAttrsB={})

innatedb = input_formats.ReadSettings(name="InnateDB", separator=";", nameColA=0, nameColB=2,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=0, 
                inFile=os.path.join(ROOT, 'data', 'innatedb.csv'),
                extraEdgeAttrs={
                    "references": (4, ":")
                },
                extraNodeAttrsA={},
                extraNodeAttrsB={})

depod = input_formats.ReadSettings(name="DEPOD", separator=";", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=0, 
                inFile=os.path.join(ROOT, 'data', 'depod.csv'),
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

pp = input_formats.ReadSettings(name="PhosphoPoint", separator=";", nameColA=1, nameColB=3,
                nameTypeA="entrez", nameTypeB="entrez",
                typeA="protein", typeB="protein", isDirected=0, header=True, 
                inFile=os.path.join(ROOT, 'data', 'phosphopoint.csv'),
                extraEdgeAttrs={"phosphopoint_category":4},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

arn = input_formats.ReadSettings(name="ARN", separator=",", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=0, 
                inFile=os.path.join(ROOT, 'data', 'arn.csv'),
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

nrf2 = input_formats.ReadSettings(name="NRF2ome", separator=",", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=0, 
                inFile=os.path.join(ROOT, 'data', 'nrf2ome.csv'),
                extraEdgeAttrs={},
                extraNodeAttrsA={},
                extraNodeAttrsB={})

netpath = input_formats.ReadSettings(name="NetPath", separator=";", nameColA=0, nameColB=1,
                nameTypeA="hgnc", nameTypeB="hgnc",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'netpath.csv'),
                extraEdgeAttrs={},
                extraNodeAttrsA={
                    "netpath_pathways": (2, ",")},
                extraNodeAttrsB={
                    "netpath_pathways": (2, ",")})

biocarta = input_formats.ReadSettings(name="BioCarta", separator=";", nameColA=0, nameColB=2,
                nameTypeA="entrez", nameTypeB="entrez",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'biocarta-pid.csv'),
                extraEdgeAttrs={},
                extraNodeAttrsA={
                    "biocarta_pathways": (4, ",")},
                extraNodeAttrsB={
                    "biocarta_pathways": (4, ",")})

nci_pid = input_formats.ReadSettings(name="NCI-PID", separator=";", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'nci-pid.csv'),
                extraEdgeAttrs={},
                extraNodeAttrsA={
                    "nci_pid_pathways": (2, ",")},
                extraNodeAttrsB={
                    "nci_pid_pathways": (2, ",")})

reactome = input_formats.ReadSettings(name="Reactome", separator=";", nameColA=0, nameColB=1,
                nameTypeA="uniprot", nameTypeB="uniprot",
                typeA="protein", typeB="protein", isDirected=1, 
                inFile=os.path.join(ROOT, 'data', 'reactome-pid.csv'),
                extraEdgeAttrs={},
                extraNodeAttrsA={
                    "reactome_pathways": (2, ",")},
                extraNodeAttrsB={
                    "reactome_pathways": (2, ",")})

gdsc_comp_target = input_formats.ReadSettings(name="GDSC", 
                separator=";", nameColA=1, nameColB=0,
                nameTypeA="pubchem", nameTypeB="genesymbol",
                typeA="drug", typeB="protein", isDirected=1, 
                inFile="gdsc.sif",
                extraEdgeAttrs={},
                extraNodeAttrsA={
                    "gene_name": 2},
                extraNodeAttrsB={})

gdsc_lst = input_formats.ReadList(name="GDSC", separator=";", nameCol=0,
                nameType="genesymbol", typ="protein",
                inFile=os.path.join(ROOT, 'data', 'gdsc.sif'),
                extraAttrs={'drugs': 2})

gdsc_lst = input_formats.ReadList(name="atg", separator=";", nameCol=0,
                nameType="genesymbol", typ="protein",
                inFile=os.path.join(ROOT, 'data', 'autophagy.list'),
                extraAttrs={'drugs': 2})

cgc = input_formats.ReadList(name="CancerGeneCensus", separator="|", nameCol=2,
                nameType="entrez", typ="protein",
                inFile=os.path.join(ROOT, 'data', 'cancer_gene_census.csv'),
                extraAttrs={})

intogene_cancer = input_formats.ReadList(name="Intogene", separator="\t", nameCol=1,
               nameType="genesymbol", typ="protein",
              inFile=os.path.join(ROOT, 'data', 'intogene_cancerdrivers.tsv'),
              extraAttrs={})

aidan_list = input_formats.ReadList(name="aidan_list", separator=";", nameCol=0,
                nameType="uniprot", typ="protein",
                inFile=os.path.join(ROOT, 'data', 'aidan_list_uniprot'),
                extraAttrs={})
