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

import os

import pypath.share.common as common
import pypath.share.cache as cache
import pypath.share.settings as settings

urls = {
    'uniprot_pdb': {
        'label': 'Getting PDB IDs of 3D structures for UniProtIDs',
        'url': 'https://ftp.uniprot.org/pub/databases/uniprot/'
            'current_release/knowledgebase/complete/docs/pdbtosp.txt'
    },
    'uniprot_basic': {
        'label': 'URL for UniProt queries',
        'url': 'https://rest.uniprot.org/uniprotkb/stream',
        'lists': 'https://legacy.uniprot.org/uploadlists/',
        'datasheet': 'https://rest.uniprot.org/uniprotkb/%s.txt',
        'history': 'https://rest.uniprot.org/unisave/%s?format=tsv&version=*',
        'deleted_sp': 'https://ftp.expasy.org/databases/uniprot/'
            'current_release/knowledgebase/complete/docs/delac_sp.txt',
        'deleted_tr': 'https://ftp.expasy.org/databases/uniprot/'
            'current_release/knowledgebase/complete/docs/delac_tr.txt.gz',
        'speindex_old': 'ftp://ftp.uniprot.org/pub/databases/uniprot/'
            'knowledgebase/docs/speindex.txt',
        'speindex': 'https://ftp.expasy.org/databases/uniprot/current_release/'
            'knowledgebase/complete/docs/speindex.txt',
        'taxids_old': 'https://www.uniprot.org/taxonomy/?query='
            '*&format=tab&force=true&columns=id&compress=yes',
        'taxids': 'https://rest.uniprot.org/taxonomy/stream?compressed=true'
            '&fields=id%2Ccommon_name%2Cscientific_name'
            '&format=tsv&query=%28%2A%29',
        'speclist': 'https://ftp.uniprot.org/pub/databases/uniprot/'
            'current_release/knowledgebase/complete/docs/speclist.txt',
    },
    'uniprot_idmapping': {
        'run': 'https://rest.uniprot.org/idmapping/run',
        'poll': 'https://rest.uniprot.org/idmapping/status/%s',
        'details': 'https://rest.uniprot.org/idmapping/details/%s',
        'fields': 'https://rest.uniprot.org/configure/idmapping/fields',
    },
    'corum': {
        'label':
        'CORUM is a database of protein complexes, downloadable in csv format',
        'url_old':
        'http://mips.helmholtz-muenchen.de/genre/proj/corum/allComplexes.csv',
        'url': 'http://mips.helmholtz-muenchen.de/corum/download/'
            'allComplexes.txt.zip',
        'url_rescued':
            'https://rescued.omnipathdb.org/CORUM_allComplexes.txt.zip',
    },
    'pfam_pdb': {
        'label': 'PDB-Pfam mapping and names of Pfam domains',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/'
            'pdb_pfam_mapping.txt',
    },
    'pfam_up': {
        'label': 'Mapping Pfam regions to UniProt',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/'
            'Pfam-A.regions.tsv.gz',
    },
    '3dcomplex_contact': {
        'label': 'This file contains the topology '
            'definition of the complexes',
        'url_old': 'http://shmoo.weizmann.ac.il/elevy/3dcomplexV4/dataV4/'
            'contactDefinition.txt',
        'url': 'http://shmoo.weizmann.ac.il/elevy/3dcomplexV6/dataV6/'
            'contactDefinition.txt',
    },
    '3dcomplex_correspondancy': {
        'label': 'This is the dictionary of chain names',
        'url_old': 'http://shmoo.weizmann.ac.il/elevy/3dcomplexV4/dataV4/'
            'pdb_chain_corresV2.txt',
        'url': 'http://shmoo.weizmann.ac.il/elevy/3dcomplexV6/dataV6/'
            'chainDefinition.txt',
    },
    'pdb_chains': {
        'label': 'Corresponding UniProt IDs and residue numbers '
            'for each chain in PDB structures',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/'
            'tsv/pdb_chain_uniprot.tsv.gz',
    },
    'complex_portal': {
        'label': 'Complexes curated by IntAct',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/'
            'current/psi25/',
    },
    'pisa_interfaces': {
        'label': 'Base URL for download interface data from PDBe PISA',
        'url': 'http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?'
    },
    'catalytic_sites': {
        'label': 'Catalytic Site Atlas',
        'url': 'http://www.ebi.ac.uk/thornton-'
            'srv/databases/CSA/downloads/CSA_2_0_121113.txt',
    },
    '3did_ddi': {
        'label': 'Domain-domain interactions derived from 3D structures',
        'url': 'http://3did.irbbarcelona.org/download/current/3did_flat.gz'
    },
    '3did_dmi': {
        'label': 'Domain-motif interactions from 3DID',
        'url': 'http://3did.irbbarcelona.org/download/current/'
            '3did_dmi_flat.gz',
    },
    'swissprot_full': {
        'label': 'Full UniProt/Swissprot database',
        'url': 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/'
            'current_release/uniprot_sprot.dat.gz',
    },
    'instruct_human': {
        'label':
        'Protein interactome networks annotated to 3D structural resolution',
        'url': 'http://instruct.yulab.org/download/sapiens.sin',
    },
    'i3d_human': {
        'label': 'Interactome3D representative dataset for human proteins',
        'url': 'http://interactome3d.irbbarcelona.org/user_data/human/'
            'download/representative/interactions.dat',
    },
    'instruct_offsets': {
        'label': 'Offsets between PDB chains and UniProt sequences',
        'url': 'http://instruct.yulab.org/download/indexing_uniprot2pdb.txt'
    },
    'compleat': {
        'label': 'Curated and inferred complexes from multiple databases',
        'old_url': 'http://www.flyrnai.org/compleat/ComplexDownload'
        '?requestType=complexDownload&org=Human',
        'rescued': 'https://rescued.omnipathdb.org/compleat.tsv',
    },
    'switches.elm': {
        'label': 'Curated data on molecular switches in cellular regulation',
        'url': 'http://switches.elm.eu.org/downloads/switches.ELM-v1.txt'
    },
    'ols': {
        'label': 'WSDL interface for the Ontology Lookup Service',
        'url_wsdl': 'http://www.ebi.ac.uk/ontology-lookup/OntologyQuery.wsdl',
        'url': 'https://www.ebi.ac.uk/ols/api/ontologies',
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
        'url_old': 'http://pdb.org/pdb/rest/das/pdb_uniprot_mapping/'
            'alignment?query=',
        'url': 'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/',
    },
    'pepcyber': {
        'label': 'MySQL injection to PEPCyber website :)',
        'url': 'http://www.pepcyber.org/PPEP/search_result.php'
            '?domain=Any&ppbd_symbol=Any&search_field=symbol&'
            'query_value=%27+OR+1&binding_sequence=&go_id=Any&'
            'Submit=Search',
        'rescued': 'https://rescued.omnipathdb.org/'
            'pepcyber_search_result.html',
        'details': 'http://www.pepcyber.org/PPEP/idetail.php?iid=%u',
        'details_rescued': 'https://rescued.omnipathdb.org/pepcyber_details/'
            'pepcyber_%u.html',
    },
    'pdzbase': {
        'label': 'Manually curated interactions of PDZ domain proteins',
        'url_rescued': 'http://rescued.omnipathdb.org/PDZBase.html',
        'url': 'http://abc.med.cornell.edu/pdzbase/allinteractions',
    },
    'pdz_details': {
        'label': 'Details of interactions in PDZbase',
        'url': 'http://abc.med.cornell.edu/pdzbase/interaction_detail/%u'
    },
    'psite_reg': {
        'label': 'PhosphoSite annotated regulatory sites',
        'url_old': 'http://www.phosphosite.org/downloads/Regulatory_sites.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'\
            'Regulatory_sites.gz',
    },
    'psite_bp': {
        'label': 'PhosphoSite kinase substrates in BioPAX format',
        'url_old': 'http://www.phosphosite.org/downloads/'\
            'Kinase_substrates.owl.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'\
            'Kinase_substrates.owl.gz',
    },
    'psite_ac': {
        'label': 'PhosphoSite acetylation sites',
        'url_old': 'http://www.phosphosite.org/downloads/'\
            'Acetylation_site_dataset.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'\
            'Acetylation_site_dataset.gz',
    },
    'psite_kin': {
        'label': 'PhosphoSite kinase-substrate interactions',
        'url_old': 'http://www.phosphosite.org/downloads/'\
            'Kinase_Substrate_Dataset.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'\
            'Kinase_Substrate_Dataset.gz',
    },
    'psite_me': {
        'label': 'PhosphoSite methylation sites',
        'url_old': 'http://www.phosphosite.org/downloads/'\
            'Methylation_site_dataset.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'
            'Methylation_site_dataset.gz',
    },
    'psite_ga': {
        'label': 'PhosphoSite O-GalNAc sites',
        'url_old': 'http://www.phosphosite.org/downloads/'\
            'O-GalNAc_site_dataset.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'
            'O-GalNAc_site_dataset.gz',
    },
    'psite_gl': {
        'label': 'PhosphoSite O-GlcNAc sites',
        'url_old': 'http://www.phosphosite.org/downloads/'\
            'O-GlcNAc_site_dataset.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'\
            'O-GlcNAc_site_dataset.gz',
    },
    'psite_p': {
        'label': 'PhosphoSite phosphorylation sites',
        'url_old': 'http://www.phosphosite.org/downloads/'\
            'Phosphorylation_site_dataset.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'\
            'Phosphorylation_site_dataset.gz',
    },
    'psite_sm': {
        'label': 'Sumoylation sites',
        'url_old': 'http://www.phosphosite.org/downloads/'\
            'Sumoylation_site_dataset.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'\
            'Sumoylation_site_dataset.gz',
    },
    'psite_ub': {
        'label': 'Ubiquitination sites',
        'url_old': 'http://www.phosphosite.org/downloads/'\
            'Ubiquitination_site_dataset.gz',
        'url': 'https://rescued.omnipathdb.org/phosphosite/'\
            'Ubiquitination_site_dataset.gz',
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
        '2009-10-22/2009-10-22-domino-full-binary.mitab26',
        'rescued': 'http://rescued.omnipathdb.org/'\
            'domino_2009-10-22-full-binary.mitab26.gz',
    },
    'hprd_all': {
        'label': 'HPRD all data in flat files',
        'url': 'http://www.hprd.org/RELEASE9/HPRD_FLAT_FILES_041310.tar.gz',
        'url_rescued': 'https://rescued.omnipathdb.org/'
            'HPRD_FLAT_FILES_041310.tar.gz',
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
        'url_old': 'http://elm.eu.org/elms/browse_instances.tsv?q=*',
        'url': 'http://elm.eu.org/instances.tsv?q=*',
    },
    'elm_class': {
        'label': 'List of ELM classes',
        'url_old': 'http://elm.eu.org/elms/browse_elms.tsv',
        'url': 'http://elm.eu.org/elms/elms_index.tsv',
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
            'http://dbptm.mbc.nctu.edu.tw/download/Nitrosylation.tgz',
        ],
        'url': 'http://dbptm.mbc.nctu.edu.tw/download/experiment/%s.txt.gz',
        'datasets': [
            'Phosphorylation',
            'Acetylation',
            'Ubiquitination',
            'Succinylation',
            'Methylation',
            'Malonylation',
            'N-linkedGlycosylation',
            'O-linkedGlycosylation',
            'Sumoylation',
            'S-nitrosylation',
            'Glutathionylation',
            'Amidation',
            'Hydroxylation',
            'PyrrolidoneCarboxylicAcid',
            'Glutarylation',
            'Palmitoylation',
            'Gamma-carboxyglutamicAcid',
            'Crotonylation',
            'Oxidation',
            'Myristoylation',
            'C-linkedGlycosylation',
            'Sulfation',
            'Formylation',
            'Citrullination',
            'GPI-anchor',
            'Nitration',
            'S-diacylglycerol',
            'Carboxylation',
            'Lipoylation',
            'Carbamidation',
            'Neddylation',
            'Pyruvate',
            'S-linkedGlycosylation',
        ],
        'old_table': 'http://rescued.omnipathdb.org/old_dbptm.tab',
    },
    'phosnw': {
        'label': 'Human kinase-substrate relationships and phosphosites',
        'url': 'http://phosphonetworks.org/download/highResolutionNetwork.csv'
    },
    'kegg_pws': {
        'label': 'KEGG Pathways and KEGG MEDICUS',
        'list_url': 'http://www.genome.jp/kegg/pathway.html',
        'kgml_url':
        'http://www.kegg.jp/kegg-bin/download?entry=%s&format=kgml',
        'biopax_l3': 'http://www.pathwaycommons.org/archives/PC2/'
            'v8/PathwayCommons.8.kegg.BIOPAX.owl.gz',
        'kgml_url_2': 'http://rest.kegg.jp/get/%s/kgml',
        'medicus': 'ftp://ftp.genome.jp/pub/kegg/medicus/network/network',
        'dbget': 'https://www.kegg.jp/dbget-bin/www_bget?%s',
        'pw_annot': 'https://www.pathwaycommons.org/archives/PC2/v12/'
            'PathwayCommons12.kegg.uniprot.gmt.gz',
    },
    'kegg_api': {
        'label': 'KEGG API',
        'url': 'https://rest.kegg.jp/%s'
    },
    'pathbank': {
        'label': 'A functional annotation and process description network '
            'database.',
        'pw_annot': 'https://pathbank.org/downloads/'
            'pathbank_all_proteins.csv.zip',
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
        'url_old': 'ftp://ftp.uniprot.org/pub/databases/uniprot/'
            'current_release/knowledgebase/complete/'
            'uniprot_sprot_varsplic.fasta.gz',
        'url': 'ftp://ftp.expasy.org/databases/uniprot/current_release/'
            'knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz',
    },
    'kinclass': {
        'label': 'Kinase families and groups',
        'url': 'http://kinase.com/static/colt/data/human/kinome/'
            'tables/Table%20S2.txt',
        'rescued': 'https://rescued.omnipathdb.org/kinase.com__Table_S2.txt',
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
        'url_old': 'ftp://ftp.uniprot.org/pub/databases/uniprot/'
            'current_release/knowledgebase/complete/docs/sec_ac.txt',
        'url': 'https://ftp.expasy.org/databases/uniprot/current_release/'
            'knowledgebase/complete/docs/sec_ac.txt',
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
        'url': 'http://www.pazar.info/tftargets/tftargets.zip',
        'url_rescued': 'http://rescued.omnipathdb.org/pazar_tftargets.zip',
    },
    'htri': {
        'label': 'TF-target gene lists from HTRI',
        'url': 'https://rescued.omnipathdb.org/HTRIdb.csv',
        'url_old1': 'http://www.lbbc.ibb.unesp.br/htri/datasHTRIdb?down=3',
        'url_old2': 'http://www.lbbc.ibb.unesp.br/htri/consulta?'\
            'type=1&all=true&down=3&iconss1.x=57&iconss1.y=48',
        'init_url': 'http://www.lbbc.ibb.unesp.br/htri/consulta?'\
            'type=1&all=true&down=3',
    },
    'oreganno_old': {
        'label': 'TF-target gene lists from ORegAnno, previous version',
        'url': 'http://www.oreganno.org/oregano/htdocs'
        '/data/oreganno_UCSC_08Nov10.txt.gz'
    },
    'oreganno': {
        'label': 'TF-target gene lists from ORegAnno',
        'url_old': 'http://py-gi1.stanford.edu:8080/oregano/htdocs/'
        'downloads/ORegAnno_Combined_2015.09.16.tsv',
        'url': 'http://www.oreganno.org/dump/'\
            'ORegAnno_Combined_2016.01.19.tsv',
    },
    'cpdb': {
        'label': 'All human interactions from ConsensusPathDB',
        'url':
        'http://cpdb.molgen.mpg.de/download/ConsensusPathDB_human_PPI.gz'
    },
    'goa': {
        'label': 'UniProt GO annotations from GOA',
        'go_url': 'http://geneontology.org/gene-associations/goa_%s.gaf.gz',
        'ebi_url_old': 'ftp://ftp.ebi.ac.uk/pub/'
            'databases/GO/goa/%s/goa_%s.gaf.gz',
        'ebi_url': 'https://ftp.ebi.ac.uk/pub/databases/GO/'
            'goa/%s/goa_%s.gaf.gz',
    },
    'quickgo_old': {
        'label': 'UniProt GO annotations from QuickGO',
        'url': 'http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&'
        'limit=-1%s&termUse=%s&tax=%u&col=proteinID,goID,goName,aspect'
    },
    'quickgo': {
        'label': 'UniProt GO annotations from QuickGO',
        'url': 'https://www.ebi.ac.uk/QuickGO/services/annotation/'\
        'downloadSearch?includeFields=goName,name&taxonId=%u'
        '&aspect=%s%s'
    },
    'quickgo_rest': {
        'label': 'REST API of the QuickGO Gene Ontology server',
        'graph': 'https://www.ebi.ac.uk/QuickGO/services/ontology/'\
            'go/terms/graph?startIds=%s,relations=%s',
        'terms': 'https://www.ebi.ac.uk/QuickGO/services/ontology/'\
            'go/terms?page=%u',
        'annot': 'https://www.ebi.ac.uk/QuickGO/services/annotation/'\
            'downloadSearch?'\
            'includeFields=goName&'\
            'includeFields=taxonName&'\
            'includeFields=name&'\
            'aspect=%s&'\
            'goUsageRelationships=%s&'\
            'taxonId=%u&'\
            'limit=100&'\
            'page=%u',
        'desc': 'https://www.ebi.ac.uk/QuickGO/services/ontology/'\
            'go/terms/%s/descendants%s',
    },
    'quickgo_desc': {
        'label': 'UniProt GO annotations from QuickGO',
        'url': 'https://www.ebi.ac.uk/QuickGO/services/annotation/'\
        'downloadSearch?includeFields=goName,name&taxonId=%u'\
        '&goUsage=descendants&goId=%s'
    },
    'goose': {
        'label': 'Legacy MySQL query interface of AmiGO',
        'url': 'http://amigo.geneontology.org/goose?query=%s'\
            '&mirror=bbop&limit=0&format=text',
    },
    'golr': {
        'label': 'AmiGO Solr service URL',
        'url': 'http://golr.berkeleybop.org/select%s'
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
    'netpath_pw': {
        'label': 'NetPath pathway page',
        'mainpage': 'http://mirror.netpath.org/index.html',
        'url': 'http://mirror.netpath.org/pathways?path_id=NetPath_%u'
    },
    'netpath_psimi': {
        'label': 'Batch download of NetPath pathways in PSI-MI format',
        'url': 'http://rescued.omnipathdb.org/NetPath_PSI-MI.zip'
    },
    'netpath_bp': {
        'label': 'Netpath pathways one by one in BioPAX level 3 format',
        'biopax_l3': 'http://www.netpath.org/data/biopax/NetPath_%u.owl'
    },
    'proteomemap': {
        'label':
        'Human Proteome Map: Mass-spec expression data in '\
            'healthy human tissues',
        'url':
        'http://www.humanproteomemap.org/Download_HPM/HPM_protein_level_'
            'expression_matrix_Kim_et_al_052914.csv'
    },
    'proteinatlas': {
        'label': 'Human Protein Atlas: Immuncytochemistry expression data in '
            'healthy human cells or cancer cells or cell lines',
        'normal': 'http://www.proteinatlas.org/download/'\
            'normal_tissue.tsv.zip',
        'pathology': 'http://www.proteinatlas.org/download/pathology.tsv.zip',
        'subcell': 'https://www.proteinatlas.org/download/'\
            'subcellular_location.tsv.zip',
        'secretome_old': 'https://stke.sciencemag.org/highwire/filestream/'\
            '216463/field_highwire_adjunct_files/1/aaz0274_Data_file_S2.xlsx',
        'secretome': 'https://www.science.org/action/downloadSupplement?'\
            'doi=10.1126%2Fscisignal.aaz0274&'\
            'file=aaz0274_data_file_s2.xlsx',
        'ig_genes': 'https://stke.sciencemag.org/highwire/filestream/'\
            '216463/field_highwire_adjunct_files/0/aaz0274_Data_file_S1.xlsx',

    },
    'lincs-compounds': {
        'label':
        'List of small molecules in LINCS with synonyms and PubChem IDs',
        'url': 'http://lincs.hms.harvard.edu/db/datasets/20000/smallmolecules'
            '?search=&output_type=.csv'
    },
    'eutils': {
        'label': 'NCBI E-Utils API',
        'esummary': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',
        'pmc-idconv': (
            'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/'
            'v1.0/?ids=%s&format=json'
        ),
    },
    'pubmed': {
        'label': 'PubMed baseurl',
        'url': 'http://www.ncbi.nlm.nih.gov/pubmed/%s'
    },
    'hpmr': {
        'label': 'Human Plasma Membrane Receptome',
        'url':
            'http://www.receptome.org/SearchDB/findGenes.asp?textName='
    },
    'tfcensus': {
        'label': 'A census of human transcription factors (Vaquerizas 2009)',
        'url': (
            r'https://static-content.springer.com/esm/'
            r'art%3A10.1038%2Fnrg2538/MediaObjects/'
            r'41576_2009_BFnrg2538_MOESM6_ESM.txt'
        )
    },
    'tfcheckp': {
        'label': 'A comprehensive list of transcription factors',
        'url':
        'http://www.tfcheckpoint.org/data/TFCheckpoint_download_180515.txt'
    },
    'gtp': {
        'label': 'Guide to Pharmacology: ligand-receptor interactions',
        'url': 'http://www.guidetopharmacology.org/DATA/%s.csv',
    },
    'giant': {
        'label': 'JSON network from GIANT by gene based query',
        'init_url': 'http://giant.princeton.edu/',
        'url':
        'http://giant.princeton.edu/contexts/integration/'\
            'query/%u/?%s&prior=0.1'
    },
    'msigdb': {
        'label': 'Molecular Signatures Database',
        'login1': 'https://www.gsea-msigdb.org/gsea/login',
        'login2': 'https://www.gsea-msigdb.org/gsea/j_spring_security_check',
        'url': 'https://www.gsea-msigdb.org/gsea/msigdb/'
            'download_file.jsp?filePath=/msigdb/release/'
            '%s.%s/%s.v%s.%s.%s.gmt',
        'coll': 'https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp',
        'url_stem': 'https://www.gsea-msigdb.org/gsea/%s',
        'one_set': 'hhttps://www.gsea-msigdb.org/gsea/msigdb/'
            'download_geneset.jsp?geneSetName=%s&fileType=txt',
        'url_all': 'https://www.gsea-msigdb.org/gsea/msigdb/'\
            'download_file.jsp?filePath=/msigdb/release/7.5.1/'\
            'msigdb_v7.5.1_files_to_download_locally.zip',
    },
    'msigdb_old': {
        'label': 'Molecular Signatures Database',
        'login1': 'http://software.broadinstitute.org/gsea/login.jsp',
        'login2':
        'http://software.broadinstitute.org/gsea/j_spring_security_check',
        'url': 'http://software.broadinstitute.org/gsea/msigdb/'
        'download_file.jsp?filePath=/resources/msigdb/7.0/%s.v7.0.%s.gmt',
        'coll': 'http://www.broadinstitute.org/gsea/msigdb/collections.jsp',
        'url_stem': 'http://www.broadinstitute.org/gsea/%s',
        'one_set':
        'http://www.broadinstitute.org/gsea/msigdb/download_geneset.jsp?'
        'geneSetName=%s&fileType=txt',
        'url_all': 'http://software.broadinstitute.org/gsea/msigdb/'\
            'download_file.jsp?filePath=/resources/msigdb/7.0/'\
            'msigdb_v7.0_files_to_download_locally.zip',
    },
    'disgenet': {
        'api_url': 'https://www.disgenet.org/api',
        'disease_id_mappings': 'https://www.disgenet.org/static/disgenet_ap1/files/downloads/disease_mappings.tsv.gz',
        'variant_gene_mappings': 'https://www.disgenet.org/static/disgenet_ap1/files/downloads/variant_to_gene_mappings.tsv.gz',
        'annotation_label': 'Disease-gene associations',
        'annotations': 'http://www.disgenet.org/static/disgenet_ap1/files/'\
            'downloads/%s_gene_disease_associations.tsv.gz',
        'annotation_datasets': ['curated', 'literature', 'befree', 'all'],
    },
    'hsn': {
        'label': 'The Wang Human Signaling Network',
        'rescued': (
            'http://rescued.omnipathdb.org/WangLab_HumanSignalingNet_v6.csv'
        ),
        'researchgate': 'https://www.researchgate.net/profile/Edwin-Wang-6/'
            'publication/264488606_Human_Signaling_Network_Version_6/data/'
            '53e102ed0cf2235f35271ea3/'
            'HuamnSignalingNet-v6.csv?origin=publication_detail',
    },
    'li2012': {
        'label': 'Human Phosphotyrosine Signaling Network',
        'url': 'http://genome.cshlp.org/content/suppl/2011/12/27/'
        'gr.128819.111.DC1/Supplementary_files_S1-S5.xls'
    },
    'trip': {
        'label': 'The TRP channel database',
        'base': 'http://www.trpchannel.org/',
        'base_rescued': 'https://rescued.omnipathdb.org/trip/main.html',
        'show': 'http://www.trpchannel.org/proteins/show?id=%s',
        'show_rescued': 'https://rescued.omnipathdb.org/trip/show-%s.html',
        'intr': 'http://www.trpchannel.org/interactions/show?trp=%s&'
        'interactor=%s',
        'url': 'http://www.trpchannel.org/20141116.csv',
        'json': 'http://www.trpchannel.org/proteins/getjson'
    },
    'signor': {
        'label': 'SIGNOR pathways',
        'list_url': 'https://signor.uniroma2.it/downloads.php',
        'base_url': 'https://signor.uniroma2.it/%s',
        'all_url': 'https://signor.uniroma2.it/DownloadServlet?action=all',
        'all_url_new': 'https://signor.uniroma2.it/download_entity.php',
        'complexes': 'https://signor.uniroma2.it/download_complexes.php',
    },
    'hiii14': {
        'label': 'Rolland et al 2014, Human Interactome II',
        'url':
        'http://www.cell.com/cms/attachment/2021150198/2041343422/mmc3.xlsx',
        'article_url': 'https://www.cell.com/cell/fulltext/'
            'S0092-8674(14)01422-6',
    },
    'kinome': {
        'label': 'List of human kinases',
        'url': 'http://kinase.com/static/colt/data/human/kinome/tables/'
            'Kincat_Hsap.08.02.xls',
        'rescued': 'https://rescued.omnipathdb.org/Kincat_Hsap.08.02.xls',
    },
    'phosphatome': {
        'label': 'List of human and other phosphatases',
        'webpage': 'http://phosphatome.net/3.0',
        'url_from_page': 'https://www.science.org/doi/suppl/10.1126/'\
            'scisignal.aag1796/suppl_file/aag1796_tables_s1_to_s23.zip',
        'url': 'https://www.science.org/action/downloadSupplement?'
            'doi=10.1126%2Fscisignal.aag1796&'
            'file=aag1796_tables_s1_to_s23.zip',
    },
    'dgidb': {
        'label': 'Druggable genes compiled from multiple resources',
        'main_url': 'http://dgidb.genome.wustl.edu/search_categories',
        'url': 'http://dgidb.genome.wustl.edu/druggable_gene_categories/%s?'
            'sources=BaderLabGenes,CarisMolecularIntelligence,'
            'FoundationOneGenes,GO,GuideToPharmacologyGenes,HopkinsGroom,'
            'MskImpact,RussLampel,dGene',
        'categories': 'https://www.dgidb.org/data/monthly_tsvs/'
            '2021-Jan/categories.tsv',
        'interactions': 'https://www.dgidb.org/data/monthly_tsvs/'
            '2022-Feb/interactions.tsv',
    },
    'acsn': {
        'label': 'Atlas of Cancer Signaling Networks',
        'biopax_l3': 'https://acsn.curie.fr/files/acsn_v1.1.owl',
        'sif': 'https://acsn.curie.fr/files/acsn_ppi.sif',
        'ppi': 'http://rescued.omnipathdb.org/acsn_ppi.txt',
        'names': 'http://rescued.omnipathdb.org/acsn_names.gmt',
    },
    'nci-pid': {
        'label': 'National Cancer Institute -- Pathway Interaction Database',
        'biopax_l3':
        'ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_3/'\
            'NCI-Nature_Curated.bp3.owl.gz'
    },
    'panther': {
        'label': 'PANTHER pathways',
        'biopax_l3':
        'ftp://ftp.pantherdb.org//pathway/current_release/BioPAX.tar.gz'
    },
    'laudanna': {
        'label':
        'Directionality and effect information inferred from multiple'\
            ' sources by Laudanna Lab',
        'sigflow':
        'http://dp.univr.it/~laudanna/LCTST/downloads/files/SignalingFlow.EA',
        'sigdir':
        'http://dp.univr.it/~laudanna/LCTST/downloads/files/'\
            'SignalingDirection.EA',
        'sigdir_rescued': 'https://rescued.omnipathdb.org/'\
            'Laudanna_SignalingDirection.EA',
        'sigflow_rescued': 'https://rescued.omnipathdb.org/'\
            'Laudanna_SignalingFlow.EA',
    },
    'biogrid': {
        'label': 'BioGRID version 2 tab format',
        'mv': 'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/'\
            'BIOGRID-MV-Physical-LATEST.tab2.zip',
        'all': 'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/'\
            'BIOGRID-ALL-LATEST.tab2.zip'
    },
    'wang': {
        'label':
        'Human Signaling Network compiled by Wang Group, version 6 (2014)',
        'old_url': 'http://www.bri.nrc.ca/wang/cancerMap/'
            'HumanSignalingNetwork_v6.csv',
        'rescued': 'https://rescued.omnipathdb.org/'
            'WangLab_HumanSignalingNetwork_v6-format2.csv',
        'cui': 'https://www.embopress.org/action/downloadSupplement?'
            'doi=10.1038%2Fmsb4100200&file=msb4100200-sup-0010.xls',
        'cui_init': 'https://www.embopress.org/doi/full/10.1038/msb4100200',
    },
    'graphviz': {
        'label': 'List of Graphviz attributes',
        'url': 'http://www.graphviz.org/doc/info/attrs.html'
    },
    'hid': {
        'label': 'CCSB Human Interactome Project',
        'lit-bm-13': 'http://interactome.dfci.harvard.edu/H_sapiens/'
            'download/Lit-BM-13.tsv',
        'hi-ii-14':
        'http://interactome.dfci.harvard.edu/H_sapiens/download/HI-II-14.psi',
        'hi-i-05':
        'http://interactome.dfci.harvard.edu/H_sapiens/download/HI-I-05.psi',
        'lit-bm-17': 'https://rescued.omnipathdb.org/LitBM-17.psi',
        'lit-bm': 'http://www.interactome-atlas.org/data/Lit-BM.tsv',
        'huri': 'http://www.interactome-atlas.org/data/HuRI.psi',
        'hi-union': 'http://www.interactome-atlas.org/data/HI-union.psi',
        'yang-2016': 'http://www.interactome-atlas.org/data/Yang-16.psi',
        'hi-i-05-pmi': 'http://www.interactome-atlas.org/data/H-I-05.psi',
        'hi-ii-14-pmi': 'http://www.interactome-atlas.org/data/HI-II-14.psi',
        'yu-2011': 'http://www.interactome-atlas.org/data/Yu-11.psi',
    },
    'ca1': {
        'label': 'Supplementary Online Materials for Ma\'ayan 2005',
        'url_old': 'http://science.sciencemag.org/highwire/filestream/586741/'
            'field_highwire_adjunct_files/1/Maayan_SOM_External_Files.zip',
        'url': 'https://www.science.org/action/downloadSupplement?'
            'doi=10.1126%2Fscience.1108876&'
            'file=maayan_som_external_files.zip',
    },
    'ccmap': {
        'label': 'Cancer Cell Map from PathwayCommons 2011 snapshot',
        'edges':
        'http://www.pathwaycommons.org/archives/PC1/last_release-2011/'\
            'tab_delim_network/by_source/cell-map-edge-attributes.txt.zip',
        'nodes':
        'http://www.pathwaycommons.org/archives/PC1/last_release-2011/'\
            'tab_delim_network/by_source/cell-map-node-attributes.txt.zip'
    },
    'pathguide': {
        'label':
        'Collection of metabolic and signaling pathway and molecular '\
            'interaction resources',
        'url':
        'http://pathguide.org/fullrecord.php?organisms=all&'\
            'availability=all&standards=all&order=alphabetic&DBID=%u'
    },
    'cgc': {
        'label': 'Cancer Gene Census: list of cancer related (driver) genes',
        'host': 'sftp-cancer.sanger.ac.uk',
        'file': '/files/grch38/cosmic/v76/cancer_gene_census.csv',
        'url_new': 'https://cancer.sanger.ac.uk/cosmic/file_download/'
            'GRCh38/cosmic/v94/cancer_gene_census.csv',
    },
    'havugimana': {
        'label': 'Census of human soluble protein complexes',
        'url': 'https://www.cell.com/cms/10.1016/j.cell.2012.08.011/'
            'attachment/852fdb4f-ef77-4001-b513-a6687e6e070f/mmc3.xls',
        'article': 'https://www.cell.com/fulltext/S0092-8674(12)01006-9',
    },
    'matrixdb': {
        'label': 'MatrixDB in house curated interactions, PSI-MI tab format',
        'old_url': 'http://matrixdb.ibcp.fr/download/matrixdb_CORE.tab.gz',
        'url': 'http://matrixdb.univ-lyon1.fr/download/matrixdb_CORE.tab.gz',
        'ecm_proteins': 'http://matrixdb.univ-lyon1.fr/'\
            'download//proteins_ECM.csv',
        'secreted_proteins': 'http://matrixdb.univ-lyon1.fr/'
            'download//proteins_Secreted.csv',
        'membrane_proteins': 'http://matrixdb.univ-lyon1.fr/'
            'download//proteins_Membrane.csv',
    },
    'innatedb': {
        'label': 'InnateDB PSI-MI tab',
        'url':
        'http://innatedb.com/download/interactions/innatedb_ppi.mitab.gz'
    },
    'dip': {
        'label': 'DIP PSI-MI tab',
        'login': 'http://dip.doe-mbi.ucla.edu/dip/Login.cgi',
        'url_old': 'http://dip.mbi.ucla.edu/dip/file?'
            'ds=current&fn=Hsapi20170205%s&ff=txt',
        'ik': 'http://dip.doe-mbi.ucla.edu/dip/DIPview.cgi?IK=%u',
        'url': 'https://dip.doe-mbi.ucla.edu/dip/File.cgi?FN=2017/tab25/'
            'Hsapi20170205.txt.gz',
        'url_rescued': 'https://rescued.omnipathdb.org/dip/'
            'Hsapi20170205%s.txt.gz',
    },
    'vaquerizas2009': {
        'label': 'A census of human transcription factors: function, '
        'expression and evolution; Supplementary Table S3',
        'url': (
            r'https://static-content.springer.com/esm/'
            r'art%3A10.1038%2Fnrg2538/MediaObjects/'
            r'41576_2009_BFnrg2538_MOESM6_ESM.txt'
        ),
    },
    'spike': {
        'label': 'SPIKE database XML',
        'url': 'http://www.cs.tau.ac.il/~spike/download/LatestSpikeDB.xml.zip'
    },
    'mppi': {
        'label': 'MIPS PPI full download',
        'url': 'http://mips.helmholtz-muenchen.de/proj/ppi/data/mppi.gz',
        'url_rescued': 'https://rescued.omnipathdb.org/mppi.gz',
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
        'mitab': 'https://ftp.ebi.ac.uk/pub/databases/intact/'
        'current/psimitab/intact.zip'
    },
    'death': {
        'label': 'DeathDomain webpage',
        'url_dead': 'http://www.deathdomain.org/proteins/show?family=%s',
        'url_alive': 'http://rescued.omnipathdb.org/deathdomain.tsv',
    },
    'string': {
        'label': 'STRING',
        'actions': 'https://stringdb-downloads.org/download/'
            'protein.actions.v12.0/%u.protein.actions.v12.0.txt.gz',
        'links': 'https://stringdb-downloads.org/download/protein.links.detailed.v12.0/%u'
            '.protein.links.detailed.v12.0.txt.gz',
        'physical_links': 'https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/%u'
            '.protein.physical.links.detailed.v12.0.txt.gz',
        'species': 'https://stringdb-downloads.org/download/species.v12.0.txt'
    },
    'wikipw': {
        'label': 'WikiPathways human biopax',
        'biopax_l3_old': 'http://wikipathways.org//wpi/cache/wikipathways_'
            'Homo_sapiens_Curation-AnalysisCollection__owl.zip',
        'biopax_l3': 'http://data.wikipathways.org/20191010/gpml/'
            'wikipathways-20191010-gpml-Homo_sapiens.zip',
    },
    'pwcommons': {
        'label': 'PathwayCommons binary SIF files',
        'url': 'https://www.pathwaycommons.org/archives/PC2/v%u/'
            'PathwayCommons%u.%s.hgnc.sif.gz',
    },
    'homologene': {
        'label': 'NCBI HomoloGene data, recent release',
        'url_original': 'https://ftp.ncbi.nih.gov/pub/HomoloGene/'
            'current/homologene.data',
        'url_rescued': 'https://rescued.omnipathdb.org/homologene.data',
    },
    'mirbase' : {
        'label': 'miRBase: miRNA main reference database',
        'aliases_old': 'ftp://mirbase.org/pub/mirbase/CURRENT/aliases.txt.gz',
        'aliases': 'https://mirbase.org/ftp/CURRENT/aliases.txt.gz',
    },
    'mir2dis': {
        'label': 'miR2Disease experimentally validated'\
            'miRNA-target interactions',
        'url': 'http://watson.compbio.iupui.edu:8080'\
            '/miR2Disease/download/miRtar.txt',
        'url_rescued': 'http://rescued.omnipathdb.org/miR2Disease_miRtar.txt',
    },
    'mirdeathdb': {
        'label': 'miRDeathDB experimentally verified'\
            'miRNA-target interactions',
        'url': 'http://www.rna-world.org/mirdeathdb/data'\
            '/miRDeathDB_all_data.txt',
        'url_rescued': 'http://rescued.omnipathdb.org/miRDeathDB_all_data.txt'
    },
    'ncrdeathdb': {
        'label': 'Successor of miRDeathDB. Apoptosis related non-coding RNA'\
            ' literature curated regulatory interactions',
        'url': 'http://www.rna-society.org/ncrdeathdb/data/'\
            'allNcRNACelldeathData.txt',
        'url_rescued': 'https://rescued.omnipathdb.org/'
            'ncRDeathDB_allNcRNACelldeathData.txt',
    },
    'mirecords': {
        'label': 'miRecords experimentally validated'\
            'miRNA-target interactions',
        'url': 'http://c1.accurascience.com/miRecords/download_data.php?v=4',
        'url_rescued': 'http://rescued.omnipathdb.org/miRecords_v4.xls',
    },
    'mirtarbase': {
        'label': 'miRTarBase experimentally validated'\
            'miRNA-target interactions',
        'strong_old': 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/'\
            '6.1/miRTarBase_SE_WR.xls',
        'strong': 'https://mirtarbase.cuhk.edu.cn/~miRTarBase/'\
            'miRTarBase_2019/cache/download/8.0/miRTarBase_SE_WR.xls',
        'all_old': 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/'\
            '6.1/miRTarBase_MTI.xlsx',
        'all': 'https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2019/'\
            'cache/download/8.0/miRTarBase_MTI.xlsx',
        'curated': 'https://mirtarbase.cuhk.edu.cn/~miRTarBase/'\
            'miRTarBase_2019/cache/download/8.0/MicroRNA_Target_Sites.xlsx',
    },
    'lncdisease': {
        'label': 'lncRNA and disease database: experimentally '\
            'verified lncRNA interactions',
        'url': 'http://210.73.221.6/files/images/ldd/data2.txt',
        'url_rescued': 'http://rescued.omnipathdb.org/LncRNADisease_data2.txt'
    },
    'lncrnadb': {
        'label': 'Literature curated lncRNA interactions',
        'url': 'http://lncrna.com/rest/all/nomenclature/species/association',
        'url_rescued': 'http://rescued.omnipathdb.org/'\
            'lncrnadb_association.xml',
    },
    'transmir': {
        'label': 'Literature curated TF-miRNA interactions',
        'url_old': 'http://www.cuilab.cn/files/images/transmir/'\
            'transmir_v1.2.txt',
        'url': 'http://www.cuilab.cn/files/images/transmir2/'
            'download/literature/hsa.tsv.gz',
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
    },
    'rami': {
        'label': 'Ligand-receptor interactions from Ramilowski 2015',
        'url_old': 'https://media.nature.com/original/nature-assets/ncomms/'\
            '2015/150722/ncomms8866/extref/ncomms8866-s3.xlsx',
        'fantom_url': 'http://fantom.gsc.riken.jp/5/suppl/'\
            'Ramilowski_et_al_2015/data/PairsLigRec.txt',
        'loc': 'http://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/'\
            'data/SubcelLoc.Ages.Proteins.txt',
        'url': r'https://static-content.springer.com/esm/'
            r'art%3A10.1038%2Fncomms8866/MediaObjects/'
            r'41467_2015_BFncomms8866_MOESM611_ESM.xlsx',
    },
    'kirouac2010': {
        'label': 'Intercellular (cytokine-receptor) interactions from'\
            'Kirouac et al 2010',
        'url_old': 'http://msb.embopress.org/content/msb/6/1/417/DC3/embed/'\
            'inline-supplementary-material-3.xls?download=true',
        'init_url': 'https://www.embopress.org/doi/full/10.1038/msb.2010.71',
        'url': r'https://www.embopress.org/action/downloadSupplement?doi='\
            r'10.1038%2Fmsb.2010.71&file=msb201071-sup-0003.xls'
    },
    'dorothea': {
        'label': 'A comprehensive resource of TF-target interactions',
        'url': 'http://saezlab.org/tfregulons/tfregulons_database_v01_'\
            '20180216__%s.tsv'
    },
    'dorothea_git': {
        'label': 'DoRothEA TF-taget database from its git repo',
        'url': 'https://github.com/saezlab/DoRothEA/raw/master/data/'
            'TFregulons/consensus/table/database_normal_20180915.csv.zip',
        'rda': 'https://github.com/saezlab/dorothea/raw/master/data/%s.rda',
    },
    'hpmri': {
        'label': 'Human Plasma Membrane Receptome "browse" view',
        'browse': 'http://www.receptome.org/Families/FamNav/'\
            'famnav.asp?undefined',
        'browse_rescued': 'https://rescued.omnipathdb.org/hpmr/famnav.html',
        'genes_rescued': 'https://rescued.omnipathdb.org/hpmr/genes.tar.gz',
    },
    'cellphonedb': {
        'label': 'Intercellular communication database from Teichmann lab',
        'interactions': 'https://www.cellphonedb.org/downloads/'\
            'interactions_cellphonedb.csv',
        'heterodimers': 'https://www.cellphonedb.org/downloads/'\
            'heterodimers.csv',
        'proteins': 'https://www.cellphonedb.org/downloads/'\
            'protein_complex_cellphonedb.csv',
    },
    'cellphonedb_git': {
        'label': 'Intercellular communication database from Teichmann lab',
        'proteins': 'https://raw.githubusercontent.com/ventolab/'\
            'cellphonedb-data/master/data/protein_input.csv',
        'complexes': 'https://raw.githubusercontent.com/ventolab/'\
            'cellphonedb-data/master/data/complex_input.csv',
        'interactions': 'https://raw.githubusercontent.com/ventolab/'\
            'cellphonedb-data/master/data/interaction_input.csv',
        'curated': 'https://raw.githubusercontent.com/ventolab/'\
            'cellphonedb-data/master/data/sources/interaction_curated.csv',
        'negative': 'https://raw.githubusercontent.com/ventolab/'\
            'cellphonedb-data/master/data/sources/excluded_interaction.csv',
    },
    'stitch': {
        'label': 'The STITCH small molecule-protein interaction database',
        'actions': 'http://stitch.embl.de/download/actions.v5.0/%u'\
            '.actions.v5.0.tsv.gz',
        'links': 'http://stitch.embl.de/download/protein_chemical.links.detailed.v5.0/%u'\
            '.protein_chemical.links.detailed.v5.0.tsv.gz'
    },
    'cspa': {
        'label': 'Cell Surface Protein Atlas',
        'url_s2': 'http://wlab.ethz.ch/cspa/data/S2_File.xlsx',
        'url_s1': 'http://wlab.ethz.ch/cspa/data/S1_File.xlsx',
    },
    'surfaceome': {
        'label': 'In silico human surfaceome',
        'url': 'http://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx',
    },
    'matrisome': {
        'label': 'MatrisomeDB 2.0: ECM proteins',
        'url_dl': 'http://matrisomeproject.mit.edu/matrisomeDB/search/?'\
            'keywords=&genesymbol=&category=1&category=2&category=5&'\
            'category=4&category=3&category=6&human=on&mouse=on&tissue=3&'\
            'tissue=4&tissue=6&tissue=5&tissue=11&tissue=12&tissue=13&'\
            'tissue=14&tissue=7&tissue=8&tissue=9&tissue=1&tissue=10&'\
            'tissue=2&download=Download',
        'url_xls': 'http://matrisomeproject.mit.edu/static/media/uploads/'\
            'Files/%s%%20in%%20silico%%20matrisome/'\
            'matrisome_%s_masterlist.xls',
        'url_rescued': 'https://rescued.omnipathdb.org/'
            'matrisome_%s_masterlist.xls',
    },
    'membranome': {
        'label': 'Database of single helix transmembrane proteins: '\
            'Membranome',
        'baseurl': 'http://lomize-group-membranome.herokuapp.com/%s%s',
    },
    'exocarta': {
        'label': 'Exosome protein, RNA and lipid database',
        'url_protein': 'http://exocarta.org/Archive/'\
            'EXOCARTA_PROTEIN_MRNA_DETAILS_5.txt',
        'url_study': 'http://exocarta.org/Archive/'
            'EXOCARTA_EXPERIMENT_DETAILS_5.txt',
    },
    'vesiclepedia': {
        'label': 'Extracellular vesicle database',
        'url_protein': 'http://microvesicles.org/Archive/'\
            'VESICLEPEDIA_PROTEIN_MRNA_DETAILS_4.1.txt',
        'url_study': 'http://microvesicles.org/Archive/'\
            'VESICLEPEDIA_EXPERIMENT_DETAILS_4.1.txt',
    },
    'locate': {
        'label': 'Locate protein localization database',
        'url': 'http://locate.imb.uq.edu.au/info_files/'\
            'LOCATE_%s_v6_20081121.xml.zip',
        'url_rescued': 'https://rescued.omnipathdb.org/'\
            'LOCATE_%s_v6_20081121.xml.zip',
    },
    'humap': {
        'label': 'Map of human protein complexes from '\
            'integration of 9,000 MS experiments',
        'humap_url': 'http://proteincomplexes.org/static/'\
            'downloads/clusters.txt',
        'humap2_url': 'http://humap2.proteincomplexes.org/static/'\
            'downloads/humap2/humap2_complexes_20200809.txt',
    },
    'adhesome': {
        'label': 'Literature curated database about adhesion related '\
            'proteins and their interactions',
        'components': 'http://www.adhesome.org/components.csv',
        'interactions': 'http://www.adhesome.org/interactions.csv',
    },
    'integrins': {
        'label': 'List of integrins from Takada et al 2007',
        'url': 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1929136/'\
            'table/T1/?report=objectonly',
    },
    'opm': {
        'label': 'Orientations of proteins in membranes database',
        'proteins': 'https://lomize-group-opm.herokuapp.com/'\
            'primary_structures?fileFormat=csv',
        'families': 'https://lomize-group-opm.herokuapp.com/'\
            'families?fileFormat=csv',
        'classes': 'https://lomize-group-opm.herokuapp.com/'\
            'classes?fileFormat=csv',
        'types': 'https://lomize-group-opm.herokuapp.com/'\
            'types?fileFormat=csv',
    },
    'topdb': {
        'label': 'Membrane protein topology database',
        'url': 'http://topdb.enzim.hu/?m=download&file=topdb_all.xml',
    },
    'hgnc': {
        'label': 'Web download interface of HGNC (genenames.org).',
        'groups': (
            'https://www.genenames.org/cgi-bin/download/custom?'
            'col=gd_app_sym&col=gd_pub_acc_ids&col=md_prot_id&'
            'col=family.name&status=Approved&hgnc_dbtag=on&'
            'order_by=gd_app_sym_sort&format=text&submit=submit'
        ),
    },
    'hippie': {
        'label': (
            'HIPPIE: human integrated protein-protein interaction resource'
        ),
        'url': (
            'http://cbdm-01.zdv.uni-mainz.de/~mschaefer/'
            'hippie/hippie_current.txt'
        ),
    },
    'cpad': {
        'label': 'Literature curated regulator-pathway-cancer relationships',
        'url': 'http://bio-bigdata.hrbmu.edu.cn/CPAD/download/CPAD-data.txt',
    },
    'compartments': {
        'label': 'Compartments subcellular localization database',
        'url': 'http://download.jensenlab.org/%s_compartment_%s_full.tsv',
    },
    'intogen': {
        'label': 'IntOGen is a database of cancer related '\
            'genes and mutations',
        'db2014_2': 'https://www.intogen.org/download?'\
            'file=intogen_cancer_drivers-2014.12.zip',
        'db2020': 'https://www.intogen.org/download?file='\
            'IntOGen-Drivers-20200201.zip',
    },
    'cancersea': {
        'label': 'Cancer Single Cell State Atlas',
        'url': 'http://biocc.hrbmu.edu.cn/CancerSEA/goDownload',
        'data_url': 'http://biocc.hrbmu.edu.cn/CancerSEA/%s',
        'rescued': 'https://rescued.omnipathdb.org/CancerSEA.html',
        'rescued_data': 'https://rescued.omnipathdb.org/cancersea/%s',
    },
    'signalink': {
        'label': 'Literature curated activity flow signaling network',
        'nodes': 'http://rescued.omnipathdb.org/SLK3_nodes_2020_03_23.csv',
        'edges': 'http://rescued.omnipathdb.org/SLK3_layer0_2020_03_23.csv',
    },
    'zhong2015': {
        'url': 'http://rescued.omnipathdb.org/zhong2015_S1.tsv',
    },
    'phosphopoint': {
        'url': 'http://rescued.omnipathdb.org/phosphopoint.csv',
    },
    'lmpid': {
        'url': 'http://rescued.omnipathdb.org/LMPID_DATA_pubmed_ref.xml',
    },
    'nrf2ome': {
        'url': 'http://rescued.omnipathdb.org/nrf2ome.csv',
    },
    'arn': {
        'url': 'http://rescued.omnipathdb.org/arn_curated.csv',
    },
    'alzpathway': {
        'url': 'http://rescued.omnipathdb.org/alzpw-ppi.csv',
    },
    'protmapper': {
        'label': 'Integration of kinase-substrate data across many '\
            'resources and text mining',
        'url': 'https://www.biorxiv.org/content/biorxiv/early/2019/11/06/'\
            '822668.1/DC1/embed/media-1.zip?download=true',
        'files': [
            'export.csv',
            'evidences.csv',
        ],
    },
    'lrdb': {
        'label': 'Curated ligand-receptor interactions',
        'url': 'https://raw.githubusercontent.com/SCA-IRCM/LRdb/'\
            'master/LRdb_122019.txt',
    },
    'baccin2019': {
        'label': 'Ligand receptor pairs in Supplementary Table 3 from'\
            'Baccin 2019',
        'url': r'https://static-content.springer.com/esm/'\
            r'art%3A10.1038%2Fs41556-019-0439-6/MediaObjects/'\
            r'41556_2019_439_MOESM3_ESM.xlsx',
    },
    'kea': {
        'label': (
            'Kinase enrichment analysis, contains many kinase related '
            'datasets compiled from various databases and the literature'
        ),
        'kinase_substrate': 'https://www.maayanlab.net/KEA2/gsl/'
            'kinase_substrate_phospho-site_level_pmid_resource_database.txt',
    },
    'iptmnet': {
        'label': 'A database of PTMs and enzyme-substrate interactions '
            'integrated from other resources and text mining',
        'ptms': 'https://research.bioinformatics.udel.edu/iptmnet_data/'
            'files/current/ptm.txt',
        'scores': 'https://research.bioinformatics.udel.edu/iptmnet_data/'
            'files/current/score.txt',
        'proteins': 'https://research.bioinformatics.udel.edu/iptmnet_data/'
            'files/current/protein.txt',
    },
    'icellnet': {
        'label': 'A literature curated collection of '
            'ligand-receptor interactions.',
        'url_old': 'https://www.biorxiv.org/content/biorxiv/early/2020/03/05/'
            '2020.03.05.976878/DC2/embed/media-2.xlsx?download=true',
        'url': 'https://raw.githubusercontent.com/soumelis-lab/ICELLNET/'
            'master/data/ICELLNETdb.csv',
    },
    'phobius': {
        'label': (
            'Plasma membrane protein orientation prediction for the '
            'human proteome'
        ),
        'url': 'http://rescued.omnipathdb.org/phobius_human_proteome.tsv',
    },
    'pro': {
        'label': 'The Protein Ontology',
        'url': 'https://proconsortium.org/download/current/pro_reasoned.obo',
        'mapping': 'https://proconsortium.org/download/current/'
            'promapping.txt',
    },
    'genecards': {
        'label': (
            'GeneCards provides a comprehensive summary about the function, '
            'lovalization, expression, variants and other properties of '
            'genes and their products',
        ),
        'url': 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s',
    },
    'almen2009': {
        'label': 'Classification of membrane proteins 10.1186/1741-7007-7-50',
        'url': (
            r'https://static-content.springer.com/esm/'
            r'art%3A10.1186%2F1741-7007-7-50/MediaObjects/'
            r'12915_2009_258_MOESM1_ESM.xls'
        ),
    },
    'ipi': {
        'label': (
            'International Protein Index, a protein reference database '
            'which has been shut down by 2012'
        ),
        'url': (
            'ftp://ftp.ebi.ac.uk/pub/databases/IPI/last_release/'
            'current/ipi.HUMAN.xrefs.gz'
        ),
    },
    'cellcellinteractions': {
        'label': 'Cell-cell Interaction Database from Bader lab',
        'url': (
            'https://baderlab.org/CellCellInteractions?'
            'action=AttachFile&do=get&target=protein_types.txt'
        ),
    },
    'italk': {
        'label': 'Intercellular communication database from iTalk',
        'url': (
            'https://github.com/Coolgenome/iTALK/blob/master/'
            'data/LR_database.rda?raw=true'
        ),
    },
    'embrace': {
        'label': 'Ligand-receptor interactions from embryonic brain '
            'cell extraction by FACS (EMBRACE)',
        'article': 'https://www.cell.com/iscience/fulltext/'
            'S2589-0042(19)30405-5',
        'url': 'https://www.cell.com/cms/10.1016/j.isci.2019.10.026/'
            'attachment/a1f66cee-15cb-4690-b95d-2480fa576b86/mmc2.xlsx',
    },
    'tcdb': {
        'label': 'Transporter classification database',
        'url_tcdb': 'http://www.tcdb.org/cgi-bin/projectv/public/tcdb.dr',
        'url_families': 'http://www.tcdb.org/cgi-bin/'
            'projectv/public/families.py',
        'url_acc2tc': 'http://www.tcdb.org/cgi-bin/'
            'projectv/public/acc2tcid.py',
    },
    'mcam': {
        'label': 'Mammalian Cell Adhesion Molecule Database',
        'url': 'https://app1.unmc.edu/mcam/exel.cfm?seq=HumanUniprot',
    },
    'gpcrdb': {
        'label': 'A comprehensive database about GPCRs',
        'families': 'https://raw.githubusercontent.com/protwis/gpcrdb_data/'
            'master/protein_data/proteins_and_families.txt',
    },
    'immunoglobe': {
        'label': 'A literature curated resource about immune interactions',
        'url': 'http://network.immunoglobe.org/networks.json',
    },
    'ensembl': {
        'label': 'Genome and transcript reference database',
        'biomart_url': 'http://www.ensembl.org/biomart/martservice?query=%s',
        'species': 'https://www.ensembl.org/info/about/species.html',
        'arraytypes': 'https://rest.ensembl.org/regulatory/species/'
            '%s/microarray?content-type=application/json',
    },
    'wojtowicz2020': {
        'label': 'Ligand-receptor interactions from a screen of human'
            'IgSF cell surface interactome.',
        'article': 'https://www.cell.com/cell/fulltext/S0092-8674(20)30933-8',
        'url': 'https://www.cell.com/cms/10.1016/j.cell.2020.07.025/'
            'attachment/656c4057-8b3b-4aa2-99ce-73278e21d93d/mmc4.xlsx',
    },
    'celltalkdb': {
        'label': 'Manually curated ligand-receptor interactions',
        'url': 'http://tcm.zju.edu.cn/celltalkdb/handler/download.php',
        'ref_url': 'http://tcm.zju.edu.cn/celltalkdb/',
        'init_url': 'http://tcm.zju.edu.cn/celltalkdb/download.php',
    },
    'cellchatdb': {
        'label': 'Manually curated ligand-receptor interactions',
        'url': 'https://raw.githubusercontent.com/jinworks/CellChat/'
            'master/data/%s.rda',
    },
    'connectomedb2020': {
        'label': 'Ligand-receptor interactions from original curation '
            'and from other databases',
        'url': 'https://raw.githubusercontent.com/asrhou/NATMI/'
            'gh-pages/_data/lsps.csv',
    },
    'talklr': {
        'label': 'Ligand-receptor interactions from other databases with '
            '187 literature references',
        'url': 'https://github.com/yuliangwang/talklr/raw/master/'
            'data/receptor_ligand.rda',
    },
    'humancellmap': {
        'label': 'A BioID screen for proximity based localization of '
            'proteins in HEK cells',
        'preys_url': 'https://humancellmap.org/resources/downloads/'
            'preys-latest.txt',
    },
    'cellcall': {
        'label': 'Database of ligand-receptor-TF pathways',
        'url': 'https://raw.githubusercontent.com/ShellyCoder/cellcall/'
            'master/inst/extdata/new_ligand_receptor_TFs%s.txt',
    },
    'cancerdrugs_db': {
        'label': 'Cancer Drug Database of licensed drugs by the Anticancer Fund',
        'url': 'https://acfdata.coworks.be/cancerdrugsdb.txt',
    },
    'biomodels': {
        'label': 'RESTful API base URL of EMBL-EBI BioModels (https://www.ebi.ac.uk/biomodels/)',
        'url': 'https://www.ebi.ac.uk/biomodels/',
    },
    'biogps': {
        'label': 'A collection of human and mouse tissue and cell line '
            'gene expression profiles',
        'url': 'http://plugins.biogps.org/download/%s',
        'datasets': {
            'human_gene_atlas_ave': 'gnf1h-gcrma.zip',
            'human_gene_atlas': 'gnf1h-gcrma-unaveraged.zip',
            'human_nci60': 'NCI60_U133A_20070815.raw.csv.zip',
            'mouse_gene_atlas_1': 'GNF1M_geneatlas_20120817.zip',
            'mouse_gene_atlas_2_ave':
                'geneatlas_MOE430_20090327.raw.avg.csv.zip',
            'mouse_gene_atlas_2': 'geneatlas_MOE430_20090327.raw.csv.zip',
            'rat_gene_atlas': 'Ratlas-gcRMA.zip',
            'human_primary_tumors':
                'primary_tumor_classification_data.txt.zip',
            'mouse_adipose_eqtl': 'adipose_MOE430_20071005.raw.csv.gz',
        },
    },
    'unichem': {
        'label': 'Mapping between drug compound IDs',
        'sources_old':
            'https://www.ebi.ac.uk/unichem/legacy/ucquery/listSources',
        'sources': 'https://www.ebi.ac.uk/unichem/api/v1/sources/',
        'mapping': 'https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/'
            'wholeSourceMapping/src_id%s/src%ssrc%s.txt.gz',
    },
    'redhuman': {
        'label': 'A modern and reduced complexity human metabolic network',
        'url': 'https://github.com/saezlab/ocean/raw/master/data/'
            'recon2_redhuman.RData',
    },
    'cellinker': {
        'label': 'Human and mouse ligand-receptor interactions',
        'lr': 'http://www.rna-society.org/cellinker/download/%s.txt',
        'complex': 'http://www.rna-society.org/cellinker/'
            'download/%s-complex.csv',
        'smol': 'http://www.rna-society.org/cellinker/download/%s-sMOL.txt',
    },
    'cellinker_rescued': {
        'label': 'Human and mouse ligand-receptor interactions',
        'lr': 'https://rescued.omnipathdb.org/cellinker/%s.txt',
        'complex': 'https://rescued.omnipathdb.org/cellinker/%s-complex.csv',
        'smol': 'https://rescued.omnipathdb.org/cellinker/%s-sMOL.txt',
    },
    'pubchem': {
        'label': 'NCBI database of bioactive compounds, substances and '
            'assays',
        'ftp': 'https://ftp.ncbi.nlm.nih.gov/pubchem/%s/Extras/%s-%s.gz',
    },
    'scconnect': {
        'label': 'Database of ligand-receptor interactions',
        'annot': 'https://raw.githubusercontent.com/JonETJakobsson/'
            'scConnect/master/scConnect/data/Gene_annotation/2020-5/'
            '%s/%ss.csv',
        'intera': 'https://github.com/JonETJakobsson/scConnect/raw/master/'
            'scConnect/data/GTP_tables_clean/2020-5/interactions.csv',
    },
    'progeny': {
        'label': 'Pathway responsive genes with weights',
        'url': 'https://github.com/saezlab/progeny/raw/master/'
            'data/model_%s_full.rda',
    },
    'celltypist': {
        'label': 'Immune cell type marker genes',
        'url': 'https://raw.githubusercontent.com/Teichlab/celltypist_wiki/'
            'main/atlases/Pan_Immune_CellTypist/v2/encyclopedia/'
            'encyclopedia_table.xlsx',
    },
    'panglaodb': {
        'label': 'Cell type marker genes',
        'url': 'https://panglaodb.se/markers/'
            'PanglaoDB_markers_27_Mar_2020.tsv.gz',
    },
    'cytosig': {
        'label': 'Cytokine perturbation footprints',
        'url': 'https://raw.githubusercontent.com/data2intelligence/'
            'CytoSig/master/CytoSig/signature.centroid',
    },
    'proteins': {
        'label': 'Key biological data from UniProt and data from Large Scale '
            'Studies (LSS) mapped to UniProt.',
        'url': 'https://www.ebi.ac.uk/proteins/api/%s',
    },
    'lambert2018': {
        'label': 'A new compendium of human transcription factors',
        's1': 'https://www.cell.com/cms/10.1016/j.cell.2018.01.029/'
            'attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx',
        'article': 'https://www.cell.com/cell/fulltext/S0092-8674(18)30106-5',
    },
    'interpro': {
        'label': 'Protein families, domains and functional sites',
        'entries': 'https://ftp.ebi.ac.uk/pub/databases/'
            'interpro/current_release/interpro.xml.gz',
        'annotations': 'https://www.ebi.ac.uk/interpro/api/entry/InterPro/'
            'protein/%s/taxonomy/uniprot/%s?page_size=%u',
        'interpro2go': 'http://www.geneontology.org/external2go/interpro2go',
    },
    'drugcentral': {
        'label': 'Drug-target interactions',
        'interactions': 'https://unmtid-shinyapps.net/download/DrugCentral'
            '/2021_09_01/drug.target.interaction.tsv.gz',
        'SMILES_InChI' : 'https://unmtid-shinyapps.net/download/DrugCentral'
            '/2021_09_01/structures.smiles.tsv',
    },
    'drugbank': {
        'label': 'DrugBank database',
        'all_structures': 'https://go.drugbank.com/releases/5-1-9/'
            'downloads/all-structure-links',
        'all_drugs': 'https://go.drugbank.com/releases/5-1-9/downloads/'
            'all-drug-links',
        'drug_target_identifiers' : 'https://go.drugbank.com/releases/'
            '5-1-9/downloads/target-all-polypeptide-ids',
        'drug_enzyme_identifiers' : 'https://go.drugbank.com/releases/5-1-9/'
            'downloads/enzyme-all-polypeptide-ids',
        'drug_carrier_identifiers' : 'https://go.drugbank.com/releases/5-1-9/'
            'downloads/carrier-all-polypeptide-ids',
        'drug_transporter_identifiers' : 'https://go.drugbank.com/releases/'
            '5-1-9/downloads/transporter-all-polypeptide-ids',
        'full_database': 'https://go.drugbank.com/releases/5-1-9'
            '/downloads/all-full-database',
    },
    'chembl': {
        'label': 'ChEMBL database',
        'url': 'https://www.ebi.ac.uk',
        'target': '/chembl/api/data/target.json?limit=1000',
        'assay' : '/chembl/api/data/assay.json?limit=1000',
        'activity' : '/chembl/api/data/activity.json?limit=1000',
        'molecule' : '/chembl/api/data/molecule.json?limit=1000',
        'document' : '/chembl/api/data/document.json?limit=1000',
        'drug_indication' : '/chembl/api/data/drug_indication.json?limit=1000',
        'mechanism' : '/chembl/api/data/mechanism.json?limit=1000',
    },
    'hpo': {
        'label': 'HPO database',
        'ontology': 'https://purl.obolibrary.org/obo/hp.obo',
        'ontology_new': 'https://purl.obolibrary.org/obo/hp.obo'
            'human-phenotype-ontology/master/hp.obo',
        'disease' : 'https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa',
        'gene' : 'https://purl.obolibrary.org/obo/hp/hpoa/'
            'genes_to_phenotype.txt',
    },
    'oma': {
        'label': 'A database for the inference of orthologs among '
            'complete genomes',
        'url': 'https://omabrowser.org/api/pairs/',
    },
    'pathophenodb': {
        'label': 'A database of pathogens and their phenotypes '
            'for diagnostic support in infections.',
        'url': 'http://patho.phenomebrowser.net/sparql/',
    },
    'clinvar': {
        'label': 'A database archiving submissions that interpret '
            'the effect of a single variant or set of variants on phenotype.',
        'url': 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/'
            'tab_delimited/variant_summary.txt.gz',
        'url_citations': 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt',
    },
    'pharos_api': {
        'label': 'Pharos GraphQL API Endpoint',
        'url': 'https://pharos-api.ncats.io/graphql',
    },
    'ctdbase': {
        'label': 'Comparative Toxicogenomics Database',
        'url': 'http://ctdbase.org/reports/%s',
    },
    'diseases': {
        'label': 'Disease-gene associations mined from literature',
        'url': 'https://download.jensenlab.org/human_disease_%s_%s.tsv',
    },
    'sider': {
        'label': (
            'A Database contains information on marketed medicines and their '
            'recorded adverse drug reactions.'
        ),
        'drug_names':
            'http://sideeffects.embl.de/media/download/drug_names.tsv',
        'drug_atcs':
            'http://sideeffects.embl.de/media/download/drug_atc.tsv',
        'meddra_all':
            'http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz',
        'meddra_freq':
            'http://sideeffects.embl.de/media/download/meddra_freq.tsv.gz',
        'meddra_tsv':
            'http://sideeffects.embl.de/media/download/meddra.tsv.gz',
    },
    'humsavar': {
        'label': (
            'A Database contains all missense variants annotated in '
            'UniProtKB/Swiss-Prot human entries'
        ),
        'url': (
            'https://ftp.uniprot.org/pub/databases/uniprot/current_release/'
            'knowledgebase/complete/docs/humsavar.txt'
        ),
    },
    'twosides': {
        'label': 'A comprehensive database drug-drug-effect relationships',
        'url': 'https://tatonettilab.org/resources/nsides/TWOSIDES.csv.gz'
    },
    'offsides': {
        'label': (
            'A database of drug side-effects that were found, '
            'but are not listed on the official FDA label'
        ),
        'url': 'https://tatonettilab.org/resources/nsides/OFFSIDES.csv.gz'
    },
    'adrecs': {
        'label': (
            'A comprehensive ADR ontology database that provides both '
            'standardization and hierarchical classification of ADR terms.'
        ),
        'drug_information':
            'http://bioinf.xmu.edu.cn/ADReCS/download/Drug_information.xlsx',
        'terminology':
            'http://www.bio-add.org/ADReCS/download/ADR_ontology.xlsx',
        'adrecs_drugs':
            'http://www.bio-add.org/ADReCS/download/Drug_ADR.txt.gz',
    },
    'ddinter': {
        'label': (
            'DDInter is a comprehensive, professional, and open-access '
            'database specific to drug-drug interactions.'
        ),
        'source': 'http://ddinter.scbdd.com/ddinter/data-source/',
        'mapping': 'http://ddinter.scbdd.com/ddinter/drug-detail/%s',
        'interaction':
            'http://ddinter.scbdd.com/ddinter/grapher-datasource/%s',
    },
    'reactome': {
        'sbml': 'https://reactome.org/download/current/homo_sapiens.3.1.sbml.tgz',
        'biopax_l3': 'https://reactome.org/download/current/biopax.zip',
        'biopax_l2': 'https://reactome.org/download/current/biopax2.zip',
        'txt': 'https://reactome.org/download/current/%s.txt',
        'label': (
            'REACTOME is an open-source, open access, manually curated '
            'and peer-reviewed pathway database.'
        ),
        'uniprot2all': 'https://reactome.org/download/current/'
            'UniProt2Reactome_All_Levels.txt',
        'chebi2all': 'https://reactome.org/download/current/'
            'ChEBI2Reactome_All_Levels.txt',
        'pathways': 'https://reactome.org/download/current/'
            'ReactomePathways.txt',
        'pathway_relations': 'https://reactome.org/download/current/'
            'ReactomePathwaysRelation.txt',
    },
    'rhea': {
        'label': 'database of metabolic reactions',
        'url': 'https://ftp.expasy.org/databases/rhea/tsv/%s.tsv',
    },
    'compath': {
        'label': (
            'A database of proposed and accepted mappings by the '
            'users/curators between a pair of pathways.'
        ),
        'url': 'https://compath.scai.fraunhofer.de/export_mappings',
    },
    'trrust': {
        'label': 'Transcriptional Regulatory Relationships Unraveled by Sentence-based Text mining',
        'scraping_url': 'https://www.grnpedia.org/trrust/data/search_list.%s.htm',
        'tsv_url': 'https://www.grnpedia.org/trrust/data/trrust_rawdata.%s.tsv'
    },
    'collectri': {
        'label': 'Comprehensive collection of transcriptional regulatory interactions',
        'url': 'https://rescued.omnipathdb.org/CollecTRI.csv',
    },
    'ramp': {
        'label': 'RaMP metabolomic pathway and metabolite identifier database',
        'url': 'https://figshare.com/ndownloader/files/38534654',
        'api': 'https://rampdb.nih.gov/api/%s',
        'github_sqlite': (
            'https://media.githubusercontent.com/media/ncats/RaMP-DB/'
            'refs/heads/main/db/RaMP_SQLite_v%s.sqlite.gz'
        ),
    },
    'expasy': {
        'label': 'Enzyme nomenclature database',
        'enzclass': 'https://ftp.expasy.org/databases/enzyme/enzclass.txt',
        'enzymes': 'https://ftp.expasy.org/databases/enzyme/enzyme.dat',
    },
    'hmdb': {
        'label': 'Human Metabolome Database',
        'metabolites': 'https://hmdb.ca/system/downloads/'
            'current/hmdb_metabolites.zip',
        'proteins': 'https://hmdb.ca/system/downloads/'
            'current/hmdb_proteins.zip',
        'svg': 'https://hmdb.ca/structures/%s/image.svg',
        'sdf': 'https://hmdb.ca/system/downloads/current/structures.zip',
    },
    'lipidmaps': {
        'label': 'A database of lipid structures',
        'sdf': 'https://lipidmaps.org/files/?file=LMSD&ext=sdf.zip',
    },
    'opentargets': {
        'label': 'Open Targets',
        'url': 'http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/'
            '22.11/output/etl/json/%s',
    },
    'gutmgene': {
        'label': 'A comprehensive database for target genes of gut microbes and microbial metabolites',
        'url_human': 'http://bio-annotation.cn/gutmgene/public/res/Gut%20Microbe%20and%20Gene-human.txt',
        'url_mouse': 'http://bio-annotation.cn/gutmgene/public/res/Gut%20Microbe%20and%20Gene-mouse.txt',
    },
    'swisslipids': {
        'label': 'A knowledge resource for lipids and their biology',
        'url': 'https://swisslipids.org/api/file.php?cas=download_files&file=%s.tsv',
    },
    'lipidmaps': {
        'label': '',
        'lmsd': 'https://lipidmaps.org/files/?file=LMSD&ext=sdf.zip',
    },
}


files = {
    'phosphosite': {
        'curated': os.path.join(
            cache.get_cachedir(),
            'phosphosite_curated.pickle',
        ),
        'noref': os.path.join(
            cache.get_cachedir(),
            'phosphosite_noref.pickle',
        )
    },
}
