"""BioGRID: archive reader, mapper and extractor on lines of the Tab 3.0 file of all organisms."""

import io
import json
import re
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import qualifiers, slot_name
from omnipath_core.keys import entity_key, relation_key
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import apo, biogrid, bto, go, hpo
from pypath.inputs_v2.parsers.biogrid import cell_line_accessions, iter_interactions
from pypath.inputs_v2.parsers.obo import parse_obo_text

HEADER = [
    '#BioGRID Interaction ID', 'Entrez Gene Interactor A',
    'Entrez Gene Interactor B', 'BioGRID ID Interactor A',
    'BioGRID ID Interactor B', 'Systematic Name Interactor A',
    'Systematic Name Interactor B', 'Official Symbol Interactor A',
    'Official Symbol Interactor B', 'Synonyms Interactor A',
    'Synonyms Interactor B', 'Experimental System',
    'Experimental System Type', 'Author', 'Publication Source',
    'Organism ID Interactor A', 'Organism ID Interactor B', 'Throughput',
    'Score', 'Modification', 'Qualifications', 'Tags', 'Source Database',
    'SWISS-PROT Accessions Interactor A', 'TREMBL Accessions Interactor A',
    'REFSEQ Accessions Interactor A', 'SWISS-PROT Accessions Interactor B',
    'TREMBL Accessions Interactor B', 'REFSEQ Accessions Interactor B',
    'Ontology Term IDs', 'Ontology Term Names', 'Ontology Term Categories',
    'Ontology Term Qualifier IDs', 'Ontology Term Qualifier Names',
    'Ontology Term Types', 'Organism Name Interactor A',
    'Organism Name Interactor B',
]
HUMAN = 'Homo sapiens'
NO_TERMS = ['-'] * 6

# Lines of BIOGRID-ALL-5.0.262.tab3.txt (the line numbers are those of the file).
# Line 2: a two-hybrid pair from a small-scale study.
TWO_HYBRID = [
    '103', '6416', '2318', '112315', '108607', '-', '-', 'MAP2K4', 'FLNC',
    'JNKK|JNKK1|MAPKK4|MEK4|MKK4|PRKMK4|SAPKK-1|SAPKK1|SEK1|SERK1|SKK1',
    'ABP-280|ABP280A|ABPA|ABPL|FLN2|MFM5|MPD4', 'Two-hybrid', 'physical',
    'Marti A (1997)', 'PUBMED:9006895', '9606', '9606', 'Low Throughput',
    '-', '-', '-', '-', 'BIOGRID', 'P45985', '-', 'NP_003001|NP_001268364',
    'Q14315', 'Q59H94', 'NP_001120959|NP_001449', *NO_TERMS, HUMAN, HUMAN,
]
# Line 888: the first gene has two Swiss-Prot entries and a systematic name.
TWO_PROTEINS = [
    '30184', '140290', '10608', '126602', '115854', 'PRED77', 'RP11-478C1.1',
    'TCP10L', 'MXD4', 'C21orf77|TCP10A-2', 'MAD4|MST149|MSTP149|bHLHc12',
    'Two-hybrid', 'physical', 'Jiang DJ (2004)', 'PUBMED:15469726', '9606',
    '9606', 'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'Q9NV44|Q8TDR4',
    '-', 'NP_653260', 'Q14582', '-', 'NP_006445', *NO_TERMS, HUMAN, HUMAN,
]
# Line 121404: a protein captured with an RNA; RNU2-1 has no protein.
PROTEIN_RNA = [
    '241375', '10915', '6066', '116120', '111993', '-', '-', 'TCERG1',
    'RNU2-1', 'CA150|TAF2S|Urn1', 'RNU2|U2', 'Affinity Capture-RNA',
    'physical', 'Lin KT (2004)', 'PUBMED:15456888', '9606', '9606',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'O14776', '-',
    'NP_006697|NP_001369477|NP_001035095', '-', '-', '-', *NO_TERMS, HUMAN,
    HUMAN,
]
# Line 121443: a kinase assay with the modification it observed.
KINASE_ASSAY = [
    '241416', '23235', '200186', '116840', '128308', '-', 'RP11-422P24.6',
    'SIK2', 'CRTC2', 'LOH11CR1I|QIK|SNF1LK2', 'TORC-2|TORC2',
    'Biochemical Activity', 'physical', 'Screaton RA (2004)',
    'PUBMED:15454081', '9606', '9606', 'Low Throughput', '-',
    'Phosphorylation', '-', '-', 'BIOGRID', 'Q9H0K1', 'A0A024R3G7',
    'NP_056006', 'Q53ET0', '-', 'NP_859066', *NO_TERMS, HUMAN, HUMAN,
]
# Line 129890: two proteins in the same fraction of a preparation.
CO_FRACTIONATION = [
    '251118', '10849', '9533', '116060', '114909', '-', 'RP3-337H4.4',
    'CD3EAP', 'POLR1C', 'ASE-1|ASE1|CAST|PAF49', 'RPA39|RPA40|RPA5|RPAC1|TCS3',
    'Co-fractionation', 'physical', 'Yamamoto K (2004)', 'PUBMED:15226435',
    '9606', '9606', 'Low Throughput', '-', '-', '-', '-', 'BIOGRID',
    'O15446', '-', 'NP_036231|NP_001284519', 'O15160', '-',
    'NP_001305805|NP_976035', *NO_TERMS, HUMAN, HUMAN,
]
# Line 282914: reported by a study of both throughputs, with a score.
BOTH_THROUGHPUTS = [
    '423426', '23239', '7398', '116843', '113241', '-', '-', 'PHLPP1', 'USP1',
    'PHLPP|PLEKHE1|SCOP', 'UBP', 'Affinity Capture-MS', 'physical',
    'Sowa ME (2009)', 'PUBMED:19615732', '9606', '9606',
    'High Throughput|Low Throughput', '5.1', '-',
    'exogenous expression of bait|reciprocal IP-MS experiment', '-',
    'BIOGRID', 'O60346', '-', 'NP_919431', 'O94782', '-',
    'NP_001017415|NP_001017416|NP_003359', 'BTO:0000007', 'HEK-293 cell',
    'cell line', 'BTO:0000414', 'epithelial cell', '-', HUMAN, HUMAN,
]
# Line 2305902: a preprint, cited by its DOI.
PREPRINT = [
    '3308836', '259307', '2734', '129252', '108996', 'UNQ636/PRO1265',
    'RP11-252A24.7', 'IL4I1', 'GLG1', 'FIG1|LAAO|LAO',
    'CFR-1|ESL-1|MG-160|MG160', 'Affinity Capture-MS', 'physical',
    'Gatineau J (2021)', 'DOI:10.1101/2021.12.23.473686', '9606', '9606',
    'High Throughput', '-', '-', '-', '-', 'BIOGRID', 'Q96RQ9', '-',
    'NP_758962|NP_001244947|NP_001244946|NP_690863|NP_001372568', 'Q92896',
    '-', 'NP_036333|NP_001139139|NP_001139138', *NO_TERMS, HUMAN, HUMAN,
]
# Line 2333622: a viral participant without Entrez Gene or protein accession.
VIRAL = [
    '3337444', '-', '3304', '4383869', '109536', '-', 'DAAP-21F2.7', 'ORF3b',
    'HSPA1B', 'SARS-CoV2 ORF3b|SARS-CoV-2 ORF3b|3b', 'HSP70-1B|HSP70-2',
    'Protein-RNA', 'physical', 'Banerjee AK (2020)', 'PUBMED:33080218',
    '2697049', '9606', 'High Throughput', '-', '-', '-', '-', 'BIOGRID', '-',
    '-', '-', 'P08107', 'A8K5I0', 'NP_005337', *NO_TERMS,
    'Severe acute respiratory syndrome coronavirus 2', HUMAN,
]
# Line 783431: a genetic interaction with phenotype terms.
SYNTHETIC_LETHALITY = [
    '1278840', '55294', '701', '120581', '107166', '-', '-', 'FBXW7', 'BUB1B',
    'AGO|CDC4|FBW6|FBW7|FBX30|FBXO30|FBXW6|SEL-10|SEL10|hAgo|hCdc4',
    'BUB1beta|BUBR1|Bub1A|MAD3L|MVA1|SSK1|hBUBR1', 'Synthetic Lethality',
    'genetic', 'Bailey ML (2015)', 'PUBMED:26354767', '9606', '9606',
    'Low Throughput', '-', '-', 'shRNA knockdown in FBW7 -/-|shRNA screen',
    '-', 'BIOGRID', 'Q969H0', 'S4R3U4|G0Z2K0',
    'NP_001013433|NP_361014|NP_060785|NP_001243998', 'O60566', '-',
    'NP_001202', 'HP:0001507', 'Growth abnormality', 'phenotype',
    'HP:0001507', 'Growth abnormality', '-', HUMAN, HUMAN,
]
# Line 166018: a fission yeast gene rescued by a human gene.
DOSAGE_RESCUE = [
    '300037', '2540618', '6047', '277144', '111974', 'SPBC3D6.11c',
    'RES4-26', 'slx8', 'RNF4', '-', 'SLX5|SNURF', 'Dosage Rescue', 'genetic',
    'Prudden J (2007)', 'PUBMED:17762865', '284812', '9606', 'Low Throughput',
    '-', '-', '-', '-', 'BIOGRID', 'P87176', '-', 'NP_595523', 'P78317', '-',
    'NP_001171939|NP_001171938|NP_002929', 'APO:0000147|APO:0000087',
    'heat sensitivity|resistance to chemicals', 'phenotype|phenotype', '-|-',
    '-|-', 'wild type|wild type', 'Schizosaccharomyces pombe (972h)', HUMAN,
]
# Line 48576: a two-hybrid pair of baker's yeast proteins.
YEAST_TWO_HYBRID = [
    '140962', '850822', '854357', '31400', '34580', 'YLR131C', 'YOR185C',
    'ACE2', 'GSP2', 'DNA-binding transcription factor ACE2|L000000020',
    'CNR2|Ran GTPase GSP2|L000000737', 'Two-hybrid', 'physical',
    'Jensen TH (2000)', 'PUBMED:11027275', '559292', '559292',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'P21192', '-',
    'NP_013232', 'P32836', '-', 'NP_014828', '-', '-', '-', '-', '-', '-',
    'Saccharomyces cerevisiae (S288c)', 'Saccharomyces cerevisiae (S288c)',
]
# Line 123858: a mouse protein captured with a human one.
MOUSE_HUMAN = [
    '243914', '15951', '672', '200523', '107140', '-', '-', 'Ifi204', 'BRCA1',
    'Ifi16|p204', 'BRCAI|BRCC1|BROVCA1|FANCS|IRIS|PNCA4|PPP1R53|PSCP|RNF53',
    'Affinity Capture-Western', 'physical', 'Aglipay JA (2003)',
    'PUBMED:14654789', '10090', '9606', 'Low Throughput', '-', '-',
    'BRCA1 (aa 502-802) interacted with the IFI16 Pyrin domain (aa 1-130);BRCA1-IFI16 role in p53-mediated transmission of DNA damage signals and apoptosis',
    '-', 'BIOGRID', 'Q08619|P15092', 'E9QKL6', 'NP_032355', 'P38398',
    'E9PFZ0|A0A024R1V0', 'NP_009229|NP_009228|NP_009230|NP_009231|NP_009225',
    '-', '-', '-', '-', '-', '-', 'Mus musculus', 'Homo sapiens',
]
# Line 295941: a yeast study with a plain score.
PLAIN_SCORE = [
    '439165', '852758', '851412', '33132', '31918', 'YGL119W', 'YDL143W',
    'COQ8', 'CCT4', 'ABC1|protein kinase COQ8|L000000010',
    'ANC2|TCP4|L000000087', 'Affinity Capture-MS', 'physical',
    'Breitkreutz A (2010)', 'PUBMED:20489023', '559292', '559292',
    'High Throughput', '0.946', '-', '-', '-', 'BIOGRID', 'P27697', '-',
    'NP_011396', 'P39078', '-', 'NP_010138', '-', '-', '-', '-', '-', '-',
    'Saccharomyces cerevisiae (S288c)', 'Saccharomyces cerevisiae (S288c)',
]
# Line 1933523: a negative score.
NEGATIVE_SCORE = [
    '2846057', '121512', '5339', '125734', '111355', '-', '-', 'FGD4', 'PLEC',
    'CMT4H|FRABP|ZFYVE6', 'EBS1|EBSO|HD1|LGMD2Q|PCN|PLEC1|PLEC1b|PLTN',
    'Affinity Capture-MS', 'physical', 'Mueller PM (2020)', 'PUBMED:32203420',
    '9606', '9606', 'High Throughput', '-0.107558193', '-',
    'prey proteins were classified into confidence-level categories of gold (highest confidence), silver or bronze based on various criteria|tagged bait proteins were affinity purified and associated proteins identified by mass spectrometry|the score provided is the CompPASS Z-score|this interaction was classified by the authors as being in the bronze list',
    '-', 'BIOGRID', 'Q96M96', 'E9PJX4|F8W1R0|B7Z493', 'NP_640334', 'Q15149',
    '-',
    'NP_000436|NP_958785|NP_958784|NP_958786|NP_958781|NP_958780|NP_958783|NP_958782',
    '-', '-', '-', '-', '-', '-', 'Homo sapiens', 'Homo sapiens',
]
# Line 1496517: a score with an exponent.
EXPONENT_SCORE = [
    '2360999', '1355', '56829', '107747', '121203', '-', 'PRO1677', 'COX15',
    'ZC3HAV1', 'CEMCOX2', 'ARTD13|FLB6421|PARP13|ZAP|ZC3H2|ZC3HDC2',
    'Affinity Capture-MS', 'physical', 'Malty RH (2017)', 'PUBMED:29128334',
    '9606', '9606', 'High Throughput', '8.90312e-05', '-',
    'HEK293 cell line (Antibody; IP/MS)|HEK293 cell line (FLAG TAGGED; AP/MS)|Normalized Spectral Abundance Factor (NSAF) scores were associated with each mitochondrial protein (MP) interaction. High confidence MP interactions were filtered at a confidence threshold of 90%.',
    '-', 'BIOGRID', 'Q7KZN9', '-',
    'NP_001307903|NP_004367|NP_510870|NP_001307905|NP_001307904', 'Q7Z2W4',
    'C9J6P4', 'NP_064504|NP_001350420|NP_078901', '-', '-', '-', '-', '-', '-',
    'Homo sapiens', 'Homo sapiens',
]
# Line 503966: co-fractionation with a score.
CO_FRACTIONATION_SCORE = [
    '740709', '7307', '11338', '113157', '116466', '-', '-', 'U2AF1', 'U2AF2',
    'FP793|RN|RNU2AF1|U2AF35|U2AFBP', 'U2AF65', 'Co-fractionation', 'physical',
    'Havugimana PC (2012)', 'PUBMED:22939629', '9606', '9606',
    'High Throughput', '0.874', '-', 'Denoised score >= 0.75', '-', 'BIOGRID',
    'Q01081', 'Q7Z780', 'NP_001020374|NP_001020375|NP_006749', 'P26368', '-',
    'NP_001012496|NP_009210', '-', '-', '-', '-', '-', '-', 'Homo sapiens',
    'Homo sapiens',
]
# Line 120537: two proteins seen in the same place.
CO_LOCALIZATION = [
    '240138', '27352', '4771', '118164', '110844', 'RP5-1042K10.9', '-',
    'SGSM3', 'NF2', 'MAP|RABGAP5|RUSC3|RUTBC3|RabGAP-5|rabGAPLP',
    'ACN|BANF|SCH', 'Co-localization', 'physical', 'Lee IK (2004)',
    'PUBMED:15541357', '9606', '9606', 'Low Throughput', '-', '-', '-', '-',
    'BIOGRID', 'Q96HU1', 'B9A6J5|B4DMS2|B4DVE3',
    'NP_001336972|NP_001336969|NP_001336970|NP_001336971|NP_001336976|NP_001336968|NP_001288778|NP_001336975|NP_001336974|NP_001336973|NP_056520',
    'P35240', 'A0A024R1J9|A0A024R1D9|A0A024R1J8|A0A024R1F6|A0A024R1I0',
    'NP_000259|NP_861971|NP_861970|NP_861546|NP_861966|NP_861967|NP_057502|NP_861968|NP_861969',
    '-', '-', '-', '-', '-', '-', 'Homo sapiens', 'Homo sapiens',
]
# Line 1448657: a fly pair that only has TrEMBL and RefSeq accessions of the first gene.
FLY_TWO_HYBRID = [
    '2311451', '117365', '34286', '72885', '60384', 'Dmel_CG1898',
    'Dmel_CG3959', 'HBS1', 'pelo', 'CG1898|DmelCG1898',
    'CG3959|DmelCG3959|ms(2)01559', 'Two-hybrid', 'physical', 'Yang F (2015)',
    'PUBMED:26124316', '7227', '7227', 'Low Throughput', '-', '-', '-', '-',
    'FLYBASE', '-', 'Q9W074', 'NP_652729|NP_001286906', 'P48612', '-',
    'NP_476982', '-', '-', '-', '-', '-', '-', 'Drosophila melanogaster',
    'Drosophila melanogaster',
]
# Line 1754637: a co-crystal structure of a viral protein with itself, from a preprint.
VIRAL_PREPRINT = [
    '2658320', '43740578', '43740578', '4383854', '4383854', 'GU280_gp01_nsp5',
    'GU280_gp01_nsp5', 'nsp5', 'nsp5',
    'ORF1ab|ORF1ab-nsp5|SARS-CoV2 nsp5|SARS-CoV-2 nsp5|Mpro|3CLpro|3C-like proteinase|R1AB_SARS2|PRO_0000449623',
    'ORF1ab|ORF1ab-nsp5|SARS-CoV2 nsp5|SARS-CoV-2 nsp5|Mpro|3CLpro|3C-like proteinase|R1AB_SARS2|PRO_0000449623',
    'Co-crystal Structure', 'physical', 'Dai W (2020)',
    'DOI:10.1101/2020.03.25.996348', '2697049', '2697049', 'Low Throughput',
    '-', '-', '-', '-', 'BIOGRID', 'P0DTD1', '-', 'YP_009725301', 'P0DTD1',
    '-', 'YP_009725301', '-', '-', '-', '-', '-', '-',
    'Severe acute respiratory syndrome coronavirus 2',
    'Severe acute respiratory syndrome coronavirus 2',
]
# Line 2752705: a viral protein displayed against a human receptor.
SURFACE_DISPLAY = [
    '3776731', '43740568', '59272', '4383848', '121864', 'GU280_gp02',
    'UNQ868/PRO1885', 'S', 'ACE2',
    'spike|SARS-CoV2 S|SARS-CoV2 spike|S protein|surface|SARS-CoV-2 spike|SARS-CoV-2 S|SPIKE_SARS2|PRO_0000449646',
    'ACEH', 'Surface Display', 'physical', 'Tsujino S (2024)',
    'PUBMED:38961765', '2697049', '9606', 'Low Throughput', '-', '-', '-', '-',
    'BIOGRID', 'P0DTC2', '-', 'YP_009724390', 'Q9BYF1', '-', 'NP_068576', '-',
    '-', '-', '-', '-', '-', 'Severe acute respiratory syndrome coronavirus 2',
    'Homo sapiens',
]
# Line 2824277: a thermal shift assay, the only one in the file.
THERMAL_SHIFT = [
    '3853326', '851055', '855341', '31604', '35475', 'YLR342W', 'YMR295C',
    'FKS1', 'YMR295C',
    'CND1|CWH53|ETG1|GSC1|PBR1|1,3-beta-D-glucan synthase|L000000445|L000000732',
    'IBI2', 'Thermal Shift Assay', 'physical', 'Li J (2025)',
    'PUBMED:41354657', '559292', '559292', 'Low Throughput', '-', '-', '-',
    '-', 'BIOGRID', 'P38631', '-', 'NP_013446', 'Q03559', '-', 'NP_014024',
    '-', '-', '-', '-', '-', '-', 'Saccharomyces cerevisiae (S288c)',
    'Saccharomyces cerevisiae (S288c)',
]
# Line 60619: a yeast dosage lethality.
DOSAGE_LETHALITY = [
    '153745', '854531', '855574', '34734', '35680', 'YOR349W', 'YNL148C',
    'CIN1', 'ALF1', 'L000000336', 'L000003527', 'Dosage Lethality', 'genetic',
    'Feierbach B (1999)', 'PUBMED:9885248', '559292', '559292',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'P40987', '-',
    'NP_014994', 'P53904', '-', 'NP_014251', 'APO:0000112', 'inviable',
    'phenotype', '-', '-', '-', 'Saccharomyces cerevisiae (S288c)',
    'Saccharomyces cerevisiae (S288c)',
]
# Line 96172: a negative score of an E-MAP.
GENETIC_NEGATIVE_SCORE = [
    '206363', '850369', '850333', '30990', '30959', 'YCR011C', 'YCL025C',
    'ADP1', 'AGP1', 'putative ATP-dependent permease ADP1|L000000049',
    'YCC5|amino acid transporter AGP1|L000003271', 'Negative Genetic',
    'genetic', 'Schuldiner M (2005)', 'PUBMED:16269340', '559292', '559292',
    'High Throughput', '-5.6431', '-',
    'An Epistatic MiniArray Profile (E-MAP) analysis was used to quantitatively score genetic interactions based on fitness defects estimated from the colony size of double versus single mutants.  Genetic interactions were considered significant if they had an S score > 2.5 for positive interactions (suppression)  and S score < -2.5 for negative interactions (synthetic sick/lethality).',
    '-', 'BIOGRID', 'P25371', '-', 'NP_009937', 'P25376', '-', 'NP_009905',
    'APO:0000063', 'colony size', 'phenotype', '-', '-', '-',
    'Saccharomyces cerevisiae (S288c)', 'Saccharomyces cerevisiae (S288c)',
]
# Line 96173: a positive score of an E-MAP.
GENETIC_PLAIN_SCORE = [
    '206364', '850369', '850513', '30990', '31115', 'YCR011C', 'YFL031W',
    'ADP1', 'HAC1', 'putative ATP-dependent permease ADP1|L000000049',
    'ERN4|IRE15|transcription factor HAC1|L000002611', 'Positive Genetic',
    'genetic', 'Schuldiner M (2005)', 'PUBMED:16269340', '559292', '559292',
    'High Throughput', '2.5219', '-',
    'An Epistatic MiniArray Profile (E-MAP) analysis was used to quantitatively score genetic interactions based on fitness defects estimated from the colony size of double versus single mutants.  Genetic interactions were considered significant if they had an S score > 2.5 for positive interactions (suppression)  and S score < -2.5 for negative interactions (synthetic sick/lethality).',
    '-', 'BIOGRID', 'P25371', '-', 'NP_009937', 'P41546', '-', 'NP_116622',
    'APO:0000063', 'colony size', 'phenotype', '-', '-', '-',
    'Saccharomyces cerevisiae (S288c)', 'Saccharomyces cerevisiae (S288c)',
]
# Line 2566213: a score with an exponent, of human genes.
GENETIC_EXPONENT_SCORE = [
    '3584384', '805', '808', '107256', '107259', '-', '-', 'CALM2', 'CALM3',
    'CAMII|PHKD|PHKD2', 'HEL-S-72|PHKD|PHKD3', 'Negative Genetic', 'genetic',
    'Ito T (2021)', 'PUBMED:34857952', '9606', '9606', 'High Throughput',
    '6.50523e-05', '-',
    'CRISPR GI screen|Cell Line:  PATU8988S_PANCREAS score (0.0368732253612743)|Cell Line:  PK1_PANCREAS score (0.000026907546703447)|Cell Line: A549_LUNG score (0.000065052302380019)|Cell Line: GI1_CENTRAL_NERVOUS_SYSTEM score (0.0408961003123878)|Cell Line: HSC5_SKIN score (0.0069166791692945)|Cell Line: IPC298_SKIN score (0.0016831373997568)|Cell Line: MEL202_UVEA score (0.0006778402486547)|Experimental Setup: Timecourse-Synthetic Lethality|GIST: A-phenotypic negative genetic interaction|Library: Digenic Paralog CRISPR library|Significance Threshold: GEMINI FDR < 0.05',
    '-', 'BIOGRID', 'P62158', 'B4DJ51',
    'NP_001292555|NP_001292554|NP_001734|NP_001292553', 'P62158',
    'B4DJ51|Q9BRL5',
    'NP_001316855|NP_001316854|NP_001316853|NP_001316852|NP_001316851|NP_001316850|NP_005175',
    '-', '-', '-', '-', '-', '-', 'Homo sapiens', 'Homo sapiens',
]

# Rows with ontology terms and modifications, copied from the file.
# Line 493797: HEK-293 cell line and a cell type.
HEK_CELL_LINE_AND_CELL_TYPE = [
    '725220', '79400', '3308', '122660', '109540', '-', '-', 'NOX5', 'HSPA4',
    '-', 'APG-2|HEL-S-5a|HS24/P52|HSPH2|RY|hsp70|hsp70RY',
    'Affinity Capture-Western', 'physical', 'Chen F (2012)', 'PUBMED:23023377',
    '9606', '9606', 'Low Throughput', '-', '-',
    'Exogenous expression of bait (HEK-293 cells)', '-', 'BIOGRID', 'Q96PH1',
    'A3QRJ0', 'NP_078781|NP_001171708|NP_001171709', 'P34932', '-',
    'NP_002145', 'BTO:0000007|BTO:0001025',
    'HEK-293 cell|peripheral blood mononuclear cell', 'cell line|cell type',
    '-|-', '-|-', '-|-', 'Homo sapiens', 'Homo sapiens',
]
# Line 647352: HeLa cell with the disease it derives from as qualifier.
HELA_WITH_DISEASE = [
    '926198', '6206', '55720', '112120', '120842', '-', '-', 'RPS12', 'TSR1',
    'S12', '-', 'Co-fractionation', 'physical', 'Kristensen AR (2012)',
    'PUBMED:22863883', '9606', '9606', 'High Throughput', '-', '-', '-', '-',
    'BIOGRID', 'P25398', '-', 'NP_001007', 'Q2NL82', '-', 'NP_060598',
    'BTO:0000567', 'HeLa cell', 'cell line', 'DOID:3702',
    'cervical adenocarcinoma', '-', 'Homo sapiens', 'Homo sapiens',
]
# Line 2757096: HEK-293 cell and kidney that BioGRID labels "phenotype", with a cell type.
HEK_AS_PHENOTYPE = [
    '3781150', '-', '550', '4383959', '107031', 'VE5_HPV31', '-', 'E5', 'AUP1',
    'E5 protein', '-', 'Affinity Capture-MS', 'physical', 'Eckhardt M (2018)',
    'PUBMED:30209081', '10585', '9606', 'High Throughput', '-', '-', '-', '-',
    'BIOGRID', 'P17385', '-', '-', 'Q9Y679', '-', 'NP_853553',
    'BTO:0000671|CL:0002518|BTO:0000007',
    'kidney|kidney epithelial cell|HEK-293 cell',
    'phenotype|cell type|phenotype', '-|-|BTO:0002733',
    '-|-|embryonic kidney cell line', '-|-|-', 'Human papillomavirus (31)',
    'Homo sapiens',
]
# Line 2815807: a cell line that BioGRID labels "phenotype".
PC3_AS_PHENOTYPE = [
    '3844416', '10481', '3500', '115744', '109721', '-', '-', 'HOXB13',
    'IGHG1', 'PSGD', '-', 'Affinity Capture-MS', 'physical',
    'Ersoy-Fazlioglu B (2025)', 'PUBMED:41343677', '9606', '9606',
    'High Throughput', '-', '-', '-', '-', 'BIOGRID', 'Q92826', 'Q4KR72',
    'NP_006352', '-', '-', '-', 'BTO:0001061', 'PC-3 cell', 'phenotype', '-',
    '-', '-', 'Homo sapiens', 'Homo sapiens',
]
# Line 631954: a chemical.
CHEMICAL = [
    '907739', '4056', '4056', '110234', '110234', '-', '-', 'LTC4S', 'LTC4S',
    '-', '-', 'Co-crystal Structure', 'physical', 'Ago H (2007)',
    'PUBMED:17632548', '9606', '9606', 'Low Throughput', '-', '-', 'trimer',
    '-', 'BIOGRID', 'Q16873', '-', 'NP_665874', 'Q16873', '-', 'NP_665874',
    'CHEBI:16856', 'glutathione', 'chemical', '-', '-', '-', 'Homo sapiens',
    'Homo sapiens',
]
# Line 469468: a biochemical activity without modification.
NO_MODIFICATION = [
    '684318', '25508', '7316', '247540', '113164', '-', '-', 'Pam', 'UBC', '-',
    'HMG20', 'Biochemical Activity', 'physical', 'Chew GH (2005)',
    'PUBMED:16098968', '10116', '9606', 'Low Throughput', '-',
    'No Modification', '-', '-', 'BIOGRID', 'P14925', '-', 'NP_037132',
    'P0CG48', '-', 'NP_066289', '-', '-', '-', '-', '-', '-',
    'Rattus norvegicus', 'Homo sapiens',
]
# Line 69214: two phenotypes of a partial rescue.
PARTIAL_RESCUE = [
    '164095', '854535', '850800', '34738', '31383', 'YOR353C', 'YLR110C',
    'SOG2', 'CCW12', '-', '-', 'Dosage Rescue', 'genetic',
    'Kurischko C (2005)', 'PUBMED:15972461', '559292', '559292',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'Q08817', '-',
    'NP_014998', 'Q12127', '-', 'NP_013211', 'APO:0000111|APO:0000049',
    'viability|morphology', 'phenotype|phenotype', '-|-', '-|-',
    'partial rescue|partial rescue', 'Saccharomyces cerevisiae (S288c)',
    'Saccharomyces cerevisiae (S288c)',
]
# Line 338437: a zebrafish anatomy term with a PATO qualifier.
QUALIFIED_ANATOMY = [
    '505707', '30465', '117510', '78646', '79558', '-', '-', 'lhx5', 'sfrp5',
    'lim5', 'D144|frp1|zgc:101849', 'Phenotypic Suppression', 'genetic',
    'Peng G (2006)', 'PUBMED:16854974', '7955', '7955', 'Low Throughput', '-',
    '-', '-', '-', 'BIOGRID', 'P52889', '-', 'NP_571293', '-', 'Q6YNR8',
    'NP_571933', 'ZFA:0000107', 'eye', 'phenotype', 'PATO:0000646',
    'malformed', '-', 'Danio rerio', 'Danio rerio',
]
# Line 336297: a term with a PATO qualifier and a type.
QUALIFIED_TYPED = [
    '503356', '386629', '30701', '88941', '78856', '-', '-', 'tob1a', 'lef1',
    'tob1', '-', 'Phenotypic Suppression', 'genetic', 'Xiong B (2006)',
    'PUBMED:16890162', '7955', '7955', 'Low Throughput', '-', '-',
    'Overexpression of lef1 alleviated the ventralizing effect of tob1a overexpression.',
    '-', 'BIOGRID', '-', 'Q7SZR1', 'NP_938177', '-', 'Q9W7C0', 'NP_571501',
    'ZFA:0001105', 'embryonic structure', 'phenotype', 'PATO:0000636',
    'ventralized', 'wild type', 'Danio rerio', 'Danio rerio',
]
# Line 1564075: a human CRISPR screen with a cell and a phenotype.
CRISPR_SCREEN = [
    '2456438', '819', '819', '107269', '107269', '-', '-', 'CAMLG', 'CAMLG',
    'CAML', 'CAML', 'Negative Genetic', 'genetic', 'Horlbeck MA (2018)',
    'PUBMED:30033366', '9606', '9606', 'High Throughput', '-5.268160661', '-',
    'CRISPR GI screen|Cell Line: K562 EFO:0002067/Jurkat EFO:0002796|Experimental Setup: Timecourse|GIST: A-phenotypic negative/positive genetic interaction|Interactions in this CRISPR interference (CRISPRi) analysis were considered to be significant when GI <= -3 (negative genetic interaction) or GI >= 3 (positive genetic interaction).|K562 cell line Replicate Average GI score = -5.268160661|Library: CRISPRi v1|Significance Threshold: (positive genetic interaction) 3<GI<-3 (negative genetic interaction)',
    '-', 'BIOGRID', 'P49069', '-', 'NP_001736', 'P49069', '-', 'NP_001736',
    'BTO:0000664|HP:0001507', 'K-562 cell|Growth abnormality',
    'cell type|phenotype', '-|-', '-|-', '-|-', 'Homo sapiens', 'Homo sapiens',
]
# Line 2943997: a type on a term that is not a phenotype.
TYPED_PROCESS = [
    '3978428', '3456', '55681', '109678', '120810', 'RP11-113D19.1', '-',
    'IFNB1', 'SCYL2', 'IFB|IFF|IFNB', 'CVAK104', 'Phenotypic Suppression',
    'genetic', 'Miyakawa K (2012)', 'PUBMED:23047923', '9606', '9606',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'P01574', '-',
    'NP_002167', 'Q6P3W7', '-', 'NP_060458', 'GO:0019076',
    'viral release from host cell', 'biological process', '-', '-',
    'wild type', 'Homo sapiens', 'Homo sapiens',
]

# Rows with ontology terms of the kinds the other rows lack.
# Line 2946264: a CRISPR activation screen with a phenotype of the Gene Ontology.
GO_PHENOTYPE = [
    '4102079', '64843', '4066', '122318', '110244', '-', '-', 'ISL2', 'LYL1',
    '-', 'bHLHa18', 'Positive Genetic', 'genetic', 'Norman TM (2019)',
    'PUBMED:31395745', '9606', '9606', 'High Throughput', '3.758534192', '-',
    'CRISPR GI Activation screen|Cell Line:K-562 (BTO:0000664/CELLOSAURUS:CVCL_0004)|Experimental Setup: In vitro proliferation|GIST: A-phenotypic positive genetic interaction|Library: CRISPRa GI library|Significance Threshold: -2.77>GI score>2.83 (Figure 1F)',
    '-', 'BIOGRID', 'Q96A47', '-', 'NP_665804', 'P12980', '-', 'NP_005574',
    'GO:0008283', 'cell proliferation', 'phenotype', '-', '-', '-',
    'Homo sapiens', 'Homo sapiens',
]
# Line 2285734: a phenotype of the Phenotype And Trait Ontology.
PATO_PHENOTYPE = [
    '3202817', '25842', '55723', '117368', '120845', 'RP3-329L24.1', '-',
    'ASF1A', 'ASF1B', 'CGI-98|CIA|HSPC146', 'CIA-II', 'Synthetic Lethality',
    'genetic', 'De Kegel B (2021)', 'PUBMED:34529928', '9606', '9606',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'Q9Y294', '-',
    'NP_054753', 'Q9NVP2', 'A0A024R7G4', 'NP_060624', 'PATO:0000169',
    'viability', 'phenotype', '-', '-', '-', 'Homo sapiens', 'Homo sapiens',
]
# Line 338001: a phenotype of the Mammalian Phenotype Ontology.
MOUSE_PHENOTYPE = [
    '505265', '13543', '13544', '199343', '199344', 'RP23-172M21.1', '-',
    'Dvl2', 'Dvl3', '-', '-', 'Synthetic Lethality', 'genetic',
    'Etheridge SL (2008)', 'PUBMED:19008950', '10090', '10090',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'Q60838', '-',
    'NP_031914', 'Q61062', 'E9Q967', 'NP_031915', 'MP:0008762',
    'embryonic lethality', 'phenotype', '-', '-', '-', 'Mus musculus',
    'Mus musculus',
]
# Line 2846060: a phenotype of the C. elegans phenotype ontology.
WORM_PHENOTYPE = [
    '3876261', '176076', '173375', '41284', '38756', 'CELE_ZK783.1',
    'CELE_F31C3.3', 'fbn-1', 'F31C3.3', '-', '-', 'Phenotypic Enhancement',
    'genetic', 'Dziurdzik SK (2026)', 'PUBMED:41746219', '6239', '6239',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', '-', 'Q23587',
    'NP_498670', '-', 'O62191', 'NP_493625', 'WBPhenotype:0000520',
    'morphology variant', 'phenotype', '-', '-', '-', 'Caenorhabditis elegans',
    'Caenorhabditis elegans',
]
# Line 1371456: a cell line that BioGRID labels phenotype next to two phenotypes.
CELL_LINE_AS_PHENOTYPE_GENETIC = [
    '2208433', '1019', '7150', '107454', '113003', '-', 'RP3-511B24.1', 'CDK4',
    'TOP1', 'CMM3|PSK-J3', 'TOPI', 'Negative Genetic', 'genetic',
    'Shen JP (2017)', 'PUBMED:28319113', '9606', '9606', 'High Throughput',
    '-', '-',
    'CRISPR GI screen|Cell Line: A-549 EFO:0001086|Experimental Setup: Timecourse|GIST: A-phenotypic negative genetic interaction|Library: Dual-guide CRISPRn library|Significance Threshold:FDR ~ 0.3',
    '-', 'BIOGRID', 'P11802', 'A0A024RBB6', 'NP_000066', 'P11387', '-',
    'NP_003277', 'PATO:0000169|BTO:0000018|HP:0001507',
    'viability|A-549 cell|Growth abnormality', 'phenotype|phenotype|phenotype',
    '-|-|-', '-|-|-', '-|-|-', 'Homo sapiens', 'Homo sapiens',
]
# Line 1371372: HeLa cell that BioGRID labels phenotype next to two phenotypes.
HELA_AS_PHENOTYPE_GENETIC = [
    '2208349', '2625', '7465', '108895', '113303', '-', '-', 'GATA3', 'WEE1',
    'HDR|HDRS', 'WEE1A|WEE1hu', 'Negative Genetic', 'genetic',
    'Shen JP (2017)', 'PUBMED:28319113', '9606', '9606', 'High Throughput',
    '-', '-',
    'CRISPR GI screen|Cell Line: HeLa EFO:0001185|Experimental Setup: Timecourse|GIST: A-phenotypic negative genetic interaction|Library: Dual-guide CRISPRn library|Significance Threshold:FDR ~ 0.3',
    '-', 'BIOGRID', 'P23771', '-', 'NP_001002295|NP_002042', 'P30291',
    'Q86V29', 'NP_003381|NP_001137448', 'PATO:0000169|HP:0001507|BTO:0000567',
    'viability|Growth abnormality|HeLa cell', 'phenotype|phenotype|phenotype',
    '-|-|-', '-|-|-', '-|-|-', 'Homo sapiens', 'Homo sapiens',
]
# Line 496654: a worm cell that BioGRID labels phenotype, as qualifier of itself.
MUSCLE_CELL_AS_PHENOTYPE = [
    '728522', '175519', '172737', '40756', '38168', 'C36E8.5', 'CELE_ZK858.4',
    'tbb-2', 'mel-26', 'CELE_C36E8.5|TUBB8|rot-2', '-',
    'Phenotypic Suppression', 'genetic', 'Wilson KJ (2012)', 'PUBMED:22621901',
    '6239', '6239', 'Low Throughput', '-', '-', 'figure 6.', '-', 'BIOGRID',
    'P52275', '-', 'NP_497806', 'Q94420', '-', 'NP_492449', 'WBbt:0006804',
    'body wall muscle cell', 'phenotype', 'WBbt:0006804',
    'body wall muscle cell', '-', 'Caenorhabditis elegans',
    'Caenorhabditis elegans',
]
# Line 648715: a chemical and two phenotypes.
GENETIC_CHEMICAL = [
    '927565', '852028', '852793', '32477', '33163', 'YDR419W', 'YGL087C',
    'RAD30', 'MMS2', 'DBH1|L000003523',
    'E2  ubiquitin-conjugating protein MMS2|L000004015',
    'Synthetic Growth Defect', 'genetic', 'Daraba A (2014)', 'PUBMED:24465179',
    '559292', '559292', 'Low Throughput', '-', '-', 'on MMS', '-', 'BIOGRID',
    'Q04049', '-', 'NP_010707', 'P53152', '-', 'NP_011428',
    'CHEBI:25255|APO:0000106|APO:0000087',
    'methyl methanesulfonate|vegetative growth|resistance to chemicals',
    'chemical|phenotype|phenotype', '-|-|-', '-|-|-', '-|-|-',
    'Saccharomyces cerevisiae (S288c)', 'Saccharomyces cerevisiae (S288c)',
]
# Line 606341: a tissue and a phenotype of mice.
GENETIC_TISSUE = [
    '872759', '11806', '11816', '198155', '198164', '-', '-', 'Apoa1', 'Apoe',
    'Alp-1|Apoa-1|Brp-14|Ltw-1|Lvtw-1|Sep-1|Sep-2|Sep2|apo-AI|apoA-I',
    'AI255918', 'Phenotypic Enhancement', 'genetic', 'Zhao Y (2005)',
    'PUBMED:15654758', '10090', '10090', 'Low Throughput', '-', '-',
    'apoA1-/- apoe-/- double mutants have much reduced LCAT activity compared to either single mutant',
    '-', 'BIOGRID', 'Q00623', '-', 'NP_033822', 'P08226', 'Q3TXU4',
    'NP_001292748|NP_001292773|NP_001292772|NP_033826',
    'BTO:0000131|MP:0000180',
    'blood plasma|abnormal circulating cholesterol level', 'tissue|phenotype',
    '-|-', '-|-', '-|-', 'Mus musculus', 'Mus musculus',
]
# Line 632703: a phenotype with two qualifiers.
QUALIFIED_APO_PHENOTYPE = [
    '908564', '2540618', '2539686', '277144', '276241', 'SPBC3D6.11c',
    'SPBC1734.06', 'slx8', 'rhp18', '-', '-', 'Phenotypic Enhancement',
    'genetic', 'Kohler JB (2013)', 'PUBMED:24265825', '284812', '284812',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'P87176', '-',
    'NP_595523', 'O74747', '-', 'NP_595423', 'APO:0000087',
    'resistance to chemicals', 'phenotype', 'APO:0000003^CHEBI:27656',
    'decreased^camptothecin', '-', 'Schizosaccharomyces pombe (972h)',
    'Schizosaccharomyces pombe (972h)',
]
# Line 2528181: a human phenotype with a PATO qualifier.
QUALIFIED_HP_PHENOTYPE = [
    '3541992', '997', '7326', '107432', '113174', '-', '-', 'CDC34', 'UBE2G1',
    'E2-CDC34|UBC3|UBCH3|UBE2R1', 'E217K|UBC7|UBE2G', 'Negative Genetic',
    'genetic', 'Hill S (2019)', 'PUBMED:31868589', '9606', '9606',
    'Low Throughput', '-', '-',
    'CRISPR GI screen|UBE2G1 is synthetic lethal in a UBE2R1/UBE2R2 background|genetic complex',
    '-', 'BIOGRID', 'P49427', 'A0A024R1Z1', 'NP_004350', 'P62253', '-',
    'NP_003333', 'HP:0001507', 'Growth abnormality', 'phenotype',
    'PATO:0000169', 'viability', '-', 'Homo sapiens', 'Homo sapiens',
]
# Line 2942386: a cell type, a cell line that BioGRID labels phenotype with a qualifier, a process and skin.
CELL_TYPE_AND_CELL_LINE_AS_PHENOTYPE_GENETIC = [
    '3976812', '5586', '25937', '111572', '117434', 'RP5-905H16.1', '-',
    'PKN2', 'WWTR1', 'PAK2|PRK2|PRKCL2|PRO2042|Pak-2|STK7', 'TAZ',
    'Phenotypic Suppression', 'genetic', 'Killarney ST (2025)',
    'PUBMED:39560431', '9606', '9606', 'Low Throughput', '-', '-', '-', '-',
    'BIOGRID', 'Q16513', '-',
    'NP_001307637|NP_001307636|NP_001307638|NP_006247', 'Q9GZV5', '-',
    'NP_056287|NP_001335291|NP_001161750|NP_001161752',
    'CL:0000148|BTO:0002131|GO:0010468|BTO:0001253',
    'melanocyte|SK-MEL-28 cell|regulation of gene expression|skin',
    'cell type|phenotype|biological process|phenotype', '-|BTO:0000849|-|-',
    '-|melanoma cell line|-|-', '-|-|wild type|-', 'Homo sapiens',
    'Homo sapiens',
]
# Line 608826: a disease and a cell line.
PHYSICAL_DISEASE = [
    '876673', '1437', '10949', '107824', '116149', '-', '-', 'CSF2', 'HNRNPA0',
    'GMCSF', 'HNRPA0', 'Affinity Capture-RNA', 'physical', 'Bollig F (2003)',
    'PUBMED:12565831', '9606', '9606', 'Low Throughput', '-', '-', '-', '-',
    'BIOGRID', 'P04141', '-', 'NP_000749', 'Q13151', '-', 'NP_006796',
    'DOID:299|BTO:0000567', 'adenocarcinoma|HeLa cell', 'disease|cell line',
    '-|-', '-|-', '-|-', 'Homo sapiens', 'Homo sapiens',
]
# Line 183492: a phenotype of the Ascomycete Phenotype Ontology on a bacterial interaction.
PHYSICAL_APO_PHENOTYPE = [
    '318348', '938167', '939930', '855787', '857505', 'BSU04720', 'BSU04710',
    'rsbW', 'rsbV', '-', '-', 'Co-purification', 'physical', 'Dufour A (1994)',
    'PUBMED:8144446', '224308', '224308', 'Low Throughput', '-', '-', '-', '-',
    'BIOGRID', 'P17904', '-', 'NP_388353', 'P17903', '-', 'NP_388352',
    'APO:0000080', 'stress resistance', 'phenotype', '-', '-', '-',
    'Bacillus subtilis (168)', 'Bacillus subtilis (168)',
]
# Line 599779: a disease, a cell type, a phenotype of the Human Phenotype Ontology and a process.
PHYSICAL_HP_PHENOTYPE = [
    '860063', '1277', '1278', '107674', '107675', '-', '-', 'COL1A1', 'COL1A2',
    'OI4', 'OI4', 'Co-purification', 'physical', 'Pace JM (2008)',
    'PUBMED:18375391', '9606', '9606', 'Low Throughput', '-', '-', '-', '-',
    'BIOGRID', 'P02452', '-', 'NP_000079', 'P08123', '-', 'NP_000080',
    'DOID:12347|BTO:0004419|HP:0003784|GO:0030199',
    'osteogenesis imperfecta|dermal fibroblast|Type 1 collagen overmodification|collagen fibril organization',
    'disease|cell type|phenotype|biological process', '-|-|-|-', '-|-|-|-',
    '-|-|-|-', 'Homo sapiens', 'Homo sapiens',
]
# Line 604232: tissues, a process and a phenotype of the Mammalian Phenotype Ontology.
PHYSICAL_MP_PHENOTYPE = [
    '870545', '29581', '24180', '248214', '246370', '-', '-', 'Htr2b',
    'Agtr1a', '-', 'AT1|AT1A|AT1R|Agtr1', 'Affinity Capture-Western',
    'physical', 'Jaffre F (2009)', 'PUBMED:19023134', '10116', '10116',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'P30994', 'G3V8A0',
    'NP_058946', 'P25095', '-', 'NP_112247',
    'UBERON:0002349|GO:0050663|MP:0001625|BTO:0003093',
    'myocardium|cytokine secretion|cardiac hypertrophy|cardiofibroblast',
    'tissue|biological process|phenotype|tissue', '-|-|-|-', '-|-|-|-',
    '-|-|-|-', 'Rattus norvegicus', 'Rattus norvegicus',
]
# Line 639973: PC-3 cell labelled cell line, with a disease as qualifier.
PC3_CELL_LINE = [
    '918532', '10399', '6194', '115671', '112108', 'HLC7', 'RP11-513M16.6',
    'GNB2L1', 'RPS6', 'Gnb2-rs1|H12.3|HLC-7|PIG21|RACK1', 'S6',
    'Co-fractionation', 'physical', 'Marques N (2013)', 'PUBMED:23455324',
    '9606', '9606', 'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'P63244',
    'E9KL35', 'NP_006089', 'P62753', 'A2A3R6', 'NP_001001', 'BTO:0001061',
    'PC-3 cell', 'cell line', 'DOID:2526', 'prostate adenocarcinoma', '-',
    'Homo sapiens', 'Homo sapiens',
]
# Line 2913183: kidney, HEK-293T cell and a kidney cell that BioGRID labels phenotype.
KIDNEY_AS_PHENOTYPE = [
    '3946388', '9889', '65095', '115219', '122394', '-', '-', 'ZBED4', 'KRI1',
    '-', '-', 'Cross-Linking-MS (XL-MS)', 'physical', 'Zhao L (2022)',
    'PUBMED:36479438', '9606', '9606', 'High Throughput', '-', '-', '-', '-',
    'BIOGRID', 'O75132', '-', 'NP_055653', 'Q8N9T8', '-', 'NP_075384',
    'BTO:0000671|BTO:0002181|CL:0002518',
    'kidney|HEK-293T cell|kidney epithelial cell',
    'phenotype|phenotype|phenotype', '-|BTO:0002733|-',
    '-|embryonic kidney cell line|-', '-|-|-', 'Homo sapiens', 'Homo sapiens',
]
# Line 2911504: a leukocyte (CL) and K-562 cell and blood that BioGRID labels phenotype.
CELL_TYPE_AS_PHENOTYPE = [
    '3944709', '6124', '539', '112044', '107021', '-', '-', 'RPL4', 'ATP5O',
    'L4', 'ATPO|OSCP', 'Cross-Linking-MS (XL-MS)', 'physical',
    'Ryl PSJ (2020)', 'PUBMED:31746214', '9606', '9606', 'High Throughput',
    '-', '-', '-', '-', 'BIOGRID', 'P36578', '-', 'NP_000959', 'P48047', '-',
    'NP_001688', 'CL:0000738|BTO:0000664|BTO:0000089',
    'leukocyte|K-562 cell|blood', 'phenotype|phenotype|phenotype',
    '-|BTO:0002580|-', '-|chronic myeloid leukemia cell line|-', '-|-|-',
    'Homo sapiens', 'Homo sapiens',
]
# Line 605988: a tissue, a development stage and a cell type with a qualifier.
TISSUE_DEVELOPMENT_AND_CELL_TYPE = [
    '872330', '857', '4846', '107305', '110909', '-', '-', 'CAV1', 'NOS3',
    'BSCL3|CGL3|LCCNS|MSTP085|PPH3|VIP21', 'ECNOS|eNOS',
    'Affinity Capture-Western', 'physical', 'Linton EA (2003)',
    'PUBMED:12852865', '9606', '9606', 'Low Throughput', '-', '-', '-', '-',
    'BIOGRID', 'Q03135', 'A0A024R757|Q2TNI1|Q59E85|Q7Z4F3|A9XTE5',
    'NP_001744|NP_001166368|NP_001166366|NP_001166367', 'P29474',
    'A0S0A6|E7ESA7', 'NP_001153581|NP_000594|NP_001153583|NP_001153582',
    'UBERON:0000319|UBERON:0001987|BTO:0005003',
    'cytotrophoblast|placenta|villous trophoblast',
    'tissue|development|cell type', '-|-|BTO:0001078', '-|-|placenta', '-|-|-',
    'Homo sapiens', 'Homo sapiens',
]
# Line 614492: BioGRID's placeholder accession as a cell type.
PLACEHOLDER_CELL_TYPE = [
    '886743', '3665', '9111', '109873', '114561', '-', '-', 'IRF7', 'NMI',
    'IRF-7H|IRF7A|IRF7B|IRF7C|IRF7H', '-', 'Affinity Capture-Western',
    'physical', 'Wang J (2013)', 'PUBMED:23956435', '9606', '9606',
    'Low Throughput', '-', '-', '-', '-', 'BIOGRID', 'Q92985',
    'A0A024RCA2|M9RSF4', 'NP_004022|NP_004020|NP_001563', 'Q13287', 'Q8WTW2',
    'NP_004679', 'BIOGRID:0000011|BTO:0002181', 'BRENDA tissue|HEK-293T cell',
    'cell type|cell type', '-|-', '-|-', '-|-', 'Homo sapiens', 'Homo sapiens',
]
# Line 1368185: BioGRID's placeholder accession as a phenotype.
PLACEHOLDER_PHENOTYPE = [
    '2205106', '7408', '87', '113251', '106602', '-', '-', 'VASP', 'ACTN1',
    '-', 'BDPLT15', 'Affinity Capture-Western', 'physical', 'Bausen M (2006)',
    'PUBMED:16376568', '9606', '9606', 'Low Throughput', '-', '-', '-', '-',
    'BIOGRID', 'P50552', 'A0A024R0V4', 'NP_003361', 'P12814', 'A0A024R694',
    'NP_001093|NP_001123477|NP_001123476', 'BIOGRID:0000011', 'BRENDA tissue',
    'phenotype', '-', '-', '-', 'Homo sapiens', 'Homo sapiens',
]

PHYSICAL = [
    TWO_HYBRID, TWO_PROTEINS, PROTEIN_RNA, KINASE_ASSAY, CO_FRACTIONATION,
    BOTH_THROUGHPUTS, PREPRINT, VIRAL,
]
GENETIC = [SYNTHETIC_LETHALITY, DOSAGE_RESCUE]
MEMBER = 'BIOGRID-ALL-5.0.262.tab3.txt'


def table(*lines, header=HEADER):
    return '\n'.join('\t'.join(line) for line in (header, *lines)) + '\n'


def archive(text, member=MEMBER):
    """The opened all-organism archive: member name -> text handle."""
    return SimpleNamespace(result={member: io.StringIO(text)})


def read(*lines, dataset='physical', header=HEADER, **kwargs):
    parser = biogrid.resource.datasets()[dataset]._raw_parser
    return list(parser(archive(table(*lines, header=header), **kwargs)))


def changed(line, columns):
    """The row with the values of some columns replaced."""
    out = list(line)
    for name, value in columns.items():
        out[HEADER.index(name)] = value
    return out


def extract(rows, dataset='physical'):
    out = SilverExtractor('biogrid', dataset)
    for index, row in enumerate(rows):
        record = biogrid.resource.datasets()[dataset].mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def interactions(out):
    """The interaction relations of an extraction, without the links of genes
    to the phenotype terms of a genetic interaction."""
    return [r for r in out.relations if r.predicate != 'associated_with']


def participants(out):
    """The entities of an extraction that are not ontology terms."""
    return [e for e in out.entities.values() if e.entity_type != 'ontology_class']


def identifiers(entity):
    return [(i.type, i.value) for i in entity.identifiers]


def values(relation, term, scope='relation'):
    return [
        a['value'] for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


def scores(relation):
    """The numbers of the score annotation, as (value, unit, field)."""
    return [
        (a['quantity']['has_numeric_value'], a['quantity']['has_unit'],
         a['quantity']['source_field'])
        for a in relation.annotations
        if a['term'] == 'biogrid:score' and a['scope'] == 'relation'
    ]


def test_module_import_defines_the_datasets_without_downloading():
    assert list(biogrid.resource.datasets()) == ['physical', 'genetic']
    assert biogrid.config.id.name == 'BIOGRID'
    downloads = {d.download for d in biogrid.resource.datasets().values()}
    assert {d.url for d in downloads} == {biogrid._URL}
    assert {d.filename for d in downloads} == {'BIOGRID-ALL-5.0.262.tab3.zip'}
    # The release is pinned: the link is the file of the release archive.
    assert biogrid._URL == (
        'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/'
        'BIOGRID-5.0.262/BIOGRID-ALL-5.0.262.tab3.zip'
    )


def test_reader_takes_the_member_of_the_release_and_splits_by_system_type():
    physical = read(*PHYSICAL, *GENETIC)
    genetic = read(*PHYSICAL, *GENETIC, dataset='genetic')
    assert [row['#BioGRID Interaction ID'] for row in physical] == [
        line[0] for line in PHYSICAL
    ]
    assert [row['Experimental System'] for row in genetic] == [
        'Synthetic Lethality', 'Dosage Rescue',
    ]
    # The 37 columns of the file, the decoded ontology terms and their three
    # groupings.
    assert len(physical[0]) == 41
    assert set(HEADER) <= set(physical[0]) and 'Ontology Terms' in physical[0]
    assert physical[0]['Official Symbol Interactor B'] == 'FLNC'
    # No organism filter: the rows are of every organism in the file.
    rows = read(TWO_HYBRID, YEAST_TWO_HYBRID, MOUSE_HUMAN, VIRAL)
    assert [(r['Organism ID Interactor A'], r['Organism ID Interactor B'])
            for r in rows] == [
        ('9606', '9606'), ('559292', '559292'), ('10090', '9606'),
        ('2697049', '9606'),
    ]


def test_reader_rejects_another_archive_or_layout():
    with pytest.raises(ValueError, match='0 members named'):
        read(TWO_HYBRID, member='BIOGRID-ORGANISM-Homo_sapiens-5.0.262.tab3.txt')
    # A member of another release is not the pinned one.
    with pytest.raises(ValueError, match='0 members named'):
        read(TWO_HYBRID, member='BIOGRID-ALL-5.0.263.tab3.txt')
    two = archive(table(TWO_HYBRID))
    two.result['BIOGRID-ALL-5.0.262.tab3.txt.bak'] = io.StringIO('')
    two.result['BIOGRID-ALL-5.0.262.tab3.txt'] = io.StringIO(table(TWO_HYBRID))
    parser = biogrid.resource.datasets()['physical']._raw_parser
    assert len(list(parser(two))) == 1
    with pytest.raises(ValueError, match='received an HTML document'):
        read(['<head>'], header=['<!DOCTYPE html>'])
    # The obsolete Tab 2.0 layout has other column names.
    with pytest.raises(ValueError, match="missing .*'SWISS-PROT Accessions"):
        read(TWO_HYBRID[:24], header=HEADER[:23] + ['Phenotypes'])
    # The score column is needed.
    with pytest.raises(ValueError, match="missing .*'Score'"):
        read(TWO_HYBRID[:18] + TWO_HYBRID[19:],
             header=HEADER[:18] + HEADER[19:])
    with pytest.raises(ValueError, match='line 2: expected 24 named columns'):
        read(TWO_HYBRID[:20])
    with pytest.raises(ValueError, match='fields beyond the header'):
        read(TWO_HYBRID + ['-'])
    with pytest.raises(ValueError, match="unknown experimental system type 'chemical'"):
        read(changed(TWO_HYBRID, {'Experimental System Type': 'chemical'}))
    assert read() == []


def test_reader_raises_on_a_value_it_does_not_know_in_any_row():
    def bad(column, value):
        return changed(TWO_HYBRID, {column: value})

    # A system of the other type, a system BioGRID may add, an organism that
    # is not a taxon, a new throughput, an identifier of another form and a
    # score that is not a number.
    unknown = [
        (bad('Experimental System', 'Synthetic Lethality'),
         "unknown physical experimental system 'Synthetic Lethality'"),
        (bad('Experimental System', 'Telepathy'),
         "unknown physical experimental system 'Telepathy'"),
        (bad('Organism ID Interactor A', '-'),
         "unknown Organism ID Interactor A '-'"),
        (bad('Organism ID Interactor B', '0'),
         "unknown Organism ID Interactor B '0'"),
        (bad('Throughput', 'Medium Throughput'),
         "unknown Throughput 'Medium Throughput'"),
        (bad('Throughput', 'High Throughput|Medium Throughput'),
         'unknown Throughput'),
        (bad('Throughput', '-'), "unknown Throughput '-'"),
        (bad('BioGRID ID Interactor B', '-'),
         "unknown BioGRID ID Interactor B '-'"),
        (bad('Entrez Gene Interactor A', 'MAP2K4'),
         "unknown Entrez Gene Interactor A 'MAP2K4'"),
        (bad('SWISS-PROT Accessions Interactor A', 'P45985|P4598'),
         'unknown SWISS-PROT Accessions Interactor A'),
        (bad('TREMBL Accessions Interactor B', 'q59h94'),
         "unknown TREMBL Accessions Interactor B 'q59h94'"),
        (bad('REFSEQ Accessions Interactor A', 'NP_003001|P45985'),
         'unknown REFSEQ Accessions Interactor A'),
        (bad('Score', 'high'), "unknown Score 'high'"),
        (bad('Score', '>5'), "unknown Score '>5'"),
        (bad('Score', ''), "unknown Score ''"),
    ]
    for line, message in unknown:
        # Either dataset reads the whole file, so each finds the row, which
        # follows a row of the other type.
        for dataset in ('physical', 'genetic'):
            with pytest.raises(ValueError, match=re.escape(message)):
                read(SYNTHETIC_LETHALITY, line, dataset=dataset)
    # The file writes no RefSeq accession as an empty cell in some rows.
    (row,) = read(bad('REFSEQ Accessions Interactor A', ''))
    assert not [i for i in identifiers(biogrid.map_physical_interaction(row).subject)
                if i[0] == Namespace.REFSEQ_PROTEIN]


def test_two_hybrid_row_keeps_identifiers_system_throughput_and_pmid():
    (row,) = read(TWO_HYBRID)
    relation = biogrid.map_physical_interaction(row)
    assert slot_name(relation.predicate) == 'directly_physically_interacts_with'
    assert identifiers(relation.subject)[:7] == [
        (Namespace.ENTREZ, '6416'),
        (Namespace.BIOGRID, '112315'),
        (Namespace.UNIPROT, 'P45985'),
        (Namespace.REFSEQ_PROTEIN, 'NP_003001'),
        (Namespace.REFSEQ_PROTEIN, 'NP_001268364'),
        (Namespace.GENESYMBOL, 'MAP2K4'),
        (Namespace.SYNONYM, 'JNKK'),
    ]
    assert (Namespace.UNIPROT_TREMBL, 'Q59H94') in identifiers(relation.object)
    assert [(i.type, i.value) for i in relation.identifiers] == [
        ('biogrid_interaction', '103'),
    ]
    # The participant is a gene record: no protein sequence is asserted.
    assert relation.subject.molecular_form is None

    out = extract([row])
    (observed,) = out.relations
    assert observed.predicate == 'directly_physically_interacts_with'
    assert values(observed, 'has_evidence_of_type') == ['MI:0018']
    assert values(observed, 'original_predicate') == ['MI:0407']
    assert values(observed, 'BIOGRID_SYSTEM:0000001') == ['Two-hybrid']
    assert values(observed, 'BIOGRID_TP:0000001') == ['Low Throughput']
    assert values(observed, 'publications') == ['PMID:9006895']
    assert scores(observed) == []
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}


def test_association_evidence_is_not_called_direct():
    capture, fraction = extract(read(BOTH_THROUGHPUTS, CO_FRACTIONATION)).relations
    assert capture.predicate == 'interacts_with'
    assert values(capture, 'has_evidence_of_type') == ['MI:0004']
    assert values(capture, 'original_predicate') == ['MI:0915']
    assert values(capture, 'BIOGRID_SYSTEM:0000001') == ['Affinity Capture-MS']
    assert values(capture, 'BIOGRID_TP:0000001') == [
        'High Throughput', 'Low Throughput',
    ]
    # BioGRID's own PSI-MI type for co-fractionation is colocalization.
    assert fraction.predicate == 'interacts_with'
    assert values(fraction, 'original_predicate') == ['MI:0403']
    assert values(fraction, 'BIOGRID_SYSTEM:0000001') == ['Co-fractionation']


def test_co_fractionation_and_co_localization_stay_interactions_of_the_participants():
    out = extract(read(CO_FRACTIONATION, CO_FRACTIONATION_SCORE, CO_LOCALIZATION))
    assert len(out.relations) == 3
    for relation, system, method in zip(
        out.relations,
        ('Co-fractionation', 'Co-fractionation', 'Co-localization'),
        ('MI:0401', 'MI:0401', 'MI:0428'),
    ):
        assert relation.predicate == 'interacts_with'
        assert values(relation, 'original_predicate') == ['MI:0403']
        assert values(relation, 'has_evidence_of_type') == [method]
        assert values(relation, 'BIOGRID_SYSTEM:0000001') == [system]
    assert {e.entity_type for e in out.entities.values()} == {'protein'}


def test_every_listed_protein_accession_and_name_is_an_identifier():
    relation = biogrid.map_physical_interaction(read(TWO_PROTEINS)[0])
    assert identifiers(relation.subject) == [
        (Namespace.ENTREZ, '140290'),
        (Namespace.BIOGRID, '126602'),
        (Namespace.UNIPROT, 'Q9NV44'),
        (Namespace.UNIPROT, 'Q8TDR4'),
        (Namespace.REFSEQ_PROTEIN, 'NP_653260'),
        (Namespace.GENESYMBOL, 'TCP10L'),
        (Namespace.SYNONYM, 'PRED77'),
        (Namespace.SYNONYM, 'C21orf77'),
        (Namespace.SYNONYM, 'TCP10A-2'),
    ]


def test_protein_rna_systems_keep_the_genes_as_participants():
    out = extract(read(PROTEIN_RNA, VIRAL))
    captured, viral = out.relations
    # The file does not say which participant is the RNA.
    assert {e.entity_type for e in out.entities.values()} == {'gene'}
    assert captured.predicate == 'interacts_with'
    assert values(captured, 'BIOGRID_SYSTEM:0000001') == ['Affinity Capture-RNA']
    assert viral.predicate == 'directly_physically_interacts_with'
    assert values(viral, 'has_evidence_of_type') == ['MI:0686']
    relation = biogrid.map_physical_interaction(read(VIRAL)[0])
    # No Entrez Gene record, no protein accession: the BioGRID record remains.
    assert identifiers(relation.subject) == [
        (Namespace.BIOGRID, '4383869'),
        (Namespace.GENESYMBOL, 'ORF3b'),
        (Namespace.SYNONYM, 'SARS-CoV2 ORF3b'),
        (Namespace.SYNONYM, 'SARS-CoV-2 ORF3b'),
        (Namespace.SYNONYM, '3b'),
    ]
    assert sorted({e.taxon for e in out.entities.values()}) == ['2697049', '9606']


def test_rows_of_other_organisms_keep_the_organism_of_each_participant():
    out = extract(read(YEAST_TWO_HYBRID, FLY_TWO_HYBRID, MOUSE_HUMAN))
    yeast, fly, mouse = out.relations
    assert yeast.predicate == 'directly_physically_interacts_with'
    assert values(yeast, 'publications') == ['PMID:11027275']
    assert values(fly, 'BIOGRID_SYSTEM:0000001') == ['Two-hybrid']
    assert values(mouse, 'original_predicate') == ['MI:0915']
    by_taxon = {}
    for key, entity in out.entities.items():
        by_taxon.setdefault(entity.taxon, []).append(entity)
    assert sorted(by_taxon) == ['10090', '559292', '7227', '9606']
    assert len(by_taxon['559292']) == 2 and len(by_taxon['7227']) == 2
    assert len(by_taxon['10090']) == 1 and len(by_taxon['9606']) == 1
    # Both organisms of a row are kept, each on its participant.
    subject = biogrid.map_physical_interaction(read(MOUSE_HUMAN)[0]).subject
    obj = biogrid.map_physical_interaction(read(MOUSE_HUMAN)[0]).object
    assert [a.value for a in subject.annotations] == ['NCBITaxon:10090']
    assert [a.value for a in obj.annotations] == ['NCBITaxon:9606']
    # Swiss-Prot and TrEMBL accessions of a fly gene without a Swiss-Prot one.
    relation = biogrid.map_physical_interaction(read(FLY_TWO_HYBRID)[0])
    assert (Namespace.UNIPROT_TREMBL, 'Q9W074') in identifiers(relation.subject)
    assert not [i for i in identifiers(relation.subject)
                if i[0] == Namespace.UNIPROT]
    assert (Namespace.UNIPROT, 'P48612') in identifiers(relation.object)


def test_preprint_is_cited_by_its_doi():
    (observed,) = extract(read(PREPRINT)).relations
    assert values(observed, 'publications') == ['doi:10.1101/2021.12.23.473686']
    assert values(observed, 'BIOGRID_TP:0000001') == ['High Throughput']
    # A preprint of a virus, with the same participant on both sides.
    (viral,) = extract(read(VIRAL_PREPRINT)).relations
    assert values(viral, 'publications') == ['doi:10.1101/2020.03.25.996348']
    assert values(viral, 'BIOGRID_SYSTEM:0000001') == ['Co-crystal Structure']
    assert values(viral, 'original_predicate') == ['MI:0914']
    assert viral.predicate == 'interacts_with'


def test_genetic_interactions_are_between_genes_with_their_own_predicate():
    rows = read(*PHYSICAL, *GENETIC, dataset='genetic')
    out = extract(rows, dataset='genetic')
    lethality, rescue = interactions(out)
    assert {r.predicate for r in interactions(out)} == {'genetically_interacts_with'}
    assert {e.entity_type for e in participants(out)} == {'gene'}
    assert values(lethality, 'has_evidence_of_type') == ['MI:0254']
    assert values(lethality, 'original_predicate') == ['MI:2370']
    assert values(lethality, 'BIOGRID_SYSTEM:0000001') == ['Synthetic Lethality']
    assert values(lethality, 'publications') == ['PMID:26354767']
    assert values(rescue, 'original_predicate') == ['MI:2376']
    assert sorted({e.taxon for e in participants(out)}) == ['284812', '9606']
    (lethal,) = interactions(extract(read(DOSAGE_LETHALITY, dataset='genetic'), 'genetic'))
    assert values(lethal, 'original_predicate') == ['MI:2377']
    assert values(lethal, 'BIOGRID_SYSTEM:0000001') == ['Dosage Lethality']


def test_score_is_a_source_attribute_in_each_form_it_takes():
    # Plain, negative and with an exponent; the rows are of two datasets.
    physical = extract(read(
        PLAIN_SCORE, NEGATIVE_SCORE, EXPONENT_SCORE, CO_FRACTIONATION_SCORE,
        BOTH_THROUGHPUTS, TWO_HYBRID,
    ))
    assert [scores(r) for r in physical.relations] == [
        [(0.946, None, 'Score')],
        [(-0.107558193, None, 'Score')],
        [(8.90312e-05, None, 'Score')],
        [(0.874, None, 'Score')],
        [(5.1, None, 'Score')],
        [],
    ]
    genetic = extract(read(
        GENETIC_PLAIN_SCORE, GENETIC_NEGATIVE_SCORE, GENETIC_EXPONENT_SCORE,
        SYNTHETIC_LETHALITY, dataset='genetic',
    ), 'genetic')
    assert [scores(r) for r in interactions(genetic)] == [
        [(2.5219, None, 'Score')],
        [(-5.6431, None, 'Score')],
        [(6.50523e-05, None, 'Score')],
        [],
    ]
    # The scale differs between publications: it is not a confidence.
    for relation in (*physical.relations, *genetic.relations):
        assert 'has_confidence_score' not in {
            a['term'] for a in relation.annotations
        }
    # The annotation keeps the position among the others and is a number.
    relation = biogrid.map_physical_interaction(read(NEGATIVE_SCORE)[0])
    (score,) = [a for a in relation.annotations if a.term == 'biogrid:score']
    assert json.loads(score.value) == {
        'has_numeric_value': -0.107558193, 'has_unit': None,
        'has_unit_prefix': None, 'has_binary_relation': None,
        'source_field': 'Score', 'comparator': None,
    }


def test_no_effect_or_phenotype_is_invented():
    out = extract(read(KINASE_ASSAY, BOTH_THROUGHPUTS))
    kinase, capture = out.relations
    # The modification of a biochemical activity row is not turned into an
    # effect, and the publication's own score is not a confidence value.
    assert kinase.predicate == 'directly_physically_interacts_with'
    for relation in (kinase, capture):
        terms = {a['term'] for a in relation.annotations}
        assert not {t for t in terms if t.endswith('_qualifier')
                    and not t.startswith('biogrid:')}
        assert 'has_confidence_score' not in terms
    genetic = extract(read(*GENETIC, dataset='genetic'), 'genetic')
    for relation in genetic.relations:
        assert 'has_phenotype' not in {a['term'] for a in relation.annotations}


# The 31 experimental systems of the file of release 5.0.262 (every one occurs
# with one system type only), by system type.
OBSERVED_SYSTEMS = {
    'physical': {
        'Affinity Capture-Luminescence', 'Affinity Capture-MS',
        'Affinity Capture-RNA', 'Affinity Capture-Western',
        'Biochemical Activity', 'Co-crystal Structure', 'Co-fractionation',
        'Co-localization', 'Co-purification', 'Cross-Linking-MS (XL-MS)',
        'FRET', 'Far Western', 'PCA', 'Protein-RNA', 'Protein-peptide',
        'Proximity Label-MS', 'Reconstituted Complex', 'Surface Display',
        'Thermal Shift Assay', 'Two-hybrid',
    },
    'genetic': {
        'Dosage Growth Defect', 'Dosage Lethality', 'Dosage Rescue',
        'Negative Genetic', 'Phenotypic Enhancement',
        'Phenotypic Suppression', 'Positive Genetic',
        'Synthetic Growth Defect', 'Synthetic Haploinsufficiency',
        'Synthetic Lethality', 'Synthetic Rescue',
    },
}


def test_every_experimental_system_of_the_release_is_mapped():
    assert {
        kind: set(systems) for kind, systems in biogrid._SYSTEMS.items()
    } == OBSERVED_SYSTEMS
    for kind, template in (('physical', TWO_HYBRID), ('genetic', SYNTHETIC_LETHALITY)):
        for name in sorted(OBSERVED_SYSTEMS[kind]):
            system = biogrid._SYSTEMS[kind][name]
            row = read(changed(template, {
                'Experimental System': name, 'Experimental System Type': kind,
            }), dataset=kind)[0]
            out = extract([row], dataset=kind)
            (observed,) = interactions(out)
            assert values(observed, 'has_evidence_of_type') == [system.detection_method]
            assert values(observed, 'original_predicate') == [system.interaction_type]
            assert values(observed, 'BIOGRID_SYSTEM:0000001') == [name]
            if kind == 'genetic':
                assert observed.predicate == 'genetically_interacts_with'
            elif system.interaction_type == 'MI:0407':
                assert observed.predicate == 'directly_physically_interacts_with'
            else:
                assert observed.predicate == 'interacts_with'
            expected = 'gene' if system.participant.__name__ == 'Gene' else 'protein'
            assert {e.entity_type for e in participants(out)} == {expected}


def test_the_rare_systems_are_mapped_from_rows_of_the_file():
    surface, shift = extract(read(SURFACE_DISPLAY, THERMAL_SHIFT)).relations
    assert values(surface, 'BIOGRID_SYSTEM:0000001') == ['Surface Display']
    assert values(surface, 'has_evidence_of_type') == ['MI:0686']
    assert values(surface, 'original_predicate') == ['MI:0914']
    assert surface.predicate == 'interacts_with'
    assert values(shift, 'BIOGRID_SYSTEM:0000001') == ['Thermal Shift Assay']
    assert values(shift, 'original_predicate') == ['MI:0914']
    assert shift.predicate == 'interacts_with'


def test_unknown_system_or_score_in_the_mapper_raises_and_nothing_is_dropped():
    (row,) = read(TWO_HYBRID)
    unknown = dict(row, **{'Experimental System': 'Telepathy'})
    # A genetic system in a row typed physical is not a known physical system.
    mixed = dict(row, **{'Experimental System': 'Synthetic Lethality'})
    for bad in (unknown, mixed):
        for mapper in (biogrid.map_physical_interaction,
                       biogrid.map_genetic_interaction):
            with pytest.raises(ValueError, match='Unknown BioGRID physical experimental system'):
                mapper(bad)
    score = dict(row, **{'Score': 'high'})
    with pytest.raises(ValueError, match="Unknown BioGRID score 'high'"):
        biogrid.map_physical_interaction(score)
    for kept in (dict(row, Score='-'), dict(row, Score='-2.5e-3')):
        assert biogrid.map_physical_interaction(kept) is not None


def test_publication_without_a_recognised_source_is_not_a_reference():
    (row,) = read(TWO_HYBRID)
    for source in ('PUBMED:unassigned', 'ISBN:0-19-850673-2', 'DOI:n/a', '-'):
        observed, = extract([dict(row, **{'Publication Source': source})]).relations
        assert not values(observed, 'publications')


# Ontology terms and modification ------------------------------------------

KINDS = (
    'biological_process', 'cell_line', 'cell_type', 'cellular_component',
    'chemical', 'development', 'disease', 'environment', 'phenotype', 'tissue',
)
TERM_ATTRIBUTES = {
    f'biogrid:{kind}{suffix}' for kind in KINDS for suffix in ('', '_qualifier')
}
TERM_COLUMNS = [
    'Modification', 'Ontology Term IDs', 'Ontology Term Names',
    'Ontology Term Categories', 'Ontology Term Qualifier IDs',
    'Ontology Term Qualifier Names', 'Ontology Term Types',
]


def mapped(line, *context, dataset='physical'):
    """The mapped relation of a row. The kind of a term can depend on the other
    rows of the file, so ``context`` rows are read with it."""
    rows = read(*context, line, dataset=dataset)
    return getattr(biogrid, f'map_{dataset}_interaction')(rows[-1])


def terms(relation):
    """The annotations of the terms and qualifiers of a mapped relation, by
    attribute, in file order."""
    out = {}
    for a in relation.annotations:
        if a.term in TERM_ATTRIBUTES:
            out.setdefault(a.term, []).append(a.value)
    return out


def annotated(relation, term):
    return [a.value for a in relation.annotations if a.term == term]


def without_terms(line):
    return changed(line, {column: '-' for column in TERM_COLUMNS})


def with_terms(line, ids, categories, qualifiers=None, types=None):
    """The row with other terms (names are not read)."""
    count = len(ids.split('|'))
    qualifiers = qualifiers or '|'.join(['-'] * count)
    types = types or '|'.join(['-'] * count)
    return changed(line, {
        'Ontology Term IDs': ids,
        'Ontology Term Names': '|'.join(['name'] * count),
        'Ontology Term Categories': categories,
        'Ontology Term Qualifier IDs': qualifiers,
        'Ontology Term Qualifier Names': (
            '|'.join('-' if q == '-' else '^'.join(['name'] * len(q.split('^')))
                     for q in qualifiers.split('|'))
        ),
        'Ontology Term Types': types,
    })


def interaction(out):
    (observed,) = interactions(out)
    return observed


def associations(out):
    """The associations of an extraction: (subject key, namespace, identifier)
    of the term it links to, with the type of the term entity checked."""
    found = []
    for r in out.relations:
        if r.predicate == 'associated_with':
            term = out.entities[r.object_entity_key]
            assert term.entity_type == 'ontology_class'
            found.append((r.subject_entity_key, term.namespace, term.identifier))
    return sorted(found)


def linked(out, namespace, *accessions):
    """The associations both participants of the interaction have to the terms."""
    observed = interaction(out)
    return sorted(
        (key, namespace, accession)
        for key in {observed.subject_entity_key, observed.object_entity_key}
        for accession in accessions
    )


def test_terms_are_annotations_by_kind_with_the_accession_as_value():
    expected = [
        (BOTH_THROUGHPUTS, {
            'biogrid:cell_line': ['BTO:0000007'],
            'biogrid:cell_line_qualifier': ['BTO:0000414'],
        }),
        (HEK_CELL_LINE_AND_CELL_TYPE, {
            'biogrid:cell_line': ['BTO:0000007'],
            'biogrid:cell_type': ['BTO:0001025'],
        }),
        (HELA_WITH_DISEASE, {
            'biogrid:cell_line': ['BTO:0000567'],
            'biogrid:cell_line_qualifier': ['DOID:3702'],
        }),
        (PC3_CELL_LINE, {
            'biogrid:cell_line': ['BTO:0001061'],
            'biogrid:cell_line_qualifier': ['DOID:2526'],
        }),
        (PHYSICAL_DISEASE, {
            'biogrid:disease': ['DOID:299'],
            'biogrid:cell_line': ['BTO:0000567'],
        }),
        (CHEMICAL, {'biogrid:chemical': ['CHEBI:16856']}),
        (TISSUE_DEVELOPMENT_AND_CELL_TYPE, {
            'biogrid:tissue': ['UBERON:0000319'],
            'biogrid:development': ['UBERON:0001987'],
            'biogrid:cell_type': ['BTO:0005003'],
            'biogrid:cell_type_qualifier': ['BTO:0001078'],
        }),
        # The terms of a physical interaction labelled phenotype that are
        # phenotypes stay annotations of the interaction; they are not linked.
        (PHYSICAL_APO_PHENOTYPE, {'biogrid:phenotype': ['APO:0000080']}),
        (PHYSICAL_HP_PHENOTYPE, {
            'biogrid:disease': ['DOID:12347'],
            'biogrid:cell_type': ['BTO:0004419'],
            'biogrid:phenotype': ['HP:0003784'],
            'biogrid:biological_process': ['GO:0030199'],
        }),
        (PHYSICAL_MP_PHENOTYPE, {
            'biogrid:tissue': ['UBERON:0002349', 'BTO:0003093'],
            'biogrid:biological_process': ['GO:0050663'],
            'biogrid:phenotype': ['MP:0001625'],
        }),
        # BioGRID's placeholder is carried as written under its category.
        (PLACEHOLDER_CELL_TYPE, {
            'biogrid:cell_type': ['BIOGRID:0000011', 'BTO:0002181'],
        }),
        (PLACEHOLDER_PHENOTYPE, {'biogrid:phenotype': ['BIOGRID:0000011']}),
        (TWO_HYBRID, {}),
        (KINASE_ASSAY, {}),
    ]
    for line, attributes in expected:
        relation = mapped(line)
        assert terms(relation) == attributes
        # Every annotation has a value, and the participants have no links.
        assert all(a.value for a in relation.annotations)
        assert relation.subject.associations is None
        assert relation.object.associations is None
        # The annotations survive the extractor.
        out = extract(read(line))
        observed = interaction(out)
        assert not associations(out)
        assert 'ontology_class' not in {e.entity_type for e in out.entities.values()}
        for attribute, accessions in attributes.items():
            assert values(observed, attribute) == accessions
        # The terms are kept apart from the other attributes of the relation.
        assert {a['term'] for a in observed.annotations} - TERM_ATTRIBUTES >= {
            'has_evidence_of_type', 'BIOGRID_SYSTEM:0000001',
        }


def test_phenotype_labelled_anatomy_and_cells_are_filed_as_what_they_are():
    # Alone, the rows have no label of a cell line: kidney and skin are
    # tissues, so are the cell lines labelled phenotype, and a CL term is a cell type.
    assert terms(mapped(HEK_AS_PHENOTYPE)) == {
        'biogrid:tissue': ['BTO:0000671', 'BTO:0000007'],
        'biogrid:cell_type': ['CL:0002518'],
        'biogrid:tissue_qualifier': ['BTO:0002733'],
    }
    assert terms(mapped(PC3_AS_PHENOTYPE)) == {'biogrid:tissue': ['BTO:0001061']}
    assert terms(mapped(KIDNEY_AS_PHENOTYPE)) == {
        'biogrid:tissue': ['BTO:0000671', 'BTO:0002181'],
        'biogrid:cell_type': ['CL:0002518'],
        'biogrid:tissue_qualifier': ['BTO:0002733'],
    }
    assert terms(mapped(CELL_TYPE_AS_PHENOTYPE)) == {
        'biogrid:cell_type': ['CL:0000738'],
        'biogrid:tissue': ['BTO:0000664', 'BTO:0000089'],
        'biogrid:tissue_qualifier': ['BTO:0002580'],
    }
    # An accession that the file labels "cell line" in any row is a cell line
    # (and its qualifier a qualifier of a cell line), whichever row comes first.
    for context in ([HEK_CELL_LINE_AND_CELL_TYPE], [HEK_CELL_LINE_AND_CELL_TYPE]):
        assert terms(mapped(HEK_AS_PHENOTYPE, *context)) == {
            'biogrid:cell_line': ['BTO:0000007'],
            'biogrid:tissue': ['BTO:0000671'],
            'biogrid:cell_type': ['CL:0002518'],
            'biogrid:cell_line_qualifier': ['BTO:0002733'],
        }
    rows = read(PC3_AS_PHENOTYPE, PC3_CELL_LINE)
    assert terms(biogrid.map_physical_interaction(rows[0])) == {
        'biogrid:cell_line': ['BTO:0001061'],
    }
    # The label of a physical row repairs a term of a genetic row, which is
    # a tissue without it. The phenotypes of that row are linked.
    cell_line = mapped(HELA_AS_PHENOTYPE_GENETIC, HELA_WITH_DISEASE, dataset='genetic')
    assert terms(cell_line) == {'biogrid:cell_line': ['BTO:0000567']}
    assert terms(mapped(HELA_AS_PHENOTYPE_GENETIC, dataset='genetic')) == {
        'biogrid:tissue': ['BTO:0000567'],
    }
    # A cell type of BioGRID's (BTO) stays a cell type, and a cell line label
    # does not make a term of another category a cell line.
    assert terms(mapped(
        with_terms(HELA_WITH_DISEASE, 'BTO:0000567', 'cell type'),
        HELA_WITH_DISEASE,
    )) == {'biogrid:cell_type': ['BTO:0000567']}
    # The anatomy of genetic rows labelled phenotype.
    assert terms(mapped(QUALIFIED_ANATOMY, dataset='genetic')) == {
        'biogrid:tissue': ['ZFA:0000107'],
        'biogrid:tissue_qualifier': ['PATO:0000646'],
    }
    assert terms(mapped(QUALIFIED_TYPED, dataset='genetic')) == {
        'biogrid:tissue': ['ZFA:0001105'],
        'biogrid:tissue_qualifier': ['PATO:0000636'],
    }
    assert terms(mapped(MUSCLE_CELL_AS_PHENOTYPE, dataset='genetic')) == {
        'biogrid:tissue': ['WBbt:0006804'],
        'biogrid:tissue_qualifier': ['WBbt:0006804'],
    }
    for line in (QUALIFIED_ANATOMY, QUALIFIED_TYPED, MUSCLE_CELL_AS_PHENOTYPE):
        out = extract(read(line, dataset='genetic'), 'genetic')
        assert not associations(out)
    out = extract(read(QUALIFIED_ANATOMY, dataset='genetic'), 'genetic')
    assert values(interaction(out), 'biogrid:tissue') == ['ZFA:0000107']


def test_phenotype_terms_of_genetic_interactions_link_both_genes_to_the_term():
    expected = [
        (DOSAGE_LETHALITY, 'apo', ['APO:0000112']),
        (DOSAGE_RESCUE, 'apo', ['APO:0000147', 'APO:0000087']),
        (PARTIAL_RESCUE, 'apo', ['APO:0000111', 'APO:0000049']),
        (SYNTHETIC_LETHALITY, 'hpo', ['HP:0001507']),
        (QUALIFIED_HP_PHENOTYPE, 'hpo', ['HP:0001507']),
        (GO_PHENOTYPE, 'go', ['GO:0008283']),
        (PATO_PHENOTYPE, 'pato', ['PATO:0000169']),
        (MOUSE_PHENOTYPE, 'mp', ['MP:0008762']),
        (WORM_PHENOTYPE, 'wbphenotype', ['WBPhenotype:0000520']),
    ]
    for line, namespace, accessions in expected:
        relation = mapped(line, dataset='genetic')
        out = extract(read(line, dataset='genetic'), 'genetic')
        assert associations(out) == linked(out, namespace, *accessions)
        # Two genes, the terms, no other entity.
        assert len(out.entities) == 2 + len(accessions)
        assert {r.predicate for r in out.relations} == {
            'genetically_interacts_with', 'associated_with',
        }
        assert len(out.relations) == 1 + 2 * len(accessions)
        # The terms are not annotations of the interaction.
        assert 'biogrid:phenotype' not in {a.term for a in relation.annotations}
        for participant in (relation.subject, relation.object):
            assert [
                (a.predicate, a.object.type, str(a.object.identifier_type),
                 a.object.identifier)
                for a in participant.associations
            ] == [
                ('associated_with', 'ontology_class', namespace, accession)
                for accession in accessions
            ]
    # Several vocabularies in one row, with a cell line labelled phenotype
    # (a tissue here) and the type of a partial rescue.
    out = extract(read(CELL_LINE_AS_PHENOTYPE_GENETIC, dataset='genetic'), 'genetic')
    assert associations(out) == sorted(
        linked(out, 'pato', 'PATO:0000169') + linked(out, 'hpo', 'HP:0001507')
    )
    assert values(interaction(out), 'biogrid:tissue') == ['BTO:0000018']
    # A chemical and phenotypes.
    out = extract(read(GENETIC_CHEMICAL, dataset='genetic'), 'genetic')
    assert associations(out) == linked(out, 'apo', 'APO:0000106', 'APO:0000087')
    assert values(interaction(out), 'biogrid:chemical') == ['CHEBI:25255']
    # Terms of other categories are not linked: a process, a tissue.
    out = extract(read(TYPED_PROCESS, dataset='genetic'), 'genetic')
    assert not associations(out)
    assert values(interaction(out), 'biogrid:biological_process') == ['GO:0019076']
    out = extract(read(GENETIC_TISSUE, dataset='genetic'), 'genetic')
    assert associations(out) == linked(out, 'mp', 'MP:0000180')
    assert values(interaction(out), 'biogrid:tissue') == ['BTO:0000131']
    # A cell type of CL, a cell line labelled phenotype (BTO) with a qualifier,
    # a process and skin: none is a phenotype.
    out = extract(
        read(CELL_TYPE_AND_CELL_LINE_AS_PHENOTYPE_GENETIC, dataset='genetic'),
        'genetic',
    )
    assert not associations(out)
    observed = interaction(out)
    assert values(observed, 'biogrid:cell_type') == ['CL:0000148']
    assert values(observed, 'biogrid:tissue') == ['BTO:0002131', 'BTO:0001253']
    assert values(observed, 'biogrid:tissue_qualifier') == ['BTO:0000849']
    assert values(observed, 'biogrid:biological_process') == ['GO:0010468']
    # A repeated accession is linked once.
    doubled = changed(DOSAGE_LETHALITY, {
        'Ontology Term IDs': 'APO:0000112|APO:0000112',
        'Ontology Term Names': 'inviable|inviable',
        'Ontology Term Categories': 'phenotype|phenotype',
        'Ontology Term Qualifier IDs': '-|-',
        'Ontology Term Qualifier Names': '-|-',
        'Ontology Term Types': '-|-',
    })
    out = extract(read(doubled, dataset='genetic'), 'genetic')
    assert associations(out) == linked(out, 'apo', 'APO:0000112')


def test_a_gene_that_interacts_with_itself_is_linked_once():
    out = extract(read(CRISPR_SCREEN, dataset='genetic'), 'genetic')
    observed = interaction(out)
    assert observed.subject_entity_key == observed.object_entity_key
    assert associations(out) == [(observed.subject_entity_key, 'hpo', 'HP:0001507')]
    assert len(out.relations) == 2
    assert values(observed, 'biogrid:cell_type') == ['BTO:0000664']
    # Two different genes of the same row are both linked.
    out = extract(read(SYNTHETIC_LETHALITY, dataset='genetic'), 'genetic')
    observed = interaction(out)
    assert observed.subject_entity_key != observed.object_entity_key
    assert len(associations(out)) == 2


def test_no_association_on_any_physical_row():
    lines = [line for line, dataset in ALL_TERM_ROWS if dataset == 'physical']
    lines += [PHYSICAL_APO_PHENOTYPE, PHYSICAL_HP_PHENOTYPE, PHYSICAL_MP_PHENOTYPE]
    # Also with a genetic phenotype term in the row of the physical dataset.
    lines.append(with_terms(TWO_HYBRID, 'APO:0000112|HP:0001507', 'phenotype|phenotype'))
    for line in lines:
        (row,) = read(line)
        record = biogrid.map_physical_interaction(row)
        assert record.subject.associations is None
        assert record.object.associations is None
        out = extract([row])
        assert not [r for r in out.relations if r.predicate == 'associated_with']
        assert {e.entity_type for e in out.entities.values()} <= {'protein', 'gene'}
    assert terms(mapped(lines[-1])) == {'biogrid:phenotype': ['APO:0000112', 'HP:0001507']}


def test_qualifiers_are_annotations_of_the_kind_of_their_term_and_never_links():
    # A qualifier of a phenotype term of a genetic row is an annotation of the
    # interaction; the term is the link.
    for line, qualifiers_, accession, namespace in (
        (QUALIFIED_HP_PHENOTYPE, ['PATO:0000169'], 'HP:0001507', 'hpo'),
        (SYNTHETIC_LETHALITY, ['HP:0001507'], 'HP:0001507', 'hpo'),
        (QUALIFIED_APO_PHENOTYPE, ['APO:0000003', 'CHEBI:27656'], 'APO:0000087', 'apo'),
    ):
        relation = mapped(line, dataset='genetic')
        assert terms(relation) == {'biogrid:phenotype_qualifier': qualifiers_}
        out = extract(read(line, dataset='genetic'), 'genetic')
        assert values(interaction(out), 'biogrid:phenotype_qualifier') == qualifiers_
        assert associations(out) == linked(out, namespace, accession)
    # The qualifier of a cell line is a qualifier of a cell line.
    assert terms(mapped(HELA_WITH_DISEASE)) == {
        'biogrid:cell_line': ['BTO:0000567'],
        'biogrid:cell_line_qualifier': ['DOID:3702'],
    }
    # Qualifiers of one term, and an identical pair is carried once.
    multiple = changed(HELA_WITH_DISEASE, {
        'Ontology Term Qualifier IDs': 'DOID:3702^BTO:0001421^DOID:3702',
        'Ontology Term Qualifier Names': 'a^b^a',
    })
    assert terms(mapped(multiple)) == {
        'biogrid:cell_line': ['BTO:0000567'],
        'biogrid:cell_line_qualifier': ['DOID:3702', 'BTO:0001421'],
    }
    assert terms(mapped(with_terms(
        HELA_WITH_DISEASE, 'BTO:0000567|BTO:0000567', 'cell line|cell line',
        'DOID:3702|DOID:3702',
    ))) == {
        'biogrid:cell_line': ['BTO:0000567'],
        'biogrid:cell_line_qualifier': ['DOID:3702'],
    }


def test_genetic_phenotypes_keep_their_type_apart():
    lethal = mapped(DOSAGE_LETHALITY, dataset='genetic')
    assert annotated(lethal, 'biogrid:phenotype_type') == []
    rescue = mapped(DOSAGE_RESCUE, dataset='genetic')
    # Both terms are "wild type": the type is one value of the relation.
    assert annotated(rescue, 'biogrid:phenotype_type') == ['wild type']
    partial = mapped(PARTIAL_RESCUE, dataset='genetic')
    assert annotated(partial, 'biogrid:phenotype_type') == ['partial rescue']
    typed = mapped(QUALIFIED_TYPED, dataset='genetic')
    assert annotated(typed, 'biogrid:phenotype_type') == ['wild type']
    # The type of a term that is not a phenotype is also in the file.
    assert annotated(mapped(TYPED_PROCESS, dataset='genetic'), 'biogrid:phenotype_type') == [
        'wild type',
    ]
    # Several types in a row give several values.
    mixed = changed(PARTIAL_RESCUE, {
        'Ontology Term Types': 'partial rescue|undetermined',
    })
    assert annotated(mapped(mixed, dataset='genetic'), 'biogrid:phenotype_type') == [
        'partial rescue', 'undetermined',
    ]
    # Nothing is said of the phenotype beyond the link and the type.
    for relation in (lethal, rescue, partial, typed):
        assert 'has_phenotype' not in {a.term for a in relation.annotations}


# The ontology modules write the identifier of a term as the accession with its
# prefix (``HP:0001507``), as BioGRID does. Stanzas copied from the ontology files.
HP_OBO = """format-version: 1.2
ontology: hp

[Term]
id: HP:0001507
name: Growth abnormality
synonym: "Abnormal growth" EXACT []
is_a: HP:0000118 ! Phenotypic abnormality
"""
GO_OBO = """format-version: 1.2
ontology: go

[Term]
id: GO:0008283
name: cell population proliferation
namespace: biological_process
synonym: "cell proliferation" EXACT []
is_a: GO:0009987 ! cellular process
"""
APO_OBO = """format-version: 1.2
ontology: apo

[Term]
id: APO:0000112
name: inviable
is_a: APO:0000001 ! phenotype
"""
BTO_OBO = """format-version: 1.2
ontology: bto

[Term]
id: BTO:0000007
name: HEK-293 cell
synonym: "293 cell" RELATED []
is_a: BTO:0000067 ! kidney cell line
"""


def ontology_term(module, text, accession):
    """The entity the ontology module makes of the term, through the extractor."""
    (record,) = [r for r in parse_obo_text(text) if r['id'] == accession]
    out = SilverExtractor(module.__name__.rpartition('.')[2], 'terms')
    out.process_record(module.terms_schema(record), record, '1', 1)
    (key,) = [k for k, e in out.entities.items() if e.identifier == accession]
    return key, out.entities[key]


def test_a_linked_phenotype_is_the_entity_of_the_ontology_module():
    for module, text, line, accession in (
        (hpo, HP_OBO, SYNTHETIC_LETHALITY, 'HP:0001507'),
        (go, GO_OBO, GO_PHENOTYPE, 'GO:0008283'),
        (apo, APO_OBO, DOSAGE_LETHALITY, 'APO:0000112'),
    ):
        term_key, term = ontology_term(module, text, accession)
        out = extract(read(line, dataset='genetic'), 'genetic')
        links = [r for r in out.relations if r.predicate == 'associated_with']
        assert len(links) == 2
        for link in links:
            target = out.entities[link.object_entity_key]
            # The same entity type, namespace and accession (prefix kept).
            assert (target.entity_type, target.namespace, target.identifier) == (
                term.entity_type, term.namespace, accession,
            )
            # The key of a reference is the key of that identity. The key of
            # the ontology's own record also hashes the names and synonyms the
            # record carries (as for UniProt's GO links and HPO's gene
            # associations), so the two are the same entity by identity only.
            assert link.object_entity_key == entity_key(
                term.entity_type, term.namespace, accession,
            )
            assert link.object_entity_key != term_key
    # The accession of a term cited as an annotation is the identifier of the
    # ontology module's term.
    _, bto_term = ontology_term(bto, BTO_OBO, 'BTO:0000007')
    relation = mapped(HEK_CELL_LINE_AND_CELL_TYPE)
    assert annotated(relation, 'biogrid:cell_line') == [bto_term.identifier]
    assert bto_term.namespace == 'bto'


# The categories and vocabularies of the file of release 5.0.262 and the kind
# of each (the physical and the genetic rows together).
OBSERVED_KINDS = {
    ('biological process', 'GO'): 'biological_process',
    ('cell line', 'BTO'): 'cell_line',
    ('cell type', 'BIOGRID'): 'cell_type',
    ('cell type', 'BTO'): 'cell_type',
    ('cell type', 'CL'): 'cell_type',
    ('cellular component', 'GO'): 'cellular_component',
    ('chemical', 'CHEBI'): 'chemical',
    ('development', 'UBERON'): 'development',
    ('disease', 'DOID'): 'disease',
    ('environment', 'PATO'): 'environment',
    ('phenotype', 'APO'): 'phenotype',
    ('phenotype', 'BIOGRID'): 'phenotype',
    ('phenotype', 'BTO'): 'tissue',
    ('phenotype', 'CL'): 'cell_type',
    ('phenotype', 'GO'): 'phenotype',
    ('phenotype', 'HP'): 'phenotype',
    ('phenotype', 'MP'): 'phenotype',
    ('phenotype', 'PATO'): 'phenotype',
    ('phenotype', 'WBPhenotype'): 'phenotype',
    ('phenotype', 'WBbt'): 'tissue',
    ('phenotype', 'ZFA'): 'tissue',
    ('tissue', 'BTO'): 'tissue',
    ('tissue', 'EMAP'): 'tissue',
    ('tissue', 'UBERON'): 'tissue',
}


def test_every_category_and_vocabulary_of_the_release_has_a_kind():
    assert {c: biogrid._KINDS[c] for c in OBSERVED_KINDS} == OBSERVED_KINDS
    # The anatomy and cell vocabularies under "phenotype" the owner named.
    assert biogrid._KINDS[('phenotype', 'UBERON')] == 'tissue'
    assert biogrid._KINDS[('phenotype', 'EMAP')] == 'tissue'
    assert set(biogrid._KINDS.values()) == set(KINDS)
    assert set(biogrid._PHENOTYPE_VOCABULARIES) == {
        'APO', 'GO', 'HP', 'MP', 'PATO', 'WBPhenotype',
    }
    assert {
        prefix for (category, prefix), kind in biogrid._KINDS.items()
        if kind == 'phenotype' and prefix != 'BIOGRID'
    } == set(biogrid._PHENOTYPE_VOCABULARIES)
    for (category, prefix), kind in OBSERVED_KINDS.items():
        accession = f'{prefix}:0000001'
        physical = mapped(with_terms(TWO_HYBRID, accession, category))
        assert terms(physical) == {f'biogrid:{kind}': [accession]}
        genetic_line = with_terms(SYNTHETIC_LETHALITY, accession, category)
        if kind == 'phenotype' and prefix == 'BIOGRID':
            # BioGRID's placeholder is no phenotype vocabulary to link.
            with pytest.raises(ValueError, match="Unknown BioGRID phenotype vocabulary 'BIOGRID'"):
                read(genetic_line, dataset='genetic')
        elif kind == 'phenotype':
            genetic = mapped(genetic_line, dataset='genetic')
            assert terms(genetic) == {}
            assert [a.object.identifier for a in genetic.subject.associations] == [accession]
        else:
            assert terms(mapped(genetic_line, dataset='genetic')) == {
                f'biogrid:{kind}': [accession],
            }


def test_modification_is_carried_as_written_and_derives_nothing():
    kinase = mapped(KINASE_ASSAY)
    assert annotated(kinase, 'biogrid:modification') == ['Phosphorylation']
    none = mapped(NO_MODIFICATION)
    # "No Modification" is a statement of the source.
    assert annotated(none, 'biogrid:modification') == ['No Modification']
    for relation in (kinase, none):
        assert relation.predicate == 'directly_physically_interacts_with'
        assert not [a for a in relation.annotations
                    if a.term in {'object_aspect_qualifier',
                                  'object_direction_qualifier'}]
    assert annotated(mapped(TWO_HYBRID), 'biogrid:modification') == []
    # They are kept by the extractor, with the type of the phenotype.
    (observed,) = extract(read(KINASE_ASSAY)).relations
    assert values(observed, 'biogrid:modification') == ['Phosphorylation']
    (observed,) = interactions(extract(read(PARTIAL_RESCUE, dataset='genetic'), 'genetic'))
    assert values(observed, 'biogrid:phenotype_type') == ['partial rescue']
    for name in sorted(biogrid._MODIFICATIONS):
        line = changed(NO_MODIFICATION, {'Modification': name})
        assert annotated(mapped(line), 'biogrid:modification') == [name]


ALL_TERM_ROWS = [
    (BOTH_THROUGHPUTS, 'physical'), (HEK_CELL_LINE_AND_CELL_TYPE, 'physical'),
    (HELA_WITH_DISEASE, 'physical'), (HEK_AS_PHENOTYPE, 'physical'),
    (PC3_AS_PHENOTYPE, 'physical'), (CHEMICAL, 'physical'),
    (KINASE_ASSAY, 'physical'), (NO_MODIFICATION, 'physical'),
    (PHYSICAL_DISEASE, 'physical'), (PHYSICAL_APO_PHENOTYPE, 'physical'),
    (PHYSICAL_HP_PHENOTYPE, 'physical'), (PHYSICAL_MP_PHENOTYPE, 'physical'),
    (TISSUE_DEVELOPMENT_AND_CELL_TYPE, 'physical'),
    (KIDNEY_AS_PHENOTYPE, 'physical'), (CELL_TYPE_AS_PHENOTYPE, 'physical'),
    (PLACEHOLDER_CELL_TYPE, 'physical'), (PLACEHOLDER_PHENOTYPE, 'physical'),
    (SYNTHETIC_LETHALITY, 'genetic'), (DOSAGE_RESCUE, 'genetic'),
    (DOSAGE_LETHALITY, 'genetic'), (PARTIAL_RESCUE, 'genetic'),
    (QUALIFIED_ANATOMY, 'genetic'), (QUALIFIED_TYPED, 'genetic'),
    (CRISPR_SCREEN, 'genetic'), (TYPED_PROCESS, 'genetic'),
    (GO_PHENOTYPE, 'genetic'), (PATO_PHENOTYPE, 'genetic'),
    (MOUSE_PHENOTYPE, 'genetic'), (WORM_PHENOTYPE, 'genetic'),
    (CELL_LINE_AS_PHENOTYPE_GENETIC, 'genetic'), (GENETIC_CHEMICAL, 'genetic'),
    (QUALIFIED_APO_PHENOTYPE, 'genetic'), (QUALIFIED_HP_PHENOTYPE, 'genetic'),
    (CELL_TYPE_AND_CELL_LINE_AS_PHENOTYPE_GENETIC, 'genetic'),
]


def test_the_terms_do_not_change_the_identity_of_the_relation():
    for line, dataset in ALL_TERM_ROWS:
        row = read(line, dataset=dataset)[0]
        blank = read(without_terms(line), dataset=dataset)[0]
        assert row['Ontology Terms'] or row['Modification'] != '-'
        assert blank['Ontology Terms'] == []
        out, out_blank = (
            extract([r], dataset) for r in (row, blank)
        )
        observed, plain = interaction(out), interaction(out_blank)
        assert observed.relation_key == plain.relation_key
        assert observed.subject_entity_key == plain.subject_entity_key
        assert observed.object_entity_key == plain.object_entity_key
        # The same two participants; the only other entities are the terms
        # that the participants are linked to.
        assert {
            k for k, e in out.entities.items() if e.entity_type != 'ontology_class'
        } == out_blank.entities.keys()
        # The relation key hashes the annotations that are Biolink qualifiers.
        # Computed here from the mapped record with its new annotations it is the
        # same: none of them is a qualifier.
        mapper = getattr(biogrid, f'map_{dataset}_interaction')
        record, record_blank = mapper(row), mapper(blank)
        assert (
            len(record.annotations) > len(record_blank.annotations)
            or record.subject.associations
        )
        full = [{'term': a.term, 'value': a.value or '', 'scope': 'relation'}
                for a in record.annotations]
        assert qualifiers(full) == ()
        subject, obj = (observed.subject_entity_key, observed.object_entity_key)
        assert relation_key(subject, observed.predicate, obj, full) == (
            plain.relation_key
        )
        # Every annotation of a term of the extracted relation has a value
        # (a score is a quantity and has none as text).
        assert all(
            a['value'] for a in observed.annotations if a['term'] in TERM_ATTRIBUTES
        )


def test_term_lists_of_unequal_length_or_malformed_accessions_raise():
    def bad(line, **columns):
        return changed(line, columns)

    for dataset in ('physical', 'genetic'):
        # The columns are read in every row of the file, so either dataset
        # finds the row, whichever type it is.
        cases = [
            (bad(SYNTHETIC_LETHALITY, **{'Ontology Term Names': 'one|two'}),
             'unequal length'),
            (bad(HEK_CELL_LINE_AND_CELL_TYPE,
                 **{'Ontology Term Categories': 'cell line'}),
             'unequal length'),
            (bad(HEK_CELL_LINE_AND_CELL_TYPE,
                 **{'Ontology Term Types': '-'}),
             'unequal length'),
            # No accession but the other columns say there is a term.
            (bad(SYNTHETIC_LETHALITY, **{'Ontology Term IDs': '-'}),
             "Unknown ontology term accession '-'"),
            (bad(HELA_WITH_DISEASE, **{'Ontology Term IDs': 'HeLa'}),
             "Unknown ontology term accession 'HeLa'"),
            (bad(HELA_WITH_DISEASE, **{'Ontology Term IDs': 'BTO:'}),
             "Unknown ontology term accession 'BTO:'"),
            (bad(HELA_WITH_DISEASE, **{'Ontology Term IDs': 'BTO:00X'}),
             "Unknown ontology term accession 'BTO:00X'"),
            (bad(HELA_WITH_DISEASE, **{'Ontology Term IDs': ':0000567'}),
             "Unknown ontology term accession ':0000567'"),
            (bad(HELA_WITH_DISEASE, **{'Ontology Term Qualifier IDs': 'cervical'}),
             "Unknown ontology term qualifier 'cervical'"),
            (bad(HELA_WITH_DISEASE,
                 **{'Ontology Term Qualifier IDs': 'DOID:3702^cervix'}),
             'do not match qualifier identifiers|Unknown ontology term qualifier'),
            (bad(HELA_WITH_DISEASE,
                 **{'Ontology Term Qualifier Names': 'a^b'}),
             'do not match qualifier identifiers'),
            (bad(HELA_WITH_DISEASE, **{'Ontology Term Qualifier Names': '-'}),
             'do not match qualifier identifiers'),
            (bad(PARTIAL_RESCUE,
                 **{'Ontology Term Types': 'partial rescue|full rescue'}),
             "Unknown ontology term type 'full rescue'"),
            (bad(PARTIAL_RESCUE,
                 **{'Ontology Term Types': 'partial rescue|Wild type'}),
             "Unknown ontology term type 'Wild type'"),
            (bad(KINASE_ASSAY, Modification='Hydroxylation'),
             "unknown Modification 'Hydroxylation'"),
            (bad(KINASE_ASSAY, Modification='phosphorylation'),
             "unknown Modification 'phosphorylation'"),
            (bad(KINASE_ASSAY, Modification=''), "unknown Modification ''"),
            # A category or a vocabulary that the table of kinds lacks.
            (with_terms(HELA_WITH_DISEASE, 'BTO:0000567', 'sample'),
             "Unknown BioGRID ontology term category 'sample' for the vocabulary 'BTO'"),
            (with_terms(HELA_WITH_DISEASE, 'WBPhenotype:0000001', 'cell line'),
             "category 'cell line' for the vocabulary 'WBPhenotype'"),
            (with_terms(HELA_WITH_DISEASE, 'NCIT:0000001', 'disease'),
             "category 'disease' for the vocabulary 'NCIT'"),
            (with_terms(HELA_WITH_DISEASE, 'DOID:0000001|SNOMED:0000001', 'disease|disease'),
             "vocabulary 'SNOMED'"),
            # A phenotype vocabulary of a genetic row that has no namespace.
            (with_terms(SYNTHETIC_LETHALITY, 'BIOGRID:0000011', 'phenotype'),
             "Unknown BioGRID phenotype vocabulary 'BIOGRID'"),
        ]
        for line, message in cases:
            with pytest.raises(ValueError, match=message):
                read(line, dataset=dataset)
    # A header without the columns of the terms is another layout.
    with pytest.raises(ValueError, match="lacks the columns .*'Ontology Term IDs'"):
        header = [c for c in HEADER if c != 'Ontology Term IDs']
        read([v for v, c in zip(TWO_HYBRID, HEADER) if c != 'Ontology Term IDs'],
             header=header)
    # An ordinary accession of a vocabulary is read, under its category.
    for accession, category in (
        ('BTO:0000007', 'cell line'), ('WBPhenotype:0000001', 'phenotype'),
        ('UBERON:0000955', 'tissue'),
    ):
        (row,) = read(with_terms(HELA_WITH_DISEASE, accession, category))
        assert row['Ontology Terms'][0]['accession'] == accession


def test_the_cell_line_table_needs_a_second_reading_of_the_file():
    class Stream(io.StringIO):
        def seekable(self):
            return False

    parser = biogrid.resource.datasets()['physical']._raw_parser
    stream = SimpleNamespace(result={MEMBER: Stream(table(HEK_AS_PHENOTYPE))})
    with pytest.raises(ValueError, match='cannot be read a second time'):
        list(parser(stream))
    # The accessions that are cell lines are those of the whole file.
    rows = read(HEK_CELL_LINE_AND_CELL_TYPE, HELA_WITH_DISEASE, KINASE_ASSAY)
    assert cell_line_accessions(rows) == {'BTO:0000007', 'BTO:0000567'}
    assert cell_line_accessions([]) == frozenset()


def test_the_parser_decodes_the_terms_of_a_row():
    (row,) = read(HEK_AS_PHENOTYPE)
    assert row['Ontology Terms'] == [
        {'accession': 'BTO:0000671', 'category': 'phenotype',
         'qualifiers': [], 'type': None},
        {'accession': 'CL:0002518', 'category': 'cell type',
         'qualifiers': [], 'type': None},
        {'accession': 'BTO:0000007', 'category': 'phenotype',
         'qualifiers': ['BTO:0002733'], 'type': None},
    ]
    # The terms by kind, with BioGRID's category repaired by the table.
    assert row['Ontology Terms By Kind'] == {
        'tissue': ['BTO:0000671', 'BTO:0000007'], 'cell_type': ['CL:0002518'],
    }
    assert row['Ontology Qualifiers By Kind'] == {'tissue': ['BTO:0002733']}
    assert row['Phenotype Terms By Vocabulary'] == {}
    (row,) = read(QUALIFIED_TYPED, dataset='genetic')
    assert row['Ontology Terms'] == [
        {'accession': 'ZFA:0001105', 'category': 'phenotype',
         'qualifiers': ['PATO:0000636'], 'type': 'wild type'},
    ]
    assert row['Ontology Terms By Kind'] == {'tissue': ['ZFA:0001105']}
    assert row['Ontology Qualifiers By Kind'] == {'tissue': ['PATO:0000636']}
    # The phenotypes of a genetic row are grouped by vocabulary, in file order.
    (row,) = read(CELL_LINE_AS_PHENOTYPE_GENETIC, dataset='genetic')
    assert row['Phenotype Terms By Vocabulary'] == {
        'PATO': ['PATO:0000169'], 'HP': ['HP:0001507'],
    }
    assert row['Ontology Terms By Kind'] == {'tissue': ['BTO:0000018']}
    # The same terms in a physical row are annotations only.
    (row,) = read(with_terms(TWO_HYBRID, 'HP:0001507|PATO:0000169', 'phenotype|phenotype'))
    assert row['Phenotype Terms By Vocabulary'] == {}
    assert row['Ontology Terms By Kind'] == {
        'phenotype': ['HP:0001507', 'PATO:0000169'],
    }
    # Several qualifiers of a term are joined by ^.
    multiple = changed(HELA_WITH_DISEASE, {
        'Ontology Term Qualifier IDs': 'DOID:3702^BTO:0001421',
        'Ontology Term Qualifier Names': 'cervical adenocarcinoma^uterine cervix',
    })
    (row,) = read(multiple)
    assert row['Ontology Terms'][0]['qualifiers'] == ['DOID:3702', 'BTO:0001421']
    assert row['Ontology Qualifiers By Kind'] == {
        'cell_line': ['DOID:3702', 'BTO:0001421'],
    }
    # No term is an empty list, also when the row is of another dataset.
    (row,) = read(TWO_HYBRID)
    assert row['Ontology Terms'] == []
    assert row['Ontology Terms By Kind'] == {}
    # Without a table of kinds the terms are decoded and not grouped.
    rows = list(iter_interactions(
        archive(table(HEK_AS_PHENOTYPE)), system_type='physical',
    ))
    assert rows[0]['Ontology Terms By Kind'] == {}
    assert len(rows[0]['Ontology Terms']) == 3


def test_the_mapper_raises_on_a_type_or_modification_it_does_not_know():
    row = read(PARTIAL_RESCUE, dataset='genetic')[0]
    bad = dict(row, **{'Ontology Terms': [dict(
        row['Ontology Terms'][0], type='full rescue',
    )]})
    with pytest.raises(ValueError, match="Unknown BioGRID ontology term type 'full rescue'"):
        biogrid.map_genetic_interaction(bad)
    (kinase,) = read(KINASE_ASSAY)
    for value in ('Hydroxylation', 'phosphorylation', '', None):
        with pytest.raises(ValueError, match='Unknown BioGRID modification'):
            biogrid.map_physical_interaction(dict(kinase, Modification=value))
    assert biogrid.map_physical_interaction(dict(kinase, Modification='-'))
