"""
Map iPTMnet enzyme-substrate-site relations to Relation records.

iPTMnet's ``ptm.txt`` lists post-translational modification sites, one row
per PTM type, source database, substrate, site and enzyme. The rows that name
an enzyme are the enzyme-substrate relations; each becomes one record, with
the confidence score of ``score.txt`` where that table has the same
substrate, site, enzyme and PTM type.

The enzyme is identified by what the row gives: a UniProt accession (also
when written as the Protein Ontology term of that entry, ``PR:Q92793``), a
Protein Ontology term, a yeast systematic name or, for the rows that give
only a label ("PKA", "CDK1"), that label. A label is never translated to an
accession. A label that is the approved gene symbol of one human gene
(``GENE_SYMBOLS``, ``GENE_SYMBOL_SPELLINGS``) gives a protein with that gene
symbol and the label as a name, and the human taxon where the row is a UniProt
annotation of a human substrate (see ``_enzyme_taxon``); the name of a kinase
family or complex (``GROUP_NAMES``) gives a protein family; any other label
stays a name. The row states the organism of the substrate only, so the
enzyme gets no taxon from it except in that one case.
"""

from __future__ import annotations

import functools
import re
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.iptmnet import (
    PRO_PROTEIN_REGEX,
    UNIPROT,
    UNIPROT_REGEX,
    enzyme_kind,
    iter_enzyme_sites,
)
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

_BASE_URL = 'https://research.bioinformatics.udel.edu/iptmnet_data/files/current/'

_ASPECT = model.GeneOrGeneProductOrChemicalEntityAspectEnum
# The PTM types of ``ptm.txt`` with the aspect of the substrate that the
# enzyme increases and the PSI-MOD name of the modified residue. PSI-MOD
# names ubiquitination and sumoylation for lysine only.
_PTM_TYPES = {
    'PHOSPHORYLATION': (_ASPECT.phosphorylation, 'phosphorylated residue'),
    'ACETYLATION': (_ASPECT.acetylation, 'acetylated residue'),
    'METHYLATION': (_ASPECT.methylation, 'methylated residue'),
    'UBIQUITINATION': (_ASPECT.ubiquitination, None),
    'SUMOYLATION': (_ASPECT.sumoylation, None),
    'O-GLYCOSYLATION': (_ASPECT.o_linked_glycosylation, 'O-glycosylated residue'),
    'N-GLYCOSYLATION': (_ASPECT.n_linked_glycosylation, 'N-glycosylated residue'),
    'MYRISTOYLATION': (_ASPECT.myristoylation, 'myristoylated residue'),
    'FARNESYLATION': (_ASPECT.farnesylation, 'farnesylated residue'),
    'CARBOXYLATION': (_ASPECT.carboxylation, 'carboxylated residue'),
}
_LYSINE_TERMS = {
    'UBIQUITINATION': 'ubiquitinylated lysine',
    'SUMOYLATION': 'sumoylated lysine',
}
# PTM types of ``ptm.txt`` that the aspect values of the pinned Biolink model
# do not name; the nearest values (glycosylation, nitrosation, hydroxylation,
# phosphorylation) are broader or another term. No relation is made of them.
PTM_TYPES_WITHOUT_ASPECT = frozenset(
    {
        'C-GLYCOSYLATION',
        'S-GLYCOSYLATION',
        'S-NITROSYLATION',
        'DIHYDROXYLATION',
        'N-PHOSPHORYLATION',
    }
)
# The values of ``source``, kept as the source writes them: UniProt,
# PhosphAT, HPRD, SGD, GlyGen, PomBase, PhosphoGRID, P3DB, IEDB, dbSNO,
# SIGNOR, neXtProt, the Protein Ontology and IntAct; ``rlim`` and ``rlim+``
# label results of the RLIMS-P text mining tool (the download does not say
# how the two differ).
SOURCE_DATABASES = frozenset(
    {
        'unip',
        'phat',
        'hprd',
        'sgd',
        'glygen',
        'pomb',
        'pgrd',
        'p3db',
        'rlim',
        'rlim+',
        'iedb',
        'sno',
        'sign',
        'npro',
        'pro',
        'intact',
    }
)
# ``organism`` is the substrate's organism, written as UniProt names it (a few
# as "<name> OX=<taxon>"). These are the names of the rows that name an
# enzyme, with the taxon of their substrates' UniProt entries (checked in
# UniProt, October 2026). A substrate of an organism not listed here gets no
# taxon.
ORGANISM_TAXA = {
    'Anaeromyxobacter sp. (strain Fw109-5)': '404589',
    'Anopheles gambiae (African malaria mosquito)': '7165',
    'Apis mellifera (Honeybee)': '7460',
    'Aplysia californica (California sea hare)': '6500',
    'Aplysia californica OX=6500': '6500',
    'Arabidopsis thaliana (Mouse-ear cress)': '3702',
    'Archaeoglobus fulgidus (strain ATCC 49558 / DSM 4304 / JCM 9628 / NBRC 100126 / VC-16)': '224325',
    'Artemia salina OX=85549': '85549',
    'Bacillus subtilis (strain 168)': '224308',
    'Bordetella bronchiseptica (strain ATCC BAA-588 / NCTC 13252 / RB50) (Alcaligenes bronchisepticus)': '257310',
    'Bos taurus (Bovine)': '9913',
    'Brassica napus (Rape)': '3708',
    'Caenorhabditis elegans': '6239',
    'Candida albicans (strain SC5314 / ATCC MYA-2876) (Yeast)': '237561',
    'Canis lupus familiaris (Dog) (Canis familiaris)': '9615',
    'Caulobacter vibrioides (strain ATCC 19089 / CB15) (Caulobacter crescentus)': '190650',
    'Chironomus tentans OX=7153': '7153',
    'Columba livia (Rock dove)': '8932',
    'Coturnix japonica (Japanese quail) (Coturnix coturnix japonica)': '93934',
    'Cricetulus griseus (Chinese hamster) (Cricetulus barabensis griseus)': '10029',
    'Cyprinus carpio (Common carp)': '7962',
    'Danio rerio (Zebrafish) (Brachydanio rerio)': '7955',
    'Dictyostelium discoideum (Social amoeba)': '44689',
    'Drosophila melanogaster (Fruit fly)': '7227',
    'Emericella nidulans (Aspergillus nidulans)': '162425',
    'Emericella nidulans (strain FGSC A4 / ATCC 38163 / CBS 112.46 / NRRL 194 / M139) (Aspergillus nidulans)': '227321',
    'Enterobacter cloacae (strain JD6301)': '1399774',
    'Enterococcus casseliflavus (Enterococcus flavescens)': '37734',
    'Enterococcus faecalis (strain ATCC 700802 / V583)': '226185',
    'Escherichia coli (strain K12)': '83333',
    'Escherichia coli O6:H1 (strain CFT073 / ATCC 700928 / UPEC)': '199310',
    'Escherichia phage T7 (Bacteriophage T7)': '10760',
    'Gallus gallus (Chicken)': '9031',
    'Giardia intestinalis (strain ATCC 50803 / WB clone C6) (Giardia lamblia)': '184922',
    'Glycine max (Soybean) (Glycine hispida)': '3847',
    'Glycine max OX=3847': '3847',
    'Helicobacter pylori (strain ATCC 700392 / 26695) (Campylobacter pylori)': '85962',
    'Helix pomatia (Roman snail) (Edible snail)': '6536',
    'Homo sapiens (Human)': '9606',
    'Hordeum vulgare (Barley)': '4513',
    'Human herpesvirus 1 (strain 17) (HHV-1) (Human herpes simplex virus 1)': '10299',
    'Kirsten murine sarcoma virus': '11808',
    'Lacticaseibacillus casei (Lactobacillus casei)': '1582',
    'Lactococcus lactis subsp. lactis (Streptococcus lactis)': '1360',
    'Leishmania mexicana': '5665',
    'Lilium longiflorum (Trumpet lily)': '4690',
    'Lysinibacillus sphaericus (Bacillus sphaericus)': '1421',
    'Macaca mulatta (Rhesus macaque)': '9544',
    'Mesocricetus auratus (Golden hamster)': '10036',
    'Mus musculus (Mouse)': '10090',
    'Mycobacterium bovis (strain ATCC BAA-935 / AF2122/97)': '233413',
    'Mycobacterium leprae (strain TN)': '272631',
    'Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv)': '83332',
    'Mycobacterium tuberculosis (strain CDC 1551 / Oshkosh)': '83331',
    'Mycolicibacterium smegmatis (strain ATCC 700084 / mc(2)155) (Mycobacterium smegmatis)': '246196',
    'Oryctolagus cuniculus (Rabbit)': '9986',
    'Oryza sativa subsp. indica (Rice)': '39946',
    'Oryza sativa subsp. japonica (Rice)': '39947',
    'Oryzias latipes (Japanese rice fish) (Japanese killifish)': '8090',
    'Ovis aries (Sheep)': '9940',
    'Patiria pectinifera OX=7594': '7594',
    'Pediculus humanus subsp. corporis (Body louse)': '121224',
    'Phaseolus vulgaris OX=3885': '3885',
    'Physarum polycephalum (Slime mold)': '5791',
    'Plasmodium berghei (strain Anka)': '5823',
    'Plasmodium falciparum (isolate 3D7)': '36329',
    'Pseudomonas aeruginosa': '287',
    'Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)': '208964',
    'Pseudomonas syringae pv. syringae': '321',
    'Rattus norvegicus (Rat)': '10116',
    'Rattus norvegicus OX=10116': '10116',
    'Rhinella marina (Cane toad) (Bufo marinus)': '8386',
    'Rous sarcoma virus subgroup A (strain Schmidt-Ruppin) (RSV-SR-A)': '269446',
    'Rous sarcoma virus subgroup C (strain Prague) (RSV-Pr-C)': '11888',
    "Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)": '559292',
    'Salmonella typhimurium (strain 14028s / SGSC 2262)': '588858',
    'Salmonella typhimurium (strain LT2 / SGSC1412 / ATCC 700720)': '99287',
    'Schizosaccharomyces pombe (strain 972 / ATCC 24843) (Fission yeast)': '284812',
    'Serratia marcescens': '615',
    'Shewanella oneidensis (strain MR-1)': '211586',
    'Sinapis alba OX=3728': '3728',
    'Solanum tuberosum (Potato)': '4113',
    'Squalus acanthias (Spiny dogfish)': '7797',
    'Staphylococcus aureus': '1280',
    'Staphylococcus aureus (strain COL)': '93062',
    'Staphylococcus aureus (strain N315)': '158879',
    'Staphylococcus aureus (strain NCTC 8325 / PS 47)': '93061',
    'Staphylococcus aureus (strain Newman)': '426430',
    'Staphylococcus carnosus': '1281',
    'Staphylococcus carnosus (strain TM300)': '396513',
    'Streptococcus pneumoniae (strain ATCC BAA-255 / R6)': '171101',
    'Streptomyces coelicolor (strain ATCC BAA-471 / A3(2) / M145)': '100226',
    'Sus scrofa (Pig)': '9823',
    'Synechococcus elongatus (strain ATCC 33912 / PCC 7942 / FACHB-805) (Anacystis nidulans R2)': '1140',
    'Synechocystis sp. (strain PCC 6803 / Kazusa)': '1111708',
    'Tetronarce californica (Pacific electric ray) (Torpedo californica)': '7787',
    'Thermosynechococcus vestitus (strain NIES-2133 / IAM M-273 / BP-1)': '197221',
    'Tribolium castaneum (Red flour beetle)': '7070',
    'Vibrio alginolyticus': '663',
    'Vibrio campbellii (strain ATCC BAA-1116)': '2902295',
    'Xenopus laevis (African clawed frog)': '8355',
    'Xenopus laevis OX=8355': '8355',
    'Zea mays (Maize)': '4577',
}

# Neither identifier system has a member in ``Namespace``. The PRO namespace is
# used for the PRO terms that are not the protein of a UniProt entry.
PRO_NAMESPACE = 'pr'
LOCUS_NAMESPACE = 'ordered_locus_name'
# Systematic names are those of Saccharomyces cerevisiae S288c.
_YEAST_TAX_ID = '559292'
# The ``note`` of a row that names an enzyme: how the experiment was done, as
# HPRD labels it ("in vitro;in vivo" gives both). The notes of the rows that
# name an enzyme are these three values, all of source ``hprd``; any other
# value stops the mapper. Evidence & Conclusion Ontology terms, as in
# ``hprd.py`` and ``phosphosite.py``.
NOTE_EVIDENCE_TYPES = {
    'in vivo': 'ECO:0000178',  # in vivo assay evidence
    'in vitro': 'ECO:0000181',  # in vitro assay evidence
}
# "S635": residue letter and position.
_SITE_REGEX = re.compile(r'^([A-Z])([1-9]\d*)$')

# Rows whose enzyme has no accession, only a label ``enzyme_genename`` (5,577
# rows, 425 labels in the October 2026 download; all but two are of source
# ``unip``: the "by <kinase>" notes of UniProt feature annotations). The
# labels are classified in the three tables below; a label in none of them
# stays a name of a generic ``NamedThing``.
#
# Approved gene symbols of exactly one human gene, as the label is written
# (HGNC complete set, https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt,
# last modified 2026-10-09, downloaded 2026-10-12). Only the approved symbol
# decides: a label that is only a previous or an alias symbol (PCAF, LKB1,
# IKKB, CDC2) is not here and stays a name. Some of these are also
# an alias of another gene (PDK1 of PDPK1, ATR of ANTXR1): the table has the
# approved symbol, and the 4 human rows of PDK1 are the pyruvate
# dehydrogenase kinase (the substrates are PDHA1 and PDHA2). BAK1, GRIK1 and
# GRIK2 are in the table as symbols, but their rows are all of Arabidopsis
# thaliana, where they are plant kinases; no row of another organism than
# human is mapped to a gene (see ``_enzyme_taxon``).
GENE_SYMBOLS = frozenset(
    """
    AAK1 ABL1 ABL2 ACVR1B ALK ALPK1 ATM ATR AURKA AURKB AURKC BAK1 BCKDK
    BLK BRAF BRSK1 BRSK2 BTK CAMK1 CAMK2A CAMK2B CARM1 CASK CDC7 CDK1
    CDK14 CDK15 CDK16 CDK2 CDK3 CDK4 CDK5 CDK6 CDK7 CDK8 CDK9 CHEK1
    CHEK2 CIT CLK1 CLK2 CLOCK CREBBP CSK CSNK1A1 CSNK1D CSNK1E DAPK1
    DAPK2 DCAF1 DMPK DYRK1A DYRK1B DYRK2 DYRK3 EEF2K EGFR EHMT1 EHMT2
    EP300 EPHA2 EPHA4 EPHB2 EZH2 FAM20C FER FES FGFR1 FGFR3 FGR FLT3 FRK
    FYN GALNT13 GRIK1 GRIK2 GRK5 GRK6 GRK7 HASPIN HCK HIPK1 HIPK2 HIPK4
    INSR IRAK1 IRAK4 ITK JAK1 JAK2 KAT2A KAT2B KAT5 KAT6A KAT8 KMT5A
    LATS1 LATS2 LCK LRRK2 LYN MAP2K1 MAP2K2 MAP2K3 MAP2K4 MAP2K6 MAP2K7
    MAP3K10 MAP3K14 MAP3K20 MAP3K5 MAP3K6 MAP3K7 MAP4K1 MAP4K4 MAPK1
    MAPK10 MAPK11 MAPK12 MAPK13 MAPK14 MAPK3 MAPK4 MAPK6 MAPK7 MAPK8
    MAPK9 MAPKAPK2 MAPKAPK3 MAPKAPK5 MARK1 MARK2 MARK3 MARK4 MATK MELK
    MET MINK1 MKNK2 MTOR NAT10 NEDD4 NEK1 NEK11 NEK2 NEK6 NEK7 NEK9 NLK
    NRK NTRK2 NUAK1 NUAK2 OXSR1 PAK1 PAK2 PAK3 PAK4 PAK5 PASK PCK1 PDK1
    PDK2 PDK3 PDK4 PDPK1 PGK1 PIK3CG PIM1 PIM2 PIM3 PINK1 PKDCC PKMYT1
    PLK1 PLK2 PLK3 PLK4 PRKDC PRKX PRMT1 PRMT5 PRMT6 PRMT8 PRMT9 PTK6
    RAF1 RET RIOK3 RIPK1 RIPK3 ROCK1 ROCK2 RPS6KA1 RPS6KA3 RPS6KA4
    RPS6KA5 RPS6KB1 SETD6 SETD7 SGK1 SGK3 SIK1 SIK2 SLK SMYD2 SPHK1 SRC
    SRPK2 STK10 STK16 STK26 STK39 SYK TAF1 TAOK1 TBK1 TEC TGFBR1 TGFBR2
    TLK2 TNIK TNK2 TRPM7 TSSK2 TSSK4 TTBK1 TTBK2 TTK TXK TYK2 UHMK1 ULK1
    VRK1 VRK2 WEE1 WEE2 WNK1 WNK4 ZAP70
    """.split()
)
# Labels that are the approved symbol of exactly one human gene in another
# letter case (same HGNC release), with that symbol. "HPr", the phosphocarrier
# protein of bacteria, matches the symbol of the human haptoglobin-related
# gene HPR in letter case only and is not here.
GENE_SYMBOL_SPELLINGS = {
    'CaMK1': 'CAMK1',
    'CaMK1D': 'CAMK1D',
    'CaMK2A': 'CAMK2A',
    'CaMK2D': 'CAMK2D',
    'CaMK2G': 'CAMK2G',
    'CaMK4': 'CAMK4',
    'CaMKK1': 'CAMKK1',
    'CaMKK2': 'CAMKK2',
    'Cdk1': 'CDK1',
    'Haspin': 'HASPIN',
    'MAPk14': 'MAPK14',
    'Pink1': 'PINK1',
    'Src': 'SRC',
    'atm': 'ATM',
    'atr': 'ATR',
    'aurka': 'AURKA',
    'cdk1': 'CDK1',
    'cdk2': 'CDK2',
    'plk1': 'PLK1',
    'wee1': 'WEE1',
    'wee2': 'WEE2',
}
# Names of kinase families, subfamilies and complexes that several genes
# encode, as the labels write them. They are not gene symbols; each is a
# protein family, identified by the name only and never expanded to members.
# "CAMK2", "CaMK2" and "CamK2" are one family written three ways; the name is
# kept as written.
GROUP_NAMES = frozenset(
    {
        'ABL',
        'AMPK',
        'CAK',
        'CAMK',
        'CAMK2',
        'CDC42BP',
        'CK1',
        'CK2',
        'CaMK',
        'CaMK2',
        'CamK2',
        'GSK3',
        'JAK',
        'MAP2K',
        'MAP3K',
        'MAPK',
        'MLCK',
        'PAK',
        'PDGFR',
        'PHK',
        'PI3K',
        'PKA',
        'PKB',
        'PKC',
        'PKD',
        'PKG',
        'RAF',
        'TORC1',
        'TORC2',
        'cak',
    }
)
# The taxon of the enzyme of a gene symbol label. The row states the
# substrate's organism, not the enzyme's. UniProt's "by <kinase>" note on an
# entry names the kinase of the same organism as the annotated protein, so for
# a row of source ``unip`` with a human substrate the gene symbol label is the
# human gene and the enzyme is human. This is the only case in which the
# label is taken as a gene symbol: HGNC symbols are human symbols, and for a
# substrate of another organism (mouse "Cdk1", yeast "CDC28") or a row of
# another source the label is not known to be the symbol of that gene, so the
# enzyme stays a name without a taxon.
_HUMAN = 'NCBITaxon:9606'
_UNIPROT_FEATURE_SOURCE = 'unip'


def _enzyme_kind(row: dict[str, Any]) -> str | None:
    """The kind of the enzyme of a row: its identifier system, or ``label``.

    ``label`` is a row whose enzyme has no accession and is named only by
    ``enzyme_genename`` ("PKA", "CK2", "CDK1"). None for a value of no known
    form, which the parser has left out.
    """
    value = row.get('enzyme_UniProtAC') or ''
    if not value:
        return 'label' if row.get('enzyme_genename') else None
    return enzyme_kind(value)


def _uniprot_enzyme(row: dict[str, Any]) -> str | None:
    """The accession of a UniProt enzyme, or the one a PRO term of a UniProt entry names."""
    kind = _enzyme_kind(row)
    if kind == 'uniprot':
        return row['enzyme_UniProtAC']
    if kind == 'pro_protein':
        return PRO_PROTEIN_REGEX.match(row['enzyme_UniProtAC'])[1]
    return None


def _enzyme_id(kind: str):
    """``enzyme_UniProtAC`` of the rows whose enzyme identifier is of this kind."""
    return lambda row: row['enzyme_UniProtAC'] if _enzyme_kind(row) == kind else None


def _label_gene_symbol(row: dict[str, Any]) -> str | None:
    """The gene symbol of a label-only row, or None: the label is a name.

    Only the rows of UniProt annotations of a human substrate, whose label is
    an approved human gene symbol (in the letter case written or another),
    have one; the comment above ``_HUMAN`` gives the reason.
    """
    if _enzyme_kind(row) != 'label':
        return None
    if row.get('source') != _UNIPROT_FEATURE_SOURCE or _substrate_taxon(row) != _HUMAN:
        return None
    label = row['enzyme_genename']
    if label in GENE_SYMBOLS:
        return label
    return GENE_SYMBOL_SPELLINGS.get(label)


def _enzyme_type(row: dict[str, Any]) -> Any:
    kind = _enzyme_kind(row)
    if kind == 'label':
        if _label_gene_symbol(row):
            return model.Protein
        if row['enzyme_genename'] in GROUP_NAMES:
            return model.ProteinFamily
        return model.NamedThing
    # A PRO term of another kind does not say what the enzyme is.
    return model.NamedThing if kind == 'pro_term' else model.Protein


def _enzyme_label(row: dict[str, Any]) -> str | None:
    """The enzyme label of a row that has no accession, as written."""
    return row['enzyme_genename'] if _enzyme_kind(row) == 'label' else None


def _enzyme_gene(row: dict[str, Any]) -> str | None:
    # For a row with an accession the gene name is the source's. A label
    # alone does not identify the enzyme unless it is a gene symbol of a
    # human enzyme; see ``_label_gene_symbol``.
    if _enzyme_kind(row) == 'label':
        return _label_gene_symbol(row)
    return (row.get('enzyme_genename') or None) if _enzyme_kind(row) is not None else None


def _enzyme_taxon(row: dict[str, Any]) -> str | None:
    # The table gives the substrate's organism only; the enzyme can be of
    # another species (a viral kinase of a human substrate). The exceptions:
    # a yeast systematic name, and a gene symbol label of a human substrate's
    # UniProt annotation (``_label_gene_symbol``).
    if _enzyme_kind(row) == 'yeast_locus':
        return f'NCBITaxon:{_YEAST_TAX_ID}'
    if _label_gene_symbol(row):
        return _HUMAN
    return None


def _substrate_taxon(row: dict[str, Any]) -> str | None:
    taxon = ORGANISM_TAXA.get(row.get('organism') or '')
    return f'NCBITaxon:{taxon}' if taxon else None


def _evidence_types(row: dict[str, Any]) -> list[str]:
    """ECO terms of the ``note`` experiment types; an unknown type is an error."""
    terms = []
    for note in (row.get('note') or '').split(';'):
        if not note:
            continue
        if note not in NOTE_EVIDENCE_TYPES:
            raise ValueError(f'Unknown iPTMnet note: {note!r}')
        if NOTE_EVIDENCE_TYPES[note] not in terms:
            terms.append(NOTE_EVIDENCE_TYPES[note])
    return terms


def _ptm_type(row: dict[str, Any]) -> str:
    value = row.get('ptm_type')
    if value not in _PTM_TYPES:
        raise ValueError(f'Unknown iPTMnet PTM type: {value!r}')
    return value


def _source_database(row: dict[str, Any]) -> str:
    value = row.get('source')
    if value not in SOURCE_DATABASES:
        raise ValueError(f'Unknown iPTMnet source database: {value!r}')
    return value


def _substrate_form(row: dict[str, Any]) -> dict | None:
    """The modified residue, a feature of the substrate in this record.

    iPTMnet gives the position in the sequence of ``substrate_UniProtAC``,
    which is an isoform accession where the source has one. A site not
    written as a letter and a number keeps its text and gets no position.
    """
    site = row.get('site')
    if not site or not row.get('substrate_UniProtAC'):
        return None
    ptm_type = _ptm_type(row)
    match = _SITE_REGEX.match(site)
    if not match:
        return {'modifications': [{'term': _PTM_TYPES[ptm_type][1], 'description': site}]}
    residue, position = match[1], int(match[2])
    term = _PTM_TYPES[ptm_type][1]
    if term is None and residue == 'K':
        term = _LYSINE_TERMS[ptm_type]
    return {
        'modifications': [
            {
                'term': term,
                'residue': residue,
                'position': position,
                'end_position': position,
                'coordinate_reference': {
                    'identifier': {'ns': 'uniprot', 'id': row['substrate_UniProtAC']},
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
            }
        ],
    }


config = ResourceConfig(
    id=ResourceCv.IPTMNET,
    name='iPTMnet',
    url='https://research.bioinformatics.udel.edu/iptmnet/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_NC_SA_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='29145615',
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.TEXT_MINED,
    ),
    summary=(
        'Enzyme-substrate relations with modification sites, PTM types, '
        'source databases and confidence scores, across species.'
    ),
    description=(
        'iPTMnet is an integrated resource for post-translational '
        'modification (PTM) networks of the Protein Information Resource. '
        'It combines PTM sites and PTM enzyme-substrate-site relations from '
        'curated databases, among them UniProt, HPRD, SIGNOR, neXtProt, '
        'PhosphoGRID and the Protein Ontology, with results of the RLIMS-P '
        'and eFIP text mining tools, for human and many other organisms. '
        'Each relation has its PTM type, site, source and publications, and '
        'a confidence score from 0 to 4 based on the quality and the number '
        'of its sources and publications.'
    ),
)

f = FieldConfig(extract={'uniprot': rf'^({UNIPROT})$'})

enzyme = EntityBuilder(
    entity_type=_enzyme_type,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=_uniprot_enzyme),
        CV(term=PRO_NAMESPACE, value=_enzyme_id('pro_term')),
        CV(term=LOCUS_NAMESPACE, value=_enzyme_id('yeast_locus')),
        CV(term=Namespace.GENESYMBOL, value=_enzyme_gene),
        CV(term=Namespace.NAME, value=_enzyme_label),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=_enzyme_taxon),
    ),
)

substrate = EntityBuilder(
    entity_type=model.Protein,
    # The modified residue is in the substrate's molecular form.
    molecular_form=_substrate_form,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f('substrate_UniProtAC', extract='uniprot')),
        CV(
            term=Namespace.GENESYMBOL,
            value=lambda row: (
                (row.get('substrate_genename') or None)
                if UNIPROT_REGEX.match(row.get('substrate_UniProtAC') or '')
                else None
            ),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=_substrate_taxon),
    ),
)

_enzyme_substrate_schema = RelationBuilder(
    subject=enzyme,
    predicate=slots.affects,
    object=substrate,
    annotations=AnnotationsBuilder(
        # The enzyme increases this modification of the substrate; the table
        # says nothing about the substrate's activity.
        CV(
            term=slots.object_direction_qualifier,
            value=model.DirectionQualifierEnum.increased,
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=lambda row: _PTM_TYPES[_ptm_type(row)][0],
        ),
        CV(
            term=slots.publications,
            value=f(
                'pmid',
                delimiter=',',
                extract=PMID_EXTRACT,
                transform=lambda pmid: f'PMID:{pmid}',
            ),
        ),
        # The database this row of the table comes from.
        CV(term=slots.supporting_data_source, value=_source_database),
        CV(term=slots.has_evidence_of_type, value=_evidence_types),
        # 0 to 4, 4 the highest confidence.
        CV(
            term=slots.has_confidence_score,
            value=lambda row: measurement(row.get('score'), source_field='score'),
        ),
    ),
)


def enzyme_substrate_schema(row: dict[str, Any]) -> Any:
    """Map a row that names an enzyme; None if its PTM type has no Biolink aspect.

    """
    if row.get('ptm_type') in PTM_TYPES_WITHOUT_ASPECT:
        return None
    return _enzyme_substrate_schema(row)


score_download = Download(
    url=_BASE_URL + 'score.txt',
    filename='score.txt',
    subfolder='iptmnet',
    large=True,
)


resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=Download(
            url=_BASE_URL + 'ptm.txt',
            filename='ptm.txt',
            subfolder='iptmnet',
            large=True,
        ),
        mapper=enzyme_substrate_schema,
        raw_parser=functools.partial(
            iter_enzyme_sites, scores_download=score_download
        ),
    ),
)
