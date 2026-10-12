"""Read the rescued dbPTM table of enzyme-named modification sites.

dbPTM integrates modification sites from other databases. Only its
releases up to dbPTM 3.0 named the enzyme of a site, taken from the
"by <enzyme>" notes of Swiss-Prot; the per-type files dbPTM publishes today
(entry name, accession, position, type, PubMed IDs, sequence window) have
no enzyme. The dbPTM 3.0 files are gone from the dbPTM hosts, so the table
is the one OmniPath rescued from them
(``rescued.omnipathdb.org/old_dbptm.tab``): human sites, tab separated, with
a header:

    substrate kinase resnum resaa start end instance      typ             source     references
    P31749    TNK2   176    Y     170   182 KATGRYYAMKILK phosphorylation Swiss-Prot 20333297

The table is a product of the legacy pypath reader, not a file of dbPTM:
``resaa`` is the seventh letter of ``instance``, ``start`` and ``end`` are
``resnum`` -6 and +6, a site modified "by autocatalysis" has the substrate's
accession as its enzyme, and several enzymes of one site are joined by ";".
Most rows have no enzyme. A site without an enzyme is not an
enzyme-substrate relation and is not read.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = (
    'substrate',
    'kinase',
    'resnum',
    'resaa',
    'start',
    'end',
    'instance',
    'typ',
    'source',
    'references',
)
UNIPROT_REGEX = re.compile(
    r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|'
    r'[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$'
)
# The length of a sequence window that has the site in its middle. In the
# shorter windows of the table the seventh letter is not the site.
WINDOW_LENGTH = 13

# The enzyme labels of the table are the names Swiss-Prot wrote after "by".
# Each of the following is the approved HGNC symbol of a protein coding gene
# (checked against the HGNC complete set of October 2026).
GENE_SYMBOLS = frozenset(
    """
    ABL1 ABL2 ACVR1B ATM ATR AURKA AURKB AURKC BTK CAMK1 CARM1 CDC7 CDK1
    CDK14 CDK2 CDK3 CDK4 CDK5 CDK6 CDK7 CDK8 CDK9 CHEK1 CHEK2 CLK1 CLK2
    CREBBP CSK CSNK1D CSNK1E DAPK1 DMPK DYRK1A DYRK2 EGFR EHMT1 EHMT2 EP300
    FER FES FGR FRK FYN GRK5 GRK6 HCK HIPK2 HIPK4 INSR ITK JAK1 JAK2 KAT2A
    KAT2B KAT8 LATS1 LATS2 LCK LYN MAP2K1 MAP2K2 MAP2K3 MAP2K4 MAP2K6 MAP2K7
    MAP3K14 MAP3K5 MAP4K1 MAPK1 MAPK10 MAPK11 MAPK12 MAPK13 MAPK14 MAPK3
    MAPK7 MAPK8 MAPK9 MAPKAPK2 MAPKAPK3 MAPKAPK5 MARK1 MARK2 MELK MET MINK1
    MKNK2 MTOR NEK11 NEK2 NEK6 NEK9 NLK NUAK1 NUAK2 PAK1 PAK2 PAK3 PAK4 PASK
    PDK4 PDPK1 PIM1 PKMYT1 PLK1 PLK2 PLK3 PLK4 PRKDC PRKX PRMT1 PRMT5 PRMT6
    PRMT8 PTK6 RET ROCK1 ROCK2 RPS6KA1 RPS6KA3 RPS6KA4 RPS6KA5 RPS6KB1 SETD7
    SGK1 SGK3 SIK1 SIK2 SMYD2 SRC SRPK2 STK10 STK16 STK24 STK3 STK4 SYK TAF1
    TAOK1 TBK1 TEC TGFBR1 TGFBR2 TLK2 TNK2 TRPM7 TTBK1 TTK TXK UHMK1 VRK1
    VRK2 WEE1 WNK1 WNK4 ZAP70 ZDHHC2 ZDHHC5
    """.split()
)
# Names of kinase families, subfamilies and holoenzymes that several genes
# encode. They are not gene symbols and stay names of a group of proteins.
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
        'CaMK2',
        'GSK3',
        'MAP2K',
        'MAP3K',
        'MAPK',
        'MLCK',
        'PAK',
        'PDGFR',
        'PHK',
        'PKA',
        'PKB',
        'PKC',
        'PKD',
        'PKG',
        'RAF',
        'dsDNA kinase',
        'dsDNAkinase',
    }
)
# Labels that name exactly one gene by its approved HGNC symbol (October 2026)
# in another letter case, with that symbol. The label is also kept as written.
GENE_SYMBOL_SPELLINGS = {
    'CaMK1': 'CAMK1',
    'CaMK1D': 'CAMK1D',
    'CaMK2D': 'CAMK2D',
    'CaMK4': 'CAMK4',
    'CaMKK1': 'CAMKK1',
    'CaMKK2': 'CAMKK2',
}
# Former symbols, aliases and protein names of one enzyme. They are kept as
# names and not replaced by the current gene symbol.
PROTEIN_NAMES = frozenset(
    {
        'BARK1',
        'CCRK',
        'CDC2',
        'FAK1',
        'FAK2',
        'IKKA',
        'IKKB',
        'IKKE',
        'LKB1',
        'MLTK',
        'P300',
        'PAK7',
        'PCAF',
        'SETD8',
        'ZIPK',
        'p300',
    }
)
# An acetyltransferase of the bacterium Yersinia; the label names no species.
BACTERIAL_NAMES = frozenset({'Yersinia yopJ'})
# One cell that joins two enzymes in a way the table's separator missed.
JOINED_LABELS = {'MAPK1 AND MAPK3': ('MAPK1', 'MAPK3')}
# Labels that do not name one enzyme: alternatives between which the source
# does not decide, and the residue of "by Tyr-kinases". No relation is made.
NOT_ONE_ENZYME = frozenset(
    {'ATM or ATR', 'MAPK1 or MAPK3', 'MAPKAPK2 or MAPKAPK3', 'Tyr'}
)
# The modification types of the rows that name an enzyme. The other types of
# the table occur only without an enzyme.
MODIFICATION_TYPES = frozenset(
    {'phosphorylation', 'methylation', 'acetylation', 'palmitoylation'}
)
# The databases dbPTM took its sites from, as the table labels them.
SOURCE_DATABASES = frozenset(
    {
        'HPRD',
        'MeMo',
        'OGlycBase6.0',
        'PhosphoELM.10011',
        'Phosphositeplus.1010730',
        'Swiss-Prot',
        'SysPTM',
        'UbiProtDB',
        'dbSNO',
    }
)


def _enzyme_fields(label: str) -> dict[str, str | None]:
    """The identifier of an enzyme label, in the one field of its kind."""
    fields: dict[str, str | None] = dict.fromkeys(
        ('enzyme_uniprot', 'enzyme_genesymbol', 'enzyme_name', 'enzyme_group')
    )
    if UNIPROT_REGEX.match(label):
        fields['enzyme_uniprot'] = label
    elif label in GENE_SYMBOLS:
        fields['enzyme_genesymbol'] = label
    elif label in GENE_SYMBOL_SPELLINGS:
        fields['enzyme_genesymbol'] = GENE_SYMBOL_SPELLINGS[label]
        fields['enzyme_name'] = label
    elif label in PROTEIN_NAMES or label in BACTERIAL_NAMES:
        fields['enzyme_name'] = label
    elif label in GROUP_NAMES:
        fields['enzyme_group'] = label
    elif label not in NOT_ONE_ENZYME:
        raise ValueError(f'Unknown dbPTM enzyme label: {label!r}')
    return fields


def iter_enzyme_sites(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per enzyme of each site that names an enzyme.

    Rows with an empty ``kinase`` cell are skipped. Each dictionary has the
    table's columns as written, and:

    - ``enzyme``: one label of the ``kinase`` cell.
    - ``enzyme_uniprot``, ``enzyme_genesymbol``, ``enzyme_name``,
      ``enzyme_group``: the label in the field of its kind; all None for a
      label that does not name one enzyme ("MAPK1 or MAPK3"). A label that
      is a gene symbol in another letter case ("CaMK1") has the symbol in
      ``enzyme_genesymbol`` and the label as written in ``enzyme_name``.

    A table without a header column of ``COLUMNS`` or without rows is an error.
    - ``position``: ``resnum`` as an integer.
    - ``residue``: ``resaa`` if ``instance`` is a 13 residue window, which has
      the site in its middle; None for the shorter windows.
    - ``pmids``: the PubMed IDs of ``references``; "-" stands for none.
    """
    for row in iter_tsv(opener, required=COLUMNS):
        if not row['kinase']:
            continue
        if row['typ'] not in MODIFICATION_TYPES:
            raise ValueError(f'Unknown dbPTM modification type with an enzyme: {row["typ"]!r}')
        if row['source'] not in SOURCE_DATABASES:
            raise ValueError(f'Unknown dbPTM source database: {row["source"]!r}')
        row['position'] = int(row['resnum'])
        row['residue'] = (
            row['resaa'] if len(row['instance']) == WINDOW_LENGTH else None
        )
        row['pmids'] = [
            ref
            for ref in map(str.strip, row['references'].split(';'))
            if is_pmid(ref)
        ]
        for cell in row['kinase'].split(';'):
            for label in JOINED_LABELS.get(cell, (cell,)):
                yield {**row, 'enzyme': label, **_enzyme_fields(label)}
