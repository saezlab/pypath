"""Read the phosphorylation data file of MIMP, the sites that name a kinase.

MIMP predicts the effect of mutations on phosphorylation from kinase
specificity models. Its phosphorylation data file is not a prediction: it
holds the experimentally identified human phosphosites MIMP combined from
PhosphoSitePlus, Phospho.ELM and HPRD, each matched to the longest CCDS
isoform of its gene, and for a part of the sites the kinases those databases
name. The specificity models were trained on these kinase annotations.

Source file (tab separated, with a header, one site per line):

    gene   gene_acc   kinase      position  flank            npmid  source
    YWHAB  NM_003404  PDKC;PRKCD  60        VVGARRSSWRVISSI  8      HPRD;PhosphoELM;PhosphoSitePlus

``gene_acc`` is the RefSeq transcript of the isoform, ``flank`` the site
with seven residues on each side, ``npmid`` the number of publications of
the site. Most sites have no kinase; a site without a kinase is not an
enzyme-substrate relation and is not read.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = ('gene', 'gene_acc', 'kinase', 'position', 'flank', 'npmid', 'source')
REFSEQ_REGEX = re.compile(r'^NM_\d+$')
# Seven residues on each side of the site; "-" pads a window at a terminus.
FLANK_REGEX = re.compile(r'^[A-Z-]{7}[A-Z][A-Z-]{7}$')
# A kinase label is the gene symbol of the kinase as it was in 2015, a few of
# them since replaced (ADRBK1, ZAK, PAK7), or one of the following.
GENE_SYMBOL_REGEX = re.compile(r'^[A-Z][A-Z0-9]*$')
# A kinase group of the Phospho.ELM vocabulary, e.g. "PKA_GROUP". The file
# does not say which kinase of the group phosphorylates the site, and the
# group is not replaced by its members.
GROUP_SUFFIX = '_GROUP'
# Kinase names that HGNC has never had as a gene symbol (checked against the
# HGNC complete set of October 2026): a fusion protein and three names of the
# source databases.
KINASE_NAMES = frozenset({'BCR/ABL', 'CCDPK', 'PDKC', 'PRK4'})
# The databases MIMP combined, as the source column labels them.
SOURCE_DATABASES = frozenset({'HPRD', 'PhosphoELM', 'PhosphoSitePlus'})


def _kinase_fields(label: str) -> dict[str, str | None]:
    """The identifier of a kinase label, in the one field of its kind."""
    fields: dict[str, str | None] = dict.fromkeys(
        ('enzyme_genesymbol', 'enzyme_name', 'enzyme_group')
    )
    if label.endswith(GROUP_SUFFIX):
        fields['enzyme_group'] = label
    elif label in KINASE_NAMES:
        fields['enzyme_name'] = label
    elif GENE_SYMBOL_REGEX.match(label):
        fields['enzyme_genesymbol'] = label
    else:
        raise ValueError(f'Unknown MIMP kinase label: {label!r}')
    return fields


def iter_kinase_sites(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per kinase of each phosphosite that names a kinase.

    Rows with an empty ``kinase`` cell are skipped. Each dictionary has the
    file's columns as written, and:

    - ``enzyme``: one label of the ``kinase`` cell.
    - ``enzyme_genesymbol``, ``enzyme_name``, ``enzyme_group``: the label in
      the field of its kind.
    - ``site_position``: ``position`` as an integer.
    - ``residue``: the middle letter of ``flank``.
    - ``sources``: the databases of ``source``.

    A file without a column of ``COLUMNS`` or without rows is an error.
    """
    for row in iter_tsv(opener, required=COLUMNS):
        if not row['kinase']:
            continue
        if (
            not row['gene']
            or not REFSEQ_REGEX.match(row['gene_acc'])
            or not row['position'].isdigit()
            or not FLANK_REGEX.match(row['flank'])
        ):
            raise ValueError(f'Malformed MIMP phosphosite: {row!r}')
        sources = row['source'].split(';')
        unknown = sorted(set(sources) - SOURCE_DATABASES)
        if unknown:
            raise ValueError(f'Unknown MIMP source database: {unknown}')
        row['site_position'] = int(row['position'])
        row['residue'] = row['flank'][7]
        row['sources'] = sources
        for label in row['kinase'].split(';'):
            yield {**row, 'enzyme': label, **_kinase_fields(label)}
