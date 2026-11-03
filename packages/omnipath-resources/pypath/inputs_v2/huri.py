"""
Parse the CCSB human binary interactome datasets and emit Relation records.

The Human Reference Interactome portal (interactome-atlas.org) offers each
dataset of the Center for Cancer Systems Biology as a PSI-MI TAB 2.7 file,
with one row per ORF pair and assay, and as a two-column table of Ensembl
gene pairs. Neither has a header line. The systematic yeast two-hybrid datasets are read from
their PSI-MI TAB files. The test space screens and Lit-BM exist only as gene
pairs. HuRI.psi also contains the pairs of the test space screens, marked as
such in each row; they are left out of the ``huri`` dataset, which is the
published HuRI map, and are the content of ``test_space_screens_2019``.

HI-union is not a dataset here: it is the union of the screens below, and
importing it would count every observation twice.
"""

from __future__ import annotations

import functools
import re
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2 import _mitab
from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.parsers.huri import MITAB_COLUMNS, iter_screen, parse_screens


config = ResourceConfig(
    id=ResourceCv.HURI,
    name='HuRI',
    full='Human Reference Interactome',
    url='https://interactome-atlas.org/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='32296183',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.LITERATURE_CURATED,
    ),
    summary=(
        'Human protein-protein interactions from the CCSB yeast two-hybrid '
        'screens (HuRI and its predecessors) and the Lit-BM literature set, '
        'as physical associations.'
    ),
    description=(
        'The Human Reference Interactome (HuRI) is a map of about 53,000 '
        'binary protein-protein interactions found by systematically '
        'screening nearly all pairs of human protein-coding genes with '
        'yeast two-hybrid assays and validating them in orthogonal assays. '
        'Its portal also distributes the earlier systematic maps of the '
        'Center for Cancer Systems Biology (H-I-05, Venkatesan-09, Yu-11, '
        'HI-II-14, Yang-16) and Lit-BM, a set of binary interactions from '
        'the literature that are supported by at least two pieces of '
        'experimental evidence. The PSI-MI files type every interaction '
        '"physical association" (MI:0915), which does not say that the '
        'proteins touch, so the module states that the proteins interact and '
        'not a direct physical interaction; the gene pair tables state no '
        'interaction type and are read the same way.'
    ),
)

_BASE_URL = 'https://interactome-atlas.org/data/'
_HUMAN = 'NCBITaxon:9606'

# PSI-MI TAB 2.7 columns, in file order (see the parser).
_MITAB_COLUMNS = MITAB_COLUMNS
# The two-column files of Ensembl gene pairs have no header line either.
_GENE_PAIR_COLUMNS = ('ensembl_gene_a', 'ensembl_gene_b')

_INTERACTOR_TYPES = {'MI:0326': Protein}
# The files give the prey the bait's accession (MI:0496); the label names the role.
_EXPERIMENTAL_ROLES = {
    'psi-mi:MI:0496(bait)': 'MI:0496',
    'psi-mi:MI:0496(prey)': 'MI:0498',
}
# HuRI.psi predates the paper and carries a placeholder. The portal's download
# table gives Luck et al., Nature 2020 as the publication of that dataset.
_UNASSIGNED_PUBLICATIONS = {'unassigned1304': 'PMID:32296183'}
_EXPERIMENTAL_ROLE = 'psi_mi:experimental_role'
# The screens of the portal that found the pair, as the numbers the file
# gives (``comment:Found in screens 4,7.``): an attribute of the observation.
_SCREEN = 'huri:screen'
_ENSEMBL_NAMESPACES = {
    'G': Namespace.ENSG,
    'T': Namespace.ENST,
    'P': Namespace.ENSP,
}
_ENSEMBL_RE = re.compile(r'ENS([GTP])\d{11}(?:\.\d+)?')


def _primary_identifier(value: object) -> tuple[Any, str] | None:
    prefix, _, accession = str(value or '').partition(':')
    if prefix == 'uniprotkb' and accession:
        return Namespace.UNIPROT, accession
    return _ensembl_identifier(value)


def _ensembl_identifier(value: object) -> tuple[Any, str] | None:
    prefix, _, accession = str(value or '').partition(':')
    match = _ENSEMBL_RE.fullmatch(accession) if prefix == 'ensembl' else None
    return (_ENSEMBL_NAMESPACES[match.group(1)], accession) if match else None


def _identifier_pairs(suffix: str):
    """The tested ORF's protein accession, then its Ensembl cross-references."""

    def _pairs(row: dict[str, Any]) -> list[tuple[Any, str]]:
        primary = _primary_identifier(row.get(f'id_{suffix}'))
        if primary is None:
            return []
        alternatives = str(row.get(f'alt_ids_{suffix}') or '').split('|')
        pairs = [primary, *filter(None, map(_ensembl_identifier, alternatives))]
        return list(dict.fromkeys(pairs))

    return _pairs


def _publication(value: str) -> str | None:
    prefix, _, accession = value.partition(':')
    if prefix == 'pubmed' and (curie := pmid_curie(accession)):
        return curie
    return _UNASSIGNED_PUBLICATIONS.get(value)


f = FieldConfig(
    extract={
        'mi': r'^psi-mi:(MI:\d{4})\(',
        'tax': r'^taxid:(\d+)\(',
        'author_score': r'^author score: ?([0-9.eE+-]+)$',
        'ensembl_gene': r'^(ENSG\d{11})$',
    },
    map={
        'interactor_type': _INTERACTOR_TYPES,
        'experimental_role': _EXPERIMENTAL_ROLES,
    },
    transform={
        'taxon': lambda v: f'NCBITaxon:{v}',
        'publication': _publication,
        'score': lambda v: measurement(v, source_field='author score'),
    },
    delimiter='|',
)


def _orf_product(suffix: str) -> EntityBuilder:
    return EntityBuilder(
        # A participant of another PSI-MI interactor type is not emitted.
        entity_type=f(
            f'interactor_type_{suffix}',
            extract='mi',
            map='interactor_type',
        ),
        identifiers=IdentifiersBuilder(
            CV.from_pairs(_identifier_pairs(suffix)),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.in_taxon,
                value=f(f'taxon_{suffix}', extract='tax', transform='taxon'),
            ),
            CV(
                term=_EXPERIMENTAL_ROLE,
                value=f(
                    f'experimental_role_{suffix}',
                    map='experimental_role',
                ),
            ),
        ),
    )


def _gene_product(suffix: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.ENSG,
                value=f(f'ensembl_gene_{suffix}', extract='ensembl_gene'),
            ),
        ),
        # The portal maps the human interactome only; ENSG is the human prefix.
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_HUMAN)),
    )


def _screens(row: dict[str, Any]) -> list[str]:
    """The numbers of the screens in which the pair was found."""
    return parse_screens(row.get('annotations_interaction'))[0]


def _is_positive(row: dict[str, Any]) -> bool:
    """False for a row the format marks as a negative result, or marks oddly."""
    return str(row.get('negative') or '-').lower() in {'-', 'false'}


def _predicate(row: dict[str, Any]):
    """The predicate that the PSI-MI interaction type of the row supports.

    Every row of the files states "physical association" (MI:0915), which gives
    ``interacts_with``; a type that states a direct interaction would give
    ``directly_physically_interacts_with``. The papers' description of the
    assay is not a statement of the row. Any other type is an error.
    """
    types = _mitab.mi_terms(row.get('interaction_type'))
    predicate = _mitab.physical_predicate(types)
    if predicate is None:
        raise ValueError(
            f'Unknown PSI-MI interaction type in a CCSB file: '
            f'{row.get("interaction_type")!r}'
        )
    return predicate


# One relation per ORF pair and assay.
screen_schema = RelationBuilder(
    subject=_orf_product('a'),
    predicate=_predicate,
    object=_orf_product('b'),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_evidence_of_type,
            value=f('detection_method', extract='mi'),
        ),
        # The interaction type the row states, kept as the source's predicate.
        CV(
            term=slots.original_predicate,
            value=_mitab.mi_terms_of('interaction_type'),
        ),
        CV(
            term=slots.publications,
            value=f('publication', transform='publication'),
        ),
        CV(
            term=slots.has_confidence_score,
            value=f('confidence', extract='author_score', transform='score'),
        ),
        CV(term=_SCREEN, value=_screens),
    ),
)


def map_screen_interaction(row: dict[str, Any]):
    # A negative result must never become an interaction.
    if not _is_positive(row):
        return None
    # A plain function, so that an unknown interaction type is raised and not
    # swallowed by the builder.
    _predicate(row)
    return screen_schema(row)


def _gene_pair_schema(predicate, publication: str | None = None):
    annotations = (
        [CV(term=slots.publications, value=publication)] if publication else []
    )
    return RelationBuilder(
        subject=_gene_product('a'),
        predicate=predicate,
        object=_gene_product('b'),
        annotations=AnnotationsBuilder(*annotations),
    )


# The gene pair tables state no interaction type for a pair, so a pair is an
# interaction and not a direct one, whatever the assays behind it.
#
# Pairs found in the yeast two-hybrid screens of the assay test space, which
# the portal's download table attributes to Luck et al., Nature 2020.
test_space_schema = _gene_pair_schema(
    slots.interacts_with,
    publication='PMID:32296183',
)
# Pairs from the literature with at least two pieces of evidence, one of them
# from a binary assay. The file names no method or publication per pair.
lit_bm_schema = _gene_pair_schema(slots.interacts_with)


def _dataset(
    filename: str, mapper, columns: tuple[str, ...], raw_parser=None
) -> Dataset:
    return Dataset(
        download=Download(
            url=_BASE_URL + filename,
            filename=filename,
            subfolder='huri',
            large=True,
        ),
        mapper=mapper,
        # The files have no header line. A row with another number of fields
        # is an error: the layout changed or the download is not the table.
        raw_parser=raw_parser
        or functools.partial(iter_tsv, fieldnames=columns, strict=True),
    )


def _screen(filename: str, exclude_test_space: bool = False) -> Dataset:
    return _dataset(
        filename,
        map_screen_interaction,
        _MITAB_COLUMNS,
        functools.partial(iter_screen, exclude_test_space=exclude_test_space),
    )


def _gene_pairs(filename: str, mapper) -> Dataset:
    return _dataset(filename, mapper, _GENE_PAIR_COLUMNS)


resource = Resource(
    config,
    # HuRI.psi also holds the pairs of the test space screens; they are the
    # content of the test_space_screens_2019 dataset, not of the HuRI map.
    huri=_screen('HuRI.psi', exclude_test_space=True),
    hi_ii_14=_screen('HI-II-14.psi'),
    yang_2016=_screen('Yang-16.psi'),
    yu_2011=_screen('Yu-11.psi'),
    venkatesan_2009=_screen('Venkatesan-09.psi'),
    hi_i_05=_screen('H-I-05.psi'),
    test_space_screens_2019=_gene_pairs(
        'Test_space_screens-19.tsv',
        test_space_schema,
    ),
    lit_bm=_gene_pairs('Lit-BM.tsv', lit_bm_schema),
)
