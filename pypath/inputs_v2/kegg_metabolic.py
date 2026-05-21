"""Parse KEGG metabolic reaction network data into reaction entities."""

from __future__ import annotations

import io
import logging

_log = logging.getLogger(__name__)

_KEGG_ORGANISM_CODES: dict[int, str] = {
    9606:  'hsa',   # Homo sapiens (human)
    10090: 'mmu',   # Mus musculus (mouse)
    10116: 'rno',   # Rattus norvegicus (rat)
    7955:  'dre',   # Danio rerio (zebrafish)
    6239:  'cel',   # Caenorhabditis elegans
    7227:  'dme',   # Drosophila melanogaster
    4932:  'sce',   # Saccharomyces cerevisiae
}

# Organism-agnostic KEGG REST endpoints (shared across all organisms).
_SHARED_URLS: dict[str, str] = {
    'link_reaction': 'https://rest.kegg.jp/link/reaction/enzyme',
    'list_compound': 'https://rest.kegg.jp/list/compound',
    'conv_chebi':    'https://rest.kegg.jp/conv/chebi/compound',
    'list_reaction': 'https://rest.kegg.jp/list/reaction',
}


def kegg_organism_code(ncbi_tax_id: int) -> str | None:
    """Return the KEGG organism code for an NCBI taxonomy ID, or None."""
    return _KEGG_ORGANISM_CODES.get(ncbi_tax_id)


# ---------------------------------------------------------------------------
# inputs_v2 schema objects (used by the omnipath-build pipeline).
# Guarded so that missing CV terms or API changes in the local pypath do
# not prevent this module from being imported in other contexts.
# ---------------------------------------------------------------------------

try:
    from pypath.internals.cv_terms import (
        BiologicalRoleCv,
        EntityTypeCv,
        IdentifierNamespaceCv,
        LicenseCV,
        ResourceCv,
        UpdateCategoryCV,
    )
    from pypath.internals.tabular_builder import (
        AnnotationsBuilder,
        CV,
        EntityBuilder,
        FieldConfig,
        IdentifiersBuilder,
        Member,
        MembershipBuilder,
    )
    from pypath.inputs_v2.base import Dataset, Resource, ResourceConfig
    from pypath.inputs_v2.parsers.kegg import _raw as _raw_parser

    config = ResourceConfig(
        id=ResourceCv.KEGG_METABOLIC,
        name='KEGG Metabolic Reactions (mouse)',
        url='https://www.genome.jp/kegg/reaction/',
        license=LicenseCV.KEGG_ACADEMIC,
        update_category=UpdateCategoryCV.REGULAR,
        pubmed='10592173',
        primary_category='metabolism',
        description=(
            'Mouse (mmu) metabolic reaction network from KEGG. '
            'Reactions are linked to their catalysing enzymes via UniProt IDs '
            'and to substrate / product metabolites via KEGG compound and ChEBI '
            'identifiers. Built from KEGG REST bulk endpoints '
            '(list, conv, link) using the dlmachine download infrastructure.'
        ),
    )

    f = FieldConfig(
        extract={
            'kegg_cpd': r'^(?:cpd:)?(.+)$',
            'chebi':    r'^(?:chebi:)?(\d+)$',
        },
    )

    reactions_schema = EntityBuilder(
        entity_type=EntityTypeCv.REACTION,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.KEGG_REACTION, value=f('reaction_id')),
        ),
        membership=MembershipBuilder(
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('uniprot_id')),
                    ),
                    annotations=AnnotationsBuilder(
                        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('taxon_id')),
                    ),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.CATALYST),
                ),
            ),
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.SMALL_MOLECULE,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('compound_kegg_id', extract='kegg_cpd')),
                        CV(term=IdentifierNamespaceCv.CHEBI,          value=f('compound_chebi',   extract='chebi')),
                        CV(term=IdentifierNamespaceCv.NAME,            value=f('compound_name')),
                    ),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.SUBSTRATE),
                ),
            ),
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.SMALL_MOLECULE,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('compound_kegg_id', extract='kegg_cpd')),
                        CV(term=IdentifierNamespaceCv.CHEBI,          value=f('compound_chebi',   extract='chebi')),
                        CV(term=IdentifierNamespaceCv.NAME,            value=f('compound_name')),
                    ),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.PRODUCT),
                ),
            ),
        ),
    )

except Exception as _e:  # pragma: no cover
    _log.debug('kegg_metabolic: schema objects unavailable (%s)', _e)
    config = None  # type: ignore[assignment]
    reactions_schema = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Simple multi-file downloader
# ---------------------------------------------------------------------------

class _MultiOpener:
    """Minimal opener whose .result is a dict of named text handles."""

    def __init__(self, result: dict) -> None:
        self.result = result


def _download_kegg_files(organism_code: str) -> _MultiOpener:
    """
    Download all six KEGG REST endpoints needed for the reaction network
    and return a :class:`_MultiOpener` whose ``.result`` dict maps each
    endpoint key to an open text handle (cached by dlmachine).

    Files are downloaded into ``<pypath_data_dir>/kegg_metabolic_{organism_code}/``
    and opened directly with plain ``open()`` to avoid the cachedir ``Opener``
    ``UnboundLocalError`` bug on plain TSV files.
    """
    from pathlib import Path
    from pypath.share.downloads import get_download_manager, _resolve_data_dir

    urls = {
        'conv_uniprot': f'https://rest.kegg.jp/conv/uniprot/{organism_code}',
        'link_enzyme':  f'https://rest.kegg.jp/link/enzyme/{organism_code}',
        **_SHARED_URLS,
    }
    subfolder = f'kegg_metabolic_{organism_code}'
    data_dir: Path = _resolve_data_dir() / subfolder
    data_dir.mkdir(parents=True, exist_ok=True)
    dm = get_download_manager()
    result: dict = {}

    for key, url in urls.items():
        filename = (
            f'kegg_{key}_{organism_code}.tsv'
            if key in ('conv_uniprot', 'link_enzyme')
            else f'kegg_{key}.tsv'
        )
        file_path = data_dir / filename
        _log.debug('[KEGG] %s → %s', key, file_path)
        dm.download(url, dest=str(file_path))
        result[key] = open(file_path, encoding='utf-8')  # noqa: SIM115

    return _MultiOpener(result)


# ---------------------------------------------------------------------------
# Lightweight Resource / Dataset wrappers
# ---------------------------------------------------------------------------

class _KeggDataset:
    """Minimal dataset whose raw() downloads KEGG REST files and yields reaction dicts."""

    def __init__(self, organism_code: str) -> None:
        self._org = organism_code

    def raw(self, **_kwargs):
        from pypath.inputs_v2.parsers.kegg import _raw
        opener = _download_kegg_files(self._org)
        yield from _raw(opener, data_type='reactions', organism=self._org)


class _KeggResource:
    """Minimal Resource wrapper with a ``reactions`` dataset attribute."""

    def __init__(self, organism_code: str) -> None:
        self.reactions = _KeggDataset(organism_code)


# ---------------------------------------------------------------------------
# Public factory
# ---------------------------------------------------------------------------

def make_kegg_resource(organism_code: str = 'mmu') -> _KeggResource:
    """
    Return a KEGG resource configured for the given KEGG organism code.

    Downloads all six KEGG REST endpoints (conv/uniprot, link/enzyme,
    link/reaction, list/compound, conv/chebi, list/reaction) via the
    dlmachine cache infrastructure.  Each file is stored under
    ``kegg_metabolic_{organism_code}/`` in the pypath data directory.

    Args:
        organism_code: KEGG organism code (e.g. ``'hsa'``, ``'mmu'``, ``'rno'``).
            Use :func:`kegg_organism_code` to convert an NCBI taxonomy ID.

    Returns:
        Object with a ``reactions`` attribute whose ``raw()`` method yields
        flat reaction dicts ready for conversion to COSMOS Interaction records.
    """
    return _KeggResource(organism_code)
