"""Parse KEGG metabolic reaction network data into reaction entities."""

from __future__ import annotations

from collections.abc import Iterable
import logging
import time
import urllib.request

_log = logging.getLogger(__name__)

_KEGG_REST = 'https://rest.kegg.jp'
# KEGG REST `get` flat-file requests are limited to 10 entries.
_KEGG_BATCH_SIZE = 10
_KEGG_REST_DELAY = 0.0

_KEGG_ORGANISM_CODES: dict[int, str] = {
    9606:  'hsa',   # Homo sapiens (human)
    10090: 'mmu',   # Mus musculus (mouse)
    10116: 'rno',   # Rattus norvegicus (rat)
    7955:  'dre',   # Danio rerio (zebrafish)
    6239:  'cel',   # Caenorhabditis elegans
    7227:  'dme',   # Drosophila melanogaster
    4932:  'sce',   # Saccharomyces cerevisiae
}
_DEFAULT_KEGG_ORGANISMS: tuple[str, ...] = tuple(_KEGG_ORGANISM_CODES.values())

# Organism-agnostic KEGG REST endpoints (shared across all organisms).
_SHARED_URLS: dict[str, str] = {
    'link_reaction': 'https://rest.kegg.jp/link/reaction/enzyme',
    'list_compound': 'https://rest.kegg.jp/list/compound',
    'conv_chebi':    'https://rest.kegg.jp/conv/chebi/compound',
    'conv_pubchem':  'https://rest.kegg.jp/conv/pubchem/compound',
    'list_reaction': 'https://rest.kegg.jp/list/reaction',
}
_SHARED_SUBFOLDER = 'kegg_shared'


def kegg_organism_code(ncbi_tax_id: int) -> str | None:
    """Return the KEGG organism code for an NCBI taxonomy ID, or None."""
    return _KEGG_ORGANISM_CODES.get(ncbi_tax_id)


def _kegg_organism_codes(
    organism: str | int | Iterable[str | int] | None = None,
) -> tuple[str, ...]:
    """Return KEGG organism codes, defaulting to all configured organisms."""
    if organism is None:
        return _DEFAULT_KEGG_ORGANISMS

    if isinstance(organism, str):
        if organism in {'all', '*'}:
            return _DEFAULT_KEGG_ORGANISMS
        return (organism,)

    if isinstance(organism, int):
        organism_code = kegg_organism_code(organism)
        if organism_code is None:
            raise ValueError(f'Unsupported KEGG organism taxon ID: {organism}')
        return (organism_code,)

    organism_codes: list[str] = []
    for item in organism:
        organism_codes.extend(_kegg_organism_codes(item))
    return tuple(dict.fromkeys(organism_codes))


# ---------------------------------------------------------------------------
# inputs_v2 schema objects (used by the omnipath-build pipeline).
# Guarded so that missing CV terms or API changes in the local pypath do
# not prevent this module from being imported in other contexts.
# ---------------------------------------------------------------------------

try:
    from pypath.inputs_v2.base import Dataset, Resource, ResourceConfig
    from pypath.internals.cv_terms import (
        BiologicalRoleCv,
        EntityTypeCv,
        IdentifierNamespaceCv,
        InteractionMetadataCv,
        LicenseCV,
        MoleculeAnnotationsCv,
        ParticipantMetadataCv,
        ResourceCv,
        UpdateCategoryCV,
    )
    from pypath.internals.tabular_builder import (
        CV,
        AnnotationsBuilder,
        EntityBuilder,
        FieldConfig,
        IdentifiersBuilder,
        MembersFromList,
        MembershipBuilder,
    )

    config = ResourceConfig(
        id=ResourceCv.KEGG_METABOLIC,
        name='KEGG Metabolic Reactions',
        url='https://www.genome.jp/kegg/reaction/',
        license=LicenseCV.KEGG_ACADEMIC,
        update_category=UpdateCategoryCV.REGULAR,
        pubmed='10592173',
        primary_category='metabolism',
        description=(
            'Biochemical reactions from KEGG REACTION, with exact KEGG '
            'compound participants parsed from reaction equations and enriched '
            'with EC, Rhea, KO, RCLASS, ChEBI, and PubChem identifiers.'
        ),
    )

    f = FieldConfig(
        extract={
            'kegg_cpd': r'^(?:cpd:)?(C\d{5})$',
            'chebi':    r'^(?:chebi:)?(\d+)$',
            'pubchem':  r'^(?:pubchem:)?(\d+)$',
            'rhea':     r'^(?:RHEA:)?(\d+)$',
        },
    )

    reactions_schema = EntityBuilder(
        entity_type=EntityTypeCv.REACTION,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.KEGG_REACTION, value=f('reaction_id')),
            CV(term=IdentifierNamespaceCv.NAME, value=f('reaction_name', delimiter=';')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=InteractionMetadataCv.CONVERSION_DIRECTION, value=f('conversion_direction')),
            CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('reaction_equation')),
            CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('reaction_definition')),
            CV(term=MoleculeAnnotationsCv.EC_NUMBER, value=f('ec_numbers', delimiter=';')),
            CV(term=IdentifierNamespaceCv.RHEA_ID, value=f('rhea_ids', delimiter=';', extract='rhea')),
            CV(term=IdentifierNamespaceCv.KEGG, value=f('ko_ids', delimiter=';')),
            CV(term=IdentifierNamespaceCv.KEGG, value=f('rclass_ids', delimiter=';')),
            CV(term=IdentifierNamespaceCv.KEGG_PATHWAY, value=f('pathway_ids', delimiter=';')),
            CV(term=MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value=f('pathway_names', delimiter=';')),
        ),
        membership=MembershipBuilder(
            MembersFromList(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('reactant_kegg_id', delimiter='||', extract='kegg_cpd', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.NAME,          value=f('reactant_name',    delimiter='||', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.CHEBI,         value=f('reactant_chebi',   delimiter='||', extract='chebi', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.PUBCHEM,       value=f('reactant_pubchem', delimiter='||', extract='pubchem', preserve_indices=True)),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.REACTANT),
                    CV(term=ParticipantMetadataCv.STOICHIOMETRY, value=f('reactant_stoichiometry', delimiter='||', preserve_indices=True)),
                ),
            ),
            MembersFromList(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('product_kegg_id', delimiter='||', extract='kegg_cpd', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.NAME,          value=f('product_name',    delimiter='||', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.CHEBI,         value=f('product_chebi',   delimiter='||', extract='chebi', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.PUBCHEM,       value=f('product_pubchem', delimiter='||', extract='pubchem', preserve_indices=True)),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.PRODUCT),
                    CV(term=ParticipantMetadataCv.STOICHIOMETRY, value=f('product_stoichiometry', delimiter='||', preserve_indices=True)),
                ),
            ),
            MembersFromList(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('uniprot_ids', delimiter=';')),
                ),
                entity_annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('taxon_id')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.CATALYST),
                ),
            ),
        ),
    )

    pathways_schema = EntityBuilder(
        entity_type=EntityTypeCv.PATHWAY,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.KEGG_PATHWAY, value=f('pathway_id')),
            CV(term=IdentifierNamespaceCv.NAME, value=f('pathway_name')),
        ),
        membership=MembershipBuilder(
            MembersFromList(
                entity_type=EntityTypeCv.REACTION,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.KEGG_REACTION, value=f('reaction_ids', delimiter=';', preserve_indices=True)),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.PATHWAY_COMPONENT),
                ),
            ),
            MembersFromList(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.KEGG, value=f('protein_member_kegg_ids', delimiter='||', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.ENTREZ, value=f('protein_member_entrez_ids', delimiter='||', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('protein_member_uniprot_ids', delimiter='||', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('protein_member_names', delimiter='||', preserve_indices=True)),
                ),
                entity_annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('taxon_id')),
                    CV(term=IdentifierNamespaceCv.KEGG_REACTION, value=f('protein_member_reaction_ids', delimiter='||', preserve_indices=True)),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.PATHWAY_COMPONENT),
                ),
            ),
            MembersFromList(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('small_molecule_member_kegg_ids', delimiter='||', extract='kegg_cpd', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('small_molecule_member_names', delimiter='||', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.CHEBI, value=f('small_molecule_member_chebi_ids', delimiter='||', extract='chebi', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.PUBCHEM, value=f('small_molecule_member_pubchem_ids', delimiter='||', extract='pubchem', preserve_indices=True)),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.PATHWAY_COMPONENT),
                ),
            ),
            MembersFromList(
                entity_type=EntityTypeCv.CV_TERM,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.KEGG, value=f('ortholog_member_kegg_ids', delimiter='||', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('ortholog_member_names', delimiter='||', preserve_indices=True)),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.PATHWAY_COMPONENT),
                ),
            ),
            MembersFromList(
                entity_type=EntityTypeCv.PATHWAY,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.KEGG_PATHWAY, value=f('linked_pathway_ids', delimiter='||', preserve_indices=True)),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('linked_pathway_names', delimiter='||', preserve_indices=True)),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=BiologicalRoleCv.PATHWAY_COMPONENT),
                ),
            ),
        ),
    )

except Exception as _e:  # pragma: no cover
    _log.debug('kegg: schema objects unavailable (%s)', _e)
    config = None  # type: ignore[assignment]
    reactions_schema = None  # type: ignore[assignment]
    pathways_schema = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Simple multi-file downloader
# ---------------------------------------------------------------------------

class _MultiOpener:
    """Minimal opener whose .result is a dict of named text handles."""

    def __init__(self, result: dict) -> None:
        self.result = result


def _reaction_ids(list_reaction_handle) -> list[str]:
    """Return bare KEGG reaction IDs from a /list/reaction handle."""
    if hasattr(list_reaction_handle, 'seek'):
        list_reaction_handle.seek(0)
    reaction_ids: list[str] = []
    for line in list_reaction_handle:
        if isinstance(line, bytes):
            line = line.decode('utf-8', 'ignore')
        if not line.strip():
            continue
        reaction = line.split('\t', 1)[0].strip()
        if reaction.startswith('rn:'):
            reaction = reaction[3:]
        reaction_ids.append(reaction)
    if hasattr(list_reaction_handle, 'seek'):
        list_reaction_handle.seek(0)
    return reaction_ids


def _completed_kegg_entries(path) -> int:
    """Count completed KEGG flat-file entries in a combined cache file."""
    if not path.exists() or path.stat().st_size == 0:
        return 0
    count = 0
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('ENTRY'):
                count += 1
    return count


def _fetch_reaction_details(
    reaction_ids: list[str],
    dest: str,
    delay: float = _KEGG_REST_DELAY,
) -> None:
    """Fetch KEGG reaction flat files in REST batches and write one combined file."""
    from pathlib import Path

    dest_path = Path(dest)
    tmp_path = dest_path.with_suffix(dest_path.suffix + '.tmp')
    expected = len(reaction_ids)

    completed = _completed_kegg_entries(dest_path)
    if completed == expected:
        return
    if completed and not tmp_path.exists():
        dest_path.replace(tmp_path)
    completed = _completed_kegg_entries(tmp_path)
    if completed > expected:
        tmp_path.unlink(missing_ok=True)
        completed = 0

    # Resume only at a batch boundary so a killed process cannot leave a
    # partial batch that shifts subsequent reaction IDs.
    if completed % _KEGG_BATCH_SIZE:
        tmp_path.unlink(missing_ok=True)
        completed = 0
    mode = 'a' if completed else 'w'
    _log.info('[KEGG] fetching reaction entries %d/%d', completed, expected)

    with open(tmp_path, mode, encoding='utf-8') as out:
        for start in range(completed, expected, _KEGG_BATCH_SIZE):
            batch = reaction_ids[start:start + _KEGG_BATCH_SIZE]
            url = f"{_KEGG_REST}/get/{'+'.join(batch)}"
            request = urllib.request.Request(
                url,
                headers={'User-Agent': 'pypath-kegg/1.0'},
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                out.write(response.read().decode('utf-8'))
            if delay:
                time.sleep(delay)
            if (start + len(batch)) % 500 == 0 or start + len(batch) >= expected:
                _log.info(
                    '[KEGG] fetched reaction entries %d/%d',
                    start + len(batch),
                    expected,
                )

    if _completed_kegg_entries(tmp_path) != expected:
        raise RuntimeError(
            f'Incomplete KEGG reaction cache: '
            f'{_completed_kegg_entries(tmp_path)}/{expected} entries'
        )
    tmp_path.replace(dest_path)


def _pathway_ids(list_pathway_handle) -> list[str]:
    """Return KEGG pathway IDs from a /list/pathway/{organism} handle."""
    if hasattr(list_pathway_handle, 'seek'):
        list_pathway_handle.seek(0)
    pathway_ids: list[str] = []
    for line in list_pathway_handle:
        if isinstance(line, bytes):
            line = line.decode('utf-8', 'ignore')
        if not line.strip():
            continue
        pathway = line.split('\t', 1)[0].strip()
        if pathway.startswith('path:'):
            pathway = pathway[5:]
        pathway_ids.append(pathway)
    if hasattr(list_pathway_handle, 'seek'):
        list_pathway_handle.seek(0)
    return pathway_ids


def _fetch_pathway_kgml(
    pathway_ids: list[str],
    dest_dir,
    *,
    force_refresh: bool = False,
    delay: float = _KEGG_REST_DELAY,
) -> None:
    """Fetch organism-specific KGML pathway files one by one."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    for index, pathway_id in enumerate(pathway_ids, start=1):
        file_path = dest_dir / f'{pathway_id}.kgml'
        if not force_refresh and file_path.exists() and file_path.stat().st_size > 0:
            continue
        url = f'{_KEGG_REST}/get/{pathway_id}/kgml'
        request = urllib.request.Request(
            url,
            headers={'User-Agent': 'pypath-kegg/1.0'},
        )
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                file_path.write_text(response.read().decode('utf-8'), encoding='utf-8')
        except Exception as exc:  # pragma: no cover - network/API defensive path
            _log.warning('[KEGG] failed to fetch KGML for %s: %s', pathway_id, exc)
        if delay:
            time.sleep(delay)
        if index % 50 == 0 or index == len(pathway_ids):
            _log.info(
                '[KEGG] fetched pathway KGML files %d/%d',
                index,
                len(pathway_ids),
            )


def _download_kegg_files(
    organism_code: str,
    *,
    force_refresh: bool = False,
) -> _MultiOpener:
    """
    Download KEGG REST endpoints needed for the reaction network
    and return a :class:`_MultiOpener` whose ``.result`` dict maps each
    endpoint key to an open text handle (cached by dlmachine).

    Organism-specific files are downloaded into
    ``<pypath_data_dir>/kegg_{organism_code}/``. Shared KEGG reaction and
    compound files are downloaded once into ``<pypath_data_dir>/kegg_shared/``.
    Files are opened directly with plain ``open()`` to avoid the cachedir
    ``Opener`` ``UnboundLocalError`` bug on plain TSV files.
    """
    from pathlib import Path

    from pypath.share.downloads import _resolve_data_dir, get_download_manager

    organism_urls = {
        'conv_uniprot': f'https://rest.kegg.jp/conv/uniprot/{organism_code}',
        'link_enzyme':  f'https://rest.kegg.jp/link/enzyme/{organism_code}',
        'list_pathway': f'https://rest.kegg.jp/list/pathway/{organism_code}',
    }
    base_dir: Path = _resolve_data_dir()
    data_dir = base_dir / f'kegg_{organism_code}'
    shared_dir = base_dir / _SHARED_SUBFOLDER
    data_dir.mkdir(parents=True, exist_ok=True)
    shared_dir.mkdir(parents=True, exist_ok=True)
    dm = get_download_manager()
    result: dict = {}

    for key, url in organism_urls.items():
        filename = f'kegg_{key}_{organism_code}.tsv'
        file_path = data_dir / filename
        _log.debug('[KEGG] %s → %s', key, file_path)
        dm.download(url, dest=str(file_path))
        result[key] = open(file_path, encoding='utf-8')  # noqa: SIM115

    for key, url in _SHARED_URLS.items():
        file_path = shared_dir / f'kegg_{key}.tsv'
        _log.debug('[KEGG] %s → %s', key, file_path)
        dm.download(url, dest=str(file_path))
        result[key] = open(file_path, encoding='utf-8')  # noqa: SIM115

    details_path = shared_dir / 'kegg_get_reaction.txt'
    reaction_ids = _reaction_ids(result['list_reaction'])
    if (
        force_refresh
        or _completed_kegg_entries(details_path) != len(reaction_ids)
    ):
        _log.info('[KEGG] fetching %d reaction flat-file entries', len(reaction_ids))
        _fetch_reaction_details(reaction_ids, str(details_path))
    result['get_reaction'] = open(details_path, encoding='utf-8')  # noqa: SIM115

    kgml_dir = data_dir / 'kgml'
    pathway_ids = _pathway_ids(result['list_pathway'])
    _log.info('[KEGG] ensuring %d pathway KGML files', len(pathway_ids))
    _fetch_pathway_kgml(pathway_ids, kgml_dir, force_refresh=force_refresh)
    result['kgml_pathways'] = {
        pathway_id: open(kgml_dir / f'{pathway_id}.kgml', encoding='utf-8')  # noqa: SIM115
        for pathway_id in pathway_ids
        if (kgml_dir / f'{pathway_id}.kgml').exists()
    }

    return _MultiOpener(result)


# ---------------------------------------------------------------------------
# Lightweight Resource / Dataset wrappers
# ---------------------------------------------------------------------------

class _KeggDataset:
    """Minimal dataset whose raw() downloads KEGG REST files and yields rows."""

    def __init__(
        self,
        organism_codes: str | int | Iterable[str | int] | None,
        data_type: str,
    ) -> None:
        self._orgs = _kegg_organism_codes(organism_codes)
        self._data_type = data_type

    def raw(self, **kwargs):
        from pypath.inputs_v2.parsers.kegg import _raw
        organism_codes = _kegg_organism_codes(kwargs.pop('organism', self._orgs))
        for organism in organism_codes:
            opener = _download_kegg_files(
                organism,
                force_refresh=bool(kwargs.get('force_refresh', False)),
            )
            yield from _raw(opener, data_type=self._data_type, organism=organism, **kwargs)


class _KeggResource:
    """Minimal Resource wrapper with ``reactions`` and ``pathways`` datasets."""

    def __init__(
        self,
        organism_codes: str | int | Iterable[str | int] | None,
    ) -> None:
        self.reactions = _KeggDataset(organism_codes, 'reactions')
        self.pathways = _KeggDataset(organism_codes, 'pathways')


# ---------------------------------------------------------------------------
# Public factory
# ---------------------------------------------------------------------------

def make_kegg_resource(
    organism_code: str | int | Iterable[str | int] | None = None,
) -> _KeggResource:
    """
    Return a KEGG resource configured for one or more KEGG organisms.

    Downloads KEGG REST endpoints for organism-specific gene mappings,
    shared compound cross-references, and batched reaction flat files via the
    dlmachine cache infrastructure. Organism-specific files are stored under
    ``kegg_{organism_code}/`` and shared files under ``kegg_shared/`` in the
    pypath data directory.

    Args:
        organism_code: KEGG organism code, NCBI taxonomy ID, iterable of either,
            ``'all'``, or ``None``. ``None`` and ``'all'`` load all configured
            organism codes. Use :func:`kegg_organism_code` to convert an NCBI
            taxonomy ID explicitly.

    Returns:
        Object with ``reactions`` and ``pathways`` attributes whose ``raw()``
        methods yield flat rows ready for conversion to entities.
    """
    return _KeggResource(organism_code)


def _default_raw(data_type: str, opener=None, **kwargs):
    from pypath.inputs_v2.parsers.kegg import _raw
    requested_organism = kwargs.pop('organism', None)
    if opener is None:
        organism_codes = _kegg_organism_codes(requested_organism)
        for organism in organism_codes:
            organism_opener = _download_kegg_files(
                organism,
                force_refresh=bool(kwargs.get('force_refresh', False)),
            )
            yield from _raw(
                organism_opener,
                data_type=data_type,
                organism=organism,
                **kwargs,
            )
        return

    organism_codes = _kegg_organism_codes(requested_organism or 'mmu')
    if len(organism_codes) != 1:
        raise ValueError('A single organism must be specified when passing opener.')
    yield from _raw(opener, data_type=data_type, organism=organism_codes[0], **kwargs)


def _default_reactions_raw(opener=None, **kwargs):
    yield from _default_raw('reactions', opener=opener, **kwargs)


def _default_pathways_raw(opener=None, **kwargs):
    yield from _default_raw('pathways', opener=opener, **kwargs)


if config is not None and reactions_schema is not None and pathways_schema is not None:
    resource = Resource(
        config,
        reactions=Dataset(
            download=None,
            mapper=reactions_schema,
            raw_parser=_default_reactions_raw,
        ),
        pathways=Dataset(
            download=None,
            mapper=pathways_schema,
            raw_parser=_default_pathways_raw,
        ),
    )
else:  # pragma: no cover
    resource = None  # type: ignore[assignment]
