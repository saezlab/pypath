"""
Parse MetaChatDB metabolite-sensor interactions and emit Relation records.

MetaChatDB is the literature-supported metabolite-sensor table of MetaChat
(human and mouse). A sensor is a receptor, a nuclear receptor or a transporter
of the metabolite. The files also carry the pairs MetaChat imported from
MEBOCOST, MRClinksDB and MetalinksDB; only the pairs MetaChatDB itself reports
are imported, the others come from those resources' own modules.
"""

from __future__ import annotations

import re
from typing import Any, Iterator

from biolink_model.datamodel.model import ChemicalEntity, Protein, slots
from omnipath_core.interaction_profiles import TRANSPORT_QUALIFIERS
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._class_terms import chemont_association, with_chemont
from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    AssociationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.METACHATDB,
    name='MetaChatDB',
    url='https://github.com/SonghaoLuo/MetaChat',
    license=LicenseCV.MIT,
    update_category=UpdateCategoryCV.IRREGULAR,
    topics=(Topic.CELL_CELL_COMMUNICATION, Topic.METABOLISM),
    content=(Content.COMPOUND_TARGET_INTERACTIONS, Content.MOLECULAR_TRANSPORT),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Literature-supported metabolite–sensor pairs: receptors, nuclear receptors and transporters.',
    description=(
        'MetaChatDB is a literature-supported database of metabolite-sensor '
        'interactions for human and mouse, compiled for the MetaChat '
        'metabolic cell-communication method.'
    ),
)

BASE_URL = (
    'https://raw.githubusercontent.com/SonghaoLuo/MetaChat/refs/heads/main/'
    'metachat/preprocessing/_data/MetaChatDB/MetaChatDB_%s.tsv'
)
TAXA = {'human': '9606', 'mouse': '10090'}
OWN_SOURCE = 'MetaChatDB'
# The values of `Sensor.Type`, one or two of them comma-separated.
TRANSPORTER = 'transporter'
SENSOR_TYPES = {'receptor', 'nuclear receptor', TRANSPORTER}
# The metabolite's ClassyFire superclass: a ChemOnt label.
CLASSYFIRE_FIELDS = ('Metabolite.Class',)
_EC = re.compile(r'^EC [\d.n-]+$')
_SECTION = re.compile(r'\s*\[(?:Cleaved into|Includes):[^\]]*\]')


def _values(value: Any, delimiter: str) -> list[str]:
    """Items of a delimited cell; the source writes an empty list as NA."""
    items = (item.strip() for item in str(value or '').split(delimiter))
    return [item for item in items if item and item != 'NA']


def iter_metachatdb(opener, **kwargs: Any) -> Iterator[dict[str, Any]]:
    """The rows MetaChatDB itself reports, one per sensor type.

    `Sensor.Type` names two types for some sensors (CD36 is a receptor and a
    transporter of fatty acids); each is a statement of its own.
    """
    for row in iter_tsv(opener, **kwargs):
        if OWN_SOURCE not in _values(row.get('Sources'), ';'):
            continue
        for sensor_type in _values(row.get('Sensor.Type'), ','):
            if sensor_type not in SENSOR_TYPES:
                raise ValueError(f'MetaChatDB: unknown sensor type {sensor_type!r}')
            yield {**row, 'sensor_type': sensor_type}


def _top_level_groups(text: str) -> tuple[str, list[str]]:
    """Split a UniProt protein name into the name and its parenthesised
    alternative names, keeping parentheses nested inside a name."""
    head, groups, depth, current = [], [], 0, []
    for char in text:
        if char == '(':
            if depth:
                current.append(char)
            depth += 1
        elif char == ')' and depth:
            depth -= 1
            if depth:
                current.append(char)
            else:
                groups.append(''.join(current).strip())
                current = []
        elif depth:
            current.append(char)
        else:
            head.append(char)
    return ''.join(head).strip(), groups


def _own_names(value: Any) -> str:
    """UniProt's protein names field without its `[Cleaved into: ...]` and
    `[Includes: ...]` sections, which name other chains or domains."""
    return _SECTION.sub('', str(value or '')).strip()


def protein_name(value: Any) -> str | None:
    """The recommended name of `Sensor.Name`, which is UniProt's protein names field."""
    name, _ = _top_level_groups(_own_names(value))
    return name or None


def protein_synonyms(value: Any) -> list[str]:
    """The alternative names of `Sensor.Name`; EC numbers are not names."""
    _, groups = _top_level_groups(_own_names(value))
    return [group for group in groups if group and not _EC.match(group)]


def _evidence(row: dict[str, Any]) -> dict[str, list[str]]:
    """`Evidences`: PubMed IDs, URLs and named source databases (Recon2)."""
    parsed: dict[str, list[str]] = {'pubmed': [], 'url': [], 'source': []}
    for item in _values(row.get('Evidences'), ';'):
        if curie := pmid_curie(item):
            parsed['pubmed'].append(curie)
        elif item.startswith(('http://', 'https://')):
            parsed['url'].append(item)
        else:
            parsed['source'].append(item)
    return parsed


def _supporting_sources(row: dict[str, Any]) -> list[str]:
    """The other databases `Sources` names, then those `Evidences` names."""
    sources = [s for s in _values(row.get('Sources'), ';') if s != OWN_SOURCE]
    return sources + [s for s in _evidence(row)['source'] if s not in sources]


def _is_transporter(row: dict[str, Any]) -> bool:
    return row.get('sensor_type') == TRANSPORTER


f = FieldConfig()


def interactions_schema(taxon_id: str) -> RelationBuilder:
    """A receptor binds its metabolite; a transporter moves it.

    Binding is metabolite `interacts_with` sensor. Transport is sensor
    `affects` metabolite with the transport qualifiers, as MEBOCOST and
    MRClinksDB state a transporter and its substrate.
    """
    metabolite = EntityBuilder(
        entity_type=ChemicalEntity,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.HMDB, value=f('HMDB.ID')),
            CV(term=Namespace.NAME, value=f('Metabolite.Name')),
        ),
        # `Metabolite.Class` is the ClassyFire superclass: a ChemOnt term.
        associations=AssociationsBuilder(chemont_association()),
    )
    sensor = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('Sensor.Gene')),
            CV(term=Namespace.NAME, value=f('Sensor.Name', transform=protein_name)),
            CV(term=Namespace.SYNONYM, value=lambda row: protein_synonyms(row.get('Sensor.Name'))),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
        ),
    )
    return RelationBuilder(
        subject=lambda row: (sensor if _is_transporter(row) else metabolite).build(row),
        predicate=lambda row: slots.affects if _is_transporter(row) else slots.interacts_with,
        object=lambda row: (metabolite if _is_transporter(row) else sensor).build(row),
        annotations=AnnotationsBuilder(
            *(
                CV(term=term, value=lambda row, value=value: value if _is_transporter(row) else None)
                for term, value in TRANSPORT_QUALIFIERS
            ),
            CV(term='metachatdb:sensor_type', value=f('sensor_type')),
            CV(term=slots.publications, value=lambda row: _evidence(row)['pubmed']),
            CV(term=slots.source_record_urls, value=lambda row: _evidence(row)['url']),
            CV(term=slots.supporting_data_source, value=_supporting_sources),
            # Context the table gives with each pair, kept as written.
            CV(
                term='metachatdb:sensor_pathway',
                value=lambda row: _values(row.get('Sensor.Pathway'), ';'),
            ),
            CV(
                term='metachatdb:metabolite_pathway',
                value=lambda row: _values(row.get('Metabolite.Pathway'), ';'),
            ),
            CV(
                term='metachatdb:metabolite_transporter',
                value=lambda row: _values(row.get('Metabolite.transporter'), ';'),
            ),
            CV(
                term='metachatdb:long_range_channel',
                value=lambda row: _values(row.get('Long.Range.Channel'), ','),
            ),
        ),
    )


def _dataset(species: str) -> Dataset:
    return Dataset(
        download=Download(
            url=BASE_URL % species,
            filename=f'MetaChatDB_{species}.tsv',
            subfolder='metachatdb',
            ext='tsv',
        ),
        mapper=interactions_schema(TAXA[species]),
        raw_parser=with_chemont(iter_metachatdb, CLASSYFIRE_FIELDS),
    )


resource = Resource(
    config,
    human=_dataset('human'),
    mouse=_dataset('mouse'),
)
