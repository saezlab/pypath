"""Source context stays on evidence without changing biological statement identity."""

import importlib
from pathlib import Path

from omnipath_resolver.resolver import EntityResolver
from omnipath_build.silver import SilverExtractor
from omnipath_build.writer import ParquetWriter
from omnipath_core.biolink import direction_sign, qualifiers
from omnipath_core.keys import relation_key
from omnipath_core.source_attributes import (
    CELLULAR_LOCATION,
    CONVERSION_DIRECTION,
    PARTICIPANT_ROLE,
    TRAIT_TYPE,
)
import pyarrow.parquet as pq
import pytest

from pypath.internals.tabular_builder import EntityBuilder, RelationBuilder

REACTION_CASES = [
    (module, mapper, row)
    for module, mappers, row in (
        (
            'rhea',
            ('reactions_schema', 'transport_reactions_schema'),
            {
                'rhea_id': '123',
                'participant_chebi': '1||2||3',
                'participant_role': 'reactant||reactant||product',
                'participant_compartment': 'in||in||out',
                'uniprot': 'P12345',
            },
        ),
        (
            'recon3d',
            ('reactions_schema', 'transport_reactions_schema'),
            {
                'bigg_reaction_id': 'RX1',
                'reactants': 'a:c:2||b:c:n',
                'products': 'a:e:(n+1)',
            },
        ),
        (
            'metatlas',
            ('reactions_schema', 'transport_reactions_schema'),
            {
                'human_gem_reaction_id': 'RX1',
                'reactants': 'a:c:2||b:c:n',
                'products': 'a:e:(n+1)',
            },
        ),
        (
            'kegg',
            ('reactions_schema',),
            {
                'reaction_id': 'R00001',
                'reactant_kegg_id': 'C00001||C00002',
                'reactant_stoichiometry': '2||n',
                'product_kegg_id': 'C00001',
                'product_stoichiometry': '(n+1)',
                'uniprot_ids': 'P12345',
            },
        ),
        (
            'reactome',
            ('reactions_schema',),
            {
                'entity_type': 'reaction',
                'reactome_stable_id': 'R-HSA-1',
                'participant_entity_type': 'chemical||chemical||protein',
                'participant_role': 'reactant||product||reactant',
                'participant_chebi': '1||2||__MISSING__',
                'participant_uniprot': '__MISSING__||__MISSING__||P12345',
                'participant_compartment': 'cytosol||extracellular||cytosol',
                'participant_stoichiometry': 'n||(n+1)||1',
            },
        ),
    )
    for mapper in mappers
]


def _mapper(module: str, name: str) -> EntityBuilder | RelationBuilder:
    return getattr(importlib.import_module('pypath.inputs_v2.' + module), name)


def _extract(
    module: str, mapper: str, row: dict, locator: str = 'reactions:1'
) -> SilverExtractor:
    result = SilverExtractor(module, mapper)
    result.process_record(_mapper(module, mapper)(row), row, locator, 1)
    return result


def _values(annotations: list[dict], term: str) -> list[str]:
    return [a['value'] for a in annotations if a['term'] == term]


@pytest.mark.parametrize('module,mapper,row', REACTION_CASES)
def test_direction_broadcasts_to_each_membership_without_changing_identity(
    module: str, mapper: str, row: dict
) -> None:
    """Broadcast source assertions without adding qualifiers or changing edge keys."""
    original = _extract(module, mapper, row)
    contextual = _extract(
        module,
        mapper,
        {
            **row,
            'direction': 'REVERSIBLE',
            'conversion_direction': 'UNSUPPORTED',
        },
    )
    baseline = [
        r
        for r in original.relations
        if r.predicate in {'has_input', 'has_output', 'enabled_by'}
    ]
    observations = [
        r
        for r in contextual.relations
        if r.predicate in {'has_input', 'has_output', 'enabled_by'}
    ]
    assert len(observations) >= 3
    assert [r.relation_key for r in observations] == [
        r.relation_key for r in baseline
    ]
    for relation in observations:
        assert set(_values(relation.annotations, CONVERSION_DIRECTION)) == {
            'REVERSIBLE',
            'UNSUPPORTED',
        }
        assert qualifiers(relation.annotations) == ()
        assert direction_sign(relation.annotations) == 0
        assert all(
            a['scope'] == 'relation'
            for a in relation.annotations
            if a['term'] == CONVERSION_DIRECTION
        )
    assert not any(
        _values(entity.annotations, CONVERSION_DIRECTION)
        for entity in contextual.entities.values()
    )
    assert all(
        not _values(r.annotations, CONVERSION_DIRECTION) for r in baseline
    )


@pytest.mark.parametrize('module,mapper,row', REACTION_CASES)
def test_direction_keeps_source_events_separate(
    module: str, mapper: str, row: dict
) -> None:
    """Keep directions attached to the row that asserted them."""
    first = _extract(
        module, mapper, dict(row, direction='LEFT-TO-RIGHT'), 'reactions:1'
    )
    second = _extract(
        module, mapper, dict(row, direction='REVERSIBLE'), 'reactions:2'
    )
    first_edges = [
        r
        for r in first.relations
        if r.predicate in {'has_input', 'has_output', 'enabled_by'}
    ]
    second_edges = [
        r
        for r in second.relations
        if r.predicate in {'has_input', 'has_output', 'enabled_by'}
    ]
    assert [r.relation_key for r in first_edges] == [
        r.relation_key for r in second_edges
    ]
    for edges, locator, expected in (
        (first_edges, 'reactions:1', 'LEFT-TO-RIGHT'),
        (second_edges, 'reactions:2', 'REVERSIBLE'),
    ):
        assert all(
            r.row_id == locator
            and _values(r.annotations, CONVERSION_DIRECTION) == [expected]
            for r in edges
        )


@pytest.mark.parametrize('missing', ['', 'invalid', 'bad-CHEBI:999-tail'])
@pytest.mark.parametrize(
    'mapper', ['reactions_schema', 'transport_reactions_schema']
)
def test_rhea_missing_identifier_preserves_compartment_pairing(
    missing: str, mapper: str
) -> None:
    """Retain the correct compartment after an omitted or invalid middle identity."""
    row = {
        'rhea_id': '123',
        'participant_chebi': f'1||{missing}||3',
        'participant_role': 'reactant||reactant||product',
        'participant_compartment': 'in||unknown||out',
        'participant_display_name': '|| ||',
    }
    result = _extract('rhea', mapper, row)
    edges = [
        r
        for r in result.relations
        if r.predicate in {'has_input', 'has_output'}
    ]
    assert [
        (
            r.predicate,
            result.entities[r.object_entity_key].identifier,
            _values(r.annotations, CELLULAR_LOCATION),
        )
        for r in edges
    ] == [('has_input', '1', ['in']), ('has_output', '3', ['out'])]
    assert not any(_values(r.annotations, CONVERSION_DIRECTION) for r in edges)


@pytest.mark.parametrize(
    'module,mapper,row',
    [
        case
        for case in REACTION_CASES
        if case[0] in {'recon3d', 'metatlas', 'kegg', 'reactome'}
    ],
)
def test_symbolic_coefficients_survive_source_context_addition(
    module: str, mapper: str, row: dict
) -> None:
    """Keep symbolic coefficient values alongside new source attributes."""
    result = _extract(module, mapper, dict(row, direction='UNSUPPORTED'))
    coefficients = {
        value
        for r in result.relations
        for value in _values(r.annotations, 'stoichiometry')
    }
    assert {'n', '(n+1)'} <= coefficients


@pytest.mark.parametrize(
    'subtype,expected', [('cancer', ['cancer']), ('', []), (None, [])]
)
def test_macdb_trait_type_is_an_ordinary_entity_attribute(
    subtype: str | None, expected: list[str]
) -> None:
    """Publish meaningful trait types and omit missing values."""
    result = _extract(
        'macdb',
        'trait_terms_schema',
        {
            'Trait_Ontology_ID': 'T1',
            'Trait_Ontology': 'Trait',
            'Trait_Type': subtype,
        },
    )
    assert len(result.entities) == 1
    assert (
        _values(next(iter(result.entities.values())).annotations, TRAIT_TYPE)
        == expected
    )


def _publish(
    path: Path, module: str, mapper: str, rows: list[dict]
) -> list[list[dict]]:
    writer = ParquetWriter(path)
    resolver = EntityResolver(
        library_dir=path / 'no-library', defer_aliases=True
    )
    try:
        for index, row in enumerate(rows):
            observations = _extract(
                module, mapper, row, f'interactions:{index}'
            )
            writer.append_observations(observations, resolver)
        paths = writer.close()[:3]
        return [pq.read_table(item).to_pylist() for item in paths]
    finally:
        resolver.close()


def test_ligand_receptor_roles_follow_canonical_flips_per_evidence(
    tmp_path: Path,
) -> None:
    """Preserve opposite roles on the same canonical edge as separate evidence."""
    rows = [
        {
            'Interaction ID': 'CDB1',
            'Species': 'human',
            'Ligand Symbols': 'ALPHA',
            'Receptor Symbols': 'BETA',
        },
        {
            'Interaction ID': 'CDB2',
            'Species': 'human',
            'Ligand Symbols': 'BETA',
            'Receptor Symbols': 'ALPHA',
        },
    ]
    entities, relations, _ = _publish(
        tmp_path, 'connectomedb', 'interactions_schema', rows
    )
    assert len(relations) == 1
    relation = relations[0]
    assert relation['sign'] == 0 and relation['is_directed'] is False
    assert relation['evidence_count'] == 2
    assert relation['relation_key'] == relation_key(
        relation['subject_entity_key'],
        'interacts_with',
        relation['object_entity_key'],
    )
    by_key = {e['entity_key']: e for e in entities}
    for evidence in relation['evidence']:
        index = int(evidence['row_id'].rsplit(':', 1)[1])
        role_values = {}
        for annotation in evidence['annotations']:
            if annotation['term'] != PARTICIPANT_ROLE:
                continue
            endpoint_key = relation[f'{annotation["scope"]}_entity_key']
            symbols = {
                i['id']
                for i in by_key[endpoint_key]['identifiers']
                if i['ns'] == 'genesymbol'
            }
            expected = rows[index][
                annotation['value'].capitalize() + ' Symbols'
            ]
            assert expected in symbols
            role_values[annotation['value']] = endpoint_key
        assert set(role_values) == {'ligand', 'receptor'}
        assert role_values['ligand'] != role_values['receptor']
    assert not any(
        _values(e['annotations'], PARTICIPANT_ROLE) for e in entities
    )


@pytest.mark.parametrize(
    'mapper', ['reactions_schema', 'transport_reactions_schema']
)
def test_rhea_short_compartment_array_does_not_broadcast(mapper: str) -> None:
    """Leave unreported compartments absent instead of copying the first value."""
    result = _extract(
        'rhea',
        mapper,
        {
            'rhea_id': '123',
            'participant_chebi': '1||2||3',
            'participant_role': 'reactant||reactant||product',
            'participant_compartment': 'in',
        },
    )
    edges = [
        r
        for r in result.relations
        if r.predicate in {'has_input', 'has_output'}
    ]
    assert [_values(r.annotations, CELLULAR_LOCATION) for r in edges] == [
        ['in'],
        [],
        [],
    ]


@pytest.mark.parametrize(
    'mapper', ['reactions_schema', 'transport_reactions_schema']
)
def test_rhea_short_role_array_does_not_invent_roles(mapper: str) -> None:
    """Do not assign an unreported role to later source participants."""
    result = _extract(
        'rhea',
        mapper,
        {
            'rhea_id': '123',
            'participant_chebi': '1||2||3',
            'participant_role': 'reactant',
            'participant_compartment': 'in||unknown||out',
        },
    )
    edges = [
        r
        for r in result.relations
        if r.predicate in {'has_input', 'has_output'}
    ]
    assert (
        len(edges) == 1
        and result.entities[edges[0].object_entity_key].identifier == '1'
    )


@pytest.mark.parametrize(
    'mapper', ['reactions_schema', 'transport_reactions_schema']
)
def test_rhea_short_identifier_array_does_not_duplicate_identity(
    mapper: str,
) -> None:
    """Do not repeat a source identity at unrelated participant positions."""
    result = _extract(
        'rhea',
        mapper,
        {
            'rhea_id': '123',
            'participant_chebi': '1',
            'participant_role': 'reactant||reactant||product',
            'participant_compartment': 'in||unknown||out',
        },
    )
    edges = [
        r
        for r in result.relations
        if r.predicate in {'has_input', 'has_output'}
    ]
    assert (
        len(edges) == 1
        and result.entities[edges[0].object_entity_key].identifier == '1'
    )


@pytest.mark.parametrize('module', ['rhea', 'recon3d'])
def test_direction_events_survive_parquet_aggregation(
    tmp_path: Path, module: str
) -> None:
    """Keep each row's direction when identical participant edges aggregate."""
    row = next(
        row
        for source, mapper, row in REACTION_CASES
        if source == module and mapper == 'reactions_schema'
    )
    rows = [
        dict(row, direction=value) for value in ('LEFT-TO-RIGHT', 'REVERSIBLE')
    ]
    entities, relations, _ = _publish(
        tmp_path, module, 'reactions_schema', rows
    )
    direct = [
        r
        for r in relations
        if r['predicate'] in {'has_input', 'has_output', 'enabled_by'}
    ]
    assert len(direct) >= 3
    for relation in direct:
        assert relation['evidence_count'] == 2
        assert relation['sign'] == 0
        assert relation['relation_key'] == relation_key(
            relation['subject_entity_key'],
            relation['predicate'],
            relation['object_entity_key'],
        )
        for evidence in relation['evidence']:
            index = int(evidence['row_id'].rsplit(':', 1)[1])
            assert _values(evidence['annotations'], CONVERSION_DIRECTION) == [
                rows[index]['direction']
            ]
            attributes = [
                a
                for a in evidence['annotations']
                if a['term'] == CONVERSION_DIRECTION
            ]
            assert all(
                a['scope'] == 'relation' and a['source'] == module
                for a in attributes
            )
    assert not any(
        _values(e['annotations'], CONVERSION_DIRECTION) for e in entities
    )
