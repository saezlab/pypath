"""
The ChEBI input keeps the role branch, and leaves application out.

The tests parse a small synthetic OBO text: the container holds no ChEBI
download. Each test states one rule.
"""

import io

import pytest

from omnipath_build.silver import SilverExtractor
from omnipath_resources.inputs_v2 import chebi
from omnipath_resources.inputs_v2.parsers import chebi as chebi_parser


_OBO = '''format-version: 1.2
ontology: chebi

[Term]
id: CHEBI:24431
name: chemical entity

[Term]
id: CHEBI:23367
name: molecular entity
is_a: CHEBI:24431 ! chemical entity

[Term]
id: CHEBI:17234
name: glucose
is_a: CHEBI:23367 ! molecular entity
property_value: chemrof:inchi_string "InChI=1S/C6H12O6/c7-1-2-3(9)4(10)5(11)6(12)8/h1,3-6,8-12H,2H2" xsd:string
property_value: chemrof:inchi_key_string "GZCGUPFRVQAUEE-UHFFFAOYSA-N" xsd:string
property_value: chemrof:generalized_empirical_formula "C6H12O6" xsd:string
relationship: RO:0000087 CHEBI:77746 ! human metabolite

[Term]
id: CHEBI:15377
name: water
is_a: CHEBI:23367 ! molecular entity
property_value: chemrof:inchi_key_string "XLYOFNOQVPJJNP-UHFFFAOYSA-N" xsd:string
relationship: RO:0000087 CHEBI:25212 ! metabolite
relationship: RO:0000087 CHEBI:23888 ! drug

[Term]
id: CHEBI:50906
name: role

[Term]
id: CHEBI:24432
name: biological role
is_a: CHEBI:50906 ! role

[Term]
id: CHEBI:25212
name: metabolite
is_a: CHEBI:24432 ! biological role

[Term]
id: CHEBI:77746
name: human metabolite
is_a: CHEBI:25212 ! metabolite

[Term]
id: CHEBI:33232
name: application
is_a: CHEBI:50906 ! role

[Term]
id: CHEBI:23888
name: drug
is_a: CHEBI:33232 ! application

[Term]
id: CHEBI:90001
name: drug metabolite
is_a: CHEBI:25212 ! metabolite
is_a: CHEBI:23888 ! drug

[Term]
id: CHEBI:90002
name: role with a structure
is_a: CHEBI:24432 ! biological role
property_value: chemrof:inchi_key_string "AAAAAAAAAAAAAA-UHFFFAOYSA-N" xsd:string
property_value: chemrof:smiles_string "C" xsd:string

[Term]
id: CHEBI:90003
name: obsolete role
is_a: CHEBI:24432 ! biological role
is_obsolete: true
'''

_ROLE_IDS = {
    'CHEBI:50906',
    'CHEBI:24432',
    'CHEBI:25212',
    'CHEBI:77746',
    'CHEBI:90001',
    'CHEBI:90002',
}
_APPLICATION_IDS = {'CHEBI:33232', 'CHEBI:23888'}

# The fields of the molecule rows that the role branch could change, as the
# parser built them before it kept roles. The full rows were checked equal.
_MOLECULES_BEFORE = {
    'CHEBI:15377': {
        'ancestor_terms': ['CHEBI:23367', 'CHEBI:24431'],
        'is_a': ['CHEBI:23367'],
        'relationships': [
            {'type': 'RO:0000087', 'target': 'CHEBI:25212', 'target_name': 'metabolite'},
            {'type': 'RO:0000087', 'target': 'CHEBI:23888', 'target_name': 'drug'},
        ],
        'inchi': '',
        'inchikey': 'XLYOFNOQVPJJNP-UHFFFAOYSA-N',
        'formula': '',
    },
    'CHEBI:17234': {
        'ancestor_terms': ['CHEBI:23367', 'CHEBI:24431'],
        'is_a': ['CHEBI:23367'],
        'relationships': [
            {'type': 'RO:0000087', 'target': 'CHEBI:77746', 'target_name': 'human metabolite'},
        ],
        'inchi': 'InChI=1S/C6H12O6/c7-1-2-3(9)4(10)5(11)6(12)8/h1,3-6,8-12H,2H2',
        'inchikey': 'GZCGUPFRVQAUEE-UHFFFAOYSA-N',
        'formula': 'C6H12O6',
    },
    'CHEBI:23367': {
        'ancestor_terms': ['CHEBI:24431'],
        'is_a': ['CHEBI:24431'],
        'relationships': [],
        'inchi': '',
        'inchikey': '',
        'formula': '',
    },
    'CHEBI:24431': {
        'ancestor_terms': [],
        'is_a': [],
        'relationships': [],
        'inchi': '',
        'inchikey': '',
        'formula': '',
    },
}


@pytest.fixture(scope = 'module')
def records():

    terms = chebi_parser._parse_terms(io.StringIO(_OBO))

    return chebi_parser._build_records(terms)


@pytest.fixture(scope = 'module')
def rows(records):

    return {row['chebi_id']: row for row in records['molecules']}


@pytest.fixture(scope = 'module')
def silver(rows):

    extractor = SilverExtractor('chebi', 'molecules')

    for row_number, row in enumerate(rows.values(), start = 1):

        record = chebi.molecules_schema(row)
        extractor.process_record(record, row, row['chebi_id'], row_number)

    return extractor


def _entity(silver, entity_key):

    entity = silver.entities[entity_key]

    return entity.entity_type, entity.namespace, entity.identifier


def _edges(silver, predicate):

    return {
        (
            _entity(silver, relation.subject_entity_key)[2],
            _entity(silver, relation.object_entity_key)[2],
        )
        for relation in silver.relations
        if relation.predicate == predicate
    }


def _chebi(term_id):

    return term_id.removeprefix('CHEBI:')


def test_roles_become_chemical_entities(rows, silver):

    assert _ROLE_IDS <= set(rows)

    entities = {
        entity.identifier: entity.entity_type
        for entity in silver.entities.values()
        if entity.evidence
    }

    for role_id in ('CHEBI:25212', 'CHEBI:77746'):

        assert entities[_chebi(role_id)] == 'chemical_entity'


def test_sub_role_is_subclass_of_role(rows, silver):

    assert rows['CHEBI:77746']['is_a'] == ['CHEBI:25212']
    assert rows['CHEBI:25212']['is_a'] == ['CHEBI:24432']

    relations = [
        relation
        for relation in silver.relations
        if relation.predicate == 'subclass_of'
    ]
    sub_role = [
        relation
        for relation in relations
        if _entity(silver, relation.subject_entity_key)[2] == '77746'
    ]

    assert [_entity(silver, r.object_entity_key) for r in sub_role] == [
        ('chemical_entity', 'chebi', '25212'),
    ]
    assert {r.statement_kind for r in sub_role} == {'ontology'}
    assert {r.source for r in sub_role} == {'chebi'}
    assert {('25212', '24432'), ('24432', '50906')} <= _edges(silver, 'subclass_of')


def test_has_role_points_to_the_sub_role_entity(silver):

    roles = _edges(silver, 'has_chemical_role')

    assert ('17234', '77746') in roles
    assert ('15377', '25212') in roles

    targets = {
        _entity(silver, relation.object_entity_key)
        for relation in silver.relations
        if relation.predicate == 'has_chemical_role'
    }

    assert ('chemical_entity', 'chebi', '77746') in targets


def test_role_never_becomes_a_structured_molecule(records, rows, silver):

    role = rows['CHEBI:90002']

    assert role['inchikey'] == role['inchi'] == role['smiles'] == ''

    ontology_ids = {term['id'] for term in records['ontology_terms']}

    assert _ROLE_IDS <= ontology_ids

    role_entities = [
        entity
        for entity in silver.entities.values()
        if entity.evidence and entity.identifier in {_chebi(i) for i in _ROLE_IDS}
    ]

    assert len(role_entities) == len(_ROLE_IDS)

    for entity in role_entities:

        assert entity.namespace == 'chebi'
        assert not {
            identifier['ns']
            for identifier in entity.identifiers
        } & {'inchikey', 'inchi', 'smiles'}

    assert not [
        id_row
        for role_id in _ROLE_IDS
        for id_row in chebi._id_translation_rows(rows[role_id])
    ]


def test_application_branch_stays_excluded(rows, records):

    assert not _APPLICATION_IDS & set(rows)
    assert not _APPLICATION_IDS & {term['id'] for term in records['ontology_terms']}
    assert 'CHEBI:90003' not in rows


def test_role_with_a_parent_below_application_keeps_its_role_parent(rows):

    assert rows['CHEBI:90001']['is_a'] == ['CHEBI:25212']


def test_molecule_records_are_unchanged(rows):

    molecules = {
        term_id: {key: rows[term_id][key] for key in expected}
        for term_id, expected in _MOLECULES_BEFORE.items()
    }

    assert molecules == _MOLECULES_BEFORE
    assert set(rows) == set(_MOLECULES_BEFORE) | _ROLE_IDS
