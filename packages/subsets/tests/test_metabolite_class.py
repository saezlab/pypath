"""
The metabolite class product over a small published release (``metabolite_class_fixture``).

One test per rule of the plan. The rule file of most tests enables all three signals: the
cascade code is the same in both stages, so signals 2 and 3 are tested before they ship.
"""

import json
import os
from pathlib import Path

import duckdb
import pytest

from omnipath_core import subset_layout as layout
from omnipath_subsets.checks import check_metabolite_class, check_release
from omnipath_subsets.metabolite_class import RuleError, read_rules, write_metabolite_class
from omnipath_subsets.parquet_products import build_release_subsets
from metabolite_class_fixture import (
    CORTISOL_17S_KEY,
    CORTISOL_KEY,
    write_release,
    write_rules,
)

GOLDEN = Path(__file__).parent / 'golden' / 'metabolite_class.json'


def build(tmp_path, rules_text = None):
    """
    Build the product over the fixture release; returns its directory and manifest.
    """

    root = tmp_path / 'data'
    version = write_release(root)
    rules = write_rules(tmp_path, rules_text) if rules_text else write_rules(tmp_path)
    manifest = layout.build_product(
        root,
        version,
        'metabolite_class',
        lambda directory, context: write_metabolite_class(
            directory,
            context,
            rules_path = rules,
        ),
    )

    return layout.subset_directory(root, version, 'metabolite_class'), manifest


def rows(directory, table = 'metabolite_class'):
    """
    The rows of one table of the product, as dicts.
    """

    con = duckdb.connect()
    cursor = con.execute(f"SELECT * FROM read_parquet('{directory / (table + '.parquet')}')")
    names = [d[0] for d in cursor.description]

    return [dict(zip(names, row)) for row in cursor.fetchall()]


def by_entity(directory):
    """
    The rows of ``metabolite_class`` by entity key.
    """

    return {row['entity_key']: row for row in rows(directory)}


@pytest.fixture
def built(tmp_path):

    return build(tmp_path)


@pytest.fixture
def table(built):

    return by_entity(built[0])


# The cascade.


def test_each_signal_sets_its_class_and_an_unclaimed_compound_has_no_row(table):

    assert table['cortisol']['metabolite_class'] == 'hmdb'
    assert table['menthol']['metabolite_class'] == 'chebi_role'
    assert table['pyruvate']['metabolite_class'] == 'metabolic_network'
    assert 'screening_compound' not in table
    assert {
        'receptor',
        'rhea_reaction',
        'metatlas_reaction',
        'hmdb_disease',
    }.isdisjoint(table)


def test_the_first_signal_wins_and_every_asserting_resource_is_named(table):

    cortisol = table['cortisol']

    assert (cortisol['metabolite_class'], cortisol['match_kind']) == ('hmdb', 'entity')
    assert cortisol['supported_by'] == ['chebi', 'hmdb']


def test_a_cross_reference_to_hmdb_from_another_resource_does_not_make_an_hmdb_molecule(table):

    assert 'xref_compound' not in table


def test_the_role_terms_themselves_are_not_carriers(table):

    assert 'role_metabolite' not in table
    assert 'role_human_metabolite' not in table


# Connectivity groups.


def test_a_stereoisomer_takes_the_class_of_its_group_through_the_member_that_gives_it(table):

    variant = table['cortisol_17s']

    assert variant['metabolite_class'] == 'hmdb'
    assert variant['match_kind'] == 'structure_key'
    assert variant['via_entity_key'] == 'cortisol'
    assert variant['supported_by'] == []


def test_a_molecule_its_group_moves_to_a_better_class_keeps_its_own_resources(table):

    ion = table['lactate_ion']

    assert (ion['metabolite_class'], ion['match_kind']) == ('hmdb', 'structure_key')
    assert ion['via_entity_key'] == 'lactate'
    assert ion['supported_by'] == ['rhea']


def test_an_entity_row_names_no_member(table):

    assert all(
        (row['via_entity_key'] is None) == (row['match_kind'] == 'entity')
        for row in table.values()
    )


def test_a_fallback_and_a_lipid_keep_their_own_class_without_an_inchikey(table):

    fallback, lipid = table['fallback'], table['lipid']

    assert (fallback['metabolite_class'], fallback['match_kind']) == (
        'metabolic_network',
        'entity',
    )
    assert (lipid['metabolite_class'], lipid['match_kind']) == ('hmdb', 'entity')
    assert fallback['inchikey'] is None
    assert lipid['inchikey'] is None


def test_no_entity_is_merged_with_another(table):

    assert len({row['entity_key'] for row in table.values()}) == len(table) == 8


# The readable columns.


def test_the_label_is_the_hmdb_name_then_the_chebi_name_then_the_release_label(table):

    assert table['cortisol']['label'] == 'Cortisol'
    assert table['menthol']['label'] == 'menthol'
    assert table['cortisol_17s']['label'] == 'Cortisol (17S)'


def test_the_inchikey_is_the_entity_identifier_and_hmdb_the_lowest_seven_digit_accession(table):

    assert table['cortisol']['inchikey'] == CORTISOL_KEY
    assert table['cortisol_17s']['inchikey'] == CORTISOL_17S_KEY
    assert table['cortisol']['hmdb'] == 'HMDB0000063'
    assert table['lactate']['hmdb'] == 'HMDB0000190'
    assert table['pyruvate']['hmdb'] is None


def test_the_hmdb_accession_is_the_record_s_own_and_not_a_lower_one_the_resolver_added(table):

    # The HMDB record of cortisol also carries HMDB0000001 from another record.
    assert table['cortisol']['hmdb'] == 'HMDB0000063'


# The report.


def report(directory, kind):
    """
    The report rows of one kind, by class.
    """

    return {
        row['metabolite_class']: row
        for row in rows(directory, 'metabolite_class_report')
        if row['kind'] == kind
    }


def test_the_report_gives_the_size_of_every_class(built):

    sizes = {name: row['n'] for name, row in report(built[0], 'class_size').items()}

    assert sizes == {'hmdb': 5, 'chebi_role': 1, 'metabolic_network': 2}


def test_the_report_counts_each_combination_of_signals(built):

    overlap = {
        tuple(json.loads(row['detail'])['classes']): row['n']
        for row in rows(built[0], 'metabolite_class_report')
        if row['kind'] == 'overlap'
    }

    assert overlap == {
        ('chebi_role', 'hmdb'): 1,
        ('hmdb',): 2,
        ('chebi_role',): 1,
        ('metabolic_network',): 3,
    }


def test_the_report_lists_the_groups_that_pulled_in_variants(built):

    shared = report(built[0], 'shared_key')
    keys = {k['key']: k for k in json.loads(shared['hmdb']['detail'])['keys']}

    assert {name: row['n'] for name, row in shared.items()} == {
        'hmdb': 2,
        'chebi_role': 0,
        'metabolic_network': 0,
    }
    assert keys[CORTISOL_KEY[:14]]['via'] == 'cortisol'
    assert keys[CORTISOL_KEY[:14]]['labels'] == ['Cortisol (17S)']


def test_the_report_counts_the_claimed_molecules_without_a_structure_key(built):

    unkeyed = report(built[0], 'no_key')

    assert {name: row['n'] for name, row in unkeyed.items()} == {
        'hmdb': 1,
        'chebi_role': 0,
        'metabolic_network': 1,
    }
    assert json.loads(unkeyed['hmdb']['detail'])['labels'] == ['PC(34:1)']


# The rule file.


def test_the_packaged_rule_file_claims_hmdb_alone_and_a_commented_signal_adds_nothing(tmp_path):

    root = tmp_path / 'data'
    version = write_release(root)
    layout.build_product(root, version, 'metabolite_class', write_metabolite_class)
    table = by_entity(layout.subset_directory(root, version, 'metabolite_class'))

    assert {row['metabolite_class'] for row in table.values()} == {'hmdb'}
    assert {'menthol', 'pyruvate', 'fallback'}.isdisjoint(table)
    # Without the reaction signal the charged lactate is an unclaimed variant.
    assert table['lactate_ion']['supported_by'] == []
    assert [s.signal for s in read_rules().sources] == ['resource']


BROKEN = {
    'unknown class': 'sources:\n  - {resource: hmdb, priority: 1, class: hmbd}\n',
    'default class': 'sources:\n  - {resource: hmdb, priority: 1, class: none}\n',
    'two signals': 'sources:\n  - {resource: hmdb, chebi_role_root: X, priority: 1, class: hmdb}\n',
    'unknown key': 'sources:\n  - {resource: hmdb, priority: 1, class: hmdb, taxon: 9606}\n',
    'shared priority': (
        'sources:\n  - {resource: hmdb, priority: 1, class: hmdb}\n'
        '  - {resource: chebi, priority: 1, class: chebi_role}\n'
    ),
    'worse class first': (
        'sources:\n  - {resource: rhea, priority: 1, class: metabolic_network}\n'
        '  - {resource: hmdb, priority: 2, class: hmdb}\n'
    ),
    'network without resources': (
        'sources:\n  - {predicates: [has_input], priority: 1, class: metabolic_network}\n'
    ),
}
CLASSES = (
    'classes:\n  - {name: hmdb, precedence: 1}\n  - {name: chebi_role, precedence: 2}\n'
    '  - {name: metabolic_network, precedence: 3}\n  - {name: none, precedence: 4}\n'
    'default: none\n'
)


@pytest.mark.parametrize('case', sorted(BROKEN))
def test_a_broken_rule_file_is_refused_before_the_build_reads_the_release(tmp_path, case):

    with pytest.raises(RuleError):

        read_rules(write_rules(tmp_path, CLASSES + BROKEN[case]))


# The product.


def test_two_builds_give_identical_tables(tmp_path):

    first, _ = build(tmp_path / 'one')
    second, _ = build(tmp_path / 'two')

    for table in ('metabolite_class', 'metabolite_class_report'):

        assert rows(first, table) == rows(second, table)


def test_every_row_carries_the_release_digest(built):

    directory, manifest = built
    digest = manifest['release_manifest_sha256']

    for table in ('metabolite_class', 'metabolite_class_report'):

        assert {row['build_id'] for row in rows(directory, table)} == {digest}

    assert manifest['row_contract'] == layout.ROW_CONTRACTS['metabolite_class']
    assert manifest['stats']['classes'] == {
        'hmdb': 5,
        'chebi_role': 1,
        'metabolic_network': 2,
    }


def test_the_release_build_writes_the_product(tmp_path):

    root = tmp_path / 'data'
    version = write_release(root)
    manifests = build_release_subsets(root, version, ['metabolite_class'])

    assert set(manifests) == {'metabolite_class'}
    assert layout.read_subset(root, version, 'metabolite_class') is not None


def test_the_checks_pass_on_a_built_product_and_find_a_broken_one(tmp_path):

    root = tmp_path / 'data'
    version = write_release(root)
    rules = write_rules(tmp_path)
    layout.build_product(
        root,
        version,
        'metabolite_class',
        lambda d, c: write_metabolite_class(d, c, rules_path = rules),
    )
    directory = layout.subset_directory(root, version, 'metabolite_class')

    assert check_release(root, version, ('metabolite_class',)) == {'metabolite_class': []}

    # The ion names a member of another class, menthol gains a row of an unknown class, and
    # cortisol is listed twice, so the size of hmdb in the report is wrong.
    broken = tmp_path / 'broken.parquet'
    members = directory / 'metabolite_class.parquet'
    con = duckdb.connect()
    con.execute(f"CREATE VIEW t AS SELECT * FROM read_parquet('{members}')")
    con.execute(
        "COPY (SELECT * REPLACE (CASE WHEN entity_key = 'lactate_ion' THEN 'menthol' "
        'ELSE via_entity_key END AS via_entity_key) FROM t '
        "UNION ALL SELECT * REPLACE ('nonsense' AS metabolite_class) FROM t "
        "WHERE entity_key = 'menthol' "
        f"UNION ALL SELECT * FROM t WHERE entity_key = 'cortisol') TO '{broken}' "
        '(FORMAT parquet)'
    )
    found = {
        f.check
        for f in check_metabolite_class(broken, directory / 'metabolite_class_report.parquet')
    }

    assert {
        'classes',
        'no_duplicates',
        'structure_key_members',
        'report_classes',
    } <= found


# The frozen rows.


def _frozen(directory):
    """
    Both tables of the product without the release digest, which changes with the release.
    """

    def plain(row):

        return {k: v for k, v in row.items() if k != 'build_id'}

    return {
        'metabolite_class': [plain(r) for r in rows(directory)],
        'metabolite_class_report': [
            plain(r)
            for r in rows(directory, 'metabolite_class_report')
        ],
    }


def test_the_product_equals_the_frozen_rows(built):

    found = _frozen(built[0])

    if os.environ.get('UPDATE_GOLDEN'):

        GOLDEN.write_text(json.dumps(found, indent = 1, sort_keys = True) + '\n')

    assert found == json.loads(GOLDEN.read_text())
