"""The per-resource MetSigDB mapping rules (cycle 010).

Pure unit tests. The rules are a frozen contract, so they are checked without a
database: `contracts/mapping-rules.md` is the reference.
"""

from __future__ import annotations

import pytest

from omnipath_build.metsigdb.mapping import RESOURCES, rule_for

# The frozen per-resource set_type matrix. A resource contributes exactly one
# semantic, and no row's semantics come from source-side free text.
FROZEN_SET_TYPES = {
    'KEGG': 'pathway',
    'Reactome': 'pathway',
    'WikiPathways': 'pathway',
    'MACdb': 'disease',
    'ClassyFire': 'chemical_class',
}

# The `data_source.name` each rule reads from.
FROZEN_SOURCES = {
    'KEGG': 'kegg',
    'Reactome': 'reactome',
    'WikiPathways': 'wikipathways',
    'MACdb': 'macdb',
    'ClassyFire': 'hmdb',
}


def test_every_contract_resource_has_a_rule():
    assert {rule.name for rule in RESOURCES} == set(FROZEN_SET_TYPES)


def test_set_type_matches_the_frozen_matrix():
    assert {rule.name: rule.set_type for rule in RESOURCES} == FROZEN_SET_TYPES


def test_source_name_matches_the_frozen_matrix():
    assert {rule.name: rule.source_name for rule in RESOURCES} == FROZEN_SOURCES


def test_rules_are_unique():
    names = [rule.name for rule in RESOURCES]
    assert len(names) == len(set(names))


def test_classyfire_declares_its_hierarchy_source():
    """ClassyFire is two sources: HMDB assigns, ChemOnt supplies the hierarchy."""
    classyfire = rule_for('ClassyFire')
    assert classyfire.hierarchy_source_name == 'chemont'
    assert all(
        rule.hierarchy_source_name is None
        for rule in RESOURCES
        if rule.name != 'ClassyFire'
    )


def test_no_rule_derives_an_organism_from_an_identifier():
    """The organism is read from the source, never guessed from the id.

    Deriving 9606 from an `R-HSA-` prefix gave the right answer for Reactome
    and null for the 38 species WikiPathways covers, because a WP id names no
    species while the source record does.
    """
    assert not any(hasattr(rule, 'organism_sql') for rule in RESOURCES)


def test_rule_for_rejects_an_unknown_resource():
    with pytest.raises(KeyError):
        rule_for('SMPDB')


def test_no_sql_file_carries_a_bare_percent_sign():
    """psycopg2 interpolates these files, comments included.

    A bare `%` reads as a malformed placeholder and fails the whole load with
    "dict is not a sequence", which names neither the file nor the character.
    The first version of this test covered the extraction files only, and the
    fourth occurrence of the bug was in the shared publication file.
    """
    import re
    from pathlib import Path

    import omnipath_build.metsigdb.build as build

    named = re.compile(r'%\([a-z_]+\)s')
    offenders = []
    for path in sorted(Path(build._SQL_DIR).glob('*.sql')):
        if '%' in named.sub('', path.read_text(encoding='utf-8')):
            offenders.append(path.name)
    assert not offenders, f'bare percent sign in: {", ".join(offenders)}'
