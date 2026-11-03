"""A malformed chemical identifier value must never become evidence at all
(spec 011 T043, spec.md US1 acceptance scenario 4).

`InChIKey=none`, `none`, `null`, `na` and similar placeholder strings are
values the live build database's own `entity` table already carries as
Standard InChIKey canonical identities today -- confirmed by
`tests/test_identity_validity.py`'s live-database check. That happens
because a malformed value, once written as evidence, can still become an
entity's own asserted canonical identity through the direct-assertion
path, not only through the resolver join T040/T041 fixed. The fix has to
sit before the value is ever written, not after.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import SimpleNamespace

import pytest

from omnipath_build.evidence_projector import (
    EvidenceProjectorBase,
    _MutableProjectionStats,
)


@dataclass
class _FakeWriter:
    rows: list[dict] = field(default_factory=list)

    def write(self, row: dict) -> None:
        self.rows.append(row)


@dataclass
class _FakeWriters:
    entity: _FakeWriter = field(default_factory=_FakeWriter)
    identifier: _FakeWriter = field(default_factory=_FakeWriter)
    entity_annotation: _FakeWriter = field(default_factory=_FakeWriter)
    relation_annotation: _FakeWriter = field(default_factory=_FakeWriter)
    annotation: _FakeWriter = field(default_factory=_FakeWriter)
    relation: _FakeWriter = field(default_factory=_FakeWriter)
    annotation_relation: _FakeWriter = field(default_factory=_FakeWriter)
    ontology_relation: _FakeWriter = field(default_factory=_FakeWriter)


def _project_one(*identifiers):
    entity = SimpleNamespace(
        type='chemical',
        identifiers=[
            SimpleNamespace(type=t, value=v) for t, v in identifiers
        ],
        annotations=[],
        associations=[],
        membership=[],
    )
    writers = _FakeWriters()
    EvidenceProjectorBase()._flatten_entity_tree(
        entity,
        source='s',
        dataset='d',
        row_id=1,
        occurrence_id='occ',
        parent_entity_evidence_id=None,
        entity_role='representative',
        writers=writers,
        seen_annotations=set(),
        stats=_MutableProjectionStats(),
    )
    return writers.identifier.rows


@pytest.mark.parametrize(
    'garbage_value',
    ['InChIKey=none', 'none', 'null', 'na', '-'],
)
def test_malformed_inchikey_value_is_never_written(garbage_value):
    rows = _project_one(('Standard Inchi Key:MI:1101', garbage_value))
    assert not rows, (
        f'{garbage_value!r} was written as evidence -- a malformed chemical '
        'value must never reach identifier_evidence at all'
    )


def test_a_valid_inchikey_is_still_written():
    rows = _project_one(
        ('Standard Inchi Key:MI:1101', 'BQJCRHHNABKAKN-KBQPJGBKSA-N'),
    )
    assert rows and rows[0]['identifier'] == 'BQJCRHHNABKAKN-KBQPJGBKSA-N'


def test_a_non_chemical_identifier_is_unaffected():
    """The reject-on-malformed rule is chemical-specific -- a non-chemical
    identifier type is never subject to it, whatever its raw form."""
    rows = _project_one(('Uniprot:MI:0486', 'not-a-real-accession'))
    assert rows and rows[0]['identifier'] == 'not-a-real-accession'
