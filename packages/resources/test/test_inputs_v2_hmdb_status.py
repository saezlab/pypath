"""
HMDB's status of a metabolite, read from the XML and emitted as an annotation.

The fixture is a small synthetic file in HMDB's schema: the container cannot
download HMDB.
"""

import io
from types import SimpleNamespace

import pytest

from omnipath_build.silver import SilverExtractor
from omnipath_resources.inputs_v2 import hmdb
from omnipath_resources.inputs_v2.parsers.hmdb import _raw


STATUSES = ('quantified', 'detected', 'expected', 'predicted')


def _metabolite(accession, status = None, extra = ''):

    status_element = '' if status is None else f'<status>{status}</status>'

    return f"""
  <metabolite>
    <version>5.0</version>
    {status_element}
    <accession>{accession}</accession>
    <name>Metabolite {accession}</name>
    <description>A metabolite.</description>
    <inchikey>BRMWTNUJHUMWMS-LURJTMIESA-N</inchikey>
    <kegg_id>C01152</kegg_id>
    <synonyms>
      <synonym>Synonym one</synonym>
      <synonym>Synonym two</synonym>
    </synonyms>
    <biological_properties>
      <cellular_locations>
        <cellular>Cytoplasm</cellular>
      </cellular_locations>
    </biological_properties>
    {extra}
  </metabolite>"""


def _opener(*metabolites):

    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<hmdb xmlns="http://www.hmdb.ca">' +
        ''.join(metabolites) +
        '\n</hmdb>\n'
    )

    return SimpleNamespace(
        result = {'hmdb_metabolites.xml': io.BytesIO(xml.encode('utf-8'))},
    )


def _rows(*metabolites):

    return list(_raw(_opener(*metabolites)))


def _status_annotations(row):

    record = hmdb.metabolites_schema(row)
    extractor = SilverExtractor(hmdb.config.id.value, 'fixture')
    extractor.process_record(record, row, 'fixture:1', 1)

    return [
        a['value']
        for entity in extractor.entities.values()
        for a in (entity.annotations or [])
        if a['term'] == hmdb.STATUS
    ]


@pytest.mark.parametrize('status', STATUSES)
def test_each_status_is_read(status):

    rows = _rows(_metabolite('HMDB0000001', status))

    assert [row['status'] for row in rows] == [status]


@pytest.mark.parametrize('status', STATUSES)
def test_each_status_is_emitted_as_the_annotation(status):

    row, = _rows(_metabolite('HMDB0000001', status))

    assert _status_annotations(row) == [status]


def test_the_annotation_term_is_hmdb_status():

    assert hmdb.STATUS == 'hmdb:status'


def test_no_status_element_gives_no_status_annotation():

    row, = _rows(_metabolite('HMDB0000001'))

    assert row['status'] == ''
    assert _status_annotations(row) == []


def test_an_empty_status_element_gives_no_status_annotation():

    row, = _rows(_metabolite('HMDB0000001', '  '))

    assert _status_annotations(row) == []


@pytest.mark.parametrize(
    'written, emitted',
    [
        ('  Quantified\n', 'quantified'),
        ('DETECTED', 'detected'),
        ('\n    Expected  ', 'expected'),
        ('Detected  and\n Quantified', 'detected and quantified'),
    ],
)
def test_whitespace_and_case_are_normalized(written, emitted):

    row, = _rows(_metabolite('HMDB0000001', written))

    assert _status_annotations(row) == [emitted]


def test_a_nested_status_is_not_the_metabolite_status():

    nested = (
        '<protein_associations><protein>'
        '<status>predicted</status>'
        '</protein></protein_associations>'
    )
    row, = _rows(_metabolite('HMDB0000001', extra = nested))

    assert _status_annotations(row) == []


def test_each_metabolite_keeps_its_own_status():

    rows = _rows(
        _metabolite('HMDB0000001', 'quantified'),
        _metabolite('HMDB0000002'),
        _metabolite('HMDB0000003', 'predicted'),
    )

    assert [(row['accession'], row['status']) for row in rows] == [
        ('HMDB0000001', 'quantified'),
        ('HMDB0000002', ''),
        ('HMDB0000003', 'predicted'),
    ]


def test_the_other_fields_of_the_record_are_unchanged():

    with_status, = _rows(_metabolite('HMDB0000001', 'quantified'))
    without_status, = _rows(_metabolite('HMDB0000001'))

    assert with_status.pop('status') == 'quantified'
    assert without_status.pop('status') == ''
    assert with_status == without_status
    assert with_status['accession'] == 'HMDB0000001'
    assert with_status['name'] == 'Metabolite HMDB0000001'
    assert with_status['synonyms'] == 'Synonym one;Synonym two'
    assert with_status['cellular_locations'] == ['Cytoplasm']


def test_the_other_annotations_of_the_entity_are_unchanged():

    with_status, = _rows(_metabolite('HMDB0000001', 'quantified'))
    without_status, = _rows(_metabolite('HMDB0000001'))

    def other_annotations(row):

        record = hmdb.metabolites_schema(row)

        return [
            (str(a.term), a.value)
            for a in record.annotations
            if str(a.term) != hmdb.STATUS
        ]

    assert other_annotations(with_status) == other_annotations(without_status)
    assert other_annotations(with_status)
