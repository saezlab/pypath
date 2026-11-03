"""Source-specific Biolink contract for the CollecTRI input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import collectri as mod
from pypath.inputs_v2.parsers.base import iter_csv

# Lines copied from CollecTRI.csv; the reference lists of the first two are
# shortened to their first and last entries.
SOURCE = """\
source,target,weight,TF.category,resources,PMID,sign.decision
MYC,TERT,1,dbTF,ExTRI;HTRI;TRRUST;TFactS;NTNU.Curated;Pavlidis2021;DoRothEA_A,10022128;9988278,PMID
AP1,JUN,1,dbTF,ExTRI;TRRUST;NTNU.Curated,10022869;9834073,PMID
SPI1,BGLAP,1,dbTF,ExTRI,10022617,default activation
SRF,MIR206,1,dbTF,GEREDB,20655308,regulon
AHR,CA9,-1,dbTF,TRRUST;GEREDB;Pavlidis2021,19154183,PMID
ETV3,FGF2,-1,dbTF,ExTRI,10913304,regulon
"""


def rows():
    return list(iter_csv(SimpleNamespace(result=io.StringIO(SOURCE))))


def extract(row):
    record = mod.resource.interactions.mapper(row)
    assert record is not None
    out = SilverExtractor('collectri', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_shared_reader_yields_the_source_columns():
    parsed = rows()
    assert len(parsed) == 6
    assert parsed[2] == {
        'source': 'SPI1',
        'target': 'BGLAP',
        'weight': '1',
        'TF.category': 'dbTF',
        'resources': 'ExTRI',
        'PMID': '10022617',
        'sign.decision': 'default activation',
    }


def test_literature_supported_activation_is_increased_expression():
    record, out = extract(rows()[0])
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'object_aspect_qualifier') == ['expression']
    assert annotations(relation, 'publications') == [
        'PMID:10022128',
        'PMID:9988278',
    ]
    assert annotations(relation, 'supporting_data_source') == [
        'ExTRI',
        'HTRI',
        'TRRUST',
        'TFactS',
        'NTNU.Curated',
        'Pavlidis2021',
        'DoRothEA_A',
    ]
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'protein',
        Namespace.GENESYMBOL,
        'MYC',
    )
    assert (target.entity_type, target.namespace, target.identifier) == (
        'gene',
        Namespace.GENESYMBOL,
        'TERT',
    )
    assert subject.taxon == target.taxon == '9606'


def test_literature_supported_repression_is_decreased_expression():
    _, out = extract(rows()[4])
    relation = out.relations[0]
    assert annotations(relation, 'object_direction_qualifier') == ['decreased']
    assert annotations(relation, 'object_aspect_qualifier') == ['expression']


@pytest.mark.parametrize('index', [2, 3, 5])
def test_default_and_regulon_inferred_signs_are_not_asserted(index):
    record, out = extract(rows()[index])
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    assert all(
        'qualifier' not in slot_name(a.term) for a in record.annotations
    )


def test_dimer_name_is_a_complex_without_invented_members():
    record, out = extract(rows()[1])
    complex_ = out.entities[out.relations[0].subject_entity_key]
    assert complex_.entity_type == 'macromolecular_complex'
    assert (complex_.namespace, complex_.identifier) == (Namespace.NAME, 'AP1')
    assert not record.subject.membership
    assert all(
        i.type != Namespace.GENESYMBOL for i in record.subject.identifiers
    )
    assert {e.entity_type for e in out.entities.values()} == {
        'macromolecular_complex',
        'gene',
    }


def test_mirna_gene_target_keeps_its_gene_symbol():
    _, out = extract(rows()[3])
    target = out.entities[out.relations[0].object_entity_key]
    assert (target.entity_type, target.namespace, target.identifier) == (
        'gene',
        Namespace.GENESYMBOL,
        'MIR206',
    )


def test_only_pubmed_ids_become_publications():
    row = dict(rows()[4], PMID='10405202;S5203;pazar_ABS;PMID:16139800')
    _, out = extract(row)
    assert annotations(out.relations[0], 'publications') == [
        'PMID:10405202',
        'PMID:16139800',
    ]


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    # The file has ``0`` and 16-digit fused numbers (``1451920411085504``).
    row = dict(rows()[4], PMID='0;0123456;123456789;1451920411085504;PMID:19154183')
    _, out = extract(row)
    assert annotations(out.relations[0], 'publications') == ['PMID:19154183']


@pytest.mark.parametrize('missing', ['source', 'target'])
def test_row_without_a_participant_emits_no_relation(missing):
    assert mod.resource.interactions.mapper(dict(rows()[4], **{missing: ''})) is None


@pytest.mark.parametrize(
    'field,value,message',
    [
        ('weight', '0', 'Unknown CollecTRI weight'),
        ('weight', '0.5', 'Unknown CollecTRI weight'),
        ('sign.decision', 'TF role', 'Unknown CollecTRI sign decision'),
        ('sign.decision', '', 'Unknown CollecTRI sign decision'),
    ],
)
def test_unknown_sign_values_are_rejected(field, value, message):
    with pytest.raises(ValueError, match=message):
        mod.resource.interactions.mapper(dict(rows()[4], **{field: value}))
