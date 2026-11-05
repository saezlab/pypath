"""SLC (Meixner et al. 2020): transporter-substrate pairs as the legacy network read them."""

import io
from types import SimpleNamespace

from omnipath_build.silver import SilverExtractor
from omnipath_core.interaction_profiles import interaction_label
import openpyxl

from omnipath_resources.inputs_v2 import slc
from omnipath_resources.inputs_v2.parsers import slc as parsers

ANNOTATION_HEADER = (
    'Entrez_ID', 'Ensembl_ID', 'HGNC_symbol', 'Substrate_class', 'Substrates', 'Substrates_PMID',
    'Coupled_ions', 'Coupled_ions_PMID', 'Transport_mechanism', 'Transport_system_PMID',
    'Subcellular_localization', 'Subcellular_localization_PMID',
)
MAPPING_HEADER = ('Substrate_name_annotation', 'chebi_id', 'chebi_term', 'synonym', 'description', 'iri')


def _opener(sheet, header, rows):
    """An opener over an .xlsx whose first sheet is a README, as in the tables."""
    book = openpyxl.Workbook()
    book.active.title = 'README'
    table = book.create_sheet(sheet)
    for row in (header, *rows):
        table.append(row)
    data = io.BytesIO()
    book.save(data)
    return SimpleNamespace(result=data.getvalue())


def _records(annotation_rows, mapping_rows):
    mapping = _opener(parsers.MAPPING_SHEET, MAPPING_HEADER, mapping_rows)
    return list(
        parsers.transport_substrates(
            _opener(parsers.ANNOTATION_SHEET, ANNOTATION_HEADER, annotation_rows),
            mapping_download=SimpleNamespace(open=lambda **_: mapping),
        )
    )


MAPPING = [
    ('L-glutamate', 'CHEBI:29985', 'L-glutamate(1-)', '', '', ''),
    ('sodium', 'CHEBI:29101', 'sodium(1+)', '', '', ''),
    ('Organic cation drugs', None, None, '', '', ''),
]


def test_substrates_and_coupled_ions_with_chebi_become_pairs():
    records = _records(
        [(
            '6505', 'ENSG00000106688', 'SLC1A1', 'Amino acid',
            'L-glutamate; Organic cation drugs; Unknown', '15369674, 11827462; 18337592',
            'sodium', '12133831 (rabbit)', 'Symporter', None, 'Plasma membrane', None,
        )],
        MAPPING,
    )
    assert [(r['genesymbol'], r['role'], r['chebi'], r['pmids']) for r in records] == [
        ('SLC1A1', 'substrate', 'CHEBI:29985', ['11827462', '15369674', '18337592']),
        ('SLC1A1', 'coupled_ion', 'CHEBI:29101', ['12133831']),
    ]
    assert records[0]['entrez'] == '6505' and records[0]['ensg'] == 'ENSG00000106688'
    assert records[0]['transport_mechanism'] == 'Symporter'
    assert records[0]['localization'] == 'Plasma membrane'


def test_orphans_and_unknown_attributes_have_no_values():
    records = _records(
        [
            ('56172', 'ENSG00000154122', 'ANKH', 'Orphan', None, None, 'Unknown', None, 'Unknown', None, 'Unknown', None),
            ('1', 'ENSG1', 'SLCX', 'Other', 'L-glutamate', None, None, None, 'Unknown', None, 'Unknown', None),
        ],
        MAPPING,
    )
    assert [(r['genesymbol'], r['transport_mechanism'], r['localization'], r['pmids']) for r in records] == [
        ('SLCX', None, None, []),
    ]


def test_pmids_keep_identifiers_and_drop_notes():
    assert parsers.pmids('21128598 9260930 11907186 11907186') == ['11907186', '21128598', '9260930']
    assert parsers.pmids('22943001  - no direct proof') == ['22943001']
    assert parsers.pmids(None) == []


def test_pair_is_a_transport_relation():
    row = {
        'genesymbol': 'SLC1A1', 'entrez': '6505', 'ensg': 'ENSG00000106688',
        'role': 'coupled_ion', 'substrate_name': 'sodium', 'chebi': 'CHEBI:29101',
        'pmids': ['12133831'], 'substrate_class': 'Amino acid',
        'transport_mechanism': 'Symporter', 'localization': 'Plasma membrane',
    }
    ex = SilverExtractor(slc.config.id.value, 'transports')
    ex.process_record(slc.transport_schema.build(row), row, 'transports:1', 1)
    (relation,) = ex.relations
    subject = ex.entities[relation.subject_entity_key]
    obj = ex.entities[relation.object_entity_key]
    assert {(str(i['ns']), i['id']) for i in subject.identifiers} == {
        ('genesymbol', 'SLC1A1'), ('entrez', '6505'), ('ensg', 'ENSG00000106688'),
    }
    assert ('chebi', '29101') in {(str(i['ns']), i['id']) for i in obj.identifiers}
    assert interaction_label(relation.predicate, relation.annotations) == 'Transports'
    assert {(a['term'], a['value']) for a in relation.annotations} >= {
        ('publications', 'PMID:12133831'),
        ('slc:role', 'coupled_ion'),
        ('slc:transport_mechanism', 'Symporter'),
        ('slc:substrate_class', 'Amino acid'),
        ('slc:subcellular_localization', 'Plasma membrane'),
    }
