"""Li 2012 phosphotyrosine circuits: readers, mappers and extractor on source rows."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import entity_type, slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import li2012
from pypath.inputs_v2.parsers import li2012 as li2012_parser

S1_TITLE = ['File S1: A literature-derived human TK circuit dataset ']
S1_HEADER = [
    'substrate family', 'protein/pY position', 'TK writers', 'SH2/PTB Readers',
    'signaling route', 'tissue specificity', 'evo-patterns', None, 'Note:',
]
# Rows copied from File S1: the circuits of three sites, then the blank row
# and the start of the reference list that follows the circuits.
S1_ROWS = [
    ['BCR', 'BCR/177', 'ABL', 'PTPN SH2', 'CTK_Cytoplasm', 'ubiquitous', 'PPP', None, 'X: unknown tissue specificity'],
    ['BCR', 'BCR/177', 'ABL', 'GRB2 SH2', 'CTK_Cytoplasm', 'ubiquitous', 'PPP'],
    ['BCR', 'BCR/177', 'FER', 'PTPN SH2', 'CTK_Cytoplasm', 'ubiquitous', 'PPP'],
    ['BCR', 'BCR/177', 'FER', 'GRB2 SH2', 'CTK_Cytoplasm', 'ubiquitous', 'PPP'],
    ['blnk', 'BLNK/72', 'INSR', 'VAV SH2', 'RTK_Membrane', 'ubiquitous', 'BVP'],
    ['blnk', 'BLNK/72', 'INSR', 'GRB2 SH2', 'RTK_Membrane', 'ubiquitous', 'BVP'],
    ['ITGB3', 'ITGB3/773', 'SRC', 'GRB2 SH2', 'CTK_Membrane', 'tissue specific', 'PPP'],
    ['ITGB3', 'ITGB3/773', 'SRC', 'SHC PTB', 'CTK_Membrane', 'tissue specific', 'PPP'],
    ['blnk', 'BLNK/72', 'SYK', 'VAV SH2', 'CTK_Membrane', 'ubiquitous', 'PVP'],
    ['blnk', 'BLNK/72', 'SYK', 'GRB2 SH2', 'CTK_Membrane', 'ubiquitous', 'PVP'],
    ['ITGB3', 'ITGB3/773', 'SYK', 'GRB2 SH2', 'CTK_Membrane', 'tissue specific', 'PPP'],
    [],
    ['List of References used for constructing the curated dataset [1-71]:'],
    [1.0, 'Manning, B.D. and L.C. Cantley, AKT/PKB signaling: navigating downstream. Cell, 2007. 129(7): p. 1261-74.'],
]
S3_TITLE = [
    'File S3: The tissue-specificity and localization of human TK families '
    'and SH2 domain-containing protein families '
]
S3_HEADER = list(li2012_parser.FAMILY_COLUMNS)
# Rows copied from File S3 (a single gene ID is a number cell in the workbook).
S3_ROWS = [
    ['ABL', 'PTK/SH2', 'ABL2,ABL1,', '27,25', 'ubiquitous', 'cytoplasm'],
    ['SHC', 'SH2/PTB', 'SHC3,SHC4,SHC2,SHC1,', '53358,399694,25759,6464', 'ubiquitous', 'cytoplasm'],
    ['APPL', 'PTB', 'APPL,APPL2', '26060, 55198', 'ubiquitous', 'cytoplasm'],
    ['GULP1', 'PTB', 'gulp1', 51454.0, 'ubiquitous', 'cytoplasm'],
]


def workbook(s1_rows=S1_ROWS, s1_header=S1_HEADER, s3_rows=S3_ROWS):
    book = openpyxl.Workbook()
    for name, rows in (
        ('File S1', [S1_TITLE, s1_header, *s1_rows]),
        ('File S3', [S3_TITLE, S3_HEADER, *s3_rows]),
    ):
        table = book.create_sheet(name)
        for row in rows:
            table.append(row)
    del book['Sheet']
    data = io.BytesIO()
    book.save(data)
    return SimpleNamespace(result=io.BytesIO(data.getvalue()))


def extract(dataset, rows):
    out = SilverExtractor('li2012', dataset)
    mapper = li2012.resource.datasets()[dataset].mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def test_circuit_reader_skips_the_reference_list_and_rejects_bad_rows():
    circuits = list(li2012_parser.iter_circuits(workbook()))
    assert len(circuits) == 11
    assert circuits[0] == {
        'substrate family': 'BCR',
        'protein/pY position': 'BCR/177',
        'TK writers': 'ABL',
        'SH2/PTB Readers': 'PTPN SH2',
        'signaling route': 'CTK_Cytoplasm',
        'tissue specificity': 'ubiquitous',
        'evo-patterns': 'PPP',
    }
    bad = [['BCR', 'BCR/Y177', 'ABL', 'PTPN SH2', 'CTK_Cytoplasm', 'ubiquitous', 'PPP']]
    with pytest.raises(ValueError, match='Malformed Li 2012 circuit'):
        list(li2012_parser.iter_circuits(workbook(s1_rows=bad)))
    header = ['substrate family', 'site', *S1_HEADER[2:]]
    with pytest.raises(ValueError, match='lacks columns'):
        list(li2012_parser.iter_circuits(workbook(s1_header=header)))


def test_circuits_split_into_distinct_kinase_site_and_site_reader_rows():
    kinase_rows = list(li2012_parser.iter_kinase_substrate(workbook()))
    assert [(r['TK writers'], r['substrate'], r['position']) for r in kinase_rows] == [
        ('ABL', 'BCR', 177),
        ('FER', 'BCR', 177),
        ('INSR', 'BLNK', 72),
        ('SRC', 'ITGB3', 773),
        ('SYK', 'BLNK', 72),
        ('SYK', 'ITGB3', 773),
    ]
    # The circuits of a pair keep their other columns as lists.
    assert kinase_rows[0]['SH2/PTB Readers'] == ['PTPN SH2', 'GRB2 SH2']
    assert kinase_rows[0]['evo-patterns'] == ['PPP']

    binding_rows = list(li2012_parser.iter_phosphosite_binding(workbook()))
    assert [
        (r['substrate'], r['position'], r['reader_family'], r['reader_domain'])
        for r in binding_rows
    ] == [
        ('BCR', 177, 'PTPN', 'SH2'),
        ('BCR', 177, 'GRB2', 'SH2'),
        ('BLNK', 72, 'VAV', 'SH2'),
        ('BLNK', 72, 'GRB2', 'SH2'),
        ('ITGB3', 773, 'GRB2', 'SH2'),
        ('ITGB3', 773, 'SHC', 'PTB'),
    ]
    assert binding_rows[3]['TK writers'] == ['INSR', 'SYK']
    assert binding_rows[3]['evo-patterns'] == ['BVP', 'PVP']


def test_kinase_family_increases_phosphorylation_of_the_tyrosine():
    rows = list(li2012_parser.iter_kinase_substrate(workbook()))
    record = li2012.resource.kinase_substrate.mapper(rows[0])
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.type == 'protein_family'
    assert record.subject.identifiers == [(Namespace.NAME, 'ABL')]
    assert record.subject.molecular_form is None
    assert record.object.type == 'protein'
    assert record.object.identifiers == [(Namespace.GENESYMBOL, 'BCR')]
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'Y',
            'position': 177,
            'end_position': 177,
            'coordinate_reference': {
                'identifier': None,
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': 'BCR/177',
        }
    ]

    out = extract('kinase_substrate', rows)
    assert len(out.relations) == 6
    assert len({r.relation_key for r in out.relations}) == 6
    for relation in out.relations:
        assert relation.predicate == 'affects'
        qualifiers = {
            a['term']: a['value'] for a in relation.annotations if a['scope'] == 'relation'
        }
        # Phosphorylation, not activity; the sheet cites no publication per row.
        assert qualifiers == {
            'object_direction_qualifier': 'increased',
            'object_aspect_qualifier': 'phosphorylation',
        }
    proteins = [e for e in out.entities.values() if e.entity_type == 'protein']
    families = [e for e in out.entities.values() if e.entity_type == 'protein_family']
    assert {e.taxon for e in proteins} == {'9606'}
    # A family is not a gene: no gene symbol and no taxon is given to it.
    assert {e.namespace for e in families} == {Namespace.NAME}
    assert {e.taxon for e in families} == {None}


def test_reader_binding_is_an_interaction_without_a_direction_or_directness():
    rows = list(li2012_parser.iter_phosphosite_binding(workbook()))
    record = li2012.resource.phosphosite_binding.mapper(rows[5])
    # The sheet gives no interaction type: not a stated direct interaction.
    assert slot_name(record.predicate) == 'interacts_with'
    assert record.subject.identifiers == [(Namespace.NAME, 'SHC')]
    # The binding domain is not part of the family: a region in the molecular
    # form made a different entity than the family of File S3.
    assert record.subject.molecular_form is None
    site = record.object.molecular_form['modifications'][0]
    assert (site['term'], site['residue'], site['position']) == (
        'phosphorylated residue',
        'Y',
        773,
    )

    out = extract('phosphosite_binding', rows)
    assert [r.predicate for r in out.relations] == ['interacts_with'] * 6
    assert len({r.relation_key for r in out.relations}) == 6
    # The domain the source names after the family is an attribute of the relation.
    assert [
        [a['value'] for a in r.annotations if a['term'] == 'li2012:binding_domain']
        for r in out.relations
    ] == [
        ['SH2'],
        ['SH2'],
        ['SH2'],
        ['SH2'],
        ['SH2'],
        ['PTB'],
    ]
    for relation in out.relations:
        assert {a['term'] for a in relation.annotations if a['scope'] == 'relation'} == {
            'li2012:binding_domain'
        }


def test_reader_family_is_the_family_of_file_s3():
    # SHC is a reader family in the circuits and a family with members in
    # File S3 (the rows are copied from the source).
    binding = extract('phosphosite_binding', list(li2012_parser.iter_phosphosite_binding(workbook())))
    families = extract('families', list(li2012_parser.iter_families(workbook())))
    shc = next(k for k, e in families.entities.items() if e.entity_type == 'protein_family' and e.identifier == 'SHC')
    readers = {r.subject_entity_key for r in binding.relations}
    assert shc in readers
    assert shc in binding.entities
    assert binding.entities[shc].molecular_form is None
    # PTPN, VAV, GRB2 and SHC are the four reader families of the sample.
    assert len(readers) == 4


def test_no_relation_without_a_participant():
    kinase_row = next(li2012_parser.iter_kinase_substrate(workbook()))
    mapper = li2012.resource.kinase_substrate.mapper
    assert mapper({**kinase_row, 'TK writers': None}) is None
    assert mapper({**kinase_row, 'substrate': None, 'position': None}) is None
    binding_row = next(li2012_parser.iter_phosphosite_binding(workbook()))
    mapper = li2012.resource.phosphosite_binding.mapper
    assert mapper({**binding_row, 'reader_family': None, 'reader_domain': None}) is None


def test_families_keep_their_members_and_gene_identifiers():
    rows = list(li2012_parser.iter_families(workbook()))
    assert rows[0]['members'] == [
        {'symbol': 'ABL2', 'entrez': '27'},
        {'symbol': 'ABL1', 'entrez': '25'},
    ]
    assert rows[2]['members'] == [
        {'symbol': 'APPL', 'entrez': '26060'},
        {'symbol': 'APPL2', 'entrez': '55198'},
    ]
    # The source's spelling of a symbol is kept; the gene ID is the accession.
    assert rows[3]['members'] == [{'symbol': 'gulp1', 'entrez': '51454'}]

    record = li2012.map_family(rows[0])
    assert entity_type(record.type) == 'protein_family'
    assert record.identifiers == [(Namespace.NAME, 'ABL')]
    assert [m.member.identifiers for m in record.membership] == [
        [(Namespace.ENTREZ, '27'), (Namespace.GENESYMBOL, 'ABL2')],
        [(Namespace.ENTREZ, '25'), (Namespace.GENESYMBOL, 'ABL1')],
    ]

    out = extract('families', rows)
    assert [r.predicate for r in out.relations] == ['has_member'] * 9
    families = {
        e.identifier: key
        for key, e in out.entities.items()
        if e.entity_type == 'protein_family'
    }
    assert set(families) == {'ABL', 'SHC', 'APPL', 'GULP1'}
    # The family of a kinase-substrate record is the same silver entity.
    kinase = extract('kinase_substrate', list(li2012_parser.iter_kinase_substrate(workbook()))[:1])
    assert families['ABL'] in kinase.entities

    misaligned = [['ABL', 'PTK/SH2', 'ABL2,ABL1,', '27', 'ubiquitous', 'cytoplasm']]
    with pytest.raises(ValueError, match='do not align'):
        list(li2012_parser.iter_families(workbook(s3_rows=misaligned)))
