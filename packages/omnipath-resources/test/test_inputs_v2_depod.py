"""DEPOD phosphatase-substrate sheets: reader, mapper and extractor on source rows."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import depod
from pypath.inputs_v2.parsers import depod as depod_parser

PROTEIN_HEADER = list(depod_parser.PROTEIN_COLUMNS)
# Rows copied from PPase_protSubtrates_201903.xls (a single PubMed ID is a
# number cell in the workbook).
PROTEIN_ROWS = [
    ['P36873', 'PPP1CC', 'O00168', 'FXYD1', 'Ser-88', 'SIRRLsTRRR-', 'In vitro', 26668322.0, 1.0],
    ['P67775', 'PPP2CA', 'O14745', 'SLC9A3R1', 'Ser-339, Ser-340', 'AHQKRsSKRAP HQKRSsKRAPQ', 'In vitro and in vivo', '22301917, 22523600', 3.0],
    ['O75365', 'PTP4A3', 'O14745', 'SLC9A3R1', 'N/A', 'N/A', 'In vivo', 25897829.0, 1.0],
    ['Q6ZVD8', 'PHLPP2', 'P05771', 'PRKCB', 'Ser-660 in Isoform Beta-II ', 'EFEGFsFVNSE', 'In vitro and in vivo', 18162466.0, 2.0],
    ['Q8TF42', 'UBASH3B', 'Q05397', 'PTK2', 'N/A', 'N/A', 'N/A', 'UniProt', 0.0],
    ['Q8IXW5', 'RPAP2', 'P24928', 'POLR2A', 'Ser-1619,Ser-1647,Ser-1661,Ser-1675,Ser-1689,Ser-1703,Ser-1717,Ser-1731,Ser-1766,Ser-1794', 'SYSPTsPSYSP', 'In vitro and in vivo', 22137580.0, 2.0],
    ['Q05209', 'PTPN12', 'P00533', 'EGFR', 'Tyr-1148', '1148 is not a Tyr in HsEGFR', 'in vivo', 21376233.0, 1.0],
    ['Q8WYL5', 'SSH1', 'Q9BR76', 'CORO1B', 'Ser-2', 29, 'In vitro and in vivo', 17350576.0, 2.0],
]
NON_PROTEIN_HEADER = list(depod_parser.NON_PROTEIN_COLUMNS)
# Rows copied from PPase_nonprotSubtrates_201903.xls.
NON_PROTEIN_ROWS = [
    ['Q96T60', 'PNKP', 'CHEBI:1359', 'C03463', "3'-phosphopolynucleotide", "3' position", 'In vitro', 10446192.0, 1.0],
    ['Q9NX62', 'IMPAD1', 'CHEBI:25448', 'NA', 'myo-inositol phosphate', 'N/A', 'N/A', 'UniProt', 1.0],
]


def workbook(sheet, header, rows):
    book = openpyxl.Workbook()
    table = book.active
    table.title = sheet
    table.append(header)
    for row in rows:
        table.append(row)
    data = io.BytesIO()
    book.save(data)
    return SimpleNamespace(result=io.BytesIO(data.getvalue()))


def protein_rows(rows=PROTEIN_ROWS, header=PROTEIN_HEADER):
    return list(depod_parser.iter_protein_substrates(workbook('PSprots', header, rows)))


def non_protein_rows():
    return list(
        depod_parser.iter_non_protein_substrates(
            workbook('PSnonprots', NON_PROTEIN_HEADER, NON_PROTEIN_ROWS)
        )
    )


def extract(dataset, rows):
    out = SilverExtractor('depod', dataset)
    mapper = depod.resource.datasets()[dataset].mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def site(out, relation):
    form = out.entities[relation.object_entity_key].molecular_form
    return form['modifications'][0] if form else None


def test_reader_yields_one_row_per_site_with_its_window():
    rows = protein_rows()
    assert [(r['site'], r['residue'], r['position'], r['window']) for r in rows[:6]] == [
        ('Ser-88', 'S', 88, 'SIRRLsTRRR-'),
        ('Ser-339', 'S', 339, 'AHQKRsSKRAP'),
        ('Ser-340', 'S', 340, 'HQKRSsKRAPQ'),
        # No site in the source.
        (None, None, None, None),
        # Isoform numbering: the text is kept, no position is read.
        ('Ser-660 in Isoform Beta-II', None, None, None),
        (None, None, None, None),
    ]
    # One repeat motif given once for ten sites is not assigned to a site.
    assert [(r['site'], r['position'], r['window']) for r in rows[6:16]] == [
        (f'Ser-{position}', position, None)
        for position in (1619, 1647, 1661, 1675, 1689, 1703, 1717, 1731, 1766, 1794)
    ]
    # A comment in the windows column becomes part of the site text, and no
    # position is read.
    assert (rows[16]['site'], rows[16]['residue'], rows[16]['position'], rows[16]['window']) == (
        'Tyr-1148 (1148 is not a Tyr in HsEGFR)',
        None,
        None,
        None,
    )
    # A number in the windows column is not a window and not a comment.
    assert (rows[17]['site'], rows[17]['position'], rows[17]['window']) == ('Ser-2', 2, None)
    assert len(rows) == 18
    assert rows[0]['pmids'] == ['26668322']
    assert rows[0]['Literature source/s'] == '26668322'
    assert rows[1]['pmids'] == ['22301917', '22523600']
    assert rows[5]['pmids'] == []


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    row = [*PROTEIN_ROWS[1]]
    row[7] = '0123456, 123456789, 22301917'
    parsed = protein_rows([row])
    assert [p['pmids'] for p in parsed] == [['22301917']] * len(parsed)
    for relation in extract('enzyme_substrate', parsed).relations:
        assert annotations(relation, 'publications') == ['PMID:22301917']


def test_reader_rejects_a_changed_header():
    header = ['Phosphatase'] + PROTEIN_HEADER[1:]
    with pytest.raises(ValueError, match='lacks columns'):
        protein_rows(header=header)


def test_phosphatase_decreases_substrate_phosphorylation_at_the_site():
    row = protein_rows()[0]
    record = depod.resource.enzyme_substrate.mapper(row)
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.identifiers == [
        (Namespace.UNIPROT, 'P36873'),
        (Namespace.NAME, 'PPP1CC'),
    ]
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'O00168'),
        (Namespace.NAME, 'FXYD1'),
    ]
    assert record.subject.molecular_form is None
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'de-phosphorylated residue',
            'residue': 'S',
            'position': 88,
            'end_position': 88,
            'coordinate_reference': {
                'identifier': {'ns': 'uniprot', 'id': 'O00168'},
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            # The site with its window as written (the window used to be
            # left out of the description).
            'description': 'Ser-88 (SIRRLsTRRR-)',
        }
    ]

    out = extract('enzyme_substrate', [row])
    (relation,) = out.relations
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'object_direction_qualifier') == ['decreased']
    assert annotations(relation, 'publications') == ['PMID:26668322']
    (score,) = [
        a['quantity'] for a in relation.annotations if a['term'] == 'has_confidence_score'
    ]
    assert score['has_numeric_value'] == 1.0
    assert score['source_field'] == 'Reliability score'
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}


def test_sites_of_one_pair_are_distinct_records():
    out = extract('enzyme_substrate', protein_rows()[1:3])
    first, second = out.relations
    assert first.subject_entity_key == second.subject_entity_key
    assert first.relation_key != second.relation_key
    assert (site(out, first)['position'], site(out, second)['position']) == (339, 340)
    assert annotations(first, 'publications') == ['PMID:22301917', 'PMID:22523600']


def test_missing_and_free_text_sites_get_no_position():
    rows = protein_rows()
    out = extract('enzyme_substrate', rows[3:6])
    no_site, free_text, uniprot_only = out.relations
    assert site(out, no_site) is None
    assert site(out, free_text) == {
        'term': 'de-phosphorylated residue',
        'residue': None,
        'position': None,
        'end_position': None,
        'coordinate_reference': None,
        'description': 'Ser-660 in Isoform Beta-II',
    }
    # "UniProt" is a data source, not a publication.
    assert annotations(uniprot_only, 'publications') == []
    assert annotations(uniprot_only, 'supporting_data_source') == ['UniProt']
    with pytest.raises(ValueError, match='Unknown DEPOD literature source'):
        depod.resource.enzyme_substrate.mapper(
            {**rows[0], 'Literature source/s': 'personal communication'}
        )


def test_dephosphorylation_is_not_an_effect_on_activity():
    out = extract('enzyme_substrate', protein_rows())
    assert len(out.relations) == 18
    for relation in out.relations:
        assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
        assert annotations(relation, 'object_direction_qualifier') == ['decreased']
        assert annotations(relation, 'causal_mechanism_qualifier') == []
    # A row without a phosphatase or substrate accession yields no relation.
    row = protein_rows()[0]
    mapper = depod.resource.enzyme_substrate.mapper
    assert mapper({**row, 'Substrate accession numbers': None, 'Substrate entry names': None}) is None
    assert mapper({**row, 'Phosphatase accession numbers': None, 'Phosphatase entry names': None}) is None


def test_non_protein_substrate_is_a_chemical_without_taxon_or_site():
    rows = non_protein_rows()
    record = depod.resource.non_protein_substrates.mapper(rows[0])
    assert record.object.identifiers == [
        (Namespace.CHEBI, '1359'),
        (Namespace.KEGG, 'C03463'),
        (Namespace.NAME, "3'-phosphopolynucleotide"),
    ]
    assert record.object.molecular_form is None

    out = extract('non_protein_substrates', rows)
    assert [r.predicate for r in out.relations] == ['affects', 'affects']
    chemicals = [e for e in out.entities.values() if e.entity_type == 'chemical_entity']
    assert len(chemicals) == 2
    assert all(e.taxon is None and e.molecular_form is None for e in chemicals)
    # "NA" is not a KEGG compound.
    second = depod.resource.non_protein_substrates.mapper(rows[1])
    assert all(i.type != Namespace.KEGG for i in second.object.identifiers)
    for relation in out.relations:
        assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
        assert annotations(relation, 'object_direction_qualifier') == ['decreased']
    assert annotations(out.relations[1], 'supporting_data_source') == ['UniProt']


def test_assay_type_is_evidence_of_type():
    # The assay type column, as the source writes it: with and without a
    # capital letter, one or both assays, ex vivo, and "N/A".
    rows = protein_rows()
    out = extract('enzyme_substrate', [rows[0], rows[1], rows[3], rows[5]])
    in_vitro, both, in_vivo, missing = [
        annotations(r, 'has_evidence_of_type') for r in out.relations
    ]
    assert in_vitro == ['ECO:0000181']
    assert sorted(both) == ['ECO:0000178', 'ECO:0000181']
    assert in_vivo == ['ECO:0000178']
    assert missing == []
    lower_case = {**rows[0], 'Assay/s to infer dephosphorylation': 'in vitro and in vivo'}
    out = extract('enzyme_substrate', [lower_case])
    assert sorted(annotations(out.relations[0], 'has_evidence_of_type')) == [
        'ECO:0000178',
        'ECO:0000181',
    ]
    non_protein = non_protein_rows()
    ex_vivo = {**non_protein[0], 'Assay/s to infer dephosphorylation': 'In vitro and ex vivo'}
    out = extract('non_protein_substrates', [non_protein[0], ex_vivo, non_protein[1]])
    assert [sorted(annotations(r, 'has_evidence_of_type')) for r in out.relations] == [
        ['ECO:0000181'],
        ['ECO:0000181', 'ECO:0005502'],
        [],
    ]
    with pytest.raises(ValueError, match='Unknown DEPOD assay type'):
        depod.resource.enzyme_substrate.mapper(
            {**rows[0], 'Assay/s to infer dephosphorylation': 'in silico'}
        )


def test_sequence_window_is_in_the_site_description_only_where_aligned():
    rows = protein_rows()
    # One window per site.
    out = extract('enzyme_substrate', rows[1:3])
    assert [site(out, r)['description'] for r in out.relations] == [
        'Ser-339 (AHQKRsSKRAP)',
        'Ser-340 (HQKRSsKRAPQ)',
    ]
    # One repeat motif for ten sites and a number are not windows of the
    # site: the description is the site alone.
    out = extract('enzyme_substrate', [rows[6], rows[17]])
    assert [site(out, r)['description'] for r in out.relations] == ['Ser-1619', 'Ser-2']
    # A site cell with free text keeps its text.
    out = extract('enzyme_substrate', [rows[16]])
    assert site(out, out.relations[0])['description'] == (
        'Tyr-1148 (1148 is not a Tyr in HsEGFR)'
    )
