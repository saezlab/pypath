"""Source-specific Biolink contract for the miRecords input."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from pypath.internals.cv_terms import LicenseCV

from pypath.inputs_v2 import mirecords as mod
from pypath.inputs_v2.parsers.mirecords import COLUMNS, iter_interactions, repair_encoding

N = len(COLUMNS)


def line(pmid, target_species, gene, refseq, mirna_species, mirna, **cells):
    """A row of miRecords_v4.xls; the cells by column position (zero based)."""
    row = [''] * N
    row[0:4] = [pmid, target_species, gene, refseq]
    row[4] = 'NA'
    row[5:7] = [mirna_species, mirna]
    for index, value in cells.items():
        row[int(index.lstrip('c'))] = value
    return row


# Rows copied from the sheet (cell values as xlrd reads them: PubMed IDs and site
# numbers are floats, a gene name that Excel turned into a date is a float).
ROWS = [
    # 0: a human pair; the miRNA level was raised, the mRNA level did not change.
    line(
        12808467.0, 'Homo sapiens', 'C21orf33', 'NM_198155', 'Homo sapiens', 'hsa-miR-23a',
        c4=1.0, c7='overexpression by mature miRNA transfection}{mutation',
        c8='luciferase', c10='{Western blot}{Western blot}', c11='{unchanged}{unchanged}',
    ),
    # 1: a miRNA name between square brackets.
    line(
        18803879.0, 'Homo sapiens', 'BCL2', 'NM_000633.2', 'Homo sapiens', '[hsa-miR-34a]',
    ),
    # 2: no species prefix in the name, and a space around the name.
    line(
        18728182.0, 'Mus musculus', 'Pten', 'NM_008960.2 ', 'Mus musculus', ' [miR-29a] ',
    ),
    # 3: the gene name is an Excel date (1-Sep-2008), the accession is kept.
    line(15685193.0, 'Homo sapiens', 39692.0, 'NM_004404', 'Homo sapiens', 'hsa-miR-373'),
    # 4: mojibake in the gene name.
    line(
        20351193.0, 'Homo sapiens', 'Œ∫B-Ras2 ', 'NM_017595.5 ', 'Homo sapiens', 'hsa-let-7a',
    ),
    # 5: a virus as the target; no gene name, the genome accession remains.
    line(
        17613256.0, 'Vesicular stomatitis Indiana virus', '', 'NC_001560',
        'Mus musculus', 'mmu-miR-93',
    ),
    # 6: species written 'human' and 'Danio rerio ' with a space.
    line(20484036.0, 'human', 'mTOR', 'NM_004958.3', 'Homo sapiens', 'hsa-miR-99a'),
    line(21148433.0, 'Danio rerio ', 'spred1', 'NM_213232.1', 'Danio rerio ', 'dre-miR-126a'),
    # 8: species cells that do not name a species.
    line(15845854.0, '[Virus]', 'HIVgene', 'NM_000001', 'Homo sapiens', 'hsa-miR-32'),
    line(21742974.0, 'Mus musculus', 'SOCS1', 'NM_009896.2', 'us musculus', 'mmu-let-7i'),
    # 10: nothing identifies the target.
    line(15854907.0, 'Ovis aries', 'unknown', 'unknown', 'Ovis aries', 'oar-miR-433'),
    # 11: a GenBank accession in the RefSeq column, and PubMed ID 0.
    line(0.0, 'Homo sapiens', 'ABC1', 'AF000673', 'Homo sapiens', 'hsa-miR-22'),
    # 12: an accession that is not valid, with a gene name.
    line(18952845.0, 'Homo sapiens', 'LAMB3', 'NM213137', 'Homo sapiens', 'hsa-miR-1'),
]


def workbook(rows=ROWS, header=COLUMNS):
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.append(list(header))
    for row in rows:
        sheet.append(row)
    content = io.BytesIO()
    book.save(content)
    content.seek(0)
    return content


def parsed(rows=ROWS, header=COLUMNS):
    return list(iter_interactions(SimpleNamespace(result=workbook(rows, header))))


def extract(row):
    record = mod.map_interaction(row)
    assert record is not None
    out = SilverExtractor('mirecords', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def test_module_has_one_dataset_and_a_noncommercial_license():
    assert list(mod.resource.datasets()) == ['interactions']
    assert mod.config.license == LicenseCV.NONCOMMERCIAL


def test_parser_reads_the_rows_by_position_and_decodes_the_cells():
    rows = parsed()
    assert len(rows) == len(ROWS)
    assert rows[0] == {
        'pmid': '12808467',
        'target_species': 'Homo sapiens',
        'target_gene': 'C21orf33',
        'target_gene_excel_date': None,
        'target_refseq': 'NM_198155',
        'site_number': 1,
        'mirna_species': 'Homo sapiens',
        'mirna': 'hsa-miR-23a',
        'mirna_bracketed': False,
        'mirna_regulation': [
            'overexpression by mature miRNA transfection',
            'mutation',
        ],
        'reporter': 'luciferase',
        'test_method': ['Western blot', 'Western blot'],
        'mrna_level': ['unchanged', 'unchanged'],
    }
    # 'NA' (not applicable) is not a value.
    assert rows[1]['site_number'] is None
    assert rows[1]['mirna_regulation'] == []


def test_parser_repairs_only_what_is_exact():
    rows = parsed()
    # Spaces around values and inside species names.
    assert rows[2]['target_refseq'] == 'NM_008960.2'
    assert rows[7]['target_species'] == rows[7]['mirna_species'] == 'Danio rerio'
    # Brackets are taken off the name and noted.
    assert (rows[1]['mirna'], rows[1]['mirna_bracketed']) == ('hsa-miR-34a', True)
    assert (rows[2]['mirna'], rows[2]['mirna_bracketed']) == ('miR-29a', True)
    # Mac Roman text that is UTF-8 is converted back; other text is not touched.
    assert rows[4]['target_gene'] == 'κB-Ras2'
    assert repair_encoding('Smad 4') == 'Smad 4'
    assert repair_encoding('Café') == 'Café'
    # An Excel date is not a gene name.
    assert rows[3]['target_gene'] is None
    assert rows[3]['target_gene_excel_date'] == 39692
    # 'unknown' is a missing value.
    assert rows[10]['target_gene'] is None
    assert rows[10]['target_refseq'] is None
    # A species that is not a name stays as it is.
    assert rows[8]['target_species'] == '[Virus]'
    assert rows[9]['mirna_species'] == 'us musculus'


def test_parser_rejects_changed_columns_and_empty_tables():
    changed = list(COLUMNS)
    changed[3] = 'Target gene_accession'
    with pytest.raises(ValueError, match='columns'):
        parsed(header=changed)
    with pytest.raises(ValueError, match='no rows'):
        parsed(rows=[])
    bad = line('abc', 'Homo sapiens', 'MYC', 'NM_002467', 'Homo sapiens', 'hsa-miR-1')
    with pytest.raises(ValueError, match='Pubmed_id'):
        parsed(rows=[bad])


def test_pair_is_a_mirna_affecting_a_gene_without_a_direction():
    relation, subject, target = extract(parsed()[0])
    assert relation.predicate == 'affects'
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'microrna',
        Namespace.NAME,
        'hsa-miR-23a',
    )
    assert (target.entity_type, target.namespace, target.identifier) == (
        'gene',
        Namespace.GENESYMBOL,
        'C21orf33',
    )
    assert (Namespace.REFSEQ, 'NM_198155') in identifiers(target)
    assert subject.taxon == target.taxon == '9606'
    assert annotations(relation, 'publications') == ['PMID:12808467']
    # The columns name no effect on the target: nothing may be invented.
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    assert not annotations(relation, 'has_confidence_level')


def test_bracketed_name_is_kept_without_brackets_and_not_translated():
    rows = parsed()
    _, subject, _ = extract(rows[1])
    assert subject.identifier == 'hsa-miR-34a'
    _, subject, _ = extract(rows[2])
    assert subject.identifier == 'miR-29a'
    assert subject.taxon == '10090'


def test_excel_date_gene_is_identified_by_its_accession_only():
    _, _, target = extract(parsed()[3])
    assert (target.namespace, target.identifier) == (Namespace.REFSEQ, 'NM_004404')
    assert all(ns != Namespace.GENESYMBOL for ns, _ in identifiers(target))
    _, _, target = extract(parsed()[5])
    assert (target.namespace, target.identifier) == (Namespace.REFSEQ, 'NC_001560')
    assert target.taxon == '11277'


def test_species_come_from_exact_names():
    rows = parsed()
    assert extract(rows[6])[2].taxon == '9606'  # 'human'
    assert extract(rows[7])[1].taxon == extract(rows[7])[2].taxon == '7955'
    # A cell that is not a species name gives no taxon, on either side.
    assert extract(rows[8])[2].taxon is None
    assert extract(rows[8])[1].taxon == '9606'
    assert extract(rows[9])[1].taxon is None
    assert extract(rows[9])[2].taxon == '10090'


def test_unknown_species_name_raises():
    row = dict(parsed()[0], target_species='Homo erectus')
    with pytest.raises(ValueError, match='species'):
        mod.map_interaction(row)


def test_row_without_gene_name_and_accession_gives_no_relation():
    rows = parsed()
    assert mod.map_interaction(rows[10]) is None
    # The Excel date alone is not enough either.
    assert mod.map_interaction(dict(rows[3], target_refseq=None)) is None
    # An accession of a kind the module does not know is not an identifier...
    assert mod.map_interaction(dict(rows[12], target_gene=None)) is None
    # ...but the gene name is one.
    _, _, target = extract(rows[12])
    assert identifiers(target) == {(Namespace.GENESYMBOL, 'LAMB3')}


@pytest.mark.parametrize('number', ['0123456', '123456789'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(number):
    relation, _, _ = extract(dict(parsed()[0], pmid=number))
    assert not annotations(relation, 'publications')


def test_genbank_accession_and_missing_publication():
    relation, _, target = extract(parsed()[11])
    assert (Namespace.GENBANK, 'AF000673') in identifiers(target)
    assert not any(ns == Namespace.REFSEQ for ns, _ in identifiers(target))
    # PubMed ID 0 is no publication.
    assert not annotations(relation, 'publications')
