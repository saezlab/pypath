"""Shared readers of standard file formats."""

import csv
import io
from types import SimpleNamespace

import openpyxl
import pytest

from pypath.inputs_v2.parsers.base import (
    content_text,
    iter_csv,
    iter_excel,
    iter_gmt,
    iter_html_links,
    iter_html_tables,
    iter_psimi_xml,
    iter_semicolon,
    iter_tsv,
    read_bytes,
    read_text,
    read_xml_root,
)


def text(content):
    """A download opened as a text handle."""
    return SimpleNamespace(result=io.StringIO(content))


def workbook(rows, title='S1'):
    book = openpyxl.Workbook()
    book.active.title = title
    for row in rows:
        book.active.append(row)
    content = io.BytesIO()
    book.save(content)
    return SimpleNamespace(result=io.BytesIO(content.getvalue()))


def test_excel_rows_start_below_the_header_and_skip_empty_rows():
    rows = [['Table S1'], [], ['Ligand', 'Receptor', 'PMID'], ['TGFB1', ' TGFBR1 ', 123], [None] * 3, ['EGF', 'EGFR', None]]
    assert list(iter_excel(workbook(rows), sheet='S1', header_row=3)) == [
        {'Ligand': 'TGFB1', 'Receptor': 'TGFBR1', 'PMID': 123},
        {'Ligand': 'EGF', 'Receptor': 'EGFR', 'PMID': None},
    ]
    assert list(iter_excel(workbook(rows[2:]))) == list(iter_excel(workbook(rows), header_row=3))


def test_excel_sheet_must_exist():
    with pytest.raises(KeyError):
        list(iter_excel(workbook([['a'], [1]]), sheet='missing'))


def test_read_bytes_joins_chunks_and_accepts_no_download():
    assert read_bytes(SimpleNamespace(result=iter([b'ab', b'c']))) == b'abc'
    assert read_bytes(SimpleNamespace(result=b'abc')) == b'abc'
    assert read_bytes(SimpleNamespace(result='abc')) == b'abc'
    assert read_bytes(None) == b''


def test_excel_repeated_header_names_need_explicit_fieldnames():
    rows = [['Gene', 'ID', 'Gene', 'ID'], ['A', 1, 'B', 2]]
    with pytest.raises(ValueError, match='repeats'):
        list(iter_excel(workbook(rows)))
    names = ('gene_a', 'id_a', 'gene_b', 'id_b')
    assert list(iter_excel(workbook(rows), fieldnames=names)) == [dict(zip(names, ['A', 1, 'B', 2]))]
    assert list(iter_excel(workbook(rows[1:]), header_row=0, fieldnames=names)) == [dict(zip(names, ['A', 1, 'B', 2]))]


def test_delimited_rows_are_keyed_by_the_header_like_dictreader():
    content = 'a,b,c\n1,"x,y",3\n\n4,5\n6,7,8,9\n , ,\n'
    rows = list(iter_csv(text(content)))
    assert rows == list(csv.DictReader(io.StringIO(content)))
    assert rows[0] == {'a': '1', 'b': 'x,y', 'c': '3'}
    # A short row has None for the columns it lacks, a long one its extra fields under None.
    assert rows[1] == {'a': '4', 'b': '5', 'c': None}
    assert rows[2] == {'a': '6', 'b': '7', 'c': '8', None: ['9']}
    assert list(iter_tsv(text('a\tb\n1\t2\n'))) == [{'a': '1', 'b': '2'}]
    assert list(iter_semicolon(text('a;b\n1;2\n'))) == [{'a': '1', 'b': '2'}]
    assert list(iter_csv(text(''))) == [] and list(iter_csv(None)) == []
    with pytest.raises(ValueError, match='HTML document'):
        list(iter_tsv(text('<!DOCTYPE html>\n<html><body>503</body></html>\n')))


def test_delimited_fieldnames_name_a_file_without_header_or_replace_a_skipped_one():
    names = ('gene', 'site')
    assert list(iter_tsv(text('ABL1\tY226\nEGFR\tY1092\n'), fieldnames=names)) == [
        {'gene': 'ABL1', 'site': 'Y226'},
        {'gene': 'EGFR', 'site': 'Y1092'},
    ]
    # A header that spans two lines is skipped and the names are given instead.
    content = 'Partner 1,Partner 2\nGene,Gene\nA,B\n'
    assert list(iter_csv(text(content), fieldnames=('gene_a', 'gene_b'), skip_rows=2)) == [
        {'gene_a': 'A', 'gene_b': 'B'}
    ]
    # Without fieldnames the row after the skipped lines is the header.
    assert list(iter_csv(text('Table 1\na,b\n1,2\n'), skip_rows=1)) == [{'a': '1', 'b': '2'}]


def test_delimited_strict_rejects_another_field_count_and_skips_blank_rows():
    names = ('x', 'y')
    content = 'h1,h2\n\nA,"b,c"\n , \n,,\nD,e\n'
    assert list(iter_csv(text(content), fieldnames=names, skip_rows=1, strict=True)) == [
        {'x': 'A', 'y': 'b,c'},
        {'x': 'D', 'y': 'e'},
    ]
    # The line number counts every line of the file, skipped and blank ones too.
    with pytest.raises(ValueError, match='Line 4 has 3 fields, expected 2'):
        list(iter_csv(text('h1,h2\n\nA,b\nC,d,e\n'), fieldnames=names, skip_rows=1, strict=True))
    with pytest.raises(ValueError, match='Line 2 has 1 fields, expected 2'):
        list(iter_tsv(text('a\tb\n<html>429 Too Many Requests</html>\n'), strict=True))
    # Not strict: the same rows pass, short and long.
    assert len(list(iter_csv(text('A,b\nC,d,e\nF\n'), fieldnames=names))) == 3


def test_delimited_quotechar_none_reads_quotes_as_data():
    content = 'id\tname\n1\t"alpha\n2\t5\'-3\' "exo"\n'
    assert list(iter_tsv(text(content), quotechar=None)) == [
        {'id': '1', 'name': '"alpha'},
        {'id': '2', 'name': '5\'-3\' "exo"'},
    ]
    # With quoting, the unbalanced quote swallows the rest of the file.
    assert len(list(iter_tsv(text(content)))) == 1


def test_delimited_repeated_header_names_are_an_error_unless_allowed():
    content = 'Kinase;Entrez GeneID;Gene;Entrez GeneID;\nAKT1;207;ESR1;2099;\n'
    with pytest.raises(ValueError, match=r"repeats column names \['Entrez GeneID'\]"):
        list(iter_semicolon(text(content)))
    # The opt-out keeps what csv.DictReader does: the last column of the name.
    assert list(iter_semicolon(text(content), allow_repeated=True)) == [
        {'Kinase': 'AKT1', 'Entrez GeneID': '2099', 'Gene': 'ESR1', '': ''}
    ]
    assert list(iter_semicolon(text(content), allow_repeated=True)) == list(
        csv.DictReader(io.StringIO(content), delimiter=';')
    )
    names = ('kinase', 'kinase_id', 'gene', 'gene_id', 'empty')
    assert list(iter_semicolon(text(content), fieldnames=names, skip_rows=1)) == [
        dict(zip(names, ['AKT1', '207', 'ESR1', '2099', '']))
    ]
    # Several unnamed columns are not a repeated name.
    assert list(iter_csv(text('a,,\n1,2,3\n'))) == [{'a': '1', '': '3'}]
    with pytest.raises(ValueError, match='repeats'):
        list(iter_csv(text('1,2\n'), fieldnames=('a', 'a')))


def test_delimited_readers_accept_binary_handles_whole_content_and_line_iterators():
    expected = [{'a': 'é', 'b': '2'}]
    assert list(iter_csv(SimpleNamespace(result=io.BytesIO('a,b\r\né,2\r\n'.encode())))) == expected
    assert list(iter_csv(SimpleNamespace(result='\ufeffa,b\né,2\n'))) == expected
    assert list(iter_csv(SimpleNamespace(result=b'\xef\xbb\xbfa,b\n\xc3\xa9,2\n'))) == expected
    assert list(iter_csv(SimpleNamespace(result=iter(['a,b\n', 'é,2\n'])))) == expected
    assert list(iter_csv(SimpleNamespace(result={'member.csv': io.StringIO('a,b\né,2\n')}))) == expected
    latin = SimpleNamespace(result=io.BytesIO('a,b\né,2\n'.encode('latin-1')))
    assert list(iter_csv(latin, encoding='latin-1')) == expected


def test_gmt_sets_have_a_name_a_description_and_members():
    lines = ['IGF1/2*\tna\tIGF1\tIGF2\t\n', '\n', 'RNF11\tna\tRNF11\n', 'EMPTY\tno members\n']
    expected = [
        {'name': 'IGF1/2*', 'description': 'na', 'members': ['IGF1', 'IGF2']},
        {'name': 'RNF11', 'description': 'na', 'members': ['RNF11']},
        {'name': 'EMPTY', 'description': 'no members', 'members': []},
    ]
    assert list(iter_gmt(SimpleNamespace(result=lines))) == expected
    assert list(iter_gmt(text(''.join(lines)))) == expected
    assert list(iter_gmt(SimpleNamespace(result=''.join(lines).encode()))) == expected
    assert list(iter_gmt(None)) == []
    with pytest.raises(ValueError, match='line 2 has no description column'):
        list(iter_gmt(SimpleNamespace(result=['A\tna\tB\n', 'IGF1\n'])))
    with pytest.raises(ValueError, match='HTML document'):
        list(iter_gmt(text('<html><body>Not found</body></html>\n')))


def test_read_text_reads_any_download_as_one_string():
    assert read_text(text('a\nb\n')) == 'a\nb\n'
    assert read_text(SimpleNamespace(result=io.BytesIO('é'.encode()))) == 'é'
    assert read_text(SimpleNamespace(result={'member': b'text', 'other': b'more'})) == 'text'
    assert read_text(SimpleNamespace(result=iter([b'a\n', 'b\n']))) == 'a\nb\n'
    assert read_text(SimpleNamespace(result='é'.encode('latin-1')), encoding='latin-1') == 'é'
    assert read_text(None) == '' and read_text(SimpleNamespace(result=None)) == ''
    assert content_text(None) == '' and content_text(b'x') == 'x' and content_text('x') == 'x'


def test_html_tables_have_rows_of_cells_with_text_links_and_content():
    html = (
        '<table class="outer"><tr><th colspan="2">Head</th></tr>'
        '<tr><td>alpha <b>beta</b><!-- note -->gamma</td>'
        "<td title='t'><a href='/x?p=1&amp;q=2' id='first'>one<br>two</a> between <a href='/y'>second</a>"
        ' a &lt; b<a name="anchor">plain</a>'
        '<table id="inner"><tr><td>nested</td></tr></table> after</td></tr>'
        '<tr></tr>'
        '<tr><td><script>ignored("<table>")</script>kept<style>td {}</style>'
        '<tr><td>unclosed</table><p>outside</p>'
    )
    outer, inner = iter_html_tables(html)
    assert outer['attrs'] == {'class': 'outer'}
    assert [[cell['text'] for cell in row] for row in outer['rows']] == [
        ['Head'],
        ['alpha betagamma', 'one two between second a < bplain after'],
        ['kept'],
        ['unclosed'],
    ]
    (head,) = outer['rows'][0]
    assert head['header'] and head['attrs'] == {'colspan': '2'} and head['content'] == ['Head']
    first, second = outer['rows'][1]
    assert not first['header'] and first['links'] == []
    # Text nodes as written, split where a tag or a comment splits them.
    assert first['content'] == ['alpha ', 'beta', 'gamma']
    assert second['attrs'] == {'title': 't'}
    assert second['links'] == ['/x?p=1&q=2', '/y']
    assert second['content'] == [
        {'href': '/x?p=1&q=2', 'text': 'onetwo', 'attrs': {'href': '/x?p=1&q=2', 'id': 'first'}},
        ' between ',
        {'href': '/y', 'text': 'second', 'attrs': {'href': '/y'}},
        ' a < b',
        'plain',
        ' after',
    ]
    assert inner['attrs'] == {'id': 'inner'}
    assert [[cell['text'] for cell in row] for row in inner['rows']] == [['nested']]
    assert set(first) == {'header', 'attrs', 'text', 'links', 'content'}
    assert list(iter_html_tables('<p>no table</p>')) == []


def test_html_links_have_an_address_text_and_attributes():
    links = list(iter_html_links('<a name="x">anchor</a><a href="a?b=1" class="c">  two\n <b>words</b> </a><a href="c">'))
    assert links == [
        {'href': 'a?b=1', 'attrs': {'href': 'a?b=1', 'class': 'c'}, 'text': 'two words'},
        {'href': 'c', 'attrs': {'href': 'c'}, 'text': ''},
    ]


def test_read_xml_root_parses_handles_content_archives_and_chunks():
    xml = '<root><item id="1"/></root>'
    for result in (io.BytesIO(xml.encode()), io.StringIO(xml), xml, xml.encode(), {'a.xml': io.BytesIO(xml.encode())}, iter([b'<root><item ', b'id="1"/></root>'])):
        root = read_xml_root(SimpleNamespace(result=result))
        assert root.tag == 'root' and root[0].get('id') == '1'
    assert read_xml_root(None) is None and read_xml_root(SimpleNamespace(result=None)) is None


PSIMI_25 = (
    '<entrySet level="2" version="5" xmlns="net:sf:psidev:mi"><entry>'
    '<source><names><shortLabel>NetPath</shortLabel></names></source>'
    '<availabilityList><availability id="1">Free  to\nuse</availability></availabilityList>'
    '<experimentList><experimentDescription id="10"><bibref><xref><primaryRef db="pubmed" id="123"/></xref></bibref>'
    '<hostOrganismList><hostOrganism ncbiTaxId="9606"><names><shortLabel>Human</shortLabel></names></hostOrganism></hostOrganismList>'
    '<interactionDetectionMethod><names><shortLabel>two  hybrid</shortLabel></names>'
    '<xref><primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0018" refType="identity"/></xref></interactionDetectionMethod>'
    '</experimentDescription></experimentList>'
    '<interactorList><interactor id="1"><names><shortLabel>TSLP</shortLabel><alias type="gene name"> TSLP1 </alias></names>'
    '<xref><primaryRef db="HPRD" id="06108"/><secondaryRef db="Entrez gene" id="85480"/></xref>'
    '<interactorType><names><shortLabel>protein</shortLabel></names><xref><primaryRef db="psi-mi" id="MI:0326"/></xref></interactorType>'
    '<organism ncbiTaxId="9606"/></interactor></interactorList>'
    '<interactionList><interaction id="7"><names><fullName> TSLP binds CRLF2 </fullName></names>'
    '<experimentList><experimentRef>10</experimentRef></experimentList>'
    '<participantList><participant id="3"><interactorRef>1</interactorRef>'
    '<experimentalRoleList><experimentalRole><names><shortLabel>bait</shortLabel></names>'
    '<xref><primaryRef db="psi-mi" id="MI:0496"/></xref></experimentalRole></experimentalRoleList></participant>'
    '<participant id="4"><interactor id="2"><names><shortLabel>CRLF2</shortLabel></names></interactor></participant>'
    '</participantList><interactionType><xref><primaryRef db="psi-mi" id="MI:0218"/></xref></interactionType>'
    '</interaction></interactionList></entry></entrySet>'
)
PSIMI_1 = (
    '<entrySet xmlns="net:sf:psidev:mi" version="1" level="1"><entry>'
    '<source release="1.0" releaseDate="2004-01-01"><names><shortLabel>MPPI</shortLabel></names>'
    '<bibref><xref><primaryRef db="pubmed" id="456"/></xref></bibref></source>'
    '<experimentList><experimentDescription id="E1"><bibref><xref><primaryRef db="pubmed" id="123"/></xref></bibref>'
    '<hostOrganism ncbiTaxId="4932"><names><shortLabel>yeast</shortLabel></names></hostOrganism>'
    '<interactionDetection><names><shortLabel>two  hybrid</shortLabel></names>'
    '<xref><primaryRef db="MI" id="0018"/></xref></interactionDetection>'
    '<confidence unit="score" value="0.9"/></experimentDescription></experimentList>'
    '<interactorList><proteinInteractor id="P1"><names><shortLabel>a</shortLabel></names>'
    '<xref><primaryRef db=" Swiss-Prot " id="P12345" version="2"/></xref><organism ncbiTaxId="10090"/></proteinInteractor></interactorList>'
    '<interactionList><interaction><experimentList><experimentRef ref="E1"/></experimentList>'
    '<participantList><proteinParticipant><proteinInteractorRef ref="P1"/><role>bait</role>'
    '<isTaggedProtein>true</isTaggedProtein></proteinParticipant>'
    '<proteinParticipant><proteinInteractor id="P2"><names><shortLabel>b</shortLabel></names></proteinInteractor>'
    '<role>prey</role></proteinParticipant></participantList></interaction></interactionList></entry></entrySet>'
)


def test_psimi_xml_level_2_5_resolves_references_and_keeps_values_as_written():
    (interaction,) = iter_psimi_xml(PSIMI_25)
    assert list(interaction) == [
        'entry_index', 'source', 'availability', 'id', 'label', 'name', 'xrefs',
        'interaction_types', 'experiments', 'participants',
    ]
    assert interaction['entry_index'] == 1 and interaction['id'] == '7'
    assert interaction['source'] == {'label': 'NetPath', 'name': None, 'xrefs': []}
    assert interaction['availability'] == ['Free to use']
    assert (interaction['label'], interaction['name']) == (None, 'TSLP binds CRLF2')
    assert [term['mi'] for term in interaction['interaction_types']] == ['MI:0218']
    (experiment,) = interaction['experiments']
    assert list(experiment) == [
        'id', 'label', 'name', 'bibliography', 'xrefs', 'host_organisms',
        'detection_method', 'participant_identification_method',
    ]
    assert experiment['bibliography'] == [
        {'primary': True, 'db': 'pubmed', 'db_ac': None, 'id': '123', 'ref_type': None, 'ref_type_ac': None}
    ]
    assert experiment['host_organisms'] == [{'taxon': '9606', 'label': 'Human', 'name': None}]
    # Level 2.5 strips a value and keeps the whitespace inside it.
    assert experiment['detection_method']['label'] == 'two  hybrid'
    assert experiment['detection_method']['mi'] == 'MI:0018'
    assert experiment['detection_method']['xrefs'][0]['db_ac'] == 'MI:0488'
    assert experiment['participant_identification_method'] is None
    first, second = interaction['participants']
    assert list(first) == ['id', 'interactor', 'biological_role', 'experimental_roles']
    assert first['id'] == '3' and first['biological_role'] is None
    assert [role['mi'] for role in first['experimental_roles']] == ['MI:0496']
    assert first['interactor']['aliases'] == [{'value': 'TSLP1', 'type': 'gene name', 'type_ac': None}]
    assert [(ref['db'], ref['id'], ref['primary']) for ref in first['interactor']['xrefs']] == [
        ('HPRD', '06108', True), ('Entrez gene', '85480', False),
    ]
    assert first['interactor']['type']['mi'] == 'MI:0326'
    assert first['interactor']['organism'] == {'taxon': '9606', 'label': None, 'name': None}
    assert second['interactor']['label'] == 'CRLF2' and second['interactor']['xrefs'] == []
    assert list(iter_psimi_xml(PSIMI_25.encode(), level='2.5')) == [interaction]


def test_psimi_xml_level_1_has_its_own_keys():
    (interaction,) = iter_psimi_xml(PSIMI_1)
    assert list(interaction) == [
        'entry_index', 'index', 'source', 'label', 'name', 'xrefs',
        'interaction_types', 'confidence', 'experiments', 'participants',
    ]
    assert interaction['index'] == 1 and interaction['confidence'] is None
    assert interaction['source'] == {
        'label': 'MPPI', 'name': None, 'xrefs': [],
        'bibliography': [{'primary': True, 'db': 'pubmed', 'id': '456', 'secondary': None, 'version': None}],
        'release': '1.0', 'release_date': '2004-01-01',
    }
    (experiment,) = interaction['experiments']
    assert experiment['id'] == 'E1'
    assert experiment['host_organisms'] == [{'taxon': '4932', 'label': 'yeast', 'name': None}]
    # Level 1 collapses whitespace inside a value and writes the accession as MI:0018.
    assert experiment['detection_method']['label'] == 'two hybrid'
    assert experiment['detection_method']['mi'] == 'MI:0018'
    assert experiment['confidence'] == {'unit': 'score', 'value': '0.9'}
    first, second = interaction['participants']
    assert list(first) == ['interactor', 'role', 'confidence', 'is_tagged', 'is_overexpressed']
    assert (first['role'], first['is_tagged'], first['is_overexpressed']) == ('bait', True, None)
    assert first['interactor'] == {
        'id': 'P1', 'label': 'a', 'name': None,
        'xrefs': [{'primary': True, 'db': 'Swiss-Prot', 'id': 'P12345', 'secondary': None, 'version': '2'}],
        'organism': {'taxon': '10090', 'label': None, 'name': None},
    }
    assert second['interactor']['id'] == 'P2' and second['role'] == 'prey'
    assert list(iter_psimi_xml(PSIMI_1, level='1')) == [interaction]


def test_psimi_xml_rejects_other_levels_and_other_documents():
    with pytest.raises(ValueError, match="Unsupported PSI-MI XML level '3'"):
        list(iter_psimi_xml(PSIMI_25.replace('level="2" version="5"', 'level="3" version="0"')))
    with pytest.raises(ValueError, match='Unsupported PSI-MI XML level None'):
        list(iter_psimi_xml('<entrySet xmlns="net:sf:psidev:mi"/>'))
    # The level a caller expects is checked, because the keys differ.
    with pytest.raises(ValueError, match='Expected PSI-MI XML level 1, found level 2.5'):
        list(iter_psimi_xml(PSIMI_25, level='1'))
    with pytest.raises(ValueError, match='Expected PSI-MI XML level 2.5, found level 1'):
        list(iter_psimi_xml(PSIMI_1, level='2.5'))
    with pytest.raises(ValueError, match='entrySet'):
        list(iter_psimi_xml('<html><body>503</body></html>'))
    with pytest.raises(ValueError, match="unknown interactor '9'"):
        list(iter_psimi_xml(PSIMI_25.replace('<interactorRef>1</interactorRef>', '<interactorRef>9</interactorRef>')))
    with pytest.raises(ValueError, match="unknown experiment 'E2'"):
        list(iter_psimi_xml(PSIMI_1.replace('<experimentRef ref="E1"/>', '<experimentRef ref="E2"/>')))


def test_required_columns_and_rows_are_checked_when_asked_for():
    import io

    from pypath.inputs_v2.parsers.base import iter_tsv

    table = lambda text: SimpleNamespace(result=io.StringIO(text))
    assert list(iter_tsv(table('a\tb\n1\t2\n'), required=('a', 'b'))) == [{'a': '1', 'b': '2'}]
    with pytest.raises(ValueError, match='lacks expected columns: c'):
        list(iter_tsv(table('a\tb\n1\t2\n'), required=('a', 'c')))
    with pytest.raises(ValueError, match='no rows'):
        list(iter_tsv(table('a\tb\n'), required=('a',)))
    assert list(iter_tsv(table('a\tb\n'))) == []
    with pytest.raises(ValueError, match='lacks expected columns: PMID'):
        list(iter_excel(workbook([['Ligand'], ['EGF']]), required=('Ligand', 'PMID')))
    with pytest.raises(ValueError, match='no rows'):
        list(iter_excel(workbook([['Ligand']]), required=('Ligand',)))


def test_pubmed_identifier_is_one_to_eight_digits_without_a_leading_zero():
    import re

    from pypath.inputs_v2._publications import PMID_EXTRACT, is_pmid, pmid_curie

    assert pmid_curie('8986718') == 'PMID:8986718'
    assert pmid_curie(' 12345 ') == 'PMID:12345'
    assert pmid_curie(26668322) == 'PMID:26668322'
    for value in ('0', '012345', '169175071', '26668322.0', 'N.N.', '', None, 26668322.0, True):
        assert pmid_curie(value) is None
        assert not is_pmid(value)
    assert re.match(PMID_EXTRACT, '7').group(1) == '7'
    assert re.match(PMID_EXTRACT, '12345678x') is None


def test_rdata_reads_named_objects_and_rows_in_slices(tmp_path, monkeypatch):
    import pandas as pd
    import rdata

    from pypath.inputs_v2.parsers import base
    from pypath.inputs_v2.parsers.base import iter_rdata

    first = pd.DataFrame({'gene': pd.Series(['EGFR', ' TP53 ', ''], dtype=object), 'weight': [1.5, -2.0, 0.25]})
    second = pd.DataFrame({'gene': pd.Series(['MYC'], dtype=object), 'weight': [3.0]})

    def opener(path):
        return SimpleNamespace(result=io.BytesIO(path.read_bytes()))

    rda = tmp_path / 'two.rda'
    rdata.write_rda(rda, {'first': first, 'second': second})
    rows = list(iter_rdata(opener(rda), name='first'))
    # Text is stripped, an empty value is None, numbers are plain values.
    assert rows == [
        {'gene': 'EGFR', 'weight': 1.5},
        {'gene': 'TP53', 'weight': -2.0},
        {'gene': None, 'weight': 0.25},
    ]
    assert [type(value) for value in rows[0].values()] == [str, float]
    assert list(iter_rdata(opener(rda), name='second')) == [{'gene': 'MYC', 'weight': 3.0}]
    # A file with two objects needs the name, and the name must be one of them.
    with pytest.raises(ValueError, match='name the object'):
        list(iter_rdata(opener(rda)))
    with pytest.raises(ValueError, match="not 'third'"):
        list(iter_rdata(opener(rda), name='third'))
    # The rows do not depend on the size of a slice.
    monkeypatch.setattr(base, 'RDATA_ROWS', 2)
    assert list(iter_rdata(opener(rda), name='first')) == rows
    # An .rds file holds one object without a name.
    rds = tmp_path / 'one.rds'
    rdata.write_rds(rds, first)
    assert list(iter_rdata(opener(rds))) == rows
    assert list(iter_rdata(opener(rds), name='first')) == rows
    assert list(iter_rdata(SimpleNamespace(result=b''))) == []

