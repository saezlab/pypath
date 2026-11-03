"""Source contract of the PDZBase interaction input."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.parsers.base import iter_html_tables
from pypath.inputs_v2.parsers.pdzbase import iter_pdzbase
from pypath.inputs_v2.pdzbase import resource


FIELDS = (
    'pdz_name', 'pdz_accession', 'pdz_domain', 'ligand_name',
    'ligand_accession', 'organism', 'reference', 'interaction_id',
)


def row(*values):
    return dict(zip(FIELDS, values, strict=True))


# Records of the PDZBase "all interactions" page.
HUMAN_ROW = row(
    'NHERF-1', 'O14745', '1', 'CFTR', 'P13569', 'human', '9613608', '31001',
)
ISOFORM_ROW = row(
    'AF-6', 'P55196', '1', 'EPHB2', 'P29323-2', 'human', '9707552', '30584',
)
# Labelled ``human``: the PDZ protein is of rat and the ligand of chicken.
MIXED_SPECIES_ROW = row(
    'MUPP1', 'O55164', '10', 'TAPP-2', 'Q9DF37', 'human', '11802782', '30968',
)
# The same pair and publication, bound by another PDZ domain of MUPP1.
OTHER_DOMAIN_ROW = dict(MIXED_SPECIES_ROW, pdz_domain='13', interaction_id='30971')


def extract(source_row):
    record = resource.datasets()['interactions'].mapper(source_row)
    assert record is not None
    out = SilverExtractor('pdzbase', 'interactions')
    out.process_record(record, source_row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def annotations(relation):
    return sorted((a['scope'], a['term'], a['value']) for a in relation.annotations)


def test_interaction_is_pdz_protein_to_ligand_with_its_publication():
    relation, pdz, ligand = extract(HUMAN_ROW)
    # The page states no interaction type per record.
    assert relation.predicate == 'interacts_with'
    assert relation.upstream_id == '31001'
    assert identifiers(pdz) == {
        (Namespace.UNIPROT, 'O14745'),
        (Namespace.NAME, 'NHERF-1'),
    }
    assert identifiers(ligand) == {
        (Namespace.UNIPROT, 'P13569'),
        (Namespace.NAME, 'CFTR'),
    }
    assert (pdz.namespace, pdz.identifier) == ('uniprot', 'O14745')
    assert {pdz.entity_type, ligand.entity_type} == {'protein'}
    assert annotations(relation) == [
        ('relation', 'publications', 'PMID:9613608'),
    ]


def test_isoform_accession_is_kept_as_written():
    _, pdz, ligand = extract(ISOFORM_ROW)
    assert (Namespace.UNIPROT, 'P29323-2') in identifiers(ligand)
    assert ligand.molecular_form['isoform_identifier'] == {
        'ns': 'uniprot',
        'id': 'P29323-2',
    }
    # An accession without an isoform number asserts no isoform.
    assert pdz.molecular_form is None


def test_organism_label_is_not_a_taxon_of_the_proteins():
    relation, pdz, ligand = extract(MIXED_SPECIES_ROW)
    assert pdz.taxon is None and ligand.taxon is None
    assert all(a['term'] != 'in_taxon' for a in relation.annotations)


def test_pdz_domain_number_is_not_turned_into_a_feature_or_a_description():
    relation, pdz, _ = extract(MIXED_SPECIES_ROW)
    assert pdz.molecular_form is None
    assert {a['term'] for a in relation.annotations} == {'publications'}
    # Two records that differ only in the domain number are the same statement.
    other, _, _ = extract(OTHER_DOMAIN_ROW)
    assert other.relation_key == relation.relation_key


@pytest.mark.parametrize('role', ['pdz', 'ligand'])
@pytest.mark.parametrize('value', ['', 'NHERF-1', 'O14745;P13569'])
def test_no_relation_without_an_accession_of_each_protein(role, value):
    mapper = resource.datasets()['interactions'].mapper
    assert mapper(dict(HUMAN_ROW, **{f'{role}_accession': value})) is None


def test_reference_that_is_not_a_pubmed_id_gives_no_publication():
    relation, _, _ = extract(dict(HUMAN_ROW, reference='unpublished'))
    assert annotations(relation) == []


@pytest.mark.parametrize('number', ['0123456', '123456789'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(number):
    relation, _, _ = extract(dict(HUMAN_ROW, reference=number))
    assert annotations(relation) == []


def test_row_without_a_record_link_is_still_an_interaction():
    relation, _, _ = extract(dict(HUMAN_ROW, interaction_id=None))
    assert relation.upstream_id == 'fixture:1'
    assert annotations(relation) == [
        ('relation', 'publications', 'PMID:9613608'),
    ]


def cell(content):
    return (
        '<td width="80" class="textblock" bgcolor="#CCFFCC" valign="bottom">\n'
        f'{content}\n</td>'
    )


def label(text):
    return cell(f'<a class="textblockbold">\n<b>\n</b> {text} </a>')


def accession(value):
    return cell(
        f'<a target="_blank" href="http://www.uniprot.org/uniprot/{value}">'
        f'{value}</a>'
    )


def record(name, pdz, domain, ligand_name, ligand, organism, pmid, number):
    reference = cell(
        '<a target="_blank" href="http://www.ncbi.nlm.nih.gov:80/entrez/'
        'query.fcgi?cmd=Retrieve&amp;db=pubmed&amp;dopt=Abstract&amp;'
        f'list_uids={pmid}">{pmid}</a>'
    )
    return (
        '<tr align="center">\n'
        + cell(f'<a class="textblockbold">{name}</a>')
        + accession(pdz) + label(domain) + label(ligand_name)
        + accession(ligand) + label(organism) + reference
        + '<td align="left" valign="middle">\n'
        f'<input type="checkbox" name="interaction-{number}" /> \n</td>\n'
        '<td align="left" valign="middle">\n'
        f'<a href="interaction_detail/{number}"><img src="details.gif"/></a>\n'
        '</td>\n</tr>'
    )


HEADERS = (
    'PDZ name', 'PDZ accession code', 'PDZ domain', 'Ligand name',
    'Ligand accession code', 'Organism', 'Reference',
)


def page(*records, headers=HEADERS):
    # The layout of the PDZBase page: the interactions table sits in a form
    # inside two layout tables, after a table of page furniture.
    header = '<tr align="center">\n' + ''.join(
        f'<td width="80" class="textblock" valign="bottom">\n<b>{name}</b>\n</td>'
        for name in headers
    ) + '</tr>'
    return (
        '<html><head><title>PDZBase web service</title></head><body>\n'
        '<TABLE BORDER="0" WIDTH="750"><TR>'
        '<TD class="heading" colspan="2">Interactions</TD></TR></TABLE>\n'
        '<table width="750"><tr><td>\n<table bgcolor="#ffffcc"><tr>'
        '<td align="left">\n<form method="post" action="alignment">\n'
        '<table width="750" cellspacing="2" cellpadding="2">\n'
        + header + '\n' + '\n'.join(records)
        + '\n</table>\n</form>\n</td></tr></table>\n</td></tr></table>\n'
        '</body></html>'
    )


def opener(text):
    return SimpleNamespace(result=text)


def test_parser_reads_the_interactions_table_of_the_page():
    html = page(
        record('NHERF-1', 'O14745', '1', 'CFTR', 'P13569', 'human', '9613608', '31001'),
        record('AF-6', 'P55196', '1', 'EPHB2', 'P29323-2', 'human', '9707552', '30584'),
    )
    assert list(iter_pdzbase(opener(html))) == [HUMAN_ROW, ISOFORM_ROW]


def test_parser_reads_bytes_and_line_iterators():
    html = page(
        record('NHERF-1', 'O14745', '1', 'CFTR', 'P13569', 'human', '9613608', '31001'),
    )
    assert list(iter_pdzbase(opener(html.encode('utf-8')))) == [HUMAN_ROW]
    assert list(iter_pdzbase(opener(iter(html.splitlines(keepends=True))))) == [
        HUMAN_ROW
    ]


def test_parser_rejects_a_page_without_the_interactions_table():
    html = page(
        record('NHERF-1', 'O14745', '1', 'CFTR', 'P13569', 'human', '9613608', '31001'),
        headers=('PDZ name', 'PDZ accession code', 'Ligand name'),
    )
    with pytest.raises(ValueError, match='no table with the columns'):
        list(iter_pdzbase(opener(html)))
    with pytest.raises(ValueError, match='no table with the columns'):
        list(iter_pdzbase(opener('<html><body>Not Found</body></html>')))


def test_parser_rejects_an_interactions_table_without_rows():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_pdzbase(opener(page())))


def test_parser_rejects_a_row_with_missing_cells():
    html = page('<tr><td>NHERF-1</td><td>O14745</td></tr>')
    with pytest.raises(ValueError, match='row 1 has 2 cells'):
        list(iter_pdzbase(opener(html)))


def test_html_table_reader_separates_nested_tables_and_keeps_links():
    html = (
        '<table id="outer"><tr><th>Name</th><td>before '
        '<table id="inner"><tr><td><a href="/x">X</a> and<br>'
        '<a href="/y">Y</a></td></tr></table> after</td></tr>'
        '<tr></tr></table><script>var t = "<table>";</script>'
    )
    outer, inner = iter_html_tables(html)
    assert outer['attrs'] == {'id': 'outer'}
    assert [[c['text'] for c in r] for r in outer['rows']] == [
        ['Name', 'before after']
    ]
    assert [c['header'] for c in outer['rows'][0]] == [True, False]
    (inner_cell,) = inner['rows'][0]
    assert inner_cell['text'] == 'X and Y'
    assert inner_cell['links'] == ['/x', '/y']
