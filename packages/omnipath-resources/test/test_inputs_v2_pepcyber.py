"""Source contract of the PepCyber:P~Pep interaction input."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import pepcyber
from pypath.inputs_v2.parsers import pepcyber as parser

HEADER = (
    '<tr><td><a href="javascript:checkAll();">All</a> <a>None</a></td>'
    '<td><a><strong>PPBD class</strong></a></td>'
    '<td><a><strong>PPBD protein</strong></a></td>'
    '<td><a><strong>Substrate protein</strong></a></td>'
    '<td><a><strong>Binding sequence</strong></a></td>'
    '<td><a><strong>Index site</strong></a></td>'
    '<td><a><strong>Evidence</strong></a></td>'
    '<td><a><strong>#Records</strong></a></td>'
    '<td><a><strong>Category</strong></a></td>'
    '<td>Details</td></tr>'
)
ENTREZ = (
    '<a href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve'
    '&dopt=full_report&list_uids={gene_id} " target="_blank">{symbol}</a>'
)


def search_row(number, domain, ppbd, substrate, sequence, site, evidence, records, category):
    link = (
        f'<a href="ppbd.php?s={ppbd}">{ppbd}</a>' if ppbd else '<a href="ppbd.php?s="></a>'
    )
    return (
        f'<tr bgcolor="#FFFFFF"><td><a name="{number}"></a>'
        f'<input type="checkbox" value="{number}"/></td>'
        f'<td>{domain}</td><td>{link}</td><td>{substrate}</td>'
        f'<td><span class="sequence">{sequence}</span></td><td>{site}</td>'
        f'<td>{evidence}</td><td>{records}</td><td>{category}</td>'
        f'<td><a href="javascript:call_url(\'idetail.php?iid={number}\',\'ul_{number}\');">'
        'Details</a> <a href="#top">Top</a></td></tr>'
        '<tr bgcolor="#FFFFFF"><td colspan="13"><div id="ul_%s"></div></td></tr>' % number
    )


# Rows of the search result page.
BAD = search_row(
    1, '14-3-3', 'LOC284100', ENTREZ.format(gene_id=12015, symbol='Bad'),
    'SPFRGRSRpSAPPNLWAA', 136, 'Inference', 1, 'A',
)
# The PPBD protein has no name.
NO_PPBD = search_row(
    1105, '14-3-3', '', ENTREZ.format(gene_id=443531, symbol='AA-NAT'),
    'SPGRQRRHpTLPANEFRC', 31, 'Structure', 1, 'A',
)
# The substrate has no link: its name is free text, here a viral protein with a
# detail page that gives its UniProt accession.
HCV = search_row(
    297, '14-3-3', 'YWHAE', '<a>HCVgp1</a>', 'GVRATRKTpSERSQPRGR', 53,
    'Inference,mutagenesis', 2, 'A',
)
# No Gene ID, and a name that is a character-set artifact: no protein.
NOT_FOUND = search_row(
    10078, 'SH2', 'PIK3R1', '<a>Î´ÕÒµ½</a>', '', 475,
    'Mutagenesis', 1, 'C',
)
# No peptide and no site: category B and C are without one.
NO_PEPTIDE = search_row(
    5282, 'SH2', 'PTPN11', ENTREZ.format(gene_id=10818, symbol='FRS2'),
    '', 'N/A', ',Inference', 3, 'B',
)
# Two phosphorylated residues in the peptide, one index site.
TWO_SITES = search_row(
    10001, 'SH2', 'PTPN11', ENTREZ.format(gene_id=887, symbol='CCKBR'),
    "LS<font color='#770000'>pY</font>TTI<font color='#770000'>pY</font>TTISTLGP",
    438, 'Mutagenesis', 1, 'C',
)
# A site, but a peptide without the phosphorylated residue marked.
UNMARKED = search_row(
    11, '14-3-3', 'SFN', ENTREZ.format(gene_id=7157, symbol='TP53'),
    'HLKSKKGQSTSRHKKLM', 376, 'Inference', 2, 'A',
)
# Interaction 270: the only evidence is a structure.
STRUCTURE = search_row(
    270, '14-3-3', 'YWHAB', ENTREZ.format(gene_id=879837, symbol='exoS'),
    'SQGLLDALDLASKPERS', 427, 'Structure', 1, 'A',
)
ROWS = (BAD, NO_PPBD, HCV, NOT_FOUND, NO_PEPTIDE, TWO_SITES, UNMARKED)


def search_page(*rows, header=HEADER):
    return (
        '<html><body><table><tr><td>PepCyber</td></tr></table>'
        f'<table><tr><td>9530 interactions found.</td></tr>{header}{"".join(rows)}'
        '</table></body></html>'
    )


def detail_page(ppbd, substrate, records):
    def protein(symbol, refseq, swissprot, description='', link='ppbd.php?s='):
        return (
            '<table><tr><td>Gene name</td><td>'
            f'<a href="{link}">{symbol}</a>'
            f'{description}</td></tr>'
            '<tr><td>Synonyms</td><td></td></tr>'
            '<tr><td>GO annotations</td><td>GO_component:</td></tr>'
            f'<tr><td>RefSeq accession</td><td><a>{refseq}</a></td></tr>'
            f'<tr><td>SwissProt accession</td><td><a>{swissprot}</a></td></tr></table>'
        )

    rows = ''.join(
        f'<tr><td>{category}</td><td>{evidence}</td><td><a>{citation}[<a href='
        f'"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed'
        f'&list_uids={pmid}">{pmid}</a>]</a></td></tr>'
        for category, evidence, citation, pmid in records
    )
    return (
        '<table><tr><td>Record details:</td></tr></table>'
        + protein(*ppbd)
        + protein(*substrate)
        + '<table><tr><td>Category</td><td>Evidence</td><td>Reference</td></tr>'
        + rows
        + '</table>'
    )


# Detail pages of interactions 1 and 297.
DETAILS = {
    1: detail_page(
        ('LOC284100', 'XP_375443', 'A2RUB5'),
        (
            'Bad', 'NP_031548', 'Q61337', '( Description: Bcl-associated death promoter. )',
            'http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&list_uids=12015',
        ),
        [
            (
                'A', 'Inference',
                'Yano S, Tokumitsu H, Soderling TR., Calcium promotes cell survival., ',
                9859994,
            )
        ],
    ),
    297: detail_page(
        ('YWHAE', 'NP_006752', 'P62258'),
        ('HCVgp1', 'NP_751919', 'P26664'),
        [
            ('A', 'Mutagenesis', 'Aoki H, Hepatitis C virus core protein., ', 10644344),
            ('A', 'Inference', 'Lee SK, Interaction of HCV core protein., ', 17141194),
        ],
    ),
}


class Opened:
    def __init__(self, text):
        self.result = text.encode('latin-1')


class Details:
    """The Downloads of the detail directory and of a detail page."""

    def __init__(self, pages):
        self.pages = pages
        self.opened = []

    def open(self, force_refresh=False, number=None):
        if number is None:
            links = ''.join(f'<a href="pepcyber_{n}.html">pepcyber_{n}.html</a>\n' for n in self.pages)
            return Opened(f'<html><body><h1>Index of /pepcyber_details/</h1><pre>{links}</pre>')
        self.opened.append(number)
        return Opened(self.pages[number])


def read(*rows, pages=DETAILS, header=HEADER):
    details = Details(pages)
    parsed = list(
        parser.iter_pepcyber(
            Opened(search_page(*rows, header=header)),
            details=details,
            detail_numbers=details,
        )
    )
    return parsed, details


def extract(parsed_row):
    record = pepcyber.resource.datasets()['interactions'].mapper(parsed_row)
    assert record is not None
    out = SilverExtractor('pepcyber', 'interactions')
    out.process_record(record, parsed_row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return [(i['ns'], i['id']) for i in entity.identifiers]


def attributes(relation):
    return sorted((a['term'], a['value']) for a in relation.annotations)


def test_search_rows_are_read_with_their_fields():
    parsed, _ = read(*ROWS)
    assert [row['number'] for row in parsed] == [1, 1105, 297, 10078, 5282, 10001, 11]
    first = parsed[0]
    assert {k: v for k, v in first.items() if k != 'detail'} == {
        'ppbd_class': '14-3-3',
        'ppbd_symbol': 'LOC284100',
        'substrate_symbol': 'Bad',
        'binding_sequence': 'SPFRGRSRpSAPPNLWAA',
        'index_site': 136,
        'evidence': ['Inference'],
        'n_records': 1,
        'category': 'A',
        'number': 1,
        'substrate_entrez': '12015',
        'phospho_residue': 'S',
    }
    by_number = {row['number']: row for row in parsed}
    # The source writes ",Inference": the empty entry is dropped.
    assert by_number[5282]['evidence'] == ['Inference']
    assert by_number[5282]['index_site'] is None
    assert by_number[5282]['binding_sequence'] is None
    assert by_number[1105]['ppbd_symbol'] is None
    assert by_number[297]['substrate_entrez'] is None
    # A peptide with two phosphorylated residues or none marked gives no residue.
    assert by_number[10001]['phospho_residue'] is None
    assert by_number[11]['phospho_residue'] is None


def test_only_the_rescued_detail_pages_are_read():
    parsed, details = read(*ROWS)
    assert sorted(details.opened) == [1, 297]
    by_number = {row['number']: row for row in parsed}
    assert by_number[11]['detail'] is None
    assert by_number[1]['detail'] == {
        'ppbd': {
            'symbol': 'LOC284100', 'description': None, 'entrez': None,
            'refseq': 'XP_375443', 'swissprot': 'A2RUB5',
        },
        'substrate': {
            'symbol': 'Bad', 'description': 'Bcl-associated death promoter.', 'entrez': '12015',
            'refseq': 'NP_031548', 'swissprot': 'Q61337',
        },
        'records': [
            {
                'category': 'A',
                'evidence': 'Inference',
                'pmid': '9859994',
                'reference': (
                    'Yano S, Tokumitsu H, Soderling TR., Calcium promotes cell survival., [9859994]'
                ),
            }
        ],
    }
    assert [r['pmid'] for r in by_number[297]['detail']['records']] == ['10644344', '17141194']


def test_reader_rejects_a_changed_page():
    with pytest.raises(ValueError, match='no table with the columns'):
        read(BAD, header='<tr><td>x</td></tr>')
    with pytest.raises(ValueError, match='listed twice'):
        read(BAD, BAD)
    bad_site = BAD.replace('>136<', '>Y136<')
    with pytest.raises(ValueError, match='index site'):
        read(bad_site)
    with pytest.raises(ValueError, match='protein table lacks'):
        parser.read_detail('<table><tr><td>Gene name</td><td>x</td></tr></table>')
    with pytest.raises(ValueError, match='lacks the two protein tables'):
        parser.read_detail('<table><tr><td>Category</td></tr></table>')
    with pytest.raises(ValueError, match='lists no pepcyber'):
        list(parser.iter_detail_numbers(Opened('<html><a href="x.html">x</a></html>')))


def test_interaction_is_ppbd_protein_to_phosphorylated_substrate():
    parsed, _ = read(BAD)
    relation, subject, obj = extract(parsed[0])
    # Category A (mediated by a PPBD) with inference as the evidence: the source
    # does not say that the two proteins touch.
    assert relation.predicate == 'interacts_with'
    assert subject.entity_type == obj.entity_type == 'protein'
    assert set(identifiers(subject)) == {
        (Namespace.UNIPROT, 'A2RUB5'),
        (Namespace.GENESYMBOL, 'LOC284100'),
        (Namespace.REFSEQ_PROTEIN, 'XP_375443'),
    }
    assert set(identifiers(obj)) == {
        (Namespace.UNIPROT, 'Q61337'),
        (Namespace.ENTREZ, '12015'),
        (Namespace.REFSEQ_PROTEIN, 'NP_031548'),
        (Namespace.GENESYMBOL, 'Bad'),
    }
    # No species: the source does not say, and its accessions are not all human.
    assert subject.taxon is None and obj.taxon is None
    assert obj.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'S',
            'position': 136,
            'end_position': 136,
            'coordinate_reference': {
                'identifier': {'ns': 'refseq_protein', 'id': 'NP_031548'},
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': None,
        }
    ]
    assert subject.molecular_form is None
    assert attributes(relation) == [
        ('pepcyber:binding_domain', '14-3-3'),
        ('pepcyber:binding_sequence', 'SPFRGRSRpSAPPNLWAA'),
        ('pepcyber:category', 'A'),
        ('pepcyber:evidence_type', 'Inference'),
        ('publications', 'PMID:9859994'),
    ]


def test_category_does_not_state_a_direct_interaction():
    parsed, _ = read(BAD, NO_PEPTIDE, TWO_SITES, UNMARKED)
    # Categories A and B (mediated by a PPBD) and C (not known) all give the
    # broad predicate when the evidence is not a structure.
    assert [r['category'] for r in parsed] == ['A', 'B', 'C', 'A']
    assert {extract(row)[0].predicate for row in parsed} == {'interacts_with'}


def test_a_structure_among_the_evidence_is_a_direct_interaction():
    parsed, _ = read(STRUCTURE, HCV)
    structure, other = parsed
    assert structure['evidence'] == ['Structure']
    assert extract(structure)[0].predicate == 'directly_physically_interacts_with'
    # The same category without a structure is not.
    assert other['category'] == 'A'
    assert extract(other)[0].predicate == 'interacts_with'


def test_unknown_category_is_an_error_not_a_default():
    (row,) = read(BAD)[0]
    with pytest.raises(ValueError, match='Unknown PepCyber category'):
        pepcyber.resource.interactions.mapper(dict(row, category='D'))


def test_no_site_without_a_marked_residue_and_a_position():
    parsed, _ = read(NO_PEPTIDE, TWO_SITES, UNMARKED)
    for row in parsed:
        _, _, obj = extract(row)
        assert obj.molecular_form is None


def test_substrate_without_a_sequence_in_a_detail_page_has_no_reference_sequence():
    page = search_row(
        12, 'SH2', 'GRB2', ENTREZ.format(gene_id=1956, symbol='EGFR'),
        'TAENAEpYLRVAPQSSEF', 1092, 'Peptide library', 1, 'A',
    )
    parsed, _ = read(page)
    _, _, obj = extract(parsed[0])
    modification = obj.molecular_form['modifications'][0]
    assert (modification['residue'], modification['position']) == ('Y', 1092)
    assert modification['coordinate_reference']['identifier'] is None


def test_substrate_known_only_by_a_detail_accession_is_kept():
    parsed, _ = read(HCV)
    relation, _, obj = extract(parsed[0])
    assert (Namespace.UNIPROT, 'P26664') in identifiers(obj)
    assert [v for t, v in attributes(relation) if t == 'pepcyber:evidence_type'] == [
        'Inference',
        'Mutagenesis',
    ]
    assert [v for t, v in attributes(relation) if t == 'publications'] == [
        'PMID:10644344',
        'PMID:17141194',
    ]


def test_no_relation_without_a_participant():
    mapper = pepcyber.resource.datasets()['interactions'].mapper
    parsed, _ = read(NO_PPBD, NOT_FOUND)
    # No PPBD protein name; a free-text substrate name with no Gene ID or accession.
    assert [mapper(row) for row in parsed] == [None, None]


def test_unknown_category_and_evidence_type_raise():
    mapper = pepcyber.resource.datasets()['interactions'].mapper
    parsed, _ = read(BAD)
    with pytest.raises(ValueError, match='Unknown PepCyber category'):
        mapper({**parsed[0], 'category': 'D'})
    with pytest.raises(ValueError, match='Unknown PepCyber evidence type'):
        mapper({**parsed[0], 'evidence': ['Western blot']})


def test_a_row_without_details_has_no_publication_and_no_invented_species():
    parsed, _ = read(NO_PEPTIDE)
    relation, subject, obj = extract(parsed[0])
    assert not [t for t, _ in attributes(relation) if t == 'publications']
    assert set(identifiers(subject)) == {(Namespace.GENESYMBOL, 'PTPN11')}
    assert set(identifiers(obj)) == {
        (Namespace.ENTREZ, '10818'),
        (Namespace.GENESYMBOL, 'FRS2'),
    }
    assert obj.taxon is None
