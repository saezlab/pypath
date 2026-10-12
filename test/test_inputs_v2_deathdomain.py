"""DeathDomain database pages: table reading and Biolink mapping."""

import functools
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor

from pypath.inputs_v2 import deathdomain as mod
from pypath.inputs_v2.parsers import deathdomain as parser

# Header lines of the four tables, as the pages write them.
HEADERS = {
    'In vitro interaction': (
        '<tr><th colspan="14" class="title"><i>In vitro</i> interaction</th></tr>'
        '<tr><th colspan="5">&nbsp;</th><th colspan="3">DD1</th><th colspan="3">DD2</th>'
        '<th rowspan="2">Reference</th></tr>'
        '<tr><th>Family</th><th>DD1</th><th></th><th>DD2</th><th>Method</th><th>Species</th>'
        '<th>Region</th><th>Expression</th><th>Species</th><th>Region</th><th>Expression</th></tr>'
    ),
    'In vivo interaction': (
        '<tr><th colspan="14" class="title"><i>In vivo</i> interaction </th></tr>'
        '<tr><th colspan="5">&nbsp;</th><th rowspan="2">Endogenous<br>expression</th>'
        '<th rowspan="2">Overexpression</th><th colspan="2">DD1</th><th colspan="2">DD2</th>'
        '<th rowspan="2">Reference</th></tr>'
        '<tr><th>Family</th><th>DD1</th><th></th><th>DD2</th><th>Method</th><th>Species</th>'
        '<th>Region</th><th>Species</th><th>Region</th></tr>'
    ),
    'Characterization': (
        '<tr><th colspan="12" class="title">Characterization</th></tr>'
        '<tr><th colspan="5">&nbsp;</th><th rowspan="2">Stoichiometry</th><th rowspan="2">Affinity</th>'
        '<th colspan="2">Protein-Protein interface</th><th rowspan="2">Reference</th></tr>'
        '<tr><th>Family</th><th>DD1</th><th></th><th>DD2</th><th>Method</th><th>Mutation</th>'
        '<th>Complex structure</th></tr>'
    ),
    'Functional Role': (
        '<tr><th colspan="12" class="title">Functional Role</th></tr>'
        '<tr><th>Family</th><th>DD1</th><th></th><th>DD2</th><th>Method</th><th>Death-related</th>'
        '<th>Death-unrelated</th><th>Reference</th></tr>'
    ),
}
# Lines of the CARD page (capture of 29 April 2016).
LINES = {
    'In vitro interaction': [
        ['CARD', 'Apaf1', 'Apaf1', 'In vitro purification protein assembly(Size-exclusion chromatography)',
         'Not specified', 'Not specified', 'Not specified', 'Not specified', 'Not specified', 'Not specified',
         '10944114'],
        ['CARD', 'Apaf1', 'Caspase9', 'GST fusion protein pull-down',
         'Not specified', '1-97', 'E.coli', 'Not specified', 'Not specified', 'In vitro translation', '9390557'],
        ['CARD', 'Caspase9', 'Apaf1', 'GST fusion protein pull-down',
         'Not specified', 'Not specified', 'In vitro translation', 'Not specified', '1-97', 'E.coli', '9390557'],
        ['CARD', 'Caspase9', 'Caspase9', 'In vitro purification protein assembly(2D gel electrophoresis)',
         'U937 lysate', '', '', 'U937 lysate', '', '', '14747474'],
    ],
    'In vivo interaction': [
        ['CARD', 'Bcl-10', 'CARD11', 'Co-immunoprecipitation', '', 'HEK293',
         'Mouse', 'Not specified', 'Human', '1-470', '16809782'],
        ['CARD', 'CARD11', 'Bcl-10', 'Co-immunoprecipitation', '', 'HEK293',
         'Human', '1-470', 'Mouse', 'Not specified', '16809782'],
        ['CARD', 'ASC', 'Caspase1', 'Co-immunoprecipitation', '', 'HEK293',
         'Human', '105-195', 'Not specified', 'Not specified', '11967258'],
    ],
    'Characterization': [
        ['CARD', 'Apaf1', 'Caspase9', 'Structure', '1:1', '',
         'D27A, E40A (Apaf1) R13A, R56A (Caspase9)', '3YGS', '10376594'],
    ],
    'Functional Role': [
        ['CARD', 'Apaf1', 'Apaf1', 'Caspase cleaved assay', 'Caspase9 activation', '', '9837928'],
    ],
}


def _line(cells):
    family, first, second, *rest, reference = cells
    return (
        f'<tr><td class="c0" align="center"><a href="/summaries/{family}">{family}</a></td>'
        f'<td class="c0"><a href="/proteins/show?name={first}&amp;family={family}">{first}</a></td>'
        f'<td class="c0"><a href="/interactions/{family}?p1={first}&amp;p2={second}" '
        'title="see this interaction only"><img alt="Link" src="/images/link.gif" /></a></td>'
        f'<td class="c0"><a href="/proteins/show?name={second}&amp;family={family}">{second}</a></td>'
        + ''.join(f'<td class="c1">{value}</td>' for value in rest)
        + f'<td class="c0"><a href="/references/show/{reference}" title="">{reference}</a></td></tr>'
    )


def _page(lines=LINES, headers=HEADERS):
    return (
        '<html><body><div id="loginpart"><table cellpadding="0"><tr><td>Login</td></tr></table></div>'
        + ''.join(
            f'<table width="100%"><tr><td>Sort by</td></tr></table>'
            f'<table width="100%" class="tab"><thead>{headers[title]}</thead><tbody>'
            + ''.join(map(_line, lines[title]))
            + '</tbody></table>'
            for title in headers
        )
        + '</body></html>'
    )


def observations():
    return list(parser.mark_second_listings(parser.iter_family_page(_page(), 'CARD')))


def extract(rows):
    out = SilverExtractor('deathdomain', 'interactions')
    for number, row in enumerate(rows, 1):
        record = mod.map_deathdomain_observation(row)
        assert record is not None
        out.process_record(record, row, f'fixture:{number}', number)
    return out


def annotations(relation, term, scope='relation'):
    return [a['value'] for a in relation.annotations if a['term'] == term and a['scope'] == scope]


def test_each_table_line_is_one_row_with_its_own_columns():
    rows = observations()
    assert [r['table'] for r in rows] == ['in_vitro'] * 4 + ['in_vivo'] * 3 + [
        'characterization', 'functional_role']
    pull_down = rows[1]
    assert pull_down == {
        'table': 'in_vitro', 'page_family': 'CARD', 'family': 'CARD',
        'protein_a': 'Apaf1', 'protein_b': 'Caspase9', 'method': 'GST fusion protein pull-down',
        'a_species': 'Not specified', 'a_region': '1-97', 'a_expression': 'E.coli',
        'b_species': 'Not specified', 'b_region': 'Not specified', 'b_expression': 'In vitro translation',
        'reference': '9390557', 'mirrored': False,
    }
    assert rows[4]['overexpression'] == 'HEK293' and rows[4]['a_species'] == 'Mouse'
    assert rows[7]['complex_structure'] == '3YGS' and rows[7]['stoichiometry'] == '1:1'
    assert rows[8]['death_related'] == 'Caspase9 activation'


def test_page_with_another_layout_is_an_error():
    with pytest.raises(ValueError, match='lacks the tables'):
        list(parser.iter_family_page('<html><body>404 Error: Page Not Found</body></html>', 'CARD'))
    short = dict(LINES, **{'Functional Role': [['CARD', 'Apaf1', 'Apaf1', 'Caspase cleaved assay', '9837928']]})
    with pytest.raises(ValueError, match='has 6 cells, expected 8'):
        list(parser.iter_family_page(_page(lines=short), 'CARD'))
    renamed = dict(HEADERS, **{'Functional Role': HEADERS['Functional Role'].replace('>DD2<', '>Partner<')})
    with pytest.raises(ValueError, match='unexpected header'):
        list(parser.iter_family_page(_page(headers=renamed), 'CARD'))
    retitled = dict(HEADERS, **{'Functional Role': HEADERS['Functional Role'].replace('Functional Role', 'Disease')})
    with pytest.raises(ValueError, match="unknown table 'Disease'"):
        list(parser.iter_family_page(_page(headers=retitled), 'CARD'))


def test_second_listing_of_an_observation_is_marked_once():
    rows = observations()
    assert [r['mirrored'] for r in rows] == [False, False, True, False, False, True, False, False, False]
    assert mod.deathdomain_exclusion(rows[2]) == mod.EXCLUDED_MIRRORED
    assert mod.map_deathdomain_observation(rows[2]) is None
    # A homomeric observation is listed once and kept.
    assert mod.deathdomain_exclusion(rows[0]) is None
    # Two identical lines are two observations, not a pair of listings.
    twice = list(parser.mark_second_listings([dict(rows[1], mirrored=None)] * 2))
    assert [r['mirrored'] for r in twice] == [False, False]


def test_interaction_keeps_names_method_and_publication():
    rows = observations()
    out = extract([rows[1]])
    (relation,) = out.relations
    assert relation.predicate == 'interacts_with'
    assert annotations(relation, 'has_evidence_of_type') == ['MI:0096']
    assert annotations(relation, 'publications') == ['PMID:9390557']
    # The tables give the database's protein names and no accession.
    assert {(e.entity_type, str(e.namespace), e.identifier, e.taxon) for e in out.entities.values()} == {
        ('protein', 'name', 'Apaf1', None), ('protein', 'name', 'Caspase9', None)}
    assert all(e.identity_scope == 'deathdomain' for e in out.entities.values())


@pytest.mark.parametrize('reference', ['0123456', '123456789', '0'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(reference):
    (relation,) = extract([dict(observations()[1], reference=reference)]).relations
    assert relation.predicate == 'interacts_with'
    assert annotations(relation, 'publications') == []


def test_species_is_stated_only_where_the_line_names_one():
    rows = observations()
    out = extract([rows[4], rows[6], rows[3]])
    cross_species, one_species, lysate = out.relations
    assert annotations(cross_species, 'in_taxon', 'subject') == ['NCBITaxon:10090']
    assert annotations(cross_species, 'in_taxon', 'object') == ['NCBITaxon:9606']
    assert annotations(one_species, 'in_taxon', 'subject') == ['NCBITaxon:9606']
    assert annotations(one_species, 'in_taxon', 'object') == []
    # "U937 lysate" is where the protein came from, not a species label.
    assert annotations(lysate, 'in_taxon', 'subject') == []
    assert lysate.subject_entity_key == lysate.object_entity_key
    assert annotations(lysate, 'has_evidence_of_type') == []
    with pytest.raises(ValueError, match='Unknown DeathDomain species'):
        mod.map_deathdomain_observation(dict(rows[4], a_species='Hamster'))


def test_functional_role_and_unknown_method_are_not_mapped_to_interactions():
    rows = observations()
    assert mod.deathdomain_exclusion(rows[8]) == mod.EXCLUDED_FUNCTIONAL_ROLE
    assert mod.map_deathdomain_observation(rows[8]) is None
    structure = mod.map_deathdomain_observation(rows[7])
    assert not [a for a in structure.annotations if a.term is mod.slots.has_evidence_of_type]
    with pytest.raises(ValueError, match='Unknown DeathDomain method'):
        mod.map_deathdomain_observation(dict(rows[1], method='Luciferase assay'))
    assert mod.map_deathdomain_observation(dict(rows[1], protein_b='')) is None


def test_parser_reads_the_four_family_pages(monkeypatch):
    opened = []

    def fake_download(family):
        def open_(force_refresh=False):
            opened.append((family, force_refresh))
            return SimpleNamespace(result=_page().replace('CARD', family))
        return SimpleNamespace(open=open_)

    monkeypatch.setattr(
        mod.resource.interactions,
        '_raw_parser',
        functools.partial(
            parser.iter_deathdomain_observations,
            family_downloads={f: fake_download(f) for f in mod.FAMILY_CAPTURES},
        ),
    )
    rows = list(mod.resource.interactions.raw(force_refresh=True))
    assert opened == [('CARD', True), ('DD', True), ('DED', True), ('PYD', True)]
    assert len(rows) == 36 and {r['page_family'] for r in rows} == set(mod.FAMILY_CAPTURES)
    assert sum(mod.deathdomain_exclusion(r) is None for r in rows) == 24


def test_dataset_reads_the_four_pages_from_its_downloads():
    raw_parser = mod.resource.interactions._raw_parser
    assert raw_parser.func is parser.iter_deathdomain_observations
    assert raw_parser.keywords == {'family_downloads': mod._family_downloads}
    assert list(mod._family_downloads) == ['CARD', 'DD', 'DED', 'PYD']
    with pytest.raises(ValueError, match='downloads of the pages'):
        list(parser.iter_deathdomain_observations(None))


@pytest.mark.parametrize(
    ('method', 'term'),
    [
        ('Co-immunoprecipitation', 'MI:0019'),
        ('In vitro purification protein assembly(immunoprecipitation)', 'MI:0019'),
        ('Yeast two-hybrid', 'MI:0018'),
        ('Two-hybrid', 'MI:0018'),
        ('GST fusion protein pull-down', 'MI:0096'),
        ('His fusion protein pull-down', 'MI:0096'),
        ('His Fusion protein pull-down', 'MI:0096'),
        ('In vitro purification protein assembly(Pull down)', 'MI:0096'),
        ('Affinity purification-mass spectrometry', 'MI:0004'),
        ('Fluorescence resonance energy transfer', 'MI:0055'),
        ('Bimolecular fluorescence complementation', 'MI:0809'),
        ('Co-Immunofluorescence staining', 'MI:0403'),
        ('Size-exclusion chromatography', 'MI:0071'),
        ('Gel filtration', 'MI:0071'),
        ('In vitro purification protein assembly(Size-exclusion chromatography)', 'MI:0071'),
        ('Surface plasmon resonance', 'MI:0107'),
        ('In vitro purification protein assembly(Surface plasmon resonance)', 'MI:0107'),
        ('Isothermal titration calorimety', 'MI:0065'),
        ('In vitro purification protein assembly(Isothermal titration calorimety)', 'MI:0065'),
        ('In vitro purification protein assembly(chemical crosslinking)', 'MI:0030'),
        ('In vitro purification protein assembly(Native gel-electrophoresis)', None),
        ('In vitro purification protein assembly(2D gel electrophoresis)', None),
        ('Structure', None),
        ('Mass spectrometry, MALS', None),
        ('Fluorescence probe labeling', None),
        ('', None),
    ],
)
def test_every_method_label_has_its_own_psi_mi_term(method, term):
    assert mod.METHOD_TERMS[method] == term
    relation = mod.map_deathdomain_observation(dict(observations()[1], method=method))
    values = [
        a.value for a in relation.annotations
        if a.term is mod.slots.has_evidence_of_type
    ]
    assert values == ([term] if term else [])


def test_the_method_crosswalk_has_no_label_outside_this_list():
    assert len(mod.METHOD_TERMS) == 26
