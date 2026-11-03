"""Source contract of the HPMR receptor interaction input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2 import hpmr as mod
from pypath.inputs_v2.parsers.base import iter_html_tables
from pypath.inputs_v2.parsers.hpmr import (
    iter_interactions,
    parse_receptor_page,
)

GENE_PAGE = 'http://www.receptome.org/SearchDB/getGenePage.asp'
PUBMED = (
    'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi'
    '?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids='
)


def partner(role, key, text):
    return (
        f"<span class='IntRow'><a href='{GENE_PAGE}?Param={key}&ProtId=1"
        f"&ProtType={role}' target='_self' title='Click to Open Details Page'>"
        f'{role}</a> {text}</span>'
    )


def non_gene(name):
    return f"<span class='IntRow'><u>non-gene ligand:</u>  {name}</span>"


def comment(number, body):
    return (
        "<span style='margin-left: 10px' class='IntRow'>"
        f"<a id='OpenComm_{number}' href='#' onClick='show()'>[read more..]</a>"
        f"<div id='Comm_{number}' style='display: none'> {body} </div></span>"
    )


def reference(pmid, label):
    return (
        "<a target='_blank' title='click to open reference in new window' "
        f"href='{PUBMED}{pmid}'>{label}</a>"
    )


def entry(*parts):
    joined = '&nbsp;<b>and</b>&nbsp;'.join(parts[:-1]) + parts[-1]
    return f"<TR BgColor='#EFEFEF' ><td>{joined}</td></TR>"


NO_COMMENT = "<span style='margin-left: 10px' class='IntRow'></span>"


def page(symbol, title, entrez, *entries):
    """A receptor page with the markup of the site, cut to the parts read."""
    table = (
        "<TABLE width='100%' class='ints_table' border=0 cellpadding='2'>"
        "<TR BgColor='#FFFFFF'><td></td></TR>" + ''.join(entries) + '</TABLE>'
        if entries
        else 'No known Receptor-Ligand interactions'
    )
    return f"""<html><head><title>HPMR - Receptor</title></head><body>
<div class="paragraph" style="border-top:none;margin-top: 1px;">
<span class="supertitle">{symbol}  -  {title}</span>
<div style="margin-top:5px"><span class="title">Synonyms:</span>KIAA1568,SAX3
<span class="title">External links:</span>
<a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=Retrieve&dopt=Graphics&list_uids={entrez}" target="_blank">Entrez Gene</a>,
<a href="http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=602431" target="_blank">Omim</a>
</div></div>
<div class="paragraph">
<span class="supertitle">
Receptor {symbol} interacts with:</span>
        <div id="GeneInts" style="display:block"><div>{table}<input type='hidden' id='intslist' value=', '></input></div>
</div></div>
</body></html>"""


# The interaction entries of four receptor pages, as the site wrote them.
ROBO2 = page(
    'ROBO2', 'ROUNDABOUT, DROSOPHILA, HOMOLOG OF, 2', '6092',
    entry(
        partner('Ligand', '9353', 'SLIT2 - SLIT, DROSOPHILA, HOMOLOG OF, 2'),
        partner('Receptor', '6382', 'SDC1 - SYNDECAN 1'),
        comment(
            1491,
            'SLIT2-mediated ROBO2 signaling restricts kidney induction to a '
            'single site ' + reference('15130495', '(Grieshammer et al. 2004)')
            + '.Syndecan is a necessary component of Slit/Robo signaling and '
            'is required in the Slit target cells '
            + reference('14761655', '(Steigemann et al. 2004)') + '.',
        ),
    ),
)
SDC1 = page(
    'SDC1', 'SYNDECAN 1', '6382',
    entry(
        partner('Ligand', '9353', 'SLIT2 - SLIT, DROSOPHILA, HOMOLOG OF, 2'),
        partner('Receptor', '10047201', 'ROBO2 - ROUNDABOUT, DROSOPHILA, HOM ...'),
        comment(
            1491,
            'SLIT2-mediated ROBO2 signaling restricts kidney induction to a '
            'single site ' + reference('15130495', '(Grieshammer et al. 2004)'),
        ),
    ),
)
GPR63 = page(
    'GPR63', 'G PROTEIN-COUPLED RECEPTOR 63', '81491',
    entry(
        non_gene('Sphingosine_1-phosphate'),
        non_gene('dioleoylphosphatidic_acid'),
        comment(
            1977,
            'Sphingosine 1-phosphate and dioleoylphosphatidic acid are low '
            'affinity agonists for the orphan receptor GPR63.\n'
            + reference('12618218', 'Cell Signal. 2003 '),
        ),
    ),
)
AVPR1B = page(
    'AVPR1B', 'ARGININE VASOPRESSIN RECEPTOR 1B', '553',
    entry(partner('Ligand', '551', 'AVP - ARG-VASOPRESSIN'), NO_COMMENT),
)
NO_INTERACTIONS = page(
    'P2RY10', 'PUTATIVE PURINERGIC RECEPTOR P2Y10', '27334'
)
ERROR_REPORT = (
    '<!DOCTYPE html><html lang="en"><head><title>TemplateDoesNotExist at '
    '/SearchDB/getGenePage.asp</title></head><body><div id="summary">'
    '<h1>TemplateDoesNotExist at /SearchDB/getGenePage.asp</h1></div>'
    '</body></html>'
)


def test_page_yields_one_record_for_each_partner_of_an_entry():
    ligand, receptor = parse_receptor_page(ROBO2)
    assert ligand == {
        'receptor': 'ROBO2',
        'receptor_title': 'ROUNDABOUT, DROSOPHILA, HOMOLOG OF, 2',
        'receptor_entrez': '6092',
        'interaction_id': '1491',
        'partner_role': 'Ligand',
        'partner': 'SLIT2',
        'partner_title': 'SLIT, DROSOPHILA, HOMOLOG OF, 2',
        'partner_page_id': '9353',
        'entry': ['Ligand SLIT2', 'Receptor SDC1'],
        'comment': (
            'SLIT2-mediated ROBO2 signaling restricts kidney induction to a '
            'single site (Grieshammer et al. 2004) .Syndecan is a necessary '
            'component of Slit/Robo signaling and is required in the Slit '
            'target cells (Steigemann et al. 2004) .'
        ),
        'references': ['15130495', '14761655'],
    }
    assert (receptor['partner_role'], receptor['partner']) == ('Receptor', 'SDC1')
    assert receptor['partner_page_id'] == '6382'
    assert receptor['references'] == ligand['references']


def test_page_reads_non_gene_ligands_and_entries_without_comment():
    first, second = parse_receptor_page(GPR63)
    assert [(r['partner_role'], r['partner']) for r in (first, second)] == [
        ('non-gene ligand', 'Sphingosine_1-phosphate'),
        ('non-gene ligand', 'dioleoylphosphatidic_acid'),
    ]
    assert first['partner_page_id'] is None
    assert first['references'] == ['12618218']
    (bare,) = parse_receptor_page(AVPR1B)
    assert (bare['partner'], bare['partner_title']) == ('AVP', 'ARG-VASOPRESSIN')
    assert bare['interaction_id'] is None
    assert bare['comment'] is None
    assert bare['references'] == []


def test_page_without_interactions_and_error_report_are_told_apart():
    assert parse_receptor_page(NO_INTERACTIONS) == []
    assert parse_receptor_page(ERROR_REPORT) is None
    with pytest.raises(ValueError, match='no interaction section'):
        parse_receptor_page('<html><head><title>Search</title></head></html>')


def test_page_with_an_unknown_layout_is_rejected():
    changed = ROBO2.replace('>Ligand</a>', '>Agonist</a>')
    with pytest.raises(ValueError, match='Unexpected link'):
        parse_receptor_page(changed)


def archive(**pages):
    return SimpleNamespace(
        result={
            name: io.BytesIO(text.encode('utf-8')) for name, text in pages.items()
        }
    )


def pages():
    return archive(
        **{
            'gene_10047201-1-Receptor.html': ROBO2,
            # A second protein record of the same gene repeats its entries.
            'gene_10047201-2-Receptor.html': ROBO2,
            'gene_6382-1-Receptor.html': SDC1,
            'gene_13540557-1-Receptor.html': GPR63,
            'gene_10092633-1-Receptor.html': NO_INTERACTIONS,
            'gene_1100166-1-Receptor.html': ERROR_REPORT,
        }
    )


def test_archive_yields_each_ligand_pair_once():
    records = list(iter_interactions(pages(), roles=mod.LIGAND_ROLES))
    assert [(r['partner'], r['receptor'], r['page']) for r in records] == [
        ('SLIT2', 'ROBO2', 'gene_10047201-1-Receptor.html'),
        ('Sphingosine_1-phosphate', 'GPR63', 'gene_13540557-1-Receptor.html'),
        ('dioleoylphosphatidic_acid', 'GPR63', 'gene_13540557-1-Receptor.html'),
        ('SLIT2', 'SDC1', 'gene_6382-1-Receptor.html'),
    ]


def test_archive_yields_the_two_views_of_a_receptor_pair_once():
    records = list(iter_interactions(pages(), roles=mod.RECEPTOR_ROLES))
    assert [(r['partner'], r['receptor']) for r in records] == [('SDC1', 'ROBO2')]


def test_archive_without_receptor_pages_or_interactions_is_rejected():
    with pytest.raises(ValueError, match='no receptor page'):
        list(iter_interactions(archive(**{'gene_1-1-Receptor.html': ERROR_REPORT})))
    with pytest.raises(ValueError, match='no interaction'):
        list(
            iter_interactions(archive(**{'gene_1-1-Receptor.html': NO_INTERACTIONS}))
        )


def extract(dataset, source_row):
    record = mod.resource.datasets()[dataset].mapper(source_row)
    assert record is not None
    out = SilverExtractor('hpmr', dataset)
    out.process_record(record, source_row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return [(i['ns'], i['id']) for i in entity.identifiers]


def annotations(relation, scope='relation'):
    return sorted(
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def rows(html, roles):
    return list(iter_interactions(archive(**{'gene_1-1-Receptor.html': html}), roles=roles))


def test_ligand_pair_is_ligand_to_receptor_with_publications_and_comment():
    (row,) = rows(ROBO2, mod.LIGAND_ROLES)
    relation, ligand, receptor = extract('interactions', row)
    assert relation.predicate == 'interacts_with'
    assert identifiers(ligand) == [(Namespace.GENESYMBOL, 'SLIT2')]
    # The extractor orders identifiers by namespace, the gene ID first.
    assert identifiers(receptor) == [
        (Namespace.ENTREZ, '6092'),
        (Namespace.GENESYMBOL, 'ROBO2'),
    ]
    assert {ligand.entity_type, receptor.entity_type} == {'protein'}
    assert {ligand.taxon, receptor.taxon} == {'9606'}
    assert annotations(relation) == [
        ('description', row['comment']),
        ('publications', 'PMID:14761655'),
        ('publications', 'PMID:15130495'),
    ]


def test_ligand_dataset_gives_ligand_and_receptor_roles():
    (row,) = rows(ROBO2, mod.LIGAND_ROLES)
    relation, ligand, receptor = extract('interactions', row)
    # The page marks the partner as a ligand and is the page of a receptor.
    assert row['partner_role'] == 'Ligand'
    assert (Namespace.GENESYMBOL, 'SLIT2') in identifiers(ligand)
    assert (Namespace.GENESYMBOL, 'ROBO2') in identifiers(receptor)
    assert [v for t, v in annotations(relation, 'subject') if t == PARTICIPANT_ROLE] == ['ligand']
    assert [v for t, v in annotations(relation, 'object') if t == PARTICIPANT_ROLE] == ['receptor']


def test_receptor_receptor_pair_has_no_ligand_or_receptor_roles():
    # The partner of the receptor dataset is marked a `Receptor`: the page does
    # not call either side a ligand, so no role is stated for either.
    (row,) = rows(ROBO2, mod.RECEPTOR_ROLES)
    relation, partner_, receptor = extract('receptor_interactions', row)
    assert row['partner_role'] == 'Receptor'
    assert (Namespace.GENESYMBOL, 'SDC1') in identifiers(partner_)
    assert all(a['term'] != PARTICIPANT_ROLE for a in relation.annotations)


def test_entry_with_a_co_receptor_asserts_no_complex():
    # HPMR lists SLIT2 and SDC1 together on the ROBO2 page without saying how
    # the three are arranged: each is a partner of ROBO2 and nothing more.
    (row,) = rows(ROBO2, mod.LIGAND_ROLES)
    record = mod.map_hpmr_interaction(row)
    out = SilverExtractor('hpmr', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    assert [r.predicate for r in out.relations] == ['interacts_with']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert 'SDC1' not in {e.identifier for e in out.entities.values()}
    assert {a['term'] for a in out.relations[0].annotations} == {
        'description',
        'in_taxon',
        'publications',
        PARTICIPANT_ROLE,
    }
    assert {
        a['term']
        for a in out.relations[0].annotations
        if a['scope'] == 'relation'
    } == {'description', 'publications'}


def test_receptor_pair_is_its_own_dataset():
    (row,) = rows(ROBO2, mod.RECEPTOR_ROLES)
    relation, partner_, receptor = extract('receptor_interactions', row)
    assert relation.predicate == 'interacts_with'
    assert identifiers(partner_) == [(Namespace.GENESYMBOL, 'SDC1')]
    assert (Namespace.GENESYMBOL, 'ROBO2') in identifiers(receptor)
    # Each dataset maps only its own kind of partner.
    assert mod.map_hpmr_interaction(row) is None
    (ligand_row,) = rows(ROBO2, mod.LIGAND_ROLES)
    assert mod.map_hpmr_receptor_interaction(ligand_row) is None


def test_non_gene_ligand_is_a_chemical_without_species():
    first, second = rows(GPR63, mod.LIGAND_ROLES)
    relation, ligand, receptor = extract('interactions', first)
    assert ligand.entity_type == 'chemical_entity'
    assert identifiers(ligand) == [(Namespace.NAME, 'Sphingosine 1-phosphate')]
    assert ligand.taxon is None
    # The page lists it under the ligands, so it is a ligand without species.
    assert annotations(relation, 'subject') == [(PARTICIPANT_ROLE, 'ligand')]
    assert (Namespace.GENESYMBOL, 'GPR63') in identifiers(receptor)
    assert ('publications', 'PMID:12618218') in annotations(relation)
    _, ligand, _ = extract('interactions', second)
    assert identifiers(ligand) == [(Namespace.NAME, 'dioleoylphosphatidic acid')]


def test_stimulus_or_unnamed_non_gene_ligand_is_not_a_chemical():
    (row,) = rows(AVPR1B, mod.LIGAND_ROLES)
    non_gene_row = dict(row, partner_role='non-gene ligand')
    assert mod.map_hpmr_interaction(dict(non_gene_row, partner='light')) is None
    assert mod.map_hpmr_interaction(dict(non_gene_row, partner='ssRNA')) is None
    assert mod.map_hpmr_interaction(dict(non_gene_row, partner=None)) is None
    with pytest.raises(ValueError, match='non-gene ligand'):
        mod.map_hpmr_interaction(dict(non_gene_row, partner='caffeine'))


def test_enzyme_named_in_place_of_its_product_is_not_a_ligand():
    (row,) = rows(AVPR1B, mod.LIGAND_ROLES)
    # HPMR gives tryptophan hydroxylase as the ligand of serotonin receptors.
    enzyme = dict(row, receptor='HTR2A', receptor_entrez='3356', partner='TPH1')
    assert mod.map_hpmr_interaction(enzyme) is None
    assert mod.map_hpmr_interaction(dict(row, partner='DUSP18')) is None
    # The same gene as a receptor's receptor partner is not affected.
    assert mod.map_hpmr_receptor_interaction(
        dict(enzyme, partner_role='Receptor')
    ) is not None


def test_pair_without_comment_or_reference_has_no_annotation():
    (row,) = rows(AVPR1B, mod.LIGAND_ROLES)
    relation, ligand, receptor = extract('interactions', row)
    assert annotations(relation) == []
    # The cut-off partner title is not an identifier.
    assert identifiers(ligand) == [(Namespace.GENESYMBOL, 'AVP')]


def test_reference_written_with_stray_characters_is_cleaned():
    (row,) = rows(AVPR1B, mod.LIGAND_ROLES)
    references = [' 12810718', 'PMID: 16284174', '9832610)', '12802337.', '9832610']
    relation, _, _ = extract('interactions', dict(row, references=references))
    assert annotations(relation) == [
        ('publications', 'PMID:12802337'),
        ('publications', 'PMID:12810718'),
        ('publications', 'PMID:16284174'),
        ('publications', 'PMID:9832610'),
    ]
    # The references are the ``list_uids`` of PubMed links, so a short number is
    # a PubMed ID (the shared rule is 1 to 8 digits); text, a leading zero and
    # nine digits are not.
    relation, _, _ = extract('interactions', dict(row, references=['95', 'n/a']))
    assert annotations(relation) == [('publications', 'PMID:95')]
    relation, _, _ = extract(
        'interactions', dict(row, references=['0123456', '123456789', 'n/a'])
    )
    assert annotations(relation) == []


def test_no_relation_without_a_participant():
    (row,) = rows(AVPR1B, mod.LIGAND_ROLES)
    assert mod.map_hpmr_interaction(dict(row, partner='')) is None
    assert mod.map_hpmr_interaction(dict(row, receptor=None)) is None


def test_general_reader_keeps_cell_content_in_order_and_nested_tables_apart():
    html = (
        "<table id='outer'><tr><th>Name</th><th colspan='2'>Link</th></tr>"
        "<tr><td>alpha <b>beta</b></td>"
        "<td><a href='/x' title='t'>first</a> between <a href='/y'>second</a>"
        "<table class='inner'><tr><td>nested</td></tr></table></td></tr>"
        '<tr></tr></table>'
    )
    outer, inner = iter_html_tables(html)
    assert outer['attrs'] == {'id': 'outer'}
    header, body = outer['rows']
    assert [(c['header'], c['text']) for c in header] == [
        (True, 'Name'),
        (True, 'Link'),
    ]
    assert header[1]['attrs'] == {'colspan': '2'}
    assert body[0]['content'] == ['alpha ', 'beta']
    assert body[1]['text'] == 'first between second'
    assert body[1]['links'] == ['/x', '/y']
    assert body[1]['content'] == [
        {'href': '/x', 'text': 'first', 'attrs': {'href': '/x', 'title': 't'}},
        ' between ',
        {'href': '/y', 'text': 'second', 'attrs': {'href': '/y'}},
    ]
    assert inner['attrs'] == {'class': 'inner'}
    assert [[c['text'] for c in row] for row in inner['rows']] == [['nested']]
