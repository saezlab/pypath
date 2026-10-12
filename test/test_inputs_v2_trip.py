"""TRIP Database: shared HTML/CSV readers, parsed rows and Biolink mapping."""

import functools
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import trip as mod
from pypath.inputs_v2.parsers import trip as parser
from pypath.inputs_v2.parsers.base import (
    iter_csv,
    iter_html_links,
    iter_html_tables,
    read_text,
)

# The five header lines and lines of the 6 August 2015 export, as written.
HEADER = '''\
,,,Screening,,,,,,Validation,,,,,,,,,,,,Characterization,,,,,,Functional consequence,,,Reference,
,,,Experimental screening,,,,,Non-experimental screening,In vitro validation,,,,,,In vivo validation,,,,,,Binding region mapping,,,,Stoichiometry ,Affinity (Kd),Post-translational modification,Subcellular trafficking,Activity,,
,,,TRP channel construct,,Interactor source,,,,Assay with recombinant proteins,,,,,,Assay with endogenous proteins,Assay with overexpressed proteins,,,,,,,,,,,,,,,
,,,,,,,,,TRP channel construct,,,Interactor construct,,,Cell or tissue,Cell or tissue,TRP channel construct,,Interactor construct,,TRP channel,,Interactor,,,,,,,,
TRP channels,Interactors,Methods,Species,Region,Species,Organ/tissue ,Sample type,,Species,Region,Expression system,Species,Region,Expression system,,,Species,Region,Species,Region,Species,Region,Species,Region,,,,,,,Overlap
'''
LINES = '''\
TRPP1,PKD1/Polycystin-1,Yeast two-hybrid,,,,,,,,,,,,,,,,,,,Human,742-968,Human,1227-1339,,,,,,9192675,
TRPC3,TRPC1,Co-immunoprecipitation,,,,,,,,,,,,,,HEK293T,Human,Full-length,Human,Full-length,,,,,,,,,,9215637,O
TRPC1,TRPC3,Co-immunoprecipitation,,,,,,,,,,,,,,HEK293T,Human,Full-length,Human,Full-length,,,,,,,,,,9215637,O
TRPC1,TRPC3,Inference,,,,,,Prediction,,,,,,,,,,,,,,,,,,,,,,9215637,O
TRPC3,IP3R1,Co-immunoprecipitation,,,,,,,,,,,,,,HEK293,Human,Full-length,Rat,Full-length,,,,,,,,,,10518223,
TRPC3,IP3R1,Patch clamp,,,,,,,,,,,,,,,,,,,,,,,,,,,Activation,10518223,
TRPC3,TrkB,Co-immunoprecipitation,,,,,,,,,,,,,Rat brain lysates,,,,,,,,,,,,,,,10677043,
TRPC4,NHERF1,Cell surface biotinylation,,,,,,,,,,,,,,,,,,,,,,,,,,Increase in plasma membrane level,,12154080,
TRPV1,PKA,In vitro PTM assay,,,,,,,,,,,,,,,,,,,,,,,,,Phosphorylation (Ser-116),,,12194871,
TRPC3,PLCƣ1,Yeast two-hybrid,,,,,,,Human,1-340,Yeast,Rat,750-1290,Yeast,,,,,,,,,,,,,,,,12437926,
TRPC6,Calmodulin,Fusion protein-pull down assay,,,,,,,Mouse,828-883,E. coli,Not specified,Full-length,In vitro translation,,,,,,,,,,,,,,,,12379180,
TRPC4,TRPC1,Patch clamp,,,,,,,,,,,,,,,,,,,,,,,,,,,New channel creation,11301024,O
TRPC1,TRPC4,Patch clamp,,,,,,,,,,,,,,,,,,,,,,,,,,,New channel creation,11301024,O
'''
EXPORT = SimpleNamespace(result='﻿' + HEADER + LINES)


def _protein_row(name, *accessions):
    cells = ''.join(
        f'<td align="center"> <a href=\'http://www.uniprot.org/uniprot/{a}\' target=\'_blank\'>'
        f'{a}<img src="/assets/external.gif" alt="External" /></a> </td>'
        if a
        else '<td align="center">  </td>'
        for a in accessions
    )
    return (
        f'<tr><td align="center" id="1"><a href="/proteins/show?id={name}">{name}</a></td>'
        f'{cells}</tr>'
    )


def _protein_table(*rows):
    return (
        '<table width="100%" class="tab"><thead><tr><th rowspan="2">Name</th>'
        '<th colspan="3">UniprotKB ID</th></tr><tr><th>Human</th><th>Mouse</th><th>Rat</th></tr>'
        f'</thead><tbody>{"".join(rows)}</tbody></table>'
    )


# Rows of the protein list as captured on 18 April 2021. The later site lost
# the accessions of names with Greek letters.
PROTEIN_LIST = (
    '<html><body><div id="loginpart"><table><tr><td>Login</td></tr></table></div>'
    + _protein_table(
        _protein_row('TRPC1', 'P48995', 'Q61056', 'Q9QX01'),
        _protein_row('TRPC2', '', 'Q9R244', 'Q9R283'),
        _protein_row('TRPC3', 'Q13507', 'Q9QZC1', 'Q9JMI9'),
        _protein_row('TRPC6', 'Q9Y210', 'Q61143', 'Q99N78'),
        _protein_row('TRPP1', 'Q13563', 'O35245', 'D4A5F1'),
    )
    + _protein_table(
        _protein_row('CaMKIIƢ', '', '', ''),
        _protein_row('Calmodulin', 'P62158', 'P62204', 'P62161'),
        _protein_row('DHC', 'Q14204', 'Q9JHU4?', 'P38650'),
        _protein_row('IP3R1', 'Q14643', 'P11881', 'P29994'),
        _protein_row('PKD1/Polycystin-1', 'P98161', 'O08852', 'Q9ERV0'),
    )
    + '</body></html>'
)
ACCESSIONS = parser.protein_accessions(PROTEIN_LIST)


def observations():
    return list(parser._observations(parser.iter_trip_table(EXPORT), ACCESSIONS))


def extract(rows):
    out = SilverExtractor('trip', 'interactions')
    for row in rows:
        record = mod.map_trip_observation(row)
        assert record is not None
        out.process_record(record, row, f'fixture:{row["line"]}', int(row['line']))
    return out


def annotations(relation, term, scope='relation'):
    return [a['value'] for a in relation.annotations if a['term'] == term and a['scope'] == scope]


def test_html_tables_reader_keeps_cells_links_and_nested_tables_apart():
    html = (
        '<table class="outer"><tr><th colspan="2">Head</th></tr>'
        '<tr><td>a<br>b &amp; c<table id="inner"><tr><td>nested</td></tr></table></td>'
        '<td title="t"><a href="/x?p=1&amp;q=2">link</a><script>ignored()</script>'
        '<tr><td>unclosed</table><p>outside</p>'
    )
    outer, inner = iter_html_tables(html)
    assert outer['attrs'] == {'class': 'outer'}
    assert [[c['text'] for c in row] for row in outer['rows']] == [
        ['Head'], ['a b & c', 'link'], ['unclosed'],
    ]
    assert outer['rows'][0][0]['header'] and outer['rows'][0][0]['attrs'] == {'colspan': '2'}
    assert outer['rows'][1][1]['links'] == ['/x?p=1&q=2']
    assert outer['rows'][1][1]['attrs'] == {'title': 't'}
    assert inner['attrs'] == {'id': 'inner'} and inner['rows'][0][0]['text'] == 'nested'
    assert list(iter_html_tables('<p>no table</p>')) == []


def test_html_links_reader():
    links = list(iter_html_links('<a name="x">anchor</a><a href="a?b=1">  two\n words </a><a href="c">'))
    assert [(link['href'], link['text']) for link in links] == [('a?b=1', 'two words'), ('c', '')]


def test_delimited_rows_reader_uses_given_names_and_rejects_other_layouts():
    opener = SimpleNamespace(result=b'\xef\xbb\xbfh1,h2\n\nA,"b,c"\n')
    assert list(iter_csv(opener, fieldnames=('x', 'y'), skip_rows=1, strict=True)) == [{'x': 'A', 'y': 'b,c'}]
    assert list(iter_csv(SimpleNamespace(result=''), fieldnames=('x',), strict=True)) == []
    with pytest.raises(ValueError, match='Line 1 has 1 fields'):
        list(iter_csv(SimpleNamespace(result='<html>429 Too Many Requests</html>\n'), fieldnames=('x', 'y'), strict=True))
    assert read_text(SimpleNamespace(result={'member': b'text'})) == 'text'
    assert read_text(None) == ''


def test_protein_list_gives_accessions_by_name_and_species():
    assert ACCESSIONS['TRPC3'] == {'9606': 'Q13507', '10090': 'Q9QZC1', '10116': 'Q9JMI9'}
    # No human entry for TRPC2, and "Q9JHU4?" is not an accession.
    assert ACCESSIONS['TRPC2'] == {'10090': 'Q9R244', '10116': 'Q9R283'}
    assert ACCESSIONS['DHC'] == {'9606': 'Q14204', '10116': 'P38650'}
    assert ACCESSIONS['CaMKIIβ'] == {}
    with pytest.raises(ValueError, match='no UniprotKB ID table'):
        parser.protein_accessions('<html><body>Too Many Requests</body></html>')


def test_one_row_per_line_with_category_species_and_accession():
    rows = observations()
    assert len(rows) == 13 and [r['line'] for r in rows][:2] == ['6', '7']
    first = rows[0]
    assert first['categories'] == ['characterization']
    assert (first['channel_name'], first['channel_species'], first['channel_uniprot']) == (
        'TRPP1', 'Human', 'Q13563')
    assert (first['interactor_name'], first['interactor_taxon'], first['interactor_uniprot']) == (
        'PKD1/Polycystin-1', '9606', 'P98161')
    assert first['binding_channel_region'] == '742-968'
    cross_species = rows[4]
    assert (cross_species['channel_uniprot'], cross_species['interactor_uniprot']) == ('Q13507', 'P29994')
    assert rows[5]['categories'] == ['functional_consequence']
    assert (rows[8]['ptm_type'], rows[8]['ptm_residue']) == ('Phosphorylation', 'Ser-116')
    # The export writes Greek letters 0x210 code points too low.
    assert rows[9]['interactor'] == 'PLCƣ1' and rows[9]['interactor_name'] == 'PLCγ1'
    assert parser.repair_greek('CaMKIIƢ') == 'CaMKIIβ'
    assert parser.repair_greek('Catenin-Ƥ2 TRPV1') == 'Catenin-δ2 TRPV1'


def test_changed_header_is_an_error():
    moved = SimpleNamespace(result=HEADER.replace('Functional consequence', 'Function') + LINES)
    with pytest.raises(ValueError, match="'Functional consequence' was expected"):
        list(parser.iter_trip_table(moved))
    with pytest.raises(ValueError, match='fields'):
        list(parser.iter_trip_table(SimpleNamespace(result='<html>503</html>')))


def test_species_is_never_guessed():
    rows = observations()
    endogenous = rows[6]  # co-immunoprecipitation from rat brain lysates
    assert endogenous['channel_species'] is None and endogenous['interactor_species'] is None
    assert endogenous['channel_uniprot'] is None and endogenous['interactor_uniprot'] is None
    out = extract([endogenous])
    assert {(str(e.namespace), e.identifier, e.taxon) for e in out.entities.values()} == {
        ('name', 'TRPC3', None), ('name', 'TrkB', None)}
    unspecified = rows[10]  # interactor species "Not specified"
    assert (unspecified['channel_uniprot'], unspecified['interactor_uniprot']) == ('Q61143', None)
    with pytest.raises(ValueError, match='Unknown TRIP species'):
        parser._species(dict(rows[0], binding_channel_species='Hamster'), 'channel')
    assert parser._species(dict(rows[0], binding_interactor_species='Xenopus'), 'interactor') == 'Xenopus'
    assert parser.SPECIES_TAXA['Xenopus'] is None


@pytest.mark.parametrize('number', ['0123456', '123456789'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(number):
    (relation,) = extract([dict(observations()[0], reference=number)]).relations
    assert annotations(relation, 'publications') == []


def test_physical_interaction_keeps_method_publication_and_source_accessions():
    out = extract(observations()[:1])
    (relation,) = out.relations
    assert relation.predicate == 'interacts_with'
    assert annotations(relation, 'has_evidence_of_type') == ['MI:0018']
    assert annotations(relation, 'publications') == ['PMID:9192675']
    assert annotations(relation, 'in_taxon', 'subject') == ['NCBITaxon:9606']
    assert not [a for a in relation.annotations if 'qualifier' in a['term']]
    entities = {e.identifier: e for e in out.entities.values()}
    assert set(entities) == {'P98161', 'Q13563'}
    assert entities['P98161'].namespace == Namespace.UNIPROT and entities['P98161'].taxon == '9606'
    assert {(i['ns'], i['id']) for i in entities['Q13563'].identifiers} == {
        ('uniprot', 'Q13563'), ('name', 'TRPP1')}
    # The interactor is the subject, the TRP channel the object.
    assert out.entities[relation.subject_entity_key].identifier == 'P98161'


def test_predictions_and_second_listings_are_excluded():
    rows = observations()
    assert [mod.trip_exclusion(r) for r in rows[1:4]] == [
        None, mod.EXCLUDED_MIRRORED, mod.EXCLUDED_NON_EXPERIMENTAL]
    assert mod.map_trip_observation(rows[2]) is None
    assert mod.map_trip_observation(rows[3]) is None
    out = extract([rows[1]])
    assert [r.predicate for r in out.relations] == ['interacts_with']
    with pytest.raises(ValueError, match='non-experimental screening'):
        mod.trip_exclusion(dict(rows[0], screening_non_experimental='Guess'))
    # Functional consequences are directed: both listings are kept.
    assert [r['mirrored'] for r in rows[11:13]] == [False, False]


def test_functional_consequences_use_an_exact_crosswalk():
    rows = observations()
    out = extract([rows[5], rows[7], rows[8], rows[11]])
    activation, trafficking, phosphorylation, new_channel = out.relations
    assert [r.predicate for r in out.relations] == ['affects'] * 4
    assert annotations(activation, 'object_aspect_qualifier') == ['activity']
    assert annotations(activation, 'object_direction_qualifier') == ['increased']
    # A change of membrane level names a compartment: no direction is stated.
    assert annotations(trafficking, 'object_aspect_qualifier') == ['localization']
    assert annotations(trafficking, 'object_direction_qualifier') == []
    assert annotations(phosphorylation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(phosphorylation, 'object_direction_qualifier') == ['increased']
    assert not [a for a in new_channel.annotations if 'qualifier' in a['term']]
    # Functional assays have no PSI-MI interaction detection term.
    assert annotations(activation, 'has_evidence_of_type') == []
    assert mod.functional_effect(dict(rows[5], activity='Inhibition')) == (
        mod._ASPECT.activity, mod._DIRECTION.decreased)
    assert mod.functional_effect(dict(rows[5], activity='Sensitization')) == (mod._ASPECT.activity, None)
    assert mod.functional_effect(dict(rows[8], ptm_type='Dephosphorylation'))[1] == mod._DIRECTION.decreased
    assert mod.functional_effect(dict(rows[8], ptm_type='Glycosylation')) == (mod._ASPECT.glycosylation, None)
    assert mod.functional_effect(rows[0]) is None


def test_unknown_values_raise_instead_of_getting_a_default():
    rows = observations()
    with pytest.raises(ValueError, match='Unknown TRIP activity'):
        mod.map_trip_observation(dict(rows[5], activity='Modulation'))
    with pytest.raises(ValueError, match='Unknown TRIP trafficking'):
        mod.map_trip_observation(dict(rows[7], trafficking='Moved'))
    with pytest.raises(ValueError, match='Unknown TRIP modification'):
        mod.map_trip_observation(dict(rows[8], ptm='Sumoylation', ptm_type='Sumoylation'))
    with pytest.raises(ValueError, match='Unknown TRIP method'):
        mod.map_trip_observation(dict(rows[0], method='Telepathy'))
    assert mod.map_trip_observation(dict(rows[0], interactor_name='')) is None


def test_parser_joins_the_protein_list_and_passes_force_refresh(monkeypatch):
    opened = []

    def open_(force_refresh=False):
        opened.append(force_refresh)
        return SimpleNamespace(result=PROTEIN_LIST)

    monkeypatch.setattr(
        mod.resource.interactions,
        '_raw_parser',
        functools.partial(
            parser.iter_trip_observations,
            proteins_download=SimpleNamespace(open=open_),
        ),
    )
    monkeypatch.setattr(mod.resource.interactions, 'download', SimpleNamespace(open=lambda **_: EXPORT))
    rows = list(mod.resource.interactions.raw(force_refresh=True))
    assert opened == [True]
    assert len(rows) == 13 and rows[0]['interactor_uniprot'] == 'P98161'
    assert sum(mod.map_trip_observation(r) is not None for r in rows) == 11


def test_dataset_reads_the_protein_list_from_its_download():
    raw_parser = mod.resource.interactions._raw_parser
    assert raw_parser.func is parser.iter_trip_observations
    assert raw_parser.keywords == {'proteins_download': mod._proteins_download}
    with pytest.raises(ValueError, match='download of the protein list'):
        list(parser.iter_trip_observations(EXPORT))


@pytest.mark.parametrize(
    ('column', 'value', 'aspect', 'direction'),
    [
        ('activity', 'Activation', 'activity', 'increased'),
        ('activity', 'Inhibition', 'activity', 'decreased'),
        ('activity', 'Sensitization', 'activity', None),
        ('activity', 'Desensitization', 'activity', None),
        ('activity', 'New channel creation', None, None),
        ('trafficking', 'Increase in plasma membrane level', 'localization', None),
        ('trafficking', 'Decrease in plasma membrane level', 'localization', None),
        ('trafficking', 'Increase in lysosomal membrane level', 'localization', None),
        ('trafficking', 'Retain in the endoplasmic reticulum', 'localization', None),
        ('trafficking', 'Internalization from membrane by ligand', 'localization', None),
        ('ptm_type', 'Phosphorylation', 'phosphorylation', 'increased'),
        ('ptm_type', 'Phosphorylayion', 'phosphorylation', 'increased'),
        ('ptm_type', 'Dephosphorylation', 'phosphorylation', 'decreased'),
        ('ptm_type', 'Ubiquitination', 'ubiquitination', 'increased'),
        ('ptm_type', 'De-ubiquitination', 'ubiquitination', 'decreased'),
        ('ptm_type', 'Glycosylation', 'glycosylation', None),
    ],
)
def test_every_functional_consequence_value_has_its_own_aspect_and_direction(
    column, value, aspect, direction
):
    row = {column: value}
    if column == 'ptm_type':
        row['ptm'] = value
    assert mod.functional_effect(row) == (
        mod._ASPECT[aspect] if aspect else None,
        mod._DIRECTION[direction] if direction else None,
    )


# Further functional consequence lines of the export (a glycosylation, a
# dephosphorylation, and a phosphorylation without a residue).
PTM_LINES = '''\
TRPV5,Klotho,In vitro PTM assay,,,,,,,,,,,,,,,,,,,,,,,,,Glycosylation (Asn-358),,,16239475,
TRPV6,PTP1B,In vivo PTM assay,,,,,,,,,,,,,,,,,,,,,,,,,Dephosphorylation (Tyr-162),,,17197020,
TRPV4,SRC,In vitro PTM assay,,,,,,,,,,,,,,,,,,,,,,,,,Phosphorylation,,,12538589,
'''


def test_modified_residue_is_a_feature_of_the_channel():
    rows = observations()
    extra = list(parser._observations(parser.iter_trip_table(SimpleNamespace(result=HEADER + PTM_LINES)), ACCESSIONS))
    pka = rows[8]
    klotho, ptp1b, src = extra
    assert (pka['ptm_amino_acid'], pka['ptm_position']) == ('Ser', 116)
    out = extract([pka])
    channel = next(e for e in out.entities.values() if e.label == 'TRPV1' or 'TRPV1' in str(e.identifiers))
    interactor = next(e for e in out.entities.values() if e is not channel)
    assert channel.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'S',
            'position': 116,
            'end_position': 116,
            'coordinate_reference': {
                'identifier': None,
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': None,
        }
    ]
    # The enzyme is not modified.
    assert interactor.molecular_form is None
    # The relation is unchanged: the interactor phosphorylates the channel.
    assert annotations(out.relations[0], 'object_aspect_qualifier') == ['phosphorylation']
    # A dephosphorylation acts on a phosphorylated residue.
    (site,) = mod.channel_form(ptp1b)['modifications']
    assert (site['term'], site['residue'], site['position']) == ('phosphorylated residue', 'Y', 162)
    # PSI-MOD has no glycosylated residue under this name: the type is kept.
    (site,) = mod.channel_form(klotho)['modifications']
    assert (site['term'], site['residue'], site['position']) == (None, 'N', 358)
    assert site['description'] == 'Glycosylation'
    # A modification without a residue states no site.
    assert src['ptm'] == 'Phosphorylation' and mod.channel_form(src) is None
    assert mod.map_trip_observation(src).object.molecular_form is None
    # Other lines have no site.
    assert mod.map_trip_observation(rows[5]).object.molecular_form is None


def test_unknown_amino_acid_raises_and_a_zero_position_is_no_site():
    pka = observations()[8]
    with pytest.raises(ValueError, match='Unknown TRIP amino acid'):
        mod.map_trip_observation(dict(pka, ptm_amino_acid='Xaa', ptm_residue='Xaa-116'))
    assert mod.channel_form(dict(pka, ptm_position=0)) is None
