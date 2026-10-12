"""NetPath PSI-MI interactions: reader, parsed rows and Biolink mapping."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import netpath as mod
from pypath.inputs_v2.parsers import netpath as parser
from pypath.inputs_v2.parsers.base import iter_psimi_xml

MI = 'xmlns="net:sf:psidev:mi"'

# Interactors as NetPath writes them: (label, full name, alias, HPRD, OMIM,
# sequence database, sequence accession, Entrez Gene, taxon, organism).
TSLP = ('TSLP', 'Thymic stromal lymphopoietin', 'Thymic stromal lymphopoietin isoform 1',
        '06108', '607003', 'Refseq', 'NP_149024.1', '85480', '10090', 'Mouse')
CRLF2 = ('CRLF2', 'Cytokine receptor like factor 2', 'TSLPR',
         '02288', '400023', 'Refseq', 'NP_071431.2', '64109', '10090', 'Mouse')
IL7R = ('IL7R', 'Interleukin 7 receptor', 'IL-7R-alpha',
        '00893', '146661', 'Refseq', 'NP_002176.2', '3575', '10090', 'Mouse')
NUMB = ('NUMB', 'Numb homolog', 'S171',
        '04767', '603728', 'Refseq', 'NP_001005743.1', '8650', '9606', 'Human')
# NetPath labels this UniProt accession "GenBank".
ITCH = ('ITCH', 'Itchy homolog E3 ubiquitin protein ligase', 'dJ468O1.1',
        '07565', '606409', 'GenBank', 'Q96J02', '83737', '9606', 'Human')
LEP = ('LEP', 'Leptin', 'Obese protein',
       '01249', '164160', 'Refseq', 'NP_000221.1', '3952', '9606', 'Human')
LEPR = ('LEPR', 'Leptin receptor', 'Leptin receptor isoform 3',
        '03001', '601007', 'Refseq', 'NP_002294.2', '3953', '9606', 'Human')
IKBKG = ('IKBKG', 'Inhibitor of kappa light polypeptide gene enhancer in B-cells, kinase gamma',
         'NEMO', '02217', '300248', 'Refseq', 'NP_003630.1', '8517', '10090', 'Mouse')

VIVO = ('in vivo', 'MI:0493')
VITRO = ('in vitro', 'MI:0492')
HUMAN = ('9606', 'Human')
MAMMALIA = ('40674', 'Mammalia')


def _term(tag, label, mi):
    return (
        f'<{tag}><names><shortLabel>{label}</shortLabel><fullName>{label}</fullName></names>'
        f'<xref><primaryRef db="psi-mi" dbAc="MI:0488" id="{mi}" refType="identity"/></xref></{tag}>'
    )


def _interactor(number, record):
    label, name, alias, hprd, omim, seq_db, seq_id, entrez, taxon, organism = record
    return (
        f'<interactor id="{number}"><names><shortLabel>{label}</shortLabel>'
        f'<fullName>{name}</fullName><alias type="gene name" typeAc="MI:0301">{alias}</alias></names>'
        f'<xref><primaryRef db="HPRD" dbAc="MI:0468" id="{hprd}" refType="identity"/>'
        f'<secondaryRef db="OMIM" dbAc="MI:0480" id="{omim}"/>'
        f'<secondaryRef db="{seq_db}" id="{seq_id}"/>'
        f'<secondaryRef db="Entrez gene" dbAc="MI:0477" id="{entrez}"/></xref>'
        + _term('interactorType', 'protein', 'MI:0326')
        + f'<organism ncbiTaxId="{taxon}"><names><shortLabel>{organism}</shortLabel></names></organism>'
        '</interactor>'
    )


def _entry(pathway, name, experiments, interactors, participants=None, extra='',
           interaction_type=('physical interaction', 'MI:0218')):
    """One NetPath ``<entry>``; participants are (interactor number, role)."""
    participants = participants or [
        (number, ('neutral', 'MI:0497')) for number in range(1, len(interactors) + 1)
    ]
    return (
        f'<entry><source><names><shortLabel>PandeyLab;Institute of Bioinformatics</shortLabel></names>'
        f'<xref><primaryRef db="NetPath" dbAc="0" id="{pathway}"/></xref></source>'
        '<experimentList>'
        + ''.join(
            f'<experimentDescription id="{number}"><bibref><xref>'
            + ''.join(
                f'<{"primaryRef" if i == 0 else "secondaryRef"} db="pubmed" id="{pmid}"/>'
                for i, pmid in enumerate(pmids)
            )
            + '</xref></bibref><hostOrganismList>'
            f'<hostOrganism ncbiTaxId="{host[0]}"><names><shortLabel>{host[1]}</shortLabel></names></hostOrganism>'
            '</hostOrganismList>' + _term('interactionDetectionMethod', *method)
            + '</experimentDescription>'
            for number, (pmids, host, method) in enumerate(experiments, 1)
        )
        + '</experimentList><interactorList>'
        + ''.join(_interactor(number, record) for number, record in enumerate(interactors, 1))
        + extra
        + f'</interactorList><interactionList><interaction id="1"><names><fullName>{name}</fullName></names>'
        '<experimentList>'
        + ''.join(f'<experimentRef>{n}</experimentRef>' for n in range(1, len(experiments) + 1))
        + '</experimentList><participantList>'
        + ''.join(
            f'<participant id="{i}"><interactorRef>{number}</interactorRef>'
            + (
                '<experimentalRoleList>' + _term('experimentalRole', *role) + '</experimentalRoleList>'
                if role
                else ''
            )
            + '</participant>'
            for i, (number, role) in enumerate(participants, 1)
        )
        + '</participantList>' + _term('interactionType', *interaction_type)
        + '</interaction></interactionList></entry>'
    )


def _file(*entries):
    return f'<?xml version="1.0" encoding="UTF-8"?><entrySet level="2" version="5" {MI}>{"".join(entries)}</entrySet>'.encode()


# NetPath_24 (TSLP), entries 1 and 3.
NETPATH_24 = _file(
    _entry(24, 'Interaction between TSLP and CRLF2',
           [(['10974032', '10881176'], HUMAN, VIVO), (['10974032', '10881176'], HUMAN, VITRO)],
           [TSLP, CRLF2]),
    _entry(24, 'Interaction between TSLP and CRLF2, IL7R',
           [(['10881176'], HUMAN, VIVO)], [TSLP, CRLF2, IL7R]),
)
# NetPath_22 (Leptin), entry 59: the export appends an interactor without
# cross-references, typed "small molecule", to this ligand-receptor interaction.
NETPATH_22 = _file(
    _entry(22, 'Interaction between LEP, LEP and LEPR, LEPR',
           [(['9218470', '15051545'], MAMMALIA, VIVO)], [LEP, LEPR],
           participants=[(1, ('neutral', 'MI:0497')), (2, ('neutral', 'MI:0497')), (3, None)],
           extra='<interactor id="3"><names><shortLabel>LEP</shortLabel><fullName>LEP</fullName></names>'
           + _term('interactorType', 'small molecule', 'MI:0328') + '</interactor>'),
)
# NetPath_3 (Notch), entry 54, and NetPath_21 (RANKL), entry 53.
NETPATH_3 = _file(
    _entry(3, 'Interaction between NUMB and ITCH', [(['12682059'], HUMAN, VIVO)], [NUMB, ITCH],
           participants=[(1, ('bait', 'MI:0496')), (2, ('prey', 'MI:0498'))]),
)
NETPATH_21 = _file(
    _entry(21, 'Interaction between IKBKG and IKBKG, IKBKG',
           [(['19830703'], ('10090', 'Mouse'), VIVO)], [IKBKG]),
)

ARCHIVE = SimpleNamespace(
    result={
        'PSI-MI/': b'',
        'PSI-MI/NetPath_24_PSIMI.xml': NETPATH_24,
        'PSI-MI/NetPath_3_PSIMI.xml': NETPATH_3,
        'PSI-MI/NetPath_22_PSIMI.xml': NETPATH_22,
        'PSI-MI/NetPath_21_PSIMI.xml': NETPATH_21,
    }
)


def extract(dataset, rows):
    out = SilverExtractor('netpath', dataset)
    mapper = mod.resource.datasets()[dataset].mapper
    for number, row in enumerate(rows, 1):
        record = mapper(row)
        assert record is not None
        out.process_record(record, row, f'fixture:{number}', number)
    return out


def annotations(relation, term, scope='relation'):
    return [a['value'] for a in relation.annotations if a['term'] == term and a['scope'] == scope]


def test_psimi_reader_resolves_compact_references():
    first, second = iter_psimi_xml(NETPATH_24)
    assert first['entry_index'] == 1 and second['entry_index'] == 2
    assert first['source']['xrefs'][0]['db'] == 'NetPath'
    assert first['name'] == 'Interaction between TSLP and CRLF2'
    assert [t['mi'] for t in first['interaction_types']] == ['MI:0218']
    assert [e['detection_method']['mi'] for e in first['experiments']] == ['MI:0493', 'MI:0492']
    assert [x['id'] for x in first['experiments'][0]['bibliography']] == ['10974032', '10881176']
    assert first['experiments'][0]['host_organisms'] == [{'taxon': '9606', 'label': 'Human', 'name': None}]
    tslp = first['participants'][0]
    assert tslp['interactor']['label'] == 'TSLP'
    assert tslp['interactor']['organism']['taxon'] == '10090'
    assert tslp['interactor']['type']['mi'] == 'MI:0326'
    assert [(x['db'], x['id'], x['primary']) for x in tslp['interactor']['xrefs']] == [
        ('HPRD', '06108', True), ('OMIM', '607003', False),
        ('Refseq', 'NP_149024.1', False), ('Entrez gene', '85480', False),
    ]
    assert [r['mi'] for r in tslp['experimental_roles']] == ['MI:0497']
    assert len(second['participants']) == 3


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    rows = [r for r in parser.iter_netpath_interactions(ARCHIVE) if r['pathway_id'] == 'NetPath_24']
    row = dict(rows[0], pubmed=['0123456', '123456789', '10974032'])
    out = extract('interactions', [row])
    assert annotations(out.relations[0], 'publications') == ['PMID:10974032']


def test_psimi_reader_reads_the_expanded_form_and_rejects_other_documents():
    expanded = (
        f'<entrySet level="2" version="5" {MI}><entry><interactionList><interaction id="7">'
        '<experimentList><experimentDescription id="1"><bibref><xref>'
        '<primaryRef db="pubmed" id="123"/></xref></bibref></experimentDescription></experimentList>'
        '<participantList><participant id="1"><interactor id="5"><names><shortLabel>A</shortLabel></names>'
        '</interactor></participant></participantList></interaction></interactionList></entry></entrySet>'
    )
    (interaction,) = iter_psimi_xml(expanded)
    assert interaction['id'] == '7'
    assert interaction['experiments'][0]['bibliography'][0]['id'] == '123'
    assert interaction['experiments'][0]['detection_method'] is None
    assert interaction['participants'][0]['interactor']['label'] == 'A'
    assert interaction['participants'][0]['experimental_roles'] == []
    with pytest.raises(ValueError, match='entrySet'):
        list(iter_psimi_xml('<html><body>503</body></html>'))
    dangling = expanded.replace(
        '<interactor id="5"><names><shortLabel>A</shortLabel></names></interactor>',
        '<interactorRef>9</interactorRef>',
    )
    with pytest.raises(ValueError, match='unknown interactor'):
        list(iter_psimi_xml(dangling))


def test_one_row_per_interaction_and_experiment_in_pathway_order():
    rows = list(parser.iter_netpath_observations(ARCHIVE))
    assert [(r['pathway_id'], r['entry_index'], r['experiment_id']) for r in rows] == [
        ('NetPath_3', 1, '1'), ('NetPath_21', 1, '1'), ('NetPath_22', 1, '1'),
        ('NetPath_24', 1, '1'), ('NetPath_24', 1, '2'), ('NetPath_24', 2, '1'),
    ]
    tslp = rows[3]
    assert tslp['pubmed'] == ['10974032', '10881176']
    assert (tslp['detection_method'], tslp['detection_method_mi']) == ('in vivo', 'MI:0493')
    assert (tslp['host_taxon'], tslp['host_organism']) == ('9606', 'Human')
    assert tslp['participants'][0]['entrez'] == ['85480']
    assert tslp['participants'][0]['taxon'] == '10090'
    assert rows[4]['detection_method_mi'] == 'MI:0492'
    # Binary and complex datasets split the observations; none is in both.
    assert len(list(parser.iter_netpath_interactions(ARCHIVE))) == 5
    assert [r['interaction_name'] for r in parser.iter_netpath_complexes(ARCHIVE)] == [
        'Interaction between TSLP and CRLF2, IL7R'
    ]


def test_archive_without_pathway_files_is_an_error():
    with pytest.raises(ValueError, match='no NetPath'):
        list(parser.iter_netpath_observations(SimpleNamespace(result={'readme.txt': b'x'})))
    with pytest.raises(ValueError, match='zip archive'):
        list(parser.iter_netpath_observations(SimpleNamespace(result='<html>503</html>')))


def test_binary_interaction_keeps_source_identifiers_method_and_pathway():
    rows = [r for r in parser.iter_netpath_interactions(ARCHIVE) if r['pathway_id'] == 'NetPath_24']
    out = extract('interactions', rows)
    # Two experiments are two observations of the same statement.
    assert [r.predicate for r in out.relations] == ['interacts_with'] * 2
    assert len({r.relation_key for r in out.relations}) == 1
    first, second = out.relations
    assert annotations(first, 'has_evidence_of_type') == ['MI:0493']
    assert annotations(second, 'has_evidence_of_type') == ['MI:0492']
    assert annotations(first, 'publications') == ['PMID:10974032', 'PMID:10881176']
    assert annotations(first, 'has_topic') == ['NetPath_24']
    assert annotations(first, 'in_taxon') == ['NCBITaxon:9606']
    assert annotations(first, 'psi_mi:experimental_role', 'subject') == ['MI:0497']
    entities = {e.identifier: e for e in out.entities.values()}
    assert set(entities) == {'85480', '64109'}
    tslp = entities['85480']
    assert (tslp.entity_type, tslp.namespace) == ('protein', Namespace.ENTREZ)
    assert {(i['ns'], i['id']) for i in tslp.identifiers} == {
        ('entrez', '85480'), ('refseq_protein', 'NP_149024.1'), ('hprd', '06108'),
        ('omim', '607003'), ('genesymbol', 'TSLP'),
        ('name', 'Thymic stromal lymphopoietin'),
        ('synonym', 'Thymic stromal lymphopoietin isoform 1'),
    }


def test_mouse_label_on_a_human_record_asserts_no_taxon():
    rows = list(parser.iter_netpath_interactions(ARCHIVE))
    relation = mod.map_netpath_interaction(next(r for r in rows if r['pathway_id'] == 'NetPath_24'))
    # NetPath gives "Mouse" with the identifiers of the human HPRD entry.
    assert not [a for a in relation.subject.annotations if a.term is mod.slots.in_taxon]
    out = extract('interactions', rows)
    taxa = {e.identifier: e.taxon for e in out.entities.values()}
    assert taxa['85480'] is None and taxa['64109'] is None
    assert taxa['8650'] == '9606' and taxa['3952'] == '9606'


def test_roles_missing_values_and_mislabelled_accessions():
    (row,) = [r for r in parser.iter_netpath_interactions(ARCHIVE) if r['pathway_id'] == 'NetPath_3']
    out = extract('interactions', [row])
    (relation,) = out.relations
    assert annotations(relation, 'psi_mi:experimental_role', 'subject') == ['MI:0496']
    assert annotations(relation, 'psi_mi:experimental_role', 'object') == ['MI:0498']
    itch = next(e for e in out.entities.values() if e.identifier == '83737')
    namespaces = {i['ns'] for i in itch.identifiers}
    # A UniProt accession under the "GenBank" label is not a GenBank accession:
    # the parser reads it by its form and the module emits it as UniProt.
    assert 'genbank' not in namespaces
    assert ('uniprot', 'Q96J02') in {(i['ns'], i['id']) for i in itch.identifiers}
    assert row['participants'][1]['uniprot'] == ['Q96J02']
    assert row['participants'][1]['genbank'] == []
    assert mod._accessions({'genbank': ['AAC24716.1', 'Q96J02', '']}, 'genbank') == ['AAC24716.1']
    assert mod._accessions({'uniprot': ['Q96J02', 'AAC24716.1']}, 'uniprot') == ['Q96J02']
    # NetPath writes a missing OMIM number as "0"; " 841" is an Entrez Gene ID.
    assert mod._accessions({'omim': ['0']}, 'omim') == []
    assert mod._accessions({'entrez': [' 841']}, 'entrez') == ['841']


def test_homo_oligomer_is_a_self_interaction():
    (row,) = [r for r in parser.iter_netpath_interactions(ARCHIVE) if r['pathway_id'] == 'NetPath_21']
    assert row['participant_count'] == 1
    out = extract('interactions', [row])
    (relation,) = out.relations
    assert relation.predicate == 'interacts_with'
    assert relation.subject_entity_key == relation.object_entity_key
    assert annotations(relation, 'in_taxon') == ['NCBITaxon:10090']


def test_participant_without_cross_reference_is_not_mapped():
    (row,) = [r for r in parser.iter_netpath_interactions(ARCHIVE) if r['pathway_id'] == 'NetPath_22']
    assert row['participant_count'] == 3
    assert [p['symbol'] for p in row['participants']] == ['LEP', 'LEPR']
    assert [(p['symbol'], p['type']) for p in row['unidentified_participants']] == [
        ('LEP', 'small molecule')
    ]
    out = extract('interactions', [row])
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.identifier for e in out.entities.values()} == {'3952', '3953'}
    assert annotations(out.relations[0], 'in_taxon') == ['NCBITaxon:40674']


def test_complex_interaction_is_a_complex_with_members_not_pairs():
    (row,) = parser.iter_netpath_complexes(ARCHIVE)
    assert mod.map_netpath_interaction(row) is None
    out = extract('complexes', [row])
    assert [r.predicate for r in out.relations] == ['has_member'] * 3
    assert not [r for r in out.relations if r.predicate == 'interacts_with']
    (complex_,) = [e for e in out.entities.values() if e.entity_type == 'macromolecular_complex']
    # No accession is invented for the complex.
    assert complex_.namespace == 'row' and not complex_.identifiers
    assert {(a['term'], a['value']) for a in complex_.annotations} == {
        ('has_evidence_of_type', 'MI:0493'), ('publications', 'PMID:10881176'),
        ('has_topic', 'NetPath_24'), ('in_taxon', 'NCBITaxon:9606'),
    }
    members = {e.identifier for e in out.entities.values() if e.entity_type == 'protein'}
    assert members == {'85480', '64109', '3575'}
    # The organism of the experiment, as for the interactions.
    assert complex_.taxon == '9606'
    assert mod.map_netpath_complex(next(parser.iter_netpath_interactions(ARCHIVE))) is None


def test_unknown_interaction_and_interactor_types_are_not_given_a_default():
    row = next(parser.iter_netpath_interactions(ARCHIVE))
    with pytest.raises(ValueError, match='Unknown NetPath interaction type'):
        mod.map_netpath_interaction(dict(row, interaction_type_mi='MI:0407'))
    chemical = dict(row['participants'][0], type='small molecule', type_mi='MI:0328')
    with pytest.raises(ValueError, match='Unknown NetPath interactor type'):
        mod.map_netpath_interaction(dict(row, participants=[chemical, row['participants'][1]]))
    catalysis = _file(_entry(4, 'x', [(['1'], HUMAN, VIVO)], [NUMB, ITCH],
                             interaction_type=('phosphorylation', 'MI:0217')))
    (row,) = parser.iter_netpath_interactions(SimpleNamespace(result={'NetPath_4_PSIMI.xml': catalysis}))
    with pytest.raises(ValueError, match="'MI:0217'"):
        mod.map_netpath_interaction(row)


def test_pathway_names_come_from_the_browse_page():
    page = SimpleNamespace(result=(
        '<html><body><a href="index.html">Home</a><table><tr>'
        '<td><a href="pathways?path_id=NetPath_24">Thymic stromal lymphopoietin (TSLP)</a></td>'
        '<td><a href="pathways?path_id=NetPath_7"> Transforming growth factor beta\n(TGF-beta) receptor </a></td>'
        '</tr></table></body></html>'
    ))
    rows = list(parser.iter_netpath_pathways(page))
    assert rows == [
        {'pathway_id': 'NetPath_24', 'name': 'Thymic stromal lymphopoietin (TSLP)'},
        {'pathway_id': 'NetPath_7', 'name': 'Transforming growth factor beta (TGF-beta) receptor'},
    ]
    out = extract('pathways', rows)
    pathway = next(e for e in out.entities.values() if e.identifier == 'NetPath_24')
    assert (pathway.entity_type, pathway.namespace) == ('pathway', 'netpath')
    assert pathway.label == 'Thymic stromal lymphopoietin (TSLP)'
    assert not out.relations
    with pytest.raises(ValueError, match='lists no pathway'):
        list(parser.iter_netpath_pathways(SimpleNamespace(result='<html>503</html>')))
