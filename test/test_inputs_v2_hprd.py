"""HPRD flat files: readers, mappers and extractor on rows copied from the archive."""

import io

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import hprd
from pypath.inputs_v2.parsers import hprd as parser

D = parser.ARCHIVE_DIRECTORY

# Rows of POST_TRANSLATIONAL_MODIFICATIONS.txt, with the line they have in the file.
PTM_ROWS = {
    # a kinase, with "None" in the experiment column
    60: '00020\tACTN1\t00020_1\tNP_001093.1\t12\tY\tPTK2\t02859\tPhosphorylation\tNone\t11369769,16291744',
    # the phosphatase of the same site
    61: '00020\tACTN1\t00020_1\tNP_001093.1\t12\tY\tPTPN1\t01477\tDephosphorylation\tin vitro;in vivo\t16291744',
    # a substrate whose sequence is given as a UniProt accession
    8064: '00925\tITPR1\t00925_1\tQ14643\t421\tS\tCDK1\t00302\tPhosphorylation\tin vitro;in vivo\t14635192',
    # a GenBank protein accession
    5247: '00558\tEPB41\t00558_1\tAAD42222.1\t660\tY\tEGFR\t00579\tPhosphorylation\tin vitro;in vivo\t1647028',
    # a cleavage, in vitro only
    303: '00024\tRFC1\t00024_1\tNP_002904.3\t1116\tD\tCASP3\t02799\tProteolytic Cleavage\tin vitro\t9351817',
    # a ubiquitinated lysine
    1726: '00182\tADRBK1\t00182_1\tNP_001610.2\t19\tK\tMDM2\t01272\tUbiquitination\tin vitro;in vivo\t17006543',
    # a modification type that has no aspect
    18566: '01859\tTP53\t01859_1\tNP_000537.3\t370\tK\tFBXO11\t06985\tNeddylation\tin vitro;in vivo\t17098746',
    # an enzyme with only an HPRD identifier, here the substrate itself
    52334: '06663\t-\t06663_1\tNP_057626.2\t178\tT\t-\t06663\tPhosphorylation\tin vitro;in vivo\t15037601,18691976,19651622',
    # the site cell "19;-"
    3291: '00302\tCDK1\t00302_2\tNP_203698.1\t19;-\tY\tCDK1\t00302\tPhosphorylation\tin vivo\t15592455,19413330',
    # a spelling of the experiment type
    15350: '01499\tPRKCB\t01499_1\tNP_002729.2\t16\tS\tPRKCB\t01499\tPhosphorylation\tinvitro\t8327493',
    # an eraser of a methyl group, found by yeast 2-hybrid
    78618: '11822\tHIST2H3C\t11822_1\tNP_066403.2\t36\tK\tKDM2A\t05741\tDemethylation\tyeast 2-hybrid\t16362057',
    # a reference that is not a PubMed ID
    9217: '01043\tRUNX1\t01043_2\tNP_001001890.1\t249\tS\tHIPK2\t06039\tPhosphorylation\tin vitro;in vivo\t169175071',
    # sites without an enzyme: an in vivo phosphorylation and a disulfide bridge
    11: '00001\tALDH1A1\t00001_1\tNP_000680.2\t267\tT\t-\t-\tPhosphorylation\tin vivo\t20166139',
    394: '00032\tACTB\t00032_1\tNP_001092.1\t285;374\tC\t-\t-\tDisulfide Bridge\tin vitro;in vivo\t7876306,8605189',
}

# Rows of BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt.
PPI_ROWS = (
    # 0: a protein with itself
    'ALDH1A1\t00001\tNP_000680.2\tALDH1A1\t00001\tNP_000680.2\tin vivo;yeast 2-hybrid\t12081471,16189514',
    # 1
    'ITGA7\t02761\tNP_001138468.1\tCHRNA1\t00007\tNP_001034612.1\tin vivo\t10910772',
    # 2, 3: one pair in both orientations, with the same evidence
    'JUN\t01302\tNP_002219.1\tSNAPC5\t16183\tNP_006040.1\t in vivo;in vitro\t15867431',
    'SNAPC5\t16183\tNP_006040.1\tJUN\t01302\tNP_002219.1\t in vivo;in vitro\t15867431',
    # 4: codes of the experiment column the archive does not explain
    'SVIL\t04992\tNP_068506.2\tSH3PXD2A\t10228\tNP_055446.2\tVV\t19109420',
    # 5: an abbreviation
    'SVIL\t04992\tNP_068506.2\tKRT18\t01020\tNP_000215.1\tY2H\t20309963',
    # 6: two PubMed IDs written together
    'MAPK1\t01496\tNP_002736.3\tIQGAP1\t04541\tNP_003861.1\tin vitro;in vivo\t1497021909368021',
    # 7: no gene symbol
    'SLC2A4\t00688\tNP_001033.1\t-\t10126\tNP_006446.1\tin vivo\t16396496',
    # 8: a GenBank accession
    'PAK1\t03995\tAAC24716.1\tERBB2\t01281\tNP_004439.2\tin vivo\t9774445',
    # 9: a trailing space in the identifier
    'HOOK1\t07428\tNP_056972.1\tHOOK3\t09699 \tNP_115786.1\tin vivo\t18799622',
)

# Rows of PROTEIN_COMPLEXES.txt.
COMPLEX_LINES = (
    'COM_1\t00011\tASCL1\tNP_004307.2\tin vivo\t8900141,8948587',
    'COM_1\t00918\tTCF3\tNP_003191.1\tin vivo\t8900141,8948587',
    'COM_1\t02809\tMEF2C\tNP_002388.2\tin vivo\t8900141,8948587',
    # a complex with a member that is no protein of the database
    'COM_1410\tNone\t-\t-\tin vivo\t10074429,10196146',
    'COM_1410\t00032\tACTB\tNP_001092.1\tin vivo\t10074429,10196146',
    'COM_1410\t03347\tCAPZA1\tNP_006126.1\tin vivo\t10074429,10196146',
    # a complex whose experiment column is empty and a member written twice
    'COM_1098\t10936\tELP3\tNP_060561.3\t\t11714725',
    'COM_1098\t15443\tELP2\tNP_060725.1\t\t11714725',
    'COM_1098\t15443\tELP2\tNP_060725.1\t\t11714725',
    # a second complex of the same three members as COM_1, from other work
    'COM_2567\t04042\tMNAT1\tNP_002422.1\tin vitro;in vivo\t8521818',
    'COM_2567\t15993\tCDK7\tNP_001790.1\tin vitro;in vivo\t8521818',
)


def archive(member, lines):
    text = '\n'.join(lines) + '\n'
    return type('Opener', (), {'result': {D + member: io.BytesIO(text.encode())}})()


def ptm(*lines, member='POST_TRANSLATIONAL_MODIFICATIONS.txt'):
    return archive(member, lines)


def ptm_rows(*keys):
    return list(parser.iter_enzyme_sites(ptm(*(PTM_ROWS[k] for k in keys))))


def all_rows(*keys):
    return list(parser.iter_modification_sites(ptm(*(PTM_ROWS[k] for k in keys))))


def interactions(lines=PPI_ROWS):
    return list(
        parser.iter_interactions(archive('BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt', lines))
    )


def complexes(lines=COMPLEX_LINES):
    return list(parser.iter_complexes(archive('PROTEIN_COMPLEXES.txt', lines)))


def extract(dataset, rows):
    out = SilverExtractor('hprd', dataset)
    mapper = getattr(hprd.resource, dataset).mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term, scope='relation'):
    return [
        a['value'] for a in relation.annotations if a['term'] == term and a['scope'] == scope
    ]


def test_importing_the_module_downloads_nothing_and_has_three_datasets():
    assert hprd.config.license.name == 'NONCOMMERCIAL'
    for name in ('enzyme_substrate', 'interactions', 'complexes'):
        assert getattr(hprd.resource, name).download.url.startswith(
            'https://rescued.omnipathdb.org/'
        )


# --- reading the modification table


def test_reader_gives_all_rows_and_the_enzyme_ones_only_when_asked():
    rows = all_rows(60, 11, 394, 52334)
    assert [r['has_enzyme'] for r in rows] == [True, False, False, True]
    assert [r['line'] for r in ptm_rows(60, 11, 394, 52334)] == [1, 4]
    first = rows[0]
    assert (first['substrate_hprd_id'], first['substrate_isoform_id']) == ('00020', '00020_1')
    assert (first['substrate_accession'], first['substrate_accession_kind']) == (
        'NP_001093.1',
        'refseq_protein',
    )
    assert (first['position'], first['second_position'], first['residue']) == (12, None, 'Y')
    assert (first['enzyme_name'], first['enzyme_hprd_id']) == ('PTK2', '02859')
    # The source's "None" states no experiment type; the cell stays as written.
    assert first['experiment_type'] == 'None'
    assert first['experiment_types'] == []
    assert first['references'] == ['11369769', '16291744']


def test_reader_decodes_site_cells_and_accession_kinds_and_keeps_in_vitro_rows():
    rows = {r['line']: r for r in all_rows(3291, 394, 8064, 5247, 52334, 303)}
    assert (rows[1]['position'], rows[1]['second_position']) == (19, None)
    assert rows[1]['site'] == '19;-'
    assert (rows[2]['position'], rows[2]['second_position']) == (285, 374)
    assert rows[3]['substrate_accession_kind'] == 'uniprot'
    assert rows[4]['substrate_accession_kind'] == 'genbank'
    # The row of an enzyme known only by its HPRD identifier.
    assert (rows[5]['enzyme_name'], rows[5]['enzyme_hprd_id']) == (None, '06663')
    assert rows[5]['substrate_gene_symbol'] is None
    # No row is dropped for its experiment type.
    assert [r['experiment_types'] for r in ptm_rows(303, 60)] == [['in vitro'], []]


def test_experiment_types_are_stripped_and_without_repeats():
    assert parser.experiment_types(' in vivo;in vitro;in vitro;in vivo') == ['in vivo', 'in vitro']
    assert parser.experiment_types('None;in vivo') == ['in vivo']
    assert parser.experiment_types('') == []
    assert parser.experiment_types('VT; VV;Y2H') == ['VT', 'VV', 'Y2H']


def test_reader_rejects_a_changed_or_wrong_file():
    with pytest.raises(ValueError, match='no member'):
        list(parser.iter_enzyme_sites(archive('OTHER.txt', [PTM_ROWS[60]])))
    with pytest.raises(ValueError, match='no rows'):
        list(parser.iter_enzyme_sites(ptm()))
    with pytest.raises(ValueError, match='fields'):
        list(parser.iter_enzyme_sites(ptm(PTM_ROWS[60] + '\textra')))
    with pytest.raises(ValueError, match='site not understood'):
        list(parser.iter_enzyme_sites(ptm(PTM_ROWS[60].replace('\t12\t', '\tY12\t'))))
    with pytest.raises(ValueError, match='accession not understood'):
        list(parser.iter_enzyme_sites(ptm(PTM_ROWS[60].replace('NP_001093.1', 'X1'))))
    with pytest.raises(ValueError, match='identifier not understood'):
        list(parser.iter_enzyme_sites(ptm(PTM_ROWS[60].replace('02859', '2859'))))
    with pytest.raises(ValueError, match='fields'):
        list(parser.iter_enzyme_sites(ptm('<html><body>503</body></html>')))


# --- enzyme-substrate relations


def test_kinase_increases_substrate_phosphorylation_at_the_site():
    (row,) = ptm_rows(8064)
    record = hprd.resource.enzyme_substrate.mapper(row)
    assert record.subject.identifiers == [('hprd', '00302'), (Namespace.GENESYMBOL, 'CDK1')]
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'Q14643'),
        ('hprd', '00925'),
        ('hprd_isoform', '00925_1'),
        (Namespace.GENESYMBOL, 'ITPR1'),
    ]
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'S',
            'position': 421,
            'end_position': 421,
            'coordinate_reference': {
                'identifier': {'ns': 'uniprot', 'id': 'Q14643'},
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': None,
        }
    ]
    out = extract('enzyme_substrate', [row])
    (relation,) = out.relations
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'has_evidence_of_type') == ['ECO:0000181', 'ECO:0000178']
    assert annotations(relation, 'publications') == ['PMID:14635192']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}


def test_refseq_isoform_is_the_sequence_the_position_refers_to():
    (row,) = ptm_rows(60)
    substrate = hprd.resource.enzyme_substrate.mapper(row).object
    assert substrate.identifiers[0] == (Namespace.REFSEQ_PROTEIN, 'NP_001093.1')
    assert substrate.molecular_form['sequence_identifiers'] == [
        {'ns': 'refseq_protein', 'id': 'NP_001093.1'}
    ]
    assert substrate.molecular_form['modifications'][0]['coordinate_reference'][
        'identifier'
    ] == {'ns': 'refseq_protein', 'id': 'NP_001093.1'}
    # GenBank and UniProt accessions are the coordinate reference only.
    (genbank,) = ptm_rows(5247)
    form = hprd.resource.enzyme_substrate.mapper(genbank).object.molecular_form
    assert form['sequence_identifiers'] is None
    assert form['modifications'][0]['coordinate_reference']['identifier'] == {
        'ns': 'genbank',
        'id': 'AAD42222.1',
    }
    # Two isoforms of one protein are two substrates.
    other = dict(row, substrate_isoform_id='00020_2', substrate_accession='NP_001094.1')
    out = extract('enzyme_substrate', [row, other])
    assert len({r.object_entity_key for r in out.relations}) == 2


def test_phosphatase_decreases_phosphorylation_and_no_activity_effect_is_given():
    kinase, phosphatase = (
        extract('enzyme_substrate', ptm_rows(k)).relations[0] for k in (60, 61)
    )
    assert annotations(kinase, 'object_direction_qualifier') == ['increased']
    assert annotations(phosphatase, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(phosphatase, 'object_direction_qualifier') == ['decreased']
    out = extract('enzyme_substrate', ptm_rows(60, 61, 303, 1726, 78618))
    terms = {a['term'] for r in out.relations for a in r.annotations}
    assert not terms & {'activity', 'activity_or_abundance', 'subject_aspect_qualifier'}
    for relation in out.relations:
        assert 'activity' not in annotations(relation, 'object_aspect_qualifier')


def test_other_modification_types_use_their_own_aspect():
    out = extract('enzyme_substrate', ptm_rows(303, 1726, 78618))
    assert [
        (annotations(r, 'object_aspect_qualifier'), annotations(r, 'object_direction_qualifier'))
        for r in out.relations
    ] == [
        (['cleavage'], ['increased']),
        (['ubiquitination'], ['increased']),
        (['methylation'], ['decreased']),
    ]
    cleavage, ubiquitin, demethylation = (
        hprd.resource.enzyme_substrate.mapper(r).object.molecular_form['modifications'][0]
        for r in ptm_rows(303, 1726, 78618)
    )
    # No PSI-MOD name for a cleavage site: HPRD's type stays as description.
    assert (cleavage['term'], cleavage['description']) == (None, 'Proteolytic Cleavage')
    assert (cleavage['residue'], cleavage['position']) == ('D', 1116)
    assert ubiquitin['term'] == 'ubiquitinylated lysine'
    assert demethylation['term'] == 'methylated residue'


def test_experiment_types_are_evidence_and_rows_are_not_filtered_by_them():
    out = extract('enzyme_substrate', ptm_rows(303, 8064, 60, 15350, 78618))
    assert [annotations(r, 'has_evidence_of_type') for r in out.relations] == [
        ['ECO:0000181'],
        ['ECO:0000181', 'ECO:0000178'],
        [],
        ['ECO:0000181'],
        ['ECO:0000068'],
    ]


def test_sites_without_an_enzyme_are_not_relations():
    assert ptm_rows(11, 394) == []
    out = extract('enzyme_substrate', all_rows(11, 394, 60))
    assert len(out.relations) == 1
    assert all_rows(11, 394)[1]['second_position'] == 374


def test_type_without_aspect_gives_no_relation_and_unknown_values_stop():
    (neddylation,) = ptm_rows(18566)
    assert hprd.resource.enzyme_substrate.mapper(neddylation) is None
    assert extract('enzyme_substrate', [neddylation]).relations == []
    (row,) = ptm_rows(60)
    with pytest.raises(ValueError, match='modification type'):
        hprd.resource.enzyme_substrate.mapper(dict(row, modification_type='Phosphorylated'))
    with pytest.raises(ValueError, match='experiment type'):
        hprd.resource.enzyme_substrate.mapper(dict(row, experiment_types=['in cellulo']))
    with pytest.raises(ValueError, match='two positions'):
        hprd.resource.enzyme_substrate.mapper(dict(row, second_position=5))


def test_enzyme_known_only_by_hprd_id_and_site_with_empty_second_position():
    (row,) = ptm_rows(52334)
    record = hprd.resource.enzyme_substrate.mapper(row)
    assert record.subject.identifiers == [('hprd', '06663')]
    assert record.object.identifiers == [
        (Namespace.REFSEQ_PROTEIN, 'NP_057626.2'),
        ('hprd', '06663'),
        ('hprd_isoform', '06663_1'),
    ]
    (row,) = ptm_rows(3291)
    assert hprd.resource.enzyme_substrate.mapper(row).object.molecular_form['modifications'][0][
        'position'
    ] == 19


def test_only_well_formed_pubmed_ids_become_publications():
    (row,) = ptm_rows(9217)
    out = extract('enzyme_substrate', [row])
    assert annotations(out.relations[0], 'publications') == []


# --- interactions


def test_interaction_is_a_physical_interaction_between_identified_proteins():
    rows = interactions()
    assert len(rows) == 10
    record = hprd.resource.interactions.mapper(rows[1])
    assert record.subject.identifiers == [
        (Namespace.REFSEQ_PROTEIN, 'NP_001138468.1'),
        ('hprd', '02761'),
        (Namespace.GENESYMBOL, 'ITGA7'),
    ]
    assert record.object.identifiers[1:] == [('hprd', '00007'), (Namespace.GENESYMBOL, 'CHRNA1')]
    out = extract('interactions', rows[1:2])
    (relation,) = out.relations
    assert relation.predicate == 'interacts_with'
    assert annotations(relation, 'has_evidence_of_type') == ['ECO:0000178']
    assert annotations(relation, 'publications') == ['PMID:10910772']
    assert {e.taxon for e in out.entities.values()} == {'9606'}
    # The source does not say the interaction is direct.
    assert relation.predicate != 'directly_physically_interacts_with'


def test_experiment_types_of_interactions():
    out = extract('interactions', interactions())
    evidence = [annotations(r, 'has_evidence_of_type') for r in out.relations]
    assert evidence[0] == ['ECO:0000178', 'ECO:0000068']
    assert evidence[2] == ['ECO:0000178', 'ECO:0000181']
    # "VV" is not given a meaning; "Y2H" is the abbreviation of yeast 2-hybrid.
    assert evidence[4] == []
    assert evidence[5] == ['ECO:0000068']
    rows = interactions()
    with pytest.raises(ValueError, match='experiment type'):
        hprd.resource.interactions.mapper(dict(rows[1], experiment_types=['in situ']))


def test_self_interaction_pair_in_both_orientations_and_odd_cells():
    rows = interactions()
    out = extract('interactions', rows)
    assert len(out.relations) == 10
    assert out.relations[0].subject_entity_key == out.relations[0].object_entity_key
    # Both orientations are rows of the source and both are emitted.
    assert rows[2]['interactor_1_hprd_id'] == rows[3]['interactor_2_hprd_id']
    assert rows[2]['interactor_2_hprd_id'] == rows[3]['interactor_1_hprd_id']
    # No gene symbol: the protein is identified by its other identifiers.
    assert rows[7]['interactor_2_gene_symbol'] is None
    assert out.relations[7].object_entity_key
    # The identifier cell is stripped.
    assert rows[9]['interactor_2_hprd_id'] == '09699'
    assert rows[8]['interactor_1_accession_kind'] == 'genbank'
    # Two PubMed IDs written together are not split or repaired.
    assert annotations(out.relations[6], 'publications') == []
    assert rows[6]['references'] == ['1497021909368021']


def test_interaction_reader_rejects_other_files():
    with pytest.raises(ValueError, match='no member'):
        list(parser.iter_interactions(ptm(PPI_ROWS[1])))
    with pytest.raises(ValueError, match='fields'):
        interactions([PPI_ROWS[1] + '\tx'])
    with pytest.raises(ValueError, match='identifier not understood'):
        interactions([PPI_ROWS[1].replace('02761', 'ITGA7')])


# --- complexes


def test_complex_has_its_members_and_is_not_expanded_to_pairs():
    rows = complexes()
    assert [r['complex_id'] for r in rows] == ['COM_1', 'COM_1410', 'COM_1098', 'COM_2567']
    out = extract('complexes', rows[:1])
    assert [r.predicate for r in out.relations] == ['has_member'] * 3
    assert not [r for r in out.relations if r.predicate == 'interacts_with']
    (cx,) = [e for e in out.entities.values() if e.entity_type == 'macromolecular_complex']
    assert cx.identifier == 'COM_1'
    members = [e for e in out.entities.values() if e.entity_type == 'protein']
    assert {i['id'] for e in members for i in e.identifiers if i['ns'] == 'hprd'} == {
        '00011',
        '00918',
        '02809',
    }
    assert {e.taxon for e in members} == {'9606'}
    assert {(a['term'], a['value']) for a in cx.annotations} == {
        ('has_evidence_of_type', 'ECO:0000178'),
        ('publications', 'PMID:8900141'),
        ('publications', 'PMID:8948587'),
    }
    # The source gives the complex no name and no stoichiometry.
    assert [i['ns'] for i in cx.identifiers] == ['hprd_complex']


def test_complex_with_an_unidentified_member_is_not_mapped():
    rows = complexes()
    assert rows[1]['unidentified_rows'] == 1
    assert [m['hprd_id'] for m in rows[1]['members']] == ['00032', '03347']
    assert hprd.resource.complexes.mapper(rows[1]) is None
    assert extract('complexes', [rows[1]]).relations == []


def test_repeated_member_rows_are_one_member_and_empty_experiment_is_none():
    row = complexes()[2]
    assert [m['hprd_id'] for m in row['members']] == ['10936', '15443']
    assert row['duplicate_rows'] == 1
    assert row['experiment_types'] == []
    out = extract('complexes', [row])
    assert len(out.relations) == 2
    assert [a for a in out.entities.values() if a.entity_type == 'macromolecular_complex']


def test_complex_with_fewer_than_two_members_is_not_mapped():
    row = complexes()[3]
    assert hprd.resource.complexes.mapper(row) is not None
    assert hprd.resource.complexes.mapper(dict(row, members=row['members'][:1])) is None


def test_complex_reader_rejects_a_split_or_inconsistent_complex():
    split = COMPLEX_LINES[:1] + COMPLEX_LINES[9:10] + COMPLEX_LINES[1:2]
    with pytest.raises(ValueError, match='not in one block'):
        complexes(split)
    other = COMPLEX_LINES[1].replace('in vivo', 'in vitro')
    with pytest.raises(ValueError, match='differs'):
        complexes([COMPLEX_LINES[0], other])
    twice = COMPLEX_LINES[0].replace('ASCL1', 'XXXX')
    with pytest.raises(ValueError, match='twice'):
        complexes([COMPLEX_LINES[0], twice])
    with pytest.raises(ValueError, match='no member'):
        list(parser.iter_complexes(ptm(COMPLEX_LINES[0])))
