"""MatrixDB: parsers, mappers and extractor on source rows."""
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import entity_type, slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import matrixdb
from pypath.inputs_v2.parsers import matrixdb as parser

HEADER = ['#ID(s) interactor A', 'ID(s) interactor B', 'Alt. ID(s) interactor A', 'Alt. ID(s) interactor B', 'Alias(es) interactor A', 'Alias(es) interactor B', 'Interaction detection method(s)', 'Publication 1st author(s)', 'Publication Identifier(s)', 'Taxid interactor A', 'Taxid interactor B', 'Interaction type(s)', 'Source database(s)', 'Interaction identifier(s)', 'Confidence value(s)', 'Expansion method(s)', 'Biological role(s) interactor A', 'Biological role(s) interactor B', 'Experimental role(s) interactor A', 'Experimental role(s) interactor B', 'Type(s) interactor A', 'Type(s) interactor B', 'Xref(s) interactor A', 'Xref(s) interactor B', 'Interaction Xref(s)', 'Annotation(s) interactor A', 'Annotation(s) interactor B', 'Interaction annotation(s)', 'Host organism(s)', 'Interaction parameter(s)', 'Creation date', 'Update date', 'Checksum(s) interactor A', 'Checksum(s) interactor B', 'Interaction Checksum(s)', 'Negative', 'Feature(s) interactor A', 'Feature(s) interactor B', 'Stoichiometry(s) interactor A', 'Stoichiometry(s) interactor B', 'Identification method participant A', 'Identification method participant B']
# Rows copied from matrixdb_CORE.tab, cell by cell and with the quoting of the
# file: some cells are wrapped in quotes with their own quotes doubled, others only
# have the quotes doubled. The file pads every line with empty cells.
FIRST = [
    '"chebi:""CHEBI:28304"""',
    'uniprotkb:P15127',
    'matrixdb:GAG_1',
    '-',
    'matrixdb:Heparin(short label)',
    '"uniprotkb:""Insr""(gene name)"',
    '"psi-mi:""MI:0400""(""affinity technology"")"',
    'Ori A et al.(2011)',
    'pubmed:21454685|imex:IM-25784',
    '-',
    'taxid:10116(Rattus norvegicus)',
    '"psi-mi:""MI:0914""(association)"',
    '"psi-mi:""MI:0917""(matrixdb)"',
    'matrixdb:A0JPN2__B2GUY2__...__Q9Z2H4_21454685_matrixdb_1|intact:EBI-15184828|imex:IM-25784-1',
    '-',
    '"psi-mi:""MI:1060""(spoke expansion)"',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0496""(bait)"',
    '"psi-mi:""MI:0498""(prey)"',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '"figure legend:""Supplemental Table 1""|curation depth:imex curation"',
    'taxid:-1(in vitro)',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    'false',
    '-',
    '-',
    '-',
    '-',
    '"psi-mi:""MI:0102""(sequence tag identification)"',
    '"psi-mi:""MI:0102""(sequence tag identification)"',
]
TRUNCATED = [
    '"chebi:""CHEBI:28304"""',
    'uniprotkb:P04157',
    'matrixdb:GAG_1',
    '-',
    'matrixdb:Heparin(short label)',
    '"uniprotkb:""Ptprc""(gene name)"',
    '"psi-mi:""MI:0400""(""affinity technology"")"',
    'Ori A et al.(2011)',
    'pubmed:21454685|imex:IM-25784',
    '-',
    'taxid:10116(Rattus norvegicus)',
    '"psi-mi:""MI:0914""(association)"',
    '"psi-mi:""MI:0917""(matrixdb)"',
    'matrixdb:A0JPN2__B2GUY2__...__Q9Z2H4_21454685_matrixdb_1|intact:EBI-15184828|imex:IM-25784-1',
    '-',
    '"psi-mi:""MI:1060""(spoke expansion)"',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0496""(bait)"',
    '"psi-mi:""MI:0498""(prey)"',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '"figure legend:""Supplemental Table 1""|curation depth:imex curation"',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
    '',
]
PFRAG = [
    '"chebi:""CHEBI:29105"""',
    'uniprotkb:P39060-PRO_0000005794',
    'matrixdb:CAT_1',
    'matrixdb:PFRAG_1_human',
    'matrixdb:Zinc(short label)',
    'matrixdb:Endostatin(short label)',
    '"psi-mi:""MI:0013""(""biophysical"")"',
    'Ding YH et al.(1998)',
    'pubmed:9724722',
    '-',
    'taxid:9606(Homo sapiens)',
    '"psi-mi:""MI:0407""(direct interaction)"',
    '"psi-mi:""MI:0917""(matrixdb)"',
    'matrixdb:CAT_1__PFRAG_1_human_9724722_MatrixDB_1',
    '-',
    '-',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0497""(neutral component)"',
    '"psi-mi:""MI:0497""(neutral component)"',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '"comment:""Binding site: His1, His3, His11, Asp76""|curation depth:mimix curation"',
    'taxid:-1(in vitro)',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    'false',
    '-',
    '-',
    '-',
    '-',
    '"psi-mi:""MI:0396""(predetermined participant)"',
    '"psi-mi:""MI:0078""(nucleotide sequence identification)"',
]
COLOC = [
    '"chebi:""CHEBI:28815"""',
    'uniprotkb:P15502',
    'matrixdb:GAG_2',
    '-',
    'matrixdb:Heparan sulfate(short label)',
    '"uniprotkb:""ELN""(gene name)"',
    '"psi-mi:""MI:0040""(""electron microscopy"")"',
    'Gheduzzi D et al.(2005)',
    'pubmed:15748998',
    '-',
    'taxid:9606(Homo sapiens)',
    '"psi-mi:""MI:0403""(colocalization)"',
    '"psi-mi:""MI:0917""(matrixdb)"',
    'matrixdb:GAG_2__P15502_15748998_MatrixDB_1',
    '-',
    '-',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0497""(neutral component)"',
    '"psi-mi:""MI:0497""(neutral component)"',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '"figure legend:""1""|curation depth:mimix curation"',
    'taxid:-1(in vitro)',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    'false',
    '-',
    '-',
    '-',
    '-',
    '"psi-mi:""MI:0422""(immunostaining)"',
    '"psi-mi:""MI:0396""(predetermined participant)"',
]
RANGE = [
    '"chebi:""CHEBI:29108"""',
    'uniprotkb:P07214',
    'matrixdb:CAT_3',
    '-',
    'matrixdb:Calcium(short label)',
    '"uniprotkb:""SPARC""(gene name)"',
    '"psi-mi:""MI:0016""(""circular dichroism"")"',
    'Maurer P et al.(1992)',
    'pubmed:1555584|imex:IM-14109',
    '-',
    'taxid:10090(Mus musculus)',
    '"psi-mi:""MI:0407""(direct interaction)"',
    '"psi-mi:""MI:0917""(matrixdb)"',
    'matrixdb:CAT_3__P07214_1555584_MatrixDB_1|imex:IM-14109-1',
    '-',
    '-',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0497""(neutral component)"',
    '"psi-mi:""MI:0497""(neutral component)"',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '"figure legend:""2""|curation depth:imex curation"',
    'taxid:-1(in vitro)',
    'kd:9e-07~3e-07(molar)',
    '-',
    '-',
    '-',
    '-',
    '-',
    'false',
    '-',
    '-',
    '-',
    '1',
    '"psi-mi:""MI:0396""(predetermined participant)"',
    '"psi-mi:""MI:0396""(predetermined participant)"',
]
KON = [
    'chebi:""CHEBI:28304"',
    'complex portal:CPX-1794',
    'matrixdb:GAG_1',
    'matrixdb:MULT_15_human',
    'matrixdb:Heparin(short label)',
    'matrixdb:Integrin-alpha5beta1(short label)',
    '"psi-mi:""MI:0107""(""surface plasmon resonance"")"',
    'Faye C et al.(2009)',
    'pubmed:19502598|imex:IM-16038',
    '-',
    'taxid:9606(Homo sapiens)',
    '"psi-mi:""MI:0407""(direct interaction)"',
    '"psi-mi:""MI:0917""(matrixdb)"',
    'matrixdb:GAG_1__MULT_15_human_19502598_MatrixDB_4|imex:IM-16038-7',
    '-',
    '-',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0499""(unspecified role)"',
    '"psi-mi:""MI:0496""(bait)"',
    '"psi-mi:""MI:0498""(prey)"',
    '-',
    '-',
    '-',
    '-',
    '-',
    '"author-name:""heparin"""',
    '"author-name:""full-length alpha5beta1 integrin ectodomain"""',
    '"figure legend:""6A, 6C""|comment:""experiment realised with the full-length alpha5beta1 integrin ectodomain""|curation depth:imex curation"',
    'taxid:-1(in vitro)',
    'kd:1.55e-08(molar)|kon:14700(per mole per second)|koff:2.28e-05(second -1)',
    '-',
    '-',
    '-',
    '-',
    '-',
    'false',
    '"biotin tag"":?-?',
    '""sufficient binding region"":1-951,1-708',
    '-',
    '-',
    'psi-mi:""MI:0396""(predetermined participant)',
    'psi-mi:""MI:0396""(predetermined participant)',
]
MULT = [
    'matrixdb:MULT_3_sheep',
    'complex portal:CPX-1876',
    '-',
    'matrixdb:MULT_78_human',
    'matrixdb:Collagen-I(short label)',
    'matrixdb:PDGF-BB(short label)',
    'psi-mi:""MI:0107""(""surface plasmon resonance"")',
    'Gohring W et al.(1998)',
    'pubmed:9692901|imex:IM-14103',
    'taxid:9940(Ovis aries)',
    'taxid:9606(Homo sapiens)',
    'psi-mi:""MI:0407""(direct interaction)',
    'psi-mi:""MI:0917""(matrixdb)',
    'matrixdb:MULT_3_sheep__MULT_78_human_9692901_MatrixDB_1|imex:IM-14103-15',
    '-',
    '-',
    'psi-mi:""MI:0499""(unspecified role)',
    'psi-mi:""MI:0499""(unspecified role)',
    'psi-mi:""MI:0498""(prey)',
    'psi-mi:""MI:0496""(bait)',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    'figure legend:""Table 2""|curation depth:imex curation',
    'taxid:-1(in vitro)',
    'kd:3.3e-08(molar)|kon:190000(per mole per second)|koff:0.0062(second -1)',
    '-',
    '-',
    '-',
    '-',
    '-',
    'false',
    '-',
    '-',
    '-',
    '-',
    'psi-mi:""MI:0396""(predetermined participant)',
    'psi-mi:""MI:0396""(predetermined participant)',
]
ENZYME = [
    'uniprotkb:P16636',
    'uniprotkb:O43897',
    '-',
    '-',
    'uniprotkb:""Lox""(gene name)',
    'uniprotkb:""Tll1""(gene name)',
    'psi-mi:""MI:0415""(""enzymatic study"")',
    'Uzel MI et al.(2001)',
    'pubmed:11313359',
    'taxid:10116(Rattus norvegicus)',
    'taxid:9606(Homo sapiens)',
    'psi-mi:""MI:0407""(direct interaction)',
    'psi-mi:""MI:0917""(matrixdb)',
    'matrixdb:O43897__P16636_11313359_MatrixDB_1',
    '-',
    '-',
    'psi-mi:""MI:0502""(enzyme target)',
    'psi-mi:""MI:0501""(enzyme)',
    'psi-mi:""MI:0497""(neutral component)',
    'psi-mi:""MI:0497""(neutral component)',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    'curation depth:mimix curation',
    'taxid:-1(in vitro)',
    '-',
    '-',
    '-',
    '-',
    '-',
    '-',
    'false',
    '-',
    '-',
    '-',
    '-',
    'psi-mi:""MI:0396""(predetermined participant)',
    'psi-mi:""MI:0816""(molecular weight estimation by staining)',
]

ROWS = [FIRST, TRUNCATED, PFRAG, COLOC, RANGE, KON, MULT, ENZYME]

ECM_CSV = (
    'protein_id,species (NCBI),go_terms,uniprot_keywords,dataset,interaction_count\n'
    'P39060,9606,Basement membrane|Extracellular region|Extracellular space,'
    'Basement membrane|Extracellular matrix|Extracellular space/secreted,Swiss-Prot,157\n'
    'P05067,9606,Extracellular region|Extracellular space,Extracellular space/secreted,Swiss-Prot,1115\n'
    'Q5ZKC9,9031,Extracellular space,,Swiss-Prot,1\n'
    'A8K2U0,9606,,Extracellular space/secreted,Swiss-Prot,37\n'
    'A0A024R1N1,9606,Extracellular space,,UniProtKB unreviewed (TrEMBL),2\n'
)


def opener(text):
    return SimpleNamespace(result=text)


def interaction_file(rows=ROWS, header=HEADER, padding=3):
    lines = ['\t'.join([*header, *[''] * padding])]
    for row in rows:
        lines.append('\t'.join([*row, *[''] * padding]))
    return '\n'.join(lines) + '\n'


def parsed(rows=ROWS):
    return list(parser.iter_interactions(opener(interaction_file(rows))))


def extract(dataset, rows):
    out = SilverExtractor('matrixdb', dataset)
    mapper = matrixdb.resource.datasets()[dataset].mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def attributes(relation, scope='relation'):
    result = {}
    for a in relation.annotations:
        if a['scope'] == scope:
            result.setdefault(a['term'], []).append(a['value'])
    return result


def quantities(relation):
    return [
        a['quantity']
        for a in relation.annotations
        if a['scope'] == 'relation' and a['term'] == 'has_quantitative_value'
    ]


def test_parser_undoes_the_two_kinds_of_quoting_of_the_file():
    rows = parsed()
    assert len(rows) == 8
    first = rows[0]
    # Wrapped in quotes with doubled quotes, and doubled quotes only.
    assert first['id_a'] == 'chebi:"CHEBI:28304"'
    assert first['detection_methods'] == 'psi-mi:"MI:0400"("affinity technology")'
    assert rows[6]['detection_methods'] == 'psi-mi:"MI:0107"("surface plasmon resonance")'
    assert rows[6]['interaction_types'] == 'psi-mi:"MI:0407"(direct interaction)'
    assert first['negative'] == 'false'


def test_the_row_that_stops_early_has_no_negative_flag_and_is_kept():
    rows = parsed()
    truncated = rows[1]
    assert truncated['negative'] == '-'
    assert truncated['host_organisms'] == '-'
    assert truncated['interaction_annotations'].startswith('figure legend')
    # Not a negative interaction: it is mapped.
    assert matrixdb.exclusion_reason(truncated) is None
    assert matrixdb.map_interaction(truncated) is not None


@pytest.mark.parametrize(
    'text',
    [
        '',
        interaction_file([]),
        interaction_file(header=[*HEADER[:41], 'Other']),
        interaction_file(header=HEADER[:41]),
        interaction_file([[*FIRST[:35], 'maybe', *FIRST[36:]]]),
    ],
)
def test_parser_rejects_files_that_are_not_the_interactions(text):
    with pytest.raises(ValueError):
        list(parser.iter_interactions(opener(text)))


def test_parser_rejects_content_after_the_42nd_cell():
    header, first, *rest = interaction_file(padding=0).split('\n')
    text = '\n'.join([header, first + '\tstray', *rest])
    with pytest.raises(ValueError, match='after the 42nd'):
        list(parser.iter_interactions(opener(text)))


def test_negative_interaction_is_excluded():
    row = {**parsed()[0], 'negative': 'true'}
    assert matrixdb.exclusion_reason(row) == 'negative result'
    assert matrixdb.map_interaction(row) is None


def test_glycosaminoglycan_and_protein_participants():
    row = parsed()[0]
    record = matrixdb.map_interaction(row)
    assert slot_name(record.predicate) == 'interacts_with'
    assert entity_type(record.subject.type) == 'chemical_entity'
    assert record.subject.identifiers == [
        (Namespace.CHEBI, '28304'),
        ('matrixdb', 'GAG_1'),
        (Namespace.NAME, 'Heparin'),
    ]
    assert entity_type(record.object.type) == 'protein'
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'P15127'),
        (Namespace.GENESYMBOL, 'Insr'),
    ]
    assert record.object.molecular_form is None
    assert record.identifiers == [
        ('matrixdb_interaction', 'A0JPN2__B2GUY2__...__Q9Z2H4_21454685_matrixdb_1'),
        (Namespace.INTACT, 'EBI-15184828'),
        (Namespace.IMEX, 'IM-25784-1'),
    ]
    out = extract('interactions', [row])
    attrs = attributes(out.relations[0])
    # The taxon of the interaction in vitro is no organism; the experiment's
    # organism is the participant's.
    assert attrs['has_evidence_of_type'] == ['MI:0400']
    assert attrs['original_predicate'] == ['MI:0914']
    assert attrs['publications'] == ['PMID:21454685']
    assert attrs['matrixdb:expansion_method'] == ['MI:1060']
    assert attributes(out.relations[0], 'subject')['psi_mi:experimental_role'] == ['MI:0496']
    assert attributes(out.relations[0], 'object')['psi_mi:experimental_role'] == ['MI:0498']
    # "unspecified role" is no role.
    assert 'psi_mi:biological_role' not in attributes(out.relations[0], 'subject')
    protein = next(e for e in out.entities.values() if e.entity_type == 'protein')
    assert protein.taxon == '10116'


def test_protein_fragment_keeps_its_chain_and_matrixdb_identifier():
    record = matrixdb.map_interaction(parsed()[2])
    assert slot_name(record.predicate) == 'directly_physically_interacts_with'
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'P39060-PRO_0000005794'),
        ('matrixdb', 'PFRAG_1_human'),
        (Namespace.NAME, 'Endostatin'),
    ]
    assert record.object.molecular_form['sequence_identifiers'] == [
        {'ns': 'uniprot', 'id': 'P39060-PRO_0000005794'}
    ]
    assert record.subject.identifiers[:2] == [
        (Namespace.CHEBI, '29105'),
        ('matrixdb', 'CAT_1'),
    ]


def test_complexes_and_multimers_are_complexes_without_members():
    record = matrixdb.map_interaction(parsed()[5])
    assert entity_type(record.object.type) == 'macromolecular_complex'
    assert record.object.identifiers == [
        (Namespace.COMPLEXPORTAL, 'CPX-1794'),
        ('matrixdb', 'MULT_15_human'),
        (Namespace.NAME, 'Integrin-alpha5beta1'),
    ]
    multimer = matrixdb.map_interaction(parsed()[6])
    assert entity_type(multimer.subject.type) == 'macromolecular_complex'
    assert multimer.subject.identifiers == [
        ('matrixdb', 'MULT_3_sheep'),
        (Namespace.NAME, 'Collagen-I'),
    ]
    out = extract('interactions', parsed())
    # A complex is never expanded into its members, and no membership is made.
    assert 'has_member' not in {r.predicate for r in out.relations}


def test_colocalization_is_not_a_physical_interaction():
    record = matrixdb.map_interaction(parsed()[3])
    assert slot_name(record.predicate) == 'colocalizes_with'


def test_kinetic_parameters_are_measurements_and_ranges_are_left_out():
    rows = parsed()
    out = extract('interactions', [rows[5], rows[4]])
    kinetics = [quantities(r) for r in out.relations]
    assert len(kinetics[0]) == 3
    assert (kinetics[0][0]['source_field'], kinetics[0][0]['has_unit']) == ('kd', 'molar')
    assert kinetics[0][0]['has_numeric_value'] == 1.55e-08
    assert kinetics[0][1]['source_field'] == 'kon'
    assert kinetics[0][2]['source_field'] == 'koff'
    # "kd:9e-07~3e-07(molar)" is an interval of two values: no measurement.
    assert kinetics[1] == []


def test_enzyme_roles_are_kept_on_the_participants():
    out = extract('interactions', [parsed()[7]])
    # Lox (P16636) is the enzyme target and Tll1 (O43897) the enzyme.
    assert attributes(out.relations[0], 'subject')['psi_mi:biological_role'] == ['MI:0502']
    assert attributes(out.relations[0], 'object')['psi_mi:biological_role'] == ['MI:0501']


def test_unknown_interactor_or_interaction_type_raises():
    row = parsed()[0]
    with pytest.raises(ValueError, match='Unknown MatrixDB biomolecule'):
        matrixdb.map_interaction({**row, 'id_a': 'matrixdb:XYZ_1'})
    with pytest.raises(ValueError, match='Unknown MatrixDB interaction type'):
        matrixdb.map_interaction({**row, 'interaction_types': 'psi-mi:"MI:0208"(genetic interaction)'})


def test_ecm_protein_parser_and_mapper():
    rows = list(parser.iter_ecm_proteins(opener(ECM_CSV)))
    assert rows[0] == {
        'accession': 'P39060',
        'taxon': '9606',
        'go_terms': ['Basement membrane', 'Extracellular region', 'Extracellular space'],
        'keywords': ['Basement membrane', 'Extracellular matrix', 'Extracellular space/secreted'],
        'dataset': 'Swiss-Prot',
        'interaction_count': 157,
    }
    assert rows[3]['go_terms'] == [] and rows[2]['keywords'] == []
    record = matrixdb.map_ecm_protein(rows[0])
    assert slot_name(record.predicate) == 'associated_with'
    assert record.subject.identifiers == [(Namespace.UNIPROT, 'P39060')]
    assert record.object.identifier == 'ecm_proteins'
    out = extract('ecm_proteins', rows)
    # The chicken protein and the unreviewed entry are not kept.
    assert len(rows) == 5 and len(out.relations) == 3
    attrs = attributes(out.relations[0])
    # The names of the criteria table become the accessions that table gives.
    assert attrs['matrixdb:go_cellular_component'] == ['GO:0005604', 'GO:0005576', 'GO:0005615']
    assert attrs['matrixdb:uniprot_keyword'] == ['KW-0084', 'KW-0272', 'KW-0964']
    # The count is a Measurement: no text value, a quantity of 157 from its column.
    assert attrs['matrixdb:interaction_count'] == ['']
    (count,) = [
        a['quantity'] for a in out.relations[0].annotations if a['term'] == 'matrixdb:interaction_count'
    ]
    assert count['has_numeric_value'] == 157.0 and count['source_field'] == 'interaction_count'
    # Every kept protein is a reviewed human entry, so the dataset is no attribute.
    assert 'matrixdb:uniprot_dataset' not in attrs
    assert {e.taxon for e in out.entities.values() if e.entity_type == 'protein'} == {'9606'}
    assert 'matrixdb:go_cellular_component' not in attributes(out.relations[2])


def test_only_human_reviewed_proteins_of_the_list_are_kept():
    rows = {row['accession']: row for row in parser.iter_ecm_proteins(opener(ECM_CSV))}
    # A reviewed protein of chicken, and an unreviewed human entry.
    assert rows['Q5ZKC9']['dataset'] == 'Swiss-Prot' and rows['Q5ZKC9']['taxon'] == '9031'
    assert matrixdb.map_ecm_protein(rows['Q5ZKC9']) is None
    assert rows['A0A024R1N1']['dataset'] == 'UniProtKB unreviewed (TrEMBL)'
    assert matrixdb.map_ecm_protein(rows['A0A024R1N1']) is None
    assert matrixdb.map_ecm_protein(rows['A8K2U0']) is not None
    # An unknown dataset is an error also for a protein that would not be kept.
    with pytest.raises(ValueError, match='UniProt dataset'):
        matrixdb.map_ecm_protein({**rows['Q5ZKC9'], 'dataset': 'PIR'})


@pytest.mark.parametrize(
    'text',
    [
        'protein_id,species (NCBI),go_terms,uniprot_keywords,dataset,interaction_count\n',
        ECM_CSV.replace('P39060', 'P39060X', 1),
        ECM_CSV + 'P05067,9606,,,Swiss-Prot,3\n',
        ECM_CSV.replace(',157', ',many'),
        ECM_CSV.replace('protein_id', 'accession'),
    ],
)
def test_ecm_parser_rejects_what_it_cannot_read(text):
    with pytest.raises(ValueError):
        list(parser.iter_ecm_proteins(opener(text)))


def test_unknown_selection_criterion_raises():
    row = next(parser.iter_ecm_proteins(opener(ECM_CSV)))
    with pytest.raises(ValueError, match='Unknown MatrixDB GO term'):
        matrixdb.map_ecm_protein({**row, 'go_terms': ['Extracellular exosome']})
    with pytest.raises(ValueError, match='Unknown MatrixDB UniProt keyword'):
        matrixdb.map_ecm_protein({**row, 'keywords': ['Secreted']})
    with pytest.raises(ValueError, match='UniProt dataset'):
        matrixdb.map_ecm_protein({**row, 'dataset': 'PIR'})


def test_list_term_is_the_object_of_the_protein_statements():
    terms = extract('terms', list(parser.iter_list_terms(opener(ECM_CSV))))
    proteins = extract('ecm_proteins', list(parser.iter_ecm_proteins(opener(ECM_CSV))))
    (term,) = terms.entities.values()
    objects = [proteins.entities[r.object_entity_key] for r in proteins.relations]
    assert {(e.entity_type, e.namespace, e.identifier) for e in objects} == {
        (term.entity_type, term.namespace, term.identifier)
    } == {('ontology_class', 'matrixdb_protein_list', 'ecm_proteins')}
    record = matrixdb.map_list_term({'list': 'ecm_proteins'})
    assert (Namespace.NAME, 'ECM and secreted proteins') in record.identifiers


def test_nothing_is_asserted_beyond_the_list():
    proteins = extract('ecm_proteins', list(parser.iter_ecm_proteins(opener(ECM_CSV))))
    # The GO terms and keywords are attributes of the statement: no entity for
    # them and no location is asserted.
    assert {e.entity_type for e in proteins.entities.values()} == {'protein', 'ontology_class'}
    assert {r.predicate for r in proteins.relations} == {'associated_with'}


def test_the_vocabulary_dataset_is_called_terms():
    assert list(matrixdb.resource.datasets()) == ['terms', 'ecm_proteins', 'interactions']
