"""CancerDrugs_DB drug-target table: parser, mapper and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import cancerdrugsdb
from pypath.inputs_v2.parsers.cancerdrugsdb import iter_drug_targets

DB = '<a href="https://go.drugbank.com/drugs/{0}">{0}</a>'
CH = '<a href="https://www.ebi.ac.uk/chembl/compound_report_card/{0}">{0}</a>'
HEADER = (
    'Product\tEMA\tFDA\tEN\tOther\tWHO\tYear\tGeneric\tDrugBank ID\tATC\tChEMBL'
    '\tIndications\tTargets\tLast Update\t'
)
# Rows copied from the source file (CRLF line ends, a tab after the last
# column); long indication and target cells are shortened.
ROWS = [
    # 0: a plain row.
    'Abemaciclib\tY\tY\tN\t\tN\t2017\tN\t' + DB.format('DB12001') + '\tL01EF03\t'
    + CH.format('CHEMBL3301610')
    + '\tAdvanced Breast Cancer; Metastatic Breast Cancer\tKRAS; ERBB2; TP53\t2025-04-02\t',
    # 1: no year, an agency under "Other".
    'Aclarubicin\tN\tN\tN\tJapan\tN\t\tN\t' + DB.format('DB11617') + '\tL01DB04\t'
    + CH.format('CHEMBL502620') + '\t\tTOP2A; CSF2\t2025-12-16\t',
    # 2: no targets.
    'Afamitresgene Autoleucel\tN\tY\tN\t\tN\t2024\t\t' + DB.format('DB18592') + '\t\t'
    + CH.format('CHEMBL4802254') + '\tUnresectable synovial sarcoma\t\t2024-10-30\t',
    # 3: two DrugBank identifiers, no targets.
    'Avutometinib, Defactinib\tN\tY\tN\t\tN\t2025\t\t' + DB.format('DB15254') + '; '
    + DB.format('DB12282') + '\t\t\tOvarian cancer\t\t2025-12-16\t',
    # 4: a class-level ATC code.
    'Imlunestrant\tN\tY\tN\t\tN\t2025\t\t' + DB.format('DB19043') + '\tL02BA\t'
    + CH.format('CHEMBL5095183') + '\tBreast cancer\tBCYRN1P2\t2025-12-16\t',
    # 5: DrugBank placeholder, an HTML character reference in the indications.
    'Nab-Paclitaxel\tY\tY\tN\t\tN\t2005\tY\t'
    '<a href="https://go.drugbank.com/drugs/Not found in DrugBank">Not found in DrugBank</a>'
    '\tL01CD01\t' + CH.format('CHEMBL428647')
    + "\tKaposi&#39;s sarcoma\tTUBA4A; TP53\t2025-04-02\t",
    # 6: a ChEMBL identifier with a space before it; several ATC codes.
    'Quizartinib\tY\tY\tN\t\tN\t2023\tN\t' + DB.format('DB12874') + '\tL01EX11; S01BA01\t'
    '<a href="https://www.ebi.ac.uk/chembl/compound_report_card/ CHEMBL576982"> CHEMBL576982</a>'
    '\tNewly Diagnosed AML\tFLT3; KIT\t2025-04-02\t',
]


def parse(rows=ROWS, header=HEADER):
    text = '\r\n'.join([header, *rows]) + '\r\n'
    return list(iter_drug_targets(SimpleNamespace(result=io.StringIO(text))))


def extract(records):
    out = SilverExtractor('cancerdrugsdb', 'drug_targets')
    for index, record in enumerate(records):
        mapped = cancerdrugsdb.resource.drug_targets.mapper(record)
        if mapped is not None:
            out.process_record(mapped, record, f'fixture:{index}', index)
    return out


def test_parser_gives_one_record_per_drug_and_target():
    records = parse()
    # Two rows have no targets and yield nothing: 3 + 2 + 1 + 2 + 2 targets.
    assert len(records) == 10
    first = records[0]
    assert first['product'] == 'Abemaciclib'
    assert first['target'] == 'KRAS'
    assert first['targets'] == ['KRAS', 'ERBB2', 'TP53']
    assert first['drugbank_ids'] == ['DB12001']
    assert first['chembl_ids'] == ['CHEMBL3301610']
    assert first['atc_codes'] == ['L01EF03']
    assert (first['ema'], first['fda'], first['en'], first['other'], first['who']) == ('Y', 'Y', 'N', '', 'N')
    assert (first['year'], first['generic'], first['last_update']) == ('2017', 'N', '2025-04-02')
    assert [r['product'] for r in records].count('Afamitresgene Autoleucel') == 0
    aclarubicin = [r for r in records if r['product'] == 'Aclarubicin']
    assert aclarubicin[0]['other'] == 'Japan' and aclarubicin[0]['year'] == ''


def test_parser_decodes_identifier_cells_and_text():
    by_product = {r['product']: r for r in parse()}
    # The placeholder for a drug DrugBank does not list is not an identifier.
    assert by_product['Nab-Paclitaxel']['drugbank_ids'] == []
    assert by_product['Nab-Paclitaxel']['indications'] == "Kaposi's sarcoma"
    # The space inside the ChEMBL link is not part of the identifier.
    assert by_product['Quizartinib']['chembl_ids'] == ['CHEMBL576982']
    assert by_product['Quizartinib']['atc_codes'] == ['L01EX11', 'S01BA01']
    # Two DrugBank identifiers of a combination row are both kept.
    two = parse([ROWS[3].replace('\t\tOvarian cancer\t\t', '\t\tOvarian cancer\tXYZ\t')])
    assert two[0]['drugbank_ids'] == ['DB15254', 'DB12282']


def test_parser_rejects_changed_table():
    with pytest.raises(ValueError, match='lacks columns'):
        parse(header=HEADER.replace('Targets', 'Genes'))
    with pytest.raises(ValueError, match='holds no drugs'):
        parse(rows=[])
    with pytest.raises(ValueError, match='unexpected DrugBank ID'):
        parse([ROWS[0].replace('DB12001</a>', 'DB1200</a>')])
    with pytest.raises(ValueError, match='not a link'):
        parse([ROWS[0].replace('<a href="https://go.drugbank.com/drugs/DB12001">DB12001</a>', 'DB12001')])
    with pytest.raises(ValueError, match='fewer fields'):
        parse(['Abemaciclib\tY\tY'])


def test_drug_is_associated_with_each_target_gene():
    record = cancerdrugsdb.resource.drug_targets.mapper(parse()[0])
    assert slot_name(record.predicate) == 'associated_with'
    assert record.subject.identifiers == [
        (Namespace.DRUGBANK, 'DB12001'),
        (Namespace.CHEMBL, 'CHEMBL3301610'),
        ('atc', 'L01EF03'),
        (Namespace.NAME, 'Abemaciclib'),
    ]
    assert record.object.identifiers == [(Namespace.GENESYMBOL, 'KRAS')]

    out = extract(parse()[:3])
    assert [r.predicate for r in out.relations] == ['associated_with'] * 3
    assert len({r.relation_key for r in out.relations}) == 3
    types = {e.entity_type for e in out.entities.values()}
    assert types == {'chemical_entity', 'protein'}
    # No species on the drug, and none stated for the targets.
    assert {e.taxon for e in out.entities.values()} == {None}


def test_identifier_systems_stay_apart_and_nothing_is_invented():
    out = extract(parse())
    chemicals = {e.label: e for e in out.entities.values() if e.entity_type == 'chemical_entity'}
    namespaces = {
        name: {str(i['ns']) for i in e.identifiers} for name, e in chemicals.items()
    }
    assert namespaces['Abemaciclib'] == {'drugbank', 'chembl', 'atc', 'name'}
    # The placeholder gives no DrugBank identifier; a class-level ATC code
    # identifies no chemical.
    assert namespaces['Nab-Paclitaxel'] == {'chembl', 'atc', 'name'}
    assert namespaces['Imlunestrant'] == {'drugbank', 'chembl', 'name'}
    assert {i['id'] for i in chemicals['Quizartinib'].identifiers if str(i['ns']) == 'atc'} == {
        'L01EX11',
        'S01BA01',
    }
    # No mechanism, direction, approval or indication is asserted.
    assert len(out.relations) == len(parse())
    for relation in out.relations:
        assert relation.predicate == 'associated_with'
        assert relation.annotations == []
    # A drug without targets has no relation and no entity.
    assert 'Afamitresgene Autoleucel' not in chemicals
    assert 'Avutometinib, Defactinib' not in chemicals
