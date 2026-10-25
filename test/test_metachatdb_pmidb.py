"""MetaChatDB and PMI-DB: parsers on small fixtures and the
Biolink statements their rows become."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import annotation_term, annotation_value
from omnipath_core.interaction_profiles import interaction_label


def opener(text):
    return SimpleNamespace(result=io.StringIO(text))


def extract(source, record):
    assert record is not None
    out = SilverExtractor(source, 'fixture')
    out.process_record(record, {}, 'fixture:1', 1)
    return out


def relation_annotations(out):
    """Term -> values of the main relation (the first one with a statement)."""
    relation = next(r for r in out.relations if r.predicate in ('interacts_with', 'affects'))
    values = {}
    for annotation in relation.annotations:
        if annotation.get('scope', 'relation') != 'relation':
            continue
        values.setdefault(annotation['term'], []).append(annotation['value'])
    return relation, values


def entity(out, key):
    return out.entities[key]


def ids(observation):
    return {(i['ns'], i['id']) for i in observation.identifiers}


# MetaChatDB ------------------------------------------------------------------

METACHAT_HEADER = (
    'Sensor.Gene\tSensor.Name\tSensor.Type\tSensor.Pathway\tHMDB.ID\t'
    'Metabolite.Name\tMetabolite.transporter\tMetabolite.Class\t'
    'Metabolite.Pathway\tLong.Range.Channel\tEvidences\tSources\n'
)


def metachat_row(**values):
    row = {
        'Sensor.Gene': 'CD36',
        'Sensor.Name': 'Platelet glycoprotein 4 (EC 1.2.3.4) (Fatty acid translocase) (FAT)',
        'Sensor.Type': 'receptor',
        'Sensor.Pathway': 'WP1: A pathway',
        'HMDB.ID': 'HMDB0001043',
        'Metabolite.Name': 'Arachidonic acid',
        'Metabolite.transporter': '',
        'Metabolite.Class': 'Lipids and lipid-like molecules',
        'Metabolite.Pathway': 'NA',
        'Long.Range.Channel': 'Blood,Urine',
        'Evidences': '123; https://en.wikipedia.org/wiki/CD36; Recon2',
        'Sources': 'MetaChatDB; MRCLinkdb',
        'chemont_ids': ['CHEMONTID:0000012'],
    }
    row.update(values)
    return row


def test_metachatdb_keeps_own_rows_and_splits_sensor_types():
    from pypath.inputs_v2.metachatdb import iter_metachatdb

    own = metachat_row(**{'Sensor.Type': 'receptor,transporter'})
    imported = metachat_row(**{'Sources': 'Mebocost; MRCLinkdb'})
    columns = METACHAT_HEADER.strip().split('\t')
    text = METACHAT_HEADER + ''.join(
        '\t'.join(row[c] for c in columns) + '\n' for row in (own, imported)
    )
    rows = list(iter_metachatdb(opener(text)))
    assert [row['sensor_type'] for row in rows] == ['receptor', 'transporter']


def test_metachatdb_rejects_an_unknown_sensor_type():
    from pypath.inputs_v2.metachatdb import iter_metachatdb

    row = metachat_row(**{'Sensor.Type': 'enzyme'})
    columns = METACHAT_HEADER.strip().split('\t')
    with pytest.raises(ValueError, match='enzyme'):
        list(iter_metachatdb(opener(METACHAT_HEADER + '\t'.join(row[c] for c in columns) + '\n')))


def test_metachatdb_receptor_binds_its_metabolite():
    from pypath.inputs_v2.metachatdb import interactions_schema

    out = extract('metachatdb', interactions_schema('9606')(metachat_row(sensor_type='receptor')))
    relation, values = relation_annotations(out)
    assert relation.predicate == 'interacts_with'
    # Unqualified, as MEBOCOST and MRClinksDB state a receptor and its ligand.
    assert interaction_label(relation.predicate, relation.annotations) is None
    chemical = entity(out, relation.subject_entity_key)
    sensor = entity(out, relation.object_entity_key)
    assert chemical.entity_type == 'chemical_entity'
    assert ('hmdb', 'HMDB0001043') in ids(chemical)
    assert sensor.entity_type == 'protein' and sensor.taxon == '9606'
    assert ('genesymbol', 'CD36') in ids(sensor)
    assert ('name', 'Platelet glycoprotein 4') in ids(sensor)
    assert ('synonym', 'Fatty acid translocase') in ids(sensor)
    assert not any('EC 1.2.3.4' in value for _, value in ids(sensor))
    assert values['publications'] == ['PMID:123']
    assert values['source_record_urls'] == ['https://en.wikipedia.org/wiki/CD36']
    assert values['supporting_data_source'] == ['MRCLinkdb', 'Recon2']
    assert values['metachatdb:sensor_type'] == ['receptor']
    assert values['metachatdb:long_range_channel'] == ['Blood', 'Urine']
    assert 'metachatdb:metabolite_pathway' not in values  # NA
    assert 'object_aspect_qualifier' not in values
    # The ClassyFire class is a ChemOnt term of the chemical.
    assert any(r.predicate == 'associated_with' for r in out.relations)


def test_metachatdb_transporter_moves_its_metabolite():
    from pypath.inputs_v2.metachatdb import interactions_schema

    out = extract('metachatdb', interactions_schema('10090')(metachat_row(sensor_type='transporter')))
    relation, values = relation_annotations(out)
    assert relation.predicate == 'affects'
    assert interaction_label(relation.predicate, relation.annotations) == 'Transports'
    assert entity(out, relation.subject_entity_key).entity_type == 'protein'
    assert entity(out, relation.subject_entity_key).taxon == '10090'
    assert entity(out, relation.object_entity_key).entity_type == 'chemical_entity'


def test_metachatdb_protein_names_leave_out_other_chains():
    from pypath.inputs_v2.metachatdb import protein_name, protein_synonyms

    value = (
        'Adhesion G protein-coupled receptor B1 (Brain-specific angiogenesis '
        'inhibitor 1) [Cleaved into: Vasculostatin-40 (Vstat40)]'
    )
    assert protein_name(value) == 'Adhesion G protein-coupled receptor B1'
    assert protein_synonyms(value) == ['Brain-specific angiogenesis inhibitor 1']
    assert protein_synonyms('Channel (Na(+) channel) (EC 7.1.1.1)') == ['Na(+) channel']


# PMI-DB ----------------------------------------------------------------------

PMIDB_INTERACTIONS = (
    'KEGG\tUniprot_KB_id\tinteraction\tprobe\tmean\tliterature\tTechnology\t'
    'Sample source\tMS Quality control method\tCandidate selection cutoff\r\n'
    'C00002\tP00350\tTRUE\tNA\tNA\tPMID:29307493\tLiP-SMap\tEscherichia coli\tDDA and DIA\tFC>2,Q<0.01\r\n'
    'C00002\tP00363\tFALSE\tNA\tNA\tPMID:29307493\tLiP-SMap\tEscherichia coli\tDDA and DIA\tFC>2,Q<0.01\r\n'
    'C00003\tQ9Y6K9\tTRUE\tUV probe\t1.5\tPMID:26524241\tLipid probes\thuman\tSILAC\tFC>=3.0\r\n'
)
PMIDB_METABOLITES = (
    "C00002\tHMDB00538\tAdenosine 5'-triphosphate\thttps://www.genome.jp/x\tNC1=NC=NC2\r\n"
    'C00003\tHMDB00902\tNAD\thttps://www.genome.jp/y\tNC(=O)C1\r\n'
)
PMIDB_PROTEINS = (
    'P00350\tgnd\tNA\thttps://www.uniprot.org/uniprot/P00350\r\n'
    'Q9Y6K9\tIKBKG,\tGO:0005515\thttps://www.uniprot.org/uniprot/Q9Y6K9\r\n'
)


@pytest.fixture
def pmidb_rows(monkeypatch):
    from pypath.inputs_v2 import pmidb

    files = {'metabolites_infor': PMIDB_METABOLITES, 'protein_infor': PMIDB_PROTEINS}
    monkeypatch.setattr(
        pmidb,
        '_download',
        lambda name: SimpleNamespace(open=lambda **_: opener(files[name])),
    )
    return list(pmidb.iter_pmidb(opener(PMIDB_INTERACTIONS)))


def test_pmidb_keeps_interacting_pairs_with_their_identifiers(pmidb_rows):
    assert [(r['KEGG'], r['Uniprot_KB_id']) for r in pmidb_rows] == [
        ('C00002', 'P00350'),
        ('C00003', 'Q9Y6K9'),
    ]
    assert pmidb_rows[0]['taxon_id'] == '562'
    assert pmidb_rows[0]['hmdb'] == 'HMDB00538'
    assert pmidb_rows[1]['gene'] == 'IKBKG,'


def test_pmidb_rejects_an_unknown_sample_source(monkeypatch):
    from pypath.inputs_v2 import pmidb

    monkeypatch.setattr(pmidb, '_download', lambda name: SimpleNamespace(open=lambda **_: opener('')))
    text = PMIDB_INTERACTIONS.replace('\thuman\t', '\trat\t')
    with pytest.raises(ValueError, match='rat'):
        list(pmidb.iter_pmidb(opener(text)))


def test_pmidb_pair_is_a_binding_with_study_context(pmidb_rows):
    from pypath.inputs_v2.pmidb import interactions_schema

    out = extract('pmidb', interactions_schema(pmidb_rows[1]))
    relation, values = relation_annotations(out)
    assert relation.predicate == 'interacts_with'
    chemical = entity(out, relation.subject_entity_key)
    protein = entity(out, relation.object_entity_key)
    assert {('kegg', 'C00003'), ('hmdb', 'HMDB00902'), ('smiles', 'NC(=O)C1'), ('name', 'NAD')} <= ids(chemical)
    assert ids(protein) == {('uniprot', 'Q9Y6K9'), ('genesymbol', 'IKBKG')}
    assert protein.taxon == '9606'
    assert values['publications'] == ['PMID:26524241']
    assert values['pmidb:technology'] == ['Lipid probes']
    assert values['pmidb:probe'] == ['UV probe']
    assert values['pmidb:candidate_selection_cutoff'] == ['FC>=3.0']


def test_pmidb_na_context_is_left_out(pmidb_rows):
    from pypath.inputs_v2.pmidb import interactions_schema

    _, values = relation_annotations(extract('pmidb', interactions_schema(pmidb_rows[0])))
    assert 'pmidb:probe' not in values and 'pmidb:mean' not in values


@pytest.mark.parametrize('module', ['metachatdb', 'pmidb'])
def test_module_imports_and_names_its_datasets(module):
    from importlib import import_module

    assert import_module('pypath.inputs_v2.' + module).resource.datasets()
