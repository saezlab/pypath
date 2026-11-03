"""kinase.com kinome: reader, mapper and extractor on rows of Kincat_Hsap.08.02.xls."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import KnowledgeOrigin

from pypath.inputs_v2 import kinasedotcom
from pypath.inputs_v2.parsers import kinasedotcom as parser

HEADER = [
    'Name', 'SK', 'Group', 'Family', 'Subfamily', 'Pseudogene?', 'Header',
    'cDNA', 'Protein', 'Kinase Domain', 'Notes', 'KD2_Sequence', 'KD2_Group',
    'KD2_Family', 'KD2_Subfamily', 'Additional Aliases', 'Entrez_GeneID',
    'RefSeq cDNA (Jan 07)', 'Refseq mRNA', 'gi_mRNA', 'Refseq protein',
    'gi_protein', 'Old_Name', 'Entrez_Symbol',
]


def kinome_row(name, sk, group, family, subfamily, pseudogene, kd2, entrez, symbol):
    row = dict.fromkeys(HEADER, '')
    row.update(
        {
            'Name': name, 'SK': sk, 'Group': group, 'Family': family,
            'Subfamily': subfamily, 'Pseudogene?': pseudogene,
            'KD2_Group': kd2[0], 'KD2_Family': kd2[1], 'KD2_Subfamily': kd2[2],
            'Entrez_GeneID': entrez, 'Entrez_Symbol': symbol,
        }
    )
    return [row[column] for column in HEADER]


# Rows of the table. A cell without a value is an empty string, a blank
# second-domain group a space, and the Entrez ID a number.
NONE = (' ', '', '')
ROWS = [
    kinome_row('AKT1', 'SK018', 'AGC', 'Akt', '', 'N', NONE, 207.0, 'AKT1'),
    kinome_row('CRIK', 'SK695', 'AGC', 'DMPK', 'CRIK', 'N', NONE, 11113.0, 'CIT'),
    kinome_row('DMPK1', 'SK111', 'AGC', 'DMPK', 'GEK', 'N', NONE, 1760.0, 'DMPK'),
    kinome_row('RSK1', 'SK336', 'AGC', 'RSK', 'RSK', 'N', ('CAMK', 'RSKb', 'RSKb'), 6196.0, 'RPS6KA2'),
    kinome_row('Obscn', 'SK601', 'CAMK', 'Trio', '', 'N', ('CAMK', 'Trio', ''), 84033.0, 'OBSCN'),
    kinome_row('SgK110', 'SK592', 'Other', 'NKF1', '', 'N', NONE, '', ''),
    kinome_row('CK2a1-rs', 'SK750', 'CMGC', 'CK2', '', 'R', NONE, '', ''),
    kinome_row('MARKps30', 'SK599', 'CAMK', 'CAMKL', 'MARK', 'Y', NONE, '', ''),
]


def workbook(rows=ROWS, header=HEADER):
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.title = 'Sheet1'
    for row in [header, *rows]:
        sheet.append(row)
    data = io.BytesIO()
    book.save(data)
    return SimpleNamespace(result=io.BytesIO(data.getvalue()))


def tree():
    return {t['name']: t for t in parser.iter_kinase_tree(workbook())}


def extract(trees):
    out = SilverExtractor('kinasedotcom', 'families')
    for index, node in enumerate(trees):
        out.process_record(kinasedotcom.map_family(node), node, f'fixture:{index}', index)
    return out


def test_reader_builds_the_group_family_subfamily_tree():
    groups = tree()
    assert list(groups) == ['AGC', 'CAMK', 'Other']
    agc = groups['AGC']
    assert [c['name'] for c in agc['children']] == ['Akt', 'DMPK', 'RSK']
    akt, dmpk, rsk = agc['children']
    # AKT1 has no subfamily: it is a member of the family.
    assert akt['members'] == [
        {'sk': 'SK018', 'name': 'AKT1', 'entrez': '207', 'symbol': 'AKT1', 'kinase_domain': '1'}
    ]
    assert [(c['level'], c['name']) for c in dmpk['children']] == [
        ('subfamily', 'CRIK'),
        ('subfamily', 'GEK'),
    ]
    # The second kinase domain places RSK1 in another group.
    assert groups['CAMK']['children'][0]['name'] == 'RSKb'
    second = groups['CAMK']['children'][0]['children'][0]
    assert (second['level'], second['name']) == ('subfamily', 'RSKb')
    assert [(m['name'], m['kinase_domain']) for m in second['members']] == [('RSK1', '2')]
    # Both domains of Obscn are in the same family.
    trio = next(c for c in groups['CAMK']['children'] if c['name'] == 'Trio')
    assert [(m['name'], m['kinase_domain']) for m in trio['members']] == [
        ('Obscn', '1'),
        ('Obscn', '2'),
    ]
    # A kinase the table gives no Entrez Gene ID is kept with its accession.
    assert groups['Other']['children'][0]['members'][0]['entrez'] is None


def test_reader_leaves_out_pseudogenes_and_rejects_bad_input():
    names = {m['name'] for g in tree().values() for m in _members(g)}
    assert 'MARKps30' not in names and 'CK2a1-rs' not in names
    assert 'CAMKL' not in {c['name'] for c in tree()['CAMK']['children']}
    bad = [kinome_row('X', 'SK1', 'AGC', 'Akt', '', '?', NONE, 1.0, 'X')]
    with pytest.raises(ValueError, match='Unknown kinase.com Pseudogene'):
        list(parser.iter_kinase_tree(workbook(bad)))
    gap = [kinome_row('X', 'SK1', 'AGC', '', 'Sub', 'N', NONE, 1.0, 'X')]
    with pytest.raises(ValueError, match='Malformed kinase.com classification'):
        list(parser.iter_kinase_tree(workbook(gap)))
    with pytest.raises(ValueError, match='Table lacks expected columns'):
        list(parser.iter_kinase_tree(workbook(header=HEADER[:-1], rows=[r[:-1] for r in ROWS])))
    with pytest.raises(ValueError, match='no rows'):
        list(parser.iter_kinase_tree(workbook(rows=[])))


def _members(node):
    yield from node['members']
    for child in node['children']:
        yield from _members(child)


def test_families_are_nested_with_their_kinases_as_members():
    out = extract(tree().values())
    families = {
        (e.namespace, e.identifier): e
        for e in out.entities.values()
        if e.entity_type == 'gene_family'
    }
    # A name is only unique within its level: the family RSK has the subfamily RSK.
    assert {'kinasedotcom_group', 'kinasedotcom_family', 'kinasedotcom_subfamily'} == {
        ns for ns, _ in families
    }
    assert ('kinasedotcom_family', 'RSK') in families
    assert ('kinasedotcom_subfamily', 'RSK') in families
    assert {e.taxon for e in families.values()} == {'9606'}

    by_key = {key: (e.namespace, e.identifier) for key, e in out.entities.items()}
    edges = {
        (by_key[r.subject_entity_key], by_key[r.object_entity_key])
        for r in out.relations
        if r.predicate == 'has_member'
    }
    assert (('kinasedotcom_group', 'AGC'), ('kinasedotcom_family', 'DMPK')) in edges
    assert (('kinasedotcom_family', 'DMPK'), ('kinasedotcom_subfamily', 'GEK')) in edges
    assert (('kinasedotcom_family', 'RSK'), ('kinasedotcom_subfamily', 'RSK')) in edges
    # The kinase is a member of the narrowest level only, and has the
    # Entrez Gene ID as identifier.
    assert (('kinasedotcom_subfamily', 'GEK'), (Namespace.ENTREZ, '1760')) in edges
    assert (('kinasedotcom_group', 'AGC'), (Namespace.ENTREZ, '1760')) not in edges
    assert (('kinasedotcom_family', 'DMPK'), (Namespace.ENTREZ, '1760')) not in edges


def test_kinase_gene_identifiers_and_domain_attribute():
    out = extract(tree().values())
    genes = {e.identifier: e for e in out.entities.values() if e.entity_type == 'gene'}
    assert {e.taxon for e in genes.values()} == {'9606'}
    crik = genes['11113']
    ids = {(i['ns'], i['id']) for i in crik.identifiers}
    # CRIK is kinase.com's name for CIT: a synonym, not the label.
    assert ids == {
        (Namespace.ENTREZ, '11113'),
        (Namespace.GENESYMBOL, 'CIT'),
        ('kinasedotcom', 'SK695'),
        (Namespace.SYNONYM, 'CRIK'),
    }
    # A kinase without an Entrez ID is identified by its kinase.com accession.
    sgk = next(e for e in out.entities.values() if e.entity_type == 'gene' and e.identifier == 'SK592')
    assert sgk.namespace == 'kinasedotcom'
    assert {(i['ns'], i['id']) for i in sgk.identifiers} == {
        ('kinasedotcom', 'SK592'),
        (Namespace.NAME, 'SgK110'),
    }
    domains = sorted(
        a['value']
        for r in out.relations
        for a in r.annotations
        if a['term'] == 'kinasedotcom:kinase_domain'
        and out.entities[r.object_entity_key].identifier == '84033'
    )
    assert domains == ['1', '2']
    # Pseudogenes are not entities.
    assert not any(e.identifier in {'SK599', 'SK750'} for e in out.entities.values())


def test_knowledge_origin_is_computationally_predicted():
    assert tuple(kinasedotcom.config.knowledge_origin) == (KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,)


def test_the_dataset_of_families_is_called_families():
    assert list(kinasedotcom.resource.datasets()) == ['families']
