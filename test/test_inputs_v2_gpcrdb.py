"""GPCRdb: reader of the outline and the family tree, mapper and extractor."""

import inspect
import io
import json
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import gpcrdb
from pypath.inputs_v2.parsers import gpcrdb as gpcrdb_parser

# The first receptor, the angiotensin AT1 receptor (rat and mouse have two
# genes, a bar separates their accessions) and the CRF1 receptor (an isoform
# in the human accession), copied from proteins_and_families.txt.
HTR1A = (
    '"gpcr","1","5-Hydroxytryptamine receptors","1","5-HT<sub>1A</sub> receptor",'
    '"","","","ADRBRL1|5-HT1A|ADRB2RL1|serotonin receptor 1A|Gpcr18|'
    '5-hydroxytryptamine (serotonin) receptor 1A, G protein-coupled","5286",'
    '"HTR1A","5-hydroxytryptamine receptor 1A","5q11.2-q13","NM_000524",'
    '"NP_000515","P08908","3350","2845","Htr1a","5-hydroxytryptamine receptor 1A",'
    '"2q16","NM_012585","NP_036717","P19327","24473","MGI:96273","Htr1a",'
    '"5-hydroxytryptamine (serotonin) receptor 1A","13 D2.1","NM_008308",'
    '"NP_032334","Q64264","15550"'
)
AT1 = (
    '"gpcr","6","Angiotensin receptors","34","AT<sub>1</sub> receptor","","","",'
    '"AG2S|AT2R1|AT2R1A","336","AGTR1","angiotensin II receptor type 1",'
    '"3q21-q25","NM_000685","NP_000676","P30556","185","2070|2071",'
    '"Agtr1a|Agtr1b","angiotensin II receptor, type 1a|angiotensin II receptor, '
    'type 1b","17q12|2q24","NM_031009|NM_030985","NP_112271|NP_112247",'
    '"P29089|P25095","81638|24180","MGI:87965|MGI:87964","Agtr1b|Agtr1a",'
    '"angiotensin II receptor, type 1b|angiotensin II receptor, type 1a",'
    '"3 7.6 cM|13 16.0 cM","NM_175086|NM_177322","NP_780295|NP_796296",'
    '"P29755|P29754","11608|11607"'
)
CRF1 = (
    '"gpcr","19","Corticotropin-releasing factor receptors","212",'
    '"CRF<sub>1</sub> receptor","","","","CRFR1|CRF-R","2357","CRHR1",'
    '"corticotropin releasing hormone receptor 1","17q21.31","NM_004382",'
    '"NP_004373","P34998-2","1394","61276","Crhr1","corticotropin releasing '
    'hormone receptor 1","10q32.1","NM_030999","NP_112261","P35353","58959",'
    '"MGI:88498","Crhr1","corticotropin releasing hormone receptor 1","11 E1",'
    '"NM_007762","NP_031788","P35347","12921"'
)
OUTLINE = f"""\
Class A (Rhodopsin) | gpcrdba
    Aminergic receptors
        5-Hydroxytryptamine receptors
            {HTR1A}
    Peptide receptors
        Angiotensin receptors
            {AT1}
Class B1 (Secretin) | gpcrdbb
    Peptide receptors
        Corticotropin-releasing factor receptors
            {CRF1}
"""
TREE = [
    {'slug': '000', 'name': 'Parent family', 'parent': None},
    {'slug': '001', 'name': 'Class A (Rhodopsin)', 'parent': {'slug': '000', 'name': 'Parent family'}},
    {'slug': '001_001', 'name': 'Aminergic receptors', 'parent': {'slug': '001', 'name': 'Class A (Rhodopsin)'}},
    {'slug': '001_001_001', 'name': '5-Hydroxytryptamine receptors', 'parent': {'slug': '001_001', 'name': 'Aminergic receptors'}},
    {'slug': '001_002', 'name': 'Peptide receptors', 'parent': {'slug': '001', 'name': 'Class A (Rhodopsin)'}},
    {'slug': '001_002_001', 'name': 'Angiotensin receptors', 'parent': {'slug': '001_002', 'name': 'Peptide receptors'}},
    {'slug': '002', 'name': 'Class B1 (Secretin)', 'parent': {'slug': '000', 'name': 'Parent family'}},
    {'slug': '002_001', 'name': 'Peptide receptors', 'parent': {'slug': '002', 'name': 'Class B1 (Secretin)'}},
    {'slug': '002_001_002', 'name': 'Corticotropin-releasing factor receptors', 'parent': {'slug': '002_001', 'name': 'Peptide receptors'}},
    # A class of the tree that the outline does not have.
    {'slug': '012', 'name': 'G-Protein', 'parent': {'slug': '000', 'name': 'Parent family'}},
]


class Tree:
    def __init__(self, nodes):
        self.nodes = nodes

    def open(self, **_kwargs):
        return SimpleNamespace(result=io.StringIO(json.dumps(self.nodes)))


def read(outline=OUTLINE, tree=TREE):
    return list(
        gpcrdb_parser.iter_families(
            SimpleNamespace(result=io.StringIO(outline)), tree=Tree(tree)
        )
    )


def extract(rows):
    out = SilverExtractor('gpcrdb', 'families')
    for index, row in enumerate(rows):
        out.process_record(gpcrdb.map_family(row), row, f'fixture:{index}', index)
    return out


def by_slug(rows):
    return {row['slug']: row for row in rows}


def test_reader_gives_each_level_with_its_slug_and_children():
    rows = by_slug(read())
    assert set(rows) == {
        '001', '001_001', '001_001_001', '001_002', '001_002_001',
        '002', '002_001', '002_001_002',
    }
    assert rows['001'] == {
        'slug': '001', 'name': 'Class A (Rhodopsin)', 'level': 'class',
        'children': [
            {'slug': '001_001', 'name': 'Aminergic receptors'},
            {'slug': '001_002', 'name': 'Peptide receptors'},
        ],
    }
    # The same name at different places of the tree has different slugs.
    assert rows['002_001']['name'] == rows['001_002']['name'] == 'Peptide receptors'
    assert rows['001_001']['children'] == [
        {'slug': '001_001_001', 'name': '5-Hydroxytryptamine receptors'}
    ]


def test_reader_gives_the_receptors_of_each_species_column():
    rows = by_slug(read())
    assert rows['001_001_001']['members'] == [
        {'accession': 'P08908', 'taxon': '9606', 'species': 'human'},
        {'accession': 'P19327', 'taxon': '10116', 'species': 'rat'},
        {'accession': 'Q64264', 'taxon': '10090', 'species': 'mouse'},
    ]
    # Two genes in a species: both accessions, each as a member.
    accessions = [
        (m['species'], m['accession']) for m in rows['001_002_001']['members']
    ]
    assert accessions == [
        ('human', 'P30556'), ('rat', 'P29089'), ('rat', 'P25095'),
        ('mouse', 'P29755'), ('mouse', 'P29754'),
    ]
    # The isoform number stays with the accession.
    assert rows['002_001_002']['members'][0]['accession'] == 'P34998-2'


def test_names_are_given_without_html_markup():
    # GPCRdb writes the subscript of "GABAB receptors" as markup, in the outline
    # and in the tree; the join to the tree is by the name as written.
    marked = 'GABA<sub>B</sub> receptors'
    tree = [
        {**node, 'name': marked} if node['slug'] == '001_001_001' else node
        for node in TREE
    ]
    outline = OUTLINE.replace('        5-Hydroxytryptamine receptors', f'        {marked}')
    rows = by_slug(read(outline, tree))
    assert rows['001_001_001']['name'] == 'GABAB receptors'
    assert rows['001_001']['children'] == [{'slug': '001_001_001', 'name': 'GABAB receptors'}]
    out = extract(list(rows.values()))
    names = {
        i['id'] for e in out.entities.values() for i in e.identifiers
        if i['ns'] == Namespace.NAME
    }
    assert 'GABAB receptors' in names and not any('<' in name for name in names)


def test_reader_rejects_a_node_that_is_not_in_the_tree_and_bad_lines():
    with pytest.raises(ValueError, match='not in the family tree'):
        read(tree=[node for node in TREE if node['slug'] != '001_002_001'])
    with pytest.raises(ValueError, match='indented by 6'):
        read(OUTLINE.replace('    Peptide receptors\n        Angiotensin', '      Peptide receptors\n        Angiotensin'))
    with pytest.raises(ValueError, match='accession'):
        read(OUTLINE.replace('"P08908"', '"5-HT1A"'))
    with pytest.raises(ValueError, match='nothing below'):
        read(OUTLINE.replace('            ' + CRF1 + '\n', ''))
    with pytest.raises(ValueError, match='empty'):
        read('')
    with pytest.raises(ValueError, match='empty'):
        read(tree=[])


def test_a_family_has_its_slug_name_and_receptors_as_members():
    rows = read()
    out = extract(rows)
    family = next(
        e for e in out.entities.values()
        if (e.namespace, e.identifier) == (gpcrdb.FAMILY, '001_001_001')
        and e.entity_type == 'protein_family'
    )
    assert {i['id'] for i in family.identifiers if i['ns'] == Namespace.NAME} == {
        '5-Hydroxytryptamine receptors'
    }
    members = {
        (
            out.entities[r.object_entity_key].entity_type,
            out.entities[r.object_entity_key].namespace,
            out.entities[r.object_entity_key].identifier,
            out.entities[r.object_entity_key].taxon,
        )
        for r in out.relations
        if r.predicate == 'has_member' and r.subject_entity_key == family.entity_key
    }
    assert members == {
        ('protein', Namespace.UNIPROT, 'P08908', '9606'),
        ('protein', Namespace.UNIPROT, 'P19327', '10116'),
        ('protein', Namespace.UNIPROT, 'Q64264', '10090'),
    }


def test_levels_are_nested_class_ligand_type_family():
    out = extract(read())
    names = {
        k: [i['id'] for i in e.identifiers if i['ns'] == Namespace.NAME][0]
        for k, e in out.entities.items()
        if e.entity_type == 'protein_family'
    }
    pairs = {
        (names[r.subject_entity_key], names[r.object_entity_key])
        for r in out.relations
        if r.predicate == 'has_member' and r.object_entity_key in names
    }
    assert pairs == {
        ('Class A (Rhodopsin)', 'Aminergic receptors'),
        ('Class A (Rhodopsin)', 'Peptide receptors'),
        ('Aminergic receptors', '5-Hydroxytryptamine receptors'),
        ('Peptide receptors', 'Angiotensin receptors'),
        ('Class B1 (Secretin)', 'Peptide receptors'),
        ('Peptide receptors', 'Corticotropin-releasing factor receptors'),
    }


def test_nothing_is_asserted_beyond_membership():
    out = extract(read())
    assert {r.predicate for r in out.relations} == {'has_member'}
    # The class name carries no letter of the numbering scheme.
    assert not any('gpcrdb' in i['id'] for e in out.entities.values() for i in e.identifiers)
    # A receptor is a protein of its species and nothing else is said of it.
    proteins = [e for e in out.entities.values() if e.entity_type == 'protein']
    assert len(proteins) == 11 and all(e.taxon for e in proteins)


def test_unknown_level_raises():
    with pytest.raises(ValueError, match='Unknown GPCRdb level'):
        gpcrdb.map_family({'slug': '1', 'name': 'x', 'level': 'order'})


def test_resource_datasets():
    assert set(gpcrdb.resource.datasets()) == {'families'}


def test_module_has_no_reading_code_and_the_dataset_binds_the_parser():
    # The family tree is read by the parser; the module only names the download.
    assert 'opener' not in inspect.getsource(gpcrdb)
    parser = gpcrdb.resource.families._raw_parser
    assert parser.func is gpcrdb_parser.iter_families
    assert parser.keywords == {'tree': gpcrdb.download_tree}
