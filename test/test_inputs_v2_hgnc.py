"""HGNC gene groups: parser, mapper and extractor on source rows."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import entity_type, slot_name
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import KnowledgeOrigin

from pypath.inputs_v2 import hgnc
from pypath.inputs_v2.parsers import hgnc as parser

# Rows copied from family.csv (the unread columns are not shown).
FAMILY = (
    '"id","abbreviation","name","external_note","pubmed_ids","desc_comment",'
    '"desc_label","desc_source","desc_go","typical_gene"\n'
    '"3","FSCN","Fascin family","NULL","21618240","Fascins are actin-binding '
    'proteins that cross-link filamentous actin.","NULL","Hashimoto et al, 2011",'
    '"NULL","FSCN1"\n'
    '"139","NULL","G protein-coupled receptors","NULL","NULL","Receptors which '
    'transduce extracellular signals across the cell membrane.",'
    '"G-protein coupled receptor","UniProt|http://purl.uniprot.org/keywords/297",'
    '"http://purl.uniprot.org/go/0004930","NULL"\n'
    '"639","","Mitochondrial respiratory chain complexes","You can view genes.",'
    '"NULL","NULL","NULL","NULL","NULL","NULL"\n'
    '"640","NULL","Mitochondrial complex I: NADH:ubiquinone oxidoreductase subunits",'
    '"Please note.","16756485,25209663","Complex I is an enzyme of the respiratory '
    'chains.","NADH dehydrogenase (ubiquinone)","Wikipedia","NULL","NULL"\n'
    '"642","UQCR","Mitochondrial complex III: ubiquinol-cytochrome c reductase '
    'complex subunits","NULL","NULL","The coenzyme Q : cytochrome c-oxidoreductase.",'
    '"NULL","NULL","NULL","NULL"\n'
)
ALIAS = (
    '"id","family_id","alias"\n'
    '"4278","3","Fascins"\n'
    '"2094","139","G protein coupled receptors"\n'
    '"2099","139","GPCRs"\n'
    '"502","642","Coenzyme Q : cytochrome c-oxidoreductase subunits"\n'
)
HIERARCHY = '"parent_fam_id","child_fam_id"\n"639","640"\n"639","642"\n'
MEMBERS = (
    '"hgnc_id","family_id"\n"11148","3"\n"3960","3"\n"3961","3"\n'
    '"2579","642"\n"7427","642"\n"99999999","642"\n'
)
SYMBOLS = (
    'HGNC ID\tApproved symbol\n'
    'HGNC:11148\tFSCN1\nHGNC:3960\tFSCN2\nHGNC:3961\tFSCN3\n'
    'HGNC:2579\tCYC1\nHGNC:7427\tMT-CYB\n'
)


class Source:
    def __init__(self, text):
        self.text = text

    def open(self, **_kwargs):
        return SimpleNamespace(result=self.text)


def groups(family=FAMILY, alias=ALIAS, hierarchy=HIERARCHY, members=MEMBERS, symbols=SYMBOLS):
    return list(
        parser.iter_groups(
            Source(family).open(),
            alias_download=Source(alias),
            hierarchy_download=Source(hierarchy),
            members_download=Source(members),
            symbols_download=Source(symbols),
        )
    )


def extract(rows):
    out = SilverExtractor('hgnc', 'families')
    for index, row in enumerate(rows):
        out.process_record(hgnc.map_group(row), row, f'fixture:{index}', index)
    return out


def test_the_dataset_of_families_is_called_families():
    assert list(hgnc.resource.datasets()) == ['families']


def test_parser_joins_the_tables():
    rows = {r['id']: r for r in groups()}
    assert list(rows) == ['3', '139', '639', '640', '642']
    fascin = rows['3']
    assert fascin['name'] == 'Fascin family'
    assert fascin['abbreviation'] == 'FSCN'
    assert fascin['aliases'] == ['Fascins']
    assert fascin['pmids'] == ['21618240']
    assert fascin['members'] == [
        {'hgnc_id': '11148', 'symbol': 'FSCN1'},
        {'hgnc_id': '3960', 'symbol': 'FSCN2'},
        {'hgnc_id': '3961', 'symbol': 'FSCN3'},
    ]
    assert rows['639']['children'] == ['640', '642']
    assert rows['639']['members'] == []
    assert rows['640']['pmids'] == ['16756485', '25209663']
    # NULL and the empty string are both no value.
    assert rows['639']['abbreviation'] is None
    assert rows['139']['abbreviation'] is None
    assert rows['3']['description'].startswith('Fascins are actin-binding')
    assert rows['139']['go'] == 'GO:0004930'
    assert rows['3']['go'] is None
    # A gene without an approved symbol keeps its HGNC id.
    assert rows['642']['members'][-1] == {'hgnc_id': '99999999', 'symbol': None}


@pytest.mark.parametrize(
    'tables',
    [
        {'family': FAMILY.split('\n')[0] + '\n'},
        {'family': ''},
        {'family': FAMILY + FAMILY.split('\n')[1] + '\n'},
        {'members': '"hgnc_id","family_id"\n"11148","1"\n'},
        {'hierarchy': '"parent_fam_id","child_fam_id"\n"639","1"\n'},
        {'members': '"hgnc_id","family_id"\n"11148","3"\n"11148","3"\n'},
        {'symbols': 'HGNC ID\tApproved symbol\nA1BG\tA1BG\n'},
        {'family': FAMILY.replace('"21618240"', '"PMID21618240"')},
        {'family': FAMILY.replace('http://purl.uniprot.org/go/0004930', 'GO:0004930')},
        {'alias': '"id","family_id"\n"1","3"\n'},
    ],
)
def test_parser_rejects_tables_it_cannot_join(tables):
    with pytest.raises(ValueError):
        groups(**tables)


def test_group_is_a_gene_family_with_genes_and_subgroups_as_members():
    rows = {r['id']: r for r in groups()}
    record = hgnc.map_group(rows['639'])
    assert entity_type(record.type) == 'gene_family'
    assert record.identifiers == [
        ('hgnc_group', '639'),
        (Namespace.NAME, 'Mitochondrial respiratory chain complexes'),
    ]
    assert [slot_name(m.predicate) for m in record.membership] == ['has_member'] * 2
    assert [m.member.identifiers for m in record.membership] == [
        [('hgnc_group', '640')],
        [('hgnc_group', '642')],
    ]

    fascin = hgnc.map_group(rows['3'])
    assert fascin.identifiers == [
        ('hgnc_group', '3'),
        (Namespace.NAME, 'Fascin family'),
        (Namespace.SYNONYM, 'FSCN'),
        (Namespace.SYNONYM, 'Fascins'),
    ]
    assert [m.member.identifiers for m in fascin.membership] == [
        [(Namespace.HGNC, '11148'), (Namespace.GENESYMBOL, 'FSCN1')],
        [(Namespace.HGNC, '3960'), (Namespace.GENESYMBOL, 'FSCN2')],
        [(Namespace.HGNC, '3961'), (Namespace.GENESYMBOL, 'FSCN3')],
    ]


def test_extractor_yields_has_member_edges_and_group_attributes():
    rows = groups()
    out = extract(rows)
    # 6 member genes (one without a symbol), 2 subgroups and the GO term.
    assert len(out.relations) == 6 + 2 + 1
    assert sum(r.predicate == 'has_member' for r in out.relations) == 8
    genes = [e for e in out.entities.values() if e.entity_type == 'gene']
    assert len(genes) == 6
    assert {e.taxon for e in genes} == {'9606'}
    assert {e.namespace for e in genes} == {Namespace.HGNC}
    gpcr = next(e for e in out.entities.values() if e.identifier == '139')
    assert gpcr.entity_type == 'gene_family'


def test_go_term_is_the_same_entity_as_the_go_ontology_term():
    out = extract(groups())
    goal = [r for r in out.relations if r.predicate == 'associated_with']
    assert len(goal) == 1
    term = out.entities[goal[0].object_entity_key]
    assert (term.entity_type, term.namespace, term.identifier) == (
        'ontology_class',
        Namespace.GO,
        'GO:0004930',
    )


def test_hierarchy_is_not_expanded_into_members_or_pairs():
    out = extract(groups())
    families = {
        key: e for key, e in out.entities.items() if e.entity_type == 'gene_family'
    }
    parent = next(k for k, e in families.items() if e.identifier == '639')
    object_types = {
        out.entities[r.object_entity_key].entity_type
        for r in out.relations
        if r.subject_entity_key == parent
    }
    # The genes of the subgroups are not repeated in the group above.
    assert object_types == {'gene_family'}
    assert {r.predicate for r in out.relations if r.subject_entity_key == parent} == {
        'has_member'
    }


def test_knowledge_origin_is_literature_curated():
    assert tuple(hgnc.config.knowledge_origin) == (KnowledgeOrigin.LITERATURE_CURATED,)
