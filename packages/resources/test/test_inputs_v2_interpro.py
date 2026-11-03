"""InterPro: readers, mappers and extractor on lines of the release files and of UniProt's table."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import interpro
from pypath.inputs_v2.parsers import interpro as parser

ENTRY_LIST = (
    'ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\n'
    'IPR000001\tDomain\tKringle\n'
    'IPR000003\tFamily\tRetinoid X receptor/HNF4\n'
    'IPR000008\tDomain\tC2 domain\n'
    'IPR014705\tDomain\tSynaptotagmin-17, C2B domain\n'
    'IPR033884\tDomain\tCalpain C2 domain\n'
    'IPR000053\tFamily\tThymidine/pyrimidine-nucleoside phosphorylase\n'
    'IPR013466\tFamily\tThymidine phosphorylase/AMP phosphorylase\n'
    'IPR017713\tFamily\tAMP phosphorylase\n'
)
TREE = (
    'IPR000008::C2 domain::\n'
    '--IPR014705::Synaptotagmin-17, C2B domain::\n'
    '--IPR033884::Calpain C2 domain::\n'
    'IPR000053::Thymidine/pyrimidine-nucleoside phosphorylase::\n'
    '--IPR013466::Thymidine phosphorylase/AMP phosphorylase::\n'
    '----IPR017713::AMP phosphorylase::\n'
    # An entry that is also at the top of the file is listed again below its parent.
    'IPR013466::Thymidine phosphorylase/AMP phosphorylase::\n'
    '--IPR017713::AMP phosphorylase::\n'
)
INTERPRO2GO = (
    '!date: 2026/08/22 19:05:25\n'
    '!Mapping of InterPro entries to GO\n'
    'InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:DNA binding ; GO:0003677\n'
    'InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:nuclear steroid receptor activity ; GO:0003707\n'
    'InterPro:IPR013466 Thymidine phosphorylase/AMP phosphorylase > GO:pentosyltransferase activity ; GO:0016763\n'
    'InterPro:IPR017713 AMP phosphorylase > GO:pyrimidine deoxyribonucleoside metabolic process ; GO:0046125\n'
)
UNIPROT = (
    'Entry\tInterPro\n'
    'P04637\tIPR008967;IPR012346;IPR057064;IPR011615;IPR040926;IPR036674;IPR010991;IPR013872;IPR002117;\n'
    'Q9Y6K9\tIPR032419;IPR021063;IPR034735;IPR051301;\n'
    'B0L3A2\t\n'
)


def entries(entry_list=ENTRY_LIST, tree=TREE):
    download = SimpleNamespace(open=lambda **_: SimpleNamespace(result=tree))
    rows = parser.iter_entries(
        SimpleNamespace(result=entry_list), tree_download=download
    )
    return list(rows)


def extract(dataset, rows):
    out = SilverExtractor('interpro', dataset)
    mapper = interpro.resource.datasets()[dataset].mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def test_entries_have_type_name_and_parents():
    rows = {r['accession']: r for r in entries()}
    assert rows['IPR000008'] == {
        'accession': 'IPR000008', 'type': 'Domain', 'name': 'C2 domain', 'parents': [],
    }
    assert rows['IPR014705']['parents'] == ['IPR000008']
    # The chain of three levels, and an entry listed twice under the same parent.
    assert rows['IPR013466']['parents'] == ['IPR000053']
    assert rows['IPR017713']['parents'] == ['IPR013466']
    assert rows['IPR000001']['parents'] == []


def test_reader_rejects_inconsistent_files():
    with pytest.raises(ValueError, match='is named'):
        entries(tree=TREE.replace('AMP phosphorylase::', 'AMP phosphorylase 2::'))
    with pytest.raises(ValueError, match='two names'):
        entries(tree=TREE.replace('AMP phosphorylase::', 'AMP phosphorylase 2::', 1))
    with pytest.raises(ValueError, match='not in entry.list'):
        entries(tree=TREE + 'IPR999999::Unlisted::\n')
    with pytest.raises(ValueError, match='Tree depth jumps'):
        entries(tree='IPR000008::C2 domain::\n----IPR014705::Synaptotagmin-17, C2B domain::\n')
    with pytest.raises(ValueError, match='Malformed line'):
        entries(tree='C2 domain\n')
    with pytest.raises(ValueError, match='Malformed InterPro entry row'):
        entries(entry_list=ENTRY_LIST.replace('IPR000001', 'IPR1'))
    with pytest.raises(ValueError, match='lacks expected columns'):
        entries(entry_list=ENTRY_LIST.replace('ENTRY_TYPE', 'TYPE'))


def test_terms_have_names_types_and_the_hierarchy():
    out = extract('terms', entries())
    named = {e.identifier: e for e in out.entities.values() if e.label}
    assert named['IPR000008'].entity_type == 'ontology_class'
    assert named['IPR000008'].namespace == 'interpro'
    assert interpro.NAMESPACE is Namespace.INTERPRO
    assert named['IPR000008'].label == 'C2 domain'
    types = {
        e.identifier: [a['value'] for a in e.annotations if a['term'] == 'interpro:type']
        for e in named.values()
    }
    assert types['IPR000003'] == ['Family'] and types['IPR000001'] == ['Domain']
    by_key = {k: e.identifier for k, e in out.entities.items()}
    edges = {(by_key[r.subject_entity_key], r.predicate, by_key[r.object_entity_key]) for r in out.relations}
    assert edges == {
        ('IPR014705', 'subclass_of', 'IPR000008'),
        ('IPR033884', 'subclass_of', 'IPR000008'),
        ('IPR013466', 'subclass_of', 'IPR000053'),
        ('IPR017713', 'subclass_of', 'IPR013466'),
    }
    # The top-level entry has no parent: the hierarchy is not inferred from names.
    assert not any(s in {'IPR000008', 'IPR000053', 'IPR000001'} for s, _, _ in edges)


def test_go_terms_are_read_from_interpro2go():
    rows = list(parser.iter_go_terms(SimpleNamespace(result=INTERPRO2GO)))
    assert rows[0] == {'accession': 'IPR000003', 'go_terms': ['GO:0003677', 'GO:0003707']}
    out = extract('go_terms', rows)
    pairs = {
        (out.entities[r.subject_entity_key].identifier, out.entities[r.object_entity_key].identifier)
        for r in out.relations
    }
    assert pairs == {
        ('IPR000003', 'GO:0003677'), ('IPR000003', 'GO:0003707'),
        ('IPR013466', 'GO:0016763'), ('IPR017713', 'GO:0046125'),
    }
    assert {out.entities[r.object_entity_key].namespace for r in out.relations} == {Namespace.GO}
    with pytest.raises(ValueError, match='Malformed line'):
        list(parser.iter_go_terms(SimpleNamespace(result='InterPro:IPR000003 x > GO:y\n')))


def test_proteins_are_associated_with_the_entries_that_match_them():
    rows = list(parser.iter_protein_entries(SimpleNamespace(result=UNIPROT)))
    assert rows[1] == {'accession': 'Q9Y6K9', 'entries': ['IPR032419', 'IPR021063', 'IPR034735', 'IPR051301']}
    assert rows[2] == {'accession': 'B0L3A2', 'entries': []}
    out = extract('protein_entries', rows)
    proteins = [e for e in out.entities.values() if e.entity_type == 'protein']
    # A protein with no InterPro entry gives no statement and no entity.
    assert sorted(e.identifier for e in proteins) == ['P04637', 'Q9Y6K9']
    assert {e.taxon for e in proteins} == {'9606'}
    assert len(out.relations) == 9 + 4
    assert {r.predicate for r in out.relations} == {'associated_with'}
    with pytest.raises(ValueError, match='Malformed UniProt InterPro row'):
        list(parser.iter_protein_entries(SimpleNamespace(result='Entry\tInterPro\nP1\tPF00001;\n')))


def test_the_object_of_a_protein_statement_is_the_term_of_the_entry():
    terms = extract('terms', entries())
    term_keys = {(e.entity_type, e.namespace, e.identifier) for e in terms.entities.values()}
    proteins = extract('protein_entries', [{'accession': 'P04637', 'entries': ['IPR000003', 'IPR000001']}])
    objects = {
        (e.entity_type, e.namespace, e.identifier)
        for e in (proteins.entities[r.object_entity_key] for r in proteins.relations)
    }
    assert objects <= term_keys
    assert len(objects) == 2


def test_the_description_says_where_the_protein_statements_come_from():
    description = interpro.config.description
    assert "cross-references of UniProt" in description
    assert 'without the positions' in description


def test_the_vocabulary_dataset_is_called_terms():
    assert list(interpro.resource.datasets()) == ['terms', 'go_terms', 'protein_entries']
