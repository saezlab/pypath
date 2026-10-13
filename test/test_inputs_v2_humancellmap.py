"""Human Cell Map: parser, mapper and extractor on lines of preys-v1.txt."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.humancellmap import LOCALIZATION_NAMESPACE, config, resource
from pypath.inputs_v2.parsers.humancellmap import iter_localizations, iter_terms

HEADER = (
    'symbol\tNMF rank\tMMF localization\tSAFE domain\tSAFE localization\tEntrez\tBioGRID\t'
    'HGNC\tRefSeq\tUniProt'
)
# Lines of preys-v1.txt. Two NMF labels, and SAFE "no prediction".
AAK1 = 'AAK1\t5\tactin cytoskeleton, cytosol\t1\tno prediction\t22848\t116520\tHGNC:19679\tNM_014911\tQ2M2I8'
# NMF rank 19 is "miscellaneous"; SAFE has a domain and a label.
AAR2 = 'AAR2\t19\tmiscellaneous\t5\tspliceosomal complex I\t25980\t117464\tHGNC:15886\tNM_015511\tQ9Y312'
# Neither method has anything to say for the SAFE column: "-" for domain and label.
ABCE1 = 'ABCE1\t19\tmiscellaneous\t-\t-\t6059\t111986\tHGNC:69\tNM_002940\tP61221'
# No UniProt accession (the last column is empty).
ATP5MF = (
    'ATP5MF-PTCD1\t4\tmitochondrial matrix\t21\tmitochondrial inner membrane, '
    'mitochondrial intermembrane space, mitochondrial matrix\t100526740\t1529301\t'
    'HGNC:38844\tNM_001198879.1\t'
)
ABCD1 = (
    'ABCD1\t13\tmitochondrial outer membrane, peroxisome\t2\tmitochondrial outer membrane, '
    'peroxisome\t215\t106717\tHGNC:61\tNM_000033\tP33897'
)


def opener(*lines, header=HEADER):
    return SimpleNamespace(result=io.StringIO('\n'.join([header, *lines]) + '\n'))


def process(dataset, row):
    out = SilverExtractor('humancellmap', dataset)
    out.process_record(resource.datasets()[dataset].mapper(row), row, 'fixture:1', 1)
    return out


def relation_of(out):
    (relation,) = out.relations
    return (
        out.entities[relation.subject_entity_key],
        relation,
        out.entities[relation.object_entity_key],
    )


def attributes(relation):
    return {a['term']: a['value'] for a in relation.annotations if a['scope'] == 'relation'}


def test_license_is_the_one_the_site_states_and_the_file_is_the_pinned_version():
    assert config.license.name == 'CC_BY_4_0'
    assert resource.localizations.download.url.endswith('/preys-v1.txt')


def test_a_prey_is_one_row_per_label_of_each_method():
    rows = list(iter_localizations(opener(AAK1, ABCD1)))
    assert [(r['uniprot'], r['method'], r['label'], r['group']) for r in rows] == [
        ('Q2M2I8', 'NMF', 'actin cytoskeleton', '5'),
        ('Q2M2I8', 'NMF', 'cytosol', '5'),
        ('P33897', 'NMF', 'mitochondrial outer membrane', '13'),
        ('P33897', 'NMF', 'peroxisome', '13'),
        ('P33897', 'SAFE', 'mitochondrial outer membrane', '2'),
        ('P33897', 'SAFE', 'peroxisome', '2'),
    ]


def test_no_prediction_dashes_and_miscellaneous_are_no_localization():
    rows = list(iter_localizations(opener(AAK1, AAR2, ABCE1)))
    labels = {(r['symbol'], r['method'], r['label']) for r in rows}
    assert ('AAK1', 'SAFE', 'no prediction') not in labels
    assert not any(label == 'miscellaneous' for _, _, label in labels)
    assert not any(symbol == 'ABCE1' for symbol, _, _ in labels)
    assert ('AAR2', 'SAFE', 'spliceosomal complex I') in labels
    assert all(label not in ('-', 'no prediction') for _, _, label in labels)


def test_statement_is_located_in_the_label_of_the_map_with_method_and_group():
    rows = list(iter_localizations(opener(AAR2)))
    protein, relation, term = relation_of(process('localizations', rows[0]))
    assert relation.predicate == 'located_in'
    assert (protein.entity_type, protein.namespace, protein.identifier, protein.taxon) == (
        'protein', 'uniprot', 'Q9Y312', '9606'
    )
    assert {(i['ns'], i['id']) for i in protein.identifiers} == {
        (Namespace.UNIPROT, 'Q9Y312'), (Namespace.ENTREZ, '25980'),
        (Namespace.HGNC, 'HGNC:15886'), (Namespace.GENESYMBOL, 'AAR2'),
    }
    assert (term.entity_type, term.namespace, term.identifier) == (
        'ontology_class', LOCALIZATION_NAMESPACE, 'spliceosomal complex I'
    )
    assert attributes(relation) == {
        'humancellmap:method': 'SAFE', 'humancellmap:safe_domain': '5',
    }


def test_nmf_statement_has_the_rank():
    rows = list(iter_localizations(opener(AAK1)))
    _, relation, _ = relation_of(process('localizations', rows[0]))
    assert attributes(relation) == {'humancellmap:method': 'NMF', 'humancellmap:nmf_rank': '5'}


def test_a_prey_without_a_uniprot_accession_is_a_gene():
    rows = list(iter_localizations(opener(ATP5MF)))
    assert len(rows) == 4
    assert all(r['uniprot'] is None for r in rows)
    subject, _, _ = relation_of(process('localizations', rows[0]))
    assert subject.entity_type == 'gene'
    assert {(i['ns'], i['id']) for i in subject.identifiers} == {
        (Namespace.ENTREZ, '100526740'), (Namespace.HGNC, 'HGNC:38844'),
        (Namespace.GENESYMBOL, 'ATP5MF-PTCD1'),
    }


def test_terms_are_the_labels_each_once_with_a_name():
    terms = list(iter_terms(opener(AAK1, AAR2, ABCD1)))
    assert {t['label']: t['methods'] for t in terms} == {
        'actin cytoskeleton': ['NMF'], 'cytosol': ['NMF'],
        'spliceosomal complex I': ['SAFE'],
        'mitochondrial outer membrane': ['NMF', 'SAFE'], 'peroxisome': ['NMF', 'SAFE'],
    }
    out = SilverExtractor('humancellmap', 'terms')
    term = next(t for t in terms if t['label'] == 'cytosol')
    out.process_record(resource.terms.mapper(term), term, 'fixture:1', 1)
    (entity,) = out.entities.values()
    assert {(i['ns'], i['id']) for i in entity.identifiers} == {
        (LOCALIZATION_NAMESPACE, 'cytosol'), (Namespace.NAME, 'cytosol'),
    }


def test_objects_have_the_type_namespace_and_identifier_of_their_terms():
    terms = list(iter_terms(opener(AAK1, AAR2, ABCD1, ATP5MF)))
    defined = {}
    for term in terms:
        out = SilverExtractor('humancellmap', 'terms')
        out.process_record(resource.terms.mapper(term), term, 'fixture:1', 1)
        (entity,) = out.entities.values()
        defined[entity.identifier] = (entity.entity_type, entity.namespace, entity.identifier)
    for row in iter_localizations(opener(AAK1, AAR2, ABCD1, ATP5MF)):
        _, _, term = relation_of(process('localizations', row))
        assert defined[term.identifier] == (term.entity_type, term.namespace, term.identifier)


def test_no_ontology_accession_is_given_to_a_label_and_the_proximity_interactions_are_not_read():
    for term in iter_terms(opener(AAK1, AAR2)):
        entity = resource.terms.mapper(term)
        assert {i.type for i in entity.identifiers} == {LOCALIZATION_NAMESPACE, Namespace.NAME}
    # Only the localization and term datasets exist: the SAINT report is not read.
    assert set(resource.datasets()) == {'terms', 'localizations'}


def test_unknown_method_raises_in_the_mapper():
    (row, *_) = iter_localizations(opener(AAK1))
    with pytest.raises(ValueError, match='Unknown Human Cell Map method'):
        resource.localizations.mapper({**row, 'method': 'PCA'})


def test_parser_raises_on_a_changed_header_and_an_empty_table():
    with pytest.raises(ValueError, match='lacks expected columns'):
        list(iter_localizations(opener(AAK1, header=HEADER.replace('UniProt', 'Uniprot'))))
    with pytest.raises(ValueError):
        list(iter_localizations(opener()))
