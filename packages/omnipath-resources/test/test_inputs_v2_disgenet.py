"""DisGeNET curated gene-disease associations: parsers, mappers and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import disgenet
from pypath.inputs_v2.parsers.disgenet import iter_associations, iter_diseases

HEADER = (
    'geneId\tgeneSymbol\tDSI\tDPI\tdiseaseId\tdiseaseName\tdiseaseType\tdiseaseClass\tdiseaseSemanticType\tscore'
    '\tEI\tYearInitial\tYearFinal\tNofPmids\tNofSnps\tsource'
)
# Lines of curated_gene_disease_associations.tsv; the numbers are padded with blanks.
LINES = [
    '        1\tA1BG\t0.857\t0.172\tC0019209\tHepatomegaly\tphenotype\tC06;C23\tFinding\t0.30\t\t2017\t2017\t   1\t   0\tCTD_human',
    '        1\tA1BG\t0.857\t0.172\tC0036341\tSchizophrenia\tdisease\tF03\tMental or Behavioral Dysfunction\t0.30\t\t2015\t2015\t   1\t   0\tCTD_human',
    '        2\tA2M\t0.564\t0.724\tC0002395\tAlzheimer\'s Disease\tdisease\tC10;F03\tDisease or Syndrome\t0.40\t0.848484848484849\t1998\t2016\t   3\t   0\tCTD_human',
    '        2\tA2M\t0.564\t0.724\tC0009375\tColonic Neoplasms\tgroup\tC04;C06\tNeoplastic Process\t0.30\t\t2004\t2004\t   1\t   0\tCTD_human',
    '       16\tAARS\t0.667\t0.345\tC2750090\tCharcot-Marie-Tooth Disease, Axonal, Type 2n\tdisease\tC10;C16\tDisease or Syndrome\t0.70\t\t2010\t2016\t   3\t   2\tCTD_human;ORPHANET;UNIPROT',
    '        2\tA2M\t0.564\t0.724\tC0011570\tMental Depression\tdisease\tF01\tMental or Behavioral Dysfunction\t0.31\t1\t1987\t2000\t   2\t   0\tPSYGENET',
]


def opener(lines=LINES, header=HEADER):
    text = '\n'.join([header, *lines]) + '\n'
    return SimpleNamespace(result=io.StringIO(text))


def rows(**kwargs):
    return list(iter_associations(opener(**kwargs)))


def extract(dataset, records):
    out = SilverExtractor('disgenet', dataset)
    mapper = getattr(disgenet.resource, dataset).mapper
    for index, record in enumerate(records):
        out.process_record(mapper(record), record, f'fixture:{index}', index)
    return out


def values(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term and a.get('scope') == 'relation']


def quantity(relation, term):
    return [a['quantity']['has_numeric_value'] for a in relation.annotations if a['term'] == term][0]


def test_parser_trims_the_padded_cells():
    records = rows()
    assert len(records) == 6
    assert records[0] == {
        'entrez_id': '1',
        'symbol': 'A1BG',
        'disease_id': 'C0019209',
        'disease_name': 'Hepatomegaly',
        'disease_type': 'phenotype',
        'dsi': '0.857',
        'dpi': '0.172',
        'score': '0.30',
        'ei': None,
        'year_initial': '2017',
        'year_final': '2017',
        'n_pmids': '1',
        'n_snps': '0',
        'sources': ['CTD_human'],
    }
    assert records[4]['sources'] == ['CTD_human', 'ORPHANET', 'UNIPROT']
    assert records[2]['ei'] == '0.848484848484849'


def test_parser_rejects_a_changed_header_and_an_empty_table():
    with pytest.raises(ValueError, match='lacks expected columns'):
        rows(header=HEADER.replace('NofPmids', 'N'))
    with pytest.raises(ValueError, match='has no rows'):
        rows(lines=[])
    with pytest.raises(ValueError, match='two names or types'):
        list(iter_diseases(opener([LINES[1], LINES[1].replace('Schizophrenia', 'Schizo')])))


def test_diseases_are_umls_terms_with_their_names():
    diseases = list(iter_diseases(opener(LINES + LINES[:2])))
    assert [d['disease_id'] for d in diseases] == [
        'C0019209', 'C0036341', 'C0002395', 'C0009375', 'C2750090', 'C0011570',
    ]
    term = disgenet.resource.terms.mapper(diseases[2])
    assert term.type == 'ontology_class'
    assert term.identifiers == [('umls', 'C0002395'), (Namespace.NAME, "Alzheimer's Disease")]
    assert term.identifiers[0][0] is Namespace.UMLS
    out = extract('terms', diseases)
    entity = [e for e in out.entities.values() if e.identifier == 'C0009375'][0]
    assert entity.namespace == 'umls' and entity.label == 'Colonic Neoplasms'
    assert [a['value'] for a in entity.annotations if a['term'] == 'disgenet:disease_type'] == ['group']


def test_gene_is_associated_with_the_disease_with_the_attributes_of_the_association():
    out = extract('associations', rows())
    assert [r.predicate for r in out.relations] == ['associated_with'] * 6
    relation = out.relations[4]
    assert values(relation, 'supporting_data_source') == ['CTD_human', 'ORPHANET', 'UNIPROT']
    assert quantity(relation, 'disgenet:score') == pytest.approx(0.70)
    assert quantity(relation, 'disgenet:dsi') == pytest.approx(0.667)
    assert quantity(relation, 'disgenet:dpi') == pytest.approx(0.345)
    assert quantity(relation, 'disgenet:publication_count') == 3
    assert quantity(relation, 'disgenet:variant_count') == 2
    assert values(relation, 'disgenet:year_initial') == ['2010']
    assert values(relation, 'disgenet:year_final') == ['2016']
    # The evidence index is absent where the file has none.
    assert [a for a in out.relations[0].annotations if a['term'] == 'disgenet:evidence_index'] == []
    assert quantity(out.relations[2], 'disgenet:evidence_index') == pytest.approx(0.848484848484849)
    genes = [e for e in out.entities.values() if e.entity_type == 'gene']
    assert len(genes) == 3
    gene = [e for e in genes if any(i['id'] == 'A2M' for i in e.identifiers)][0]
    assert {(i['ns'], i['id']) for i in gene.identifiers} == {(Namespace.ENTREZ, '2'), (Namespace.GENESYMBOL, 'A2M')}
    assert gene.taxon == '9606'


def test_object_is_the_entity_of_the_terms_dataset():
    term = disgenet.resource.terms.mapper(list(iter_diseases(opener()))[0])
    statement = disgenet.resource.associations.mapper(rows()[0])
    assert statement.object.type == term.type
    assert statement.object.identifiers[0] == term.identifiers[0] == ('umls', 'C0019209')


def test_the_vocabulary_dataset_is_called_terms():
    assert set(disgenet.resource.datasets()) == {'terms', 'associations'}


def test_unknown_values_raise_and_nothing_is_inferred():
    with pytest.raises(ValueError, match='Unknown DisGeNET source database'):
        disgenet.resource.associations.mapper(rows()[0] | {'sources': ['NEWDB']})
    with pytest.raises(ValueError, match='Unknown DisGeNET disease type'):
        disgenet.resource.associations.mapper(rows()[0] | {'disease_type': 'syndrome'})
    # The score is the source's ranking, not a confidence, and no causal statement is made.
    out = extract('associations', rows())
    assert not any(a['term'] == 'has_confidence_score' for r in out.relations for a in r.annotations)
    assert {r.predicate for r in out.relations} == {'associated_with'}
