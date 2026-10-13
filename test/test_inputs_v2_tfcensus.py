"""TF census: reader, mappers and extractor on lines of Supplementary Table 2 (Vaquerizas 2009)."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import tfcensus
from pypath.inputs_v2.parsers import tfcensus as parser

PREAMBLE = (
    'A census of human transcription factors: function, expression and evolution\r\n'
    'Juan M Vaquerizas, Sarah K Kummerfeld, Sarah A Teichmann, Nicholas M Luscombe\r\n'
    '\r\n'
    'Supplementary Table 2\r\n'
    "List of TF-encoding loci classified as 'a', 'b', 'c', 'x', or selected from other "
    "databases (classified as 'other'; see Supplementary Material).\r\n"
    '\r\n'
    '--Table start--\r\n'
)
HEADER = (
    'Class\tEnsembl ID\tIPI ID\tInterpro DBD\tInterpro DNA-binding family\t'
    'HGNC symbol\tTissue-specificity\r\n'
)
# Lines of the table: a gene with a family and no domain, one with domains and a
# tissue, one with three tissues, a gene without a symbol, a class x gene with
# no tissue column at all, and classes b, c and other.
LINES = [
    'a\tENSG00000001167\tIPI00218470;IPI00333568\t\tIPR001289\tNFYA\t',
    'a\tENSG00000005073\tIPI00010754;IPI00807582\tIPR001356;IPR009057;IPR012287\t\tHOXA11\tuterus',
    'a\tENSG00000006468\tIPI00221065;IPI00306430\tIPR000418;IPR011991\t\tETV1\tadrenal.cortex;whole.brain',
    'a\tENSG00000108452\tIPI00012394\tIPR007087;IPR015880\t\t\t',
    'b\tENSG00000005801\tIPI00442264\tIPR007086;IPR007087;IPR015880\t\tZNF195\tgeneral',
    'c\tENSG00000007545\tIPI00658167\tIPR001005;IPR014778\t\tCRAMP1L\t',
    'other\tENSG00000029363\tIPI00006079;IPI00413672\t\t\tBCLAF1\tsmooth.muscle',
    'x\tENSG00000003756\tIPI00005036\tIPR007087;IPR015880\t\tRBM5',
]


def census(lines=LINES, preamble=PREAMBLE, header=HEADER):
    text = preamble + header + ''.join(line + '\r\n' for line in lines)
    return SimpleNamespace(result=text)


def extract(dataset, rows):
    out = SilverExtractor('tfcensus', dataset)
    mapper = tfcensus.resource.datasets()[dataset].mapper
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def test_reader_skips_the_preamble_and_decodes_the_cells():
    rows = list(parser.iter_census(census()))
    assert len(rows) == 8
    assert rows[0] == {
        'census_class': 'a',
        'ensembl': 'ENSG00000001167',
        'ipi': ['IPI00218470', 'IPI00333568'],
        'interpro_dbd': [],
        'interpro_family': ['IPR001289'],
        'symbol': 'NFYA',
        'tissues': [],
    }
    assert rows[2]['tissues'] == ['adrenal.cortex', 'whole.brain']
    assert rows[3]['symbol'] is None
    # A line without the last column has no tissues.
    assert rows[7]['census_class'] == 'x' and rows[7]['tissues'] == []


def test_reader_rejects_a_changed_file():
    with pytest.raises(ValueError, match='Table start'):
        list(parser.iter_census(census(preamble='Title\r\n')))
    with pytest.raises(ValueError, match='lacks expected columns'):
        list(parser.iter_census(census(header=HEADER.replace('Class', 'Category'))))
    with pytest.raises(ValueError, match='no rows'):
        list(parser.iter_census(census(lines=[])))
    with pytest.raises(ValueError, match='Malformed TF census row'):
        list(parser.iter_census(census(lines=['a\t\tIPI00218470\t\t\tNFYA\t'])))


def test_genes_are_associated_with_the_list_with_the_class_as_confidence():
    rows = list(parser.iter_census(census()))
    record = tfcensus.map_transcription_factor(rows[1])
    assert record.predicate.curie == 'biolink:associated_with'
    assert record.subject.type.class_name == 'gene'
    assert [(i.type, i.value) for i in record.subject.identifiers] == [
        (Namespace.ENSG, 'ENSG00000005073'),
        (Namespace.GENESYMBOL, 'HOXA11'),
    ]

    out = extract('transcription_factors', rows)
    # Class x is the census saying a gene is unlikely to be a TF: no statement.
    assert len(out.relations) == 7
    genes = {e.identifier: e for e in out.entities.values() if e.entity_type == 'gene'}
    assert 'ENSG00000003756' not in {i['id'] for e in genes.values() for i in e.identifiers}
    levels = {}
    for r in out.relations:
        subject = out.entities[r.subject_entity_key]
        ids = {i['id'] for i in subject.identifiers}
        confidence = [a['value'] for a in r.annotations if a['term'] == 'has_confidence_level']
        tissues = [a['value'] for a in r.annotations if a['term'] == 'tfcensus:tissue']
        levels[next(i for i in ids if i.startswith('ENSG'))] = (confidence, tissues)
    assert levels['ENSG00000006468'] == (['a'], ['adrenal.cortex', 'whole.brain'])
    assert levels['ENSG00000001167'] == (['a'], [])
    assert levels['ENSG00000007545'] == (['c'], [])
    assert levels['ENSG00000029363'] == (['other'], ['smooth.muscle'])
    # The gene without a symbol is identified by its Ensembl gene ID.
    assert any(
        e.namespace == Namespace.ENSG and e.identifier == 'ENSG00000108452'
        for e in out.entities.values()
    )
    assert {e.taxon for e in genes.values()} == {'9606'}


def test_the_object_of_a_statement_is_the_term_of_the_list():
    rows = list(parser.iter_census(census()))
    terms = extract('terms', parser.iter_census_list(census()))
    (term,) = terms.entities.values()
    out = extract('transcription_factors', rows)
    objects = {
        (out.entities[r.object_entity_key].entity_type,
         out.entities[r.object_entity_key].namespace,
         out.entities[r.object_entity_key].identifier)
        for r in out.relations
    }
    assert objects == {(term.entity_type, term.namespace, term.identifier)}
    assert term.entity_type == 'ontology_class'
    assert (term.namespace, term.identifier) == ('tfcensus_list', 'human_transcription_factors')
    assert term.label == 'Human transcription factors (TF census)'
    description = next(a['value'] for a in term.annotations if a['term'] == 'description')
    # The classes as the Supplementary Material of the paper defines them.
    assert 'class a: probable TFs with experimental evidence of regulatory function' in description
    assert 'class b: probable TFs with a protein domain arrangement equivalent' in description
    assert 'class c: possible TFs with non-promiscuous InterPro DNA-binding domains' in description
    assert 'undefined DNA-binding domains' in description
    assert 'unlikely TFs' in description
    # Class b is not about orthologues (a second-hand reading of the table).
    assert 'orthologous' not in description
    assert 'class c is a possible TF' in description


def test_interpro_entries_are_associated_by_column_and_class_x_is_left_out():
    rows = list(parser.iter_census_interpro(census()))
    assert len([r for r in rows if r['census_class'] == 'x']) == 2
    record = tfcensus.map_interpro_entry(rows[0])
    assert (record.object.identifier_type, record.object.identifier) == ('interpro', 'IPR001289')
    assert record.object.identifier_type is Namespace.INTERPRO
    out = extract('interpro_entries', rows)
    pairs = [
        (
            next(i['id'] for i in out.entities[r.subject_entity_key].identifiers if i['id'].startswith('ENSG')),
            out.entities[r.object_entity_key].identifier,
            [a['value'] for a in r.annotations if a['term'] == 'tfcensus:interpro_column'],
        )
        for r in out.relations
    ]
    assert ('ENSG00000001167', 'IPR001289', ['DNA-binding family']) in pairs
    assert ('ENSG00000005073', 'IPR012287', ['DNA-binding domain']) in pairs
    assert all(subject != 'ENSG00000003756' for subject, _, _ in pairs)
    assert len(pairs) == 1 + 3 + 2 + 2 + 3 + 2
    assert {out.entities[r.object_entity_key].namespace for r in out.relations} == {'interpro'}
    bad = dict(rows[0], interpro='IPR1')
    with pytest.raises(ValueError, match='Malformed InterPro ID'):
        tfcensus.map_interpro_entry(bad)


def test_an_unknown_class_raises():
    row = dict(list(parser.iter_census(census()))[0], census_class='d')
    with pytest.raises(ValueError, match='Unknown TF census class'):
        tfcensus.map_transcription_factor(row)


def test_the_vocabulary_dataset_is_called_terms():
    assert list(tfcensus.resource.datasets()) == ['terms', 'transcription_factors', 'interpro_entries']
