"""CancerSEA functional states and their gene sets: parsers, mappers and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import cancersea
from pypath.inputs_v2.parsers.cancersea import iter_gene_sets, iter_states

# Rows of the states table of the download page. The page gives 73 and 137
# genes; the gene sets here are cut to three, so the counts are changed.
ROW = (
    '<tr><td>{state}</td><td>{description}</td><td>{count}</td>'
    '<td><a id="sigId-1" class="sourceDetail">{sources}</a></td>'
    '<td><a href="download/signature/{file}.txt"><i class="step fi-download size-16"></i></a></td></tr>'
)
STATES = [
    dict(
        state='Angiogenesis', file='Angiogenesis', count=3, sources='15614775,GO,MsigDB',
        description='Angiogenesis ensures that cancer cells receive continuous supplies of oxygen and other nutrients.',
    ),
    dict(
        state='Cell Cycle', file='Cell_Cycle', count=3, sources='28102192,27124452,29472726,27806376,cyclebase',
        description='Cell cycle,a critical process to ensure correct cell division,lies at the heart of cancer.',
    ),
]
HEADER = (
    '<table><thead><tr><th>Signature</th><th>Description</th><th>No. genes</th><th>Source</th>'
    '<th>Download</th></tr></thead><tbody>'
)
# Another table of the page: the expression datasets.
OTHER = (
    '<table id="downloadTable"><thead><tr><th>Database ID</th><th>Dataset Name</th></tr></thead>'
    '<tbody><tr><td>EXP0001</td><td>Chen W. Genome Biol. 2018 (Rh41)</td></tr></tbody></table>'
)
GENE_SETS = {
    'Angiogenesis': 'EnsembleID\tGeneName\nENSG00000139567\tACVRL1\nENSG00000101384\tJAG1\nENSG00000154188\tANGPT1\n',
    'Cell_Cycle': 'EnsembleID\tGeneName\nENSG00000011426\tANLN\nENSG00000143401\tANP32E\nENSG00000156802\tATAD2\n',
}


def page(states=STATES, header=HEADER):
    html = '<html><body>' + OTHER + header + ''.join(ROW.format(**s) for s in states) + '</tbody></table></body></html>'
    return SimpleNamespace(result=html)


class Fake:
    def __init__(self, text):
        self.text = text

    def open(self, **_kwargs):
        return SimpleNamespace(result=io.StringIO(self.text))


def downloads(sets=GENE_SETS):
    return {name: Fake(text) for name, text in sets.items()}


def genes(**kwargs):
    return list(iter_gene_sets(page(), downloads=kwargs.get('downloads', downloads())))


def extract(dataset, records):
    out = SilverExtractor('cancersea', dataset)
    mapper = getattr(cancersea.resource, dataset).mapper
    for index, record in enumerate(records):
        out.process_record(mapper(record), record, f'fixture:{index}', index)
    return out


def test_states_come_from_the_table_with_the_right_header():
    states = list(iter_states(page()))
    assert [s['state'] for s in states] == ['Angiogenesis', 'Cell Cycle']
    assert states[1]['description'].startswith('Cell cycle,a critical process')
    assert states[1]['sources'] == '28102192,27124452,29472726,27806376,cyclebase'
    assert states[0]['n_genes'] == 3


def test_gene_sets_are_checked_against_the_page():
    records = genes()
    assert len(records) == 6
    assert records[0] == {
        'state': 'Angiogenesis',
        'description': STATES[0]['description'],
        'ensembl_id': 'ENSG00000139567',
        'symbol': 'ACVRL1',
    }
    with pytest.raises(ValueError, match='no table of states'):
        list(iter_states(page(header=HEADER.replace('Signature', 'Name'))))
    with pytest.raises(ValueError, match='the page says 3 genes, the file has 2'):
        genes(downloads=downloads(GENE_SETS | {'Cell_Cycle': GENE_SETS['Cell_Cycle'].rsplit('\n', 2)[0] + '\n'}))
    with pytest.raises(ValueError, match='differ'):
        genes(downloads=downloads({'Angiogenesis': GENE_SETS['Angiogenesis']}))
    with pytest.raises(ValueError, match='lacks expected columns'):
        genes(downloads=downloads(GENE_SETS | {'Cell_Cycle': GENE_SETS['Cell_Cycle'].replace('GeneName', 'Gene')}))


def test_states_are_terms_with_their_description():
    state = cancersea.resource.terms.mapper(list(iter_states(page()))[1])
    assert state.type == 'ontology_class'
    assert state.identifiers == [('cancersea_state', 'Cell Cycle'), (Namespace.NAME, 'Cell Cycle')]
    out = extract('terms', list(iter_states(page())))
    entity = [e for e in out.entities.values() if e.identifier == 'Cell Cycle'][0]
    assert entity.namespace == 'cancersea_state'
    assert set(cancersea.resource.datasets()) == {'terms', 'genes'}
    assert [a['value'] for a in entity.annotations if a['term'] == 'description'][0].startswith('Cell cycle')


def test_gene_is_associated_with_the_state_it_is_a_member_of():
    out = extract('genes', genes())
    assert [r.predicate for r in out.relations] == ['associated_with'] * 6
    gene = [e for e in out.entities.values() if e.entity_type == 'gene' and any(i['id'] == 'ENSG00000101384' for i in e.identifiers)][0]
    assert {(i['ns'], i['id']) for i in gene.identifiers} == {
        (Namespace.ENSG, 'ENSG00000101384'), (Namespace.GENESYMBOL, 'JAG1'),
    }
    assert gene.taxon == '9606'
    # The statement's object is the term of the ``terms`` dataset.
    state = cancersea.resource.terms.mapper(list(iter_states(page()))[0])
    statement = cancersea.resource.genes.mapper(genes()[0])
    assert statement.object.type == state.type
    assert statement.object.identifiers[0] == state.identifiers[0]


def test_no_attribute_or_other_predicate_is_invented():
    out = extract('genes', genes())
    assert {r.predicate for r in out.relations} == {'associated_with'}
    # Membership of the gene set is the whole statement.
    assert [a['term'] for r in out.relations for a in r.annotations if a.get('scope') == 'relation'] == []
    assert len({r.relation_key for r in out.relations}) == 6


def test_the_gene_sets_are_curated_not_predicted_from_expression():
    # The help page: the 14 signatures were collected from databases and the
    # literature; the download page lists their sources.
    from omnipath_core.resource_metadata import KnowledgeOrigin

    config = cancersea.resource.config
    assert tuple(config.knowledge_origin) == (KnowledgeOrigin.LITERATURE_CURATED,)
    assert 'single-cell expression data' not in config.summary
