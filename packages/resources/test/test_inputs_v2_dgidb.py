"""DGIdb gene categories and drug-gene interactions: parsers, mappers and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import dgidb
from pypath.inputs_v2.parsers.dgidb import iter_categories, iter_category_terms, iter_interactions

PREAMBLE = '# Data version: 2026-09\n# DGIdb version: v.5.0.13\n'
CATEGORY_HEADER = 'gene_claim_name\tgene_category_name\tsource_db_name\tsource_db_version'
GENE_HEADER = 'gene_claim_name\tnomenclature\tconcept_id\tgene_name\tsource_db_name\tsource_db_version'
INTERACTION_HEADER = (
    'gene_claim_name\tgene_concept_id\tgene_name\tdrug_claim_name\tdrug_concept_id\tdrug_name\t'
    'drug_is_approved\tdrug_is_immunotherapy\tdrug_is_antineoplastic\tinteraction_source_db_name\t'
    'interaction_source_db_version\tinteraction_types\tinteraction_score\tdrug_specificity_score\t'
    'gene_specificity_score\tevidence_score'
)

# Lines of categories.tsv and genes.tsv (release 2026-09).
CATEGORIES = [
    'UNIPROT:O00748\tDRUGGABLE GENOME\tHopkinsGroom\t11-Sep-2012',
    'LXRA\tNUCLEAR HORMONE RECEPTOR\tBaderLab\tFebruary 2014',
    'RORB\tKINASE\tPharos\t01-September-2026',
    # A claim that genes.tsv maps to no gene.
    'NOSUCHCLAIM\tKINASE\tPharos\t01-September-2026',
]
GENES = [
    'UNIPROT:O00748\tUniProtKB ID\thgnc:1864\tCES2\tHopkinsGroom\t11-Sep-2012',
    'LXRA\tNCBI Gene Name\thgnc:7966\tNR1H3\tBaderLab\tFebruary 2014',
    'RORB\tGene Symbol\thgnc:10259\tRORB\tPharos\t01-September-2026',
    # The claim of another source is another claim.
    'LXRA\tGene Symbol\thgnc:1\tOTHER\tPharos\t01-September-2026',
]
# Lines of interactions.tsv; the lines of drugs DGIdb did not map end after
# the source columns.
INTERACTIONS = {
    'inhibitor': 'ABL1\thgnc:76\tABL1\tIMATINIB\trxcui:282388\tIMATINIB\ttrue\tfalse\ttrue\tTALC\t12-May-2016\tinhibitor\t0.2024766898977484\t0.04891704649042844\t0.11497735249122597\t36',
    'two types': 'IUPHAR.RECEPTOR:627\thgnc:8910\tPGR\tIUPHAR.LIGAND:2805\trxcui:6964\tMIFEPRISTONE\ttrue\tfalse\tfalse\tGuideToPharmacology\t2026.2\tinhibitor|agonist\t0.09394609658678636\t0.13102780309936188\t0.14339872052276498\t5',
    'antisense': 'RAF1\thgnc:9829\tRAF1\tISIS 5132\tncit:C2721\tC-RAF ANTISENSE OLIGONUCLEOTIDE ISIS 5132\tfalse\tfalse\ttrue\tTALC\t12-May-2016\tantisense oligonucleotide\t0.8214514831377951\t3.668778486782133\t0.22390326537765057\t1',
    'binder': 'IFNA1\thgnc:5417\tIFNA1\tFARALIMOMAB\tchembl:CHEMBL2108407\tFARALIMOMAB\tfalse\tfalse\tfalse\tChEMBL\t37\tbinder\t2.128306115402469\t1.8343892433910665\t1.1602260115023713\t1',
    'positive modulator': 'HBB\thgnc:4827\tHBB\tVOXELOTOR\trxcui:2265678\tVOXELOTOR\ttrue\tfalse\tfalse\tChEMBL\t37\tpositive modulator\t0.7094353718008232\t0.6114630811303555\t0.3867420038341237\t3',
    'antibody': 'AOC3\thgnc:550\tAOC3\tIUPHAR.LIGAND:9097\tiuphar.ligand:9097\tTIMOLUMAB\tfalse\tfalse\tfalse\tGuideToPharmacology\t2026.2\tantibody\t5.202526059872703\t3.668778486782133\t1.4180540140584537\t1',
    # No type (the cell is two quote characters) and a concept without a name.
    'untyped placeholder': 'STS\thgnc:11425\tSTS\t17ALPHA-3\'-BROMOBENZYL ESTRADIOL\tchembl:CHEMBL1627429\tCHEMBL:CHEMBL1627429\tfalse\tfalse\tfalse\tDTC\t2020-09-02\t""\t6.243031271847244\t3.668778486782133\t0.8508324084350722\t2',
    # A drug that DGIdb did not map: the line ends after the source columns.
    'unmapped drug': 'RIBONUCLEOSIDE-DIPHOSPHATE REDUCTASE M2\thgnc:10452\tRRM2\tGemcitabine prodrug\t\t\t\t\t\tTTD\t2020.06.01\t""',
    # A gene that DGIdb did not map.
    'unmapped gene': '8231511\t\t\tPIPERONYL BUTOXIDE\tchembl:CHEMBL1201131\tPIPERONYL BUTOXIDE\ttrue\tfalse\tfalse\tChEMBL\t37\tinhibitor',
}


class Fake:
    """A download of the lines given."""

    def __init__(self, header, lines):
        self.text = PREAMBLE + '\n'.join([header, *lines]) + '\n'

    def open(self, **_kwargs):
        return SimpleNamespace(result=io.StringIO(self.text))


def parse_categories(lines=CATEGORIES, genes=GENES):
    opener = Fake(CATEGORY_HEADER, lines).open()
    return list(iter_categories(opener, genes_download=Fake(GENE_HEADER, genes)))


def parse_interactions(names=None, header=INTERACTION_HEADER):
    lines = [INTERACTIONS[name] for name in names or INTERACTIONS]
    return list(iter_interactions(Fake(header, lines).open()))


def interaction(name):
    return parse_interactions([name])


def extract(dataset, records):
    out = SilverExtractor('dgidb', dataset)
    mapper = getattr(dgidb.resource, dataset).mapper
    for index, record in enumerate(records):
        mapped = mapper(record)
        if mapped is not None:
            out.process_record(mapped, record, f'fixture:{index}', index)
    return out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term and a.get('scope') == 'relation']


def test_category_claims_are_joined_to_the_gene_concepts():
    records = parse_categories()
    assert len(records) == 4
    # The claim is joined on its name and its source database.
    assert (records[0]['gene_concept_id'], records[0]['gene_name']) == ('hgnc:1864', 'CES2')
    assert (records[1]['gene_concept_id'], records[1]['gene_name']) == ('hgnc:7966', 'NR1H3')
    assert (records[2]['gene_concept_id'], records[2]['category'], records[2]['source']) == (
        'hgnc:10259', 'KINASE', 'Pharos',
    )
    assert records[3]['gene_concept_id'] is None


def test_parsers_reject_a_changed_header_and_an_empty_table():
    with pytest.raises(ValueError, match='lacks expected columns'):
        list(
            iter_categories(
                Fake(CATEGORY_HEADER.replace('gene_category_name', 'category'), CATEGORIES).open(),
                genes_download=Fake(GENE_HEADER, GENES),
            )
        )
    with pytest.raises(ValueError, match='has no rows'):
        parse_categories(lines=[])
    with pytest.raises(ValueError, match='lacks expected columns'):
        parse_interactions(header=INTERACTION_HEADER.replace('interaction_types', 'types'))
    with pytest.raises(ValueError, match='two concepts'):
        parse_categories(genes=[*GENES, 'RORB\tGene Symbol\thgnc:2\tOTHER\tPharos\tv'])


def test_category_terms_are_listed_once_and_match_the_statements_object():
    terms = list(iter_category_terms(Fake(CATEGORY_HEADER, CATEGORIES + CATEGORIES).open()))
    assert [t['category'] for t in terms] == ['DRUGGABLE GENOME', 'NUCLEAR HORMONE RECEPTOR', 'KINASE']

    term = dgidb.resource.terms.mapper(terms[1])
    assert term.type == 'ontology_class'
    assert term.identifiers == [
        ('dgidb_category', 'NUCLEAR HORMONE RECEPTOR'),
        (Namespace.NAME, 'NUCLEAR HORMONE RECEPTOR'),
    ]
    statement = dgidb.resource.categories.mapper(parse_categories()[1])
    assert statement.object.type == term.type
    assert statement.object.identifiers[0] == term.identifiers[0]


def test_category_statement_is_an_association_with_its_source():
    out = extract('categories', parse_categories())
    # The claim that maps to no gene gives nothing.
    assert len(out.relations) == 3
    relation = out.relations[1]
    assert relation.predicate == 'associated_with'
    assert annotations(relation, 'supporting_data_source') == ['BaderLab']
    genes = {e.identifier: e for e in out.entities.values() if e.entity_type == 'gene'}
    assert set(genes) == {'1864', '7966', '10259'}
    assert {(i['ns'], i['id']) for i in genes['7966'].identifiers} == {
        (Namespace.HGNC, '7966'), (Namespace.GENESYMBOL, 'NR1H3'),
    }
    assert genes['7966'].taxon == '9606'
    category = [e for e in out.entities.values() if e.entity_type == 'ontology_class']
    assert {(e.namespace, e.identifier) for e in category} == {
        ('dgidb_category', 'DRUGGABLE GENOME'),
        ('dgidb_category', 'NUCLEAR HORMONE RECEPTOR'),
        ('dgidb_category', 'KINASE'),
    }
    # The category is not identified by its name alone.
    assert all(Namespace.NAME not in {i['ns'] for i in e.identifiers[:1]} for e in category)


def test_unknown_category_source_raises():
    bad = parse_categories()[0] | {'source': 'NewDatabase'}
    with pytest.raises(ValueError, match='Unknown DGIdb source'):
        dgidb.resource.categories.mapper(bad)


def test_interaction_types_are_split_and_empty_cells_are_none():
    records = parse_interactions()
    # The claim of two types is two records, so one record more than lines.
    assert len(records) == len(INTERACTIONS) + 1
    assert [r['interaction_type'] for r in records if r['drug_claim_name'] == 'IUPHAR.LIGAND:2805'] == [
        'inhibitor', 'agonist',
    ]
    placeholder = interaction('untyped placeholder')[0]
    assert placeholder['interaction_type'] is None
    unmapped = interaction('unmapped drug')[0]
    assert unmapped['drug_concept_id'] is None and unmapped['interaction_score'] is None
    assert unmapped['source'] == 'TTD' and unmapped['source_version'] == '2020.06.01'


@pytest.mark.parametrize(
    ('name', 'predicate', 'direction', 'mechanism'),
    [
        ('inhibitor', 'affects', 'decreased', 'inhibition'),
        ('antisense', 'affects', 'decreased', 'antisense_oligonucleotide_inhibition'),
        ('positive modulator', 'affects', 'increased', 'positive_modulation'),
        ('binder', 'interacts_with', None, 'binding'),
        # An antibody is stated to bind its target and nothing else.
        ('antibody', 'interacts_with', None, None),
        ('untyped placeholder', 'interacts_with', None, None),
        ('unmapped drug', 'interacts_with', None, None),
    ],
)
def test_interaction_type_crosswalk(name, predicate, direction, mechanism):
    out = extract('interactions', interaction(name))
    relation = out.relations[0]
    assert relation.predicate == predicate
    assert annotations(relation, 'object_direction_qualifier') == ([direction] if direction else [])
    assert annotations(relation, 'causal_mechanism_qualifier') == ([mechanism] if mechanism else [])


def test_both_types_of_a_claim_are_statements():
    out = extract('interactions', interaction('two types'))
    assert [r.predicate for r in out.relations] == ['affects', 'affects']
    assert sorted(annotations(r, 'object_direction_qualifier')[0] for r in out.relations) == [
        'decreased', 'increased',
    ]
    assert len({r.relation_key for r in out.relations}) == 2


def test_interaction_entities_and_attributes():
    out = extract('interactions', interaction('inhibitor') + interaction('untyped placeholder'))
    relation = out.relations[0]
    assert annotations(relation, 'supporting_data_source') == ['TALC']
    assert annotations(relation, 'dgidb:drug_is_approved') == ['true']
    assert annotations(relation, 'dgidb:drug_is_antineoplastic') == ['true']
    scores = [a for a in relation.annotations if a['term'] == 'dgidb:interaction_score']
    assert scores[0]['quantity']['has_numeric_value'] == pytest.approx(0.2024766898977484)
    drugs = {e.label: e for e in out.entities.values() if e.entity_type == 'chemical_entity'}
    assert {(i['ns'], i['id']) for i in drugs['IMATINIB'].identifiers} == {
        ('rxcui', '282388'), (Namespace.NAME, 'IMATINIB'), (Namespace.SYNONYM, 'IMATINIB'),
    }
    # The concept without a name has the name the source gave, and the ChEMBL
    # accession without its prefix.
    chembl = drugs["17ALPHA-3'-BROMOBENZYL ESTRADIOL"]
    assert {(i['ns'], i['id']) for i in chembl.identifiers} == {
        (Namespace.CHEMBL, 'CHEMBL1627429'), (Namespace.NAME, "17ALPHA-3'-BROMOBENZYL ESTRADIOL"),
    }
    assert not any(i['id'].startswith('CHEMBL:') for e in drugs.values() for i in e.identifiers)


def test_flags_go_through_the_true_false_crosswalk():
    out = extract('interactions', interaction('binder') + interaction('unmapped drug'))
    flags = ('dgidb:drug_is_approved', 'dgidb:drug_is_immunotherapy', 'dgidb:drug_is_antineoplastic')
    binder, unmapped = out.relations
    assert [annotations(binder, flag) for flag in flags] == [['false'], ['false'], ['false']]
    # A drug DGIdb did not map has no flags, and none is invented.
    assert [annotations(unmapped, flag) for flag in flags] == [[], [], []]
    assert dgidb.YES_NO == {'true': 'true', 'false': 'false'}


def test_drug_namespaces_are_namespace_members():
    out = extract('interactions', interaction('antisense'))
    drug = [e for e in out.entities.values() if e.entity_type == 'chemical_entity'][0]
    assert (Namespace.NCIT, 'C2721') in {(i['ns'], i['id']) for i in drug.identifiers}
    assert dgidb.DRUG_NAMESPACES['rxcui'] is Namespace.RXCUI
    assert dgidb.DRUG_NAMESPACES['drugsatfda.anda'] is Namespace.DRUGSATFDA_ANDA


def test_the_downloads_are_pinned_to_a_release():
    assert dgidb.RELEASE == '2026-09'
    for download in (dgidb.download_categories, dgidb.download_genes, dgidb.download_interactions):
        assert f'/releases/download/{dgidb.RELEASE}/' in download.url
        assert 'latest' not in download.url
    assert dgidb.RELEASE in dgidb.config.summary and dgidb.RELEASE in dgidb.config.description


def test_unmapped_gene_gives_no_statement_and_unknown_values_raise():
    gene_less = interaction('unmapped gene')[0]
    assert dgidb.resource.interactions.mapper(gene_less) is None
    good = interaction('inhibitor')[0]
    with pytest.raises(ValueError, match='Unknown DGIdb interaction type'):
        dgidb.resource.interactions.mapper(good | {'interaction_type': 'enhancer'})
    with pytest.raises(ValueError, match='Unknown DGIdb source'):
        dgidb.resource.interactions.mapper(good | {'source': 'NewDatabase'})
    with pytest.raises(ValueError, match='Unknown DGIdb concept'):
        dgidb.resource.interactions.mapper(good | {'drug_concept_id': 'pubchem:1'})
    with pytest.raises(ValueError, match='Unknown DGIdb concept'):
        dgidb.resource.interactions.mapper(good | {'gene_concept_id': 'foo:1'})
    with pytest.raises(ValueError, match='Unknown DGIdb value'):
        dgidb.resource.interactions.mapper(good | {'drug_is_approved': 'maybe'})
    with pytest.raises(ValueError, match='Unknown DGIdb value'):
        dgidb.resource.interactions.mapper(good | {'drug_is_antineoplastic': 'True'})


def test_no_direction_is_invented():
    out = extract('interactions', parse_interactions())
    for relation in out.relations:
        if relation.predicate == 'interacts_with':
            assert annotations(relation, 'object_direction_qualifier') == []
