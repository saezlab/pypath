"""IntOGen driver genes and cancer types: parsers, mappers and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import intogen
from pypath.inputs_v2.parsers.intogen import iter_cancer_types, iter_drivers

DRIVER_HEADER = (
    'SYMBOL\tTRANSCRIPT\tCOHORT\tCANCER_TYPE\tMETHODS\tMUTATIONS\tSAMPLES\t%_SAMPLES_COHORT\tQVALUE_COMBINATION\tROLE'
    '\tDOMAINS\t2D_CLUSTERS\t3D_CLUSTERS\tEXCESS_MIS\tEXCESS_NON\tEXCESS_SPL\tIS_DRIVER\tTOTAL_SAMPLES'
)
# Lines of Compendium_Cancer_Genes.tsv (release 2024.09.20).
DRIVERS = [
    'ABCC4\tENST00000645237\tCBIOP_WXS_ANGS_UNTREAT_2020\tANGS\toncodriveclustl,hotmaps\t4\t4.0\t0.1176470588235294\t1.833927821088493e-05\tAct\t\t1250:1250\t\t0.9464674697705246\t0.0\t0.0\tTrue\t34',
    'ABL1\tENST00000318560\tTCGA_WXS_UCEC\tUCEC\tcombination\t13\t10.0\t0.0194174757281553\t0.0011081336865465\tLoF\t\t\t\t0.5635097412592688\t0.7703173826926656\t0.0\tTrue\t515',
    'AKAP9\tENST00000356239\tCBIOP_WXS_BCC_UNIGE_2016_UNTREAT\tBCC\tmutpanning\t30\t20.0\t0.2631578947368421\t0.3246789478921122\tambiguous\t\t\t\t0.0\t0.0\t0.0\tTrue\t76',
]
COHORT_HEADER = (
    'COHORT\tCANCER\tCANCER_NAME\tSOURCE\tPLATFORM\tREFERENCE\tTYPE\tTREATED\tAGE\tSAMPLES\tMUTATIONS\tCOHORT_NICK'
    '\tCOHORT_NAME\tDRIVERS'
)
COHORTS = [
    'TCGA_WXS_ACC\tACC\tAdrenocortical Carcinoma\tTCGA/PanCatAtlas, phs000178\tWXS\tPMID:24071849\tPrimary\tUntreated\tAdult\t91\t11660\tTCGA_WXS_ACC\tAdrenocortical Carcinoma from TCGA/PanCatAtlas, phs000178\t15.0',
    'STJUDE_WGS_D_ACC_2018\tACC\tAdrenocortical Carcinoma\tSTJUDE\tWGS\tWEB: https://www.stjude.cloud/\tPrimary\tUntreated\tAdult\t20\t17560\tSTJUDE_WGS_D_ACC_2018\tAdrenocortical Carcinoma from STJUDE\t3.0',
    'TCGA_WXS_UCEC\tUCEC\tUterine Corpus Endometrial Carcinoma\tTCGA/PanCatAtlas, phs000178\tWXS\tPMID:24071849\tPrimary\tUntreated\tAdult\t515\t1\tTCGA_WXS_UCEC\tx\t1.0',
]


def archive(name, header, lines):
    # The archive's members are in a folder named after the release.
    text = '\n'.join([header, *lines]) + '\n'
    return SimpleNamespace(result={'2024-06-18_IntOGen/': None, f'2024-06-18_IntOGen/{name}': io.StringIO(text)})


def drivers(lines=DRIVERS, header=DRIVER_HEADER):
    return list(iter_drivers(archive('Compendium_Cancer_Genes.tsv', header, lines)))


def cancer_types(lines=COHORTS, header=COHORT_HEADER):
    return list(iter_cancer_types(archive('cohorts.tsv', header, lines)))


def extract(dataset, records):
    out = SilverExtractor('intogen', dataset)
    mapper = getattr(intogen.resource, dataset).mapper
    for index, record in enumerate(records):
        mapped = mapper(record)
        if mapped is not None:
            out.process_record(mapped, record, f'fixture:{index}', index)
    return out


def values(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term and a.get('scope') == 'relation']


def test_parser_gives_one_record_per_gene_and_cohort():
    records = drivers()
    assert len(records) == 3
    assert records[0] == {
        'symbol': 'ABCC4',
        'transcript': 'ENST00000645237',
        'cohort': 'CBIOP_WXS_ANGS_UNTREAT_2020',
        'cancer_type': 'ANGS',
        'methods': ['oncodriveclustl', 'hotmaps'],
        'mutations': '4',
        'samples': '4.0',
        'qvalue': '1.833927821088493e-05',
        'role': 'Act',
    }
    assert records[1]['methods'] == ['combination']


def test_parser_rejects_a_changed_archive_and_table():
    with pytest.raises(ValueError, match='lacks expected columns'):
        drivers(header=DRIVER_HEADER.replace('ROLE', 'MODE'))
    with pytest.raises(ValueError, match='has no rows'):
        drivers(lines=[])
    with pytest.raises(ValueError, match='IS_DRIVER'):
        drivers([DRIVERS[0].replace('\tTrue\t', '\tFalse\t')])
    with pytest.raises(ValueError, match='has no cohorts.tsv'):
        list(iter_cancer_types(archive('other.tsv', COHORT_HEADER, COHORTS)))
    with pytest.raises(ValueError, match='two names'):
        cancer_types([COHORTS[0], COHORTS[0].replace('Adrenocortical Carcinoma', 'Adrenal')])


def test_cancer_types_are_listed_once_with_their_names():
    types = cancer_types()
    assert types == [
        {'cancer_type': 'ACC', 'name': 'Adrenocortical Carcinoma'},
        {'cancer_type': 'UCEC', 'name': 'Uterine Corpus Endometrial Carcinoma'},
    ]
    term = intogen.resource.terms.mapper(types[0])
    assert term.type == 'ontology_class'
    assert term.identifiers == [('oncotree', 'ACC'), (Namespace.NAME, 'Adrenocortical Carcinoma')]
    assert term.identifiers[0][0] is Namespace.ONCOTREE
    assert set(intogen.resource.datasets()) == {'terms', 'drivers'}


def test_driver_gene_is_associated_with_the_cancer_type_of_its_cohort():
    out = extract('drivers', drivers())
    assert [r.predicate for r in out.relations] == ['associated_with'] * 3
    first = out.relations[0]
    assert values(first, 'intogen:cohort') == ['CBIOP_WXS_ANGS_UNTREAT_2020']
    assert values(first, 'intogen:role') == ['Act']
    assert values(first, 'intogen:method') == ['oncodriveclustl', 'hotmaps']
    assert values(first, 'intogen:transcript') == ['ENST00000645237']
    q = [a for a in first.annotations if a['term'] == 'intogen:qvalue_combination'][0]
    assert q['quantity']['has_numeric_value'] == pytest.approx(1.833927821088493e-05)
    # The role the source calls ambiguous stays as written.
    assert values(out.relations[2], 'intogen:role') == ['ambiguous']
    genes = {e.identifier: e for e in out.entities.values() if e.entity_type == 'gene'}
    assert set(genes) == {'ABCC4', 'ABL1', 'AKAP9'}
    assert genes['ABL1'].taxon == '9606'


def test_object_is_the_same_entity_as_the_term():
    term = intogen.resource.terms.mapper(cancer_types()[1])
    statement = intogen.resource.drivers.mapper(drivers()[1])
    assert statement.object.type == term.type
    assert statement.object.identifiers[0] == term.identifiers[0] == ('oncotree', 'UCEC')


def test_unknown_values_raise_and_nothing_else_is_inferred():
    with pytest.raises(ValueError, match='Unknown IntOGen role'):
        intogen.resource.drivers.mapper(drivers()[0] | {'role': 'Activating'})
    with pytest.raises(ValueError, match='Unknown IntOGen method'):
        intogen.resource.drivers.mapper(drivers()[0] | {'methods': ['newmethod']})
    out = extract('drivers', drivers())
    # No causal predicate: the catalog says the gene is a driver of the cancer type,
    # not that it affects it, and the mode of action is a value of the source.
    assert {r.predicate for r in out.relations} == {'associated_with'}
    assert not any(a['term'].endswith('_qualifier') for r in out.relations for a in r.annotations)
