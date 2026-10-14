"""Vesiclepedia: parser, mapper and extractor on lines of the release 5.1 tables."""

import inspect
import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import _vesicles, exocarta, vesiclepedia
from pypath.inputs_v2.parsers.exocarta import VESICLEPEDIA, iter_detections, iter_vesicle_types
from pypath.inputs_v2.vesiclepedia import TAXA, VESICLE_TYPE_NAMESPACE, resource

# Lines of VESICLEPEDIA_EXPERIMENT_DETAILS_5.1.txt.
EXPERIMENT_HEADER = (
    'EXPERIMENT ID\tPUBMED ID\tSPECIES\tEXPERIMENT DESCRIPTION\tSAMPLE\tSAMPLE SOURCE\t'
    'SAMPLE NAME\tIDENTIFICATIONS\tMETHODS\tYEAR\tISOLATION METHOD\tVESICLE TYPE'
)
EXPERIMENTS = [
    # Two vesicle types, separated by a slash.
    '376\t18346024\tHomo sapiens\tThe proteome of red cell membranes and vesicles during '
    'storage in blood bank conditions.\tRed blood cells\tcl\t'
    'Microvesicle, Storage day 21-Red blood cells\tProtein\tMass spectrometry\t2008\t-\t'
    'Microvesicles/Nanovesicles',
    # A trailing space in the species.
    '185\t22116803\tEquus caballus \tCell-Secreted Vesicles in Equine Ovarian Follicular Fluid '
    'Contain miRNAs and Proteins: A Possible New Form of Cell Communication Within the '
    'Ovarian Follicle\tFollicular fluid\tts\tFollicular fluid - Ovarian - Young mare\t'
    'Protein|miRNA\tWestern blotting|RT-PCR|Mass spectrometry [LTQ]\t2011\t'
    'Differential centrifugation|ExoQuick\tExosomes',
    # A trailing space in the vesicle type; isolation method "-" is none.
    '1499\t33982667\tHomo sapiens\tExosomal circRNA 0001445 promotes glioma progression '
    'through miRNA-127-5p/SNX5 pathway\tGlioblastoma cells\tcl\tA172\tProtein\t'
    'Western blotting\t2021\tDifferential centrifugation|Ultracentrifugation\tExosomes ',
    '694\t28842968\tHomo sapiens\tInsulin Mediated Activation of PI3K/Akt Signalling Pathway '
    'Modifies the Proteomic Cargo of Extracellular Vesicles.\tColorectal cancer cells\tcl\t'
    'LIM1215- insulin induced\tProtein\tWestern blotting\t2017\t'
    'Differential centrifugation|Ultracentrifugation\tMicrovesicles',
    '317\t19266099\tHomo sapiens\tTransfer of microRNAs by embryonic stem cell microvesicles\t'
    'Embryonic stem cells\tts\tEmbryonic stem cells\tmRNA|miRNA|snRNA\tRT-PCR\t2009\t'
    'Differential centrifugation\tMicrovesicles',
]
# Lines of VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt.
DETECTION_HEADER = 'CONTENT ID\tCONTENT TYPE\tENTREZ GENE ID\tGENE SYMBOL\tSPECIES\tEXPERIMENT ID\tMETHODS'
A1BG = '45836\tprotein\t1\tA1BG\tHomo sapiens\t376\tMass spectrometry'
HSPA8 = '18091\tprotein\t100009679\tHSPA8\tEquus caballus \t185\tWestern blotting'
# The content type "protein " with a trailing space, "mrna" and "snrna".
ACTB = '492992\tprotein \t60\tACTB\tHomo sapiens\t694\tWestern blotting'
ACTB_MRNA = '41438\tmrna\t60\tACTB\tHomo sapiens\t317\tRT-PCR'
SNRNA = '41439\tsnrna\t26826\tRNU6-2\tHomo sapiens\t317\tRT-PCR'
DETECTIONS = [A1BG, HSPA8, ACTB, ACTB_MRNA, SNRNA]


def opener(header, *lines):
    return SimpleNamespace(result=io.StringIO('\n'.join([header, *lines]) + '\n'))


class Fake:
    def __init__(self, text):
        self.text = text

    def open(self, force_refresh=False, **_):
        return SimpleNamespace(result=io.StringIO(self.text))


def rows(detections=DETECTIONS, experiments=EXPERIMENTS):
    return list(
        iter_detections(
            opener(DETECTION_HEADER, *detections),
            experiment_download=Fake('\n'.join([EXPERIMENT_HEADER, *experiments]) + '\n'),
            database=VESICLEPEDIA,
        )
    )


def process(dataset, row):
    out = SilverExtractor('vesiclepedia', dataset)
    out.process_record(resource.datasets()[dataset].mapper(row), row, 'fixture:1', 1)
    return out


def relation_of(out):
    (relation,) = out.relations
    return (
        out.entities[relation.subject_entity_key],
        relation,
        out.entities[relation.object_entity_key],
    )


def attributes(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term and a['scope'] == 'relation']


def test_import_defines_the_datasets_without_downloading():
    assert set(resource.datasets()) == {'terms', 'detections'}
    assert 'web.archive.org' in vesiclepedia.detection_download.url
    assert vesiclepedia.config.license.name == 'CC_BY_4_0'


def test_a_molecule_of_an_experiment_with_two_vesicle_types_is_two_rows():
    parsed = rows()
    assert [(r['symbol'], r['vesicle_type']) for r in parsed[:2]] == [
        ('A1BG', 'Microvesicles'), ('A1BG', 'Nanovesicles'),
    ]
    out = SilverExtractor('vesiclepedia', 'detections')
    for index, row in enumerate(parsed[:2]):
        out.process_record(resource.detections.mapper(row), row, f'fixture:{index}', index)
    objects = {out.entities[r.object_entity_key].identifier for r in out.relations}
    assert objects == {'Microvesicles', 'Nanovesicles'}
    assert {r.predicate for r in out.relations} == {'located_in'}


def test_spaces_around_values_are_removed():
    parsed = rows()
    horse = next(r for r in parsed if r['symbol'] == 'HSPA8')
    assert horse['species'] == 'Equus caballus'
    gene, relation, _ = relation_of(process('detections', horse))
    assert gene.taxon == '9796'
    # "Exosomes " of experiment 1499 is "Exosomes".
    experiments = list(iter_vesicle_types(opener(EXPERIMENT_HEADER, *EXPERIMENTS), database=VESICLEPEDIA))
    assert [t['name'] for t in experiments].count('Exosomes') == 1


def test_the_content_types_are_normalized_by_exact_case_insensitive_names():
    parsed = rows()
    types = {r['symbol']: r['content_type'] for r in parsed if r['symbol'] in ('ACTB', 'RNU6-2')}
    assert types == {'ACTB': 'mRNA', 'RNU6-2': 'snRNA'}
    actb_protein = next(r for r in parsed if r['symbol'] == 'ACTB' and r['experiment'] == '694')
    assert actb_protein['content_type'] == 'protein'


def test_statement_has_the_evidence_of_the_experiment_and_the_species_of_the_row():
    parsed = rows()
    row = next(r for r in parsed if r['symbol'] == 'ACTB' and r['experiment'] == '694')
    gene, relation, vesicle = relation_of(process('detections', row))
    assert relation.predicate == 'located_in'
    assert (gene.entity_type, gene.namespace, gene.identifier, gene.taxon) == (
        'gene', 'entrez', '60', '9606'
    )
    assert (vesicle.entity_type, vesicle.namespace, vesicle.identifier) == (
        'ontology_class', VESICLE_TYPE_NAMESPACE, 'Microvesicles'
    )
    assert attributes(relation, 'publications') == ['PMID:28842968']
    assert attributes(relation, 'vesiclepedia:sample') == ['Colorectal cancer cells']
    assert attributes(relation, 'vesiclepedia:sample_name') == ['LIM1215- insulin induced']
    assert attributes(relation, 'vesiclepedia:detection_method') == ['Western blotting']
    # The prefix names the database: nothing is attributed to ExoCarta.
    assert not [a for a in relation.annotations if a['term'].startswith('exocarta')]


def test_an_isolation_method_of_a_dash_is_none():
    parsed = rows()
    a1bg = parsed[0]
    assert a1bg['isolation_methods'] == []
    _, relation, _ = relation_of(process('detections', a1bg))
    assert attributes(relation, 'vesiclepedia:isolation_method') == []


def test_vesicle_types_are_terms_with_names_each_once_and_the_objects_match():
    terms = list(iter_vesicle_types(opener(EXPERIMENT_HEADER, *EXPERIMENTS), database=VESICLEPEDIA))
    assert [t['name'] for t in terms] == ['Microvesicles', 'Nanovesicles', 'Exosomes']
    defined = {}
    for term in terms:
        out = SilverExtractor('vesiclepedia', 'terms')
        out.process_record(resource.terms.mapper(term), term, 'fixture:1', 1)
        (entity,) = out.entities.values()
        assert {(i['ns'], i['id']) for i in entity.identifiers} == {
            (VESICLE_TYPE_NAMESPACE, term['name']), (Namespace.NAME, term['name']),
        }
        defined[entity.identifier] = (entity.entity_type, entity.namespace, entity.identifier)
    for row in rows():
        _, _, vesicle = relation_of(process('detections', row))
        assert defined[vesicle.identifier] == (
            vesicle.entity_type, vesicle.namespace, vesicle.identifier
        )


def test_no_hierarchy_between_vesicle_types_is_asserted():
    # "Microvesicles/Nanovesicles" is two types, not a type and a subtype.
    entity = resource.terms.mapper({'name': 'Nanovesicles'})
    assert entity.ontology_relations is None


def test_species_cover_every_species_of_the_molecule_table_and_nothing_is_guessed():
    assert TAXA['Homo sapiens'] == 'NCBITaxon:9606'
    assert TAXA['Salmo salar'] == 'NCBITaxon:8030'
    with pytest.raises(ValueError, match='Unknown vesiclepedia species'):
        resource.detections.mapper({**rows()[0], 'species': 'Ixodes ricinus'})


@pytest.mark.parametrize(
    'detections',
    [
        ['1\tprotein\t1\tA1BG\tHomo sapiens\t99999\tMass spectrometry'],  # unknown experiment
        ['1\tmetabolite\t1\tA1BG\tHomo sapiens\t376\tMass spectrometry'],  # unknown content type
        [],
    ],
)
def test_parser_raises_on_what_it_does_not_know(detections):
    with pytest.raises(ValueError):
        rows(detections=detections)


def test_parser_requires_the_vesicle_type_column():
    header = EXPERIMENT_HEADER.replace('\tVESICLE TYPE', '')
    with pytest.raises(ValueError, match='lacks expected columns'):
        list(iter_vesicle_types(opener(header, EXPERIMENTS[0].rsplit('\t', 1)[0]), database=VESICLEPEDIA))


def test_the_mapping_is_in_a_shared_private_module_and_no_resource_imports_the_other():
    assert vesiclepedia.map_detection is _vesicles.map_detection is exocarta.map_detection
    assert vesiclepedia.map_vesicle_type is _vesicles.map_vesicle_type is exocarta.map_vesicle_type
    assert 'inputs_v2.exocarta' not in inspect.getsource(vesiclepedia)
    assert 'inputs_v2.vesiclepedia' not in inspect.getsource(exocarta)


def test_a_dash_is_no_gene_symbol():
    # Vesiclepedia writes "-" as the symbol of genes without one (19,280 rows of
    # release 5.1, such as "71693  protein  7996741  -  Acinetobacter baumannii");
    # here the ACTB row with its symbol replaced by the dash.
    (row,) = [r for r in rows(detections=[ACTB.replace('\tACTB\t', '\t-\t')]) if r['experiment'] == '694']
    assert row['symbol'] is None and row['entrez'] == '60'
    gene, _, _ = relation_of(process('detections', row))
    assert [(i['ns'], i['id']) for i in gene.identifiers] == [('entrez', '60')]
    assert gene.label != '-'

