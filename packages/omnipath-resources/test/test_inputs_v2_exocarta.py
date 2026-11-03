"""ExoCarta: parser, mapper and extractor on lines of the release 5 tables."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import exocarta
from pypath.inputs_v2.exocarta import TAXA, VESICLE_TYPE_NAMESPACE, resource
from pypath.inputs_v2.parsers.exocarta import (
    EXOCARTA,
    iter_detections,
    iter_vesicle_types,
)

# Lines of EXOCARTA_EXPERIMENT_DETAILS_5.txt.
EXPERIMENT_HEADER = (
    'EXPERIMENT ID\tPUBMED ID\tSPECIES\tEXPERIMENT DESCRIPTION\tSAMPLE\tSAMPLE SOURCE\t'
    'SAMPLE NAME\tIDENTIFICATIONS\tMETHODS\tYEAR\tISOLATION METHOD'
)
EXPERIMENTS = [
    # A cell line sample of mouse dendritic cells.
    '3\t10545503\tMus musculus\tMolecular characterization of dendritic cell-derived exosomes. '
    'Selective accumulation of the heat shock protein hsc73.\tDendritic cells\tcl\tD1\tProtein\t'
    'Mass spectrometry [MALDI TOF]|Western blotting\t1999\t'
    'Differential centrifugation|Sucrose density gradient',
    # Two species, two sample sources.
    '15\t17486113\tMus musculus|Homo sapiens\tExosome-mediated transfer of mRNAs and microRNAs '
    'is a novel mechanism of genetic exchange between cells.\tMast cells\tcl|ts\t'
    'MC9|Bone marrow-derived mast cells|HMC-1\tProtein|mRNA|miRNA\t'
    'Mass spectrometry [MALDI TOF]|Western blotting|Microarray|miRCURY LNA Array\t2007\t'
    'Filtration|Ultracentrifugation|Sucrose density gradient',
    # No publication (PubMed ID 0).
    '126\t0\tHomo sapiens\tMesenchymal Stem Cell Exosomes: The Future MSC-based Therapy?\t'
    'Mesenchymal stem cells\tcl\thuES9.E1\tProtein\t'
    'Western blotting|Antobody array|Mass spectrometry\t2011\tHPLC',
    # No sample source (0).
    '202\t23035643\tDrosophila melanogaster\tDrosophila S2 cells secrete Wingless on exosome-like '
    'vesicles but the Wingless gradient forms independently of exosomes.\tDrosophila S2 cells\t0\t'
    'Normal-Drosophila S2 tub-Wg\tProtein\tMass spectrometry\t2012\t'
    'Filtration|Ultracentrifugation|Sucrose density gradient',
    '222\t23663360\tHomo sapiens\tCharacterization of human plasma-derived exosomal RNAs by deep '
    'sequencing\tPlasma\tts\tPlasma - Normal\tmiRNA|tRNA|snRNA|snoRNA|rRNA|piRNA|lncRNA\t'
    'RT-PCR\t2013\tExoquick',
]
# Lines of EXOCARTA_PROTEIN_MRNA_DETAILS_5.txt.
DETECTION_HEADER = 'CONTENT ID\tCONTENT TYPE\tENTREZ GENE ID\tGENE SYMBOL\tSPECIES\tEXPERIMENT ID\tMETHODS'
ITGAM = '4\tprotein\t16409\tItgam\tMus musculus\t3\tMass spectrometry|Western blotting'
# The content type "mrna" in lower case.
ACP1 = '3693\tmrna\t11431\tAcp1\tMus musculus\t15\tMicroarray'
A2M = '13638\tprotein\t2\tA2M\tHomo sapiens\t126\tMass spectrometry [LTQ-FT Ultra]'
# The species is the genus only.
CG30382 = '36131\tprotein\t246582\tCG30382\tDrosophila\t202\tMass spectrometry'
SNORA45A = '48273\tsnoRNA\t619562\tSNORA45A\tHomo sapiens\t222\tRT-PCR'
DETECTIONS = [ITGAM, ACP1, A2M, CG30382, SNORA45A]


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
            database=EXOCARTA,
        )
    )


def process(dataset, row):
    out = SilverExtractor('exocarta', dataset)
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


def test_import_defines_the_datasets_and_the_archive_urls_without_downloading():
    assert set(resource.datasets()) == {'terms', 'detections'}
    assert 'web.archive.org' in exocarta.detection_download.url
    assert exocarta.config.license.name == 'CC_BY_4_0'


def test_a_detection_has_the_fields_of_its_experiment_and_the_exosome_type():
    itgam = rows()[0]
    assert itgam == {
        'entrez': '16409', 'symbol': 'Itgam', 'species': 'Mus musculus', 'content_type': 'protein',
        'methods': ['Mass spectrometry', 'Western blotting'], 'vesicle_type': 'Exosomes',
        'experiment': '3', 'pmid': '10545503', 'sample': 'Dendritic cells',
        'sample_name': 'D1', 'sample_source': ['cl'],
        'isolation_methods': ['Differential centrifugation', 'Sucrose density gradient'],
    }


def test_statement_is_located_in_exosomes_with_the_evidence_of_the_experiment():
    gene, relation, vesicle = relation_of(process('detections', rows()[0]))
    assert relation.predicate == 'located_in'
    assert (gene.entity_type, gene.namespace, gene.identifier, gene.taxon) == (
        'gene', 'entrez', '16409', '10090'
    )
    assert {(i['ns'], i['id']) for i in gene.identifiers} == {
        (Namespace.ENTREZ, '16409'), (Namespace.GENESYMBOL, 'Itgam'),
    }
    assert (vesicle.entity_type, vesicle.namespace, vesicle.identifier) == (
        'ontology_class', VESICLE_TYPE_NAMESPACE, 'Exosomes'
    )
    assert attributes(relation, 'publications') == ['PMID:10545503']
    assert attributes(relation, 'exocarta:content_type') == ['protein']
    assert attributes(relation, 'exocarta:experiment') == ['3']
    assert attributes(relation, 'exocarta:sample') == ['Dendritic cells']
    assert attributes(relation, 'exocarta:sample_name') == ['D1']
    assert attributes(relation, 'exocarta:sample_source') == ['cl']
    assert attributes(relation, 'exocarta:isolation_method') == [
        'Differential centrifugation', 'Sucrose density gradient',
    ]
    assert attributes(relation, 'exocarta:detection_method') == [
        'Mass spectrometry', 'Western blotting',
    ]


def test_the_species_is_that_of_the_molecule_not_of_the_experiment():
    # Experiment 15 used mouse and human cells; the molecule row says mouse.
    acp1 = rows()[1]
    assert acp1['species'] == 'Mus musculus'
    gene, relation, _ = relation_of(process('detections', acp1))
    assert gene.taxon == '10090'
    assert attributes(relation, 'exocarta:sample_source') == ['cl', 'ts']
    assert attributes(relation, 'exocarta:content_type') == ['mRNA']


def test_a_publication_of_zero_and_a_sample_source_of_zero_are_none():
    a2m, cg = rows()[2], rows()[3]
    assert a2m['pmid'] is None
    _, relation, _ = relation_of(process('detections', a2m))
    assert attributes(relation, 'publications') == []
    assert cg['sample_source'] == []
    _, relation, _ = relation_of(process('detections', cg))
    assert attributes(relation, 'exocarta:sample_source') == []


def test_the_genus_drosophila_is_the_genus_taxon_and_nothing_is_guessed():
    assert TAXA['Drosophila'] == 'NCBITaxon:7215'
    gene, _, _ = relation_of(process('detections', rows()[3]))
    assert gene.taxon == '7215'


def test_rna_types_are_kept_as_the_content_type_of_the_gene_statement():
    snora = rows()[4]
    assert snora['content_type'] == 'snoRNA'
    gene, relation, _ = relation_of(process('detections', snora))
    # The subject is the gene; the molecule detected is an attribute, not a protein.
    assert gene.entity_type == 'gene'
    assert attributes(relation, 'exocarta:content_type') == ['snoRNA']


def test_the_vesicle_type_term_has_a_name_and_the_object_matches_it():
    (term,) = list(iter_vesicle_types(opener(EXPERIMENT_HEADER, *EXPERIMENTS), database=EXOCARTA))
    assert term == {'name': 'Exosomes'}
    out = SilverExtractor('exocarta', 'terms')
    out.process_record(resource.terms.mapper(term), term, 'fixture:1', 1)
    (entity,) = out.entities.values()
    assert {(i['ns'], i['id']) for i in entity.identifiers} == {
        (VESICLE_TYPE_NAMESPACE, 'Exosomes'), (Namespace.NAME, 'Exosomes'),
    }
    _, _, vesicle = relation_of(process('detections', rows()[0]))
    assert (entity.entity_type, entity.namespace, entity.identifier) == (
        vesicle.entity_type, vesicle.namespace, vesicle.identifier
    )


def test_nothing_is_asserted_beyond_exosomes_and_no_protein_entity_is_made():
    for row in rows():
        out = process('detections', row)
        assert {e.entity_type for e in out.entities.values()} == {'gene', 'ontology_class'}


@pytest.mark.parametrize(
    'detections',
    [
        ['4\tprotein\t16409\tItgam\tMus musculus\t999\tMass spectrometry'],  # unknown experiment
        ['4\tlipid\t16409\tItgam\tMus musculus\t3\tMass spectrometry'],  # unknown content type
        [],  # no rows
    ],
)
def test_parser_raises_on_what_it_does_not_know(detections):
    with pytest.raises(ValueError):
        rows(detections=detections)


def test_parser_raises_on_a_repeated_experiment_and_a_changed_header():
    with pytest.raises(ValueError, match='listed twice'):
        rows(experiments=[EXPERIMENTS[0], EXPERIMENTS[0]])
    with pytest.raises(ValueError, match='lacks expected columns'):
        list(
            iter_detections(
                opener(DETECTION_HEADER.replace('SPECIES', 'ORGANISM'), ITGAM),
                experiment_download=Fake('\n'.join([EXPERIMENT_HEADER, *EXPERIMENTS]) + '\n'),
                database=EXOCARTA,
            )
        )


def test_an_unknown_sample_source_raises():
    # The tables write cl, ts and fl (and 0 for none); nothing else is known.
    changed = EXPERIMENTS[0].replace('\tcl\tD1\t', '\tzz\tD1\t')
    with pytest.raises(ValueError, match='Unknown sample source'):
        rows(experiments=[changed, *EXPERIMENTS[1:]])


def test_unknown_species_raises_in_the_mapper():
    with pytest.raises(ValueError, match='Unknown exocarta species'):
        resource.detections.mapper({**rows()[0], 'species': 'Pan troglodytes'})


def test_the_release_read_is_stated_in_the_summary_and_the_description():
    # Release 5 of 2015 from the Internet Archive; release 6 is not reachable.
    assert 'release 5' in exocarta.config.summary and '2015' in exocarta.config.summary
    assert 'release 5' in exocarta.config.description
    assert all(
        'web.archive.org' in d.url and d.url.endswith('_5.txt')
        for d in (exocarta.experiment_download, exocarta.detection_download)
    )
