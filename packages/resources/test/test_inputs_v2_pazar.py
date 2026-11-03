"""Source-specific Biolink contract for the PAZAR input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import pazar as mod
from pypath.inputs_v2.parsers.pazar import COLUMNS, iter_tf_targets

# Lines copied from the files of tftargets.zip (tab separated, no header),
# with the spaces at the end of a cell as they are in the file.
LINES = {
    'chip': 'TF0001150\tENST00000264010\tCTCF_HUMAN\tGS0367173\tENSG00000230368\t1\t803451\t812283\tHomo sapiens\tCTCF_Cui\t21795385\tPROTEIN BINDING ASSAY::CHROMATIN IMMUNOPRECIPITATION (CHIP)',
    'tag': 'TF0001127\tENSMUST00000039174\tPRDM1_MOUSE\tGS0014801\tENSMUSG00000021239\t12\t85910790\t85936406\tMus musculus\tCanEuCre MaxiPromoter \t20463215\tPROTEIN BINDING ASSAY :: CHROMATIN IMMUNOPRECIPITATION WITH TAG SEQUENCING (CHIP-TS)',
    'emsa': 'TF0000841\tENSMUST00000043819\tSIX3_MOUSE\tGS0014253\tENSG00000164690\t7\t155592680\t155604967\tHomo sapiens\tTFe\t18836447\tELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::SUPERSHIFT ',
    'no_assay': 'TF0000324\tENST00000265354\tSRF_HUMAN\tGS0000616\tENSG00000143632\t1\t229566992\t229569845\tHomo sapiens\tABS\t12760745\tNA',
    'predicted': 'TF0001144\tENST00000263734\tHIFALPHA\tGS0014831\tENSG00000134333\t11\t18415935\t18429972\tHomo sapiens\tOrtiz_Barahona\t20061373\tCOMPUTATIONAL PREDICTION',
    'conserved': 'TF0000586\tENST00000361642\tZEB1_HUMAN\tGS0001209\tENSG00000184363\t11\t392614\t404908\tHomo sapiens\tORegAnno\t17391671\tSEQUENCE CONSERVATION::ORTHOLOGOUS GENE CONSERVATION',
    'curated': 'TF0000699\tENST00000397062\tNFE2L2(NRF2)\tGS0013552\tENSMUSG00000003849\t8\t109912125\t109927106\tMus musculus\tAREs\t12816537\tEXPERT CURATED::LITERATURE DERIVED',
    'no_pmid': 'TF0000696\tENST00000316099\tHNF4A\tGS0013435\tENSG00000135100\t12\t121416346\t121441746\tHomo sapiens\tHNF4\t\tEXPERT CURATED::LITERATURE DERIVED',
    'gene_space': 'TF0000212\tENST00000310015\tSP3_HUMAN\tGS0000425\tENSG00000117983 \t11\t1244296\t1284402\tHomo sapiens\tMUC5AC\t12077147\tELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::SUPERSHIFT ',
    'worm': 'TF0000584\tC33D3.1\tC33D3.1_WORM\tGS0001206\tY66H1B.4\tIV\t379900\t383185\tCaenorhabditis elegans\tORegAnno\t15734735\tPROTEIN BINDING ASSAY::YEAST 1-HYBRID ASSAY',
    'fly_upper': 'TF0000941\tFBGN0015664\tDREF_Drosophila melanogaster\tGS0014589\tFBGN0039044\t-\t-\t-\tDrosophila melanogaster\tTFe\t20101238\tPROTEIN BINDING ASSAY::CHROMATIN IMMUNOPRECIPITATION (CHIP)',
    'tf_digits': 'TF0001032\t00000416339\tPAX6\tGS0014683\tENSG00000189056\t7\t103112231\t103629963\tHomo sapiens\tPax6_Aniridia\t17014839\tEMSA',
    'tf_gene': 'TF0000916\tENSRNOG00000004642\tRAT_LMX1A\tGS0014442\tENSRNOG00000012052\t1\t258001134\t258001684\tRattus norvegicus\tDiabetes\t1358758\tELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::DIRECT GEL SHIFT  ',
    'pmid_in_assay': 'TF0000812\tENST00000343882\tFOXO3_HUMAN\tGS0014040\tENSG00000153094\t2\t111876955\t111924587\tHomo sapiens\tTFe\t15688014\t15688014',
    'pmid_nine_digits': 'TF0000479\tENSMUST00000020156\tMYB_MOUSE\tGS0000961\tENSMUSG00000030102\t6\t108163112\t108501103\tMus musculus\tpleiades genes\t117363689\tPROTEIN BINDING ASSAY::CHROMATIN IMMUNOPRECIPITATION (CHIP)',
    'reporter': 'TF0000561\tENSGALT00000010031\tMYOD1_CHICKEN\tGS0001120\tENSGALG00000002907\t7\t2758273\t2775953\tGallus gallus\tORegAnno\t1374396\tREPORTER GENE ASSAY::CHLORAMPHENICOL ACETYLTRANSFERASE (CAT) ASSAY',
    'northern': 'TF0000526\tENST00000327443\tSP1_HUMAN\tGS0001072\tENSG00000135517\t12\t56843286\t56848426\tHomo sapiens\tORegAnno\t9421492\tRNA EXPRESSION ASSAY::NORTHERN BLOT',
    'target_transcript': 'TF0000999\tENSMUST00000068264\tMEIS1\tGS0014640\tENSMUST00000111087\t-\t-\t-\tMus musculus\tTFe\t12183364\tELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::SUPERSHIFT ',
}


def parse(*keys):
    text = ''.join(LINES[key] + '\n' for key in keys)
    opener = SimpleNamespace(result={'pazar_fixture.csv': io.StringIO(text)})
    return list(iter_tf_targets(opener))


def row(key):
    return parse(key)[0]


def extract(key):
    parsed = row(key)
    record = mod.interactions_schema(parsed)
    assert record is not None
    out = SilverExtractor('pazar', 'interactions')
    out.process_record(record, parsed, 'fixture:1', 1)
    return record, out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_module_has_one_dataset():
    assert list(mod.resource.datasets()) == ['interactions']


def test_parser_names_the_columns_and_cleans_the_cells():
    parsed = row('tag')
    assert tuple(parsed) == COLUMNS
    assert parsed['project'] == 'CanEuCre MaxiPromoter'
    assert parsed['method'] == (
        'PROTEIN BINDING ASSAY::CHROMATIN IMMUNOPRECIPITATION WITH TAG '
        'SEQUENCING (CHIP-TS)'
    )
    assert row('emsa')['method'] == (
        'ELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::SUPERSHIFT'
    )
    assert row('gene_space')['gene'] == 'ENSG00000117983'
    assert row('no_pmid')['pmid'] == ''
    assert row('no_assay')['method'] == 'NA'


def test_parser_writes_flybase_identifiers_as_flybase_does():
    parsed = row('fly_upper')
    assert parsed['tf_transcript'] == 'FBgn0015664'
    assert parsed['gene'] == 'FBgn0039044'


def test_parser_reads_every_file_of_the_archive():
    text = LINES['chip'] + '\n'
    opener = SimpleNamespace(
        result={
            'a.csv': io.StringIO(text),
            'empty.csv': io.StringIO(''),
            'b.csv': io.StringIO(text + text),
        }
    )
    assert len(list(iter_tf_targets(opener))) == 3


@pytest.mark.parametrize(
    'result',
    [
        {},
        {'empty.csv': io.StringIO('')},
        {'short.csv': io.StringIO('TF0000001\tENST00000000001\tX\n')},
        None,
    ],
)
def test_parser_rejects_an_empty_archive_or_a_wrong_table(result):
    with pytest.raises(ValueError):
        list(iter_tf_targets(SimpleNamespace(result=result)))


def test_chip_line_is_binding_with_the_assay_and_project():
    record, out = extract('chip')
    relation = out.relations[0]
    # A ChIP line records binding, which shows no effect on expression.
    assert relation.predicate == 'interacts_with'
    assert annotations(relation, 'pazar:method') == [
        'PROTEIN BINDING ASSAY::CHROMATIN IMMUNOPRECIPITATION (CHIP)'
    ]
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    assert all(
        'qualifier' not in str(a.term) for a in record.annotations
    )
    assert annotations(relation, 'publications') == ['PMID:21795385']
    assert annotations(relation, 'has_evidence_of_type') == ['ECO:0000226']
    assert annotations(relation, 'supporting_data_source') == ['CTCF_Cui']
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'protein',
        Namespace.ENST,
        'ENST00000264010',
    )
    assert (target.entity_type, target.namespace, target.identifier) == (
        'gene',
        Namespace.ENSG,
        'ENSG00000230368',
    )
    assert target.taxon == '9606'


@pytest.mark.parametrize(
    'key,predicate',
    [
        # Binding assays, predictions, conserved sites and records of sites
        # with a curation label or no assay: the factor binds a site.
        ('chip', 'interacts_with'),
        ('tag', 'interacts_with'),
        ('emsa', 'interacts_with'),
        ('worm', 'interacts_with'),
        ('tf_digits', 'interacts_with'),
        ('predicted', 'interacts_with'),
        ('conserved', 'interacts_with'),
        ('curated', 'interacts_with'),
        ('no_assay', 'interacts_with'),
        ('pmid_in_assay', 'interacts_with'),
        # Reporter gene and expression assays: an effect on expression, with
        # no sign in the file.
        ('reporter', 'affects'),
        ('northern', 'affects'),
    ],
)
def test_predicate_follows_the_assay_of_the_line(key, predicate):
    record, out = extract(key)
    assert out.relations[0].predicate == predicate
    # Neither kind of line has a direction or an aspect.
    assert all('qualifier' not in str(a.term) for a in record.annotations)


def test_the_assay_is_carried_as_written_except_where_the_line_names_none():
    _, out = extract('reporter')
    assert annotations(out.relations[0], 'pazar:method') == [
        'REPORTER GENE ASSAY::CHLORAMPHENICOL ACETYLTRANSFERASE (CAT) ASSAY'
    ]
    for key in ('no_assay', 'pmid_in_assay'):
        _, out = extract(key)
        assert not annotations(out.relations[0], 'pazar:method')


def test_the_transcription_factor_has_no_species_and_keeps_its_name():
    record, out = extract('chip')
    subject = out.entities[out.relations[0].subject_entity_key]
    assert subject.taxon is None
    assert [(i.type, i.value) for i in record.subject.identifiers] == [
        (Namespace.ENST, 'ENST00000264010'),
        (Namespace.NAME, 'CTCF_HUMAN'),
    ]


def test_species_is_that_of_the_target_gene():
    # A human transcription factor on a mouse gene.
    _, out = extract('curated')
    relation = out.relations[0]
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert subject.identifier == 'ENST00000397062'
    assert subject.taxon is None
    assert (target.namespace, target.identifier, target.taxon) == (
        Namespace.ENSG,
        'ENSMUSG00000003849',
        '10090',
    )


@pytest.mark.parametrize(
    'key,evidence',
    [
        ('chip', ['ECO:0000226']),
        ('tag', ['ECO:0000226']),
        ('emsa', ['ECO:0000096']),
        ('predicted', ['ECO:0007672']),
        ('conserved', ['ECO:0000201']),
        ('worm', ['ECO:0000066']),
        ('tf_digits', ['ECO:0000096']),
        # No term, or no assay given.
        ('no_assay', []),
        ('curated', []),
        ('pmid_in_assay', []),
    ],
)
def test_assays_have_an_exact_evidence_term_or_none(key, evidence):
    _, out = extract(key)
    assert annotations(out.relations[0], 'has_evidence_of_type') == evidence


def test_unknown_assay_is_rejected():
    with pytest.raises(ValueError, match='Unknown PAZAR assay'):
        mod.interactions_schema(dict(row('chip'), method='NEW ASSAY::X'))


def test_unknown_species_is_rejected():
    with pytest.raises(ValueError, match='Unknown PAZAR species'):
        mod.interactions_schema(dict(row('chip'), species='Pan troglodytes'))


def test_flybase_and_worm_identifiers_keep_their_own_namespace():
    _, out = extract('fly_upper')
    relation = out.relations[0]
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert (subject.namespace, subject.identifier) == (
        Namespace.FLYBASE,
        'FBgn0015664',
    )
    assert (target.namespace, target.identifier, target.taxon) == (
        Namespace.FLYBASE,
        'FBgn0039044',
        '7227',
    )
    _, out = extract('worm')
    relation = out.relations[0]
    target = out.entities[relation.object_entity_key]
    assert (target.namespace, target.identifier, target.taxon) == (
        Namespace.NAME,
        'Y66H1B.4',
        '6239',
    )


def test_transcription_factor_given_by_gene_keeps_the_gene_namespace():
    _, out = extract('tf_gene')
    subject = out.entities[out.relations[0].subject_entity_key]
    assert (subject.namespace, subject.identifier) == (
        Namespace.ENSG,
        'ENSRNOG00000004642',
    )


def test_identifier_without_its_prefix_is_not_an_identifier():
    record, _ = extract('tf_digits')
    assert [(i.type, i.value) for i in record.subject.identifiers] == [
        (Namespace.NAME, 'PAX6')
    ]


def test_missing_and_malformed_publications_are_left_out():
    _, out = extract('no_pmid')
    assert not annotations(out.relations[0], 'publications')
    # Nine digits cannot be a PubMed ID.
    _, out = extract('pmid_nine_digits')
    assert not annotations(out.relations[0], 'publications')
    # Nor can a number with a leading zero.
    parsed = dict(row('pmid_in_assay'), pmid='0123456')
    out = SilverExtractor('pazar', 'interactions')
    out.process_record(mod.interactions_schema(parsed), parsed, 'fixture:1', 1)
    assert not annotations(out.relations[0], 'publications')


def test_target_column_with_a_transcript_gives_no_relation():
    assert mod.interactions_schema(row('target_transcript')) is None


def test_line_without_any_factor_identifier_gives_no_relation():
    parsed = dict(row('chip'), tf_transcript='', tf_name='')
    assert mod.interactions_schema(parsed) is None


def test_repeated_lines_stay_repeated_observations():
    parsed = parse('chip', 'chip')
    out = SilverExtractor('pazar', 'interactions')
    for number, line in enumerate(parsed):
        out.process_record(
            mod.interactions_schema(line), line, f'fixture:{number}', number
        )
    assert len(out.relations) == 2
