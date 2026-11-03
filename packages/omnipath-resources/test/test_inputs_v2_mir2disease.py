"""Source-specific Biolink contract for the miR2Disease input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import mir2disease as mod
from pypath.inputs_v2.parsers.mir2disease import iter_targets

# The start of miRtar.txt and lines from further down (ISO-8859-1, CRLF). The
# title of the MET line is broken over two lines in the source, and the title
# of the first LATS2 line is closed by a line with only a quote.
PREAMBLE = (
    'The following miRNA-target pairs, validated by luciferase reporter '
    'experiments, were extracted manually from the corresponding reference. '
    'A few target genes not documented here may be available in "TarBase" and '
    '"miRecords" database system.\r\n'
    '\r\n'
    'miRNA\tValidated target\tPub Date\tReference\r\n'
)
LINES = [
    # 0
    "hsa-let-7a\tKRAS\t2008\tA SNP in a let-7 microRNA complementary site in "
    "the KRAS 3' untranslated region increases non-small cell lung cancer risk.",
    # 1: a protein name instead of a gene symbol.
    'hsa-let-7a\tcaspase-3\t2008\tLet-7a microRNA suppresses therapeutics-induced '
    'cancer cell death by targeting caspase-3.',
    # 2, 3: spaces after the miRNA name.
    'hsa-let-7d \tBDNF\t2009\tmicroRNAs miR-124, let-7d and miR-181a regulate '
    'Cocaine-induced Plasticity',
    'hsa-let-7d \tD3R\t2009\tmicroRNAs miR-124, let-7d and miR-181a regulate '
    'Cocaine-induced Plasticity',
    # 4: a space after the target name.
    'hsa-miR-106a\tIL-10 \t2009\tPosttranscriptional regulation of interleukin-10 '
    'expression by hsa-miR-106a',
    # 5, 6: the same target from two miRNA forms.
    'hsa-miR-17-3p\tT\xa6\xc2RII\t2007\tSynergistic action of the microRNA-17 '
    'polycistron and Myc in aggressive cancer development.',
    'hsa-miR-17-5p\tT\xa6\xc2RII\t2007\tSynergistic action of the microRNA-17 '
    'polycistron and Myc in aggressive cancer development.',
    # 7, 8: a title broken over two lines.
    'hsa-miR-206\tMET\t2009\t"The muscle-specific microRNA miR-206 blocks human '
    'rhabdomyosarcoma growth in',
    'xenotransplanted mice by promoting myogenic differentiation."',
    # 9, 10: a title closed on a line of its own.
    'hsa-miR-373\tLATS2\t2009\t"MicroRNA-373 (miR-373) post-transcriptionally '
    'regulates large tumor suppressor, homolog 2 (LATS2) and stimulates '
    'proliferation in human esophageal cancer.',
    '"',
    'hsa-miR-373\tLATS2\t2006\tA genetic screen implicates miRNA-372 and '
    'miRNA-373 as oncogenes in testicular germ cell tumors.',
]
SOURCE = PREAMBLE + '\r\n'.join(LINES) + '\r\n'


def rows(text=SOURCE):
    return list(iter_targets(SimpleNamespace(result=io.StringIO(text))))


def extract(row):
    record = mod.resource.mirna_targets.mapper(row)
    assert record is not None
    out = SilverExtractor('mir2disease', 'mirna_targets')
    out.process_record(record, row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        record,
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_parser_skips_the_preamble_and_names_the_fields():
    parsed = rows()
    assert len(parsed) == 10
    assert parsed[0] == {
        'mirna': 'hsa-let-7a',
        'mirna_species': 'hsa',
        'target': 'KRAS',
        'year': '2008',
        'reference': (
            "A SNP in a let-7 microRNA complementary site in the KRAS 3' "
            'untranslated region increases non-small cell lung cancer risk.'
        ),
    }


def test_parser_strips_spaces_around_names():
    parsed = rows()
    assert parsed[2]['mirna'] == 'hsa-let-7d'
    assert parsed[4]['target'] == 'IL-10'


def test_parser_joins_a_title_broken_over_lines():
    parsed = rows()
    assert parsed[7]['target'] == 'MET'
    assert parsed[7]['reference'] == (
        'The muscle-specific microRNA miR-206 blocks human rhabdomyosarcoma '
        'growth in xenotransplanted mice by promoting myogenic differentiation.'
    )
    assert parsed[8]['reference'].startswith('MicroRNA-373 (miR-373)')
    assert parsed[8]['reference'].endswith('human esophageal cancer.')
    assert parsed[9]['year'] == '2006'


def test_parser_rejects_changed_header_and_empty_table():
    with pytest.raises(ValueError, match='columns'):
        rows(PREAMBLE.replace('Validated target', 'Target') + LINES[0] + '\r\n')
    with pytest.raises(ValueError, match='no rows'):
        rows(PREAMBLE)
    with pytest.raises(ValueError, match='year'):
        rows(PREAMBLE + 'hsa-let-7a\tKRAS\tlast year\tA title\r\n')


def test_pair_is_a_human_mirna_affecting_a_target_without_a_direction():
    record, relation, subject, target = extract(rows()[0])
    assert relation.predicate == 'affects'
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'microrna',
        Namespace.NAME,
        'hsa-let-7a',
    )
    assert subject.taxon == '9606'
    assert (target.entity_type, target.namespace, target.identifier) == (
        'gene',
        Namespace.NAME,
        'KRAS',
    )
    # The table names no organism for the target.
    assert not target.taxon
    # No direction, no aspect, and no publication: the table gives titles only.
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    assert not annotations(relation, 'publications')


def test_target_name_is_not_taken_for_a_gene_symbol():
    _, _, _, target = extract(rows()[1])
    assert target.identifier == 'caspase-3'
    assert all(i['ns'] != Namespace.GENESYMBOL for i in target.identifiers)


def test_mirna_forms_stay_separate():
    _, _, first, first_target = extract(rows()[5])
    _, _, second, second_target = extract(rows()[6])
    assert {first.identifier, second.identifier} == {
        'hsa-miR-17-3p',
        'hsa-miR-17-5p',
    }
    assert first_target.identifier == second_target.identifier == 'T\xa6\xc2RII'


def test_unknown_organism_code_raises():
    row = dict(rows()[0], mirna='mmu-miR-1', mirna_species='mmu')
    with pytest.raises(ValueError, match='organism code'):
        mod.resource.mirna_targets.mapper(row)


def test_no_relation_without_a_target_name():
    row = dict(rows()[0], target='')
    assert mod.resource.mirna_targets.mapper(row) is None
