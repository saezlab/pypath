"""Source-specific Biolink contract for the miRTarBase input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import mirtarbase as mod
from pypath.inputs_v2.parsers.mirtarbase import COLUMNS, iter_interactions

HEADER = ','.join(COLUMNS)

# Lines copied from miRTarBase_SE_WR.csv and miRTarBase_WE_Other.csv, which
# start with a byte order mark and end lines with CRLF.
LINES = [
    # 0: a functional block, human miRNA and gene.
    'MIRT000002,hsa-miR-20a-5p,hsa,HIF1A,3091.0,hsa,'
    'Luciferase reporter assay//Western blot//Northern blot,Functional MTI,18632605.0',
    # 1: a negative result for the same pair and publication.
    'MIRT000002,hsa-miR-20a-5p,hsa,HIF1A,3091.0,hsa,'
    'Luciferase reporter assay//Western blot//Northern blot,Non-Functional MTI,18632605.0',
    # 2: weak support.
    'MIRT003133,mmu-miR-122-5p,mmu,Aldoa,11674.0,mmu,qRT-PCR,Functional MTI (Weak),16459310.0',
    # 3: weak negative result.
    'MIRT001961,hsa-miR-155-5p,hsa,RHOA,387.0,hsa,qRT-PCR,Non-Functional MTI (Weak),18794355.0',
    # 4: a viral miRNA with a human target.
    'MIRT003063,kshv-miR-K12-5-3p,kshv,BCLAF1,9774.0,hsa,Luciferase reporter assay,'
    'Functional MTI,19098914.0',
    # 5: the miRNA name is human, the organism column is mouse.
    'MIRT001190,hsa-miR-21-5p,mmu,Pten,19211.0,mmu,qRT-PCR//Western blot,'
    'Functional MTI,25990308.0',
    # 6: a space before the miRNA name.
    'MIRT755383, hsa-miR-320b,hsa,MYC,4609.0,hsa,Western blotting//qRT-PCR,'
    'Functional MTI,38260139.0',
    # 7: the publication ID is 0.
    'MIRT004724,hsa-miR-124-3p,hsa,VIM,7431.0,hsa,qRT-PCR//Western blotting,'
    'Functional MTI,0.0',
    # 8: a virus as the target organism.
    'MIRT004803,mmu-miR-24-3p,mmu,VSIVgp5,1489835.0,vsv,Luciferase reporter assay,'
    'Functional MTI,17613256.0',
]
SOURCE = '﻿' + '\r\n'.join([HEADER, *LINES]) + '\r\n'


def rows(text=SOURCE):
    return list(iter_interactions(SimpleNamespace(result=io.StringIO(text))))


def extract(row):
    record = mod.map_interaction(row)
    assert record is not None
    out = SilverExtractor('mirtarbase', 'strong_evidence')
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


def test_module_has_one_dataset_per_evidence_table():
    assert list(mod.resource.datasets()) == [
        'strong_evidence',
        'weak_evidence_other',
    ]


def test_parser_names_the_columns_and_reads_the_float_ids():
    parsed = rows()
    assert len(parsed) == 9
    assert parsed[0] == {
        'mirtarbase_id': 'MIRT000002',
        'mirna': 'hsa-miR-20a-5p',
        'mirna_species': 'hsa',
        'target_gene': 'HIF1A',
        'entrez_id': '3091',
        'target_species': 'hsa',
        'experiments': [
            'Luciferase reporter assay',
            'Western blot',
            'Northern blot',
        ],
        'support_type': 'Functional MTI',
        'pmid': '18632605',
    }


def test_parser_strips_spaces_around_the_mirna_name():
    assert rows()[6]['mirna'] == 'hsa-miR-320b'


def test_parser_keeps_a_zero_publication_id_as_written():
    assert rows()[7]['pmid'] == '0'


def test_parser_rejects_changed_columns_and_empty_tables():
    changed = SOURCE.replace('Support Type', 'Support')
    with pytest.raises(ValueError, match='columns'):
        rows(changed)
    with pytest.raises(ValueError, match='no rows'):
        rows(HEADER + '\r\n')
    with pytest.raises(ValueError, match='not a number'):
        rows(HEADER + '\r\nMIRT1,a,hsa,B,x,hsa,qRT-PCR,Functional MTI,1\r\n')


def test_functional_block_is_a_mirna_affecting_a_gene_without_a_direction():
    parsed = rows()
    record, relation, subject, target = extract(parsed[0])
    assert relation.predicate == 'affects'
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'microrna',
        Namespace.NAME,
        'hsa-miR-20a-5p',
    )
    # The Entrez Gene ID is the primary key of the gene; the symbol follows.
    assert (target.entity_type, target.namespace, target.identifier) == (
        'gene',
        Namespace.ENTREZ,
        '3091',
    )
    assert subject.taxon == target.taxon == '9606'
    assert {
        (i['ns'], i['id']) for i in target.identifiers
    } >= {(Namespace.GENESYMBOL, 'HIF1A'), (Namespace.ENTREZ, '3091')}
    assert annotations(relation, 'publications') == ['PMID:18632605']
    assert annotations(relation, 'has_confidence_level') == ['Functional MTI']
    # The tables state no effect: no direction, no aspect, no repression.
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')


def test_weak_support_keeps_the_source_label():
    _, relation, _, _ = extract(rows()[2])
    assert annotations(relation, 'has_confidence_level') == [
        'Functional MTI (Weak)'
    ]


@pytest.mark.parametrize('index', [1, 3])
def test_non_functional_block_produces_no_relation(index):
    parsed = rows()
    assert parsed[index]['support_type'].startswith('Non-Functional')
    assert mod.map_interaction(parsed[index]) is None


def test_unknown_support_type_raises():
    row = dict(rows()[0], support_type='Probably functional')
    with pytest.raises(ValueError, match='support type'):
        mod.map_interaction(row)


def test_organism_of_each_side_comes_from_its_own_column():
    parsed = rows()
    _, _, subject, target = extract(parsed[4])
    assert subject.taxon == '37296'
    assert target.taxon == '9606'
    # The name says human, the organism column says mouse; the column wins.
    _, _, subject, target = extract(parsed[5])
    assert subject.identifier == 'hsa-miR-21-5p'
    assert subject.taxon == target.taxon == '10090'
    _, _, subject, target = extract(parsed[8])
    assert (subject.taxon, target.taxon) == ('10090', '11277')


def test_unknown_organism_code_raises():
    row = dict(rows()[0], target_species='xyz')
    with pytest.raises(ValueError, match='organism code'):
        mod.map_interaction(row)


def test_publication_id_that_is_not_a_pubmed_id_is_left_out():
    _, relation, _, _ = extract(rows()[7])
    assert not annotations(relation, 'publications')


@pytest.mark.parametrize('number', ['0', '1', '0123456', '123456789'])
def test_placeholder_leading_zero_and_nine_digit_numbers_are_not_publications(
    number,
):
    # The file writes 0 and 1 where a block has no publication; this module
    # keeps the seven or eight digits of the files' intact IDs.
    _, relation, _, _ = extract(dict(rows()[0], pmid=number))
    assert not annotations(relation, 'publications')


def test_mirna_name_is_not_translated_or_merged():
    names = set()
    for row in (rows()[0], rows()[6]):
        _, _, subject, _ = extract(row)
        names.add(subject.identifier)
        assert subject.namespace == Namespace.NAME
    assert names == {'hsa-miR-20a-5p', 'hsa-miR-320b'}


def test_no_relation_without_a_target_gene():
    row = dict(rows()[0], target_gene='', entrez_id=None)
    assert mod.map_interaction(row) is None
