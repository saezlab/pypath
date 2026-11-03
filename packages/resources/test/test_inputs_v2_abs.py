"""Source contract of the ABS binding site input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.abs import resource

# Lines of abs.gff: a site in an extracted promoter, a site in a chicken
# promoter, and the sites of one factor on the opposite strand in a human and a
# mouse promoter (the positions of the first have a space in the file).
LINES = (
    'extracted04\tA0001\tSRF\t262\t271\t.\t+\t.\t# HS:  TCCTTCTTTG\t\n'
    'M13631\tA0001\tSRF\t320\t329\t.\t+\t.\t# GG:  CCTTATACGG\t\n'
    'M11844\tA0016\tHNF3\t 399\t410\t.\t-\t.\t# HS:  AAGTCAATAATC\t\n'
    'M19524\tA0016\tHNF3\t 397\t408\t.\t-\t.\t# MM:  AAGTCAATAATC\t\n'
)


def rows(text=LINES):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(resource.binding_sites._raw_parser(opener))


def extract(row):
    record = resource.binding_sites.mapper(row)
    assert record is not None
    out = SilverExtractor('abs', 'binding_sites')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out.relations[0]


def terms(relation, term, scope=None):
    return [
        a['value']
        for a in relation.annotations
        if a['term'] == term and (scope is None or a['scope'] == scope)
    ]


def test_reader_decodes_the_lines():
    parsed = rows()
    assert parsed[0] == {
        'sequence': 'extracted04',
        'entry': 'A0001',
        'factor': 'SRF',
        'start': 262,
        'end': 271,
        'strand': '+',
        'species': 'HS',
        'site': 'TCCTTCTTTG',
    }
    assert parsed[2]['start'] == 399 and parsed[2]['strand'] == '-'


def test_reader_rejects_malformed_lines_and_an_empty_file():
    with pytest.raises(ValueError, match='fields'):
        rows('M13631\tA0001\tSRF\t320\t329\t.\t+\t.\t# GG:  CCTTATACGG\n')
    with pytest.raises(ValueError, match='not a binding site'):
        rows('M13631\tA0001\tSRF\t320\t329\t0.5\t+\t.\t# GG:  CCTTATACGG\t\n')
    with pytest.raises(ValueError, match='positions'):
        rows('M13631\tA0001\tSRF\t329\t320\t.\t+\t.\t# GG:  CCTTATACGG\t\n')
    with pytest.raises(ValueError, match='no rows'):
        rows('')


def test_the_relation_is_a_bare_interaction_of_the_factor_with_the_promoter():
    record, relation = extract(rows()[0])
    assert relation.predicate == 'interacts_with'
    # Binding is observed, not an effect on expression.
    assert terms(relation, 'object_direction_qualifier') == []
    assert terms(relation, 'object_aspect_qualifier') == []
    assert terms(relation, 'publications') == []


def test_the_species_is_that_of_the_promoter_and_the_factor_has_none():
    for row, taxon in zip(rows(), ['9606', '9031', '9606', '10090']):
        _, relation = extract(row)
        assert terms(relation, 'in_taxon', 'object') == [f'NCBITaxon:{taxon}']
        assert terms(relation, 'in_taxon', 'subject') == []


def test_the_factor_is_known_by_its_name_and_not_as_a_gene_symbol():
    record, _ = extract(rows()[0])
    assert {(i.type, i.value) for i in record.subject.identifiers} == {(Namespace.NAME, 'SRF')}


def test_the_promoter_is_a_genbank_record_or_known_only_to_abs():
    extracted, _ = extract(rows()[0])
    genbank, _ = extract(rows()[1])
    assert {(i.type, i.value) for i in extracted.object.identifiers} == {
        ('abs_promoter', 'extracted04')
    }
    assert {(i.type, i.value) for i in genbank.object.identifiers} == {
        (Namespace.GENBANK, 'M13631'),
        ('abs_promoter', 'M13631'),
    }


def test_an_unknown_species_code_raises():
    with pytest.raises(ValueError, match='Unknown ABS species'):
        resource.binding_sites.mapper({**rows()[0], 'species': 'XX'})
