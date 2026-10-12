"""Source-specific Biolink contract for the HTRIdb input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import htri as mod
from pypath.inputs_v2.parsers import htri as parser

# Lines copied from HTRIdb.csv, header included.
HEADER = 'OID;GENEID_TF;SYMBOL_TF;GENEID_TG;SYMBOL_TG;TECHNIQUE;PUBMED_ID;\n'
LINES = (
    '28627;1044;CDX1;1824;DSC2;Chromatin Immunoprecipitation;18819935;\n'
    '28627;1044;CDX1;1824;DSC2;Electrophoretic Mobility Shift Assay;18819935;\n'
    '40521;367;AR;13;AADAC;Chromatin Immunoprecipitation coupled with microarray;20610535;\n'
    '28612;1869;E2F1;652;BMP4;CpG Chromatin Immunoprecipitation;12629508;\n'
    '29384;466;ATF1;8370;HIST2H4A;Dnase I Footprinting;7779811;\n'
    '32330;1869;E2F1;1869;E2F1;Chromatin Immunoprecipitation;14576826;\n'
)


def rows(text=HEADER + LINES):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(mod.resource.interactions._raw_parser(opener))


def extract(row):
    record = mod.resource.interactions.mapper(row)
    assert record is not None
    out = SilverExtractor('htri', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_reader_keys_lines_by_the_header_and_drops_the_trailing_column():
    parsed = rows()
    assert len(parsed) == 6
    assert parsed[0] == {
        'OID': '28627',
        'GENEID_TF': '1044',
        'SYMBOL_TF': 'CDX1',
        'GENEID_TG': '1824',
        'SYMBOL_TG': 'DSC2',
        'TECHNIQUE': 'Chromatin Immunoprecipitation',
        'PUBMED_ID': '18819935',
    }


def test_reader_rejects_an_empty_file():
    with pytest.raises(ValueError, match='no rows'):
        rows(HEADER)


def test_reader_rejects_a_missing_column():
    with pytest.raises(ValueError, match='TECHNIQUE'):
        rows('OID;GENEID_TF;SYMBOL_TF;GENEID_TG;SYMBOL_TG;PUBMED_ID;\n1;2;A;3;B;1234567;\n')


def test_reader_rejects_a_line_with_another_number_of_fields():
    with pytest.raises(ValueError, match='fields'):
        rows(HEADER + '28627;1044;CDX1;1824;DSC2;Chromatin Immunoprecipitation;\n')


@pytest.mark.parametrize(
    'pmid,publications',
    [
        ('18819935', ['PMID:18819935']),
        # The shared rule: 1 to 8 digits, no leading zero.
        ('953', ['PMID:953']),
        ('0123456', []),
        ('123456789', []),
    ],
)
def test_pubmed_identifier_follows_the_shared_rule(pmid, publications):
    _, out = extract(dict(rows()[0], PUBMED_ID=pmid))
    assert annotations(out.relations[0], 'publications') == publications


def test_binding_observation_is_interacts_with_and_has_no_effect():
    _, out = extract(rows()[0])
    relation = out.relations[0]
    assert relation.predicate == 'interacts_with'
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    assert annotations(relation, 'publications') == ['PMID:18819935']
    assert annotations(relation, 'has_evidence_of_type') == ['MI:0402']
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    # The extractor chooses the canonical identifier; the symbol and the
    # Entrez Gene ID of the source are both kept.
    def ids(entity):
        return {(i['ns'], i['id']) for i in entity.identifiers}

    assert subject.entity_type == 'protein'
    assert ids(subject) == {
        (Namespace.GENESYMBOL, 'CDX1'),
        (Namespace.ENTREZ, '1044'),
    }
    assert target.entity_type == 'gene'
    assert ids(target) == {
        (Namespace.GENESYMBOL, 'DSC2'),
        (Namespace.ENTREZ, '1824'),
    }
    assert subject.taxon == target.taxon == '9606'


@pytest.mark.parametrize(
    'index,term',
    [(0, 'MI:0402'), (1, 'MI:0413'), (2, 'MI:0225'), (4, 'MI:0606')],
)
def test_technique_with_an_exact_term_is_an_evidence_type(index, term):
    _, out = extract(rows()[index])
    assert annotations(out.relations[0], 'has_evidence_of_type') == [term]


def test_technique_without_an_exact_term_has_no_evidence_type():
    _, out = extract(rows()[3])
    assert not annotations(out.relations[0], 'has_evidence_of_type')


def test_each_observation_is_its_own_relation():
    out = SilverExtractor('htri', 'interactions')
    for number, row in enumerate(rows()[:2], 1):
        out.process_record(
            mod.resource.interactions.mapper(row), row, f'fixture:{number}', number
        )
    assert len(out.relations) == 2
    assert len(out.entities) == 2


def test_self_regulation_keeps_protein_and_gene_apart():
    _, out = extract(rows()[5])
    relation = out.relations[0]
    assert relation.subject_entity_key != relation.object_entity_key
    assert {e.entity_type for e in out.entities.values()} == {'protein', 'gene'}


def test_unknown_technique_is_rejected():
    with pytest.raises(ValueError, match='Unknown HTRIdb technique'):
        mod.resource.interactions.mapper(dict(rows()[0], TECHNIQUE='Western'))


@pytest.mark.parametrize('missing', ['SYMBOL_TF', 'SYMBOL_TG'])
def test_row_without_a_symbol_and_id_emits_no_relation(missing):
    role = missing.split('_')[1]
    row = dict(rows()[0], **{missing: '', f'GENEID_{role}': ''})
    assert mod.resource.interactions.mapper(row) is None


def test_parser_module_exposes_the_columns():
    assert parser.COLUMNS[0] == 'OID'
