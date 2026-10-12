"""HIPPIE: reader, mapper and extractor on lines of hippie_current.txt (release 2.4)."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import hippie
from pypath.inputs_v2.parsers.hippie import iter_hippie, parse_details

# Lines copied from hippie_current.txt: a pair with several techniques,
# publications and source databases; one with the placeholder PubMed ID
# 2147483647; one with entry names of two UniProt entries of a gene; one found in
# another species; one from STRING with no technique; one with no entry name
# for the first protein; and the line without a name for the second protein, with
# a score of 0 and no evidence.
LINES = (
    'AL1A1_HUMAN\t216\tAL1A1_HUMAN\t216\t0.76\t'
    'experiments:in vivo,Two-hybrid;pmids:12081471,16189514,25416956;'
    'sources:HPRD,BioGRID,IntAct,MINT,I2D,Rual05\n'
    'CBP_HUMAN\t1387\tJUN_HUMAN\t3725\t0.9\t'
    'experiments:in vitro,in vivo,Reconstituted Complex;'
    'pmids:9786917,2138276,2147483647;sources:HPRD,BioGRID,I2D,Lim06\n'
    'HSP71_HUMAN,HS71B_HUMAN\t3303\tAPOB_HUMAN\t338\t0.62\t'
    'experiments:in vitro;pmids:9252351,16713569;sources:HPRD,Lim06\n'
    'LAP2_HUMAN\t55914\tERBB2_HUMAN\t2064\t0.97\t'
    'experiments:in vitro,x-ray crystallography;pmids:11278603,10878805;'
    'species:Mus musculus (Mouse);sources:HPRD,BioGRID,HomoMINT,I2D\n'
    'AKT1_HUMAN\t207\tCADH5_HUMAN\t1003\t0.49\tpmids:39232006;sources:STRING,IntAct\n'
    '\t64283\tRHOA_HUMAN\t387\t0.84\t'
    'experiments:in vitro,Two-hybrid;pmids:9199174,32203420;sources:HPRD,BioGRID\n'
    'LYN_HUMAN\t4067\t\t100124696\t0.0\t\n'
)


def read(text=LINES):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(hippie.resource.interactions._raw_parser(opener))


def extract(rows):
    out = SilverExtractor('hippie', 'interactions')
    for index, row in enumerate(rows):
        record = hippie.resource.interactions.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def values(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def identifiers(entity):
    return [(i.type, i.value) for i in entity.identifiers]


def test_module_import_defines_the_dataset_without_downloading():
    assert list(hippie.resource.datasets()) == ['interactions']
    download = hippie.resource.interactions.download
    assert download.url.endswith('/hippie/hippie_current.txt')
    assert hippie.config.id.name == 'HIPPIE'


def test_reader_unpacks_the_columns_and_the_evidence_cell():
    rows = read()
    assert len(rows) == 7
    assert rows[0] == {
        'entry_names_a': ['AL1A1_HUMAN'],
        'entrez_a': '216',
        'entry_names_b': ['AL1A1_HUMAN'],
        'entrez_b': '216',
        'score': '0.76',
        'experiments': ['in vivo', 'Two-hybrid'],
        'pmids': ['12081471', '16189514', '25416956'],
        'sources': ['HPRD', 'BioGRID', 'IntAct', 'MINT', 'I2D', 'Rual05'],
        'species': [],
    }
    assert rows[2]['entry_names_a'] == ['HSP71_HUMAN', 'HS71B_HUMAN']
    assert rows[3]['species'] == ['Mus musculus (Mouse)']
    assert rows[5]['entry_names_a'] == []
    assert rows[6]['entry_names_b'] == []
    assert rows[6]['experiments'] == rows[6]['pmids'] == rows[6]['sources'] == []


def test_reader_rejects_a_layout_it_does_not_know():
    with pytest.raises(ValueError, match='unexpected evidence item'):
        parse_details('techniques:in vivo')
    with pytest.raises(ValueError, match='unexpected evidence item'):
        parse_details('pmids:1;pmids:2')
    with pytest.raises(ValueError, match='fields'):
        read('A_HUMAN\t1\tB_HUMAN\t2\t0.5\n')
    with pytest.raises(ValueError, match='Entrez'):
        read('A_HUMAN\tx\tB_HUMAN\t2\t0.5\tpmids:1\n')
    with pytest.raises(ValueError, match='score'):
        read('A_HUMAN\t1\tB_HUMAN\t2\t1.5\tpmids:1\n')
    with pytest.raises(ValueError, match='entry name'):
        read('P12345\t1\tB_HUMAN\t2\t0.5\tpmids:1\n')
    with pytest.raises(ValueError, match='no rows'):
        list(iter_hippie(SimpleNamespace(result=io.StringIO('\n'))))


def test_pair_is_a_plain_interaction_with_score_publications_and_sources():
    out = extract(read()[:1])
    (relation,) = out.relations
    # The table mixes direct and indirect evidence: only the broad predicate.
    assert relation.predicate == 'interacts_with'
    assert values(relation, 'publications') == [
        'PMID:12081471', 'PMID:16189514', 'PMID:25416956',
    ]
    (score,) = [
        a['quantity'] for a in relation.annotations
        if a['term'] == 'has_confidence_score'
    ]
    assert score['has_numeric_value'] == pytest.approx(0.76)
    # Every database and study the row names, under the table's own labels.
    assert set(values(relation, 'supporting_data_source')) == {
        'BioGRID', 'HPRD', 'I2D', 'IntAct', 'MINT', 'Rual05',
    }
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}


def test_participant_has_the_entry_names_and_the_entrez_gene():
    rows = read()
    record = hippie.resource.interactions.mapper(rows[2])
    assert identifiers(record.subject) == [
        (Namespace.UNIPROT_ENTRY, 'HSP71_HUMAN'),
        (Namespace.UNIPROT_ENTRY, 'HS71B_HUMAN'),
        (Namespace.ENTREZ, '3303'),
    ]
    # A protein without an entry name keeps the Entrez Gene ID alone.
    record = hippie.resource.interactions.mapper(rows[5])
    assert identifiers(record.subject) == [(Namespace.ENTREZ, '64283')]
    assert identifiers(record.object) == [
        (Namespace.UNIPROT_ENTRY, 'RHOA_HUMAN'),
        (Namespace.ENTREZ, '387'),
    ]


def test_placeholder_pubmed_id_is_not_a_publication():
    (relation,) = extract(read()[1:2]).relations
    assert values(relation, 'publications') == ['PMID:9786917', 'PMID:2138276']


def test_techniques_are_not_mapped_to_evidence_types():
    # The techniques are free text of the source, not PSI-MI terms.
    (relation,) = extract(read()[:1]).relations
    terms = {a['term'] for a in relation.annotations}
    assert 'has_evidence_of_type' not in terms
    assert 'original_predicate' not in terms


def test_species_is_kept_as_the_source_writes_it():
    out = extract(read()[3:5])
    other_species, string_pair = out.relations
    assert values(other_species, 'hippie:conserved_in_species') == ['Mus musculus (Mouse)']
    assert values(string_pair, 'hippie:conserved_in_species') == []


def test_every_named_source_is_carried_and_an_unknown_one_raises():
    (relation,) = extract(read()[2:3]).relations
    assert set(values(relation, 'supporting_data_source')) == {'HPRD', 'Lim06'}
    (string_pair,) = extract(read()[4:5]).relations
    assert set(values(string_pair, 'supporting_data_source')) == {'IntAct', 'STRING'}
    row = dict(read()[0], sources=['NotADatabase'])
    with pytest.raises(ValueError, match='Unknown HIPPIE source database'):
        hippie.resource.interactions.mapper(row)


def test_low_scores_are_kept_and_no_threshold_is_applied():
    out = extract(read())
    assert len(out.relations) == 7
    scores = sorted(
        a['quantity']['has_numeric_value']
        for r in out.relations
        for a in r.annotations
        if a['term'] == 'has_confidence_score'
    )
    assert scores[:2] == [0.0, 0.49]


def test_row_without_evidence_gives_a_relation_without_publications_or_sources():
    (relation,) = extract(read()[6:7]).relations
    assert relation.predicate == 'interacts_with'
    assert values(relation, 'publications') == []
    assert values(relation, 'supporting_data_source') == []
