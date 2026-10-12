"""DIP interactions: parser, mapper and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import dip
from pypath.inputs_v2.parsers.dip import COLUMNS, iter_dip

HEADER = 'ID interactor A\tID interactor B\tAlt. ID interactor A\tAlt. ID interactor B\tAlias(es) interactor A\tAlias(es) interactor B\tInteraction detection method(s)\tPublication 1st author(s)\tPublication Identifier(s)\tTaxid interactor A\tTaxid interactor B\tInteraction type(s)\tSource database(s)\tInteraction identifier(s)\tConfidence value(s)\tProcessing Status\t\n'
# Rows copied from Hsapi20170205.txt, in this order: DIP-42E, a self
# interaction with direct interaction and no other type; DIP-172E, a mouse and
# a human protein with the obsolete "physical interaction"; DIP-398E, a record
# whose only type is the detection method "genetic interference"; DIP-10965E, a
# publication written "pubmed:rtd:9262406" and a UniProt isoform; DIP-401E, a
# RefSeq entry that holds the word "refseq"; DIP-593E, genetic interference
# together with physical interaction; DIP-98E, a participant with only a DIP
# node; DIP-188401E, a synthetic peptide (taxid -2); DIP-44789E, a high
# throughput interaction without a quality status; DIP-192027E, a UniProt
# accession written "Q13936,159", a physical association and direct interactions.
ROWS = (
    'DIP-617N|refseq:NP_000607|uniprotkb:P01730\tDIP-617N|refseq:NP_000607|uniprotkb:P01730\t-\t-\t-\t-\tMI:0114(x-ray crystallography)|MI:0031(protein cross-linking with a bifunctional reagent)\t-\tpubmed:9168119|pubmed:DIP-17838S|pubmed:9168119|pubmed:DIP-17838S\ttaxid:9606(Homo sapiens)\ttaxid:9606(Homo sapiens)\tMI:0407(direct interaction)|MI:0407(direct interaction)\tMI:0465(dip)\tDIP-42E\tdip-quality-status:core\tdip:0004(small scale)|dip:0004(small scale)\t\t-\n',
    'DIP-493N|uniprotkb:Q62311\tDIP-147N|uniprotkb:P21675\t -\t-\t-\t-\tMI:0019(coimmunoprecipitation)|MI:0404(comigration in non denaturing gel electrophoresis)|MI:0045(experimental interaction detection)\t-\tpubmed:7667268|pubmed:DIP-231S|pubmed:7667268|pubmed:DIP-231S|pubmed:7667268|pubmed:DIP-231S\ttaxid:10090(Mus musculus)\ttaxid:9606(Homo sapiens)\tMI:0218(physical interaction)|MI:0218(physical interaction)|MI:0218(physical interaction)\tMI:0465(dip)\tDIP-172E\tdip-quality-status:core\tdip:0002(small scale)|dip:0002(small scale)|dip:0002(small scale)\t\t-\n',
    'DIP-189N|refseq:NP_005179|uniprotkb:P22681\tDIP-1061N|refseq:NP_005419|uniprotkb:P15498\t -\t-\t-\t-\tMI:0254(genetic interference)\t-\tpubmed:9200440|pubmed:DIP-279S\ttaxid:9606(Homo sapiens)\ttaxid:9606(Homo sapiens)\tMI:0254(genetic interference)\tMI:0465(dip)\tDIP-398E\tdip-quality-status:core\tdip:0002(small scale)\t\t-\n',
    'DIP-6081N|refseq:NP_004299|uniprotkb:Q07960\tDIP-6082N|uniprotkb:P60953-1\t -\t-\t-\t-\tMI:0114(x-ray crystallography)\t-\tpubmed:rtd:9262406|pubmed:DIP-1612S\ttaxid:9606(Homo sapiens)\ttaxid:9606(Homo sapiens)\tMI:0407(direct interaction)\tMI:0465(dip)\tDIP-10965E\tdip-quality-status:core\tdip:0002(small scale)\t\t-\n',
    'DIP-189N|refseq:NP_005179|uniprotkb:P22681\tDIP-199N|refseq:refseq|uniprotkb:P46108\t -\t-\t-\t-\tMI:0045(experimental interaction detection)\t-\tpubmed:9178909|pubmed:DIP-185S\ttaxid:9606(Homo sapiens)\ttaxid:9606(Homo sapiens)\tMI:0218(physical interaction)\tMI:0465(dip)\tDIP-401E\tdip-quality-status:core\tdip:0002(small scale)\t\t-\n',
    'DIP-133N|uniprotkb:P23458\tDIP-57N|refseq:NP_000620|uniprotkb:P17181\t -\t-\t-\t-\tMI:0254(genetic interference)|MI:0401(biochemical)|MI:0018(two hybrid)\t-\tpubmed:9334213|pubmed:DIP-88S|pubmed:9334213|pubmed:DIP-88S|pubmed:9334213|pubmed:DIP-88S\ttaxid:9606(Homo sapiens)\ttaxid:9606(Homo sapiens)\tMI:0254(genetic interference)|MI:0218(physical interaction)|MI:0218(physical interaction)\tMI:0465(dip)\tDIP-593E\tdip-quality-status:core\tdip:0002(small scale)|dip:0002(small scale)|dip:0002(small scale)\t\t-\n',
    'DIP-879N\tDIP-1037N|uniprotkb:P07996\t -\t-\t-\t-\tMI:0045(experimental interaction detection)\t-\tpubmed:9169439|pubmed:DIP-53S\ttaxid:9606(Homo sapiens)\ttaxid:9606(Homo sapiens)\tMI:0218(physical interaction)\tMI:0465(dip)\tDIP-98E\tdip-quality-status:core\tdip:0002(small scale)\t\t-\n',
    'DIP-29390N\tDIP-199N|refseq:refseq|uniprotkb:P46108\t -\t-\t-\t-\tMI:0077(nuclear magnetic resonance)|MI:0077(nuclear magnetic resonance)\t-\tpubmed:17515907|pubmed:DIP-4548S|pubmed:17515907|pubmed:DIP-4548S\ttaxid:-2(chemical synthesis)\ttaxid:9606(Homo sapiens)\tMI:0407(direct interaction)|MI:0407(direct interaction)\tMI:0465(dip)\tDIP-188401E\tdip-quality-status:core\tdip:0002(small scale)|dip:0002(small scale)\t\t-\n',
    'DIP-6205N|refseq:NP_002257|uniprotkb:P52292\tDIP-27532N|uniprotkb:Q15653\t -\t-\t-\t-\tMI:0676(tandem affinity purification)\t-\tpubmed:14743216|pubmed:DIP-3931S\ttaxid:9606(Homo sapiens)\ttaxid:9606(Homo sapiens)\tMI:0218(physical interaction)\tMI:0465(dip)\tDIP-44789E\tdip-quality-status:\tdip:0005(high throughput)\t\t-\n',
    'DIP-60501N|uniprotkb:Q13936,159\tDIP-31794N|refseq:NP_001734|uniprotkb:P62158\t -\t-\t-\t-\tMI:0096(pull down)|MI:0226(ion exchange chromatography)|MI:0065(isothermal titration calorimetry)|MI:0114(x-ray crystallography)|MI:0065(isothermal titration calorimetry)\t-\tpubmed:16299511|pubmed:DIP-16276S|pubmed:16299511|pubmed:DIP-16276S|pubmed:16299511|pubmed:DIP-16276S|pubmed:16299511|pubmed:DIP-16276S|pubmed:16299511|pubmed:DIP-16276S\ttaxid:9606(Homo sapiens)\ttaxid:9606(Homo sapiens)\tMI:0915(physical association)|MI:0407(direct interaction)|MI:0407(direct interaction)|MI:0407(direct interaction)|MI:0407(direct interaction)\tMI:0465(dip)\tDIP-192027E\tdip-quality-status:core\tdip:0004(small scale)|dip:0004(small scale)|dip:0004(small scale)|dip:0004(small scale)|dip:0004(small scale)\t\t-\n',
)
SOURCE = HEADER + ''.join(ROWS)
(SELF, MOUSE, GENETIC, RTD, REFSEQ, MIXED, DIP_ONLY, PEPTIDE, HIGH, MALFORMED) = range(10)


def opener(text=SOURCE):
    return SimpleNamespace(result=io.StringIO(text))


def parse(text=SOURCE):
    return list(iter_dip(opener(text)))


def extract(rows):
    out = SilverExtractor('dip', 'interactions')
    for index, row in enumerate(rows):
        record = dip.resource.interactions.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term, scope='relation'):
    return [
        a['value']
        for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


def test_parser_reads_the_named_columns_and_the_field_beyond_the_header():
    rows = parse()
    assert len(rows) == 10
    assert [name for name, _ in COLUMNS][-1] == ''
    first = rows[SELF]
    assert first['id_a'] == 'DIP-617N|refseq:NP_000607|uniprotkb:P01730'
    assert first['interaction_ids'] == 'DIP-42E'
    assert first['confidence'] == 'dip-quality-status:core'
    assert (first['quality_status'], first['evidence_count']) == ('core', 2)
    # The leading space of the second alternative-identifier column is gone.
    assert rows[MOUSE]['alt_ids_a'] == '-'
    assert 'unnamed' not in first
    assert rows[HIGH]['quality_status'] is None
    assert rows[MALFORMED]['evidence_count'] == 5


def test_parser_rejects_a_changed_file():
    with pytest.raises(ValueError, match='columns'):
        parse(HEADER.replace('Taxid interactor A', 'Taxid A') + ROWS[0])
    with pytest.raises(ValueError, match='HTML'):
        parse('<html>\n<head><title>404 Not Found</title></head>\n')
    with pytest.raises(ValueError, match='no rows'):
        parse(HEADER)
    with pytest.raises(ValueError, match='beyond the header'):
        parse(HEADER + ROWS[0].rstrip('\n') + '\textra\n')
    with pytest.raises(ValueError, match='beyond the header'):
        parse(HEADER + ROWS[0].replace('\t\t-\n', '\n'))
    with pytest.raises(ValueError, match='Unknown DIP quality status'):
        parse(HEADER + ROWS[0].replace('status:core', 'status:gold'))
    with pytest.raises(ValueError, match='quality status not understood'):
        parse(HEADER + ROWS[0].replace('dip-quality-status:core', 'core'))
    # One method fewer than types and processing statuses.
    with pytest.raises(ValueError, match='evidence columns'):
        parse(HEADER + ROWS[0].replace('MI:0114(x-ray crystallography)|', ''))


def test_participants_keep_their_accessions_and_species():
    record = dip.resource.interactions.mapper(parse()[MOUSE])
    assert slot_name(record.predicate) == 'interacts_with'
    assert record.subject.identifiers == [
        (Namespace.UNIPROT, 'Q62311'),
        (Namespace.DIP, 'DIP-493N'),
    ]
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'P21675'),
        (Namespace.DIP, 'DIP-147N'),
    ]
    out = extract(parse()[MOUSE : MOUSE + 1])
    (relation,) = out.relations
    assert annotations(relation, 'in_taxon', 'subject') == ['NCBITaxon:10090']
    assert annotations(relation, 'in_taxon', 'object') == ['NCBITaxon:9606']
    # The species of each participant is kept: the pair is not filtered to human.
    assert {e.taxon for e in out.entities.values()} == {'10090', '9606'}
    assert {e.entity_type for e in out.entities.values()} == {'protein'}


def test_interaction_is_one_relation_with_its_evidence():
    out = extract(parse()[SELF : SELF + 2])
    first, second = out.relations
    assert first.predicate == 'directly_physically_interacts_with'
    assert first.upstream_id == 'DIP-42E'
    assert annotations(first, 'has_evidence_of_type') == ['MI:0114', 'MI:0031']
    # The type written in the file is kept next to the predicate.
    assert annotations(first, 'original_predicate') == ['MI:0407']
    # The placeholder "pubmed:DIP-17838S" is not a publication.
    assert annotations(first, 'publications') == ['PMID:9168119']
    assert annotations(first, 'has_confidence_level') == ['core']
    assert annotations(first, 'dip:processing_status') == ['small scale (dip:0004)']
    # The self interaction is one protein.
    assert first.subject_entity_key == first.object_entity_key
    assert [e.identifier for e in out.entities.values()][:1] == ['P01730']
    # The obsolete "physical interaction" says only that the proteins interact.
    assert second.predicate == 'interacts_with'
    assert annotations(second, 'original_predicate') == ['MI:0218']
    assert annotations(second, 'publications') == ['PMID:7667268']
    assert annotations(second, 'dip:processing_status') == ['small scale (dip:0002)']


def test_core_status_marks_the_core_rows_and_other_rows_have_none():
    out = extract(parse())
    levels = [annotations(r, 'has_confidence_level') for r in out.relations]
    assert levels == [['core']] * 8 + [[]] + [['core']]
    # High throughput is kept, with the code of the file.
    assert annotations(out.relations[HIGH], 'dip:processing_status') == [
        'high throughput (dip:0005)'
    ]
    assert annotations(out.relations[HIGH], 'has_confidence_level') == []


def test_genetic_interference_is_a_detection_method_and_gives_the_broad_predicate():
    rows = parse()
    for index in (GENETIC, MIXED):
        record = dip.resource.interactions.mapper(rows[index])
        assert slot_name(record.predicate) == 'interacts_with'
    out = extract([rows[GENETIC], rows[MIXED]])
    assert annotations(out.relations[0], 'original_predicate') == ['MI:0254']
    assert annotations(out.relations[1], 'original_predicate') == ['MI:0254', 'MI:0218']
    # A type that is neither physical nor in the broad set gives no relation.
    other = dict(rows[SELF], interaction_types='MI:0403(colocalization)|MI:0403(colocalization)')
    assert dip.exclusion_reason(other) == 'unknown interaction type'
    assert dip.resource.interactions.mapper(other) is None


def test_isoform_is_kept_and_a_rtd_reference_is_not_a_pmid():
    out = extract([parse()[RTD]])
    (relation,) = out.relations
    assert annotations(relation, 'publications') == []
    substrate = out.entities[relation.object_entity_key]
    assert (substrate.namespace, substrate.identifier) == (Namespace.UNIPROT, 'P60953-1')
    assert substrate.molecular_form['isoform_identifier'] == {
        'ns': 'uniprot',
        'id': 'P60953-1',
    }
    assert out.entities[relation.subject_entity_key].molecular_form is None


def test_malformed_and_foreign_accessions_are_left_out():
    rows = parse()
    record = dip.resource.interactions.mapper(rows[REFSEQ])
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'P46108'),
        (Namespace.DIP, 'DIP-199N'),
    ]
    record = dip.resource.interactions.mapper(rows[MALFORMED])
    assert record.subject.identifiers == [(Namespace.DIP, 'DIP-60501N')]
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'P62158'),
        (Namespace.REFSEQ_PROTEIN, 'NP_001734'),
        (Namespace.DIP, 'DIP-31794N'),
    ]
    # Direct interaction among the types wins over physical association.
    assert slot_name(record.predicate) == 'directly_physically_interacts_with'


def test_participant_with_only_a_dip_node_and_without_a_species_is_kept():
    out = extract([parse()[DIP_ONLY], parse()[PEPTIDE]])
    peptide_side = out.entities[out.relations[1].subject_entity_key]
    assert (peptide_side.namespace, peptide_side.identifier) == (Namespace.DIP, 'DIP-29390N')
    # taxid -2 is chemical synthesis, not an organism.
    assert peptide_side.taxon is None
    assert annotations(out.relations[1], 'in_taxon', 'subject') == []
    assert annotations(out.relations[1], 'in_taxon', 'object') == ['NCBITaxon:9606']
    first = out.entities[out.relations[0].subject_entity_key]
    assert (first.namespace, first.identifier) == (Namespace.DIP, 'DIP-879N')


def test_unknown_processing_status_raises_and_nothing_is_invented():
    row = dict(parse()[SELF], processing_status='dip:0009(curated)|dip:0009(curated)')
    with pytest.raises(ValueError, match='processing status'):
        dip.resource.interactions.mapper(row)
    out = extract(parse())
    assert len(out.relations) == 10
    for relation in out.relations:
        # No sign, direction or aspect is stated by the file.
        assert not annotations(relation, 'object_direction_qualifier')
        assert not annotations(relation, 'object_aspect_qualifier')
        assert not [a for a in relation.annotations if a['value'].startswith('DIP-')]
    # The source database column names DIP itself: it is not an annotation.
    assert not annotations(out.relations[0], 'supporting_data_source')
