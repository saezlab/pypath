"""MPPI: the PSI-MI XML level 1 reader, the row reader and the mapper on source records."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.mppi import accession, resource
from pypath.inputs_v2.parsers.base import iter_psimi_xml
from pypath.inputs_v2.parsers.mppi import iter_interactions

HEADER = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<entrySet xmlns="net:sf:psidev:mi" version="1" level="1"><entry>'
    '<source><names><shortLabel>MIPS</shortLabel>'
    '<fullName>Munich Information Center for Protein Sequences</fullName></names>'
    '<bibref><xref><primaryRef db="Pubmed" id="9399795" /></xref></bibref></source>'
)
FOOTER = '</entry></entrySet>'


def reference(db, accession_id):
    """A ``primaryRef``; the source leaves the ``id`` out of some of them."""
    identifier = '' if accession_id is None else f' id="{accession_id}"'
    return f'<primaryRef db="{db}"{identifier} />'


def interactor(label, name, db, accession_id, taxon):
    return (
        '<proteinParticipant><proteinInteractor>'
        f'<names><shortLabel>{label}</shortLabel><fullName>{name}</fullName></names>'
        f'<xref>{reference(db, accession_id)}</xref>'
        f'<organism ncbiTaxId="{taxon}" /></proteinInteractor>'
        '<confidence unit="boolean" value="1" />'
        '<isTaggedProtein>false</isTaggedProtein>'
        '<isOverexpressedProtein>false</isOverexpressedProtein>'
        '</proteinParticipant>'
    )


def interaction(pubmed, method, method_id, first, second):
    """An ``<interaction>`` in the layout of the MPPI file."""
    return (
        '<interaction><experimentList><experimentDescription>'
        f'<bibref><xref>{reference("PubMed", pubmed)}</xref></bibref>'
        f'<interactionDetection><names><shortLabel>{method}</shortLabel></names>'
        f'<xref><primaryRef db="MI" id="{method_id}" /></xref></interactionDetection>'
        '</experimentDescription></experimentList>'
        f'<participantList>{interactor(*first)}{interactor(*second)}</participantList>'
        '</interaction>'
    )


def document(*interactions):
    return f'{HEADER}<interactionList>{"".join(interactions)}</interactionList>{FOOTER}'


# Interaction 1 of the file.
S100 = interaction(
    '9867828', 'two hybrid', '0018',
    ('1', 'S100A8;CAGA;MRP8; calgranulin A (migration inhibitory factor-related protein 8)',
     'SP', 'P27005', '10090'),
    ('2', 'S100A9;CAGB;MRP14; calgranulin B (migration inhibitory factor-related protein 14)',
     'SP', 'P31725', '10090'),
)
# Interaction 1747: two PubMed IDs in one cell, a name whose gene names are not a list.
GRIP1 = interaction(
    '12458226; 9069286', 'coip: coimmunoprecipitation', '0019',
    ('944', 'GRIP1; Glutamate receptor-interacting protein 1', 'SP', 'Q925T6', '10090'),
    ('945', 'GRIA2 GRIA2; Glutamate receptor 2', 'SP', 'P23819', '10090'),
)
# Interaction 255: a homodimer; gene names separated by ``or``.
HNF4 = interaction(
    '9343418', 'coip: coimmunoprecipitation', '0019',
    *2 * [(
        '165',
        'HNF4A or NR2A1 or TCF14 or HNF4 or HNF-4; hepatocyte nuclear factor 4-alpha; '
        'transcription factor HNF-4',
        'SP', 'P22449', '10116',
    )],
)
# Interaction 5: both proteins are named by GenBank nucleotide accessions.
TRP3 = interaction(
    '10611319', 'copurification', '0025',
    ('5', 'TRP3; calcium influx channel protein', 'GB', 'U47050', '9606'),
    ('6', 'IP3R-3; inositol 1,4,5-trisphosphate receptor, type 3', 'GB', 'D26351', '9606'),
)
# Interaction 742: the first protein has taxon 0; a mouse and a human protein.
CAV1 = interaction(
    '9647788', 'coip: coimmunoprecipitation', '0019',
    ('499', 'CAV1; Caveolin-1', 'SP', 'Q03135', '0'),
    ('389', 'ARHA;ARH12;RHOA;RHO12; transforming protein RhoA', 'SP', 'P06749', '9606'),
)
# Interaction 600: the second reference names ERK1 or another protein, and was cut.
ERK = interaction(
    '10194762', 'copurification', '0025',
    ('425', 'STAT5A;STAT5; signal transducer and activator of transcription 5A',
     'SP', 'P42229', '9606'),
    ('426',
     'MAPK3;PRKM3;ERK1 or MAPK1;PRKM1;ERK2; extracellular signal-regulated kinase 1 or 2',
     'SP', 'P27361 (ERK1) or sp', '9606'),
)
# Interaction 20: no PubMed ID, and a protein of an unnamed mammal without accession.
FOS = interaction(
    None, 'far western blotting', '0047',
    ('33', 'c-fos; fos cellular oncogene', 'SP', 'P12841', '10116'),
    ('35', 'TFIIF-30; general transcription factor TFIIE-34', 'TREMBL', None, '40674'),
)
# Interaction 72: the accession is a dash.
CD100 = interaction(
    '8955171', 'coip: coimmunoprecipitation', '0019',
    *2 * [('57', '; surface glycoprotein CD100', 'TREMBL', '-', '9606')],
)
# Interaction 107: a name without gene names.
CEBP = interaction(
    '8946919', 'coip: coimmunoprecipitation', '0019',
    ('86', '; CCAAT/enhancer-binding protein (C/EBP)', 'SP', 'P28033', '10090'),
    ('85', 'RB1;RB-1; retinoblastoma-associated protein (Rb protein)', 'SP', 'P13405', '10090'),
)


def opener(*interactions):
    return SimpleNamespace(result=io.BytesIO(document(*interactions).encode('utf-8')))


def rows(*interactions):
    return list(iter_interactions(opener(*interactions)))


def extract(source):
    (row,) = rows(source)
    record = resource.datasets()['interactions'].mapper(row)
    assert record is not None
    out = SilverExtractor('mppi', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def annotations(relation, scope='relation'):
    return sorted(
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def test_reader_yields_the_source_values_of_an_interaction():
    (record,) = iter_psimi_xml(document(S100))
    assert record['index'] == 1
    assert record['source']['label'] == 'MIPS'
    assert record['source']['bibliography'][0]['id'] == '9399795'
    (experiment,) = record['experiments']
    assert experiment['bibliography'] == [
        {'primary': True, 'db': 'PubMed', 'id': '9867828', 'secondary': None, 'version': None}
    ]
    assert experiment['detection_method']['label'] == 'two hybrid'
    assert experiment['detection_method']['mi'] == 'MI:0018'
    first, second = record['participants']
    assert first['interactor']['label'] == '1'
    assert first['interactor']['xrefs'][0]['db'] == 'SP'
    assert first['interactor']['xrefs'][0]['id'] == 'P27005'
    assert first['interactor']['organism'] == {'taxon': '10090', 'label': None, 'name': None}
    assert first['confidence'] == {'unit': 'boolean', 'value': '1'}
    assert (first['is_tagged'], first['is_overexpressed']) == (False, False)
    assert second['interactor']['xrefs'][0]['id'] == 'P31725'
    # Elements the file does not have are absent, not defaulted.
    assert first['role'] is None
    assert record['interaction_types'] == []
    assert experiment['host_organisms'] == []


def test_reader_resolves_references_to_the_entry_lists():
    compact = (
        '<entrySet xmlns="net:sf:psidev:mi" version="1" level="1"><entry>'
        '<experimentList><experimentDescription id="E1">'
        '<bibref><xref><primaryRef db="pubmed" id="11283351"/></xref></bibref>'
        '<hostOrganism ncbiTaxId="4932"><names><shortLabel>yeast</shortLabel></names>'
        '</hostOrganism>'
        '<interactionDetection><names><shortLabel>two hybrid</shortLabel></names>'
        '<xref><primaryRef db="psi-mi" id="MI:0018"/></xref></interactionDetection>'
        '</experimentDescription></experimentList>'
        '<interactorList>'
        '<proteinInteractor id="P1"><names><shortLabel>a</shortLabel></names>'
        '<xref><primaryRef db="uniprot" id="P12345"/>'
        '<secondaryRef db="sgd" id="S0001"/></xref></proteinInteractor>'
        '<proteinInteractor id="P2"><names><shortLabel>b</shortLabel></names>'
        '</proteinInteractor></interactorList>'
        '<interactionList><interaction><experimentList><experimentRef ref="E1"/>'
        '</experimentList><participantList>'
        '<proteinParticipant><proteinInteractorRef ref="P1"/><role>bait</role>'
        '</proteinParticipant>'
        '<proteinParticipant><proteinInteractorRef ref="P2"/><role>prey</role>'
        '</proteinParticipant>'
        '</participantList></interaction></interactionList></entry></entrySet>'
    )
    (record,) = iter_psimi_xml(compact)
    (experiment,) = record['experiments']
    assert experiment['id'] == 'E1'
    assert experiment['detection_method']['mi'] == 'MI:0018'
    assert experiment['host_organisms'] == [
        {'taxon': '4932', 'label': 'yeast', 'name': None}
    ]
    first, second = record['participants']
    assert [ref['id'] for ref in first['interactor']['xrefs']] == ['P12345', 'S0001']
    assert [ref['primary'] for ref in first['interactor']['xrefs']] == [True, False]
    assert (first['role'], second['role']) == ('bait', 'prey')
    assert second['interactor']['xrefs'] == []
    with pytest.raises(ValueError, match="unknown interactor 'P3'"):
        list(iter_psimi_xml(compact.replace('ref="P2"', 'ref="P3"')))


def test_reader_rejects_other_documents():
    with pytest.raises(ValueError, match='level 1'):
        list(iter_psimi_xml(document(S100).replace('level="1"', 'level="2"'), level='1'))
    with pytest.raises(ValueError, match='entrySet'):
        list(iter_psimi_xml('<html><body>Not found</body></html>'))


def test_rows_keep_the_references_as_written():
    (row,) = rows(ERK)
    assert row['b_db'] == 'SP'
    assert row['b_id'] == 'P27361 (ERK1) or sp'
    assert row['pubmed'] == '10194762'
    assert (row['method_label'], row['method_mi']) == ('copurification', 'MI:0025')
    (row,) = rows(FOS)
    assert row['pubmed'] is None
    assert row['b_id'] is None
    assert row['b_taxon'] == '40674'


def test_row_reader_rejects_an_empty_or_changed_download():
    with pytest.raises(ValueError, match='empty'):
        list(iter_interactions(SimpleNamespace(result=io.BytesIO(b''))))
    with pytest.raises(ValueError, match='no interactions'):
        list(iter_interactions(opener()))
    two_experiments = S100.replace(
        '</experimentDescription>',
        '</experimentDescription><experimentDescription></experimentDescription>',
    )
    with pytest.raises(ValueError, match='2 experiments'):
        rows(two_experiments)


def test_interaction_is_a_physical_interaction_with_method_and_publication():
    relation, first, second = extract(S100)
    assert relation.predicate == 'interacts_with'
    assert identifiers(first) == {
        (Namespace.UNIPROT, 'P27005'),
        (Namespace.NAME, 'calgranulin A (migration inhibitory factor-related protein 8)'),
        (Namespace.SYNONYM, 'S100A8'),
        (Namespace.SYNONYM, 'CAGA'),
        (Namespace.SYNONYM, 'MRP8'),
    }
    assert (second.namespace, second.identifier) == ('uniprot', 'P31725')
    assert {first.entity_type, second.entity_type} == {'protein'}
    assert {first.taxon, second.taxon} == {'10090'}
    assert annotations(relation) == [
        ('has_evidence_of_type', 'MI:0018'),
        ('publications', 'PMID:9867828'),
    ]
    # No direction, effect, score or role is asserted.
    assert {a['term'] for a in relation.annotations} == {
        'has_evidence_of_type',
        'in_taxon',
        'publications',
    }


def test_every_pubmed_id_of_a_cell_is_a_publication():
    relation, _, second = extract(GRIP1)
    assert annotations(relation) == [
        ('has_evidence_of_type', 'MI:0019'),
        ('publications', 'PMID:12458226'),
        ('publications', 'PMID:9069286'),
    ]
    # ``GRIA2 GRIA2`` is not a list of gene names; only the protein name is kept.
    assert identifiers(second) == {
        (Namespace.UNIPROT, 'P23819'),
        (Namespace.NAME, 'Glutamate receptor 2'),
    }


def test_gene_names_separated_by_or_are_synonyms():
    relation, first, second = extract(HNF4)
    assert relation.subject_entity_key == relation.object_entity_key
    assert identifiers(first) == {
        (Namespace.UNIPROT, 'P22449'),
        (Namespace.NAME, 'hepatocyte nuclear factor 4-alpha; transcription factor HNF-4'),
        (Namespace.SYNONYM, 'HNF4A'),
        (Namespace.SYNONYM, 'NR2A1'),
        (Namespace.SYNONYM, 'TCF14'),
        (Namespace.SYNONYM, 'HNF4'),
        (Namespace.SYNONYM, 'HNF-4'),
    }


def test_name_without_gene_names_gives_no_synonym():
    _, first, _ = extract(CEBP)
    assert identifiers(first) == {
        (Namespace.UNIPROT, 'P28033'),
        (Namespace.NAME, 'CCAAT/enhancer-binding protein (C/EBP)'),
    }


def test_genbank_accessions_keep_their_namespace():
    _, first, second = extract(TRP3)
    assert (Namespace.GENBANK, 'U47050') in identifiers(first)
    assert (Namespace.GENBANK, 'D26351') in identifiers(second)
    assert not {ns for ns, _ in identifiers(first) | identifiers(second)} & {Namespace.UNIPROT}


def test_species_come_from_each_participant_and_taxon_zero_is_no_species():
    relation, first, second = extract(CAV1)
    assert first.taxon is None
    assert second.taxon == '9606'
    assert annotations(relation, 'subject') == []
    assert annotations(relation, 'object') == [('in_taxon', 'NCBITaxon:9606')]


def test_cut_reference_naming_alternatives_is_not_an_accession():
    (row,) = rows(ERK)
    assert accession(row, 'a') == (Namespace.UNIPROT, 'P42229')
    # The leading accession is one of two proteins the source could not tell apart.
    assert accession(row, 'b') is None
    assert resource.datasets()['interactions'].mapper(row) is None


@pytest.mark.parametrize('source', [FOS, CD100])
def test_no_relation_for_a_participant_without_accession(source):
    (row,) = rows(source)
    assert resource.datasets()['interactions'].mapper(row) is None


def test_accession_must_belong_to_the_named_database():
    (row,) = rows(S100)
    assert accession(dict(row, a_db='GB', a_id='Q9JLM9'), 'a') is None
    assert accession(dict(row, a_db='EMBL'), 'a') is None
    # The repaired forms are UniProt's: they do not hold for another database.
    assert accession(dict(row, a_db='GB', a_id='P29354_2'), 'a') is None
    assert accession(dict(row, a_db='GB', a_id='DLK_MOUSE'), 'a') is None


# Interaction 311: the first protein is written as an isoform, ``P29354_2``.
GRB2_ISOFORM = interaction(
    '9299436', 'two hybrid', '0018',
    ('197', 'GRB2 or ASH; growth factor receptor-bound protein 2 (isoform GRB3-3); '
     'GRB2 adapter protein (isoform GRB3-3); SH2/SH3 adapter GRB2 (isoform GRB3-3)',
     'SP', 'P29354_2', '9606'),
    ('198', 'ADA; Adenosine deaminase(ADA)', 'SP', 'P00813', '9606'),
)
# Interaction 1697: a homodimer written with a UniProt entry name.
DLK = interaction(
    '11689349', 'two hybrid', '0018',
    *2 * [('917', 'DLK1 or DLK or PREF1 or SCP-1; Delta-like protein, splice variant',
           'SP', 'DLK_MOUSE', '10090')],
)


def test_isoform_written_with_an_underscore_is_the_uniprot_isoform():
    (row,) = rows(GRB2_ISOFORM)
    assert accession(row, 'a') == (Namespace.UNIPROT, 'P29354-2')
    assert accession(row, 'b') == (Namespace.UNIPROT, 'P00813')
    relation, first, second = extract(GRB2_ISOFORM)
    assert (first.namespace, first.identifier) == ('uniprot', 'P29354-2')
    assert first.molecular_form['isoform_identifier'] == {'ns': 'uniprot', 'id': 'P29354-2'}
    assert (second.namespace, second.identifier) == ('uniprot', 'P00813')


def test_entry_name_is_a_uniprot_entry_identifier():
    (row,) = rows(DLK)
    assert accession(row, 'a') == (Namespace.UNIPROT_ENTRY, 'DLK_MOUSE')
    relation, first, _ = extract(DLK)
    assert (Namespace.UNIPROT_ENTRY, 'DLK_MOUSE') in identifiers(first)
    assert relation.subject_entity_key == relation.object_entity_key


@pytest.mark.parametrize(
    'value',
    ['o54880', 'P04901/sp', '11167', 'Q14803 or sp', 'P16092 (mouse) or sp', 'P29354_', 'P29354_x'],
)
def test_other_forms_of_a_reference_are_still_not_repaired(value):
    (row,) = rows(S100)
    assert accession(dict(row, a_id=value, a_db='TREMBL'), 'a') is None
    assert accession(dict(row, a_id=value), 'a') is None


def reversed_interaction(source):
    start = source.index('<participantList>') + len('<participantList>')
    end = source.index('</participantList>')
    body = source[start:end]
    split = body.index('<proteinParticipant>', 1)
    return source[:start] + body[split:] + body[:split] + source[end:]


def mapped(*sources):
    data = rows(*sources)
    mapper = resource.datasets()['interactions'].mapper
    return [row for row in data if mapper(row) is not None], data


def test_second_row_of_a_pair_listed_in_both_orders_is_a_mirror():
    # Interactions 1 and 2 of the file: the same experiment, proteins exchanged.
    kept, data = mapped(S100, reversed_interaction(S100))
    assert [row['mirror'] for row in data] == [False, True]
    assert [row['a_id'] for row in kept] == ['P27005']
    out = SilverExtractor('mppi', 'interactions')
    for index, row in enumerate(data):
        record = resource.datasets()['interactions'].mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    assert len(out.relations) == 1


def test_rows_that_are_not_mirrors_are_kept():
    other_method = S100.replace('id="0018"', 'id="0019"')
    other_paper = S100.replace('9867828', '1')
    different_taxon = reversed_interaction(S100).replace('10090', '9606', 1)
    # Same direction twice (interactions 28 and 29 of the file), another method,
    # another publication, and a participant whose own fields differ.
    kept, data = mapped(
        S100, S100, reversed_interaction(other_method), reversed_interaction(other_paper),
        different_taxon,
    )
    assert [row['mirror'] for row in data] == [False] * 5
    assert len(kept) == 5


def test_row_with_the_same_protein_twice_is_never_a_mirror():
    _, data = mapped(HNF4, HNF4)
    assert [row['mirror'] for row in data] == [False, False]


def test_malformed_pubmed_id_and_unknown_method_are_left_out():
    (row,) = rows(S100)
    mapper = resource.datasets()['interactions'].mapper
    out = SilverExtractor('mppi', 'interactions')
    out.process_record(
        mapper(dict(row, pubmed='n/a;123456789', method_mi=None)), row, 'fixture:1', 1
    )
    assert annotations(out.relations[0]) == []
