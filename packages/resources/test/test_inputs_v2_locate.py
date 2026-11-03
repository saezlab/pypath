"""LOCATE: streaming reader, mapper and extractor on records of the release 6 human file."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import locate
from pypath.inputs_v2.locate import CLASS_NAMESPACE, ECO_OF_GO_EVIDENCE_CODE, resource
from pypath.inputs_v2.parsers.locate import iter_localizations, iter_protein_classes


def protein(organism, cls, source_id, source_name, accn):
    return (
        f'<protein><organism>{organism}</organism><class>{cls}</class><synonyms />'
        f'<source source_id="{source_id}"><source_name>{source_name}</source_name>'
        f'<accn>{accn}</accn></source><protein_sequence>MEGK</protein_sequence></protein>'
    )


def xref(source_id, source_name, accn):
    return (
        f'<xref><source source_id="{source_id}"><source_name>{source_name}</source_name>'
        f'<accn>{accn}</accn></source></xref>'
    )


def external(goid, label, source_id, source_name, accn, evidence=''):
    return (
        '<reference><pmid>0</pmid><source />'
        f'<evidence>{evidence}</evidence><locations><location goid="{goid}"><tier1>{label}</tier1>'
        f'</location></locations><source source_id="{source_id}"><source_name>{source_name}'
        f'</source_name><accn>{accn}</accn></source></reference>'
    )


def prediction(source_id, method, location, goid, evaluation):
    goid_element = f'<goid>{goid}</goid>' if goid else ''
    return (
        f'<source source_id="{source_id}"><method>{method}</method><location>{location}'
        f'</location>{goid_element}<evaluation>{evaluation}</evaluation></source>'
    )


def literature(organism, goid, tiers, notes, pmid):
    labels = ''.join(f'<tier{i}>{t}</tier{i}>' for i, t in enumerate(tiers, 1))
    return (
        f'<reference><author>x</author><title>x</title><citation>x</citation>'
        f'<organism>{organism}</organism><locations><location goid="{goid}">{labels}</location>'
        f'</locations><notes>{notes}</notes><date_analyzed>2006-07-07</date_analyzed>'
        f'<source source_id="3"><source_name>PubMed</source_name><accn>{pmid}</accn></source></reference>'
    )


# Records of LOCATE_human_v6_20081121.xml, with the sequence, transcript, motif and memo
# elements left out (the reader does not read them).
CALNEXIN = (
    '<LOCATE_protein uid="6000002">'
    + protein('Human', 'typeI', 23, 'RefSeq Protein', 'NP_001737')
    + '<externalannot>'
    + external('GO:0005764', 'Lysosomes', 54, 'HPRD', '00252')
    + external('GO:0042470', 'melanosome', 51, 'UniProt/SPTrEMBL', 'P27824')
    + '</externalannot><scl_prediction>'
    + prediction(33, 'CELLO', 'endoplasmic reticulum', 'GO:0005783', '2.64')
    + prediction(38, 'MultiLoc', 'plasma membrane', 'GO:0005886', '0.01')
    # A method that gives the location a score of 0 does not predict it.
    + prediction(35, 'Proteome_Analyst', 'extracellular region', 'GO:0005576', '0.00')
    + '</scl_prediction><xrefs>'
    + xref(6, 'GNF SymAtlas', 821)
    + xref(17, 'Entrez Gene', 821)
    + xref(45, 'HGNC', 1473)
    + xref(47, 'UniProtKB-SwissProt', 'P27824')
    + xref(51, 'UniProt/SPTrEMBL', 'CALX_HUMAN')
    + xref(49, 'Ensembl-Peptide Human', 'ENSP00000247461')
    + xref(21, 'Entrez Protein', 'AAH03552')
    + '</xrefs></LOCATE_protein>'
)
# A curated localization with two tiers, whose goid lists two terms.
NUMA1 = (
    '<LOCATE_protein uid="6000029">'
    + protein('Human', 'cytoplasmic', 21, 'Entrez Protein', 'AAI03766')
    + '<literature>'
    + literature('Human', 'GO:0005634;GO:0005730', ['Nucleus', 'Nucleolus'], '', 16381850)
    + '</literature><xrefs>'
    + xref(23, 'RefSeq Protein', 'NP_006176')
    + xref(47, 'UniProtKB-SwissProt', 'Q14980')
    + xref(49, 'Ensembl-Peptide Human', 'ENSP00000352675')
    + xref(51, 'UniProt/SPTrEMBL', 'NUMA1_HUMAN')
    + '</xrefs></LOCATE_protein>'
)
# An external annotation with a GO evidence code (IDA), from Ensembl.
SDHA = (
    '<LOCATE_protein uid="6004326">'
    + protein('Human', 'cytoplasmic', 49, 'Ensembl-Peptide Human', 'ENSP00000332417')
    + '<externalannot>'
    + external('GO:0005739', 'mitochondrion', 49, 'Ensembl-Peptide Human', 'ENSP00000332417', 'IDA')
    + external('GO:0005739', 'mitochondrion', 54, 'HPRD', '02914')
    + '</externalannot><xrefs>'
    + xref(17, 'Entrez Gene', 6389)
    + xref(23, 'RefSeq Protein', 'NP_004159')
    + xref(51, 'UniProt/SPTrEMBL', 'Q05BB6')
    + '</xrefs></LOCATE_protein>'
)
# A reference with no GO term, and one whose GO terms are in the other order than its names.
BECN1 = (
    '<LOCATE_protein uid="6011390">'
    + protein('Human', 'cytoplasmic', 21, 'Entrez Protein', 'AAH10276')
    + '<literature>'
    + literature('Human', '', ['Golgi Apparatus', 'Golgi Trans Face'], '', 11306555)
    + literature('Human', 'GO:0005634;GO:0005737', ['Cytoplasm', 'Nucleus'], '', 11309306)
    + '</literature></LOCATE_protein>'
)
# "No prediction" has no GO term; WoLFPSORT gives a combined location with two terms.
ENSP = (
    '<LOCATE_protein uid="6000062">'
    + protein('Human', 'cytoplasmic', 49, 'Ensembl-Peptide Human', 'ENSP00000371780')
    + '<scl_prediction>'
    + prediction(35, 'Proteome_Analyst', 'No prediction', '', '0.00')
    + prediction(37, 'WoLFPSORT', 'cytoplasm_nucleus', 'GO:0005737;GO:0005634', '7.50')
    + '</scl_prediction><xrefs>'
    + xref(17, 'Entrez Gene', 728409)
    + '</xrefs></LOCATE_protein>'
)
# A human protein with a curated localization from a mouse study, with notes.
LAMA = (
    '<LOCATE_protein uid="6000698">'
    + protein('Human', 'secretome', 21, 'Entrez Protein', 'AAH66552')
    + '<literature>'
    + literature(
        'Mouse', 'GO:0005576', ['Extracellular'],
        'Extracellular; found in the basement membranes (major component).', 9151674,
    )
    + '</literature><xrefs>'
    + xref(23, 'RefSeq Protein', 'NP_001098677')
    + xref(23, 'RefSeq Protein', 'NP_002281')
    + xref(49, 'Ensembl-Peptide Human', 'ENSP00000374114')
    + '</xrefs></LOCATE_protein>'
)
RECORDS = [CALNEXIN, NUMA1, SDHA, BECN1, ENSP, LAMA]
# A record of LOCATE_mouse_v6_20081121.xml with no UniProt accession and an MGI id
# (written as digits in the file), with an MGI symbol that the reader does not use.
UPK3A = (
    '<LOCATE_protein uid="6000004">'
    + protein('Mouse', 'typeI', 21, 'Entrez Protein', 'AAF34681')
    + '<literature>'
    + literature('Cow', 'GO:0016324;GO:0005886', ['Plasma Membrane', 'Apical Plasma Membrane'], '', 12475947)
    + '</literature><scl_prediction>'
    + prediction(33, 'CELLO', 'extracellular region', 'GO:0005576', '3.29')
    + '</scl_prediction><xrefs>'
    + xref(4, 'MGI ID', 98914)
    + xref(5, 'MGI Symbol', 'Upk3a')
    + xref(23, 'RefSeq Protein', 'NP_075967')
    + xref(49, 'Ensembl-Peptide Mouse', 'ENSMUSP00000023070')
    + '</xrefs></LOCATE_protein>'
)


def document(*records):
    return (
        '<?xml version="1.0"?>\n<LOCATE_doc xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        + ''.join(records)
        + '</LOCATE_doc>'
    )


def opener(*records):
    return SimpleNamespace(result={'LOCATE.xml': io.BytesIO(document(*records).encode())})


def rows(*records):
    return list(iter_localizations(opener(*(records or RECORDS))))


def process(dataset, row):
    out = SilverExtractor('locate', dataset)
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


def of_kind(parsed, kind):
    return [r for r in parsed if r['kind'] == kind]


def test_import_defines_the_datasets_and_the_two_sources_without_downloading():
    assert set(resource.datasets()) == {'terms', 'human', 'mouse'}
    assert locate.human_download.url.startswith('https://rescued.omnipathdb.org/')
    assert 'web.archive.org' in locate.mouse_download.url
    assert locate.config.license.name == 'NAR_OPEN_ACCESS'


def test_an_external_annotation_names_the_database_and_the_record_not_the_empty_source():
    parsed = of_kind(rows(CALNEXIN), 'external annotation')
    assert [(r['go'], r['database'], r['record'], r['evidence_code']) for r in parsed] == [
        ('GO:0005764', 'HPRD', '00252', None),
        ('GO:0042470', 'UniProt/SPTrEMBL', 'P27824', None),
    ]
    _, relation, term = relation_of(process('human', parsed[0]))
    assert relation.predicate == 'located_in'
    assert (term.entity_type, term.namespace, term.identifier) == ('ontology_class', 'go', 'GO:0005764')
    assert attributes(relation, 'locate:evidence_type') == ['external annotation']
    assert attributes(relation, 'supporting_data_source') == ['HPRD']
    assert attributes(relation, 'locate:source_record') == ['HPRD:00252']


def test_a_go_evidence_code_is_kept_as_written():
    (ida, hprd) = of_kind(rows(SDHA), 'external annotation')
    assert ida['evidence_code'] == 'IDA' and hprd['evidence_code'] is None
    _, relation, _ = relation_of(process('human', ida))
    assert attributes(relation, 'locate:evidence_code') == ['IDA']
    assert attributes(relation, 'supporting_data_source') == ['Ensembl-Peptide Human']
    # The ECO term of the code, from the official GAF-ECO mapping.
    assert attributes(relation, 'has_evidence_of_type') == ['ECO:0000314']
    _, relation, _ = relation_of(process('human', hprd))
    assert attributes(relation, 'has_evidence_of_type') == []


def test_every_go_evidence_code_has_its_eco_term_and_an_unknown_code_raises():
    # The seven codes of the human file.
    assert {
        code: ECO_OF_GO_EVIDENCE_CODE[code] for code in ('IEA', 'TAS', 'IDA', 'NAS', 'ISS', 'IEP', 'IC')
    } == {
        'IEA': 'ECO:0000501', 'TAS': 'ECO:0000304', 'IDA': 'ECO:0000314', 'NAS': 'ECO:0000303',
        'ISS': 'ECO:0000250', 'IEP': 'ECO:0000270', 'IC': 'ECO:0000305',
    }
    (row, *_) = of_kind(rows(SDHA), 'external annotation')
    with pytest.raises(ValueError, match='Unknown LOCATE evidence code'):
        resource.human.mapper({**row, 'evidence_code': 'XYZ'})


def test_a_literature_reference_with_two_go_terms_is_two_statements_with_the_publication():
    parsed = of_kind(rows(NUMA1), 'literature')
    assert [(r['go'], r['pmid'], r['study_organism']) for r in parsed] == [
        ('GO:0005634', '16381850', 'Human'), ('GO:0005730', '16381850', 'Human'),
    ]
    _, relation, term = relation_of(process('human', parsed[1]))
    assert relation.predicate == 'located_in'
    assert term.identifier == 'GO:0005730'
    assert attributes(relation, 'locate:evidence_type') == ['literature']
    assert attributes(relation, 'publications') == ['PMID:16381850']
    assert attributes(relation, 'locate:study_organism') == ['Human']


def test_notes_and_the_organism_of_the_study_are_kept_and_differ_from_the_protein_s():
    (row,) = of_kind(rows(LAMA), 'literature')
    protein_, relation, _ = relation_of(process('human', row))
    assert protein_.taxon == '9606'
    assert attributes(relation, 'locate:study_organism') == ['Mouse']
    assert attributes(relation, 'description') == [
        'Extracellular; found in the basement membranes (major component).'
    ]


def test_a_reference_without_a_go_term_has_no_statement_and_no_tier_is_assigned():
    parsed = of_kind(rows(BECN1), 'literature')
    # The Golgi reference has an empty goid; the second reference lists the GO
    # terms in the other order than its names (Cytoplasm, Nucleus).
    assert [(r['go'], r['pmid']) for r in parsed] == [
        ('GO:0005634', '11309306'), ('GO:0005737', '11309306'),
    ]
    assert not [a for r in parsed for a in r if a.startswith('tier')]


def test_a_prediction_has_the_method_label_and_score_and_a_zero_score_is_none():
    parsed = of_kind(rows(CALNEXIN), 'prediction')
    assert [(r['method'], r['go'], r['score']) for r in parsed] == [
        ('CELLO', 'GO:0005783', '2.64'), ('MultiLoc', 'GO:0005886', '0.01'),
    ]
    _, relation, term = relation_of(process('human', parsed[1]))
    assert term.identifier == 'GO:0005886'
    assert attributes(relation, 'locate:evidence_type') == ['prediction']
    assert attributes(relation, 'locate:prediction_method') == ['MultiLoc']
    assert attributes(relation, 'locate:predicted_location') == ['plasma membrane']
    (score,) = [a for a in relation.annotations if a['term'] == 'locate:prediction_score']
    # A number, with the method whose scale it is on, not the text of the file.
    assert score['quantity']['has_numeric_value'] == pytest.approx(0.01)
    assert score['value'] != '0.01'
    # Predictions have no confidence score: the scales differ by method.
    assert attributes(relation, 'has_confidence_score') == []


def test_no_prediction_is_none_and_a_combined_location_is_one_statement_per_term():
    parsed = of_kind(rows(ENSP), 'prediction')
    assert [(r['go'], r['location']) for r in parsed] == [
        ('GO:0005737', 'cytoplasm_nucleus'), ('GO:0005634', 'cytoplasm_nucleus'),
    ]
    assert all(r['method'] == 'WoLFPSORT' and r['score'] == '7.50' for r in parsed)


def test_protein_is_identified_by_the_accessions_of_the_file_and_not_translated():
    parsed = rows()
    ids = {r['ids']['locate']: r['ids'] for r in parsed}
    assert ids['6000002'] == {
        'locate': '6000002', 'uniprot': ['P27824'], 'uniprot_entry': ['CALX_HUMAN'],
        'entrez': ['821'], 'ensp': ['ENSP00000247461'], 'refseq_protein': ['NP_001737'],
        'genbank_protein': [], 'hgnc': ['1473'], 'mgi': [],
    }
    # The SPTrEMBL value of NUMA1 is an entry name; the GenBank accession is the protein's own.
    assert ids['6000029']['uniprot'] == ['Q14980']
    assert ids['6000029']['genbank_protein'] == ['AAI03766']
    # A protein with no UniProt or Entrez accession keeps what it has.
    assert ids['6011390'] == {
        'locate': '6011390', 'uniprot': [], 'uniprot_entry': [], 'entrez': [], 'ensp': [],
        'refseq_protein': [], 'genbank_protein': ['AAH10276'], 'hgnc': [], 'mgi': [],
    }
    assert ids['6004326']['uniprot'] == ['Q05BB6']
    # The GenBank accessions of the xrefs are left out.
    assert ids['6000002']['genbank_protein'] == []


def test_subject_has_the_identifiers_and_the_taxon():
    (row, *_) = of_kind(rows(CALNEXIN), 'prediction')
    subject, _, _ = relation_of(process('human', row))
    assert (subject.entity_type, subject.namespace, subject.identifier, subject.taxon) == (
        'protein', 'uniprot', 'P27824', '9606'
    )
    assert {(i['ns'], i['id']) for i in subject.identifiers} == {
        (Namespace.UNIPROT, 'P27824'), (Namespace.UNIPROT_ENTRY, 'CALX_HUMAN'),
        (Namespace.ENTREZ, '821'), (Namespace.ENSP, 'ENSP00000247461'),
        (Namespace.REFSEQ_PROTEIN, 'NP_001737'), ('locate', '6000002'),
        (Namespace.HGNC, '1473'),
    }
    # The HGNC id is a gene-level identifier of the record: after all the others,
    # LOCATE's own included, so that it does not become the leading identifier.
    mapped = resource.human.mapper(row)
    assert [i.type for i in mapped.subject.identifiers][-2:] == ['locate', Namespace.HGNC]
    assert mapped.subject.identifiers[0].type == Namespace.UNIPROT


def test_the_class_of_a_protein_is_an_association_to_a_term_of_locate():
    (row,) = of_kind(rows(CALNEXIN), 'protein class')
    subject, relation, term = relation_of(process('human', row))
    assert relation.predicate == 'associated_with'
    assert (term.entity_type, term.namespace, term.identifier) == (
        'ontology_class', CLASS_NAMESPACE, 'typeI'
    )
    assert subject.identifier == 'P27824'
    assert attributes(relation, 'locate:evidence_type') == []


def test_classes_are_terms_with_names_and_the_objects_match_them():
    terms = list(iter_protein_classes(opener(*RECORDS)))
    assert [t['name'] for t in terms] == ['typeI', 'cytoplasmic', 'secretome']
    defined = {}
    for term in terms:
        out = SilverExtractor('locate', 'terms')
        out.process_record(resource.terms.mapper(term), term, 'fixture:1', 1)
        (entity,) = out.entities.values()
        assert {(i['ns'], i['id']) for i in entity.identifiers} == {
            (CLASS_NAMESPACE, term['name']), (Namespace.NAME, term['name']),
        }
        defined[entity.identifier] = (entity.entity_type, entity.namespace, entity.identifier)
    for row in of_kind(rows(), 'protein class'):
        _, _, term = relation_of(process('human', row))
        assert defined[term.identifier] == (term.entity_type, term.namespace, term.identifier)


def test_go_terms_are_not_defined_here_and_no_location_name_is_made_a_term():
    for row in rows():
        if row['kind'] == 'protein class':
            continue
        _, _, term = relation_of(process('human', row))
        assert term.namespace == Namespace.GO


def test_the_images_of_the_experiments_are_not_a_source_of_statements():
    # The release has empty experimental_data elements: no experimental kind exists.
    assert {r['kind'] for r in rows()} == {
        'literature', 'external annotation', 'prediction', 'protein class',
    }


@pytest.mark.parametrize(
    'record',
    [
        NUMA1.replace('GO:0005730', 'GO:5730'),  # not a GO identifier
        CALNEXIN.replace('<evaluation>2.64', '<evaluation>high'),  # a score that is no number
        NUMA1.replace('<organism>Human</organism><class>', '<organism></organism><class>', 1),
        CALNEXIN.replace('P27824</accn></source></xref><xref><source source_id="51">', '$$</accn></source></xref><xref><source source_id="51">'),
    ],
)
def test_reader_raises_on_what_it_does_not_know(record):
    with pytest.raises(ValueError):
        rows(record)


def test_unknown_values_raise_in_the_mapper():
    (row, *_) = of_kind(rows(CALNEXIN), 'prediction')
    for change in ({'method': 'Other'}, {'organism': 'Rat'}):
        with pytest.raises(ValueError, match='Unknown LOCATE'):
            resource.human.mapper({**row, **change})
    (row, *_) = of_kind(rows(CALNEXIN), 'external annotation')
    with pytest.raises(ValueError, match='Unknown LOCATE database'):
        resource.human.mapper({**row, 'database': 'Elsewhere'})
    (row,) = of_kind(rows(CALNEXIN), 'protein class')
    with pytest.raises(ValueError, match='Unknown LOCATE class'):
        resource.human.mapper({**row, 'class': 'typeIII'})


def test_reader_raises_on_a_file_that_is_not_locate_and_on_an_empty_one():
    with pytest.raises(ValueError, match='Not a LOCATE file'):
        list(iter_localizations(SimpleNamespace(result={'x': io.BytesIO(b'<other></other>')})))
    with pytest.raises(ValueError, match='no protein'):
        list(iter_localizations(SimpleNamespace(result={'x': io.BytesIO(document().encode())})))


def test_the_mgi_id_is_a_secondary_identifier_written_with_its_prefix():
    parsed = rows(UPK3A)
    assert {r['ids']['mgi'] and r['ids']['mgi'][0] for r in parsed} == {'MGI:98914'}
    assert all(r['ids']['uniprot'] == [] and r['ids']['hgnc'] == [] for r in parsed)
    # The MGI symbol is not an identifier of the record.
    assert all(r['ids']['genbank_protein'] == ['AAF34681'] for r in parsed)
    (row, *_) = of_kind(parsed, 'prediction')
    mapped = resource.mouse.mapper(row)
    # The MGI id comes after all the others, LOCATE's own included.
    assert [(i.type, i.value) for i in mapped.subject.identifiers] == [
        (Namespace.ENSP, 'ENSMUSP00000023070'), (Namespace.REFSEQ_PROTEIN, 'NP_075967'),
        (Namespace.GENBANK, 'AAF34681'), ('locate', '6000004'), (Namespace.MGI, 'MGI:98914'),
    ]
    subject, _, _ = relation_of(process('mouse', row))
    assert (Namespace.MGI, 'MGI:98914') in {(i['ns'], i['id']) for i in subject.identifiers}
    assert subject.taxon == '10090'


def test_an_hgnc_or_mgi_id_that_is_not_digits_raises():
    broken = CALNEXIN.replace(xref(45, 'HGNC', 1473), xref(45, 'HGNC', 'HGNC:1473'))
    assert broken != CALNEXIN
    with pytest.raises(ValueError, match='HGNC id'):
        rows(broken)
    with pytest.raises(ValueError, match='MGI ID id'):
        rows(UPK3A.replace(xref(4, 'MGI ID', 98914), xref(4, 'MGI ID', 'x98914')))


def test_a_protein_with_two_hgnc_ids_keeps_both_and_a_repeat_once():
    two = CALNEXIN.replace(
        xref(45, 'HGNC', 1473), xref(45, 'HGNC', 1473) + xref(45, 'HGNC', 1473) + xref(45, 'HGNC', 99999)
    )
    (row, *_) = rows(two)
    assert row['ids']['hgnc'] == ['1473', '99999']
