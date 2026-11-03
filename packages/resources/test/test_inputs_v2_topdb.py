"""TOPDB: reader of the XML, mappers and extractor on entries of the real file."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import topdb
from pypath.inputs_v2.parsers import topdb as topdb_parser

NS = 'xmlns="https://topdb.unitmp.org"'
# Entries of topdb_all.xml with the sequence, the cross-references and the
# experiments left out; the elements the reader uses are as in the file, with
# all the regions (the reader checks Numtm against the Membrane regions).
ADA2C = f"""\
<TOPDB {NS} ID="ADA2C_HUMAN" type="Alpha_Polytopic">
  <Name>Alpha-2C adrenergic receptor</Name>
  <CrossRef><PDB ID="6kuw"><Resolution>2.8</Resolution></PDB></CrossRef>
  <Sequence Len="462" Checksum="11d31e46a0fcd7f1b057f3401151a7a2"><Seq>MASPALAAAL</Seq></Sequence>
  <Membranes>
    <Membrane>Cell membrane;Cytoplasm;Extracellular</Membrane>
  </Membranes>
  <Topology>
    <Reliability>70.6</Reliability>
    <Numtm Count="6"/>
    <Regions>
      <Region Begin="1" End="53" Loc="Outside"/>
      <Region Begin="54" End="75" Loc="Membrane"/>
      <Region Begin="76" End="87" Loc="Inside"/>
      <Region Begin="88" End="105" Loc="Membrane"/>
      <Region Begin="106" End="125" Loc="Outside"/>
      <Region Begin="126" End="147" Loc="Membrane"/>
      <Region Begin="148" End="167" Loc="Inside"/>
      <Region Begin="168" End="185" Loc="Membrane"/>
      <Region Begin="186" End="210" Loc="Outside"/>
      <Region Begin="211" End="231" Loc="Membrane"/>
      <Region Begin="232" End="322" Loc="Inside"/>
      <Region Begin="323" End="334" Loc="Membrane"/>
      <Region Begin="335" End="462" Loc="Outside"/>
    </Regions>
  </Topology>
  <Experiments/>
</TOPDB>"""
# Two membranes (the file writes them in one or two Membranes elements).
ANTR2 = f"""\
<TOPDB {NS} ID="ANTR2_HUMAN" type="Alpha_Bitopic">
  <Name>Anthrax toxin receptor 2 </Name>
  <Membranes>
    <Membrane>Cell membrane;Cytoplasm;Extracellular</Membrane>
    <Membrane>Endoplasmic reticulum membrane;Cytoplasm;Cysternal space</Membrane>
  </Membranes>
  <Topology>
    <Reliability>86.4</Reliability>
    <Numtm Count="1"/>
    <Regions>
      <Region Begin="1" End="33" Loc="Signal"/>
      <Region Begin="34" End="318" Loc="Outside"/>
      <Region Begin="319" End="340" Loc="Membrane"/>
      <Region Begin="341" End="489" Loc="Inside"/>
    </Regions>
  </Topology>
  <Experiments/>
</TOPDB>"""
CHAIN = f"""\
<TOPDB {NS} ID="2k3c_A" type="Alpha_Bitopic">
  <Name>2k3c_A</Name>
  <CrossRef><PDB ID="2k3c"/></CrossRef>
  <Topology>
    <Reliability>83.34</Reliability>
    <Numtm Count="1"/>
    <Regions>
      <Region Begin="1" End="1" Loc="Outside"/>
      <Region Begin="2" End="26" Loc="Membrane"/>
      <Region Begin="27" End="31" Loc="Inside"/>
    </Regions>
  </Topology>
  <Experiments/>
</TOPDB>"""
BETA = f"""\
<TOPDB {NS} ID="ALKL_PSEOL" type="Beta_Polytopic">
  <Name>Outer membrane protein AlkL</Name>
  <Membranes><Membrane>Cell outer membrane;Periplasmic space;Extracellular</Membrane></Membranes>
  <Topology>
    <Reliability>95.98</Reliability>
    <Numtm Count="8"/>
    <Regions>
      <Region Begin="1" End="27" Loc="Signal"/>
      <Region Begin="28" End="41" Loc="Periplasm"/>
      <Region Begin="42" End="53" Loc="Membrane"/>
      <Region Begin="54" End="73" Loc="Outside"/>
      <Region Begin="74" End="86" Loc="Membrane"/>
      <Region Begin="87" End="92" Loc="Periplasm"/>
      <Region Begin="93" End="102" Loc="Membrane"/>
      <Region Begin="103" End="121" Loc="Outside"/>
      <Region Begin="122" End="132" Loc="Membrane"/>
      <Region Begin="133" End="141" Loc="Periplasm"/>
      <Region Begin="142" End="154" Loc="Membrane"/>
      <Region Begin="155" End="166" Loc="Outside"/>
      <Region Begin="167" End="179" Loc="Membrane"/>
      <Region Begin="180" End="187" Loc="Periplasm"/>
      <Region Begin="188" End="195" Loc="Membrane"/>
      <Region Begin="196" End="216" Loc="Outside"/>
      <Region Begin="217" End="226" Loc="Membrane"/>
      <Region Begin="227" End="230" Loc="Periplasm"/>
    </Regions>
  </Topology>
  <Experiments/>
</TOPDB>"""


def document(*entries):
    body = '\n'.join(entries)
    return SimpleNamespace(result=io.BytesIO(f'<TOPDBs>\n{body}\n</TOPDBs>\n'.encode()))


def lines(*entries):
    """The same, as the iterator of lines that a download gives."""
    data = io.BytesIO(document(*entries).result.getvalue())
    return SimpleNamespace(result=iter(data.readlines()))


def entries(*xml):
    return list(topdb_parser.iter_entries(document(*xml)))


def extract(mapper, rows, dataset):
    out = SilverExtractor('topdb', dataset)
    for index, row in enumerate(rows):
        record = mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def test_reader_gives_the_type_name_and_membranes_of_each_entry():
    rows = entries(ADA2C, ANTR2, CHAIN)
    assert rows[0] == {
        'id': 'ADA2C_HUMAN', 'type': 'Alpha_Polytopic',
        'name': 'Alpha-2C adrenergic receptor', 'kind': 'uniprot_entry',
        'membranes': [
            {'name': 'Cell membrane', 'side_1': 'Cytoplasm', 'side_2': 'Extracellular'}
        ],
        'segments': 6, 'reliability': '70.6',
    }
    # The protein name is stripped; the entry has two membranes.
    assert rows[1]['name'] == 'Anthrax toxin receptor 2'
    assert [m['name'] for m in rows[1]['membranes']] == [
        'Cell membrane', 'Endoplasmic reticulum membrane',
    ]
    # An entry of a structure without a UniProt entry has no membrane.
    assert (rows[2]['kind'], rows[2]['membranes']) == ('pdb_chain', [])


def test_reader_reads_the_lines_of_a_download():
    rows = list(topdb_parser.iter_entries(lines(ADA2C, BETA)))
    assert [r['id'] for r in rows] == ['ADA2C_HUMAN', 'ALKL_PSEOL']


def test_a_repeated_identical_entry_is_read_once_and_a_different_one_is_an_error():
    assert [r['id'] for r in entries(ADA2C, ADA2C, BETA)] == ['ADA2C_HUMAN', 'ALKL_PSEOL']
    other = ADA2C.replace('Alpha_Polytopic', 'Alpha_Bitopic')
    with pytest.raises(ValueError, match='repeated with other content'):
        entries(ADA2C, other)


def test_reader_rejects_malformed_input():
    with pytest.raises(ValueError, match='no entries'):
        entries()
    with pytest.raises(ValueError, match='name;side;side'):
        entries(ADA2C.replace('Cell membrane;Cytoplasm;Extracellular', 'Cell membrane;Cytoplasm'))
    with pytest.raises(ValueError, match='entry name or a PDB chain'):
        entries(ADA2C.replace('ADA2C_HUMAN', 'ada2c'))


def test_reader_gives_the_segment_count_and_the_reliability_of_each_entry():
    rows = {r['id']: r for r in entries(ADA2C, ANTR2, CHAIN, BETA)}
    assert {k: (r['segments'], r['reliability']) for k, r in rows.items()} == {
        'ADA2C_HUMAN': (6, '70.6'),
        'ANTR2_HUMAN': (1, '86.4'),
        '2k3c_A': (1, '83.34'),
        'ALKL_PSEOL': (8, '95.98'),
    }


def test_reader_raises_when_numtm_is_not_the_number_of_membrane_regions():
    with pytest.raises(ValueError, match='Numtm 7 but 6 Membrane regions'):
        entries(ADA2C.replace('<Numtm Count="6"/>', '<Numtm Count="7"/>'))
    # One region less is also a mismatch, whichever the count says.
    with pytest.raises(ValueError, match='Numtm 6 but 5 Membrane regions'):
        entries(ADA2C.replace('<Region Begin="323" End="334" Loc="Membrane"/>', ''))
    # The entry of a PDB structure is checked as well.
    with pytest.raises(ValueError, match='Numtm 2 but 1 Membrane regions'):
        entries(CHAIN.replace('<Numtm Count="1"/>', '<Numtm Count="2"/>'))


def test_reader_raises_on_a_missing_or_invalid_topology():
    topology = ADA2C[ADA2C.index('<Topology>'):ADA2C.index('</Topology>') + len('</Topology>')]
    with pytest.raises(ValueError, match='has no topology'):
        entries(ADA2C.replace(topology, ''))
    with pytest.raises(ValueError, match='no Numtm or no Reliability'):
        entries(ADA2C.replace('<Reliability>70.6</Reliability>', ''))
    with pytest.raises(ValueError, match='no valid Numtm Count or Reliability'):
        entries(ADA2C.replace('<Reliability>70.6</Reliability>', '<Reliability>high</Reliability>'))
    with pytest.raises(ValueError, match='no valid Numtm Count or Reliability'):
        entries(ADA2C.replace('<Numtm Count="6"/>', '<Numtm/>'))
    with pytest.raises(ValueError, match='Reliability 170.6'):
        entries(ADA2C.replace('<Reliability>70.6</Reliability>', '<Reliability>170.6</Reliability>'))


def test_term_reader_gives_membranes_and_types_once():
    rows = list(topdb_parser.iter_terms(document(ADA2C, ANTR2, ADA2C, BETA)))
    assert [(r['kind'], r['id']) for r in rows] == [
        ('membrane', 'Cell membrane'),
        ('membrane', 'Cell outer membrane'),
        ('membrane', 'Endoplasmic reticulum membrane'),
        ('type', 'Alpha_Bitopic'),
        ('type', 'Alpha_Polytopic'),
        ('type', 'Beta_Polytopic'),
    ]
    with pytest.raises(ValueError, match='two side definitions'):
        list(topdb_parser.iter_terms(document(
            ADA2C, ANTR2.replace('Cell membrane;Cytoplasm;Extracellular', 'Cell membrane;Cytoplasm;Periplasmic space'),
        )))


def test_terms_have_their_namespace_name_and_sides():
    rows = list(topdb_parser.iter_terms(document(ADA2C, BETA)))
    out = extract(topdb.map_term, rows, 'terms')
    by_id = {(e.namespace, e.identifier): e for e in out.entities.values()}
    membrane = by_id[(topdb.MEMBRANE, 'Cell outer membrane')]
    assert membrane.entity_type == 'ontology_class'
    assert {i['id'] for i in membrane.identifiers if i['ns'] == Namespace.NAME} == {
        'Cell outer membrane'
    }
    assert {a['term']: a['value'] for a in membrane.annotations} == {
        'topdb:side_1': 'Periplasmic space',
        'topdb:side_2': 'Extracellular',
    }
    kind = by_id[(topdb.TYPE, 'Beta_Polytopic')]
    assert kind.entity_type == 'ontology_class' and not kind.annotations
    with pytest.raises(ValueError, match='Unknown TOPDB term kind'):
        topdb.map_term({'kind': 'segment', 'id': 'x', 'name': 'x'})


def test_a_protein_is_located_in_its_membranes_and_has_its_type():
    rows = entries(ANTR2)
    out = extract(topdb.map_protein, rows, 'proteins')
    subject = next(e for e in out.entities.values() if e.entity_type == 'protein')
    assert (subject.namespace, subject.identifier) == (
        Namespace.UNIPROT_ENTRY, 'ANTR2_HUMAN'
    )
    # No species is asserted: the entry name's suffix is not a taxon.
    assert subject.taxon is None
    statements = {
        (
            r.predicate,
            out.entities[r.object_entity_key].entity_type,
            out.entities[r.object_entity_key].namespace,
            out.entities[r.object_entity_key].identifier,
        )
        for r in out.relations
    }
    assert statements == {
        ('located_in', 'ontology_class', topdb.MEMBRANE, 'Cell membrane'),
        ('located_in', 'ontology_class', topdb.MEMBRANE, 'Endoplasmic reticulum membrane'),
        ('associated_with', 'ontology_class', topdb.TYPE, 'Alpha_Bitopic'),
    }


def test_the_type_statement_carries_the_segment_count_and_the_reliability():
    out = extract(topdb.map_protein, entries(ADA2C, ANTR2), 'proteins')
    attributes = {}
    for relation in out.relations:
        subject = out.entities[relation.subject_entity_key].identifier
        obj = out.entities[relation.object_entity_key]
        attributes[(subject, relation.predicate, obj.identifier)] = {
            a['term']: a for a in relation.annotations
        }
    type_statement = attributes[('ADA2C_HUMAN', 'associated_with', 'Alpha_Polytopic')]
    assert set(type_statement) == {'topdb:transmembrane_segments', 'topdb:reliability'}
    # A count, a measurement without a unit (not the text '6').
    segments = type_statement['topdb:transmembrane_segments']
    assert segments['value'] == ''
    assert segments['quantity']['has_numeric_value'] == 6
    assert segments['quantity']['has_unit'] is None
    assert segments['quantity']['source_field'] == 'Numtm Count'
    reliability = type_statement['topdb:reliability']
    # A measurement without a unit, the value as given.
    assert reliability['value'] == ''
    assert reliability['quantity']['has_numeric_value'] == 70.6
    assert reliability['quantity']['has_unit'] is None
    assert reliability['quantity']['source_field'] == 'Reliability'
    other = attributes[('ANTR2_HUMAN', 'associated_with', 'Alpha_Bitopic')]
    assert other['topdb:transmembrane_segments']['quantity']['has_numeric_value'] == 1
    assert other['topdb:reliability']['quantity']['has_numeric_value'] == 86.4


def test_the_membrane_statements_are_unchanged_and_the_protein_is_one_entity():
    out = extract(topdb.map_protein, entries(ANTR2), 'proteins')
    membrane_statements = [r for r in out.relations if r.predicate == 'located_in']
    assert {out.entities[r.object_entity_key].identifier for r in membrane_statements} == {
        'Cell membrane', 'Endoplasmic reticulum membrane',
    }
    # The membranes are associations of the protein: nothing is attached.
    assert all(not r.annotations for r in membrane_statements)
    # The type statement and the membrane statements share one subject: the
    # same identifiers, no molecular form.
    assert len({r.subject_entity_key for r in out.relations}) == 1
    proteins = [e for e in out.entities.values() if e.entity_type == 'protein']
    assert len(proteins) == 1 and not proteins[0].molecular_form
    assert {
        (i['ns'], i['id']) for i in proteins[0].identifiers
    } == {
        (Namespace.UNIPROT_ENTRY, 'ANTR2_HUMAN'),
        (Namespace.NAME, 'Anthrax toxin receptor 2'),
    }


def test_a_protein_without_a_membrane_has_only_its_type():
    rows = [dict(entries(BETA)[0], membranes=[])]
    out = extract(topdb.map_protein, rows, 'proteins')
    assert {r.predicate for r in out.relations} == {'associated_with'}


def test_statement_objects_are_the_entities_of_the_terms_dataset():
    xml = (ADA2C, ANTR2, BETA)
    term_out = extract(topdb.map_term, list(topdb_parser.iter_terms(document(*xml))), 'terms')
    protein_out = extract(topdb.map_protein, entries(*xml), 'proteins')
    term_keys = {
        (e.entity_type, e.namespace, e.identifier) for e in term_out.entities.values()
    }
    objects = {
        (o.entity_type, o.namespace, o.identifier)
        for o in (
            protein_out.entities[r.object_entity_key] for r in protein_out.relations
        )
    }
    assert objects and objects <= term_keys


def test_entries_of_a_pdb_structure_give_no_statement_and_no_compartment_is_a_membrane():
    rows = entries(CHAIN, ADA2C)
    assert topdb.exclusion_reason(rows[0]) == 'no UniProt entry name'
    assert topdb.map_protein(rows[0]) is None
    out = extract(topdb.map_protein, rows, 'proteins')
    assert len([e for e in out.entities.values() if e.entity_type == 'protein']) == 1
    # "Cytoplasm" and "Extracellular" are the sides of the membrane, not terms.
    names = {
        e.identifier for e in out.entities.values() if e.entity_type == 'ontology_class'
    }
    assert names == {'Cell membrane', 'Alpha_Polytopic'}


def test_resource_datasets():
    assert set(topdb.resource.datasets()) == {'terms', 'proteins'}
