"""LMPID motif mediated interactions: parser, mapper and extractor on source records."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import lmpid
from pypath.inputs_v2.parsers.lmpid import iter_lmpid

# The first two records of LMPID_DATA_pubmed_ref.xml, with its line endings.
# Tests append further records to them.
SOURCE = b"""<?xml version="1.0"?>\r
<RECORDSET name="LMPID">\r
<RECORD>\r
<Bait_UniProt_ID>P54253</Bait_UniProt_ID>\r
<Motif_Instance>RRWSAP</Motif_Instance>\r
<Sequence_Position>772-777</Sequence_Position>\r
<Interacting_Domain>14-3-3 domain</Interacting_Domain>\r
<Prey_UniProt_ID>P92177</Prey_UniProt_ID>\r
<References>20037628, 12757707</References>\r
</RECORD>\r
<RECORD>\r
<Bait_UniProt_ID>P54253</Bait_UniProt_ID>\r
<Motif_Instance>RRWSAP</Motif_Instance>\r
<Sequence_Position>772-777</Sequence_Position>\r
<Interacting_Domain>14-3-3 domain</Interacting_Domain>\r
<Prey_UniProt_ID>P62258</Prey_UniProt_ID>\r
<References>20037628, 12757707</References>\r
</RECORD>\r
%s</RECORDSET>\r
"""


def parse(extra=b''):
    opener = SimpleNamespace(result=io.BytesIO(SOURCE % extra))
    return list(iter_lmpid(opener))


def record(bait, motif, position, domain, prey, references):
    return (
        '<RECORD>\r\n'
        f'<Bait_UniProt_ID>{bait}</Bait_UniProt_ID>\r\n'
        f'<Motif_Instance>{motif}</Motif_Instance>\r\n'
        f'<Sequence_Position>{position}</Sequence_Position>\r\n'
        f'<Interacting_Domain>{domain}</Interacting_Domain>\r\n'
        f'<Prey_UniProt_ID>{prey}</Prey_UniProt_ID>\r\n'
        f'<References>{references}</References>\r\n'
        '</RECORD>\r\n'
    ).encode()


def extract(rows):
    out = SilverExtractor('lmpid', 'interactions')
    for index, row in enumerate(rows):
        relation = lmpid.resource.interactions.mapper(row)
        if relation is not None:
            out.process_record(relation, row, f'fixture:{index}', index)
    return out


def values(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_parser_reads_motif_position_domain_and_references():
    rows = parse()
    assert rows == [
        {
            'bait_uniprot': 'P54253',
            'prey_uniprot': prey,
            'motif_instance': 'RRWSAP',
            'sequence_position': '772-777',
            'interacting_domain': '14-3-3 domain',
            'references': '20037628, 12757707',
            'motif_start': 772,
            'motif_end': 777,
            'pmids': ['20037628', '12757707'],
        }
        for prey in ('P92177', 'P62258')
    ]


def test_parser_normalises_whitespace_unknown_domain_and_odd_values():
    rows = parse(
        # Source records: an accession with a trailing space, an unknown domain.
        record('P17612 ', 'FSEF', '348-351', 'PIF-binding pocket', 'O15530', '10698939')
        + record('Q13522', 'RKIQF', '8-12', 'Unknown', 'P62136', '8611507')
        # Not in this export: a date where a range was (as in the website's
        # own export), and a reference that is not a PubMed ID.
        + record('Q12986', 'FKFNTDAAEFIPQ', '23-Nov', 'MLLE domain', 'Q13310', 'see text, 15663938 ')
    )
    assert rows[2]['bait_uniprot'] == 'P17612'
    assert rows[3]['interacting_domain'] is None
    assert rows[3]['pmids'] == ['8611507']
    assert (rows[4]['sequence_position'], rows[4]['motif_start'], rows[4]['motif_end']) == (
        '23-Nov', None, None,
    )
    assert rows[4]['references'] == 'see text, 15663938'
    assert rows[4]['pmids'] == ['15663938']


def test_parser_rejects_a_document_that_is_not_the_export():
    opener = SimpleNamespace(result=io.BytesIO(b'<html><body>Not found</body></html>'))
    with pytest.raises(ValueError, match='no RECORD elements'):
        list(iter_lmpid(opener))
    assert list(iter_lmpid(SimpleNamespace(result=io.BytesIO(b'')))) == []


def test_interaction_links_the_two_proteins_with_its_publications():
    out = extract(parse())
    first, second = out.relations
    # The export has no interaction type per record: an interaction, not a
    # stated direct one.
    assert first.predicate == 'interacts_with'
    assert values(first, 'publications') == ['PMID:20037628', 'PMID:12757707']
    # The same bait in both records, a different prey.
    assert first.subject_entity_key == second.subject_entity_key
    assert first.object_entity_key != second.object_entity_key
    # The bait carries the motif region and each prey its domain region, so
    # the entities are those proteins in those forms.
    assert sorted(
        (e.entity_type, e.namespace, e.identifier) for e in out.entities.values()
    ) == [
        ('protein', Namespace.UNIPROT, 'P54253'),
        ('protein', Namespace.UNIPROT, 'P62258'),
        ('protein', Namespace.UNIPROT, 'P92177'),
    ]


def region(participant):
    (region,) = participant.molecular_form['regions']
    return region


def test_motif_instance_is_a_region_of_the_bait_and_the_domain_one_of_the_prey():
    relation = lmpid.resource.interactions.mapper(parse()[0])
    assert region(relation.subject) == {
        'type': 'linear motif',
        'identifier': None,
        'position': 772,
        'end_position': 777,
        'coordinate_reference': {
            'identifier': None,
            'coordinate_system': 'protein',
            'position_base': 1,
        },
        'description': 'RRWSAP',
    }
    assert region(relation.object) == {
        'type': 'interacting domain',
        'identifier': None,
        'position': None,
        'end_position': None,
        'coordinate_reference': None,
        'description': '14-3-3 domain',
    }
    out = extract(parse())
    forms = {
        e.identifier: e.molecular_form['regions'][0]['description']
        for e in out.entities.values()
    }
    assert forms == {'P54253': 'RRWSAP', 'P62258': '14-3-3 domain', 'P92177': '14-3-3 domain'}


def test_records_that_differ_only_by_their_motif_are_different_relations():
    # Source records: one bait and prey with two motif instances.
    rows = parse(
        record('Q9BYG3', 'RRQKL', '10-14', 'Unknown', 'Q9Y266', '15040001')
        + record('Q9BYG3', 'KKRKR', '230-234', 'Unknown', 'Q9Y266', '15040001')
    )
    out = extract(rows[2:])
    first, second = out.relations
    assert first.relation_key != second.relation_key
    assert first.subject_entity_key != second.subject_entity_key
    # A prey whose domain is unknown has no region.
    assert out.entities[first.object_entity_key].molecular_form is None


def test_motif_without_a_readable_position_keeps_its_residues_only():
    rows = parse(record('Q12986', 'FKFNTDAAEFIPQ', '23-Nov', 'MLLE domain', 'Q13310', '15663938'))
    form = region(lmpid.resource.interactions.mapper(rows[2]).subject)
    assert form['description'] == 'FKFNTDAAEFIPQ'
    assert (form['position'], form['end_position'], form['coordinate_reference']) == (None, None, None)


def test_no_organism_is_asserted():
    relation = lmpid.resource.interactions.mapper(parse()[0])
    # The export names no organism; P92177 is a fly protein, P54253 a human one.
    assert not relation.subject.annotations
    assert not relation.object.annotations
    assert {slot_name(a.term) for a in relation.annotations} == {'publications'}
    (observed, _) = extract(parse()).relations
    assert {a['term'] for a in observed.annotations} == {'publications'}
    assert all(e.taxon is None for e in extract(parse()).entities.values())


def test_self_interaction_is_kept_and_records_without_a_partner_are_not():
    # A source record: PDZ domain and motif in the same protein.
    rows = parse(record('Q8R4T5', 'EEESQL', '389-394', 'PDZ domain', 'Q8R4T5', '17396155'))
    out = extract(rows[2:])
    (relation,) = out.relations
    # The motif and the domain are in different forms of the one protein.
    ends = (relation.subject_entity_key, relation.object_entity_key)
    assert {out.entities[key].identifier for key in ends} == {'Q8R4T5'}

    mapper = lmpid.resource.interactions.mapper
    row = rows[0]
    assert mapper(dict(row, prey_uniprot=None)) is None
    assert mapper(dict(row, bait_uniprot=None)) is None
    # A gene name or an identifier of another kind is not taken for an accession.
    assert mapper(dict(row, prey_uniprot='YWHAE')) is None
    assert mapper(dict(row, prey_uniprot='ENSP00000264335')) is None
    # Without a PubMed ID the interaction stands, with no publication made up.
    assert not mapper(dict(row, pmids=[])).annotations
