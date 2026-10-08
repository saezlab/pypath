"""OBO tag values: qualifiers and comments are not part of the identifier."""
from pypath.inputs_v2.parsers.obo import parse_obo_text

STANZA = '''
[Term]
id: MONDO:0005015
name: diabetes mellitus
alt_id: MONDO:0000000 {source="x"}
is_a: MONDO:0000001 {source="EFO:0005932", source="NCIT:C2991"} ! disease
xref: DOID:9351 {source="MONDO:equivalentTo"}
relationship: disease_has_location UBERON:0000468 {source="y"} ! multicellular organism
'''


def test_qualifiers_do_not_become_part_of_identifiers():
    term, = parse_obo_text(STANZA)
    assert term['is_a'] == ['MONDO:0000001']
    assert term['xrefs'] == ['DOID:9351']
    assert term['alt_ids'] == ['MONDO:0000000']
    assert term['relationships'][0]['target'] == 'UBERON:0000468'
