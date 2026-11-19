"""A KGML diagram's heading is not a pathway name and not a link to another pathway."""

import io

from omnipath_resources.inputs_v2.parsers import kegg

KGML = """\
<?xml version="1.0"?>
<pathway name="path:hsa00410" org="hsa" number="00410" title="beta-Alanine metabolism">
    <entry id="1" name="path:hsa00410" type="map">
        <graphics name="TITLE:beta-Alanine metabolism" type="roundrectangle"/>
    </entry>
    <entry id="2" name="path:hsa00640" type="map">
        <graphics name="Propanoate metabolism" type="roundrectangle"/>
    </entry>
</pathway>
"""


def parse():
    return kegg._parse_kgml_documents(
        {'hsa00410': io.StringIO(KGML)},
        compound_id_map={},
        uniprot_map={},
        taxon_id='9606',
    )['hsa00410']


def test_the_diagram_heading_is_not_a_linked_pathway():
    linked = parse()['linked_pathway_members']
    assert [member['pathway_id'] for member in linked] == ['hsa00640']


def test_a_linked_pathway_keeps_its_name():
    (member,) = parse()['linked_pathway_members']
    assert member['pathway_name'] == 'Propanoate metabolism'


def test_the_pathway_is_named_by_its_title_attribute():
    assert parse()['pathway_name'] == 'beta-Alanine metabolism'
