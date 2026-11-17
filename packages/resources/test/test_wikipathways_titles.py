"""A WikiPathways file names its pathway: title, organism and id come from the pathway itself."""

from rdflib import Graph

from omnipath_resources.inputs_v2.parsers import wikipathways

# A pathway file lists the pathway and, as data nodes, other pathways it refers to; both
# are typed wp:Pathway. The referenced one has no title and comes first.
PREFIXES = """
@prefix wp: <http://vocabularies.wikipathways.org/wp#> .
@prefix dc: <http://purl.org/dc/elements/1.1/> .
@prefix dcterms: <http://purl.org/dc/terms/> .
@prefix obo: <http://purl.obolibrary.org/obo/> .
@prefix wpid: <https://identifiers.org/wikipathways/> .
"""

OWN = """
wpid:WP106_r12345 a wp:Pathway ;
    dc:title "Alanine and aspartate metabolism"@en ;
    wp:organismName "Homo sapiens" ;
    wp:organism obo:NCBITaxon_9606 .
"""
REFERENCED = """
<https://identifiers.org/wikipathways/WP1> a wp:Pathway .
<https://identifiers.org/wikipathways/WP2> a wp:Pathway .
"""


def graph(*parts):
    g = Graph()
    g.parse(data=PREFIXES + "".join(parts), format="ttl")
    return g


def test_the_pathway_is_the_one_that_has_a_title_whatever_the_order_of_the_file():
    for parts in ((REFERENCED, OWN), (OWN, REFERENCED)):
        record = wikipathways._extract_pathway_record(graph(*parts))
        assert record["pathway_id"] == "WP106"
        assert record["pathway_version_id"] == "WP106_r12345"
        assert record["title"] == "Alanine and aspartate metabolism"
        assert record["organism_name"] == "Homo sapiens"
        assert record["taxon_id"] == "9606"


def test_a_file_without_any_titled_pathway_still_gives_a_stable_record():
    first = wikipathways._extract_pathway_record(graph(REFERENCED))
    second = wikipathways._extract_pathway_record(graph(*reversed([REFERENCED])))
    assert first == second
    assert first["pathway_id"] in {"WP1", "WP2"}


def test_a_file_without_a_pathway_gives_none():
    assert wikipathways._extract_pathway_record(graph("wpid:x a wp:Interaction .")) is None
