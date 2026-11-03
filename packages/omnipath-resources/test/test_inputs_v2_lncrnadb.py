"""Source contract of the lncRNAdb association input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.lncrnadb import resource
from pypath.inputs_v2.parsers.lncrnadb import _pmid_list, iter_associations

# Entries of the lncRNAdb export (the descriptions are shortened).
XML = """<?xml version="1.0" ?>
<Search_Results>
  <Results id="0">
    <nomenclature>
      <Name>G22</Name>
      <Alias>None known</Alias>
      <Biotype>unassigned</Biotype>
    </nomenclature>
    <species>
      <Entry Species="Loris tardigradus (Slender loris)">N/A</Entry>
      <Entry Species="Otolemur garnettii (Garnett's greater bushbaby)">N/A</Entry>
    </species>
    <association>
      <association id="0">
        <ComponentType>Protein</ComponentType>
        <ComponentID>PABP</ComponentID>
        <description>G22 binds poly(A) binding protein</description>
        <pubmedID>17175535</pubmedID>
      </association>
    </association>
  </Results>
  <Results id="1">
    <nomenclature>
      <Name>TUSC7</Name>
      <Alias>LINC00902, LSAMP-AS1, ENSG00000243197</Alias>
      <Biotype>unassigned</Biotype>
    </nomenclature>
    <species>
      <Entry Species="Homo sapiens (Human)">hg19 chr3:116,709,788-116,717,040</Entry>
    </species>
    <association/>
  </Results>
  <Results id="2">
    <nomenclature>
      <Name>PTENP1</Name>
      <Alias>PTEN pseudogene 1; PTENpg1</Alias>
      <Biotype>unassigned</Biotype>
    </nomenclature>
    <species>
      <Entry Species="Homo sapiens (Human)">N/A</Entry>
    </species>
    <association>
      <association id="0">
        <ComponentType>Transcript</ComponentType>
        <ComponentID>miR-20a </ComponentID>
        <description>PTENP1 acts as a decoy for miR-20a</description>
        <pubmedID>20577206</pubmedID>
      </association>
      <association id="1">
        <ComponentType>PROTEIN</ComponentType>
        <ComponentID>EZH2</ComponentID>
        <description>NULL</description>
        <pubmedID>NULL</pubmedID>
      </association>
      <association id="2">
        <ComponentType>RNA</ComponentType>
        <ComponentID>miR-21</ComponentID>
        <description>Validated by luciferase assay</description>
        <pubmedID>Validated</pubmedID>
      </association>
      <association id="3">
        <ComponentType>Gene</ComponentType>
        <ComponentID>PTEN</ComponentID>
        <description>Two studies</description>
        <pubmedID>21172659, 21874018</pubmedID>
      </association>
      <association id="4">
        <ComponentType>Protein</ComponentType>
        <ComponentID>p68</ComponentID>
        <description>confirmed by Western</description>
        <pubmedID>:22665494</pubmedID>
      </association>
      <association id="5">
        <ComponentType></ComponentType>
        <ComponentID></ComponentID>
        <description></description>
        <pubmedID>None</pubmedID>
      </association>
    </association>
  </Results>
  <Results id="3">
    <nomenclature>
      <Name>Prion-associated RNAs</Name>
      <Alias></Alias>
      <Biotype>unassigned</Biotype>
    </nomenclature>
    <species>
      <Entry Species="Human herpesvirus">N/A</Entry>
    </species>
    <association>
      <association id="0">
        <ComponentType>ncRNA</ComponentType>
        <ComponentID>miR-715</ComponentID>
        <description>Functions in cis</description>
        <pubmedID>Confirmed</pubmedID>
      </association>
    </association>
  </Results>
</Search_Results>
"""


def source(text=XML):
    return SimpleNamespace(result=io.StringIO(text))


def rows():
    return list(iter_associations(source()))


def extract(row):
    record = resource.datasets()['associations'].mapper(row)
    assert record is not None
    out = SilverExtractor('lncrnadb', 'associations')
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
        (a['term'], a['value']) for a in relation.annotations if a['scope'] == scope
    )


def by_partner(partner):
    (row,) = [r for r in rows() if r['partner'] == partner]
    return row


def test_reader_yields_one_row_per_association_and_none_for_empty_entries():
    assert [r['lncrna'] for r in rows()] == [
        'G22', 'PTENP1', 'PTENP1', 'PTENP1', 'PTENP1', 'PTENP1', 'PTENP1',
        'Prion-associated RNAs',
    ]


def test_association_is_an_interaction_from_the_lncrna_to_the_partner():
    relation, lncrna, partner = extract(by_partner('miR-20a'))
    assert relation.predicate == 'interacts_with'
    assert lncrna.entity_type == 'noncoding_rna_product'
    assert identifiers(lncrna) == {
        (Namespace.NAME, 'PTENP1'),
        (Namespace.SYNONYM, 'PTEN pseudogene 1'),
        (Namespace.SYNONYM, 'PTENpg1'),
    }
    assert partner.entity_type == 'transcript'
    assert identifiers(partner) == {(Namespace.NAME, 'miR-20a')}
    assert annotations(relation) == [
        ('description', 'PTENP1 acts as a decoy for miR-20a'),
        ('publications', 'PMID:20577206'),
    ]


def test_partner_type_is_the_components_own_type():
    assert extract(by_partner('EZH2'))[2].entity_type == 'protein'
    assert extract(by_partner('PTEN'))[2].entity_type == 'gene'
    assert extract(by_partner('miR-21'))[2].entity_type == 'rna_product'
    assert extract(by_partner('miR-715'))[2].entity_type == 'noncoding_rna_product'


def test_lncrna_of_one_species_has_its_taxon_and_of_several_has_none():
    assert extract(by_partner('miR-20a'))[1].taxon == '9606'
    assert extract(by_partner('PABP'))[1].taxon is None
    # The partner's species is not stated.
    assert extract(by_partner('miR-20a'))[2].taxon is None


def test_group_name_of_viruses_is_no_taxon():
    assert extract(by_partner('miR-715'))[1].taxon is None


def test_publication_cell_is_read_for_pubmed_ids_only():
    assert annotations(extract(by_partner('PTEN'))[0]) == [
        ('description', 'Two studies'),
        ('publications', 'PMID:21172659'),
        ('publications', 'PMID:21874018'),
    ]
    assert [t for t, _ in annotations(extract(by_partner('p68'))[0])] == [
        'description', 'publications',
    ]
    assert ('publications', 'PMID:22665494') in annotations(extract(by_partner('p68'))[0])
    # Status words in the PubMed ID cell are not publications, and ``NULL`` is no description.
    assert annotations(extract(by_partner('EZH2'))[0]) == []
    assert by_partner('miR-21')['pmids'] == []
    assert [t for t, _ in annotations(extract(by_partner('miR-21'))[0])] == ['description']


def test_placeholder_association_gives_no_relation():
    (placeholder,) = [r for r in rows() if r['partner'] is None]
    assert placeholder['partner_type'] is None
    assert resource.datasets()['associations'].mapper(placeholder) is None


def test_aliases_that_say_none_are_no_identifiers():
    _, lncrna, _ = extract(by_partner('PABP'))
    assert identifiers(lncrna) == {(Namespace.NAME, 'G22')}


def test_unknown_partner_type_and_species_are_errors():
    mapper = resource.datasets()['associations'].mapper
    with pytest.raises(ValueError, match='Unknown partner type'):
        mapper(dict(by_partner('PABP'), partner_type='Lipid'))
    with pytest.raises(ValueError, match='Unknown species'):
        mapper(dict(by_partner('miR-20a'), species=['Equus ferus']))


def test_reader_rejects_other_documents_and_empty_exports():
    with pytest.raises(ValueError, match='Search_Results'):
        list(iter_associations(source('<html><body>for sale</body></html>')))
    with pytest.raises(ValueError, match='no entries'):
        list(iter_associations(source('<Search_Results></Search_Results>')))
    empty = (
        '<Search_Results><Results><nomenclature><Name>X</Name></nomenclature>'
        '<species/><association/></Results></Search_Results>'
    )
    with pytest.raises(ValueError, match='no associations'):
        list(iter_associations(source(empty)))


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    assert _pmid_list('21172659, 21874018') == ['21172659', '21874018']
    assert _pmid_list(':22665494') == ['22665494']
    assert _pmid_list('953') == ['953']
    assert _pmid_list('0123456, 123456789') == []
    assert _pmid_list('Validated') == []
