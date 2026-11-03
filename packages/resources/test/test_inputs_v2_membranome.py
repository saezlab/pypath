"""Membranome: parser, mapper and extractor on rows of its API tables."""

import io
import json
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.membranome import MEMBRANE_NAMESPACE, TAXA, config, resource
from pypath.inputs_v2.parsers.membranome import iter_localizations, iter_membranes

# Membranes 1, 2, 4, 5, 6 and 16 of the response of /membranes (it has 24; the
# fixture's total is its own number of objects).
MEMBRANES = [
    {'id': 1, 'ordering': 1.0, 'name': 'Plasma membrane', 'topology_in': 'cytoplasmic side',
     'topology_out': 'extracellular side', 'abbreviation': 'Plasma memb.'},
    {'id': 2, 'ordering': 2.0, 'name': 'Endoplasmic reticulum membrane',
     'topology_in': 'cytoplasmic side', 'topology_out': 'lumenal side', 'abbreviation': 'ER memb.'},
    {'id': 4, 'ordering': 4.0, 'name': 'Endosome membrane', 'topology_in': 'cytoplasmic side',
     'topology_out': 'endosome lumenal side', 'abbreviation': 'Endosome memb.'},
    {'id': 5, 'ordering': 5.0, 'name': 'Golgi apparatus membrane',
     'topology_in': 'cytoplasmic side', 'topology_out': 'Golgi lumenal side',
     'abbreviation': 'Golgi memb.'},
    {'id': 6, 'ordering': 6.0, 'name': 'Lysosomal membrane', 'topology_in': 'cytoplasmic side',
     'topology_out': 'lysosome lumenal side', 'abbreviation': 'Lysosome memb.'},
    {'id': 16, 'ordering': 16.0, 'name': 'Membrane', 'topology_in': 'inner side',
     'topology_out': 'outer side', 'abbreviation': 'Memb.'},
]
MEMBRANES_JSON = json.dumps({'total_objects': len(MEMBRANES), 'objects': MEMBRANES})

# Rows of the proteins table (CSV); only the columns the reader uses are kept.
PROTEIN_HEADER = 'id,name,uniprot_id,uniprotcode,genename,topology_show_in,membrane_id,species_name_cache'
PROTEINS = [
    # Starts outside (flag false); in the plasma membrane.
    '209,Activin receptor type-1,ACVR1_HUMAN,Q04771,ACVR1,false,1,Homo sapiens',
    # Four more membranes besides the plasma membrane.
    '97,"HLA class II histocompatibility antigen, DQ alpha 2 chain",DQA2_HUMAN,P01906,HLA-DQA2,'
    'false,1,Homo sapiens',
    # The database lists the two chains as two entries of one UniProt entry.
    '7564,N-acetylglucosamine-1-phosphotransferase subunits alpha/beta,GNPTA_HUMAN_alpha,Q3T906,'
    'GNPTAB,false,5,Homo sapiens',
    # Flag true: the topology is on the inner side.
    '121,Deacylase AADAC,AAAD_HUMAN,P22760,AADAC,true,2,Homo sapiens',
    # No topology.
    '332,APRG1 tumor suppressor candidate,APRG1_HUMAN,Q8IVJ8,APRG1,,16,Homo sapiens',
]
# Rows of protein_membranes (CSV, called protein_localizations on the download page).
LOCALIZATION_HEADER = 'id,membrane_name_cache,protein_uniprot_id_cache,pubmed,main,membrane_id,protein_id'
LOCALIZATIONS = [
    '22604,Plasma membrane,ACVR1_HUMAN,,true,1,209',
    '44,Endoplasmic reticulum membrane,DQA2_HUMAN,22407913,false,2,97',
    '45,Golgi apparatus membrane,DQA2_HUMAN,22407913,false,5,97',
    '46,Endosome membrane,DQA2_HUMAN,22407913,false,4,97',
    '47,Lysosomal membrane,DQA2_HUMAN,22407913,false,6,97',
    '22484,Plasma membrane,DQA2_HUMAN,,true,1,97',
    '8289,Golgi apparatus membrane,GNPTA_HUMAN_alpha,16120602,true,5,7564',
    '26621,Endoplasmic reticulum membrane,AAAD_HUMAN,19339378,true,2,121',
    '26705,Membrane,APRG1_HUMAN,,true,16,332',
]


class Fake:
    """A download whose opener holds a string."""

    def __init__(self, text):
        self.text = text

    def open(self, force_refresh=False, **_):
        return SimpleNamespace(result=io.StringIO(self.text))


def rows(localizations=LOCALIZATIONS, proteins=PROTEINS, membranes=MEMBRANES_JSON):
    return list(
        iter_localizations(
            SimpleNamespace(result=io.StringIO('\n'.join([LOCALIZATION_HEADER, *localizations]) + '\n')),
            membrane_download=Fake(membranes),
            protein_download=Fake('\n'.join([PROTEIN_HEADER, *proteins]) + '\n'),
        )
    )


def process(dataset, row):
    out = SilverExtractor('membranome', dataset)
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


def test_license_is_cc_by_3_as_the_site_links_it():
    assert config.license.name == 'CC_BY_3_0'


def test_a_protein_has_one_row_for_each_membrane_the_database_lists():
    parsed = rows()
    assert [(r['accession'], r['membrane']) for r in parsed] == [
        ('Q04771', 'Plasma membrane'),
        ('P01906', 'Endoplasmic reticulum membrane'),
        ('P01906', 'Golgi apparatus membrane'),
        ('P01906', 'Endosome membrane'),
        ('P01906', 'Lysosomal membrane'),
        ('P01906', 'Plasma membrane'),
        ('Q3T906', 'Golgi apparatus membrane'),
        ('P22760', 'Endoplasmic reticulum membrane'),
        ('Q8IVJ8', 'Membrane'),
    ]


def test_the_side_is_given_only_on_the_membrane_of_the_protein_table():
    parsed = {(r['accession'], r['membrane']): r for r in rows()}
    # Flag false: the membrane's topology_out side; flag true: its topology_in side.
    assert parsed[('Q04771', 'Plasma membrane')]['side'] == 'extracellular side'
    assert parsed[('P22760', 'Endoplasmic reticulum membrane')]['side'] == 'cytoplasmic side'
    assert parsed[('Q3T906', 'Golgi apparatus membrane')]['side'] == 'Golgi lumenal side'
    # Other membranes of the protein have no side, and neither has a protein without topology.
    assert parsed[('P01906', 'Plasma membrane')]['side'] == 'extracellular side'
    for membrane in ('Endoplasmic reticulum membrane', 'Golgi apparatus membrane',
                     'Endosome membrane', 'Lysosomal membrane'):
        assert parsed[('P01906', membrane)]['side'] is None
    assert parsed[('Q8IVJ8', 'Membrane')]['side'] is None
    assert parsed[('Q8IVJ8', 'Membrane')]['topology_flag'] is None


def test_statement_is_located_in_the_membrane_with_side_publication_and_main_flag():
    parsed = {(r['accession'], r['membrane']): r for r in rows()}
    protein, relation, membrane = relation_of(
        process('localizations', parsed[('P22760', 'Endoplasmic reticulum membrane')])
    )
    assert relation.predicate == 'located_in'
    assert (protein.entity_type, protein.namespace, protein.identifier, protein.taxon) == (
        'protein', 'uniprot', 'P22760', '9606'
    )
    assert {(i['ns'], i['id']) for i in protein.identifiers} == {
        (Namespace.UNIPROT, 'P22760'), (Namespace.UNIPROT_ENTRY, 'AAAD_HUMAN'),
        (Namespace.GENESYMBOL, 'AADAC'), (Namespace.NAME, 'Deacylase AADAC'),
    }
    assert (membrane.entity_type, membrane.namespace, membrane.identifier) == (
        'ontology_class', MEMBRANE_NAMESPACE, 'Endoplasmic reticulum membrane'
    )
    assert attributes(relation, 'membranome:topology_side') == ['cytoplasmic side']
    assert attributes(relation, 'publications') == ['PMID:19339378']
    assert attributes(relation, 'membranome:main_membrane') == ['true']


def test_a_secondary_membrane_is_marked_so_and_a_missing_publication_is_none():
    parsed = {(r['accession'], r['membrane']): r for r in rows()}
    _, relation, _ = relation_of(process('localizations', parsed[('P01906', 'Lysosomal membrane')]))
    assert attributes(relation, 'membranome:main_membrane') == ['false']
    assert attributes(relation, 'publications') == ['PMID:22407913']
    assert attributes(relation, 'membranome:topology_side') == []
    _, relation, _ = relation_of(process('localizations', parsed[('Q04771', 'Plasma membrane')]))
    assert attributes(relation, 'publications') == []


def test_an_entry_name_with_a_chain_suffix_is_not_a_uniprot_entry_name():
    parsed = {(r['accession'], r['membrane']): r for r in rows()}
    row = parsed[('Q3T906', 'Golgi apparatus membrane')]
    assert row['entry_name'] is None
    protein, _, _ = relation_of(process('localizations', row))
    assert Namespace.UNIPROT_ENTRY not in {i['ns'] for i in protein.identifiers}


def test_the_vocabulary_dataset_is_called_terms():
    assert set(resource.datasets()) == {'terms', 'localizations'}


def test_membranes_are_terms_with_names_and_abbreviations():
    membranes = list(iter_membranes(SimpleNamespace(result=io.StringIO(MEMBRANES_JSON))))
    assert [m['name'] for m in membranes][:2] == ['Plasma membrane', 'Endoplasmic reticulum membrane']
    out = SilverExtractor('membranome', 'terms')
    out.process_record(resource.terms.mapper(membranes[1]), membranes[1], 'fixture:1', 1)
    (entity,) = out.entities.values()
    assert (entity.entity_type, entity.namespace, entity.identifier) == (
        'ontology_class', MEMBRANE_NAMESPACE, 'Endoplasmic reticulum membrane'
    )
    assert {(i['ns'], i['id']) for i in entity.identifiers} == {
        (MEMBRANE_NAMESPACE, 'Endoplasmic reticulum membrane'),
        (Namespace.NAME, 'Endoplasmic reticulum membrane'),
        (Namespace.SYNONYM, 'ER memb.'),
    }


def test_objects_have_the_type_namespace_and_identifier_of_their_terms():
    defined = {}
    for membrane in iter_membranes(SimpleNamespace(result=io.StringIO(MEMBRANES_JSON))):
        out = SilverExtractor('membranome', 'terms')
        out.process_record(resource.terms.mapper(membrane), membrane, 'fixture:1', 1)
        (entity,) = out.entities.values()
        defined[entity.identifier] = (entity.entity_type, entity.namespace, entity.identifier)
    for row in rows():
        _, _, membrane = relation_of(process('localizations', row))
        assert defined[membrane.identifier] == (
            membrane.entity_type, membrane.namespace, membrane.identifier
        )


def test_no_membrane_is_asserted_for_a_protein_the_table_does_not_list_for_it():
    # Only the rows of protein_membranes are statements: the plasma membrane of
    # DQA2 comes from its own row, and no other protein gets it by default.
    parsed = rows()
    assert sum(1 for r in parsed if r['membrane'] == 'Plasma membrane') == 2


def test_species_are_the_eighteen_of_the_database():
    assert len(TAXA) == 18
    assert TAXA['Homo sapiens'] == 'NCBITaxon:9606'
    row = {**rows()[0], 'species': 'Pan troglodytes'}
    with pytest.raises(ValueError, match='Unknown Membranome species'):
        resource.localizations.mapper(row)


@pytest.mark.parametrize(
    'kwargs',
    [
        {'localizations': ['1,Plasma membrane,NOPE_HUMAN,,true,1,209']},  # the protein differs
        {'localizations': ['1,Lysosomal membrane,ACVR1_HUMAN,,true,1,209']},  # the membrane differs
        {'localizations': ['1,Plasma membrane,ACVR1_HUMAN,,true,1,9999']},  # unknown protein
        {'localizations': ['1,Plasma membrane,ACVR1_HUMAN,,maybe,1,209']},  # unknown flag
        {'proteins': [PROTEINS[0].replace(',false,', ',perhaps,')]},  # unknown topology flag
        {'localizations': []},  # no rows
    ],
)
def test_parser_raises_on_what_does_not_agree(kwargs):
    with pytest.raises(ValueError):
        rows(**kwargs)


def test_parser_raises_when_the_membranes_are_not_all_there():
    short = json.dumps({'total_objects': 24, 'objects': MEMBRANES})
    with pytest.raises(ValueError, match='all membranes'):
        rows(membranes=short)
