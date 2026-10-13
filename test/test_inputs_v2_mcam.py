"""Source contract of the MCAM input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import KnowledgeOrigin

from pypath.inputs_v2.mcam import LIST_ID, LIST_NAMESPACE, config, resource
from pypath.inputs_v2.parsers.mcam import iter_list, iter_proteins

HEADER = 'ACC\tENTNAME\tFIELD2\tHUNIPORTID\tLENGTH\tMATCHEDQUERYID\tNF50\tNF90\tORG\tTAXGRP'
# Rows of the human list of MCAM. The file ends its lines with a carriage return
# and its rows with a tab.
ADAM8 = 'P78325\tADAM 8 precursor; CD156; Cell surface antigen MS2; CD156a antigen; A disintegrin and metalloproteinase domain 8\tADAM8\tADAM8_HUMAN\t824\tADAM8_HUMAN, \tUniRef50_P78325\tUniRef90_P78325\tHomo sapiens\tEuk/mammal\t'
AMGO1 = 'Q86WK6; Q8IW71; Q9ULQ7\tAmphoterin-induced protein 1 precursor; AMIGO-1; Alivin-2\tAMGO1\tAMGO1_HUMAN\t493\tAMGO1_HUMAN, \tUniRef50_Q86WK6\tUniRef90_Q86WK6\tHomo sapiens\tEuk/mammal\t'
C1QR1 = 'Q9NPY3; O00274\tComplement component C1q receptor precursor; C1qRp; C1qR; C1qR(p); Complement component 1 q subcomponent receptor 1; C1q/MBL/SPA receptor; CDw93; CD93 antigen\tC1QR1\tC1QR1_HUMAN\t652\tC1QR1_HUMAN, \tUniRef50_Q9NPY3\tUniRef90_Q9NPY3\tHomo sapiens\tEuk/mammal\t'


def opener(*lines, header=HEADER):
    return SimpleNamespace(result=io.StringIO('\r'.join([header, *lines, ''])))


def proteins(*lines, header=HEADER):
    return list(iter_proteins(opener(*lines, header=header)))


def extract(record, dataset='proteins'):
    out = SilverExtractor('mcam', dataset)
    out.process_record(resource.datasets()[dataset].mapper(record), record, 'fixture:1', 1)
    return out


def test_parser_reads_a_file_with_carriage_returns_and_a_trailing_tab():
    assert proteins(ADAM8, AMGO1, C1QR1) == [
        {'accessions': ['P78325'], 'entry_name': 'ADAM8_HUMAN', 'organism': 'Homo sapiens'},
        {
            'accessions': ['Q86WK6', 'Q8IW71', 'Q9ULQ7'],
            'entry_name': 'AMGO1_HUMAN',
            'organism': 'Homo sapiens',
        },
        {'accessions': ['Q9NPY3', 'O00274'], 'entry_name': 'C1QR1_HUMAN', 'organism': 'Homo sapiens'},
    ]


def test_protein_is_associated_with_the_list_by_an_association():
    out = extract(proteins(AMGO1)[0])
    (relation,) = out.relations
    assert relation.predicate == 'associated_with'
    protein = out.entities[relation.subject_entity_key]
    term = out.entities[relation.object_entity_key]
    assert (protein.entity_type, protein.namespace, protein.identifier, protein.taxon) == (
        'protein', 'uniprot', 'Q86WK6', '9606'
    )
    # The first accession of the entry is the primary one; the others are secondary.
    assert {(i['ns'], i['id']) for i in protein.identifiers} == {
        (Namespace.UNIPROT, 'Q86WK6'),
        (Namespace.UNIPROT_SEC, 'Q8IW71'),
        (Namespace.UNIPROT_SEC, 'Q9ULQ7'),
        (Namespace.UNIPROT_ENTRY, 'AMGO1_HUMAN'),
    }
    assert (term.entity_type, term.namespace, term.identifier) == (
        'ontology_class', LIST_NAMESPACE, LIST_ID
    )


def test_the_vocabulary_dataset_is_called_terms_and_the_list_is_mcam_list():
    assert 'terms' in resource.datasets()
    assert LIST_NAMESPACE == 'mcam_list'
    out = extract(proteins(ADAM8)[0])
    (relation,) = out.relations
    assert out.entities[relation.object_entity_key].namespace == 'mcam_list'


def test_the_list_is_a_term_with_a_name():
    (record,) = iter_list(opener(ADAM8, AMGO1))
    out = extract(record, 'terms')
    (term,) = out.entities.values()
    assert (term.entity_type, term.namespace, term.identifier) == (
        'ontology_class', LIST_NAMESPACE, LIST_ID
    )
    assert (Namespace.NAME, 'MCAM human cell adhesion molecule') in {
        (i['ns'], i['id']) for i in term.identifiers
    }
    assert not out.relations


def test_the_object_of_a_statement_is_the_term():
    protein_out = extract(proteins(ADAM8)[0])
    (relation,) = protein_out.relations
    obj = protein_out.entities[relation.object_entity_key]
    term_out = extract(next(iter_list(opener(ADAM8))), 'terms')
    (term,) = term_out.entities.values()
    assert (obj.entity_type, obj.namespace, obj.identifier) == (
        term.entity_type, term.namespace, term.identifier
    )


def test_no_species_other_than_human_and_no_attribute_beyond_the_list():
    out = extract(proteins(ADAM8)[0])
    (relation,) = out.relations
    assert {a['term'] for a in relation.annotations} <= {'in_taxon'}
    with pytest.raises(ValueError, match='organism'):
        proteins(ADAM8.replace('Homo sapiens', 'Mus musculus'))


def test_parser_raises_on_a_bad_header_a_repeated_entry_or_an_empty_file():
    with pytest.raises(ValueError, match='lacks expected columns'):
        proteins(ADAM8, header=HEADER.replace('ORG', 'SPECIES'))
    with pytest.raises(ValueError, match='occurs twice'):
        proteins(ADAM8, ADAM8)
    with pytest.raises(ValueError):
        proteins()
    with pytest.raises(ValueError, match='unexpected fields'):
        proteins(ADAM8 + 'x')
    with pytest.raises(ValueError, match='Empty accession'):
        proteins(ADAM8.replace('P78325', 'P78325; '))


def test_knowledge_origin_is_literature_curated():
    assert tuple(config.knowledge_origin) == (KnowledgeOrigin.LITERATURE_CURATED,)
