"""Source contract of the AlzPathway complexes input."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.alzpathway import MEMBER_TYPES, resource
from pypath.inputs_v2.parsers.alzpathway import (
    iter_celldesigner_complexes,
    iter_complexes,
)

_NAMESPACES = (
    'xmlns="http://www.sbml.org/sbml/level2/version4" '
    'xmlns:celldesigner="http://www.sbml.org/2001/ns/celldesigner" '
    'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" '
    'xmlns:bqbiol="http://biomodels.net/biology-qualifiers/"'
)


def _rdf(about, *resources):
    items = ''.join(f'<rdf:li rdf:resource="{r}"/>' for r in resources)
    return (
        f'<rdf:RDF><rdf:Description rdf:about="#{about}"><bqbiol:isVersionOf>'
        f'<rdf:Bag>{items}</rdf:Bag></bqbiol:isVersionOf></rdf:Description></rdf:RDF>'
    )


def _identity(cell_designer_class, reference=None, name=None, state=''):
    tag = {'PROTEIN': 'proteinReference', 'RNA': 'rnaReference'}.get(cell_designer_class)
    body = (
        f'<celldesigner:{tag}>{reference}</celldesigner:{tag}>'
        if reference
        else f'<celldesigner:name>{name}</celldesigner:name>'
    )
    return (
        '<celldesigner:speciesIdentity>'
        f'<celldesigner:class>{cell_designer_class}</celldesigner:class>{body}{state}'
        '</celldesigner:speciesIdentity>'
    )


def included(species_id, name, complex_id, cell_designer_class, reference=None, state=''):
    return (
        f'<celldesigner:species id="{species_id}" name="{name}"><celldesigner:annotation>'
        f'<celldesigner:complexSpecies>{complex_id}</celldesigner:complexSpecies>'
        f'{_identity(cell_designer_class, reference, name, state)}'
        '</celldesigner:annotation></celldesigner:species>'
    )


def species(species_id, name, compartment, cell_designer_class, reference=None, xrefs=()):
    return (
        f'<species id="{species_id}" name="{name}" compartment="{compartment}"><annotation>'
        f'<celldesigner:extension>{_identity(cell_designer_class, reference, name)}'
        f'</celldesigner:extension>{_rdf(species_id, *xrefs) if xrefs else ""}'
        '</annotation></species>'
    )


def reaction(reaction_id, reaction_type, products, xrefs):
    references = ''.join(f'<speciesReference species="{p}"/>' for p in products)
    return (
        f'<reaction id="{reaction_id}"><annotation><celldesigner:extension>'
        f'<celldesigner:reactionType>{reaction_type}</celldesigner:reactionType>'
        f'</celldesigner:extension>{_rdf(reaction_id, *xrefs)}</annotation>'
        f'<listOfProducts>{references}</listOfProducts></reaction>'
    )


_HOMODIMER = '<celldesigner:state><celldesigner:homodimer>5</celldesigner:homodimer></celldesigner:state>'
_PHOSPHORYLATED = (
    '<celldesigner:state><celldesigner:listOfModifications>'
    '<celldesigner:modification residue="rs1" state="phosphorylated"/>'
    '</celldesigner:listOfModifications></celldesigner:state>'
)

# Elements of alzpathway_sbml_map.xml (Additional file 1 of PMID 22647208),
# without their layout and notes. The CD3Z complex is shortened to two of its
# four proteins.
DOCUMENT = f"""<?xml version="1.0" encoding="UTF-8"?>
<sbml {_NAMESPACES} level="2" version="4">
<model id="AlzPathway_20120201">
<annotation><celldesigner:extension>
<celldesigner:listOfIncludedSpecies>
{included('s1028', 'ACh', 's3177', 'SIMPLE_MOLECULE')}
{included('s2839', 'CHRM', 's3177', 'PROTEIN', 'pr688')}
{included('s1282', 'C1Q', 's3219', 'PROTEIN', 'pr628')}
{included('s1283', 'C3B', 's3219', 'PROTEIN', 'pr631')}
{included('s1284', 'AmyloidB1-42', 's3219', 'PROTEIN', 'pr803', _HOMODIMER)}
{included('s2268', 'RXR', 's3295', 'PROTEIN', 'pr143')}
{included('s2407', 'PPARG', 's3295', 'PROTEIN', 'pr758')}
{included('s2594', 'unknown PPARG ligand', 's3295', 'UNKNOWN')}
{included('s2858', 'CD3Z complex', 's2451', 'COMPLEX')}
{included('s2458', 'CD247', 's2858', 'PROTEIN', 'pr850')}
{included('s2457', 'CD3E', 's2858', 'PROTEIN', 'pr852')}
{included('s408', 'APP', 's3324', 'PROTEIN', 'pr1')}
{included('s410', 'APPBP1', 's3324', 'PROTEIN', 'pr618')}
{included('s1696', 'APC', 's3259', 'PROTEIN', 'pr73')}
{included('s2878', 'GSK3B', 's3259', 'PROTEIN', 'pr55', _PHOSPHORYLATED)}
{included('s1457', 'APP', 's3233', 'RNA', 'rn14')}
{included('s1458', 'ACO1', 's3233', 'PROTEIN', 'pr673')}
</celldesigner:listOfIncludedSpecies>
<celldesigner:listOfProteins>
<celldesigner:protein id="pr688" name="CHRM" type="RECEPTOR"/>
<celldesigner:protein id="pr1" name="APP" type="GENERIC"/>
</celldesigner:listOfProteins>
</celldesigner:extension></annotation>
<listOfCompartments>
<compartment id="default"/>
<compartment id="c1" name="neuron"/>
<compartment id="c3" name="nucleus"/>
</listOfCompartments>
<listOfSpecies>
{species('s1025', 'ACh', 'default', 'SIMPLE_MOLECULE', xrefs=['urn:miriam:pubchem.substance:7885682'])}
{species('s1023', 'ACh', 'c1', 'SIMPLE_MOLECULE', xrefs=['urn:miriam:pubchem.substance:7885682'])}
{species('s1018', 'CHRM', 'c1', 'PROTEIN', 'pr688', ['urn:miriam:uniprot:P11229', 'urn:miriam:uniprot:?query=chrm'])}
{species('s3485', 'C1Q', 'default', 'PROTEIN', 'pr628', ['urn:miriam:uniprot:?query=c1q'])}
{species('s3489', 'C3B', 'default', 'PROTEIN', 'pr631', ['urn:miriam:pubchem.substance:192108'])}
{species('s1280', 'AmyloidB1-42', 'default', 'PROTEIN', 'pr803', ['urn:miriam:pubchem.substance:49661775'])}
{species('s1274', 'PPARG', 'c1', 'PROTEIN', 'pr758', ['urn:miriam:uniprot:P37231'])}
{species('s3798', 'PPARG', 'c3', 'PROTEIN', 'pr758', ['urn:miriam:uniprot:P37231'])}
{species('s1', 'APP', 'c1', 'PROTEIN', 'pr1', ['urn:miriam:uniprot:P05067', 'urn:miriam:uniprot:Q13564'])}
{species('s3177', 'ACh-CHRM', 'c1', 'COMPLEX')}
{species('s3219', 'C1Q-AmyloidB1-42-C3B', 'default', 'COMPLEX')}
{species('s3295', 'PPARG-RXR-unkown PPARG ligand', 'c3', 'COMPLEX')}
{species('s2451', 'ITAM', 'c1', 'COMPLEX')}
{species('s3324', 'APPBP1-APP', 'c1', 'COMPLEX')}
{species('s3259', 'APC-GSK3B-AXIN', 'c1', 'COMPLEX')}
{species('s3233', 'APP-ACO1', 'c1', 'COMPLEX')}
</listOfSpecies>
<listOfReactions>
{reaction('re322', 'HETERODIMER_ASSOCIATION', ['s3177'], ['urn:miriam:pubmed:17404842'])}
{reaction('re427', 'HETERODIMER_ASSOCIATION', ['s3219'], ['urn:miriam:pubmed:17900710'])}
{reaction('re785', 'HETERODIMER_ASSOCIATION', ['s3295', 's2258'], ['urn:miriam:kegg.pathway:hsa05010', 'urn:miriam:pubmed:18067990'])}
{reaction('re1', 'DISSOCIATION', ['s6', 's2410'], ['urn:miriam:pubmed:19524503'])}
</listOfReactions>
</model>
</sbml>
"""


def document(text=DOCUMENT):
    return SimpleNamespace(result=text.encode('utf-8'))


ROWS = {row['name']: row for row in iter_complexes(document())}


def extract(name):
    row = ROWS[name]
    record = resource.datasets()['complexes'].mapper(row)
    assert record is not None
    out = SilverExtractor('alzpathway', 'complexes')
    out.process_record(record, row, 'fixture:1', 1)
    (complex_key,) = [
        key
        for key, entity in out.entities.items()
        if (Namespace.NAME, name) in identifiers(entity)
    ]
    return out, complex_key


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def members(out, entity_key):
    return [
        (out.entities[relation.object_entity_key], relation)
        for relation in out.relations
        if relation.predicate == 'has_member'
        and relation.subject_entity_key == entity_key
    ]


def test_reader_gives_complexes_with_components_and_forming_reactions():
    row = ROWS['ACh-CHRM']
    assert row['species_id'] == 's3177'
    assert row['compartment'] == 'neuron'
    assert row['producing_reactions'] == [
        {
            'id': 're322',
            'type': 'HETERODIMER_ASSOCIATION',
            'xrefs': {'pubmed': ['17404842']},
        }
    ]
    assert [(m['name'], m['class'], m['reference']) for m in row['members']] == [
        ('ACh', 'SIMPLE_MOLECULE', None),
        ('CHRM', 'PROTEIN', 'pr688'),
    ]
    # A component has the cross-references of the species of the same molecule.
    assert row['members'][0]['xrefs'] == {'pubchem.substance': ['7885682']}
    assert row['members'][1]['xrefs'] == {'uniprot': ['P11229', '?query=chrm']}


def test_reader_reads_copies_modifications_and_nested_complexes():
    amyloid = ROWS['C1Q-AmyloidB1-42-C3B']['members'][2]
    assert (amyloid['name'], amyloid['homodimer']) == ('AmyloidB1-42', 5)
    gsk3b = ROWS['APC-GSK3B-AXIN']['members'][1]
    assert gsk3b['modifications'] == ['phosphorylated']
    (nested,) = ROWS['ITAM']['members']
    assert nested['class'] == 'COMPLEX'
    assert [m['name'] for m in nested['members']] == ['CD247', 'CD3E']
    assert ROWS['ITAM']['producing_reactions'] == []


def test_complex_has_its_components_and_the_publication_of_its_formation():
    out, key = extract('ACh-CHRM')
    complex_entity = out.entities[key]
    assert complex_entity.entity_type == 'macromolecular_complex'
    assert identifiers(complex_entity) == {(Namespace.NAME, 'ACh-CHRM')}
    assert [(a['term'], a['value']) for a in complex_entity.annotations] == [
        ('publications', 'PMID:17404842')
    ]
    found = {
        entity.entity_type: (identifiers(entity), entity.taxon)
        for entity, _ in members(out, key)
    }
    assert found == {
        'chemical_entity': (
            {(Namespace.PUBCHEM_SUBSTANCE, '7885682'), (Namespace.NAME, 'ACh')},
            None,
        ),
        # The search query in the annotation is not an accession.
        'protein': (
            {(Namespace.UNIPROT, 'P11229'), (Namespace.GENESYMBOL, 'CHRM')},
            '9606',
        ),
    }
    # The complex and its two components: no pairs between the components.
    assert [r.predicate for r in out.relations] == ['has_member', 'has_member']


def test_molecule_with_several_accessions_gets_none():
    out, key = extract('APPBP1-APP')
    found = {frozenset(identifiers(entity)) for entity, _ in members(out, key)}
    assert found == {
        frozenset({(Namespace.GENESYMBOL, 'APP')}),
        frozenset({(Namespace.GENESYMBOL, 'APPBP1')}),
    }
    assert not out.entities[key].annotations


def test_name_that_is_not_written_like_a_symbol_is_a_name():
    out, key = extract('C1Q-AmyloidB1-42-C3B')
    by_name = {
        next(iter(identifiers(entity)))[1]: (identifiers(entity), relation)
        for entity, relation in members(out, key)
    }
    assert by_name['AmyloidB1-42'][0] == {(Namespace.NAME, 'AmyloidB1-42')}
    assert by_name['C1Q'][0] == {(Namespace.GENESYMBOL, 'C1Q')}
    # A PubChem annotation of a protein is not an identifier of the protein.
    assert by_name['C3B'][0] == {(Namespace.GENESYMBOL, 'C3B')}
    # Five copies of the peptide are drawn; the other components have no number.
    assert [
        (a['term'], a['value'])
        for a in by_name['AmyloidB1-42'][1].annotations
        if a['scope'] == 'relation'
    ] == [('stoichiometry', '5')]
    assert not [
        a for a in by_name['C1Q'][1].annotations if a['term'] == 'stoichiometry'
    ]


def test_placeholder_component_is_not_a_member():
    out, key = extract('PPARG-RXR-unkown PPARG ligand')
    names = {value for entity, _ in members(out, key) for _, value in identifiers(entity)}
    assert names == {'P37231', 'PPARG', 'RXR'}
    # Only the publication of the reaction, not its KEGG pathway.
    assert [(a['term'], a['value']) for a in out.entities[key].annotations] == [
        ('publications', 'PMID:18067990')
    ]


def test_nested_complex_stays_a_complex_with_its_own_members():
    out, key = extract('ITAM')
    ((nested, relation),) = members(out, key)
    assert nested.entity_type == 'macromolecular_complex'
    assert identifiers(nested) == {(Namespace.NAME, 'CD3Z complex')}
    inner = {
        next(iter(identifiers(entity)))
        for entity, _ in members(out, relation.object_entity_key)
    }
    assert inner == {(Namespace.GENESYMBOL, 'CD247'), (Namespace.GENESYMBOL, 'CD3E')}
    # The proteins are members of the inner complex only.
    assert len(members(out, key)) == 1


def test_rna_component_is_an_rna():
    out, key = extract('APP-ACO1')
    types = {
        next(iter(identifiers(entity)))[1]: entity.entity_type
        for entity, _ in members(out, key)
    }
    assert types == {'APP': 'rna_product', 'ACO1': 'protein'}


def test_modification_state_and_compartment_are_not_asserted():
    out, key = extract('APC-GSK3B-AXIN')
    terms = {a['term'] for r in out.relations for a in r.annotations}
    terms |= {a['term'] for e in out.entities.values() for a in e.annotations}
    assert terms <= {'in_taxon'}
    assert all(e.molecular_form is None for e in out.entities.values())


def test_complex_of_placeholders_only_gives_no_record():
    row = dict(
        ROWS['ITAM'],
        members=[{'name': 'Unknown receptor', 'class': 'UNKNOWN', 'members': []}],
    )
    assert resource.datasets()['complexes'].mapper(row) is None


def test_unknown_component_class_is_an_error():
    assert 'PHENOTYPE' not in MEMBER_TYPES
    row = dict(
        ROWS['ITAM'],
        members=[{'name': 'apoptosis', 'class': 'PHENOTYPE', 'members': []}],
    )
    with pytest.raises(ValueError, match='PHENOTYPE'):
        resource.datasets()['complexes'].mapper(row)


def test_reader_rejects_a_document_that_is_not_a_celldesigner_map():
    with pytest.raises(ValueError, match='CellDesigner'):
        list(iter_celldesigner_complexes(document('<html><body/></html>')))


def test_map_without_complexes_is_rejected():
    empty = (
        f'<sbml {_NAMESPACES}><model><annotation><celldesigner:extension/>'
        '</annotation><listOfSpecies/></model></sbml>'
    )
    with pytest.raises(ValueError, match='no complex'):
        list(iter_complexes(document(empty)))
