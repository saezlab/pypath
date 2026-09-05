"""Native pathway/resource semantics exercised across the extraction boundary."""

import importlib

import pytest

from omnipath_build.silver import SilverExtractor


def extract(module, schema, row):
    record = getattr(
        importlib.import_module('pypath.inputs_v2.' + module), schema
    )(row)
    assert record is not None
    output = SilverExtractor(module, schema)
    output.process_record(record, row, 'fixture:1', 1)
    return output


@pytest.mark.parametrize(
    'module,schema,row',
    [
        ('brenda', 'schema', {'UniProt': ['P12345'], 'EC': '1.1.1.1'}),
        (
            'corum',
            'complexes_schema',
            {
                'ComplexID': '12',
                'ComplexName': 'Test complex',
                'subunits(UniProt IDs)': 'P12345;P23456',
                'Organism': 'Human',
                'GO ID': 'GO:0005737',
                'PubMed ID': '123',
            },
        ),
        (
            'intact',
            'interactions_schema',
            {
                '#ID(s) interactor A': 'uniprotkb:P12345',
                'ID(s) interactor B': 'chebi:CHEBI:15377',
                'Type(s) interactor A': 'psi-mi:"MI:0326"(protein)',
                'Type(s) interactor B': 'psi-mi:"MI:0328"(small molecule)',
                'Interaction type(s)': 'psi-mi:"MI:0915"(physical association)',
                'Taxid interactor A': 'taxid:9606(human)',
                'Publication Identifier(s)': 'pubmed:123',
            },
        ),
        (
            'kegg',
            'reactions_schema',
            {
                'reaction_id': 'R00001',
                'reaction_name': 'Test reaction',
                'reactant_kegg_id': 'C00001',
                'reactant_stoichiometry': '2',
                'product_kegg_id': 'C00002',
                'product_stoichiometry': '1',
                'uniprot_ids': 'P12345',
                'taxon_id': '9606',
            },
        ),
        (
            'macdb',
            'association_schema',
            {
                'pubchem_CID': '123',
                'trait_id': '1',
                'trait_name': 'Cancer',
                'pubmed_id': '123',
                'case_control_p-value': '0.01',
                'case_concentration': '1.5',
            },
        ),
        (
            'pfocr',
            'association_to_entity',
            {
                'data_type': 'gene',
                'identifier': '123',
                'pfocr_id': 'PMC123__fig1',
                'taxonomy_id': '9606',
                'caption': 'A figure',
            },
        ),
        (
            'reactome',
            'reactions_schema',
            {
                'entity_type': 'reaction',
                'reactome_stable_id': 'R-HSA-123',
                'participant_entity_type': 'chemical||protein',
                'participant_chebi': '15377||__MISSING__',
                'participant_uniprot': '__MISSING__||P12345',
                'participant_role': 'reactant||product',
                'participant_ncbi_tax_id': '__MISSING__||9606',
            },
        ),
        (
            'recon3d',
            'reactions_schema',
            {
                'bigg_reaction_id': 'R1',
                'reactants': 'h2o:c:2',
                'products': 'h:e:1',
                'reactant_chebi': '15377',
                'product_chebi': '15378',
                'enzyme_entrez': '123',
            },
        ),
        (
            'rhea',
            'reactions_schema',
            {
                'rhea_id': '123',
                'participant_chebi': 'CHEBI:15377||CHEBI:15378',
                'participant_role': 'reactant||product',
                'uniprot': 'P12345',
                'ec': '1.1.1.1',
                'pubmed': '123',
            },
        ),
        (
            'slctables',
            'schema',
            {
                'SLC name': 'SLC1A1',
                'Substrates': 'L-Glu',
                'Tissue and cellular expression': 'brain',
            },
        ),
        (
            'tcdb',
            '_transport_schema',
            {
                'transporter_uniprot': 'P12345',
                'substrate_id': 'CHEBI:15377',
                'tcid': '1.A.1.1.1',
            },
        ),
        (
            'wikipathways',
            'interactions_schema',
            {
                'source_entity_type': 'protein',
                'target_entity_type': 'protein',
                'source_uniprot': 'P12345',
                'target_uniprot': 'P23456',
                'interaction_types': 'Stimulation',
                'taxon_id': '9606',
            },
        ),
    ],
)
def test_native_resource_fixture(module, schema, row):
    output = extract(module, schema, row)
    assert output.entities
    for entity in output.entities.values():
        assert not entity.entity_type.startswith(('OM:', 'MI:'))
    for relation in output.relations:
        assert not relation.predicate.startswith(('OM:', 'MI:'))


def test_reaction_roles_are_edges_not_generic_membership():
    out = extract(
        'rhea',
        'reactions_schema',
        {
            'rhea_id': '123',
            'participant_chebi': '1||2',
            'participant_role': 'reactant||product',
            'uniprot': 'P12345',
        },
    )
    assert {r.predicate for r in out.relations} == {
        'has_input',
        'has_output',
        'enabled_by',
    }


def test_macdb_numeric_fields_are_not_narrative_or_causal():
    out = extract(
        'macdb',
        'association_schema',
        {
            'pubchem_CID': '1',
            'trait_id': '1',
            'case_concentration': '123',
            'case_control_p-value': '0.01',
        },
    )
    relation = out.relations[0]
    assert relation.predicate == 'associated_with'
    assert not any(a['term'] == 'description' for a in relation.annotations)
    assert any(
        a['term'] == 'p_value' and a['quantity']['has_numeric_value'] == 0.01
        for a in relation.annotations
    )


def test_unknown_direction_not_guessed():
    from pypath.inputs_v2 import reactome, wikipathways

    assert reactome._control_effect({'control_type': 'ACTIVATION-LIKE'}) is None
    assert (
        wikipathways._interaction_direction(
            {'interaction_types': 'Stimulation;Inhibition'}
        )
        is None
    )


def test_pfocr_mentions_and_invalid_identifier():
    from pypath.inputs_v2 import pfocr

    out = extract(
        'pfocr',
        'association_to_entity',
        {'identifier': '12', 'pfocr_id': 'PMC123__fig1', 'data_type': 'gene'},
    )
    assert out.relations[0].predicate == 'mentions'
    assert (
        pfocr.association_to_entity(
            {'identifier': 'unknown', 'pfocr_id': 'F1', 'data_type': 'chemical'}
        )
        is None
    )


def test_gene_classification_not_inferred_as_location():
    out = extract(
        'slctables',
        'schema',
        {
            'SLC name': 'SLC1A1',
            'Tissue and cellular expression': 'brain',
            'Substrates': 'L-Glu',
        },
    )
    assert not out.relations
    assert not any(
        a['term'] == 'description'
        for e in out.entities.values()
        for a in e.annotations
    )


def test_unknown_reaction_role_is_not_partonomy():
    out = extract(
        'rhea',
        'reactions_schema',
        {
            'rhea_id': '123',
            'participant_chebi': '1',
            'participant_role': 'unrecognized',
        },
    )
    assert not out.relations


def test_kegg_pubchem_ids_are_substances_not_compounds():
    out = extract(
        'kegg',
        'reactions_schema',
        {
            'reaction_id': 'R1',
            'reactant_kegg_id': 'C00001',
            'reactant_pubchem': '123',
        },
    )
    assert any(
        any(
            i['ns'] == 'pubchem_substance' and i['id'] == '123'
            for i in e.identifiers
        )
        for e in out.entities.values()
    )
    assert not any(
        any(i['ns'] == 'pubchem' and i['id'] == '123' for i in e.identifiers)
        for e in out.entities.values()
    )


def test_brenda_obsolete_term_is_excluded():
    from pypath.inputs_v2.brenda import enzyme_ontology_schema

    assert (
        enzyme_ontology_schema(
            {'id': 'EC:1.1.1.1', 'name': 'Old term', 'is_obsolete': 'true'}
        )
        is None
    )


def test_wikipathways_missing_taxonomy_and_version_identifier():
    out = extract(
        'wikipathways',
        'pathways_schema',
        {
            'pathway_id': 'WP1',
            'pathway_version_id': 'WP1_r1',
            'title': 'A pathway',
        },
    )
    assert len(out.entities) == 1
    out = extract(
        'wikipathways',
        'interactions_schema',
        {
            'source_entity_type': 'protein',
            'target_entity_type': 'protein',
            'source_label': 'A',
            'target_label': 'B',
            'interaction_types': 'Interaction',
            'taxon_id': '',
        },
    )
    assert len(out.entities) == 2
    assert all(e.taxon is None for e in out.entities.values())


def test_reactome_parser_uses_source_type_tokens():
    from rdflib import Graph, RDF, URIRef, Literal
    from pypath.inputs_v2.parsers import reactome

    graph = Graph()
    reaction = URIRef('urn:test:reaction')
    graph.add((reaction, RDF.type, reactome.BP.BiochemicalReaction))
    graph.add((reaction, reactome.BP.displayName, Literal('Source reaction')))
    rows = list(reactome._iterate_reactions(graph, {}, {}, {}))
    assert len(rows) == 1
    assert rows[0]['entity_type'] == 'reaction'


def test_missing_middle_reaction_role_preserves_alignment():
    out = extract(
        'rhea',
        'reactions_schema',
        {
            'rhea_id': '123',
            'participant_chebi': '1||2||3',
            'participant_role': 'reactant||unknown||product',
        },
    )
    edges = [
        (r.predicate, out.entities[r.object_entity_key].identifier)
        for r in out.relations
    ]
    assert edges == [('has_input', '1'), ('has_output', '3')]


def test_macdb_measurements_preserve_arm_and_comparator_without_inventing_units():
    out = extract(
        'macdb',
        'association_schema',
        {
            'pubchem_CID': '1',
            'trait_id': '1',
            'case_concentration': '<=1.5',
            'control_concentration': '2.5',
            'case_confidence_interval': '1-2',
        },
    )
    quantities = {
        a['quantity']['source_field']: a['quantity']
        for a in out.relations[0].annotations
        if a.get('quantity')
    }
    assert quantities['case_concentration']['has_numeric_value'] == 1.5
    assert quantities['case_concentration']['comparator'] == '<='
    assert quantities['case_concentration']['has_unit'] is None
    assert quantities['control_concentration']['has_numeric_value'] == 2.5
    assert 'case_confidence_interval' not in quantities


def test_intact_malformed_accessions_not_reassigned_to_genbank():
    from pypath.inputs_v2.intact import _parse_identifier_pairs

    assert (
        _parse_identifier_pairs(
            'ensembl:not-an-accession|refseq:unknown|entrezgene/locuslink:ABC'
        )
        == []
    )


def test_names_and_synonyms_are_supplemental_identifiers():
    from pypath.inputs_v2 import reactome

    record = reactome.reactions_schema(
        {
            'entity_type': 'reaction',
            'reactome_stable_id': 'R-HSA-123',
            'display_name': 'Example reaction',
            'synonyms': 'Alias A;Alias B',
        }
    )
    assert record.identifiers[0].type == 'reactome'
    pairs = {(str(i.type), i.value) for i in record.identifiers}
    assert ('name', 'Example reaction') in pairs
    assert ('synonym', 'Alias A') in pairs
    assert ('synonym', 'Alias B') in pairs
    out = extract(
        'reactome',
        'reactions_schema',
        {
            'entity_type': 'reaction',
            'reactome_stable_id': 'R-HSA-123',
            'display_name': 'Example reaction',
            'synonyms': 'Alias A;Alias B',
        },
    )
    assert not any(
        a['term'] in {'name', 'synonym'}
        for e in out.entities.values()
        for a in e.annotations
    )


def test_pfocr_caption_is_name_identifier_after_figure_id():
    from pypath.inputs_v2.pfocr import association_to_entity

    record = association_to_entity(
        {'identifier': '12', 'pfocr_id': 'PMC123__fig1', 'caption': 'Caption'}
    )
    assert record.subject.identifiers[0].type == 'pfocr'
    assert any(
        i.type == 'name' and i.value == 'Caption'
        for i in record.subject.identifiers
    )
    assert all(
        str(a.term) not in {'name', 'synonym'}
        for a in record.subject.annotations
    )
