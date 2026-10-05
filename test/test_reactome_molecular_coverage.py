"""BioPAX physical state stays with its own member or controller occurrence."""

import json

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF
from pypath.inputs_v2 import reactome
from pypath.inputs_v2.parsers import reactome as parser

BP = parser.BP


def uri(name):
    return URIRef('https://reactome.org/fixture/' + name)


def add_product(graph, name, accession, *, rna=False, sequence=None):
    physical, reference, xref = (
        uri(name + suffix) for suffix in ('', '-ref', '-xref')
    )
    graph.add((physical, RDF.type, BP.Rna if rna else BP.Protein))
    graph.add((physical, BP.entityReference, reference))
    graph.add((physical, BP.displayName, Literal(name)))
    graph.add(
        (reference, RDF.type, BP.RnaReference if rna else BP.ProteinReference)
    )
    graph.add((reference, BP.xref, xref))
    graph.add((xref, RDF.type, BP.UnificationXref))
    graph.add((xref, BP.db, Literal('RefSeq' if rna else 'UniProt')))
    graph.add((xref, BP.id, Literal(accession)))
    if sequence:
        graph.add((reference, BP.sequence, Literal(sequence)))
    return physical


def add_feature(
    graph,
    physical,
    name,
    *,
    present=True,
    modification=True,
    position=15,
    exact=True,
):
    feature, vocabulary, location = (
        uri(name + suffix) for suffix in ('', '-vocab', '-site')
    )
    graph.add((physical, BP.feature if present else BP.notFeature, feature))
    graph.add(
        (
            feature,
            RDF.type,
            BP.ModificationFeature if modification else BP.FragmentFeature,
        )
    )
    graph.add((feature, BP.featureLocation, location))
    graph.add((location, BP.sequencePosition, Literal(position)))
    graph.add(
        (
            location,
            BP.positionStatus,
            Literal('EQUAL' if exact else 'LESS-THAN'),
        )
    )
    if modification:
        graph.add((feature, BP.modificationType, vocabulary))
        graph.add((vocabulary, BP['term'], Literal('phosphorylated residue')))
    return feature


def indexes(graph):
    xrefs = parser._build_xref_cache(graph, BP)
    return xrefs, parser._load_entity_reference_index(graph, xrefs)


def add_control(graph, physical):
    control, reaction = uri('control'), uri('reaction')
    graph.add((control, RDF.type, BP.Catalysis))
    graph.add((control, BP.controller, physical))
    graph.add((control, BP.controlled, reaction))
    graph.add((control, BP.controlType, Literal('ACTIVATION')))
    graph.add((reaction, RDF.type, BP.BiochemicalReaction))
    graph.add((reaction, BP.displayName, Literal('source reaction')))
    return control, reaction


def test_controller_state_survives_control_parser_and_mapper():
    graph = Graph()
    physical = add_product(graph, 'controller', 'P04637-2', sequence='MAAK')
    add_feature(graph, physical, 'phosphosite')
    add_control(graph, physical)
    xrefs, references = indexes(graph)
    (row,) = parser._iterate_controls(graph, xrefs, references, {})
    entity = reactome.controller_builder(row)
    assert entity.molecular_form['isoform_identifier']['id'] == 'P04637-2'
    (feature,) = entity.molecular_form['modifications']
    assert feature['position'] == 15
    assert (
        feature['coordinate_reference']['identifier']['ns']
        == 'protein_sequence_sha256'
    )


def test_group_members_keep_distinct_states_without_parent_broadcast():
    graph = Graph()
    group = uri('group')
    graph.add((group, RDF.type, BP.PhysicalEntity))
    graph.add((group, BP.displayName, Literal('alternative states')))
    for number in (2, 3):
        physical = add_product(graph, str(number), 'P04637-' + str(number))
        add_feature(graph, physical, 'feature-' + str(number), position=number)
        graph.add((group, BP.memberPhysicalEntity, physical))
    add_control(graph, group)
    xrefs, references = indexes(graph)
    (row,) = parser._iterate_control_groups(graph, xrefs, references, {})
    mapped = reactome.control_groups_schema(row)
    assert mapped.molecular_form is None
    assert [
        member.member.molecular_form['isoform_identifier']['id']
        for member in mapped.membership
    ] == ['P04637-2', 'P04637-3']
    assert [
        member.member.molecular_form['modifications'][0]['position']
        for member in mapped.membership
    ] == [2, 3]
    (independent,) = parser._iterate_physical_groups(graph, xrefs, references)
    assert len(reactome.control_groups_schema(independent).membership) == 2


def test_negative_and_nonmodification_features_preserve_source_context_only():
    graph = Graph()
    physical = add_product(graph, 'p', 'P04637-2')
    add_feature(graph, physical, 'negative', present=False, exact=False)
    add_feature(graph, physical, 'fragment', modification=False)
    xrefs, references = indexes(graph)
    participant = parser._extract_participant_data(
        graph, physical, 'reactant', references, xrefs, {}
    )
    assert participant['molecular_form']['modifications'] is None
    contexts = participant['feature_context']
    assert {feature['present'] for feature in contexts} == {True, False}
    negative = next(feature for feature in contexts if not feature['present'])
    assert negative['locations'][0]['properties'][str(BP.positionStatus)] == [
        'LESS-THAN'
    ]
    row = parser._flatten_participants([participant])
    row.update(
        {
            'reactome_stable_id': 'R-HSA-1',
            'display_name': 'reaction',
            'entity_type': 'reaction',
        }
    )
    mapped = reactome.reactions_schema(row)
    context = json.loads(
        next(
            annotation.value
            for annotation in mapped.membership[0].annotations
            if annotation.term == 'biopax:feature_context'
        )
    )
    assert len(context) == 2


def test_transcript_reference_and_translation_are_not_dropped():
    graph = Graph()
    transcript = add_product(
        graph, 'rna', 'NM_000546.6', rna=True, sequence='ACGU'
    )
    protein = add_product(graph, 'protein', 'P04637')
    reaction = uri('translation')
    graph.add((reaction, RDF.type, BP.BiochemicalReaction))
    graph.add((reaction, BP.left, transcript))
    graph.add((reaction, BP.right, protein))
    graph.add((reaction, BP.displayName, Literal('translation')))
    xrefs, references = indexes(graph)
    (row,) = parser._iterate_reactions(graph, xrefs, references, {})
    mapped = reactome.reactions_schema(row)
    transcripts = [
        member.member
        for member in mapped.membership
        if member.member.molecular_form
        and any(
            identifier['id'] == 'NM_000546.6'
            for identifier in member.member.molecular_form[
                'sequence_identifiers'
            ]
            or []
        )
    ]
    assert len(transcripts) == 1
    assert str(transcripts[0].type) == 'rna_product'


def test_template_reaction_retains_transcript_and_product_form():
    graph = Graph()
    transcript = add_product(graph, 'template-rna', 'NM_000546.6', rna=True)
    protein = add_product(graph, 'product-protein', 'P04637-2')
    reaction = uri('template-reaction')
    graph.add((reaction, RDF.type, BP.TemplateReaction))
    graph.add((reaction, BP.template, transcript))
    graph.add((reaction, BP.product, protein))
    graph.add((reaction, BP.displayName, Literal('source translation')))
    graph.add((reaction, BP.templateDirection, Literal('FORWARD')))
    xrefs, references = indexes(graph)
    row, = parser._iterate_reactions(graph, xrefs, references, {})
    record = reactome.reactions_schema(row)
    assert len(record.membership) == 2
    assert record.membership[0].member.molecular_form['sequence_identifiers'][0]['id'] == 'NM_000546.6'
    assert record.membership[1].member.molecular_form['isoform_identifier']['id'] == 'P04637-2'
    assert any(annotation.term == 'biopax:templateDirection' and annotation.value == 'FORWARD' for annotation in record.annotations)


def test_dependent_controllers_are_one_source_group_with_separate_forms():
    graph = Graph()
    first = add_product(graph, 'first-controller', 'P04637-2')
    second = add_product(graph, 'second-controller', 'P04637-3')
    control, reaction = add_control(graph, first)
    graph.add((control, BP.controller, second))
    other_reaction = uri('second-reaction')
    graph.add((other_reaction, RDF.type, BP.BiochemicalReaction))
    graph.add((other_reaction, BP.displayName, Literal('another source reaction')))
    graph.add((control, BP.controlled, other_reaction))
    xrefs, references = indexes(graph)
    rows = list(parser._iterate_controls(graph, xrefs, references, {}))
    assert len(rows) == 2
    assert {row['controller_entity_type'] for row in rows} == {'logical_control_set'}
    assert reactome.controller_builder(rows[0]).molecular_form is None
    groups = list(parser._iterate_physical_groups(graph, xrefs, references))
    assert len(groups) == 1
    group = reactome.control_groups_schema(groups[0])
    assert {member.member.molecular_form['isoform_identifier']['id'] for member in group.membership} == {'P04637-2', 'P04637-3'}
    assert json.loads(groups[0]['controller_control_set'])['logic'] == 'AND'


def test_complex_components_keep_their_own_controller_forms():
    graph = Graph()
    complex_ = uri('complex')
    graph.add((complex_, RDF.type, BP.Complex))
    graph.add((complex_, BP.displayName, Literal('explicit assembly')))
    for number in (2, 3):
        member = add_product(graph, str(number), 'P04637-' + str(number))
        graph.add((complex_, BP.component, member))
    add_control(graph, complex_)
    xrefs, references = indexes(graph)
    row, = parser._iterate_control_groups(graph, xrefs, references, {})
    mapped = reactome.control_groups_schema(row)
    assert mapped.type == 'macromolecular_complex'
    assert mapped.molecular_form is None
    assert len(mapped.membership) == 2
