"""Source contract of the PanglaoDB cell type marker input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor

from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2.panglaodb import (
    CELL_TYPE,
    GERM_LAYER,
    ORGAN,
    exclusion_reason,
    map_marker,
    map_term,
    resource,
)
from pypath.inputs_v2.parsers.panglaodb import iter_markers, iter_terms

HEADER = (
    'species\tofficial gene symbol\tcell type\tnicknames\tubiquitousness index\t'
    'product description\tgene type\tcanonical marker\tgerm layer\torgan\t'
    'sensitivity_human\tsensitivity_mouse\tspecificity_human\tspecificity_mouse'
)
# Lines of PanglaoDB_markers_27_Mar_2020.tsv.
LINES = [
    'Mm Hs\tCTRB1\tAcinar cells\tCTRB\t0.017\tchymotrypsinogen B1\tprotein-coding gene\t1\tEndoderm\tPancreas\t1.0\t0.957143\t0.000628931\t0.0159201',
    'Mm\tTRY4\tAcinar cells\tNA\t0.007\ttrypsin 4\tprotein coding gene\t1\tEndoderm\tPancreas\tNA\tNA\tNA\tNA',
    'Hs\tCTRB2\tAcinar cells\tNA\t0.001\tchymotrypsinogen B2\tprotein-coding gene\t1\tEndoderm\tPancreas\t0.916667\t0.0\t0.00125786\t0.0',
    'Mm Hs\tRSPH1\tEpendymal cells\tFLJ32753|RSP44|RSPH10A|CILD24|TSGA2\t0.021\tradial spoke head component 1\tprotein-coding gene\t0\tEctoderm\tBrain\t1.0\t0.0\t0.0144382\t0.0',
    'Mm Hs\tITGAM\tMyeloid-derived suppressor cells\tCD11b|CR3A|CD11B\t0.025\tintegrin subunit alpha M\tprotein-coding gene\tNA\tmesoderm\tImmune system\t0.0\t0.0\t0.00313283\t0.0302855',
    'Mm Hs\tALCAM \tCardiomyocytes\tCD166|MEMD\t0.0\tactivated leukocyte cell adhesion molecule\tprotein-coding gene\t1\tMesoderm\tHeart\tNA\tNA\tNA\tNA',
    'Mm Hs\tDRD2\tGlomus cells\tNA\t0.0\tdopamine receptor D2\tprotein-coding gene\t1\tNA\tNA\t0.0\t0.0\t0.0\t0.000577478',
    'None\tH2-T10\tT cells\tNA\t0.002\tNA\tNA\tNA\tMesoderm\tImmune system\t0.0\t0.0233645\t0.0\t0.00149632',
    '4\tDEPP1\tAlpha cells\tDEPP|Fseg|C10orf10\t0.0\tDEPP1, autophagy regulator\tprotein-coding gene\t1\tEndoderm\tPancreas\tNA\tNA\tNA\tNA',
    'Mm Hs\tXCL1\tDecidual cells\tLPTN|ATAC|SCM-1a|SCM-1|LTN|SCYC1\t0.006\tX-C motif chemokine ligand 1\tprotein-coding gene\tNA\tNA\tReproductive\t1.0\t0.0\t0.00971482\t0.0',
]


def table(*lines, header=HEADER):
    return SimpleNamespace(result=io.StringIO('\n'.join([header, *lines]) + '\n'))


def markers(*lines, header=HEADER):
    return list(iter_markers(table(*lines, header=header)))


def by_taxon(rows):
    return {row['taxon']: row for row in rows}


def extract(row, dataset='markers'):
    record = resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('panglaodb', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return out


def relation_of(row):
    out = extract(row)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def annotations(relation, scope='relation'):
    return {
        a['term']: a['quantity']['has_numeric_value'] if a['quantity'] else a['value']
        for a in relation.annotations
        if a['scope'] == scope
    }


def test_row_with_two_species_gives_a_statement_for_each():
    rows = markers(LINES[0])
    assert sorted(row['taxon'] for row in rows) == ['10090', '9606']
    human, mouse = by_taxon(rows)['9606'], by_taxon(rows)['10090']
    assert (human['sensitivity'], human['specificity']) == ('1.0', '0.000628931')
    assert (mouse['sensitivity'], mouse['specificity']) == ('0.957143', '0.0159201')
    assert human['species'] == mouse['species'] == 'Mm Hs'


def test_each_species_statement_carries_the_values_of_its_species():
    rows = by_taxon(markers(LINES[0]))
    human, mouse = rows['9606'], rows['10090']
    rel_h, gene_h, cell = relation_of(human)
    rel_m, gene_m, _ = relation_of(mouse)
    assert rel_h.predicate == rel_m.predicate == 'expressed_in'
    assert (gene_h.entity_type, gene_h.namespace, gene_h.identifier) == ('gene', 'genesymbol', 'CTRB1')
    assert gene_h.taxon == '9606' and gene_m.taxon == '10090'
    # One row, two genes: the human symbol as written, the mouse one in mouse form.
    assert (gene_m.entity_type, gene_m.namespace, gene_m.identifier) == ('gene', 'genesymbol', 'Ctrb1')
    assert (cell.entity_type, cell.namespace, cell.identifier) == (
        'ontology_class',
        CELL_TYPE,
        'Acinar cells',
    )
    assert annotations(rel_h)['panglaodb:sensitivity'] == 1.0
    assert annotations(rel_m)['panglaodb:sensitivity'] == 0.957143
    assert annotations(rel_h)['panglaodb:specificity'] == 0.000628931
    assert annotations(rel_m)['panglaodb:specificity'] == 0.0159201
    # The ubiquitousness index belongs to the gene and is the same.
    assert annotations(rel_h)['panglaodb:ubiquitousness_index'] == 0.017
    assert annotations(rel_m)['panglaodb:ubiquitousness_index'] == 0.017
    assert annotations(rel_h)['panglaodb:canonical_marker'] == 'true'


def test_mouse_only_marker_is_not_asserted_for_human():
    (row,) = markers(LINES[1])
    assert row['taxon'] == '10090'
    # The table has no values for the species of a marker that is not stated for it.
    assert row['sensitivity'] is None and row['specificity'] is None
    relation, gene, _ = relation_of(row)
    assert gene.taxon == '10090'
    assert 'panglaodb:sensitivity' not in annotations(relation)
    assert 'panglaodb:specificity' not in annotations(relation)
    # The table writes TRY4; the mouse gene symbol is Try4.
    assert row['gene'] == 'Try4'
    assert gene.identifier == 'Try4'
    (human_only,) = markers(LINES[2])
    assert human_only['taxon'] == '9606'
    assert relation_of(human_only)[1].taxon == '9606'


def test_canonical_marker_flag():
    (zero,) = markers(LINES[3])[:1]
    assert zero['canonical_marker'] == '0'
    assert annotations(relation_of(zero)[0])['panglaodb:canonical_marker'] == 'false'
    (na,) = markers(LINES[4])[:1]
    assert na['canonical_marker'] is None
    # NA is no statement about the canonical status.
    assert 'panglaodb:canonical_marker' not in annotations(relation_of(na)[0])
    with pytest.raises(ValueError, match='canonical marker flag'):
        map_marker(dict(zero, canonical_marker='2'))


def test_missing_values_give_no_attribute():
    rows = by_taxon(markers(LINES[5]))
    human = rows['9606']
    assert human['sensitivity'] is None
    relation, gene, _ = relation_of(human)
    assert set(annotations(relation)) == {
        'panglaodb:ubiquitousness_index',
        'panglaodb:canonical_marker',
    }
    # The trailing space of the symbol is a defect of the table.
    assert gene.identifier == 'ALCAM'


def test_germ_layer_written_in_lower_case_is_repaired():
    rows = markers(LINES[4])
    assert {r['germ_layer'] for r in rows} == {'Mesoderm'}


def test_rows_without_a_species_are_not_mapped():
    for line in (LINES[7], LINES[8]):
        (row,) = markers(line)
        assert row['taxon'] is None
        assert exclusion_reason(row) == 'species not stated'
        assert map_marker(row) is None
    # They are not asserted for human or mouse.
    assert [r['taxon'] for r in markers(*LINES[7:9])] == [None, None]


def test_unknown_species_flag_raises():
    with pytest.raises(ValueError, match='Unknown PanglaoDB species'):
        markers(LINES[0].replace('Mm Hs', 'Hs Rn', 1))


def test_parser_rejects_a_changed_or_empty_table():
    with pytest.raises(ValueError, match='lacks expected columns'):
        markers(LINES[0], header=HEADER.replace('organ', 'tissue'))
    with pytest.raises(ValueError, match='no rows'):
        markers()
    with pytest.raises(ValueError, match='no rows'):
        list(iter_terms(table()))


def test_terms_have_names_and_hierarchy():
    rows = list(iter_terms(table(*LINES)))
    by_kind = {}
    for row in rows:
        by_kind.setdefault(row['kind'], []).append(row)
    assert {r['id'] for r in by_kind['germ_layer']} == {'Endoderm', 'Ectoderm', 'Mesoderm'}
    assert {r['id'] for r in by_kind['organ']} == {
        'Pancreas', 'Brain', 'Immune system', 'Heart', 'Reproductive',
    }
    acinar = next(r for r in by_kind['cell_type'] if r['id'] == 'Acinar cells')
    assert (acinar['organ'], acinar['germ_layer']) == ('Pancreas', 'Endoderm')
    out = extract(acinar, 'terms')
    (cell,) = (e for e in out.entities.values() if e.label)
    assert (cell.entity_type, cell.namespace, cell.identifier, cell.label) == (
        'ontology_class', CELL_TYPE, 'Acinar cells', 'Acinar cells',
    )
    edges = {
        (r.predicate, out.entities[r.object_entity_key].namespace, out.entities[r.object_entity_key].identifier)
        for r in out.relations
    }
    assert edges == {
        ('part_of', ORGAN, 'Pancreas'),
        ('associated_with', GERM_LAYER, 'Endoderm'),
    }
    organ = extract(next(r for r in by_kind['organ'] if r['id'] == 'Pancreas'), 'terms')
    assert organ.relations == []


def test_cell_type_is_part_of_its_organ_and_not_a_subclass_of_it():
    # Acinar cells are a part of the pancreas, not a kind of pancreas.
    acinar = next(
        r for r in iter_terms(table(*LINES)) if r['kind'] == 'cell_type' and r['id'] == 'Acinar cells'
    )
    out = extract(acinar, 'terms')
    predicates = {r.predicate for r in out.relations}
    assert 'subclass_of' not in predicates
    (organ_edge,) = (
        r for r in out.relations if out.entities[r.object_entity_key].namespace == ORGAN
    )
    assert organ_edge.predicate == 'part_of'


def test_cell_type_without_organ_or_germ_layer_has_no_parent():
    rows = {r['id']: r for r in iter_terms(table(*LINES)) if r['kind'] == 'cell_type'}
    glomus = rows['Glomus cells']
    assert (glomus['organ'], glomus['germ_layer']) == (None, None)
    assert extract(glomus, 'terms').relations == []
    assert {r['kind'] for r in iter_terms(table(LINES[6]))} == {'cell_type'}


def test_cell_type_with_two_organs_raises():
    other = LINES[0].replace('Pancreas', 'Liver')
    with pytest.raises(ValueError, match='two organs or germ layers'):
        list(iter_terms(table(LINES[0], other)))


def test_statement_object_is_the_entity_of_the_terms_dataset():
    (row,) = markers(LINES[2])
    _, _, cell = relation_of(row)
    terms = {r['id']: r for r in iter_terms(table(*LINES)) if r['kind'] == 'cell_type'}
    out = extract(terms['Acinar cells'], 'terms')
    (term,) = (e for e in out.entities.values() if e.label)
    assert (cell.entity_type, cell.namespace, cell.identifier) == (
        term.entity_type, term.namespace, term.identifier,
    )


def test_unknown_term_kind_raises():
    with pytest.raises(ValueError, match='term kind'):
        map_term({'kind': 'tissue', 'id': 'Lung'})


def test_nothing_else_is_asserted():
    relation, _, _ = relation_of(markers(LINES[0])[0])
    assert {a['term'] for a in relation.annotations} == {
        'panglaodb:sensitivity',
        'panglaodb:specificity',
        'panglaodb:ubiquitousness_index',
        'panglaodb:canonical_marker',
        'in_taxon',
    }
    # Organ and germ layer are properties of the cell type term, not of the statement.
    assert 'Pancreas' not in {a['value'] for a in relation.annotations}


def test_resource_declares_its_datasets_and_license():
    from pypath.inputs_v2 import panglaodb

    assert panglaodb.config.id == ResourceCv.PANGLAODB
    assert panglaodb.config.license == LicenseCV.UNSPECIFIED
    assert set(resource.datasets()) == {'terms', 'markers'}


def test_mouse_symbol_form():
    from pypath.inputs_v2.parsers.panglaodb import mouse_symbol

    assert mouse_symbol('TRY4') == 'Try4'
    assert mouse_symbol('CD36') == 'Cd36'
    # Each part between hyphens: the H2 genes keep their capital.
    assert mouse_symbol('H2-M2') == 'H2-M2'

