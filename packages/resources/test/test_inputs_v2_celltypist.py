"""Source contract of the CellTypist encyclopedia marker input."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor

from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2.celltypist import map_marker, resource
from pypath.inputs_v2.parsers.celltypist import COLUMNS, iter_markers

# Rows of encyclopedia_table.xlsx (worksheet Main).
ROWS = [
    (
        'B cells', 'Follicular B cells',
        'resting mature B lymphocytes found in the primary and secondary lymphoid follicles and participating in T cell-dependent immune responses',
        'CL:0000843', 'Intestine', 'Smillie et al. 2019', 'CXCR5, TNFRSF13B, CD22',
        'MT-RNR2, IGHA1, JCHAIN, GAS5, RPL41, BX640514.2, AC108879.1, AL512306.2, SYCP1, AL731577.2',
    ),
    (
        'B cells', 'Proliferative germinal center B cells',
        'proliferating germinal center B cells',
        'CL:0000844', 'Intestine, Lung-draining lymph node, Mesenteric lymph node, Spleen',
        'Dominguez Conde et al. 2022, Madissoon et al. 2020', 'MKI67, SUGCT, AICDA',
        'PCLAF, OSTN-AS1, GLYATL2, AC005614.1, TRIM71, TRBV11-3, TRBV10-2, AC007368.1, CTAG2, AD000090.1',
    ),
    (
        'B cells', 'Germinal center B cells',
        'proliferating mature B cells that undergo somatic hypermutation and class-switch recombination in secondary lymphoid organs',
        'CL:0000844', 'Intestine, Lung-draining lymph node, Mesenteric lymph node, Spleen',
        'Dominguez Conde et al. 2022, Smillie et al. 2019', 'POU2AF1, CD40, SUGCT',
        'AC022126.1, CTAG2, LINC02029, OR2A25, RGS13, SELE, KIAA0087, TRBV11-3, AL136456.1, AC007878.1',
    ),
    (
        'Macrophages', 'Erythrophagocytic macrophages',
        'erythrophagocytosis-related macrophages including red pulp macrophages and Kupffer cells that mainly populate the spleen and liver',
        'CL:0000235', 'Bone marrow, Intestine, Liver, Lung, Lung-draining lymph node, Mesenteric lymph node, Omentum, Spleen',
        'Dominguez Conde et al. 2022', 'CD5L, SPIC, SCL40A1',
        'SELENOP, IFI30, CCL24, AC110995.1, AL512306.2, AP001059.6, AL109615.2, C1QC, AC022126.1, OR2A25',
    ),
]


def workbook(*rows, columns=COLUMNS, sheet='Main'):
    book = openpyxl.Workbook()
    book.active.title = sheet
    book.active.append(list(columns))
    for row in rows:
        book.active.append(list(row))
    handle = io.BytesIO()
    book.save(handle)
    return SimpleNamespace(result=io.BytesIO(handle.getvalue()))


def markers(*rows, **kwargs):
    return list(iter_markers(workbook(*rows, **kwargs)))


def extract(row):
    record = map_marker(row)
    assert record is not None
    out = SilverExtractor('celltypist', 'markers')
    out.process_record(record, row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def values(relation, term, scope='relation'):
    return [a['value'] for a in relation.annotations if a['term'] == term and a['scope'] == scope]


def test_parser_gives_a_row_for_each_gene_of_the_two_columns():
    rows = markers(ROWS[0])
    assert [r['gene'] for r in rows if r['marker_source'] == 'curated_marker'] == [
        'CXCR5', 'TNFRSF13B', 'CD22',
    ]
    assert len([r for r in rows if r['marker_source'] == 'model_top_gene']) == 10
    first = rows[0]
    assert first['cell_type'] == 'Follicular B cells'
    assert first['cell_type_group'] == 'B cells'
    assert first['cell_ontology_id'] == 'CL:0000843'
    assert first['tissues'] == ['Intestine']
    assert first['datasets'] == ['Smillie et al. 2019']


def test_cell_type_is_the_cell_ontology_term():
    row = markers(ROWS[0])[0]
    relation, gene, cell = extract(row)
    assert (cell.entity_type, cell.namespace, cell.identifier) == ('ontology_class', 'cl', 'CL:0000843')
    assert (gene.entity_type, gene.namespace, gene.identifier, gene.taxon) == (
        'gene', 'genesymbol', 'CXCR5', '9606',
    )
    # The term is the Cell Ontology's own: no name of the table is given to it.
    assert cell.label is None
    assert relation.predicate == 'expressed_in'


def test_curated_markers_are_expressed_in_and_model_genes_are_associated():
    rows = markers(ROWS[0])
    predicates = {r['marker_source']: extract(r)[0].predicate for r in rows}
    assert predicates == {'curated_marker': 'expressed_in', 'model_top_gene': 'associated_with'}
    model = next(r for r in rows if r['gene'] == 'MT-RNR2')
    relation, _, _ = extract(model)
    assert values(relation, 'celltypist:marker_source') == ['model_top_gene']


def test_labels_tissues_and_datasets_are_attributes_of_the_statement():
    rows = markers(ROWS[1])
    relation, _, _ = extract(rows[0])
    assert values(relation, 'celltypist:cell_type') == ['Proliferative germinal center B cells']
    assert values(relation, 'celltypist:cell_type_group') == ['B cells']
    assert values(relation, 'celltypist:tissue') == [
        'Intestine', 'Lung-draining lymph node', 'Mesenteric lymph node', 'Spleen',
    ]
    assert values(relation, 'celltypist:dataset') == [
        'Dominguez Conde et al. 2022', 'Madissoon et al. 2020',
    ]
    assert values(relation, 'celltypist:cell_type_description') == ['proliferating germinal center B cells']


def test_two_cell_types_of_one_ontology_term_keep_their_labels():
    first = markers(ROWS[1])[0]
    second = markers(ROWS[2])[0]
    r1, _, c1 = extract(first)
    r2, _, c2 = extract(second)
    assert c1.identifier == c2.identifier == 'CL:0000844'
    assert values(r1, 'celltypist:cell_type') != values(r2, 'celltypist:cell_type')
    # SUGCT is a curated marker of both and gives a statement for each.
    sugct = [r for r in markers(ROWS[1], ROWS[2]) if r['gene'] == 'SUGCT']
    assert len(sugct) == 2


def test_gene_symbols_are_kept_as_written():
    # SCL40A1 is the table's spelling of SLC40A1.
    rows = markers(ROWS[3])
    assert 'SCL40A1' in {r['gene'] for r in rows if r['marker_source'] == 'curated_marker'}
    assert 'SLC40A1' not in {r['gene'] for r in rows if r['marker_source'] == 'curated_marker'}


def test_cell_type_without_a_cell_ontology_id_raises():
    for missing in (None, '', 'B cell', 'CL:843'):
        with pytest.raises(ValueError, match='Cell Ontology ID'):
            markers(ROWS[0][:3] + (missing,) + ROWS[0][4:])


def test_parser_rejects_a_changed_or_empty_table():
    with pytest.raises(ValueError, match='lacks expected columns'):
        markers(ROWS[0], columns=[c.replace('Curated markers', 'Markers') for c in COLUMNS])
    with pytest.raises(ValueError, match='no rows'):
        markers()


def test_unknown_marker_source_raises():
    with pytest.raises(ValueError, match='marker source'):
        map_marker(dict(markers(ROWS[0])[0], marker_source='other'))


def test_nothing_else_is_asserted():
    relation, _, _ = extract(markers(ROWS[0])[0])
    assert {a['term'] for a in relation.annotations} == {
        'celltypist:marker_source',
        'celltypist:cell_type',
        'celltypist:cell_type_group',
        'celltypist:cell_type_description',
        'celltypist:tissue',
        'celltypist:dataset',
        'in_taxon',
    }


def test_resource_declares_its_dataset_and_license():
    from pypath.inputs_v2 import celltypist

    assert celltypist.config.id == ResourceCv.CELLTYPIST
    assert celltypist.config.license == LicenseCV.CC_BY_NC_ND_4_0
    assert set(resource.datasets()) == {'markers'}
