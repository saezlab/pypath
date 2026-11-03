"""Source contract of the talklr ligand-receptor input."""

import io
from types import SimpleNamespace

import pandas as pd
import pytest
import rdata
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2.parsers.talklr import COLUMNS, iter_pairs
from pypath.inputs_v2.talklr import resource


def row(*values):
    return dict(zip(COLUMNS, values, strict=True))


# Rows of the ``receptor_ligand`` data frame in receptor_ligand.rda.
DATABASE_ROW = row(
    'A2M_LRP1', 'A2M', 'alpha-2-macroglobulin', 'LRP1',
    'low density lipoprotein receptor-related protein 1', None, 'HPMR', None,
    'HPRD', 'STRING.binding', 'STRING.experiment', 'A2M', 'LRP1', None,
    'known', 'literature supported',
)
# Supported by the curators' own reading only; the second PMID has a leading
# space in the source.
MANUAL_PMID_ROW = row(
    'FAT4_DCHS1', 'FAT4', 'FAT atypical cadherin 4', 'DCHS1',
    'dachsous cadherin-related 1', None, None, None, None, 'STRING.binding',
    None, None, None, '24056717, 21303848', 'novel', 'literature supported',
)
# The PMID cell holds a fragment of a PMID.
FRAGMENT_PMID_ROW = row(
    'IL12A_IL12RB2', 'IL12A',
    'interleukin 12A (natural killer cell stimulatory factor 1, cytotoxic '
    'lymphocyte maturation factor 1, p35)',
    'IL12RB2', 'interleukin 12 receptor, beta 2', 'DLRP', None, None, None,
    'STRING.binding', None, 'IL12A', 'IL12RB2', '943050', 'known',
    'literature supported',
)
PUTATIVE_ROW = row(
    'ANGPT1_ITGB1', 'ANGPT1', 'angiopoietin 1', 'ITGB1',
    'integrin, beta 1 (fibronectin receptor, beta polypeptide, antigen CD29 '
    'includes MDF2, MSK12)',
    None, None, None, None, 'STRING.binding', None, 'ANGPT1', 'ITGB1', None,
    'novel', 'putative',
)
# The pair name has the approved spelling; the symbol column is upper case.
UPPER_CASE_SYMBOL_ROW = row(
    'C1orf112_PLSCR4', 'C1ORF112', 'chromosome 1 open reading frame 112',
    'PLSCR4', 'phospholipid scramblase 4', None, None, None, None,
    'STRING.binding', None, None, None, None, 'novel', 'putative',
)
# The HPMR columns hold that database's symbol, here an older alias.
HPMR_ALIAS_ROW = row(
    'ADM_GPR182', 'ADM', 'adrenomedullin', 'GPR182',
    'G protein-coupled receptor 182', None, 'HPMR', None, None, None, None,
    'ADM', 'ADMR', None, 'known', 'literature supported',
)


def extract(source_row, dataset='interactions'):
    record = resource.datasets()[dataset].mapper(source_row)
    assert record is not None
    out = SilverExtractor('talklr', dataset)
    out.process_record(record, source_row, 'fixture:1', 1)
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
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def test_datasets_separate_literature_supported_and_putative_pairs():
    assert set(resource.datasets()) == {
        'interactions',
        'putative_interactions',
    }


def test_pair_is_ligand_to_receptor_with_its_databases():
    relation, ligand, receptor = extract(DATABASE_ROW)
    assert relation.predicate == 'interacts_with'
    assert relation.upstream_id == 'A2M_LRP1'
    assert identifiers(ligand) == {
        (Namespace.GENESYMBOL, 'A2M'),
        (Namespace.NAME, 'alpha-2-macroglobulin'),
    }
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'LRP1'),
        (
            Namespace.NAME,
            'low density lipoprotein receptor-related protein 1',
        ),
    }
    assert (ligand.namespace, ligand.identifier) == ('genesymbol', 'A2M')
    assert {ligand.entity_type, receptor.entity_type} == {'protein'}
    assert {ligand.taxon, receptor.taxon} == {'9606'}
    # The two STRING columns name one database.
    assert annotations(relation) == [
        ('supporting_data_source', 'HPMR'),
        ('supporting_data_source', 'HPRD'),
        ('supporting_data_source', 'STRING'),
    ]


def test_manual_pmids_are_the_publications_of_the_pair():
    relation, _, _ = extract(MANUAL_PMID_ROW)
    assert annotations(relation) == [
        ('publications', 'PMID:21303848'),
        ('publications', 'PMID:24056717'),
        ('supporting_data_source', 'STRING'),
    ]


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    source_row = dict(MANUAL_PMID_ROW, **{'PMID.Manual': '0123456, 123456789, 21303848'})
    relation, _, _ = extract(source_row)
    assert [v for t, v in annotations(relation) if t == 'publications'] == [
        'PMID:21303848'
    ]


def test_fragment_of_a_pmid_is_not_a_publication():
    relation, _, _ = extract(FRAGMENT_PMID_ROW)
    assert annotations(relation) == [
        ('supporting_data_source', 'DLRP'),
        ('supporting_data_source', 'STRING'),
    ]


def test_putative_pair_maps_like_any_other_in_its_own_dataset():
    relation, ligand, receptor = extract(PUTATIVE_ROW, 'putative_interactions')
    assert relation.predicate == 'interacts_with'
    assert relation.dataset == 'putative_interactions'
    assert (ligand.identifier, receptor.identifier) == ('ANGPT1', 'ITGB1')
    assert annotations(relation) == [('supporting_data_source', 'STRING')]


def test_symbol_is_kept_as_the_source_spells_it():
    relation, ligand, _ = extract(UPPER_CASE_SYMBOL_ROW, 'putative_interactions')
    assert relation.upstream_id == 'C1orf112_PLSCR4'
    assert (Namespace.GENESYMBOL, 'C1ORF112') in identifiers(ligand)
    assert (Namespace.GENESYMBOL, 'C1orf112') not in identifiers(ligand)


def test_pair_asserts_no_direction_role_or_curation_status():
    relation, _, _ = extract(DATABASE_ROW)
    # Databases and the participants' species only: ``known``, ``novel``,
    # ``literature supported`` and ``putative`` have no faithful slot, and a
    # ligand-receptor pair states no effect on the receptor.
    assert {a['term'] for a in relation.annotations} == {
        'in_taxon',
        'supporting_data_source',
        PARTICIPANT_ROLE,
    }
    assert {
        a['term'] for a in relation.annotations if a['scope'] == 'relation'
    } == {'supporting_data_source'}
    # The roles are annotations of the participants, not of the relation.
    assert annotations(relation, 'subject') == [
        (PARTICIPANT_ROLE, 'ligand'),
        ('in_taxon', 'NCBITaxon:9606'),
    ]
    assert annotations(relation, 'object') == [
        (PARTICIPANT_ROLE, 'receptor'),
        ('in_taxon', 'NCBITaxon:9606'),
    ]


@pytest.mark.parametrize('dataset,source_row', [
    ('interactions', DATABASE_ROW),
    ('putative_interactions', PUTATIVE_ROW),
])
def test_ligand_and_receptor_columns_give_the_participant_roles(dataset, source_row):
    relation, ligand, receptor = extract(source_row, dataset)
    assert (Namespace.GENESYMBOL, source_row['Ligand.ApprovedSymbol']) in identifiers(ligand)
    assert (Namespace.GENESYMBOL, source_row['Receptor.ApprovedSymbol']) in identifiers(receptor)
    assert [v for t, v in annotations(relation, 'subject') if t == PARTICIPANT_ROLE] == ['ligand']
    assert [v for t, v in annotations(relation, 'object') if t == PARTICIPANT_ROLE] == ['receptor']


def test_hpmr_alias_is_not_an_identifier():
    _, _, receptor = extract(HPMR_ALIAS_ROW)
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'GPR182'),
        (Namespace.NAME, 'G protein-coupled receptor 182'),
    }


def test_reader_rejects_an_unknown_database_flag_instead_of_dropping_it(tmp_path):
    # A flag that is not the column's name used to vanish from the relation.
    for column, flag in (('HPRD', 'yes'), ('STRING.binding', 'STRING')):
        rows = [list(dict(DATABASE_ROW, **{column: flag}).values())]
        opener = rda(tmp_path, rows, list(COLUMNS))
        with pytest.raises(ValueError, match='Unknown talklr database flag'):
            list(iter_pairs(opener, evidence='literature supported'))
    # The flags of the real rows pass.
    rows = [list(DATABASE_ROW.values())]
    opener = rda(tmp_path, rows, list(COLUMNS))
    assert list(iter_pairs(opener, evidence='literature supported')) == [
        DATABASE_ROW
    ]


@pytest.mark.parametrize(
    ('column', 'database'),
    [
        ('DLRP', 'DLRP'),
        ('HPMR', 'HPMR'),
        ('IUPHAR', 'IUPHAR'),
        ('HPRD', 'HPRD'),
        ('STRING.binding', 'STRING'),
        ('STRING.experiment', 'STRING'),
    ],
)
def test_each_database_column_names_its_database(column, database):
    flags = {c: None for c in ('DLRP', 'HPMR', 'IUPHAR', 'HPRD', 'STRING.binding', 'STRING.experiment')}
    relation, _, _ = extract({**DATABASE_ROW, **flags, column: column})
    assert annotations(relation) == [('supporting_data_source', database)]


@pytest.mark.parametrize('role', ['Ligand', 'Receptor'])
def test_no_relation_without_an_approved_symbol(role):
    # The name alone must not become a human protein.
    mapper = resource.datasets()['interactions'].mapper
    assert mapper(dict(DATABASE_ROW, **{f'{role}.ApprovedSymbol': None})) is None


def rda(tmp_path, rows, columns, name='receptor_ligand'):
    """An opener over an .rda file holding one named data frame."""
    path = tmp_path / 'receptor_ligand.rda'
    frame = pd.DataFrame(rows, columns=columns, dtype=object)
    rdata.write_rda(path, {name: frame})
    return SimpleNamespace(result=io.BytesIO(path.read_bytes()))


def test_reader_selects_the_rows_of_one_verdict(tmp_path):
    rows = [DATABASE_ROW, PUTATIVE_ROW, MANUAL_PMID_ROW]
    values = [list(r.values()) for r in rows]
    opener = rda(tmp_path, values, list(COLUMNS))
    assert list(iter_pairs(opener, evidence='literature supported')) == [
        DATABASE_ROW,
        MANUAL_PMID_ROW,
    ]
    opener = rda(tmp_path, values, list(COLUMNS))
    assert list(iter_pairs(opener, evidence='putative')) == [PUTATIVE_ROW]


def test_reader_rejects_an_unknown_verdict(tmp_path):
    rows = [list(dict(DATABASE_ROW, **{'Pair.Evidence': 'EXCLUDED'}).values())]
    opener = rda(tmp_path, rows, list(COLUMNS))
    with pytest.raises(ValueError, match='EXCLUDED'):
        list(iter_pairs(opener, evidence='literature supported'))


def test_reader_rejects_a_table_with_changed_columns(tmp_path):
    opener = rda(tmp_path, [['A2M', 'LRP1']], ['ligand', 'receptor'])
    with pytest.raises(ValueError, match='Pair.Name'):
        list(iter_pairs(opener, evidence='putative'))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_pairs(SimpleNamespace(result=None), evidence='putative'))
