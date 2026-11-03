"""Source contract of the LRdb ligand-receptor input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2.lrdb import DATABASES, resource
from pypath.inputs_v2.parsers.lrdb import COLUMNS, SOURCES, iter_pairs


def row(*values):
    return dict(zip(COLUMNS, values, strict=True))


# Rows of LRdb_122019.txt. Long alternative-name cells are shortened.
DATABASE_ROW = row(
    'A2M', 'LRP1', 'Alpha-2-macroglobulin',
    'Prolow-density lipoprotein receptor-related protein 1', 'CPAMD5,FWP007',
    'A2MR,APR', 'Alpha-2-M', 'LRP-1 || Alpha-2-macroglobulin receptor',
    'fantom5,HPMR,HPRD', '10652313', '', '', '',
)
LITERATURE_ROW = row(
    'FYN', 'SPN', 'FYN Proto-Oncogene, Src Family Tyrosine Kinase',
    'Sialophorin', 'SLK, SYN, P59-FYN', 'SPN, LSN,GPL115',
    'Proto-Oncogene C-Fyn|| Src/Yes-Related Novel', 'Leukosialin',
    'literature', '8910342', '', 'T cells', 'In vitro',
)
# The PMID cell was cut at 11 characters in the source.
TRUNCATED_PMID_ROW = row(
    'ADAM15', 'ITGA5',
    'Disintegrin and metalloproteinase domain-containing protein 15',
    'Integrin alpha-5', 'MDC15', 'FNRA', 'ADAM 15', 'CD49 antigen-like family member E',
    'fantom5,HPMR', '9914169,952', '', '', '',
)
NO_PMID_ROW = row(
    'ACTR2', 'LDLR', 'Actin-related protein 2',
    'Low-density lipoprotein receptor', 'ARP2', '', 'Actin-like protein 2',
    'LDL receptor', 'reactome', '', '', '', '',
)
CLEAVED_NAME_ROW = row(
    'GCG', 'GCGR',
    'Glucagon [Cleaved into: Glicentin; Glicentin-related polypeptide',
    'Glucagon receptor', '', '', 'GRPP); Oxyntomoduli || OXM', 'GL-R',
    'IUPHAR,fantom5,HPMR,HPRD', '12529935,4305077,8382034', '', '', '',
)
PARENTHESIS_PMID_ROW = row(
    'APOE', 'LRP5', 'Apolipoprotein E',
    'Low-density lipoprotein receptor-related protein 5', '', 'LR3,LRP7',
    'Apo-E', 'LRP-5', 'fantom5,HPMR', '9832610)', '', '', '',
)
SEMICOLON_SYNONYM_ROW = row(
    'CCL3L1', 'ACKR2', 'C-C motif chemokine 3-like 1',
    'Atypical chemokine receptor 2', 'D17S1718,G0S19-2,SCYA3L1;,CCL3L3',
    'CCBP2,CCR10,CMKBR9,D6', 'G0/G1 switch regulatory protein 19-2',
    'C-C chemokine receptor D6', 'fantom5,literature,HPRD', '10699158', '',
    'T lymphocyte', '',
)


def extract(source_row):
    record = resource.datasets()['interactions'].mapper(source_row)
    assert record is not None
    out = SilverExtractor('lrdb', 'interactions')
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


def test_pair_is_ligand_to_receptor_with_databases_and_publication():
    relation, ligand, receptor = extract(DATABASE_ROW)
    assert relation.predicate == 'interacts_with'
    assert identifiers(ligand) == {
        (Namespace.GENESYMBOL, 'A2M'),
        (Namespace.NAME, 'Alpha-2-macroglobulin'),
        (Namespace.SYNONYM, 'CPAMD5'),
        (Namespace.SYNONYM, 'FWP007'),
    }
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'LRP1'),
        (
            Namespace.NAME,
            'Prolow-density lipoprotein receptor-related protein 1',
        ),
        (Namespace.SYNONYM, 'A2MR'),
        (Namespace.SYNONYM, 'APR'),
    }
    assert (ligand.namespace, ligand.identifier) == ('genesymbol', 'A2M')
    assert {ligand.entity_type, receptor.entity_type} == {'protein'}
    assert {ligand.taxon, receptor.taxon} == {'9606'}
    assert annotations(relation) == [
        ('publications', 'PMID:10652313'),
        ('supporting_data_source', 'HPMR'),
        ('supporting_data_source', 'HPRD'),
        ('supporting_data_source', 'fantom5'),
    ]


def test_literature_is_not_a_supporting_database():
    relation, ligand, receptor = extract(LITERATURE_ROW)
    # The authors' own curation is carried by the PMID; the remark is prose.
    assert annotations(relation) == [
        ('description', 'In vitro'),
        ('publications', 'PMID:8910342'),
    ]
    assert (Namespace.SYNONYM, 'P59-FYN') in identifiers(ligand)
    assert (Namespace.SYNONYM, 'GPL115') in identifiers(receptor)


def test_pair_asserts_no_direction_role_or_cell_type():
    relation, _, _ = extract(LITERATURE_ROW)
    assert {a['term'] for a in relation.annotations} == {
        'description',
        'in_taxon',
        'publications',
        PARTICIPANT_ROLE,
    }
    assert {
        a['term'] for a in relation.annotations if a['scope'] == 'relation'
    } == {'description', 'publications'}
    # The receptor cell type (``T cells``) has no faithful slot.
    assert all('T cells' not in a['value'] for a in relation.annotations)
    # The roles are annotations of the participants, not of the relation.
    assert annotations(relation, 'subject') == [
        (PARTICIPANT_ROLE, 'ligand'),
        ('in_taxon', 'NCBITaxon:9606'),
    ]
    assert annotations(relation, 'object') == [
        (PARTICIPANT_ROLE, 'receptor'),
        ('in_taxon', 'NCBITaxon:9606'),
    ]


def test_ligand_and_receptor_columns_give_the_participant_roles():
    relation, ligand, receptor = extract(DATABASE_ROW)
    assert (Namespace.GENESYMBOL, 'A2M') in identifiers(ligand)
    assert (Namespace.GENESYMBOL, 'LRP1') in identifiers(receptor)
    assert [v for t, v in annotations(relation, 'subject') if t == PARTICIPANT_ROLE] == ['ligand']
    assert [v for t, v in annotations(relation, 'object') if t == PARTICIPANT_ROLE] == ['receptor']


def test_fragment_of_a_truncated_pmid_is_not_a_publication():
    relation, _, _ = extract(TRUNCATED_PMID_ROW)
    assert [v for t, v in annotations(relation) if t == 'publications'] == [
        'PMID:9914169'
    ]


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    source_row = dict(TRUNCATED_PMID_ROW, PMIDs='0123456,123456789,9914169')
    relation, _, _ = extract(source_row)
    assert [v for t, v in annotations(relation) if t == 'publications'] == [
        'PMID:9914169'
    ]


def test_pair_without_pmids_has_no_publication():
    relation, _, receptor = extract(NO_PMID_ROW)
    assert annotations(relation) == [('supporting_data_source', 'reactome')]
    # An empty synonym cell adds nothing.
    assert identifiers(receptor) == {
        (Namespace.GENESYMBOL, 'LDLR'),
        (Namespace.NAME, 'Low-density lipoprotein receptor'),
    }


def test_name_drops_the_truncated_cleavage_product_list():
    _, ligand, _ = extract(CLEAVED_NAME_ROW)
    assert identifiers(ligand) == {
        (Namespace.GENESYMBOL, 'GCG'),
        (Namespace.NAME, 'Glucagon'),
    }


def test_stray_characters_around_identifiers_are_removed():
    relation, _, _ = extract(PARENTHESIS_PMID_ROW)
    assert ('publications', 'PMID:9832610') in annotations(relation)
    relation, ligand, _ = extract(SEMICOLON_SYNONYM_ROW)
    assert {v for ns, v in identifiers(ligand) if ns == Namespace.SYNONYM} == {
        'D17S1718',
        'G0S19-2',
        'SCYA3L1',
        'CCL3L3',
    }
    # A pair from databases and the authors' curation keeps only the databases.
    assert [
        v for t, v in annotations(relation) if t == 'supporting_data_source'
    ] == ['HPRD', 'fantom5']


def test_remark_with_repeated_csv_quoting_is_read_as_prose():
    quotes = '"' * 32
    remark = f'{quotes}Crosstalk{quotes} between galectin-1, CD45 and CD3'
    relation, _, _ = extract(dict(LITERATURE_ROW, remarks=remark))
    assert ('description', '"Crosstalk" between galectin-1, CD45 and CD3') in (
        annotations(relation)
    )


def test_alternative_names_are_not_identifiers():
    _, ligand, receptor = extract(DATABASE_ROW)
    values = {value for _, value in identifiers(ligand) | identifiers(receptor)}
    assert 'Alpha-2-M' not in values
    assert 'Alpha-2-macroglobulin receptor' not in values


def test_reader_rejects_an_unknown_source_value_instead_of_dropping_it():
    # An unknown database used to vanish from the mapped relation.
    header = '\t'.join(COLUMNS)
    line = '\t'.join(dict(DATABASE_ROW, source='fantom5,predicted').values())
    with pytest.raises(ValueError, match="Unknown LRdb source 'predicted'"):
        list(iter_pairs(table([header, line])))
    # An empty cell is no known value either.
    line = '\t'.join(dict(DATABASE_ROW, source='').values())
    with pytest.raises(ValueError, match="Unknown LRdb source ''"):
        list(iter_pairs(table([header, line])))


def test_known_sources_are_the_databases_plus_literature():
    assert set(SOURCES) == set(DATABASES) | {'literature'}


def test_every_database_of_the_source_is_kept_as_the_source_spells_it():
    # The seven values of the ``source`` column that name a database; the
    # eighth value of the column, ``literature``, does not.
    source = 'fantom5,HPRD,HPMR,IUPHAR,reactome,uniprot,cellsignal.com,literature'
    relation, _, _ = extract(dict(DATABASE_ROW, source=source))
    assert sorted(
        v for t, v in annotations(relation) if t == 'supporting_data_source'
    ) == sorted(
        ['fantom5', 'HPRD', 'HPMR', 'IUPHAR', 'reactome', 'uniprot', 'cellsignal.com']
    )


@pytest.mark.parametrize('role', ['ligand', 'receptor'])
def test_no_relation_without_a_participant_symbol(role):
    # The name and synonyms alone do not identify the participant.
    mapper = resource.datasets()['interactions'].mapper
    assert mapper(dict(DATABASE_ROW, **{role: ''})) is None


def table(lines):
    return SimpleNamespace(result=io.StringIO('\n'.join(lines) + '\n'))


def test_reader_reads_the_quoted_empty_cells_of_the_source():
    header = '\t'.join(COLUMNS)
    line = (
        'AANAT\tMTNR1A\tSerotonin N-acetyltransferase\t'
        'Melatonin receptor type 1A\tSNAT\t""\tSerotonin acetylase\t'
        'Mel-1A-R || Mel1a receptor\tfantom5,HPMR\t12943195\t""\t""\t""'
    )
    assert list(iter_pairs(table([header, line]))) == [
        row(
            'AANAT', 'MTNR1A', 'Serotonin N-acetyltransferase',
            'Melatonin receptor type 1A', 'SNAT', '', 'Serotonin acetylase',
            'Mel-1A-R || Mel1a receptor', 'fantom5,HPMR', '12943195', '', '', '',
        )
    ]


def test_reader_rejects_a_table_with_changed_columns():
    with pytest.raises(ValueError, match='ligand.name'):
        list(iter_pairs(table(['ligand\treceptor\tsource', 'A2M\tLRP1\tHPRD'])))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_pairs(table(['\t'.join(COLUMNS)])))
