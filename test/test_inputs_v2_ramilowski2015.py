"""Source contract of the Ramilowski 2015 ligand-receptor input."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2 import ramilowski2015 as mod
from pypath.inputs_v2.parsers import ramilowski2015 as parser


# Rows of Supplementary Data 2, as the parser yields them.
LITERATURE = {
    'Pair.Name': 'A2M_LRP1',
    'Ligand.ApprovedSymbol': 'A2M',
    'Ligand.Name': 'alpha-2-macroglobulin',
    'Receptor.ApprovedSymbol': 'LRP1',
    'Receptor.Name': 'low density lipoprotein receptor-related protein 1',
    'DLRP': None,
    'HPMR': 'HPMR',
    'IUPHAR': None,
    'HPRD': 'HPRD',
    'STRING.binding': 'STRING.binding',
    'STRING.experiment': 'STRING.experiment',
    'HPMR.Ligand': 'A2M',
    'HPMR.Receptor': 'LRP1',
    'PMID.Manual': None,
    'Pair.Source': 'known',
    'Pair.Evidence': 'literature supported',
    'HPMR.PMID': '1702392',
    'IUPHAR.PMID': None,
    'HPRD.PMID': '10652313, 12194978',
}
STRAY_BRACKET = {
    'Pair.Name': 'APOE_LRP5',
    'Ligand.ApprovedSymbol': 'APOE',
    'Ligand.Name': 'apolipoprotein E',
    'Receptor.ApprovedSymbol': 'LRP5',
    'Receptor.Name': 'low density lipoprotein receptor-related protein 5',
    'HPMR': 'HPMR',
    'STRING.binding': 'STRING.binding',
    'Pair.Source': 'known',
    'Pair.Evidence': 'literature supported',
    'HPMR.PMID': '9832610)',
}
NO_PUBLICATION = {
    'Pair.Name': 'AREGB_EGFR',
    'Ligand.ApprovedSymbol': 'AREGB',
    'Ligand.Name': 'amphiregulin B',
    'Receptor.ApprovedSymbol': 'EGFR',
    'Receptor.Name': 'epidermal growth factor receptor',
    'IUPHAR': 'IUPHAR',
    'Pair.Source': 'known',
    'Pair.Evidence': 'literature supported',
}
PUTATIVE = {
    'Pair.Name': 'ANGPT1_ITGB1',
    'Ligand.ApprovedSymbol': 'ANGPT1',
    'Ligand.Name': 'angiopoietin 1',
    'Receptor.ApprovedSymbol': 'ITGB1',
    'Receptor.Name': 'integrin, beta 1 (fibronectin receptor, beta polypeptide, antigen CD29 includes MDF2, MSK12)',
    'STRING.binding': 'STRING.binding',
    'Pair.Source': 'novel',
    'Pair.Evidence': 'putative',
}
EXCLUDED = {
    'Pair.Name': 'ALOX5AP_ALOX5',
    'Ligand.ApprovedSymbol': 'ALOX5AP',
    'Ligand.Name': 'arachidonate 5-lipoxygenase-activating protein',
    'Receptor.ApprovedSymbol': 'ALOX5',
    'Receptor.Name': 'arachidonate 5-lipoxygenase',
    'IUPHAR': 'IUPHAR',
    'HPRD': 'HPRD',
    'STRING.binding': 'STRING.binding',
    'Pair.Source': 'known',
    'Pair.Evidence': 'EXCLUDED not receptor',
    'HPRD.PMID': '9445303',
}


def extract(row):
    record = mod.resource.datasets()['interactions'].mapper(row)
    assert record is not None
    out = SilverExtractor('ramilowski2015', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def publications(relation):
    return [
        a.value
        for a in relation.annotations
        if slot_name(a.term) == 'publications'
    ]


def data_sources(relation):
    return [
        a.value
        for a in relation.annotations
        if slot_name(a.term) == 'supporting_data_source'
    ]


def test_literature_supported_pair_reaches_silver_extractor():
    record, out = extract(LITERATURE)
    assert slot_name(record.predicate) == 'interacts_with'
    # The ligand is the subject, the receptor the object.
    assert [(i.type, i.value) for i in record.subject.identifiers] == [
        (Namespace.GENESYMBOL, 'A2M'),
        (Namespace.NAME, 'alpha-2-macroglobulin'),
    ]
    assert record.object.identifiers[0].value == 'LRP1'
    # The HPMR sheet's publication first, then the HPRD sheet's.
    assert publications(record) == [
        'PMID:1702392',
        'PMID:10652313',
        'PMID:12194978',
    ]
    assert [r.predicate for r in out.relations] == ['interacts_with']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}
    assert {e.namespace for e in out.entities.values()} == {
        Namespace.GENESYMBOL
    }
    assert {e.identifier for e in out.entities.values()} == {'A2M', 'LRP1'}


def test_ligand_and_receptor_columns_give_the_participant_roles():
    record, out = extract(LITERATURE)
    assert [
        a.value for a in record.subject.annotations if a.term == PARTICIPANT_ROLE
    ] == ['ligand']
    assert [
        a.value for a in record.object.annotations if a.term == PARTICIPANT_ROLE
    ] == ['receptor']
    # The roles are not on the relation, and after the extractor they are
    # annotations of the participant occurrences, not of the entities.
    assert all(a.term != PARTICIPANT_ROLE for a in record.annotations)
    (relation,) = out.relations
    assert sorted(
        (a['scope'], a['value'])
        for a in relation.annotations
        if a['term'] == PARTICIPANT_ROLE
    ) == [('object', 'receptor'), ('subject', 'ligand')]
    assert len(out.entities) == 2
    # A putative pair is not emitted, so it gets no roles either.
    assert mod.map_ramilowski_interaction(PUTATIVE) is None


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    row = dict(LITERATURE, **{'HPMR.PMID': '0123456, 123456789', 'HPRD.PMID': '953'})
    record, _ = extract(row)
    assert publications(record) == ['PMID:953']


def test_database_flags_are_supporting_data_sources():
    record, out = extract(LITERATURE)
    # Flags HPMR, HPRD, STRING.binding and STRING.experiment: the two STRING
    # columns name STRING once; the empty DLRP and IUPHAR cells name nothing.
    assert data_sources(record) == ['HPMR', 'HPRD', 'STRING']
    assert {slot_name(a.term) for a in record.annotations} == {
        'in_taxon',
        'supporting_data_source',
        'publications',
    }
    # The participants have their species and their role (see the next test).
    for entity in (record.subject, record.object):
        assert {
            slot_name(a.term)
            for a in entity.annotations
            if a.term != PARTICIPANT_ROLE
        } == {'in_taxon'}
    assert all(not e.annotations for e in out.entities.values())
    # Each database column maps to the identifier lrdb.py and talklr.py use.
    row = dict(
        LITERATURE,
        DLRP='DLRP',
        IUPHAR='IUPHAR',
        **{'STRING.binding': None},
    )
    assert data_sources(mod.map_ramilowski_interaction(row)) == [
        'DLRP',
        'HPMR',
        'IUPHAR',
        'HPRD',
        'STRING',
    ]
    assert data_sources(mod.map_ramilowski_interaction(NO_PUBLICATION)) == [
        'IUPHAR'
    ]


def test_unknown_database_flag_is_rejected_by_the_parser():
    content = io.BytesIO()
    book = openpyxl.Workbook()
    pairs = book.active
    pairs.title = 'All.Pairs'
    header = ['Pair.Name', *parser.DATABASE_COLUMNS]
    pairs.append(header)
    pairs.append(['A2M_LRP1', 'DLRP', None, None, 'HPRD', None, None])
    pairs.append(['A2M_LRP2', None, 'Hpmr', None, None, None, None])
    for name in parser.PMID_SHEETS:
        book.create_sheet(name).append(['Pair.Name', 'PMID'])
    book.save(content)
    content.seek(0)
    with pytest.raises(ValueError, match="Unknown Ramilowski 2015 database flag 'Hpmr'"):
        list(parser.iter_pairs(SimpleNamespace(result=content)))


def test_stray_bracket_in_a_pubmed_cell_is_read_and_garbage_is_rejected():
    assert publications(mod.map_ramilowski_interaction(STRAY_BRACKET)) == [
        'PMID:9832610'
    ]
    with pytest.raises(ValueError, match='Unknown Ramilowski 2015 PubMed ID'):
        mod.map_ramilowski_interaction(
            dict(LITERATURE, **{'HPMR.PMID': 'see text'})
        )


def test_publications_are_the_union_of_the_columns_each_once_manual_first():
    row = dict(
        LITERATURE,
        **{
            'PMID.Manual': '555',
            'HPMR.PMID': '9632725',
            'IUPHAR.PMID': '555, 777',
            'HPRD.PMID': '9632725, 11446740',
        },
    )
    assert publications(mod.map_ramilowski_interaction(row)) == [
        'PMID:555',
        'PMID:9632725',
        'PMID:777',
        'PMID:11446740',
    ]


def test_pair_without_a_listed_publication_gets_none():
    record, _ = extract(NO_PUBLICATION)
    assert publications(record) == []


def test_putative_and_excluded_pairs_are_not_emitted():
    assert mod.map_ramilowski_interaction(PUTATIVE) is None
    assert mod.map_ramilowski_interaction(EXCLUDED) is None
    for evidence in ('EXCLUDED', 'EXCLUDED not ligand'):
        row = dict(LITERATURE, **{'Pair.Evidence': evidence})
        assert mod.map_ramilowski_interaction(row) is None


def test_unknown_or_missing_evidence_is_rejected_not_defaulted():
    for evidence in ('predicted', None):
        row = dict(LITERATURE, **{'Pair.Evidence': evidence})
        with pytest.raises(
            ValueError, match='Unknown Ramilowski 2015 pair evidence'
        ):
            mod.map_ramilowski_interaction(row)


def test_missing_participant_emits_no_relation():
    row = dict(LITERATURE, **{'Receptor.ApprovedSymbol': None})
    assert mod.map_ramilowski_interaction(row) is None


def test_parser_joins_the_database_sheets_to_all_pairs():
    book = openpyxl.Workbook()
    pairs = book.active
    pairs.title = 'All.Pairs'
    header = [
        'Pair.Name', 'Ligand.ApprovedSymbol', 'Ligand.Name',
        'Receptor.ApprovedSymbol', 'Receptor.Name', 'DLRP', 'HPMR', 'IUPHAR',
        'HPRD', 'STRING.binding', 'STRING.experiment', 'HPMR.Ligand',
        'HPMR.Receptor', 'PMID.Manual', 'Pair.Source', 'Pair.Evidence',
    ]  # fmt: skip
    pairs.append(header)
    for row in (LITERATURE, PUTATIVE):
        pairs.append([row.get(column) or '' for column in header])
    hpmr = book.create_sheet('HPMR')
    hpmr.append(['Pair.Name', 'Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol', 'Ligand.HPMRSymbol', 'Receptor.HPMRSymbol', 'PMID'])  # fmt: skip
    hpmr.append(['A2M_LRP1', 'A2M', 'LRP1', 'A2M', 'LRP1', '1702392'])
    iuphar = book.create_sheet('IUPHAR')
    iuphar.append(['Pair.Name', 'Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol', 'Ligand.IUPHARSymbol', 'Receptor.IUPHARSymbol', 'PMID'])  # fmt: skip
    hprd = book.create_sheet('HPRD')
    hprd.append(['Pair.Name', 'Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol', 'GeneA.HPRDID', 'GeneA.RefSeqID', 'GeneB.HPRDID', 'GeneB.RefSeqID', 'Evidence', 'PMID'])  # fmt: skip
    hprd.append(['A2M_LRP1', 'A2M', 'LRP1', '00138', 'NP_002323.2', '00072', 'NP_000005.2', 'in vitro;in vivo', '10652313, 12194978'])  # fmt: skip
    content = io.BytesIO()
    book.save(content)
    content.seek(0)

    rows = list(parser.iter_pairs(SimpleNamespace(result=content)))

    assert rows[0] == LITERATURE
    assert rows[1]['Pair.Evidence'] == 'putative'
    assert [rows[1][column] for column in mod.PMID_COLUMNS] == [None] * 4
    assert [mod.map_ramilowski_interaction(row) is None for row in rows] == [
        False,
        True,
    ]
