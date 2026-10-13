"""Source contract of the CytoSig signature matrix input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor

from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import cytosig
from pypath.inputs_v2.cytosig import CYTOKINE, resource
from pypath.inputs_v2.parsers.cytosig import iter_cytokines, iter_signature

# The header and two rows of signature.centroid. The header has a name for each
# of the 43 cytokines and none for the column of gene symbols.
HEADER = 'Activin A\tBDNF\tBMP2\tBMP4\tBMP6\tCD40L\tCXCL12\tEGF\tFGF2\tGCSF\tGDF11\tGMCSF\tHGF\tIFN1\tIFNG\tIFNL\tIL10\tIL12\tIL13\tIL15\tIL17A\tIL1A\tIL1B\tIL2\tIL21\tIL22\tIL27\tIL3\tIL36\tIL4\tIL6\tLIF\tLTA\tMCSF\tNO\tOSM\tTGFB1\tTGFB3\tTNFA\tTRAIL\tTWEAK\tVEGFA\tWNT3A'
A4GALT = 'A4GALT\t-0.13549461211883163\t0.0\t-0.018501576481109616\t-0.03535826335658038\t0.01875768900000008\t0.06445235315889075\t0.004607347687224817\t0.0056964250938285055\t-0.04146817049103613\t0.07630236329272932\t0.1075145164614808\t-0.008029096006877179\t-0.03392385304266998\t-0.018075781499999888\t-0.03236106310465814\t0.179488786880039\t-0.035369340741539634\t0.00553631686211009\t-0.03032148382604038\t0.0\t-0.014042077125878015\t0.0011099999999999999\t0.2520213521212398\t0.03294329122038975\t-0.025046585418646183\t-0.03246967684222213\t-0.04290723580037126\t0.08025736797699956\t0.03917059728067773\t0.08590155277261212\t0.008466926651672145\t0.172277014442499\t0.5679445170294432\t0.1545099490556249\t0.05509863941238447\t0.1377930801900562\t-0.01238777154971782\t0.30077661400000105\t0.12234390776620752\t-0.03813458291730365\t0.019689339824527245\t-0.10236111650157209\t0.013742799664961058'
CSF1 = 'CSF1\t-0.0720037426000808\t0.18514274673056508\t0.05822585112355583\t0.05409874927939701\t-0.0011397589999999624\t0.036936780000001335\t0.035821347999999816\t-0.005171796352185342\t0.11802007643082402\t0.06084478421195616\t-0.13642768832241006\t0.14225601810300836\t-0.15046749453382444\t0.065393458523755\t0.2599335942500005\t0.5003696266154027\t-0.05799730433903649\t0.034045196706600216\t0.05486330666221484\t0.4325000280958258\t0.003913298999999704\t0.05242999999999999\t0.6973188760000006\t0.2504349230629419\t0.020300215638479554\t0.04694139431424027\t0.11965287098388755\t0.3697618450201876\t0.5411462684286665\t0.014246470916852612\t-0.015507190505675172\t0.22106999999999946\t-0.08381328865263654\t-0.10217488591081914\t0.061896658641874325\t0.29332729260045953\t-0.12582559467952414\t-0.009536632499999698\t0.23115505343366644\t0.2047169767785597\t0.10472908439724993\t-0.03233137562280364\t-0.04448457966146624'

def matrix(*lines, newline='\n'):
    return SimpleNamespace(result=io.StringIO(newline.join(lines) + newline))


def cells(*lines, **kwargs):
    return list(iter_signature(matrix(*lines, **kwargs)))


def extract(row, dataset='signature'):
    record = resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('cytosig', dataset)
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


def response(relation):
    (quantity,) = (a['quantity'] for a in relation.annotations if a['term'] == 'cytosig:response')
    return quantity['has_numeric_value']


def test_parser_gives_a_row_for_each_cell_with_the_gene_column_unnamed_in_the_header():
    rows = cells(HEADER, A4GALT)
    assert len(rows) == 43
    assert rows[0] == {'gene': 'A4GALT', 'cytokine': 'Activin A', 'value': '-0.13549461211883163'}
    assert rows[1] == {'gene': 'A4GALT', 'cytokine': 'BDNF', 'value': '0.0'}
    assert rows[-1]['cytokine'] == 'WNT3A'
    assert {r['gene'] for r in rows} == {'A4GALT'}
    assert len(cells(HEADER, A4GALT, CSF1)) == 86


def test_parser_reads_a_file_with_crlf_line_ends():
    assert len(cells(HEADER, A4GALT, newline='\r\n')) == 43


def test_parser_rejects_malformed_matrices():
    with pytest.raises(ValueError, match='no rows'):
        cells(HEADER)
    with pytest.raises(ValueError, match='fields'):
        cells(HEADER, A4GALT.rsplit('\t', 1)[0])
    with pytest.raises(ValueError, match='fields'):
        cells(HEADER, A4GALT + '\t0.1')
    with pytest.raises(ValueError, match='empty cell'):
        cells(HEADER, A4GALT.replace('\t0.0\t', '\t\t', 1))
    with pytest.raises(ValueError, match='repeats'):
        cells(HEADER.replace('BDNF', 'Activin A'), A4GALT)


def test_gene_is_associated_with_the_cytokine_with_the_cell_value():
    row = cells(HEADER, A4GALT)[0]
    relation, gene, cytokine = relation_of(row)
    assert relation.predicate == 'associated_with'
    assert (gene.entity_type, gene.namespace, gene.identifier, gene.taxon) == (
        'gene', 'genesymbol', 'A4GALT', '9606',
    )
    assert (cytokine.entity_type, cytokine.namespace, cytokine.identifier) == (
        'ontology_class', CYTOKINE, 'Activin A',
    )
    assert response(relation) == -0.13549461211883163


def test_cytokine_label_is_not_turned_into_a_protein():
    # IL12 is the label of a cytokine of two chains; NO is nitric oxide. The
    # statement names the label, not a gene.
    rows = {r['cytokine']: r for r in cells(HEADER, CSF1)}
    for label in ('IL12', 'NO', 'IFN1', 'MCSF', 'Activin A'):
        _, _, cytokine = relation_of(rows[label])
        assert (cytokine.entity_type, cytokine.namespace, cytokine.identifier) == (
            'ontology_class', CYTOKINE, label,
        )
    # CSF1 is a gene row of the matrix and the MCSF column is a different thing.
    relation, gene, cytokine = relation_of(rows['MCSF'])
    assert gene.identifier == 'CSF1' and cytokine.identifier == 'MCSF'


def test_zero_cells_are_statements_too():
    rows = cells(HEADER, A4GALT)
    zeros = [r for r in rows if float(r['value']) == 0.0]
    assert zeros
    relation, _, _ = relation_of(zeros[0])
    assert response(relation) == 0.0


def test_cell_that_is_not_a_number_raises():
    row = cells(HEADER, A4GALT)[0]
    for bad in ('nan', 'inf', 'x', ''):
        with pytest.raises(ValueError, match='not a finite number'):
            cytosig.map_signature(dict(row, value=bad))


def test_terms_dataset_has_the_columns_as_named_entities():
    rows = list(iter_cytokines(matrix(HEADER, A4GALT)))
    assert len(rows) == 43
    assert rows[0] == {'cytokine': 'Activin A'}
    out = extract(rows[17], 'terms')
    (term,) = (e for e in out.entities.values() if e.label)
    assert (term.entity_type, term.namespace, term.identifier, term.label) == (
        'ontology_class', CYTOKINE, 'IL12', 'IL12',
    )
    _, _, cytokine = relation_of(cells(HEADER, A4GALT)[17])
    assert (cytokine.entity_type, cytokine.namespace, cytokine.identifier) == (
        term.entity_type, term.namespace, term.identifier,
    )
    with pytest.raises(ValueError, match='header'):
        list(iter_cytokines(matrix('')))


def test_nothing_else_is_asserted():
    relation, _, _ = relation_of(cells(HEADER, A4GALT)[0])
    assert {a['term'] for a in relation.annotations} == {'cytosig:response', 'in_taxon'}


def test_resource_declares_its_datasets_and_license():
    assert cytosig.config.id == ResourceCv.CYTOSIG
    assert cytosig.config.license == LicenseCV.GPL_3_0
    assert set(resource.datasets()) == {'terms', 'signature'}
