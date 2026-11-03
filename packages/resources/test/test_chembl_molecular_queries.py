"""Execute both ChEMBL query paths on assay-linked variants and typed components."""

import json
import sqlite3

import duckdb
import pytest
from pypath.inputs_v2 import chembl
from pypath.inputs_v2.parsers.chembl import PARQUET_QUERIES, SQLITE_QUERIES

SCHEMA = {
    'activities': 'activity_id INTEGER, molregno INTEGER, assay_id INTEGER, doc_id INTEGER, standard_type TEXT, standard_relation TEXT, standard_value DOUBLE, standard_units TEXT, pchembl_value DOUBLE, data_validity_comment TEXT, action_type TEXT',
    'molecule_dictionary': 'molregno INTEGER, chembl_id TEXT',
    'compound_structures': 'molregno INTEGER, canonical_smiles TEXT, standard_inchi TEXT, standard_inchi_key TEXT',
    'assays': 'assay_id INTEGER, tid INTEGER, doc_id INTEGER, chembl_id TEXT, description TEXT, assay_type TEXT, variant_id INTEGER, assay_tax_id INTEGER, confidence_score INTEGER, assay_category TEXT, assay_subcellular_fraction TEXT, assay_tissue TEXT, assay_cell_type TEXT',
    'variant_sequences': 'variant_id INTEGER, mutation TEXT, accession TEXT, version INTEGER, isoform INTEGER, sequence TEXT',
    'target_dictionary': 'tid INTEGER, target_type TEXT, pref_name TEXT, tax_id INTEGER, organism TEXT, chembl_id TEXT',
    'target_components': 'tid INTEGER, component_id INTEGER',
    'component_sequences': 'component_id INTEGER, component_type TEXT, db_source TEXT, accession TEXT, description TEXT, sequence TEXT',
    'docs': 'doc_id INTEGER, chembl_id TEXT, pubmed_id TEXT, doi TEXT',
    'assay_parameters': 'assay_id INTEGER, standard_type TEXT, standard_value DOUBLE, standard_units TEXT',
    'drug_mechanism': 'molregno INTEGER, tid INTEGER, action_type TEXT, mechanism_of_action TEXT, mechanism_comment TEXT, selectivity_comment TEXT, binding_site_comment TEXT, direct_interaction INTEGER, molecular_mechanism INTEGER, disease_efficacy INTEGER',
    'action_type': 'action_type TEXT, description TEXT, parent_type TEXT',
}


def make_db(engine):
    con = (
        sqlite3.connect(':memory:')
        if engine == 'sqlite'
        else duckdb.connect(':memory:')
    )
    for table, columns in SCHEMA.items():
        con.execute(f'CREATE TABLE {table} ({columns})')
    con.execute(
        "INSERT INTO target_dictionary VALUES (1, 'SINGLE PROTEIN', 'target', 9606, 'Human', 'CHEMBLT1')"
    )
    con.execute(
        "INSERT INTO component_sequences VALUES (1, 'PROTEIN', 'SWISS-PROT', 'P16455', 'description, comma', 'MASG')"
    )
    con.execute('INSERT INTO target_components VALUES (1, 1)')
    con.execute("INSERT INTO molecule_dictionary VALUES (1, 'CHEMBL1')")
    con.execute(
        "INSERT INTO variant_sequences VALUES (12, 'A2G/S3D', 'P16455', 2, 1, 'MGDG'), (13, 'A2V', 'P16455', 2, 1, 'MVSG')"
    )
    for assay, variant in ((1, 12), (2, None), (3, 13)):
        con.execute(
            'INSERT INTO assays (assay_id, tid, chembl_id, description, variant_id) VALUES (?, 1, ?, ?, ?)',
            (assay, 'ASSAY' + str(assay), 'source assay', variant),
        )
        con.execute(
            'INSERT INTO activities (activity_id, molregno, assay_id, pchembl_value) VALUES (?, 1, ?, 8)',
            (assay, assay),
        )
    return con


def rows(con, query):
    cursor = con.execute(query)
    return [
        dict(zip([column[0] for column in cursor.description], row))
        for row in cursor.fetchall()
    ]


@pytest.mark.parametrize('engine', ['sqlite', 'duckdb'])
def test_activity_variants_join_by_assay_without_target_fanout(engine):
    con = make_db(engine)
    try:
        queries = SQLITE_QUERIES if engine == 'sqlite' else PARQUET_QUERIES
        activities = sorted(
            rows(con, queries['activities']), key=lambda row: row['activity_id']
        )
        assert len(activities) == 3
        assert [row['variant_id'] for row in activities] == [12, None, 13]
        forms = [
            chembl.activities_schema(row).object.molecular_form
            for row in activities
        ]
        assert [feature['alternate'] for feature in forms[0]['variants']] == [
            'G',
            'D',
        ]
        assert forms[1] is None
        assert forms[2]['variants'][0]['alternate'] == 'V'
        assays = rows(con, queries['assays'])
        assert len(assays) == 3
        variant = chembl.assay_variants_schema(
            {**assays[0], '_variant_observation': True}
        )
        assert variant.molecular_form['variants']
        (target,) = rows(con, queries['targets'])
        (component,) = json.loads(target['component_records'])
        assert component['description'] == 'description, comma'
        assert component['sequence'] == 'MASG'
        assert (
            chembl.targets_schema(target)
            .membership[0]
            .member.molecular_form['sequence_identifiers']
        )
    finally:
        con.close()
