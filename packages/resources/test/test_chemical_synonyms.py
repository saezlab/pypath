"""RaMP and ChEMBL give a compound's synonyms as synonyms: names that are not its preferred one."""

import importlib
from pathlib import Path
import sqlite3
import sys
import types

import duckdb
import pytest

try:
    import omnipath_resolver  # noqa: F401
except ImportError:  # the native extension is not built here: the Python half is enough
    stub = types.ModuleType("omnipath_resolver._omnipath_resolver")
    stub.resolve_precomputed_batch = stub.resolve_molecular_batch = lambda *a, **k: None
    sys.modules["omnipath_resolver._omnipath_resolver"] = stub
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "resolver"))

from omnipath_build.silver import SilverExtractor
from omnipath_resources.inputs_v2.parsers import chembl as chembl_queries

SEPARATOR = chr(31)


def identifiers(module, mapper, row):
    module = importlib.import_module("omnipath_resources.inputs_v2." + module)
    record = getattr(module, mapper)(row)
    assert record is not None
    extractor = SilverExtractor(module.config.id.value, "fixture")
    extractor.process_record(record, row, "fixture:1", 1)
    (entity,) = extractor.entities.values()
    return entity, {(i["ns"], i["identifier"] if "identifier" in i else i["id"]) for i in entity.identifiers}


# ----------------------------------------------------------------- RaMP


def test_a_ramp_synonym_is_a_synonym_of_its_analyte():
    entity, found = identifiers(
        "rampdb", "analytesynonym_schema", {"rampId": "RAMP_C_000000001", "Synonym": "ascorbic acid"}
    )
    assert ("synonym", "ascorbic acid") in found
    assert ("ramp", "RAMP_C_000000001") in found
    assert "chemical" in str(entity.entity_type).lower()


def test_a_ramp_synonym_of_a_gene_stays_on_the_gene():
    entity, found = identifiers(
        "rampdb", "analytesynonym_schema", {"rampId": "RAMP_G_000000001", "Synonym": "EGFR"}
    )
    assert ("synonym", "EGFR") in found
    assert "gene" in str(entity.entity_type).lower()


def test_the_synonyms_of_ramp_are_read_as_a_table_of_their_own():
    module = importlib.import_module("omnipath_resources.inputs_v2.rampdb")
    assert "analytesynonym" in module.table_names
    assert "analytesynonym" in module.schemas


# ----------------------------------------------------------------- ChEMBL


def test_the_synonyms_of_a_chembl_molecule_are_synonyms_next_to_its_preferred_name():
    _, found = identifiers(
        "chembl",
        "molecules_schema",
        {
            "chembl_id": "CHEMBL25",
            "pref_name": "ASPIRIN",
            "molecule_type": "Small molecule",
            "synonyms": SEPARATOR.join(["Acetylsalicylic acid", "Ecotrin"]),
        },
    )
    assert {("synonym", "Acetylsalicylic acid"), ("synonym", "Ecotrin"), ("name", "ASPIRIN")} <= found


def test_a_molecule_without_synonyms_still_projects():
    _, found = identifiers(
        "chembl",
        "molecules_schema",
        {"chembl_id": "CHEMBL1", "pref_name": "X", "molecule_type": "Small molecule", "synonyms": None},
    )
    assert not [i for i in found if i[0] == "synonym"]


# The molecules a build keeps are those with a measured activity or a mechanism; each query
# below is run on the same small database.
MOLECULES = """
CREATE TABLE molecule_dictionary (molregno INTEGER, chembl_id TEXT, pref_name TEXT,
  molecule_type TEXT, max_phase REAL, therapeutic_flag INTEGER, withdrawn_flag INTEGER,
  polymer_flag INTEGER, inorganic_flag INTEGER, natural_product INTEGER);
CREATE TABLE compound_structures (molregno INTEGER, canonical_smiles TEXT, standard_inchi TEXT,
  standard_inchi_key TEXT);
CREATE TABLE compound_properties (molregno INTEGER, full_mwt REAL, alogp REAL);
CREATE TABLE activities (molregno INTEGER, pchembl_value REAL);
CREATE TABLE drug_mechanism (molregno INTEGER);
CREATE TABLE molecule_synonyms (molregno INTEGER, synonyms TEXT, syn_type TEXT);
INSERT INTO molecule_dictionary VALUES (1, 'CHEMBL25', 'ASPIRIN', 'Small molecule', 4, 1, 0, 0, 0, 0);
INSERT INTO molecule_dictionary VALUES (2, 'CHEMBL2', 'BARE', 'Small molecule', 0, 0, 0, 0, 0, 0);
INSERT INTO molecule_dictionary VALUES (3, 'CHEMBL3', 'UNUSED', 'Small molecule', 0, 0, 0, 0, 0, 0);
INSERT INTO activities VALUES (1, 7.5);
INSERT INTO activities VALUES (2, 6.0);
INSERT INTO molecule_synonyms VALUES (1, 'Acetylsalicylic acid', 'INN');
INSERT INTO molecule_synonyms VALUES (1, 'Ecotrin', 'TRADE_NAME');
INSERT INTO molecule_synonyms VALUES (1, 'Ecotrin', 'OTHER');
INSERT INTO molecule_synonyms VALUES (1, '', 'OTHER');
INSERT INTO molecule_synonyms VALUES (3, 'Ignored', 'OTHER');
"""


def synonym_sets(rows):
    return {
        row["chembl_id"]: set(filter(None, (row["synonyms"] or "").split(SEPARATOR)))
        for row in rows
    }


def test_the_sqlite_query_carries_the_distinct_synonyms_of_each_molecule():
    connection = sqlite3.connect(":memory:")
    connection.row_factory = sqlite3.Row
    connection.executescript(MOLECULES)
    query = chembl_queries.SQLITE_QUERIES["molecules"]
    found = synonym_sets(connection.execute(query).fetchall())
    assert found == {"CHEMBL25": {"Acetylsalicylic acid", "Ecotrin"}, "CHEMBL2": set()}


def test_the_prepared_duckdb_query_carries_the_same_synonyms():
    connection = duckdb.connect()
    connection.execute("CREATE SCHEMA s")
    sqlite_tables = sqlite3.connect(":memory:")
    sqlite_tables.executescript(MOLECULES)
    for table in (
        "molecule_dictionary", "compound_structures", "compound_properties", "activities",
        "drug_mechanism", "molecule_synonyms",
    ):
        columns = [row[1] for row in sqlite_tables.execute(f"PRAGMA table_info({table})")]
        types_ = {"INTEGER": "BIGINT", "REAL": "DOUBLE", "TEXT": "VARCHAR"}
        definition = ", ".join(
            f"{row[1]} {types_[row[2]]}" for row in sqlite_tables.execute(f"PRAGMA table_info({table})")
        )
        connection.execute(f"CREATE TABLE s.{table} ({definition})")
        rows = sqlite_tables.execute(f"SELECT * FROM {table}").fetchall()
        if rows:
            marks = ", ".join("?" for _ in columns)
            connection.executemany(f"INSERT INTO s.{table} VALUES ({marks})", rows)
    # The prepared cache copies these tables with the build's own selections; the two that
    # only decide which molecules are kept are copied whole.
    for table in ("molecule_dictionary", "compound_structures", "compound_properties", "molecule_synonyms"):
        connection.execute(f"CREATE TABLE {table} AS {chembl_queries.DUCKDB_TABLES[table]}")
    for table in ("activities", "drug_mechanism"):
        connection.execute(f"CREATE TABLE {table} AS SELECT * FROM s.{table}")
    cursor = connection.execute(chembl_queries.PARQUET_QUERIES["molecules"])
    names = [column[0] for column in cursor.description]
    found = synonym_sets([dict(zip(names, row)) for row in cursor.fetchall()])
    assert found == {"CHEMBL25": {"Acetylsalicylic acid", "Ecotrin"}, "CHEMBL2": set()}


def test_a_cache_prepared_before_the_synonyms_is_not_reused():
    assert chembl_queries.CHEMBL_PARQUET_CACHE_VERSION >= 10
