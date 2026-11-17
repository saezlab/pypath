"""The Parquet subsets of a release are loaded into the release's PostgreSQL schema."""

import hashlib
import json

from psycopg2 import sql
import pytest

from omnipath_core.fixtures import write_dataset, metsigdb_key
from omnipath_postgres import loader
from omnipath_postgres import subset_loading
from subset_fixtures import MEMBERSHIP_ROWS, VERSION, canonical_digest, entity_uuid


def load(release, root, products=("exclusion_metabolites", "metsigdb")):
    _, schema, dsn = release
    return loader.finish_release(dsn, schema=schema, products=products, subsets_root=root)


def rows(release, query, *args):
    connection, schema, _ = release
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
        cur.execute(query, args)
        return cur.fetchall()


def test_the_membership_table_holds_the_subsets_rows_under_the_loaded_entity_ids(release, root):
    load(release, root)
    assert rows(release, "SELECT count(*) FROM metsigdb_membership") == [(MEMBERSHIP_ROWS,)]
    found = rows(
        release,
        "SELECT metabolite_entity_id::text, set_entity_id::text, has_structure, exclusion_classes, "
        "set_context, organism FROM metsigdb_membership "
        "WHERE resource = 'Reactome' AND set_source_id = 'R-HSA-70171' "
        "AND metabolite_label = 'ATP'",
    )
    assert found == [
        (
            entity_uuid(metsigdb_key("atp")),
            entity_uuid(metsigdb_key("R-HSA-70171")),
            True,
            ["cofactors", "cofactors_extended"],
            None,
            9606,
        )
    ]


def test_json_columns_arrive_as_json(release, root):
    load(release, root)
    (context,) = rows(
        release,
        "SELECT set_context FROM metsigdb_membership WHERE resource = 'ClassyFire' LIMIT 1",
    )[0]
    assert context == {"assignment": "direct", "depth": 0}
    (record,) = rows(release, "SELECT provenance_record FROM metsigdb_membership LIMIT 1")[0]
    assert record["row_id"] == "1"


def test_the_exclusion_tables_are_loaded_with_their_members_view(release, root):
    load(release, root)
    assert rows(release, "SELECT count(*) FROM exclusion_metabolite") == [(4,)]
    assert rows(release, "SELECT count(*) FROM exclusion_metabolite_report")[0][0] > 0
    assert rows(
        release,
        "SELECT exclusion_class, label FROM exclusion_metabolite_members WHERE label = 'atp' "
        "ORDER BY 1",
    ) == [("cofactors", "atp"), ("cofactors_extended", "atp")]


def test_the_filter_and_identifier_indexes_are_built(release, root):
    load(release, root)
    names = {
        name
        for (name,) in rows(
            release,
            "SELECT indexname FROM pg_indexes WHERE schemaname = current_schema()",
        )
    }
    assert {
        "metsigdb_membership_pkey",
        "metsigdb_membership_metabolite_idx",
        "metsigdb_membership_set_type_idx",
        "metsigdb_membership_structure_idx",
        "metsigdb_membership_chebi_idx",
        "metsigdb_membership_hmdb_idx",
        "metsigdb_membership_inchikey_idx",
        "metsigdb_membership_kegg_idx",
        "metsigdb_membership_pubchem_idx",
        "exclusion_metabolite_class_entity_idx",
    } <= names


def test_each_product_leaves_its_metadata_and_phase(release, root):
    load(release, root)
    found = rows(
        release,
        "SELECT product, release_id, manifest_sha256, stats FROM subset_build_metadata "
        "ORDER BY product",
    )
    assert [r[0] for r in found] == ["exclusion_metabolites", "metsigdb"]
    subset = root / "subsets" / VERSION / "metsigdb" / "subset_manifest.json"
    metsigdb = found[1]
    assert metsigdb[1] == VERSION
    assert metsigdb[3]["subset_manifest_sha256"] == hashlib.sha256(subset.read_bytes()).hexdigest()
    assert metsigdb[3]["rows"] == {"metsigdb_membership": MEMBERSHIP_ROWS}
    assert {r[0] for r in rows(release, "SELECT phase FROM parquet_phase")} >= {
        "exclusion_metabolites",
        "metsigdb",
    }


def test_loading_again_replaces_the_tables_with_the_same_rows(release, root):
    load(release, root)
    load(release, root, products=("metsigdb",))
    assert rows(release, "SELECT count(*) FROM metsigdb_membership") == [(MEMBERSHIP_ROWS,)]


def test_a_release_without_the_subset_is_refused_with_the_way_to_build_it(release, tmp_path):
    empty = write_dataset(tmp_path / "other")
    from omnipath_api.release_database import publish_newest
    import shutil

    publish_newest(empty, VERSION)
    shutil.copy(empty / "releases" / f"{VERSION}.json", empty / "releases" / "copy.json")
    # The loaded release is this one, but no subset was built for it.
    connection, schema, _ = release
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
        manifest = json.loads((empty / "releases" / f"{VERSION}.json").read_text())
        cur.execute(
            "UPDATE parquet_release SET manifest_sha256 = %s, manifest = %s",
            [canonical_digest(manifest), json.dumps(manifest)],
        )
    connection.commit()
    with pytest.raises(ValueError, match="omnipath-subsets build"):
        load(release, empty)


def test_a_subset_of_another_release_manifest_is_refused(release, root):
    # The release is republished with another manifest after its subsets were built, and the
    # schema is loaded from the new one.
    manifest_path = root / "releases" / f"{VERSION}.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["created_at"] = "another day"
    manifest_path.write_text(json.dumps(manifest))
    connection, schema, _ = release
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
        cur.execute("UPDATE parquet_release SET manifest_sha256 = %s", [canonical_digest(manifest)])
    connection.commit()
    with pytest.raises(ValueError, match="release manifest"):
        load(release, root)
    assert rows(release, "SELECT to_regclass('metsigdb_membership')") == [(None,)]


def test_a_loaded_release_that_is_not_the_published_one_is_refused(release, root):
    connection, schema, _ = release
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
        cur.execute("UPDATE parquet_release SET manifest_sha256 = %s", ["0" * 64])
    connection.commit()
    with pytest.raises(ValueError, match="loaded release"):
        load(release, root)


def test_a_subset_naming_an_entity_the_release_lacks_is_refused_and_leaves_nothing(release, root):
    connection, schema, _ = release
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
        cur.execute("DELETE FROM entity WHERE label = 'glucose'")
    connection.commit()
    with pytest.raises(ValueError, match="entities"):
        load(release, root)
    assert rows(release, "SELECT to_regclass('metsigdb_membership')") == [(None,)]
    assert rows(release, "SELECT count(*) FROM parquet_phase WHERE phase = 'metsigdb'") == [(0,)]


def test_the_data_root_is_needed_to_load_what_has_not_been_loaded(release):
    with pytest.raises(ValueError, match="data root"):
        loader.finish_release(release[2], schema=release[1], products=("metsigdb",))


def test_the_release_is_complete_once_loaded_and_built_products_are_all_done(release, root):
    load(release, root)
    # Network views and COSMOS are built in PostgreSQL and not run here: not complete.
    assert rows(release, "SELECT status FROM parquet_release")[0][0] == "metsigdb"
    assert set(subset_loading.LOADED_PRODUCTS) == {"exclusion_metabolites", "metsigdb"}
