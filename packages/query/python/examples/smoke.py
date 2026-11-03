"""Opt-in end-to-end check against the live API and file host."""

import argparse
import tempfile

from omnipath_client import connect, init


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--api-url", default="https://omnipath-metabo-dev.schaul.click/api")
    parser.add_argument("--release", default="latest")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="omnipath-client-smoke-") as directory:
        options = dict(data_dir=directory, api_url=args.api_url, release=args.release)
        init("signor", progress=print, **options)
        init("signor", **options)
        print("Second init: nothing to download", flush=True)
        with connect(directory) as db:
            version = db.sql("SELECT version FROM resources").fetchone()[0]
            rows = db.sql(
                """SELECT relation_key, subject_entity_key FROM signor.relation
                WHERE taxon = '9606' ORDER BY relation_key LIMIT 5"""
            ).fetchall()
            assert len(rows) == 5
            print(f"SQL: 5 rows from SIGNOR {version}", flush=True)
            joined = db.sql(
                """SELECT r.relation_key FROM relation r JOIN entity e
                ON r.subject_entity_key = e.entity_key AND r._resource = e._resource
                WHERE r.relation_key = ?""",
                [rows[0][0]],
            )
            assert joined.fetchall() == [(rows[0][0],)]
            print("Parameterized SQL join: passed", flush=True)
            found = db.search("TP53", limit=5).fetchall()
            assert found
            relations = db.relations("TP53").count("*").fetchone()[0]
            assert relations
            assert db.entity(found[0][0])["entity"]["label"]
            print(f"Prebuilt queries: TP53 has {relations} SIGNOR relations", flush=True)
        init("signor", materialize=True, **options)
        with connect(directory) as db:
            assert db.relations("TP53").count("*").fetchone()[0] == relations
            print("Native tables: same answer", flush=True)


if __name__ == "__main__":
    main()
