import duckdb
import pytest
from omnipath_resolver.canonical.identifiers import normalize_id, normalize_id_sql


@pytest.mark.parametrize(
    "value,expected",
    [
        ("HMDB95306", "HMDB0095306"),
        ("hmdb:00000606", "HMDB0000606"),
        ("HMDB0002111", "HMDB0002111"),
        ("12345678", "HMDB12345678"),
        ("HMDBbad", "HMDBbad"),
    ],
)
def test_hmdb_python_sql_agree(value, expected):
    assert normalize_id("hmdb", value) == expected
    with duckdb.connect() as c:
        assert (
            c.execute(
                "SELECT " + normalize_id_sql("'hmdb'", "v") + " FROM (SELECT ? AS v)", [value]
            ).fetchone()[0]
            == expected
        )
