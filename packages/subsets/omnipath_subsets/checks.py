"""Checks of the subsets of a release: what must hold of any build, relative to its release.

Nothing here is a frozen count. The tables are compared with each other (a set's size with
its members, a class's size in the report with its members), with the rules of their
contract (the set type of a resource, the prefix of a MACdb set) and with the release they
were built from (every member is a chemical entity of it, the stamp is its digest).

A finding is an ``error`` (the table breaks its contract) or a ``warning`` (a gap worth
knowing, such as sets a resource leaves unnamed). ``check_release`` reads the verified
subsets of a release; ``check_metsigdb`` and ``check_exclusion`` take table files, which is
what a test of a broken table needs.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import duckdb

from omnipath_core.subset_layout import PRODUCTS, read_subset
from omnipath_subsets.exclusion_metabolites.resolve import CHEMICAL_ENTITY_TYPES
from omnipath_subsets.exclusion_metabolites.seeds import CLASSES
from omnipath_subsets.metsigdb.mapping import KEGG_OVERVIEW_MAP_NUMBERS, RESOURCES
from omnipath_subsets.release_source import open_release

SAMPLE = 5

SET_TYPES = {rule.name: rule.set_type for rule in RESOURCES}
# The resources whose sets carry a name in every release; the others are named only where
# their source names them, so an unnamed set is a gap to report and not a broken table.
NAMED_RESOURCES = ("Reactome", "MACdb", "ClassyFire")
MANDATORY = (
    "resource",
    "set_source_id",
    "metabolite_entity_key",
    "metabolite_label",
    "metabolite_entity_type",
    "set_type",
    "set_size",
    "provenance_source",
    "build_id",
)
SET_SUB_TYPES = {
    "cancer",
    "phenotype",
    "medical intervention",
    "gene abnormality",
    "genotype",
    "overview_map",
    "metabolic_map",
}


@dataclass(frozen=True)
class Finding:
    check: str
    level: str
    message: str
    count: int = 0
    sample: tuple = ()

    def as_dict(self) -> dict:
        return asdict(self)


def _literal(value) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _quoted_list(values) -> str:
    return "[" + ", ".join(_literal(v) for v in values) + "]"


class _Table:
    """Run queries over one table file and turn the rows they find into findings."""

    def __init__(self, path: Path, name: str):
        self.con = duckdb.connect()
        self.con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet({_literal(path)})")
        self.findings: list[Finding] = []

    def expect_none(self, check: str, message: str, query: str, *, level: str = "error") -> None:
        """The rows of ``query`` break the check; ``query`` selects what to show."""
        rows = self.con.execute(query).fetchall()
        if rows:
            self.findings.append(
                Finding(check, level, message, len(rows), tuple(rows[:SAMPLE]))
            )


def check_metsigdb(
    table: str | Path,
    *,
    release: duckdb.DuckDBPyConnection | None = None,
    build_id: str | None = None,
) -> list[Finding]:
    """The membership table against its contract, and against the release if given."""
    t = _Table(Path(table), "m")
    t.expect_none(
        "mandatory_fields",
        "rows lack a mandatory field",
        f"SELECT resource, set_source_id, metabolite_entity_key FROM m WHERE {' OR '.join(f'{c} IS NULL' for c in MANDATORY)}",
    )
    t.expect_none(
        "no_duplicates",
        "a metabolite is in a set more than once",
        "SELECT resource, set_source_id, metabolite_entity_key, count(*) FROM m "
        "GROUP BY 1, 2, 3 HAVING count(*) > 1",
    )
    cases = " ".join(
        f"WHEN resource = {_literal(name)} THEN {_literal(kind)}" for name, kind in SET_TYPES.items()
    )
    t.expect_none(
        "set_types",
        "a set's type is not the one its resource contributes, or its resource is unknown",
        f"SELECT DISTINCT resource, set_type FROM m "
        f"WHERE set_type IS DISTINCT FROM CASE {cases} END",
    )
    t.expect_none(
        "set_size",
        "a set's size is not its number of members",
        "SELECT resource, set_source_id, any_value(set_size), count(*) FROM m "
        "GROUP BY 1, 2 HAVING count(DISTINCT set_size) > 1 OR any_value(set_size) <> count(*)",
    )
    t.expect_none(
        "set_size_positive",
        "a set has no members",
        "SELECT DISTINCT resource, set_source_id FROM m WHERE set_size < 1",
    )
    if build_id is not None:
        t.expect_none(
            "build_stamp",
            "rows are stamped with another release",
            f"SELECT DISTINCT build_id FROM m WHERE build_id IS DISTINCT FROM {_literal(build_id)}",
        )
    t.expect_none(
        "set_labels",
        "sets of a resource that names every set have no name",
        "SELECT DISTINCT resource, set_source_id FROM m WHERE set_label IS NULL "
        f"AND resource IN {tuple(NAMED_RESOURCES)!r}",
    )
    t.expect_none(
        "set_labels",
        "sets have no name in this release",
        "SELECT resource, count(DISTINCT set_source_id) FROM m WHERE set_label IS NULL "
        f"AND resource NOT IN {tuple(NAMED_RESOURCES)!r} GROUP BY 1",
        level="warning",
    )
    t.expect_none(
        "chemont_prefix",
        "a ClassyFire set is not a ChemOnt term",
        "SELECT DISTINCT set_source_id FROM m WHERE resource = 'ClassyFire' "
        "AND set_source_id NOT LIKE 'CHEMONTID:%'",
    )
    t.expect_none(
        "macdb_prefix",
        "a MACdb set id lacks its prefix and can be read as a ChEBI or PubChem number",
        "SELECT DISTINCT set_source_id FROM m WHERE resource = 'MACdb' "
        "AND set_source_id NOT LIKE 'MACDB:%'",
    )
    t.expect_none(
        "organism",
        "a reference map has an organism, or an organism is not a taxon number",
        "SELECT DISTINCT resource, set_source_id, organism FROM m WHERE "
        "(resource = 'KEGG' AND regexp_matches(set_source_id, '^rn[0-9]{5}$') AND organism IS NOT NULL) "
        "OR organism <= 0",
    )
    numbers = _quoted_list(KEGG_OVERVIEW_MAP_NUMBERS)
    t.expect_none(
        "overview_maps",
        "a KEGG map's sub-type disagrees with its number",
        "SELECT DISTINCT set_source_id, set_sub_type FROM m WHERE resource = 'KEGG' AND "
        f"set_sub_type IS DISTINCT FROM CASE WHEN list_contains({numbers}, "
        "regexp_extract(set_source_id, '^[a-z]+([0-9]{5})$', 1)) THEN 'overview_map' "
        "ELSE 'metabolic_map' END",
    )
    # A trait type that MACdb adds is a gap to know about, not a broken table: the API's
    # `set_sub_type` filter refuses a value outside its vocabulary until it is added there.
    t.expect_none(
        "set_sub_types",
        "a sub-type is outside the vocabulary the API filters on",
        "SELECT DISTINCT set_sub_type FROM m WHERE set_sub_type IS NOT NULL "
        f"AND set_sub_type NOT IN {tuple(sorted(SET_SUB_TYPES))!r}",
        level="warning",
    )
    t.expect_none(
        "exclusion_classes",
        "an exclusion class is outside the vocabulary",
        "SELECT DISTINCT c FROM (SELECT unnest(exclusion_classes) AS c FROM m) "
        f"WHERE c NOT IN {tuple(CLASSES)!r}",
    )
    t.expect_none(
        "has_structure",
        "the structure flag disagrees with the identifiers",
        "SELECT metabolite_entity_key, has_structure FROM m WHERE has_structure IS DISTINCT FROM "
        "(inchikey IS NOT NULL OR smiles IS NOT NULL OR metabolite_structure_key IS NOT NULL) "
        "GROUP BY ALL",
    )
    if release is not None:
        _entities(t, release, "m", "metabolite_entity_key", "metabolite_entity_type")
    return t.findings


def _entities(t: _Table, release, table: str, column: str, type_column: str | None = None) -> None:
    """Every entity of the table is a chemical entity of the release."""
    keys = [row[0] for row in t.con.execute(f"SELECT DISTINCT {column} FROM {table}").fetchall()]
    if not keys:
        return
    found = dict(
        release.execute(
            "SELECT entity_key, entity_type FROM release_entity WHERE list_contains(?, entity_key)",
            [keys],
        ).fetchall()
    )
    missing = [k for k in keys if k not in found]
    other = [k for k in keys if k in found and found[k] not in CHEMICAL_ENTITY_TYPES]
    if missing:
        t.findings.append(
            Finding(
                "canonical_entity",
                "error",
                "entities the release does not have",
                len(missing),
                tuple(missing[:SAMPLE]),
            )
        )
    if other:
        t.findings.append(
            Finding(
                "canonical_entity",
                "error",
                "entities that are not chemical entities of the release",
                len(other),
                tuple(other[:SAMPLE]),
            )
        )


def check_exclusion(
    members: str | Path,
    report: str | Path,
    *,
    release: duckdb.DuckDBPyConnection | None = None,
    build_id: str | None = None,
) -> list[Finding]:
    """The exclusion tables against their contract, and against the release if given."""
    t = _Table(Path(members), "m")
    t.con.execute(f"CREATE VIEW r AS SELECT * FROM read_parquet({_literal(Path(report))})")
    t.expect_none(
        "classes",
        "a member has no class, or one outside the vocabulary",
        "SELECT DISTINCT exclusion_class FROM m "
        f"WHERE exclusion_class IS NULL OR exclusion_class NOT IN {tuple(CLASSES)!r}",
    )
    t.expect_none(
        "no_duplicates",
        "an entity is listed twice in a class by one source",
        "SELECT entity_key, exclusion_class, assignment_source, count(*) FROM m "
        "GROUP BY 1, 2, 3 HAVING count(*) > 1",
    )
    t.expect_none(
        "structure_key_members",
        "a structure-key row does not name the member it shares its key with",
        "SELECT entity_key, exclusion_class FROM m WHERE match_kind = 'structure_key' "
        "AND (via_entity_key IS NULL OR via_entity_key = entity_key)",
    )
    t.expect_none(
        "origin",
        "a row's origin disagrees with its source",
        "SELECT DISTINCT assignment_source, origin FROM m WHERE origin IS DISTINCT FROM CASE "
        "WHEN assignment_source IN ('metaproviz_seed', 'project_seed', 'override_add') "
        "THEN 'curated' ELSE 'derived' END",
    )
    t.expect_none(
        "report_classes",
        "the report's size of a class is missing or is not the class's number of members",
        "SELECT c.cls, r.resolved, s.n FROM unnest("
        f"{_quoted_list(CLASSES)}) AS c(cls) "
        "LEFT JOIN (SELECT exclusion_class, resolved FROM r WHERE kind = 'class_size') r "
        "ON r.exclusion_class = c.cls "
        "LEFT JOIN (SELECT exclusion_class, count(DISTINCT entity_key) AS n FROM m GROUP BY 1) s "
        "ON s.exclusion_class = c.cls "
        "WHERE r.resolved IS NULL OR r.resolved IS DISTINCT FROM coalesce(s.n, 0)",
    )
    if build_id is not None:
        t.expect_none(
            "build_stamp",
            "rows are stamped with another release",
            f"SELECT DISTINCT build_id FROM m WHERE build_id IS DISTINCT FROM {_literal(build_id)}",
        )
    if release is not None:
        _entities(t, release, "m", "entity_key")
    return t.findings


def check_release(
    data_root: str | Path,
    release: str,
    products: tuple[str, ...] | None = None,
) -> dict[str, list[Finding]]:
    """The checks of each present subset of ``release``, read through the layout's checks."""
    root = Path(data_root)
    results: dict[str, list[Finding]] = {}
    connection = None
    try:
        for product in products or PRODUCTS:
            subset = read_subset(root, release, product)
            if subset is None:
                continue
            if connection is None:
                connection = open_release(root, release)
            build_id = subset.manifest["release_manifest_sha256"]
            if product == "metsigdb":
                results[product] = check_metsigdb(
                    subset.directory / "metsigdb_membership.parquet",
                    release=connection,
                    build_id=build_id,
                )
            else:
                results[product] = check_exclusion(
                    subset.directory / "exclusion_metabolite.parquet",
                    subset.directory / "exclusion_metabolite_report.parquet",
                    release=connection,
                    build_id=build_id,
                )
    finally:
        if connection is not None:
            connection.close()
    return results
