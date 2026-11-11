"""Run a query of step files (``sql/queries/<name>.sql``) on a release database.

A query is a sequence of SQL statements. Each step's result is a temporary table, a list
or a value that later steps use: lists and values go into their SQL as literals, which
is what lets DuckDB read rows through the id indexes and skip row groups. The API and
the Python and R clients run the same files (the R runner is ``r/R/steps.R``).

File format::

    -- query: name
    -- param: name TYPE [= default]          (VARCHAR, INTEGER, BIGINT, DOUBLE, BOOLEAN,
                                              JSON, or a list such as VARCHAR[])
    -- output: table <step> | json <step>    (a table step; a value step holding JSON)
    -- example: {"name": value}
    -- step: name AS table|list|value|rows [when <SQL expression>] [cache]
    <SQL>
    -- steps: fragment                       (the steps of sql/fragments/<fragment>.sql)
    -- repeat: var = <SQL> next <SQL> until <value step> max N
    <steps>
    -- end

Placeholders: ``{{name}}`` a literal (a list is ``[...]``), ``{{name*}}`` a list's items
for ``IN (...)`` (a subquery without rows when empty), ``{{name+}}`` the same one item
at a time (the statement is the ``UNION ALL`` of one statement per item, up to
``EACH_MAX`` items, more are spread like ``*``; DuckDB reads one key of a sorted table
through its zonemaps, several in one ``IN`` by whole row groups), ``{{#name}}`` a list's length,
``{{name[offset:count]}}`` a slice spread like ``*`` (offset and count are params or
numbers), ``{{table:name}}`` a rows step as a table (``VALUES``), and
``{{include:fragment table}}`` the SQL of ``sql/fragments/<fragment>.sql`` reading
``table`` (its ``{{input}}``): pieces several queries share, such as the entity
summaries. A table step is ``__name`` in later SQL. ``error('invalid: ...')`` in SQL
raises ``ValueError``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from importlib.resources import files
import json
import re
import time
from typing import Any, Callable

import duckdb

_PLACEHOLDER = re.compile(r"\{\{\s*(#?)([A-Za-z_][A-Za-z0-9_:]*)(\[[^\]]*\])?([*+]?)\s*\}\}")
_EACH = re.compile(r"\{\{\s*([A-Za-z_][A-Za-z0-9_]*)\+\s*\}\}")
EACH_MAX = 64
_STEP = re.compile(
    r"--\s*step:\s*(\w+)\s+AS\s+(table|list|value|rows)(?:\s+when\s+(.*?))?(\s+cache)?\s*$", re.I
)
_REPEAT = re.compile(r"--\s*repeat:\s*(\w+)\s*=\s*(.+?)\s+next\s+(.+?)\s+until\s+(\w+)\s+max\s+(\d+)\s*$", re.I)
_TYPES = {"VARCHAR", "INTEGER", "BIGINT", "DOUBLE", "BOOLEAN", "JSON"}
_INCLUDE = re.compile(r"\{\{\s*include:(\w+)\s+(\w+)\s*\}\}")


@lru_cache(maxsize=None)
def fragment(name: str) -> str:
    path = files("omnipath_client").joinpath("sql", "fragments", f"{name}.sql")
    if not path.is_file():
        raise KeyError(f"Unknown fragment: {name}")
    return path.read_text()


def _uncommented(text: str) -> str:
    return "\n".join(line for line in text.splitlines() if not line.strip().startswith("--")).strip()


def expand(sql: str, fragments: Callable[[str], str] = fragment) -> str:
    """The SQL with its fragments included (each reading the table given)."""
    for _ in range(8):
        expanded = _INCLUDE.sub(
            lambda m: "(" + _uncommented(fragments(m[1])).replace("{{input}}", m[2]) + ")", sql
        )
        if expanded == sql:
            return sql
        sql = expanded
    raise RecursionError("Fragments include each other too deeply")


class QueryError(ValueError):
    """A query's own check failed: invalid input (``error('invalid: ...')``)."""


@dataclass
class Param:
    name: str
    type: str
    default: Any = None

    @property
    def is_list(self) -> bool:
        return self.type.endswith("[]")


@dataclass
class Step:
    name: str
    kind: str
    sql: str
    when: str | None = None
    cache: bool = False


@dataclass
class Repeat:
    var: str
    init: str
    next: str
    until: str
    max: int
    steps: list = field(default_factory=list)


@dataclass
class Query:
    name: str
    params: dict[str, Param]
    outputs: dict[str, str]
    steps: list
    examples: list[dict[str, Any]]
    budget_ms: int | None = None


def include_steps(text: str, fragments: Callable[[str], str] = None) -> str:
    """The text with each ``-- steps: name`` line replaced by the steps of the fragment
    ``name``: step sequences several queries share."""
    fragments = fragments or fragment
    for _ in range(8):
        lines, changed = [], False
        for line in text.splitlines():
            m = re.match(r"\s*--\s*steps:\s*(\w+)\s*$", line)
            if m:
                lines.extend(fragments(m[1]).splitlines())
                changed = True
            else:
                lines.append(line)
        text = "\n".join(lines) + "\n"
        if not changed:
            return text
    raise RecursionError("Step fragments include each other too deeply")


def parse(text: str, name: str = "query", fragments: Callable[[str], str] = None) -> Query:
    """A step file's header and steps."""
    text = include_steps(text, fragments)
    params: dict[str, Param] = {}
    outputs: dict[str, str] = {}
    examples: list[dict[str, Any]] = []
    budget = None
    body: list = []
    stack = [body]
    current: Step | None = None

    def close():
        nonlocal current
        if current is not None:
            current.sql = current.sql.strip()
            current = None

    for line in text.splitlines():
        stripped = line.strip()
        if m := _STEP.match(stripped):
            close()
            current = Step(m[1], m[2].lower(), "", m[3].strip() if m[3] else None, bool(m[4]))
            stack[-1].append(current)
            continue
        if m := _REPEAT.match(stripped):
            close()
            loop = Repeat(m[1], m[2], m[3], m[4], int(m[5]))
            stack[-1].append(loop)
            stack.append(loop.steps)
            continue
        if re.match(r"--\s*end\s*$", stripped):
            close()
            if len(stack) == 1:
                raise SyntaxError(f"{name}: '-- end' without '-- repeat'")
            stack.pop()
            continue
        if current is None:
            if m := re.match(r"--\s*query:\s*(\w+)", stripped):
                name = m[1]
            elif m := re.match(r"--\s*param:\s*(\w+)\s+([A-Z]+(?:\[\])?)(?:\s*=\s*(.+))?$", stripped):
                if m[2].rstrip("[]") not in _TYPES:
                    raise SyntaxError(f"{name}: unknown type {m[2]}")
                default = json.loads(m[3]) if m[3] is not None and m[3] != "NULL" else None
                params[m[1]] = Param(m[1], m[2], default)
            elif m := re.match(r"--\s*output:\s*(table|json)\s+(\w+)", stripped):
                outputs[m[1]] = m[2]
            elif m := re.match(r"--\s*example:\s*(\{.*\})\s*$", stripped):
                examples.append(json.loads(m[1]))
            elif m := re.match(r"--\s*budget_ms:\s*(\d+)", stripped):
                budget = int(m[1])
            continue
        if not stripped.startswith("--"):  # comments are not part of the SQL
            current.sql += line + "\n"
    close()
    if len(stack) != 1:
        raise SyntaxError(f"{name}: '-- repeat' without '-- end'")
    if not outputs:
        raise SyntaxError(f"{name}: no '-- output:'")
    return Query(name, params, outputs, body, examples, budget)


@lru_cache(maxsize=None)
def load(name: str) -> Query:
    """The query ``name`` from the package's ``sql/queries``."""
    path = files("omnipath_client").joinpath("sql", "queries", f"{name}.sql")
    if not path.is_file():
        raise KeyError(f"Unknown query: {name}")
    return parse(path.read_text(), name)


def names() -> list[str]:
    folder = files("omnipath_client").joinpath("sql", "queries")
    return sorted(p.name[:-4] for p in folder.iterdir() if p.name.endswith(".sql"))


def literal(value: Any, type_: str | None = None) -> str:
    """A SQL literal: NULL, a boolean, a number, a string, a list or (JSON) a dict.

    Values are written into the SQL rather than bound: DuckDB plans a statement without
    the values of its parameters, so it could neither use an index nor skip row groups.
    """
    if value is None:
        return f"NULL::{type_}" if type_ else "NULL"
    if type_ == "JSON":
        # JSON text as DuckDB returns it stays as it is; other values are encoded.
        text = value if isinstance(value, str) else json.dumps(value, sort_keys=True, separators=(",", ":"))
        return literal(text) + "::JSON"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(int(value))
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    if isinstance(value, (list, tuple)):
        item = type_[:-2] if type_ and type_.endswith("[]") else None
        if not value:
            return f"[]::{type_}" if type_ else "[]"
        return "[" + ", ".join(literal(v, item) for v in value) + "]"
    if isinstance(value, dict):
        return literal(json.dumps(value, sort_keys=True, separators=(",", ":"))) + "::JSON"
    raise TypeError(f"Unsupported value: {value!r}")


def _coerce(param: Param, value: Any) -> Any:
    if value is None:
        return None
    if param.is_list:
        if isinstance(value, (str, int, float)) and not isinstance(value, bool):
            value = [value]
        return [_coerce(Param(param.name, param.type[:-2]), v) for v in value]
    kind = param.type
    try:
        if kind == "VARCHAR":
            if not isinstance(value, str):
                raise TypeError
            return value
        if kind in ("INTEGER", "BIGINT"):
            if isinstance(value, bool) or not isinstance(value, (int, str)):
                raise TypeError
            return int(value)
        if kind == "DOUBLE":
            return float(value)
        if kind == "BOOLEAN":
            if not isinstance(value, bool):
                raise TypeError
            return value
        return value  # JSON
    except (TypeError, ValueError):
        raise TypeError(f"{param.name} must be {kind}, not {value!r}") from None


def bind(query: Query, params: dict[str, Any]) -> dict[str, Any]:
    """The query's parameters with defaults, checked and coerced."""
    unknown = set(params) - set(query.params)
    if unknown:
        raise TypeError(f"{query.name}: unknown parameter(s) {', '.join(sorted(unknown))}")
    return {
        name: _coerce(p, params[name]) if name in params else p.default
        for name, p in query.params.items()
    }


class Run:
    """One execution of a query on one connection."""

    def __init__(self, con, query: Query, values: dict[str, Any], cache, scope, fragments=fragment):
        self.con, self.query, self.cache, self.scope = con, query, cache, scope
        self.fragments = fragments
        self.trace: list | None = None
        self.values = values  # params, list and value steps, loop variables
        self.types = {name: p.type for name, p in query.params.items()}
        self.rows: dict[str, tuple[list[str], list[tuple]]] = {}
        self.tables: list[str] = []

    def render(self, sql: str, empty: bool = False) -> str:
        """The SQL with its placeholders filled in; ``empty`` (a skipped step) spreads
        every list as empty, so that its SQL stays short."""
        sql = expand(sql, self.fragments)
        each = None if empty else _EACH.search(sql)
        if each:
            items = list(self.values.get(each[1]) or [])
            if 0 < len(items) <= EACH_MAX:
                return " UNION ALL ".join(
                    f"SELECT * FROM ({self.render(sql[:each.start()] + literal(v) + sql[each.end():])})"
                    for v in items)

        def replace(m):
            count, name, slice_, spread = m[1], m[2], m[3], m[4]
            if name.startswith("table:"):
                return self._rows_table(name[6:])
            if name not in self.values:
                raise KeyError(f"{self.query.name}: unknown placeholder {name}")
            value = self.values[name]
            if count:
                return str(len(value or []))
            if slice_:
                start, length = (self._number(x) for x in slice_[1:-1].split(":"))
                value, spread = (value or [])[start : start + length], "*"
            if spread:
                items = [] if empty else list(value or [])
                # Empty: a subquery without rows, which DuckDB answers without reading the
                # table (it reads every row group for IN (NULL)).
                return ", ".join(literal(v) for v in items) if items else "SELECT NULL WHERE false"
            return literal(value, self.types.get(name))

        return _PLACEHOLDER.sub(replace, sql)

    def _number(self, text: str) -> int:
        text = text.strip()
        return int(self.values[text]) if text in self.values else int(text)

    def _rows_table(self, name: str) -> str:
        columns, rows = self.rows[name]
        if not rows:
            return "(SELECT " + ", ".join(f"NULL AS {c}" for c in columns) + " WHERE false)"
        values = ", ".join("(" + ", ".join(literal(v) for v in row) + ")" for row in rows)
        return f"(SELECT * FROM (VALUES {values}) AS v({', '.join(columns)}))"

    def execute(self, sql: str):
        try:
            return self.con.execute(sql)
        except duckdb.Error as exc:
            message = str(exc)
            if "invalid: " in message:
                raise QueryError(message.split("invalid: ", 1)[1].strip()) from None
            raise

    def scalar(self, expression: str) -> Any:
        row = self.execute(f"SELECT ({self.render(expression)})").fetchone()
        return row[0] if row else None

    def steps(self, steps) -> None:
        for item in steps:
            if isinstance(item, Repeat):
                self.values[item.var] = self.scalar(item.init)
                for round_ in range(item.max):
                    self.steps(item.steps)
                    if self.values.get(item.until) or round_ == item.max - 1:
                        break
                    self.values[item.var] = self.scalar(item.next)
            else:
                self.step(item)

    def step(self, step: Step) -> None:
        started = time.perf_counter()
        try:
            self._step(step)
        finally:
            if self.trace is not None:
                self.trace.append((step.name, time.perf_counter() - started))

    def _step(self, step: Step) -> None:
        skip = step.when is not None and not self.scalar(step.when)
        sql = self.render(step.sql, empty=skip)
        if skip:
            sql = f"SELECT * FROM ({sql}) WHERE false"
        if step.kind == "table":
            target = f"__{step.name}"
            if step.name not in self.tables:
                self.tables.append(step.name)
            if step.cache and self.cache is not None and not skip:
                table = self.cache((self.scope, self.query.name, step.name, sql),
                                   lambda: self.execute(sql).to_arrow_table())
                self.con.register("__cached", table)
                try:
                    self.execute(f"CREATE OR REPLACE TEMP TABLE {target} AS SELECT * FROM __cached")
                finally:
                    self.con.unregister("__cached")
            else:
                self.execute(f"CREATE OR REPLACE TEMP TABLE {target} AS {sql}")
            return

        def compute():
            cursor = self.execute(sql)
            columns = [c[0] for c in cursor.description]
            return columns, cursor.fetchall(), str(cursor.description[0][1])

        if step.cache and self.cache is not None and not skip:
            columns, rows, type_ = self.cache((self.scope, self.query.name, step.name, sql), compute)
        else:
            columns, rows, type_ = compute()
        # The column's type, so that an empty list or NULL is a typed literal later.
        if "(" not in type_ and step.kind in ("list", "value"):
            self.types[step.name] = type_ + "[]" if step.kind == "list" else type_
        if step.kind == "list":
            self.values[step.name] = [row[0] for row in rows]
        elif step.kind == "value":
            self.values[step.name] = rows[0][0] if rows else None
        else:
            self.rows[step.name] = (columns, rows)
            self.values[step.name] = rows

    def drop(self) -> None:
        for name in self.tables:
            try:
                self.con.execute(f"DROP TABLE IF EXISTS __{name}")
            except duckdb.Error:
                pass


@lru_cache(maxsize=None)
def _shared_macros() -> tuple[str, str]:
    """(marker, SQL) of the macros the queries share (sql/shape.sql) as temporary
    macros: each connection gets the code's own, whatever its database stores."""
    import hashlib

    from .contracts import CONTRACT

    text = files("omnipath_client").joinpath("sql", "shape.sql").read_text()
    tables = json.dumps(CONTRACT.get("query_tables", {}), sort_keys=True)
    marker = "op_layer_" + hashlib.sha256((text + tables).encode()).hexdigest()[:12]
    return marker, re.sub(r"CREATE OR REPLACE MACRO", "CREATE OR REPLACE TEMP MACRO", text)


def prepare(con: duckdb.DuckDBPyConnection) -> None:
    """Define the shared macros on the connection, once per version of them."""
    marker, sql = _shared_macros()
    try:
        con.execute(f"SELECT {marker}()")  # defined: a trivial statement
    except duckdb.CatalogException:
        from .contracts import CONTRACT

        for table, spec in CONTRACT.get("query_tables", {}).items():
            con.execute(f"CREATE OR REPLACE TEMP TABLE {table} ({', '.join(spec['columns'])})")
            if spec["rows"]:
                values = ", ".join("(" + ", ".join(literal(v) for v in row) + ")" for row in spec["rows"])
                con.execute(f"INSERT INTO {table} VALUES {values}")
        con.execute(sql)
        con.execute(f"CREATE OR REPLACE TEMP MACRO {marker}() AS 1")


def apply_settings(con: duckdb.DuckDBPyConnection) -> None:
    """The contract's settings, which every connection that runs queries uses."""
    from .contracts import CONTRACT

    for name, value in CONTRACT["settings"].items():
        con.execute(f"SET {name} = {literal(value)}")


def run(
    con: duckdb.DuckDBPyConnection,
    name: str | Query,
    params: dict[str, Any] | None = None,
    *,
    output: str = "table",
    cache: Callable[[tuple, Callable[[], Any]], Any] | None = None,
    scope: Any = None,
    fragments: Callable[[str], str] = fragment,
    trace: list | None = None,
    macros: bool = True,
):
    """Run a query: an Arrow table for ``output="table"``, a JSON string (or None) for
    ``output="json"``.

    ``cache(key, compute)`` keeps the results of steps marked ``cache``; ``scope`` (the
    database's identity) is part of their keys. ``trace``, a list, gets (step, seconds)
    for each step run. ``macros=False`` leaves out the shared macros (a database without
    the release tables).
    """
    query = load(name) if isinstance(name, str) else name
    if output not in query.outputs:
        raise ValueError(f"{query.name} has no {output} output")
    if macros:
        prepare(con)
    run_ = Run(con, query, bind(query, params or {}), cache, scope, fragments)
    run_.trace = trace
    try:
        run_.steps(query.steps)
        target = query.outputs[output]
        if output == "json":
            return run_.values.get(target)
        return run_.execute(f"SELECT * FROM __{target}").to_arrow_table()
    finally:
        run_.drop()
