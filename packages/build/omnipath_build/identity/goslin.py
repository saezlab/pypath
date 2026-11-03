"""Goslin 2 normalized shorthand identifiers of the chemical hubs (hub-index ``goslin`` stage).

IDs are local level-qualified descriptors, not accession numbers issued by Goslin.
Only the most specific successfully parsed name level on each record is admitted.

In a general compound database (ChEBI, HMDB, RefMet) a name can parse as a lipid
by its shape alone: water's synonym "HOH" reads as the hydroxy fatty acid
"FA 18:5;OH", and an anion named "-oate" as a wax ester "WE 18:2/0:0". A name the
source writes in lipid notation (a chain such as ``16:0``) is the source's own
lipid statement and stands. A parse of any other name counts only where it
agrees with the structure the record states (a formula, or InChI's formula
layer): the same composition, hydrogen aside, since charge states differ in it.
The lipid databases' names are lipid names by construction and are not checked.
A wax ester without its alcohol ("WE 6:0/0:0") is Goslin's reading of an
"-oate" name in every hub, and is no claim anywhere.
"""

from __future__ import annotations

import concurrent.futures
import contextlib
import json
import re
import time
from collections import Counter
from functools import lru_cache
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from omnipath_resolver.goslin_cache import NormalizationCache

from omnipath_resolver.goslin import normalize

BATCH_NAMES = 128
SOURCES = ("swisslipids", "lipidmaps", "hmdb", "chebi", "refmet")
LIPID_DATABASES = frozenset({"swisslipids", "lipidmaps"})
_ELEMENT = re.compile(r"(\d*)([A-Z][a-z]?)(\d*)")
# A wax ester whose alcohol chain is empty: Goslin's reading of an "-oate" name.
EMPTY_ALCOHOL = r"^[a-z_]+:WE .*/0:0$"
# A chain written as carbons:double bonds: the name is lipid notation.
LIPID_NOTATION = r"\d+:\d+"
# Formula text no composition can be read from: R groups, repeat units, wildcards.
_UNREADABLE = re.compile(r"[^A-Za-z0-9.+\-]|R|X")
SCHEMA = pa.schema(
    [(n, pa.string()) for n in ("name", "goslin", "level", "status")]
    + [("specificity", pa.int32())]
)


def quote(value):
    return "'" + str(value).replace("'", "''") + "'"


def parse_batch(names):
    return [normalize(name) for name in names]


def heavy_atoms(formula):
    """Element counts without hydrogen, one Counter per ``.``-separated component.

    None where the formula cannot be read (R groups, repeat units, charges in
    brackets): a structure that cannot be compared vetoes nothing.
    """
    text = re.sub(r"[+-]\d*$", "", (formula or "").strip())
    if not text or _UNREADABLE.search(text):
        return None
    components = []
    for part in text.split("."):
        multiplier, body = re.match(r"(\d*)(.*)", part).groups()
        counts = Counter()
        for prefix, element, count in _ELEMENT.findall(body):
            if prefix:
                return None
            if element != "H":
                counts[element] += int(count or 1)
        if not counts and body.strip("H0123456789"):
            return None
        components.extend([counts] * int(multiplier or 1))
    return components or None


@lru_cache(maxsize=1 << 18)
def lipid_atoms(goslin):
    """Heavy-atom composition of a Goslin claim ``level:name``; None if it has none."""
    from pygoslin.parser.Parser import LipidParser

    try:
        formula = _parser(LipidParser).parse(goslin.split(":", 1)[1]).get_sum_formula()
    except Exception:
        return None
    atoms = heavy_atoms(formula)
    return atoms[0] if atoms and len(atoms) == 1 else None


_PARSERS = {}


def _parser(kind):
    if kind not in _PARSERS:
        _PARSERS[kind] = kind()
    return _PARSERS[kind]


def contradicts(goslin, formula):
    """The record's own structure has another composition than the lipid named."""
    record = heavy_atoms(formula)
    lipid = lipid_atoms(goslin)
    if record is None or lipid is None:
        return False
    return lipid != sum(record, Counter()) and all(lipid != part for part in record)


def _record_structures(c, hubs):
    """``(record_id, formula)`` of the general databases' records that state one."""
    pieces = [
        f"""SELECT {quote(hub + ":")} || CASE WHEN {quote(hub)}='chebi' THEN
              'CHEBI:' || regexp_replace(upper(trim(hub_id)), '^CHEBI:', '') ELSE trim(hub_id) END record_id,
            coalesce(min(source_id) FILTER (WHERE source_type='formula'),
                     min(regexp_extract(source_id, '^InChI=1S?/([^/]+)', 1))
                       FILTER (WHERE source_type='inchi')) formula
            FROM read_parquet({quote(path)})
            WHERE source_type IN ('formula', 'inchi') AND hub_id IS NOT NULL AND source_id IS NOT NULL
            GROUP BY 1"""
        for hub, path in hubs.items()
        if hub in SOURCES and hub not in LIPID_DATABASES
    ]
    if not pieces:
        return "SELECT NULL::VARCHAR record_id, NULL::VARCHAR formula WHERE false"
    return " UNION ALL ".join(pieces)


def build(hubs, output, workers=6, memory="4GB", cache_dir=None):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    c = duckdb.connect()
    c.execute(f"SET memory_limit={quote(memory)}")
    c.execute(f"SET threads={int(workers)}")
    c.execute("SET preserve_insertion_order=false")
    c.execute(f"SET temp_directory={quote(output / 'spill')}")
    pieces = []
    for hub, path in hubs.items():
        if hub not in SOURCES:
            continue
        # Preserve export shorthand columns and systematic fatty-acyl names.
        # Goslin 2 can parse names without a carbon:double-bond token; the parser
        # decides eligibility, without manual chemical-name translations.
        pieces.append(f"""SELECT {quote(hub)} hub,hub_id,
            {quote(hub + ":")} || CASE WHEN {quote(hub)}='chebi' THEN
            'CHEBI:' || regexp_replace(upper(trim(hub_id)), '^CHEBI:', '') ELSE trim(hub_id) END record_id,
            source_type name_type,source_id source_name,trim(part.name) AS name FROM read_parquet({quote(path)}),
            UNNEST(CASE WHEN {quote(hub)}='swisslipids' AND source_type IN ('synonym','lipid_shorthand') THEN string_split(source_id,' | ')
                        WHEN {quote(hub)}='lipidmaps' AND source_type='synonym' THEN string_split(source_id,'; ')
                        ELSE [source_id] END) part(name)
            WHERE source_type IN ('name','synonym','lipid_shorthand','systematic_name') AND hub_id IS NOT NULL
            AND trim(hub_id)<>'' AND source_id IS NOT NULL
            AND trim(part.name)<>''""")
    raw = (
        " UNION ALL ".join(pieces)
        if pieces
        else "SELECT NULL::VARCHAR hub,NULL::VARCHAR hub_id,NULL::VARCHAR record_id,NULL::VARCHAR name_type,NULL::VARCHAR source_name,NULL::VARCHAR AS name WHERE false"
    )
    c.execute(
        f"COPY (SELECT DISTINCT * FROM ({raw})) TO {quote(output / 'names.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
    )
    c.execute(
        f"COPY (SELECT DISTINCT name FROM read_parquet({quote(output / 'names.parquet')})) TO {quote(output / 'unique-names.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
    )
    total = pq.ParquetFile(output / "unique-names.parquet").metadata.num_rows
    cache = NormalizationCache(cache_dir or output.parent / ".goslin-cache", normalize)
    try:
        cached = 0
        with (
            pq.ParquetWriter(output / "cached.parquet", SCHEMA, compression="zstd") as hits,
            pq.ParquetWriter(
                output / "pending-names.parquet",
                pa.schema([("name", pa.string())]),
                compression="zstd",
            ) as misses,
        ):
            for batch in pq.ParquetFile(output / "unique-names.parquet").iter_batches(
                batch_size=4096
            ):
                names = batch.column("name").to_pylist()
                found = cache.lookup(names)
                seen = {row["name"] for row in found}
                if found:
                    hits.write_table(pa.Table.from_pylist(found, schema=SCHEMA))
                    cached += len(found)
                missing = [name for name in names if name not in seen]
                if missing:
                    misses.write_table(pa.table({"name": missing}))
        print(
            json.dumps(
                dict(
                    event="goslin_cache",
                    cached_names=cached,
                    new_names=total - cached,
                    cache=str(cache.path),
                )
            ),
            flush=True,
        )
        print(
            json.dumps(dict(event="goslin_start", total_names=total, workers=workers)),
            flush=True,
        )
        start = time.monotonic()
        done = 0
        tick = start
        with pq.ParquetWriter(output / "normalized.parquet", SCHEMA, compression="zstd") as writer:
            for batch in pq.ParquetFile(output / "cached.parquet").iter_batches():
                writer.write_batch(batch)
            # A single batch gains nothing from worker processes, which would each
            # re-import the parser; parse it here (same function, same results).
            single_batch = (
                pq.ParquetFile(output / "pending-names.parquet").metadata.num_rows <= BATCH_NAMES
            )
            with (
                contextlib.nullcontext()
                if single_batch
                else concurrent.futures.ProcessPoolExecutor(max_workers=workers)
            ) as pool:
                pending = set()

                def drain(block=False):
                    nonlocal pending, done, tick
                    ready, pending = concurrent.futures.wait(
                        pending,
                        timeout=1 if block else 0,
                        return_when=concurrent.futures.FIRST_COMPLETED,
                    )
                    for future in ready:
                        results = future.result()
                        cache.store(results)
                        writer.write_table(pa.Table.from_pylist(results, schema=SCHEMA))
                        done += len(results)
                    if time.monotonic() - tick >= 10:
                        print(
                            json.dumps(
                                dict(
                                    event="goslin_progress",
                                    parsed_names=done,
                                    total_names=total,
                                    cached_names=cached,
                                    names_per_second=round(
                                        done / max(time.monotonic() - start, 0.001), 1
                                    ),
                                    seconds=round(time.monotonic() - start),
                                )
                            ),
                            flush=True,
                        )
                        tick = time.monotonic()

                for b in pq.ParquetFile(output / "pending-names.parquet").iter_batches(
                    batch_size=BATCH_NAMES
                ):
                    if single_batch:
                        results = parse_batch(b.column("name").to_pylist())
                        cache.store(results)
                        writer.write_table(pa.Table.from_pylist(results, schema=SCHEMA))
                        done += len(results)
                        continue
                    while len(pending) >= workers * 2:
                        drain(True)
                    pending.add(pool.submit(parse_batch, b.column("name").to_pylist()))
                    drain()
                while pending:
                    drain(True)
        c.execute(f"""CREATE TEMP TABLE parsed AS SELECT n.*,p.goslin,p.level,p.specificity
            FROM read_parquet({quote(output / "names.parquet")}) n
            JOIN read_parquet({quote(output / "normalized.parquet")}) p USING(name)
            WHERE p.status='parsed' AND NOT regexp_matches(p.goslin, {quote(EMPTY_ALCOHOL)})""")
        c.execute(f"CREATE TEMP TABLE structures AS {_record_structures(c, hubs)}")
        # A claim is checked when no name of the record written in lipid
        # notation supports it.
        checked = c.execute(
            f"""SELECT p.record_id, p.goslin, min(s.formula) FROM parsed p
            JOIN structures s USING (record_id)
            WHERE p.hub NOT IN ({",".join(quote(h) for h in sorted(LIPID_DATABASES))})
              AND s.formula IS NOT NULL
            GROUP BY 1, 2
            HAVING NOT bool_or(regexp_matches(p.name, {quote(LIPID_NOTATION)}))"""
        ).fetchall()
        rejected = [(record, goslin) for record, goslin, formula in checked if contradicts(goslin, formula)]
        c.execute("CREATE TEMP TABLE rejected(record_id VARCHAR, goslin VARCHAR)")
        if rejected:
            c.executemany("INSERT INTO rejected VALUES (?, ?)", rejected)
        # Keeping only the most specific aliases prevents a record's generic synonym
        # from bridging otherwise different lipid descriptions or precision levels.
        c.execute(f"""COPY (SELECT * EXCLUDE (specificity) FROM parsed ANTI JOIN rejected USING (record_id, goslin)
            QUALIFY specificity=max(specificity) OVER(PARTITION BY record_id))
            TO {quote(output / "claims.parquet")} (FORMAT PARQUET, COMPRESSION ZSTD)""")
        counts = c.execute(f"""SELECT hub,count(DISTINCT record_id) records,count(DISTINCT goslin) identifiers
            FROM read_parquet({quote(output / "claims.parquet")}) GROUP BY hub""").fetchall()
        summary = dict(
            sources=counts,
            unique_names=total,
            parsed_names=done,
            cached_names=cached,
            structure_checked=len(checked),
            structure_rejected=len(rejected),
            seconds=round(time.monotonic() - start, 2),
            parse_status=c.execute(
                f"SELECT status,count(*) FROM read_parquet({quote(output / 'normalized.parquet')}) GROUP BY status"
            ).fetchall(),
        )
        (output / "summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(dict(event="goslin_complete", **summary)), flush=True)
    finally:
        cache.close()
        c.close()
