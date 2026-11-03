"""Per-hub index (spec 3a): independent of the identity rules, rebuilt only when the hub changes.

    <output-root>/<hub>/<sha12 of the hub parquet>/
      manifest.json  hub, input {sha256, bytes}, code_sha256, counts, timings
      records.parquet  one row per record, sorted by local_id:
                       local_id, record_id, taxon (NULL unknown), anchor ('inchikey:K' | 'uniprot:X' |
                       NULL; only when anchor_count = 1), anchor_count, reviewed
      by_id/part=XX/   lookup rows (XX = md5(identifier)[:2]) sorted by (ns, identifier, local_id):
                       ns, identifier, local_id, tag, taxon, anchor, anchor_count
      by_record/part=XX/ every normalized hub row (XX = md5(record_id)[:2]), sorted by local_id:
                       local_id, source_type, value
      xrefs.parquet    (addition to the spec) the rows of by_record whose type is a hub name
                       (cross-references and gene links), sorted by local_id: local_id, ns,
                       identifier; the identity build reads only this, never by_id or by_record.

by_id tags: native (the record's own id), claim, secondary (uniprot-sec value repeated under ns
uniprot), symbol_synonym (genesymbol-syn value repeated under genesymbol), version_stripped
(refseq/genbank accession without its version). Names, synonyms, shorthand, smiles, inchi and
formula are only in by_record. Chemical goslin hubs also get `goslin` rows (ns goslin, identifier
`<level>:<name>`, the record's most specific parsed name at species level or finer).
"""

from __future__ import annotations

import json
import os
import shutil
import time
from pathlib import Path

from .common import (
    CHEMICAL,
    EMPTY_KEYS,
    GOSLIN_HUBS,
    HUBS,
    ONTOLOGY,
    INCHIKEY_RE,
    PARTS,
    Context,
    code_sha256,
    files,
    finish_part,
    log,
    norm,
    part_dirs,
    part_done,
    quote,
    scan,
    sha256_file,
    sql_list,
)

CODE_FILES = ("hubindex.py", "common.py")
NOT_IDENTIFIERS = (
    "name",
    "synonym",
    "lipid_shorthand",
    "systematic_name",
    "smiles",
    "inchi",
    "formula",
)
PRODUCT_RE = "(AP|NP|XP|YP|WP|ZP)_[0-9]+"
GENBANK_RE = "[A-Z]{3}[0-9]{5,}"
RNA_RE = "(NM|NR|XM|XR)_[0-9]+"
VERSION = r"\.[0-9]+$"


def hub_input(path: Path, root: Path):
    """sha256/bytes of the hub file, cached by size and mtime."""
    stat = path.stat()
    cache_file = root / "sha-cache.json"
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
    key = [str(path.resolve()), stat.st_size, stat.st_mtime_ns]
    if cache.get("key") != key:
        cache = dict(key=key, sha256=sha256_file(path), bytes=stat.st_size)
        root.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(cache))
    return dict(sha256=cache["sha256"], bytes=cache["bytes"])


# --------------------------------------------------------------------------- stages
def stage_rows(ctx, directory, hub, path):
    c = ctx.connect("rows")
    ident = norm(quote(hub), "hub_id")
    value = norm("source_type", "source_id")
    out = directory / "rows"
    shutil.rmtree(out, ignore_errors=True)
    ctx.copy_partitioned(
        c,
        f"""SELECT {ident} local_id,source_type,{value} AS value,taxonomy_id taxon,
          substr(md5({quote(hub + ":")} || {ident}),1,2) part FROM {scan(path)}
        WHERE hub_id IS NOT NULL AND trim(hub_id)<>'' AND source_id IS NOT NULL AND trim(source_id)<>''""",
        out,
        "part",
        "rows-",
    )
    n = c.execute(f"SELECT count(*) FROM {files(str(out / '*' / '*.parquet'))}").fetchone()[0]
    return dict(rows=n)


def stage_goslin(ctx, directory, hub, path):
    """Goslin names of the record's most specific parsed level (persistent parse cache)."""
    from .goslin import build

    build({hub: str(path)}, directory, ctx.threads, ctx.memory, str(ctx.goslin_cache))
    shutil.rmtree(directory / "spill", ignore_errors=True)
    c = ctx.connect("goslin")
    ctx.copy(
        c,
        f"""SELECT DISTINCT {norm(quote(hub), "hub_id")} local_id,goslin,level,
          substr(md5({quote(hub + ":")} || {norm(quote(hub), "hub_id")}),1,2) part
        FROM read_parquet({quote(directory / "claims.parquet")})""",
        directory / "names.parquet",
    )
    n = c.execute(
        f"SELECT count(*),count(DISTINCT local_id) FROM read_parquet({quote(directory / 'names.parquet')})"
    ).fetchone()
    levels = c.execute(
        f"SELECT level,count(*) FROM read_parquet({quote(directory / 'names.parquet')}) GROUP BY 1"
    ).fetchall()
    return dict(names=n[0], records=n[1], levels=dict(levels))


def stage_parts(ctx, directory, hub):
    """Per part: records, by_record (final), xrefs and the unsorted by_id rows."""
    work = directory.parent
    rows = work / "rows" / "rows"
    c = ctx.connect("parts")
    names = work / "goslin" / "names.parquet"
    c.execute(
        "CREATE TEMP TABLE gn AS "
        + (
            f"SELECT local_id,goslin,part FROM read_parquet({quote(names)})"
            if names.exists()
            else "SELECT NULL::VARCHAR local_id,NULL::VARCHAR goslin,NULL::VARCHAR part WHERE false"
        )
    )
    chemical = hub in CHEMICAL
    valid_key = f"regexp_matches(value,'{INCHIKEY_RE}') AND value NOT IN {sql_list(EMPTY_KEYS)}"
    if chemical:
        key = f"source_type='inchikey' AND {valid_key}"
        count_sql = f"count(DISTINCT CASE WHEN {key} THEN value END)::INTEGER"
        anchor_sql = f"min(CASE WHEN {key} THEN 'inchikey:' || value END)"
    elif hub in ("uniprot", "rhea"):  # the record's own accession / master reaction
        own = f"source_type={quote(hub)} AND value=local_id"
        count_sql = f"coalesce(bool_or({own})::INTEGER,0)"
        anchor_sql = f"CASE WHEN bool_or({own}) THEN {quote(hub + ':')} || local_id END"
    else:
        count_sql, anchor_sql = "0", "NULL::VARCHAR"
    for p in PARTS:
        if part_done(directory, p):
            continue
        sources = part_dirs(rows, [p])
        if not sources:
            continue
        c.execute(
            f"""CREATE OR REPLACE TEMP TABLE r AS SELECT local_id,source_type,value,taxon FROM {files(*sources)}
            UNION ALL SELECT local_id,'goslin',goslin,'0' FROM gn WHERE part={quote(p)}"""
        )
        c.execute(
            f"""CREATE OR REPLACE TEMP TABLE rec AS SELECT local_id,{quote(hub + ":")} || local_id record_id,
              nullif(min(taxon) FILTER (WHERE taxon NOT IN ('','0')),'') taxon,
              coalesce(bool_or(source_type='uniprot_entry' AND split_part(value,'_',1)<>local_id),false) reviewed,
              {count_sql} anchor_count,
              {anchor_sql} anchor_min
            FROM r GROUP BY local_id"""
        )
        ctx.copy(
            c,
            """SELECT local_id,record_id,taxon,CASE WHEN anchor_count=1 THEN anchor_min END anchor,anchor_count,reviewed
            FROM rec""",
            directory / "records" / f"{p}.parquet",
        )
        ctx.copy(
            c,
            "SELECT DISTINCT local_id,source_type,value FROM r ORDER BY local_id,source_type,value",
            ctx.out / "by_record" / f"part={p}" / "data.parquet",
            ", ROW_GROUP_SIZE 65536",
        )
        ctx.copy(
            c,
            f"""SELECT DISTINCT local_id,source_type ns,value identifier FROM r
            WHERE source_type IN {sql_list(HUBS)} AND NOT (source_type={quote(hub)} AND value=local_id)""",
            directory / "xrefs" / f"{p}.parquet",
        )
        ctx.copy_partitioned(
            c,
            rf"""SELECT x.ns,x.identifier,x.local_id,x.tag,rec.taxon,
              CASE WHEN rec.anchor_count=1 THEN rec.anchor_min END anchor,rec.anchor_count,
              substr(md5(x.identifier),1,2) part
            FROM (SELECT DISTINCT * FROM (
                SELECT local_id,source_type ns,value identifier,
                  CASE WHEN source_type={quote(hub)} AND value=local_id THEN 'native' ELSE 'claim' END tag
                  FROM r WHERE source_type NOT IN {sql_list(NOT_IDENTIFIERS)}
                UNION ALL SELECT local_id,'uniprot',value,'secondary' FROM r WHERE source_type='uniprot-sec'
                UNION ALL SELECT local_id,'genesymbol',value,'symbol_synonym' FROM r WHERE source_type='genesymbol-syn'
                UNION ALL SELECT local_id,source_type,regexp_replace(value,'{VERSION}',''),'version_stripped' FROM r
                  WHERE regexp_replace(value,'{VERSION}','')<>value AND (
                    (source_type='refseq_protein' AND regexp_full_match(regexp_replace(value,'{VERSION}',''),'{PRODUCT_RE}'))
                    OR (source_type='genbank' AND regexp_full_match(regexp_replace(value,'{VERSION}',''),'{GENBANK_RE}'))
                    OR (source_type='refseq' AND regexp_full_match(regexp_replace(value,'{VERSION}',''),'{RNA_RE}'))))) x
              JOIN rec USING(local_id)""",
            directory / "by_id_raw",
            "part",
            f"p{p}-",
        )
        finish_part(directory, p)
    return {}


def stage_by_id(ctx, directory):
    """Sort and de-duplicate the by_id rows of every identifier partition."""
    raw = directory.parent / "parts" / "by_id_raw"
    c = ctx.connect("by_id")
    rows = 0
    for p in PARTS:
        source = raw / f"part={p}"
        if not source.is_dir() or part_done(directory, p):
            continue
        ctx.copy(
            c,
            f"""SELECT DISTINCT ns,identifier,local_id,tag,taxon,anchor,anchor_count
            FROM {files(str(source / "*.parquet"))} ORDER BY ns,identifier,local_id,tag""",
            ctx.out / "by_id" / f"part={p}" / "data.parquet",
            ", ROW_GROUP_SIZE 65536",
        )
        finish_part(directory, p)
    rows = c.execute(
        f"SELECT count(*) FROM {files(str(ctx.out / 'by_id' / '*' / 'data.parquet'))}"
    ).fetchone()[0]
    return dict(rows=rows)


def stage_final(ctx, directory):
    parts = directory.parent / "parts"
    c = ctx.connect("final")
    n = ctx.copy(
        c,
        f"SELECT * FROM {files(str(parts / 'records' / '*.parquet'))} ORDER BY local_id",
        ctx.out / "records.parquet",
        ", ROW_GROUP_SIZE 65536",
    )
    x = ctx.copy(
        c,
        f"SELECT * FROM {files(str(parts / 'xrefs' / '*.parquet'))} ORDER BY local_id,ns,identifier",
        ctx.out / "xrefs.parquet",
        ", ROW_GROUP_SIZE 65536",
    )
    anchored = c.execute(
        f"SELECT anchor_count,count(*) FROM read_parquet({quote(ctx.out / 'records.parquet')}) GROUP BY 1 ORDER BY 1"
    ).fetchall()
    br = c.execute(
        f"SELECT count(*) FROM {files(str(ctx.out / 'by_record' / '*' / 'data.parquet'))}"
    ).fetchone()[0]
    return dict(
        records=n, xrefs=x, by_record=br, records_by_anchor_count={str(k): v for k, v in anchored}
    )


# ----------------------------------------------------------------------- entrypoint
def build_hub_index(
    hub,
    hubs_dir,
    output_root,
    memory="7GB",
    threads=6,
    goslin_cache=None,
    min_free_gib=50,
    keep_work=False,
):
    if hub not in HUBS + ONTOLOGY:
        raise ValueError(f"Unknown hub {hub}")
    path = Path(hubs_dir) / f"{hub}.parquet"
    root = Path(output_root) / hub
    started = time.monotonic()
    t = time.monotonic()
    source = hub_input(path, root)
    sha_seconds = round(time.monotonic() - t, 2)
    code = code_sha256(*CODE_FILES)
    final = root / source["sha256"][:12]
    if (final / "manifest.json").is_file():
        log("hub_index_exists", hub=hub, path=str(final))
        return json.loads((final / "manifest.json").read_text())
    building = root / (source["sha256"][:12] + ".building")
    building.mkdir(parents=True, exist_ok=True)
    ctx = Context(building, memory, threads, goslin_cache, min_free_gib)
    ctx.timings["input_sha256"] = sha_seconds
    log("hub_index_start", hub=hub, input=source, code_sha256=code)
    w = ctx.work
    ctx.stage("rows", [w / "rows" / "rows"], lambda c, d: stage_rows(c, d, hub, path))
    if hub in GOSLIN_HUBS:
        ctx.stage(
            "goslin", [w / "goslin" / "names.parquet"], lambda c, d: stage_goslin(c, d, hub, path)
        )
    ctx.stage("parts", [w / "parts" / "records"], lambda c, d: stage_parts(c, d, hub))
    ctx.stage("by_id", [building / "by_id"], stage_by_id)
    ctx.stage("final", [building / "records.parquet", building / "xrefs.parquet"], stage_final)
    sizes = {}
    for name in ("records.parquet", "xrefs.parquet", "by_id", "by_record"):
        p = building / name
        sizes[name] = (
            p.stat().st_size if p.is_file() else sum(f.stat().st_size for f in p.rglob("*.parquet"))
        )
    manifest = dict(
        hub=hub,
        input=source,
        code_sha256=code,
        counts=ctx.info,
        output_bytes=sizes,
        timings=dict(ctx.timings, total=round(time.monotonic() - started, 2)),
    )
    if not keep_work:
        shutil.rmtree(w, ignore_errors=True)
    (building / "manifest.json").write_text(json.dumps(manifest, indent=2))
    os.replace(building, final)
    log("hub_index_done", hub=hub, path=str(final), seconds=manifest["timings"]["total"])
    return manifest
