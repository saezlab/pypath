"""MetalinksDB from the Parquet release with DuckDB, reproducing the PostgreSQL steps.

Reproduces, for the resources MetalinksDB reads:

1. The annotation linking of omnipath_postgres.projection: an evidence row's own
   annotations, plus each statement annotation whose source and dataset match the
   evidence row (or every annotation of a synthetic evidence row).
2. The interaction record of relational/db/derived_tables.py
   (_stage_interaction_class_evidence, _stage_interaction_record,
   _populate_interaction_fact_resource): class precedence per canonical relation
   over every resource of the release, direction, sign, references, affinity,
   pChEMBL, score, curation flags, organisms.
3. The preset recipe (omnipath_subsets.network_views.METALINKSDB) and the fold of
   main's reference fixture (tests/fixtures/collapse.py).

Then compares the record and the collapsed preset with the rows exported from
the PostgreSQL database loaded from the same pins.

Usage: metalinks_faithful.py <pins.json> <postgres_record.csv.gz> <out_dir>
"""
import json
import sys
import time
from pathlib import Path

import duckdb
import yaml
from omnipath_postgres.relational.classify import interaction_class as ic_module
from omnipath_postgres.relational.db import derived_tables as dt
from omnipath_subsets.network_views import METALINKSDB

ROOT = "/root/data/release-20261007-native/resources"
PINS = json.load(open(sys.argv[1]))["resources"]
PG_CSV = sys.argv[2]
OUT = Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)

# The recipe, read from the preset itself.
components = METALINKSDB.composition["components"]
CHEMBL_FILTERS = components[0]["parameters"]["filters"]
CURATED_FILTERS = components[1]["parameters"]["filters"]
REACTION_FILTERS = components[2]["parameters"]["filters"]
EXCLUDED = next(s["resources"] for s in METALINKSDB.composition["steps"] if s["operation"] == "exclude")
SCOPE = sorted(
    {*CHEMBL_FILTERS["resources"], *CURATED_FILTERS["resources"], *REACTION_FILTERS["resources"]}
    - set(EXCLUDED)
) 
SCOPE = [s for s in SCOPE if s in PINS]  # the pins may predate a resource the preset names

# The constants of the PostgreSQL derive, imported rather than restated.
PREDICATE_CLASS = yaml.safe_load(open(ic_module._RULES_PATH))["predicates"]
ANNOTATION_CLASS_TERMS = dt._ANNOTATION_CLASS_TERMS
DIRECTED_PREDICATES = dt._DIRECTED_PREDICATES
DIRECTED_CLASSES = dt._DIRECTED_CLASSES
AFFINITY_TERMS = dt._AFFINITY_TERMS
NUMERIC = dt._NUMERIC_VALUE
CURATION_TERM = dt._CURATION_TERM
PRESENCE = dt._CURATION_PRESENCE_TERMS
CLAIM_DATASET = "omnipath:published_statement"


def lit(value):
    return "'" + str(value).replace("'", "''") + "'"


def lits(values):
    return "(" + ",".join(lit(v) for v in values) + ")"


def pq(source, table):
    return f"read_parquet('{ROOT}/{source}/{PINS[source]}/{table}.parquet')"


con = duckdb.connect(str(OUT / "work.duckdb"))
con.sql("SET threads = 3; SET memory_limit = '8GB'; SET preserve_insertion_order = false")
timings = {}


def step(name, sql):
    start = time.monotonic()
    con.sql(sql)
    timings[name] = round(time.monotonic() - start, 1)


con.sql("""CREATE OR REPLACE MACRO ap_uuid(value) AS (
    (substr(md5(value),1,8)||'-'||substr(md5(value),9,4)||'-'||substr(md5(value),13,4)||'-'||
     substr(md5(value),17,4)||'-'||substr(md5(value),21,12))::UUID)""")
con.sql("CREATE OR REPLACE MACRO entity_uuid(k) AS ap_uuid(to_json(list_value('published-entity', k::VARCHAR)))")
con.sql("CREATE OR REPLACE MACRO norm_scope(s) AS CASE WHEN s = 'evidence' THEN 'relation' ELSE coalesce(s, 'relation') END")
con.sql(
    "CREATE OR REPLACE MACRO ann_value(v, q) AS "
    "CASE WHEN q.has_numeric_value IS NOT NULL THEN q.has_numeric_value::VARCHAR ELSE v END"
)


def linked_annotations(resources, terms, name):
    """Annotations linked to evidence rows, as the projection links them.

    Columns: resource, relation_id, ordinal, source, dataset, term, value, unit,
    q (quantity struct or NULL), scope.
    """
    parts = []
    for s in resources:
        term_filter = f"term IN {lits(terms)}" if terms else "true"
        parts.append(f"""
        WITH rel AS (SELECT relation_id FROM {pq(s, 'relation')} WHERE statement_kind = 'relation'),
        ev AS (
          SELECT e.relation_id, e.ordinal, coalesce(e.source, '{s}') AS source,
                 coalesce(e.dataset, '{CLAIM_DATASET}') AS dataset, false AS synthetic, e.annotations
          FROM {pq(s, 'relation_evidence')} e
          UNION ALL
          SELECT r.relation_id, -1, '{s}', '{CLAIM_DATASET}', true, NULL
          FROM rel r ANTI JOIN {pq(s, 'relation_evidence')} e USING (relation_id)
        ),
        own AS (
          SELECT relation_id, ordinal, source, dataset, a.term AS term,
                 ann_value(a.value, a.quantity) AS value, a.quantity.has_unit AS unit,
                 a.quantity AS q, norm_scope(a.scope) AS scope
          FROM (SELECT relation_id, ordinal, source, dataset, unnest(annotations) AS a FROM ev WHERE annotations IS NOT NULL)
          WHERE a.term IS NOT NULL AND a.{term_filter}
        ),
        stmt AS (
          SELECT relation_id, term, value AS raw_value,
                 CASE WHEN quantity_has_numeric_value IS NULL AND quantity_has_unit IS NULL AND quantity_has_unit_prefix IS NULL
                       AND quantity_has_binary_relation IS NULL AND quantity_source_field IS NULL AND quantity_comparator IS NULL
                      THEN NULL ELSE struct_pack(has_numeric_value := quantity_has_numeric_value, has_unit := quantity_has_unit,
                           has_unit_prefix := quantity_has_unit_prefix, has_binary_relation := quantity_has_binary_relation,
                           source_field := quantity_source_field, comparator := quantity_comparator) END AS q,
                 source AS a_source, dataset AS a_dataset, norm_scope(scope) AS scope
          FROM {pq(s, 'relation_annotation')} WHERE term IS NOT NULL AND {term_filter}
        )
        SELECT '{s}' AS resource, relation_id, ordinal, source, dataset, term, value, unit, q, scope FROM own
        UNION ALL
        SELECT '{s}', e.relation_id, e.ordinal, e.source, e.dataset, a.term, ann_value(a.raw_value, a.q),
               a.q.has_unit, a.q, a.scope
        FROM stmt a JOIN ev e ON e.relation_id = a.relation_id
         AND (e.synthetic OR ((a.a_source IS NULL OR a.a_source = e.source)
                              AND (a.a_dataset IS NULL OR a.a_dataset = e.dataset)))""")
    step(name, f"CREATE OR REPLACE TABLE {name} AS " + " UNION ALL ".join(f"SELECT * FROM ({p})" for p in parts))


def statements(resources, name):
    parts = [
        f"SELECT '{s}' AS resource, relation_id, subject_entity_key, predicate, object_entity_key, "
        f"subject_type, object_type FROM {pq(s, 'relation')} WHERE statement_kind = 'relation'"
        for s in resources
    ]
    step(name, f"CREATE OR REPLACE TABLE {name} AS " + " UNION ALL ".join(parts))


ALL = sorted(PINS)

# 1. The class evidence of every resource of the release: participant roles,
#    transport, and the annotation terms that set a class. The class belongs to
#    the canonical relation, so a resource outside the scope can set it.
class_terms = ["connectomedb:participant_role", "object_aspect_qualifier"] + sorted(
    {t for t, *_ in ANNOTATION_CLASS_TERMS}
)
class_resources = []
for s in ALL:
    try:
        n = con.sql(
            f"SELECT count(*) FROM (SELECT term FROM {pq(s, 'relation_annotation')} WHERE term IN {lits(class_terms)} "
            f"UNION ALL SELECT a.term FROM (SELECT unnest(annotations) a FROM {pq(s, 'relation_evidence')}) "
            f"WHERE a.term IN {lits(class_terms)})"
        ).fetchone()[0]
    except duckdb.IOException:
        continue
    if n:
        class_resources.append(s)
linked_annotations(class_resources, class_terms, "class_ann")
statements(sorted(set(class_resources) | set(SCOPE)), "class_stmt")
annotation_class = ",".join(
    f"({lit(t)},{lit(v)},{lit(c)},{p})" for t, v, c, p in ANNOTATION_CLASS_TERMS
)
step("class_per_triple", f"""
CREATE OR REPLACE TABLE triple_class AS
WITH per_evidence AS (
  SELECT st.subject_entity_key sk, st.predicate pred, st.object_entity_key ok,
         a.resource, a.relation_id, a.ordinal,
         bool_or(a.scope = 'subject' AND a.term = 'connectomedb:participant_role' AND a.value = 'ligand') sl,
         bool_or(a.scope = 'subject' AND a.term = 'connectomedb:participant_role' AND a.value = 'receptor') sr,
         bool_or(a.scope = 'object' AND a.term = 'connectomedb:participant_role' AND a.value = 'ligand') ol,
         bool_or(a.scope = 'object' AND a.term = 'connectomedb:participant_role' AND a.value = 'receptor') orr,
         bool_or(a.scope = 'relation' AND a.term = 'object_aspect_qualifier' AND a.value = 'transport') tr
  FROM class_ann a JOIN class_stmt st USING (resource, relation_id)
  WHERE a.term IN ('connectomedb:participant_role', 'object_aspect_qualifier')
  GROUP BY ALL
),
roles AS (
  SELECT sk, pred, ok, bool_or((sl AND orr) OR (sr AND ol)) lr, bool_or(tr) tr FROM per_evidence GROUP BY ALL
),
annotated AS (
  SELECT st.subject_entity_key sk, st.predicate pred, st.object_entity_key ok,
         arg_min(m.cls, m.prec) cls
  FROM class_ann a JOIN class_stmt st USING (resource, relation_id)
  JOIN (VALUES {annotation_class}) m(term, val, cls, prec) ON m.term = a.term AND m.val = a.value
  GROUP BY ALL
)
SELECT coalesce(r.sk, n.sk) sk, coalesce(r.pred, n.pred) pred, coalesce(r.ok, n.ok) ok,
       r.lr, r.tr, n.cls annotation_class
FROM roles r FULL JOIN annotated n USING (sk, pred, ok)
""")

# 2. The scope's statements and the annotations the record reads.
record_terms = [
    "object_direction_qualifier", "publications", "has_quantitative_value", "has_confidence_score",
    CURATION_TERM, *PRESENCE, *AFFINITY_TERMS,
]
statements(SCOPE, "stmt")
linked_annotations(SCOPE, record_terms, "rec_ann")
step("evidence_rows", f"""
CREATE OR REPLACE TABLE ev AS
SELECT s.resource, e.relation_id, e.ordinal, coalesce(e.source, s.resource) AS source,
       coalesce(e.dataset, '{CLAIM_DATASET}') AS dataset
FROM ({' UNION ALL '.join(f"SELECT '{s}' resource, relation_id, ordinal, source, dataset FROM {pq(s, 'relation_evidence')}" for s in SCOPE)}) e
JOIN (SELECT DISTINCT resource FROM stmt) s USING (resource)
UNION ALL
SELECT st.resource, st.relation_id, -1, st.resource, '{CLAIM_DATASET}'
FROM stmt st ANTI JOIN ({' UNION ALL '.join(f"SELECT '{s}' resource, relation_id FROM {pq(s, 'relation_evidence')}" for s in SCOPE)}) e
  USING (resource, relation_id)
""")
pred_map = ",".join(f"({lit(p)},{lit(c)})" for p, c in PREDICATE_CLASS.items())
step("record_evidence", f"""
CREATE OR REPLACE TABLE rec_ev AS
WITH sign AS (
  SELECT resource, relation_id, ordinal,
         nullif(bool_or(value IN ('increased','upregulated')), false) stim,
         nullif(bool_or(value IN ('decreased','downregulated')), false) inh
  FROM rec_ann WHERE term = 'object_direction_qualifier' AND scope = 'relation' GROUP BY ALL
),
classed AS (
  SELECT st.*, ev.ordinal, ev.source, ev.dataset,
    coalesce(CASE WHEN tc.lr THEN 'ligand_receptor' WHEN tc.tr THEN 'transport' END,
             tc.annotation_class,
             nullif(coalesce(pm.cls, 'other'), 'other'),
             'other') AS interaction_class
  FROM stmt st
  JOIN ev USING (resource, relation_id)
  LEFT JOIN triple_class tc ON tc.sk = st.subject_entity_key AND tc.pred = st.predicate AND tc.ok = st.object_entity_key
  LEFT JOIN (VALUES {pred_map}) pm(pred, cls) ON pm.pred = st.predicate
)
SELECT c.resource, c.relation_id, c.ordinal, c.source, c.dataset, c.predicate,
  CASE WHEN c.predicate = 'enabled_by' THEN c.object_entity_key ELSE c.subject_entity_key END sk,
  CASE WHEN c.predicate = 'enabled_by' THEN c.subject_entity_key ELSE c.object_entity_key END ok,
  CASE WHEN c.predicate = 'enabled_by' THEN c.object_type ELSE c.subject_type END subject_type,
  CASE WHEN c.predicate = 'enabled_by' THEN c.subject_type ELSE c.object_type END object_type,
  c.interaction_class,
  CASE WHEN c.predicate IN {lits(DIRECTED_PREDICATES)} THEN true
       WHEN c.interaction_class IN {lits(DIRECTED_CLASSES)} THEN true END AS is_directed,
  s.stim, s.inh
FROM classed c LEFT JOIN sign s USING (resource, relation_id, ordinal)
""")
affinity_case = f"""CASE WHEN a.term IN {lits(AFFINITY_TERMS)} AND regexp_full_match(a.value, '{NUMERIC[1:-1]}')
   AND coalesce(a.q.comparator,'=')='=' AND coalesce(a.q.has_binary_relation,'equal_to')='equal_to'
  THEN a.value::DOUBLE * CASE
    WHEN coalesce(a.q.has_unit,a.unit) IN ('nM','nanomolar') THEN 1
    WHEN coalesce(a.q.has_unit,a.unit) IN ('uM','µM','micromolar') THEN 1000
    WHEN coalesce(a.q.has_unit,a.unit) IN ('mM','millimolar') THEN 1000000
    WHEN coalesce(a.q.has_unit,a.unit) IN ('M','molar','UO:0000061') THEN
      CASE coalesce(a.q.has_unit_prefix,'') WHEN 'nano' THEN 1 WHEN 'micro' THEN 1000
        WHEN 'milli' THEN 1000000 WHEN '' THEN 1000000000 END
    WHEN a.q IS NULL AND a.unit IS NULL THEN 1 END END"""
presence = " ".join(f"WHEN a.term = {lit(t)} THEN {lit(f)}" for t, f in PRESENCE.items())
step("record", f"""
CREATE OR REPLACE TABLE pq_record AS
WITH ann AS (
  SELECT e.sk, e.ok, e.interaction_class, e.source, e.is_directed, e.stim, e.inh, e.dataset, a.term, a.value, a.unit, a.q
  FROM rec_ev e LEFT JOIN rec_ann a USING (resource, relation_id, ordinal)
)
SELECT entity_uuid(sk) subject_entity_id, entity_uuid(ok) object_entity_id, any_value(sk) sk, any_value(ok) ok,
  interaction_class, source, is_directed, stim is_stimulation, inh is_inhibition,
  list_sort(list_distinct(list(substr(value, 6)) FILTER (WHERE term='publications' AND regexp_full_match(value, '(?i)PMID:[0-9]+')))) reference_pubmed_ids,
  list_sort(list_distinct(list(substr(value, 5)) FILTER (WHERE term='publications' AND regexp_matches(value, '(?i)^doi:')))) reference_dois,
  min({affinity_case}) affinity,
  max(CASE WHEN term='has_quantitative_value' AND a_source_field='pchembl_value' AND regexp_full_match(value, '{NUMERIC[1:-1]}') THEN value::DOUBLE END) pchembl,
  max(CASE WHEN source='stitch' AND term='has_confidence_score' AND a_source_field='combined_score' AND regexp_full_match(value, '{NUMERIC[1:-1]}') THEN value::DOUBLE END) score,
  list_sort(list_distinct(list(CASE WHEN source='chembl' AND dataset='mechanisms' THEN 'mechanism_of_action'
       {presence} ELSE value END)
     FILTER (WHERE (source='chembl' AND dataset='mechanisms')
                OR ((term={lit(CURATION_TERM)} OR term IN {lits(PRESENCE)}) AND value IS NOT NULL AND value <> ''))))
     curation_flags
FROM (SELECT *, a.q.source_field a_source_field FROM ann a) a
GROUP BY subject_entity_id, object_entity_id, interaction_class, source, is_directed, stim, inh
""")
# Entity types and organisms. The canonical taxon of a published entity is its
# taxon when every pinned occurrence agrees.
step("entity_info", f"""
CREATE OR REPLACE TABLE ent AS
SELECT entity_key, any_value(entity_type) entity_type,
  CASE WHEN min(try_cast(taxon AS BIGINT)) = max(try_cast(taxon AS BIGINT)) THEN min(try_cast(taxon AS BIGINT)) END taxonomy_id
FROM ({' UNION ALL '.join(f"SELECT entity_key, entity_type, taxon FROM {pq(s, 'entity')}" for s in ALL)})
WHERE entity_key IN (SELECT sk FROM pq_record UNION SELECT ok FROM pq_record)
GROUP BY entity_key
""")
step("source_taxon", f"""
CREATE OR REPLACE TABLE source_taxon AS
SELECT resource AS source, entity_key,
  CASE WHEN count(DISTINCT try_cast(taxon AS BIGINT)) = 1 THEN min(try_cast(taxon AS BIGINT)) END taxonomy_id
FROM ({' UNION ALL '.join(f"SELECT '{s}' resource, entity_key, taxon FROM {pq(s, 'entity')}" for s in SCOPE)})
WHERE try_cast(taxon AS BIGINT) IS NOT NULL
GROUP BY ALL
""")
step("record_final", """
CREATE OR REPLACE TABLE pq_rec AS
SELECT r.*, se.entity_type subject_type, oe.entity_type object_type,
       coalesce(se.taxonomy_id, sst.taxonomy_id) subject_organism, coalesce(oe.taxonomy_id, ost.taxonomy_id) object_organism
FROM pq_record r LEFT JOIN ent se ON se.entity_key = r.sk LEFT JOIN ent oe ON oe.entity_key = r.ok
LEFT JOIN source_taxon sst ON sst.source = r.source AND sst.entity_key = r.sk
LEFT JOIN source_taxon ost ON ost.source = r.source AND ost.entity_key = r.ok
""")

# 3. The PostgreSQL record, as exported.
step("load_postgres", f"""
CREATE OR REPLACE TABLE pg_rec AS
SELECT subject_entity_id::UUID subject_entity_id, object_entity_id::UUID object_entity_id, interaction_class, source,
  subject_type, object_type, is_directed, is_stimulation, is_inhibition, subject_organism, object_organism,
  affinity, pchembl, score,
  CAST(curation_flags AS VARCHAR[]) curation_flags, CAST(reference_pubmed_ids AS VARCHAR[]) reference_pubmed_ids,
  CAST(reference_dois AS VARCHAR[]) reference_dois
FROM read_csv('{PG_CSV}', skip=1, header=true, auto_detect=true,
  types={{'is_directed':'BOOLEAN','is_stimulation':'BOOLEAN','is_inhibition':'BOOLEAN','subject_organism':'BIGINT',
          'object_organism':'BIGINT','affinity':'DOUBLE','pchembl':'DOUBLE','score':'DOUBLE',
          'curation_flags':'JSON','reference_pubmed_ids':'JSON','reference_dois':'JSON'}})
""")
# The PostgreSQL export holds every scope record row except ChEMBL's, of which
# it holds those with the mechanism flag. Restrict the Parquet side the same way.
step("restrict", """
CREATE OR REPLACE TABLE pq_cmp AS SELECT * FROM pq_rec
WHERE source <> 'chembl' OR list_contains(curation_flags, 'mechanism_of_action')
""")

KEY = "subject_entity_id, object_entity_id, interaction_class, source, is_directed, is_stimulation, is_inhibition"
ATTRS = ["subject_type", "object_type", "subject_organism", "object_organism", "affinity", "pchembl", "score",
         "curation_flags", "reference_pubmed_ids", "reference_dois"]


def norm(col, side):
    if col in ("curation_flags", "reference_pubmed_ids", "reference_dois"):
        return f"coalesce(list_sort({side}.{col}), [])"
    if col in ("affinity", "pchembl", "score"):
        return f"round({side}.{col}, 6)"
    return f"{side}.{col}"


report = {"scope": SCOPE, "class_resources": len(class_resources)}
q = lambda sql: con.sql(sql).fetchall()
report["record_rows"] = {"parquet": q("SELECT count(*) FROM pq_cmp")[0][0], "postgres": q("SELECT count(*) FROM pg_rec")[0][0]}
# The signature columns are nullable, so the key is compared NULL-safely.
null_safe_key = " AND ".join(f"p.{k} IS NOT DISTINCT FROM g.{k}" for k in KEY.split(", "))
con.sql(f"""CREATE OR REPLACE TABLE rec_cmp AS
SELECT p.source p_source, g.source g_source, p.interaction_class p_class, g.interaction_class g_class,
  {', '.join(f'{norm(c, "p")} p_{c}, {norm(c, "g")} g_{c}' for c in ATTRS)}
FROM pq_cmp p FULL JOIN pg_rec g ON {null_safe_key}""")
report["record_match"] = q("""SELECT
  count(*) FILTER (WHERE p_source IS NOT NULL AND g_source IS NOT NULL) both_sides,
  count(*) FILTER (WHERE g_source IS NULL) parquet_only,
  count(*) FILTER (WHERE p_source IS NULL) postgres_only FROM rec_cmp""")[0]
report["record_only_by_source"] = q("""SELECT coalesce(p_source, g_source), coalesce(p_class, g_class),
  count(*) FILTER (WHERE g_source IS NULL), count(*) FILTER (WHERE p_source IS NULL)
  FROM rec_cmp WHERE p_source IS NULL OR g_source IS NULL GROUP BY ALL ORDER BY 3 DESC, 4 DESC LIMIT 40""")
report["record_attribute_mismatches"] = {
    c: q(f"SELECT count(*) FROM rec_cmp WHERE p_source IS NOT NULL AND g_source IS NOT NULL AND p_{c} IS DISTINCT FROM g_{c}")[0][0]
    for c in ATTRS
}
report["record_attribute_examples"] = {
    c: q(f"SELECT p_source, p_{c}, g_{c} FROM rec_cmp WHERE p_source IS NOT NULL AND g_source IS NOT NULL AND p_{c} IS DISTINCT FROM g_{c} LIMIT 4")
    for c in ATTRS
}

# 4. The preset on each side: components, exclusion, then main's fold.
chem, partners = CURATED_FILTERS["entity_types"], CURATED_FILTERS["partner_types"]
gate = (f"((subject_type IN {lits(chem)} AND object_type IN {lits(partners)}) OR "
        f"(object_type IN {lits(chem)} AND subject_type IN {lits(partners)}))")


def preset(record, name):
    con.sql(f"""CREATE OR REPLACE TABLE {name}_rows AS
      SELECT * FROM {record} WHERE source IN {lits(CHEMBL_FILTERS['resources'])} AND {gate}
        AND list_has_all(coalesce(curation_flags, []), {list(CHEMBL_FILTERS['curation_flags'])})
      UNION ALL BY NAME
      SELECT * FROM {record} WHERE source IN {lits(CURATED_FILTERS['resources'])} AND {gate}
      UNION ALL BY NAME
      SELECT * FROM {record} WHERE source IN {lits(REACTION_FILTERS['resources'])} AND {gate}
        AND interaction_class IN {lits(REACTION_FILTERS['interaction_classes'])}""")
    con.sql(f"DELETE FROM {name}_rows WHERE source IN {lits(EXCLUDED)}")
    con.sql(f"""CREATE OR REPLACE TABLE {name} AS
      SELECT subject_entity_id, object_entity_id, interaction_class,
        min(subject_organism) subject_organism, min(object_organism) object_organism,
        bool_or(is_directed) is_directed, bool_or(is_stimulation) is_stimulation, bool_or(is_inhibition) is_inhibition,
        count(DISTINCT source) FILTER (WHERE is_stimulation IS NOT NULL OR is_inhibition IS NOT NULL) sign_source_count,
        count(DISTINCT source) FILTER (WHERE is_directed IS NOT NULL) direction_source_count,
        round(min(affinity), 6) affinity, round(max(pchembl), 6) pchembl, round(max(score), 6) score,
        list_sort(list_distinct(list(source))) sources, count(DISTINCT source) source_count,
        list_sort(list_distinct(flatten(list(coalesce(curation_flags, []))))) curation_flags,
        list_sort(list_distinct(flatten(list(coalesce(reference_pubmed_ids, []))))) reference_pubmed_ids,
        list_sort(list_distinct(flatten(list(coalesce(reference_dois, []))))) reference_dois
      FROM {name}_rows GROUP BY ALL""")


start = time.monotonic()
preset("pq_rec", "pq_ml")
timings["preset_parquet"] = round(time.monotonic() - start, 1)
preset("pg_rec", "pg_ml")
FOLD_KEY = "subject_entity_id, object_entity_id, interaction_class"
FOLD_ATTRS = ["subject_organism", "object_organism", "is_directed", "is_stimulation", "is_inhibition",
              "sign_source_count", "direction_source_count", "affinity", "pchembl", "score", "sources",
              "source_count", "curation_flags", "reference_pubmed_ids", "reference_dois"]
con.sql(f"""CREATE OR REPLACE TABLE ml_cmp AS SELECT p.sources p_present, g.sources g_present,
  {', '.join(f'p.{c} p_{c}, g.{c} g_{c}' for c in FOLD_ATTRS)}
  FROM pq_ml p FULL JOIN pg_ml g USING ({FOLD_KEY})""")
report["metalinksdb_pairs"] = {"parquet": q("SELECT count(*) FROM pq_ml")[0][0], "postgres": q("SELECT count(*) FROM pg_ml")[0][0]}
report["metalinksdb_match"] = q("""SELECT count(*) FILTER (WHERE p_present IS NOT NULL AND g_present IS NOT NULL),
  count(*) FILTER (WHERE g_present IS NULL), count(*) FILTER (WHERE p_present IS NULL) FROM ml_cmp""")[0]
report["metalinksdb_attribute_mismatches"] = {
    c: q(f"SELECT count(*) FROM ml_cmp WHERE p_present IS NOT NULL AND g_present IS NOT NULL AND p_{c} IS DISTINCT FROM g_{c}")[0][0]
    for c in FOLD_ATTRS
}
report["metalinksdb_identical_rows"] = q(
    "SELECT count(*) FROM ml_cmp WHERE p_present IS NOT NULL AND g_present IS NOT NULL AND "
    + " AND ".join(f"p_{c} IS NOT DISTINCT FROM g_{c}" for c in FOLD_ATTRS)
)[0][0]
report["timings_seconds"] = timings
json.dump(report, open(OUT / "report.json", "w"), indent=1, default=str)
print(json.dumps(report, indent=1, default=str)[:12000])
