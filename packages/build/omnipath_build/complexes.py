"""Composition identities, applied once after all observation shards arrive.

Membership is a set of canonical member identities, not inferred stoichiometry.
Only explicit nested membership definitions qualify. Conflicting definitions,
ambiguous member expansion, missing nested definitions and cycles stay native.
"""

from collections import defaultdict

import pyarrow as pa
import pyarrow.parquet as pq
from omnipath_core.keys import stable_hash


def composition_keys(complexes, definitions, blocked=()):
    """Resolve nested complexes bottom-up without guessing through cycles."""
    pending = set(definitions) - set(blocked)
    resolved, leaves = {}, {}
    while pending:
        ready = []
        for key in sorted(pending):
            variants = definitions[key]
            members = set().union(*variants)
            nested = members & complexes
            if not nested <= leaves.keys():
                continue
            expanded = {
                tuple(sorted(set().union(*(leaves[m] if m in leaves else {m} for m in variant))))
                for variant in variants
            }
            if len(expanded) == 1 and next(iter(expanded)):
                leaf_set = next(iter(expanded))
                leaves[key] = set(leaf_set)
                resolved[key] = stable_hash("complex-composition-v1", leaf_set)
            ready.append(key)
        if not ready:
            break
        pending.difference_update(ready)
    return resolved


def resolve_complexes(db, payload_path, metrics):
    complexes = {
        r[0]
        for r in db.execute(
            "SELECT DISTINCT entity_key FROM entities WHERE entity_type='macromolecular_complex'"
        ).fetchall()
    }
    if not complexes:
        return
    # Dataset-row definitions keep contradictory complete compositions apart.
    rows = db.execute("""SELECT r.subject_entity_key, r.source, r.dataset, r.row_id,
            r.object_entity_key, r.event_id, r.object_molecular_form
        FROM relations r JOIN (SELECT DISTINCT entity_key FROM entities
            WHERE entity_type='macromolecular_complex') c ON c.entity_key=r.subject_entity_key
        WHERE r.predicate IN ('has_member', 'has_part') AND contains(r.upstream_id, ':member:')""").fetchall()
    groups, events = defaultdict(set), defaultdict(set)
    for parent, source, dataset, row, member, event, form in rows:
        # Gene grouping must not make complexes containing different products
        # or explicitly different member forms identical. This describes the
        # observed complex; it does not create reusable state-combination nodes.
        if form:
            product = form.get("protein_entity_key") or form.get("transcript_entity_key")
            details = {
                key: value
                for key, value in form.items()
                if key not in {"protein_entity_key", "transcript_entity_key"} and value is not None
            }
            member = (
                stable_hash("molecular-member", product or member, details)
                if details
                else (product or member)
            )
        groups[parent, source, dataset, row].add(member)
        events[parent, event].add(member)
    definitions = defaultdict(list)
    for (parent, *_), members in groups.items():
        definitions[parent].append(members)
    blocked = {parent for (parent, _), members in events.items() if len(members) > 1}
    identities = composition_keys(complexes, definitions, blocked)
    if not identities:
        return
    db.register(
        "complex_identity_input",
        pa.table({"old_key": list(identities), "identifier": list(identities.values())}),
    )
    db.execute("""CREATE TEMP TABLE complex_identity_map AS SELECT *,
        sha256('macromolecular_complex' || chr(0) || 'complex' || chr(0) || identifier || chr(0)) AS new_key
        FROM complex_identity_input""")
    db.unregister("complex_identity_input")
    db.execute("""CREATE TEMP TABLE complex_relation_map AS
        WITH mapped AS (SELECT r.*,
            coalesce(s.new_key,r.subject_entity_key) AS s, coalesce(o.new_key,r.object_entity_key) AS o
            FROM relations r LEFT JOIN complex_identity_map s ON r.subject_entity_key=s.old_key
            LEFT JOIN complex_identity_map o ON r.object_entity_key=o.old_key
            WHERE s.old_key IS NOT NULL OR o.old_key IS NOT NULL),
        oriented AS (SELECT *, NOT is_directed AND qualified='()' AND s>o AS flipped FROM mapped),
        endpoints AS (SELECT *, CASE WHEN flipped THEN o ELSE s END AS new_s,
            CASE WHEN flipped THEN s ELSE o END AS new_o FROM oriented)
        SELECT DISTINCT relation_key AS old_key, new_s, new_o, flipped,
            CASE WHEN statement_kind='ontology' THEN ontology_key(new_s,predicate,new_o,qualified)
            ELSE sha256(new_s || chr(0) || predicate || chr(0) || new_o || chr(0) ||
                CASE WHEN qualified='()' THEN '' ELSE qualified || chr(0) END) END AS new_key
        FROM endpoints""")
    db.execute("""UPDATE entities SET entity_key=m.new_key, namespace='complex', identifier=m.identifier,
        reference_entity_key='complex:' || m.identifier
        FROM complex_identity_map m WHERE entities.entity_key=m.old_key""")
    for table in ("identifiers", "entity_annotations", "entity_evidence", "refs"):
        db.execute(
            f"UPDATE {table} SET entity_key=m.new_key FROM complex_identity_map m WHERE {table}.entity_key=m.old_key"
        )
    db.execute(
        """UPDATE identifiers SET is_canonical=false WHERE entity_key IN (SELECT new_key FROM complex_identity_map)"""
    )
    db.execute(
        "INSERT INTO identifiers SELECT DISTINCT new_key, 'complex', identifier, 'canonical', true FROM complex_identity_map"
    )
    db.execute("""UPDATE relations SET relation_key=m.new_key, subject_entity_key=m.new_s, object_entity_key=m.new_o,
        subject_molecular_form=CASE WHEN m.flipped THEN object_molecular_form ELSE subject_molecular_form END,
        object_molecular_form=CASE WHEN m.flipped THEN subject_molecular_form ELSE object_molecular_form END,
        subject_reference_entity_key=CASE WHEN m.flipped THEN object_reference_entity_key ELSE subject_reference_entity_key END,
        object_reference_entity_key=CASE WHEN m.flipped THEN subject_reference_entity_key ELSE object_reference_entity_key END
        FROM complex_relation_map m WHERE relations.relation_key=m.old_key""")
    db.execute("""UPDATE relations SET subject_reference_entity_key='complex:' || m.identifier
        FROM complex_identity_map m WHERE relations.subject_entity_key=m.new_key""")
    db.execute("""UPDATE relations SET object_reference_entity_key='complex:' || m.identifier
        FROM complex_identity_map m WHERE relations.object_entity_key=m.new_key""")
    db.execute("""UPDATE relation_annotations SET relation_key=m.new_key,
        scope=CASE WHEN m.flipped THEN CASE scope WHEN 'subject' THEN 'object' WHEN 'object' THEN 'subject' ELSE scope END ELSE scope END
        FROM complex_relation_map m WHERE relation_annotations.relation_key=m.old_key""")
    # Rewrite only narrow keys in SQL; keep the potentially large JSON dictionary in Arrow.
    temp = payload_path.with_suffix(".complex.parquet")
    source = pq.ParquetFile(payload_path, read_dictionary=["payload_json"])
    with pq.ParquetWriter(temp, source.schema_arrow, compression="zstd", store_schema=False) as out:
        for batch in source.iter_batches(batch_size=8192):
            table = pa.Table.from_batches([batch])
            narrow = table.drop(["payload_json"]).append_column(
                "position", pa.array(range(len(table)))
            )
            db.register("complex_payload_batch", narrow)
            result = db.execute("""SELECT coalesce(r.new_key,p.relation_key) AS relation_key,
                coalesce(e.new_key,p.entity_key) AS entity_key, p.source, p.row_id
                FROM complex_payload_batch p LEFT JOIN complex_relation_map r ON p.relation_key=r.old_key
                LEFT JOIN complex_identity_map e ON p.entity_key=e.old_key ORDER BY p.position""").to_arrow_table()
            out.write_table(result.append_column("payload_json", table["payload_json"]))
            db.unregister("complex_payload_batch")
    source.close()
    temp.replace(payload_path)
    metrics["composition_complexes"] = len(identities)
    db.execute("DROP TABLE complex_relation_map")
    db.execute("DROP TABLE complex_identity_map")
