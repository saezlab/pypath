-- Every veto, with what it removed. A veto that removed a derived or key-matched
-- row is the case a reader needs to see; one that removed nothing is here too,
-- because it may be stale.
INSERT INTO exclusion_metabolite_report
  (kind, exclusion_class, namespace, resolved, detail, build_id)
SELECT 'veto_applied',
       v.exclusion_class,
       v.namespace,
       count(r.entity_key)::INTEGER,
       json_object(
         'identifier', v.identifier,
         'note', v.note,
         'resolved', v.resolved,
         'entities', v.entities::JSON,
         'removed', coalesce(
           to_json(list(json_object('entity', r.entity_key,
                                    'assignment_source', r.assignment_source,
                                    'match_kind', r.match_kind)
                        ORDER BY r.entity_key, r.assignment_source, r.match_kind)
                   FILTER (WHERE r.entity_key IS NOT NULL)), '[]'::JSON)
       )::VARCHAR,
       $build_id
FROM exclusion_veto_listing v
LEFT JOIN exclusion_veto_removed r
  ON r.veto_class = v.exclusion_class AND r.namespace = v.namespace
 AND r.identifier = v.identifier
GROUP BY v.exclusion_class, v.namespace, v.identifier, v.note, v.resolved, v.entities;
