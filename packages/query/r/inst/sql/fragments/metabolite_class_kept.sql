-- Whether the relation (or endpoint table) row {{input}} passes the metabolite_classes
-- filter: one of its ends in the buckets named. Needs the relation_candidates steps.
(op_unset({{metabolite_classes}})
    OR {{input}}.relation_n IN (SELECT relation_n FROM __metabolite_class_relations))
