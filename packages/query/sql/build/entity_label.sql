-- The lowercase label of each entity that has one, in the order a substring search ranks
-- them (shorter labels, more relations, then the label and key): rank_n numbers that
-- order and the table is stored in it, so that a search for the entities whose label
-- contains some text stops at the first it needs, and reads one narrow column instead of
-- lowering every label.
SELECT row_number() OVER (ORDER BY length(label), relation_count DESC, label, entity_key)::INTEGER AS rank_n,
    entity_n, lower(label) AS term
FROM omnipath.entity
WHERE label IS NOT NULL
