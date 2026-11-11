-- Every entity key numbered in key order (entity_n), relation ends without an entity row
-- included: the release's tables name an entity by 4 bytes instead of its 64-character
-- key.
SELECT entity_key, row_number() OVER (ORDER BY entity_key)::INTEGER AS entity_n
FROM (SELECT entity_key FROM work.entity WHERE entity_key IS NOT NULL
      UNION SELECT entity_key FROM work.endpoint_stubs)
