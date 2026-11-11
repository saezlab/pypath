-- Every relation key numbered in key order (relation_n).
SELECT relation_key, row_number() OVER (ORDER BY relation_key)::INTEGER AS relation_n
FROM (SELECT DISTINCT relation_key FROM work.relation WHERE relation_key IS NOT NULL)
