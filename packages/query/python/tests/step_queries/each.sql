-- query: each
-- param: keys VARCHAR[] = NULL
-- output: json response
-- example: {"keys": ["c", "a", "zz"]}
-- example: {"keys": []}
-- example: {"keys": ["k0", "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11", "k12", "k13", "k14", "k15", "k16", "k17", "k18", "k19", "k20", "k21", "k22", "k23", "k24", "k25", "k26", "k27", "k28", "k29", "k30", "k31", "k32", "k33", "k34", "k35", "k36", "k37", "k38", "k39", "k40", "k41", "k42", "k43", "k44", "k45", "k46", "k47", "k48", "k49", "k50", "k51", "k52", "k53", "k54", "k55", "k56", "k57", "k58", "k59", "k60", "k61", "k62", "k63", "k64", "a", "f"]}

-- step: rows AS table
SELECT key, label FROM item WHERE key IN ({{keys+}})

-- step: response AS value
SELECT to_json({'labels': coalesce(list(label ORDER BY key), [])}) FROM __rows
