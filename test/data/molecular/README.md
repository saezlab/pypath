# Molecular source fixtures

`uniprot-P04637.tsv` is the official UniProt REST TSV for P04637, retrieved
5 October 2026. Its source is https://rest.uniprot.org/uniprotkb/P04637?format=tsv
with accession, sequence, sequence_version, and explicit feature fields selected.
UniProt data are CC BY 4.0. Tests use this export to validate source serialization;
synthetic fixtures cover difficult cases without asserting that they occur here.

BRENDA tests use synthetic records following the official 2026.1 README
(https://www.brenda-enzymes.org/download.php, dl-readme). EN and PM are engineering
and posttranslation modification; #protein# identifiers are EC-local. ChEMBL
synthetic SQL fixtures follow the release-36 official schema, where variants link
to assays and sequences are representative reconstructions.

`bindingdb-202610-header-row.tsv` contains the official header and first row
from https://www.bindingdb.org/rwd/bind/downloads/BindingDB_All_202610_tsv.zip,
retrieved 5 October 2026 (CC BY 3.0). It verifies the actual chain sequence
column spelling against both parser paths.
