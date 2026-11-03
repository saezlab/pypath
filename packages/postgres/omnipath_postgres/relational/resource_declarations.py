"""Main's explicit namespace authority declarations, vendored as metadata.

Extracted from ResourceConfig.mints at the pinned main pypath gitlink. No
resource module is imported or executed; no authority is inferred from data.
"""

SOURCE_COMMIT = "33f37fbaab59993d24f5c32bb4b3e7587085bcf4"
SOURCE_BLOBS = {
    "pypath/inputs_v2/chebi.py": "3c022fb464f571e8a240009438649728efd8a36c",
    "pypath/inputs_v2/chembl.py": "e7907856f653626e458eaf296109da796338cda2",
    "pypath/inputs_v2/hmdb.py": "cf0b49b1239cfe815d9f4c90d7777aa7d38e5b67",
    "pypath/inputs_v2/kegg.py": "93b678afa903fd3d9189eab0e189ad1249aedbcd",
    "pypath/inputs_v2/lipidmaps.py": "1c87cbba0dca805e2963dad894d8214171bf7bd9",
    "pypath/inputs_v2/swisslipids.py": "502e86b8b4bd74a1208c1ca58a25d15ea2bdde77",
}
IDENTIFIER_MINTS = {
    "chebi": ("Chebi:MI:0474",),
    "chembl": ("Chembl Compound:MI:0967",),
    "hmdb": ("Hmdb:OM:0004",),
    "kegg": ("Kegg Compound:MI:2012", "Kegg Reaction:MI:2013"),
    "lipidmaps": ("Lipidmaps:OM:0003",),
    "swisslipids": ("Swisslipids:OM:0009",),
}
