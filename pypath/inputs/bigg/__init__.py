"""
BiGG Models metabolite cross-reference data.

Provides functions to download and parse the BiGG universal metabolite
cross-reference table, which maps BiGG metabolite IDs to external
databases (ChEBI, HMDB, KEGG, MetaNetX, LipidMaps, etc.).

The data is sourced from
``http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt``
(a ~2 MB TSV covering all 9,090 universal metabolites across 85+ models).
"""

from pypath.inputs.bigg._metabolites import *  # noqa: F401, F403
