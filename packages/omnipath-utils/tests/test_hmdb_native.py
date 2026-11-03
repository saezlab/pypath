"""The native HMDB load: required cross-references, not a best-effort extra.

HMDB carries the richest per-metabolite cross-reference block in the field
(ChEBI, KEGG, PubChem, DrugBank, FooDB), and cycle 005 recorded it as
unreachable behind Cloudflare. The archive is now mirrored on the project's
own host (``pypath.inputs_v2.hmdb``, ``rescued.omnipathdb.org``), so the load
must actually reach every one of those namespaces, not just the narrow
synonym-to-ChEBI slice the legacy path produced -- and a failure must fail
the build, not warn and continue quietly.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

REQUIRED_TARGET_TYPES = {'chebi', 'kegg', 'pubchem', 'drugbank', 'foodb'}


def _builder():
    """A DatabaseBuilder with no database connection.

    ``DatabaseBuilder.__init__`` connects and calls ``ensure_schema`` on
    construction, which this test has no business doing -- it exercises the
    row-shaping logic in ``_long_hmdb``, not the database. Only the
    attributes ``_effective_limit`` reads are set.
    """
    from omnipath_utils.db._build import DatabaseBuilder

    builder = object.__new__(DatabaseBuilder)
    builder._max_records = None
    builder._pubchem_max_records = None
    return builder


_HMDB_ROW = {
    'accession': 'HMDB0000001',
    'name': '1-Methylhistidine',
    'synonyms': 'Pi-methylhistidine;1-MHis',
    'chebi_id': '50599',
    'pubchem_compound_id': '92105',
    'kegg_id': 'C01152',
    'drugbank_id': 'DB04151',
    'foodb_id': 'FDB012154',
    'inchikey': 'BRMWTNUJHUMWMS-LURJTMIESA-N',
    'inchi': None,
    'smiles': None,
}


def test_long_hmdb_produces_every_required_cross_reference():
    """The native load reaches ChEBI, KEGG, PubChem, DrugBank and FooDB.

    Not just the synonym-to-ChEBI slice the legacy ``pypath.inputs.hmdb``
    path produced -- every namespace HMDB's own record carries.
    """
    builder = _builder()
    captured_target_types: set[str] = set()

    def fake_populate_long_slice(data, source_type, target_type, backend_name):
        if data:
            captured_target_types.add(target_type)
        return len(data)

    with (
        patch(
            'omnipath_utils.mapping.backends._inputs_v2_adapter.raw_rows',
            return_value=[_HMDB_ROW],
        ),
        patch.object(
            builder,
            '_populate_long_slice',
            side_effect=fake_populate_long_slice,
        ),
    ):
        builder._long_hmdb()

    missing = REQUIRED_TARGET_TYPES - captured_target_types
    assert not missing, (
        f'_long_hmdb never wrote a slice targeting {missing}. The native '
        'load must reach every cross-reference HMDB publishes, not a '
        'narrow subset.'
    )


def test_hmdb_load_is_required_not_best_effort():
    """A failing HMDB load must fail the build, not warn and continue.

    ``_populate_chemical_long`` currently wraps every per-resource builder
    in one shared try/except that only logs a warning. HMDB is the richest
    cross-reference source in the chemical layer -- silently continuing
    without it is exactly the kind of gap cycle 005 already got away with
    once.
    """
    from omnipath_utils.db._build import DatabaseBuilder

    builder = _builder()

    with (
        patch.object(builder, '_long_chebi'),
        patch.object(builder, '_long_hmdb', side_effect=RuntimeError('mirror unreachable')),
        patch.object(builder, '_long_chembl'),
        patch.object(builder, '_long_kegg'),
        patch.object(builder, '_long_ramp'),
        patch.object(DatabaseBuilder, '_record_long_rollups', lambda self: None),
        pytest.raises(RuntimeError, match='mirror unreachable'),
    ):
        builder._populate_chemical_long()
