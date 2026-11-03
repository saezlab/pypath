"""Tests for the structure-bearing chemical backends + export-resolver (Milestone I).

The adapter/registration/policy tests are fast (no DB, no downloads). The
export-resolver protein projection runs only when OMNIPATH_UTILS_DB_URL points at
a built utils DB.
"""

from __future__ import annotations

import os

import pytest


# ── inputs_v2 adapter: the pure projection logic ────────────────────────────


def test_build_mapping_projects_columns():
    from omnipath_utils.mapping.backends._inputs_v2_adapter import build_mapping

    rows = [
        {'chebi_id': 'CHEBI:15377', 'inchikey': 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'},
        {'chebi_id': 'CHEBI:16236', 'inchikey': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'},
        {'chebi_id': 'CHEBI:NOSTRUCT', 'inchikey': None},  # dropped
    ]
    mapping = build_mapping(rows, 'chebi_id', 'inchikey')
    assert mapping == {
        'CHEBI:15377': {'XLYOFNOQVPJJNP-UHFFFAOYSA-N'},
        'CHEBI:16236': {'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'},
    }


def test_build_mapping_handles_list_cells():
    from omnipath_utils.mapping.backends._inputs_v2_adapter import build_mapping

    rows = [{'lm_id': 'LMFA', 'xref': ['A', 'B', '']}]
    assert build_mapping(rows, 'lm_id', 'xref') == {'LMFA': {'A', 'B'}}


# ── backend registration + id_types.yaml column wiring ──────────────────────


@pytest.mark.parametrize(
    'backend_name, yaml_key, module, dataset',
    [
        ('chebi', 'chebi', 'chebi', 'molecules'),
        ('chembl', 'chembl', 'chembl', 'molecules'),
        ('lipidmaps', 'lipidmaps', 'lipidmaps', 'lipids'),
        ('swisslipids', 'swisslipids', 'swisslipids', 'lipids'),
    ],
)
def test_adapter_backends_registered(backend_name, yaml_key, module, dataset):
    from omnipath_utils.mapping.backends import get_backend

    backend = get_backend(backend_name)
    assert backend is not None
    assert backend.yaml_key == yaml_key
    assert backend.resource_module == module
    assert backend.dataset == dataset


def test_pubchem_backend_registered():
    from omnipath_utils.mapping.backends import get_backend

    assert get_backend('pubchem') is not None


@pytest.mark.parametrize(
    'backend, structure_col',
    [
        ('chebi', 'inchikey'),
        ('chembl', 'standard_inchi_key'),
        ('lipidmaps', 'INCHI_KEY'),
        ('pubchem', 'inchikey'),
    ],
)
def test_inchikey_wired_for_backends(backend, structure_col):
    """id_types.yaml resolves inchikey to each backend's raw column."""
    from omnipath_utils.mapping._id_types import IdTypeRegistry

    reg = IdTypeRegistry.get()
    assert reg.backend_column('inchikey', backend) == structure_col


def test_chebi_source_column_wired():
    from omnipath_utils.mapping._id_types import IdTypeRegistry

    reg = IdTypeRegistry.get()
    assert reg.backend_column('chebi', 'chebi') == 'chebi_id'
    assert reg.backend_column('smiles', 'chebi') == 'smiles'


# ── export-resolver: policy + projection ────────────────────────────────────


def test_load_policy_and_accepted_sources():
    from pathlib import Path

    from omnipath_utils.db import _resolver_export as rx

    policy_path = (
        Path(rx.__file__).parents[1] / 'data' / 'resolver_policy.yaml'
    )
    policy = rx.load_policy(str(policy_path))
    assert 'chemical' in policy and 'protein' in policy
    # Canonical inchikey is excluded; candidate_only (swisslipids) not in accept.
    sources = rx._accepted_source_types(policy, 'chemical', 'inchikey')
    assert 'chebi' in sources and 'pubchem' in sources
    assert 'swisslipids' not in sources


def test_canonical_type_resolution():
    from omnipath_utils.db import _resolver_export as rx

    assert rx._canonical_type('chemical', 'Standard InChI Key') == 'inchikey'
    assert rx._canonical_type('protein', None) == 'uniprot'
    with pytest.raises(ValueError):
        rx._canonical_type('chemical', 'bogus-type')


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_DB_URL'),
    reason='OMNIPATH_UTILS_DB_URL not set; export-resolver needs a built utils DB',
)
def test_export_resolver_protein(tmp_path):
    from pathlib import Path

    from omnipath_utils.db import _resolver_export as rx

    policy_path = Path(rx.__file__).parents[1] / 'data' / 'resolver_policy.yaml'
    stats = rx.export_resolver(
        family='protein',
        policy_path=str(policy_path),
        output_dir=str(tmp_path),
        taxa=[9606],
        max_records=500,
    )
    assert stats.canonical_type == 'uniprot'
    assert stats.files and stats.rows > 0
    # Parquet schema: source_type, source_id, canonical_target, taxonomy_id (text).
    import pyarrow.parquet as pq

    table = pq.read_table(stats.files[0])
    assert table.column_names == [
        'source_type',
        'source_id',
        'canonical_target',
        'taxonomy_id',
    ]
    assert table.schema.field('taxonomy_id').type == __import__('pyarrow').string()
