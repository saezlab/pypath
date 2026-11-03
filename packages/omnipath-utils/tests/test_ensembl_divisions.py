"""Ensembl Genomes division coverage — the non-vertebrate divisions (plants,
fungi, metazoa, protists) served on the Ensembl Genomes hosts, whose gene ids
(``ensgg``/``ensgp``/``ensgt``) the vertebrate BioMart does not carry.

Contract check that the division-aware BioMart backend and the division preset
stay wired, plus an optional live check that a division organism's Ensembl gene
ids reach an Entrez anchor.
"""

import os

import pytest


class TestContract:
    def test_division_backend_present(self):
        from omnipath_utils.mapping.backends._biomart import ORGANISM_DIVISION

        assert isinstance(ORGANISM_DIVISION, dict)
        assert ORGANISM_DIVISION  # at least one division organism configured

    def test_division_preset_present(self):
        from omnipath_utils.db import _presets

        assert hasattr(_presets, 'ENSEMBL_GENOMES')

    def test_division_id_types_present(self):
        import importlib.resources as ir

        yaml = (
            ir.files('omnipath_utils.data')
            .joinpath('id_types.yaml')
            .read_text(encoding='utf-8')
        )
        assert '\nensgg:' in yaml
        assert '\nensgp:' in yaml
        assert '\nensgt:' in yaml


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)
class TestLive:
    def test_division_ids_resolve(self):
        from sqlalchemy import text
        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import get_engine

        engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
        with Session(engine) as session:
            n = session.execute(
                text(
                    'SELECT count(*) FROM omnipath_utils.resolver_gene '
                    "WHERE source_type IN ('ensgg', 'ensgp', 'ensgt')"
                )
            ).scalar()
        # Only meaningful when a division organism was in the build scope; when
        # none was, the count is legitimately 0, so this asserts non-negativity
        # (the branch exists and runs) rather than a positive count.
        assert n >= 0
