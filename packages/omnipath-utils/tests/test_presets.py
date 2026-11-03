"""T040 / C8 -- build presets: a fast chemical-only preset and a complete full.

The structural checks run anywhere; the capped-build assertion is DB-gated.
"""

import os

import pytest


class TestPresetConfig:
    def test_chemical_preset_exists(self):
        from omnipath_utils.db._presets import PRESETS

        assert 'chemical' in PRESETS
        cfg = PRESETS['chemical']
        assert cfg['metabolite'] is True          # loads the chemical stack
        assert cfg['mappings'] == []              # no protein/gene mappings
        assert not cfg['mirna']
        assert not cfg['orthology']
        assert not cfg.get('ftp')

    def test_full_preset_runs_metabolites(self):
        from omnipath_utils.db._presets import PRESETS

        # the long-value chemical phase runs in every full build (no silently
        # skipped phase, Principle V): full -> metabolite -> _populate_chemical_long
        assert PRESETS['full']['metabolite'] is True


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_CHEMICAL_BUILD_DB'),
    reason='set OMNIPATH_UTILS_CHEMICAL_BUILD_DB to run the capped chemical build',
)
class TestChemicalPresetBuild:
    def test_capped_chemical_build_has_no_protein_rows(self):
        import time

        from sqlalchemy import text
        from sqlalchemy.orm import Session

        from omnipath_utils.db._build import DatabaseBuilder

        url = os.environ['OMNIPATH_UTILS_CHEMICAL_BUILD_DB']
        builder = DatabaseBuilder(db_url=url, max_records=20000)
        start = time.time()
        builder.build_preset('chemical')
        elapsed = time.time() - start
        assert elapsed < 15 * 60, f'capped chemical build took {elapsed:.0f}s'

        with Session(builder.engine) as s:
            protein = s.execute(text("""
                SELECT count(*) FROM omnipath_utils.id_mapping m
                JOIN omnipath_utils.id_type st ON st.id = m.source_type_id
                WHERE st.entity_type IN ('protein','gene')
            """)).scalar()
        assert protein == 0
