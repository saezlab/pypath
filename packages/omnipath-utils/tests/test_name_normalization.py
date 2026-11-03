"""spec 011 T071 -- one shared name-normalization function.

Before this, ``_query.py``'s ``_lookup_key`` and ``_build.py``'s
``_long_rows`` each wrote their own ``.strip().lower()`` for a name-like
value -- harmless while identical, but two places to drift. Both now call
``normalize_name``.
"""

from omnipath_utils.mapping._id_types import normalize_name


def test_case_folded():
    assert normalize_name('Taurine') == 'taurine'
    assert normalize_name('TAURINE') == 'taurine'


def test_whitespace_folded():
    assert normalize_name('  GABA  ') == 'gaba'


def test_query_and_build_sides_share_the_function():
    import inspect

    from omnipath_utils.db import _query, _build

    assert 'normalize_name' in inspect.getsource(_query._lookup_key)
    assert 'normalize_name' in inspect.getsource(
        _build.DatabaseBuilder._long_rows
    )
