"""The lipid nomenclature grammar binding (spec 011 R9) -- an optional extra
(``omnipath-utils[lipid]``, pygoslin). Every function here degrades to
``None``/``False`` rather than raising when the grammar is unavailable,
mirroring :mod:`_chemistry`'s convention (research R14's capability-
degradation invariant, extended to this second optional toolkit).

**Chain counts come from the parse object, never from counting ``C:D``
tokens in a rendered string** (research R9) -- text-token counting is unsound
in two documented ways: a one-chain class (e.g. ``MG 18:1``) renders with
empty-position-slot tokens (``MG 18:1/0:0/0:0``, 3 tokens for 1 real chain),
and a sterol ester's backbone renders its own pseudo ``C:D`` token
(``CE 18:1`` -> ``SE 27:1/18:1``, 2 tokens for 1 real chain). Both are fixed
by reading ``lipid.fa_list`` directly: an empty/backbone slot has
``num_carbon == 0`` and is excluded; a genuinely named chain (whether an
individually resolved one, or a same-slot combined sum standing in for
several unresolved chains -- the ``TG 58:12_20:0`` case) has
``num_carbon > 0`` and counts. This also naturally produces the "partially
specified" level research R9 introduces (between ``species`` and
``molecular_species``): a lump-sum-plus-named-chain name like
``TG 58:12_20:0`` parses as 2 populated ``fa_list`` slots out of 3 possible
(``chains_listed=2 < chains_possible=3``), matching R9's own flagship
examples (``TG 58:12_20:0``, ``CL 54:6_18:2``) exactly.
"""

from __future__ import annotations

import re
from functools import cache

__all__ = [
    'lipid_nomenclature_available',
    'parse_lipid_name',
]

#: A C:D chain token -- present in every lipid shorthand, absent from
#: ordinary chemical names (cholesterol, glucose, ...). Used only as a cheap
#: prefilter before invoking the parser, never to count chains.
_CHAIN_TOKEN_RE = re.compile(r'\d+:\d+')

# Goslin's own LipidLevel names (lowercased) map straight through except
# SPECIES, which this project splits in two: a name with zero individually
# listed chains is 'species'; one with some-but-not-all listed is its own
# 'partially_specified' level (research R9) -- never collapsed back to
# 'species'. pygoslin has no such level itself; it is synthesised here from
# chains_listed vs. chains_possible, per data-model.md section 7's own
# stated invariant.
_PARTIALLY_SPECIFIED = 'partially_specified'


@cache
def lipid_nomenclature_available() -> bool:
    """Whether ``pygoslin`` is importable in this process. Cached -- a
    static fact about the running environment.
    """

    try:
        import pygoslin  # noqa: F401
    except ImportError:
        return False
    return True


def _goslin_version() -> str:
    try:
        import importlib.metadata

        return importlib.metadata.version('pygoslin')
    except Exception:  # noqa: BLE001
        return ''


_parser = None  # lazy singleton -- LipidParser construction is not free


def _get_parser():
    global _parser
    if _parser is None:
        from pygoslin.parser.Parser import LipidParser

        _parser = LipidParser()
    return _parser


def parse_lipid_name(name: str) -> dict | None:
    """Parse one lipid shorthand name -> the standardized identity fields, or
    ``None`` when the grammar is unavailable, the name carries no chain
    token at all (not lipid-shaped), or the grammar cannot parse it.

    Returns a dict matching ``lipid_name_node``'s columns (data-model.md
    section 7): ``lipid_name`` (the canonical rendered name -- the key),
    ``lipid_level``, ``chains_possible``, ``chains_listed``,
    ``lipid_category``, ``lipid_class``, ``total_carbon``, ``total_db``,
    ``sum_formula``, ``parser_version``.
    """

    if not lipid_nomenclature_available():
        return None
    if not name or not _CHAIN_TOKEN_RE.search(name):
        return None

    try:
        parsed = _get_parser().parse(name)
        canonical = parsed.get_lipid_string()
    except Exception:  # noqa: BLE001 -- pygoslin raises a family of parse errors
        return None

    lipid = parsed.lipid
    info = getattr(lipid, 'info', None)
    if info is None:
        return None

    level = getattr(info, 'level', None)
    level_slug = level.name.lower() if level is not None else None
    chains_possible = getattr(info, 'poss_fa', None)

    fa_list = getattr(lipid, 'fa_list', None) or []
    chains_listed = sum(1 for fa in fa_list if getattr(fa, 'num_carbon', 0) > 0)

    # The "partially specified" level (R9): fewer chains individually listed
    # than the class holds, but more than zero -- its own level, between
    # species (0 listed) and a fully resolved molecular_species+ (all listed).
    # Overrides whatever pygoslin itself called this parse (it has no such
    # level in its own vocabulary and typically calls this molecular_species).
    if (
        isinstance(chains_possible, int)
        and chains_possible > 0
        and 0 < chains_listed < chains_possible
    ):
        level_slug = _PARTIALLY_SPECIFIED

    def _try(fn):
        try:
            return fn()
        except Exception:  # noqa: BLE001
            return None

    category = _try(lambda: lipid.headgroup.lipid_category.name)
    lipid_class = _try(lambda: lipid.headgroup.get_lipid_string())
    sum_formula = _try(parsed.get_sum_formula)

    from pygoslin.domain.LipidLevel import LipidLevel

    species = _try(lambda: parsed.get_lipid_string(LipidLevel.SPECIES))
    total_carbon = total_db = None
    if species:
        chain = _CHAIN_TOKEN_RE.search(species)
        if chain:
            total_carbon, total_db = (int(x) for x in chain.group().split(':'))

    return {
        'lipid_name': canonical,
        'lipid_level': level_slug,
        'chains_possible': chains_possible,
        'chains_listed': chains_listed,
        'lipid_category': category,
        'lipid_class': lipid_class,
        'total_carbon': total_carbon,
        'total_db': total_db,
        'sum_formula': sum_formula,
        'parser_version': _goslin_version(),
    }
