"""Fill ``data_source_license`` from pypath's resource metadata.

The build database recorded no license anywhere. `data_source` held
``(source_id, name)`` and nothing more, so "everything usable commercially" was
a question no query could ask. This module answers it by loading the terms
pypath already maintains next to its resource definitions.

**Where the terms come from — two sources, one precedence.**

``pypath/resources/data/resources.json`` names a license per resource — a
*name*, such as ``CC BY 4.0``. It is hand-curated one resource at a time, with
the grounds cited in its commit messages, so it is the original license
curation rather than a derived artifact. It **wins wherever it has an entry**.

It covers pypath's *interaction* resources only, and the build has since loaded
chemistry, ontology and single-cell resources it never described. Those are
reached through the second source: every ``inputs_v2`` resource declares a
:class:`pypath.internals.cv_terms.LicenseCV` term — ``license`` is a required
field on ``ResourceConfig`` — and ``sync_resources_table`` already writes that
term's accession into ``resources.license``. The CV term **fills the gaps** and
never overrides a curated entry, because for a resource ``resources.json`` does
not describe there is no competing curation to lose.

Either way the name is resolved against ``pypath/data/licenses/*.json`` through
:class:`pypath.resources.licenses.Licenses`, which is where the ordinal levels
live: ``purpose``, ``sharing`` and ``attrib``, each a point on a scale from the
smallest to the greatest freedom (``pypath/internals/license.py``). The levels
are read off that model rather than transcribed here, so a license whose terms
are corrected in pypath is corrected in the build at the next load.

**Where the two disagree, the load says so.** Precedence makes the behaviour
defined; it does not make the disagreement resolved. A resource whose curated
entry names one concrete license and whose CV term names another is logged at
WARNING and reported in :class:`LicenseTableStats.conflicts`, because that is a
curation question for a person to settle and not one for this module to decide
quietly.

**The name mapping** is the awkward part on the ``resources.json`` side and is
not invented here. Its keys are display names such as ``CellChatDB``, while
``data_source.name`` holds slugs such as ``cellchat``. The build already
resolves the two through the three-name model
(``pypath.inputs_v2.resource_names``), whose filter index maps every
{slug, slugified short name, slugified synonym} onto one canonical slug, and
whose registry maps that slug back to the ``resources.json`` key. The candidate
names a source offers the index are its own slug and the short/full/synonym
names ``sync_resources_table`` already resolved into the ``resources`` table —
the same mapping, read rather than rebuilt.

**Unknown terms are exclusions, never defaults.** A source the mapping
does not reach gets a row all the same, carrying ``is_known = false`` and no
levels at all. Nothing here invents a level for it: a NULL level is silently
dropped by a range predicate, which is the right answer for the wrong reason,
so the exclusion is carried by ``is_known`` and never by the numbers. A license
that maps but whose level pypath cannot place on its scale is treated the same
way — pypath scores an unplaceable level 99, which every ``>=`` filter admits,
and admitting it is exactly the failure mode this table cannot detect
downstream. ``LicenseCV.UNSPECIFIED`` says the same thing in the vocabulary's
own words and is treated the same: it is the absence of terms, never a
permissive default, and it never turns a resource known.

Call it once per build, and **after** ``sync_resources_table``: it reads both
the resolved resource names and the declared CV term out of the ``resources``
table that step fills. Run before it, and every resource falls back to the
curated entry alone.

.. code-block:: python

    license_stats = sync_data_source_licenses(conn, schema=args.schema)
"""

from __future__ import annotations

import json
import logging
import warnings
from functools import lru_cache
from dataclasses import field, dataclass

from psycopg2 import sql
import psycopg2.extensions

logger = logging.getLogger(__name__)

# pypath scores a level it cannot place on its ordinal scale as 99, which sits
# above every real level and would therefore pass every `>=` filter.
UNPLACEABLE_LEVEL = 99

# The vocabulary's own word for "no terms recorded".
UNSPECIFIED_TERM = 'OM:0599'

# `LicenseCV` accession -> the key of the record in pypath's license database
# (`pypath/data/licenses/*.json`, keyed by each file's `name`). This is an
# identity mapping between two spellings of the same license and nothing more:
# the levels are still read off the pypath record, never written here.
#
# Four members are deliberately absent, because pypath has no license record
# for them and inventing one would be inventing terms:
#
#   OM:0505 ACADEMIC_FREE — "free for academic use"; `afl_v3.json` is the
#           Academic Free License v3.0, a software license, not this.
#   OM:0512 BIGG          — BiGG Models terms, no record.
#   OM:0515 PUBLIC        — public domain; CC0 is one instrument of it, not a
#           synonym for it.
#   OM:0599 UNSPECIFIED   — the vocabulary's own word for unknown.
#
# A resource carrying one of these stays `is_known = false`. The fix belongs in
# `pypath/data/licenses/`, not here.
CV_LICENSE_NAMES = {
    'OM:0501': 'CC BY 4.0',
    'OM:0502': 'CC0 1.0',
    'OM:0503': 'GPLv3',
    'OM:0504': 'MIT',
    'OM:0506': 'HPO',
    'OM:0507': 'GPLv2',
    'OM:0508': 'KEGG',
    'OM:0509': 'CC BY-SA 3.0',
    'OM:0510': 'CC BY-SA 4.0',
    'OM:0511': 'CC BY-NC 4.0',
    'OM:0513': 'BSD',
    'OM:0514': 'CC BY 3.0',
}


@dataclass(frozen=True)
class LicenseTableStats:
    """What one ``data_source_license`` load recorded.

    ``unmapped`` names the sources that got no license, because an unmapped
    source is a finding to surface rather than a number to round up: it is
    excluded from every license-filtered result, and the caller deserves to
    know which resources it loses.

    ``conflicts`` names the sources whose curated ``resources.json`` entry and
    whose ``inputs_v2`` CV term both name a concrete license and the two are
    not the same one, as ``(source, curated name, CV name)``. The curated entry
    was stored; the disagreement is a curation finding and is reported rather
    than resolved here.
    """

    sources: int = 0
    known: int = 0
    unmapped: tuple[str, ...] = field(default_factory=tuple)
    from_curation: int = 0
    from_cv: int = 0
    conflicts: tuple[tuple[str, str, str], ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class ResourceLicense:
    """One resource's ``data_source_license`` terms, as ordinal levels."""

    license_name: str | None = None
    license_full_name: str | None = None
    license_url: str | None = None
    purpose_level: int | None = None
    sharing_level: int | None = None
    attrib_level: int | None = None
    is_known: bool = False
    resource_key: str | None = None
    matched_on: str | None = None
    # Which of the two sources supplied the terms: ``curation`` for the
    # hand-curated ``resources.json`` entry, ``cv`` for the ``LicenseCV`` term
    # the resource declares in ``inputs_v2``.
    origin: str | None = None


_UNKNOWN = ResourceLicense()


@lru_cache(maxsize=1)
def _resources_json() -> dict:
    """The authoritative resource metadata, the file the controller reads."""

    try:
        from importlib import resources as importlib_resources

        path = importlib_resources.files('pypath.resources.data') / 'resources.json'
        with path.open('r', encoding='utf-8') as handle:
            return json.load(handle)
    except Exception:  # noqa: BLE001 - metadata missing or unreadable
        logger.warning('pypath resources.json is unavailable; no license loaded')
        return {}


@lru_cache(maxsize=1)
def _license_db() -> object | None:
    """The pypath license database: license name → :class:`License`."""

    try:
        from pypath.resources.licenses import Licenses

        return Licenses()
    except Exception:  # noqa: BLE001 - pypath without the license database
        logger.warning('pypath license database is unavailable')
        return None


@lru_cache(maxsize=1)
def _cv_terms() -> dict[str, object]:
    """The license vocabulary, as ``accession -> LicenseCV`` member."""

    try:
        from pypath.internals.cv_terms import LicenseCV
    except ImportError:  # pypath without the CV terms (older pin)
        logger.warning('pypath LicenseCV is unavailable; no CV term read')
        return {}

    return {str(term.value): term for term in LicenseCV}


@lru_cache(maxsize=1)
def _name_mapping() -> tuple[dict[str, str], dict[str, str]]:
    """The name mapping, as ``(filter index, canonical slug → JSON key)``.

    ``build_filter_index`` already folds slug, short name and synonyms onto one
    canonical slug; the registry carries that slug back to the self-spelled
    short name, which *is* the ``resources.json`` key.
    """

    try:
        from pypath.inputs_v2.resource_names import (
            resource_registry,
            build_filter_index,
        )
    except ImportError:  # pypath without the 3-name model (older pin)
        logger.warning('pypath resource_names is unavailable; no license loaded')
        return {}, {}

    registry = resource_registry()
    return build_filter_index(), {
        slug: names.short for slug, names in registry.items()
    }


def _slugify(name: str) -> str:
    from pypath.inputs_v2.resource_names import slugify

    return slugify(name)


def _levels(license_obj: object) -> tuple[int, int, int] | None:
    """The three ordinal levels, or ``None`` when pypath cannot place one.

    ``LicenseFeature.level_to_int`` warns and returns 99 for a level outside
    its vocabulary. 99 is above every real level, so storing it would admit the
    resource under every filter — the opposite of what an unrecognised license
    means.
    """

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            levels = (
                int(license_obj.purpose),
                int(license_obj.sharing),
                int(license_obj.attrib),
            )
        except (AttributeError, TypeError, ValueError):
            return None
    return None if UNPLACEABLE_LEVEL in levels else levels


def _from_record(
    name: str,
    record: object | None,
    *,
    origin: str,
    resource_key: str | None = None,
    matched_on: str | None = None,
) -> ResourceLicense:
    """One license record, as a row: with levels when pypath can place them.

    A record whose levels the model cannot place keeps its name for
    attribution and stays unknown. An unplaceable level scores 99, which sits
    above every real level and would pass every ``>=`` filter — the opposite
    of what an unrecognised license means.
    """

    levels = _levels(record) if record is not None else None
    if levels is None:
        logger.warning(
            'license %r of resource %r carries no usable levels',
            name,
            resource_key or matched_on,
        )
        return ResourceLicense(
            license_name=name,
            resource_key=resource_key,
            matched_on=matched_on,
            origin=origin,
        )
    purpose, sharing, attrib = levels
    return ResourceLicense(
        license_name=str(record),
        license_full_name=getattr(record, 'full_name', None),
        license_url=getattr(record, 'url', None),
        purpose_level=purpose,
        sharing_level=sharing,
        attrib_level=attrib,
        is_known=True,
        resource_key=resource_key,
        matched_on=matched_on,
        origin=origin,
    )


def resolve_curated_license(candidates: list[str | None]) -> ResourceLicense:
    """The hand-curated ``resources.json`` entry, if one of ``candidates`` hits.

    ``candidates`` are the names one source answers to, most canonical first:
    its build slug, then the short/full/synonym names the three-name model
    resolved for it. Returns the unknown license when none of them reaches an entry.
    """

    index, key_of = _name_mapping()
    if not index:
        return _UNKNOWN
    metadata = _resources_json()
    db = _license_db()

    for candidate in candidates:
        if not candidate:
            continue
        slug = index.get(_slugify(candidate))
        if not slug:
            continue
        key = key_of.get(slug)
        name = (metadata.get(key) or {}).get('license') if key else None
        if not isinstance(name, str):
            continue
        record = db[name] if db is not None else None
        return _from_record(
            name,
            record,
            origin='curation',
            resource_key=key,
            matched_on=candidate,
        )

    return _UNKNOWN


def resolve_cv_license(accession: str | None) -> ResourceLicense:
    """The ``LicenseCV`` term the resource declares in ``inputs_v2``.

    ``accession`` is the OM term ``sync_resources_table`` copied out of
    ``ResourceConfig.license`` into ``resources.license``. ``UNSPECIFIED``, and
    any member pypath has no license record for, resolve to unknown — with the
    accession and the term's own words kept for attribution, so that a row a
    filter must exclude still says what it was excluded for.
    """

    if not accession:
        return _UNKNOWN

    accession = accession.strip()
    term = _cv_terms().get(accession)
    name = CV_LICENSE_NAMES.get(accession)

    if name is None:
        if term is None:
            logger.warning('unknown license CV accession %r', accession)
            return _UNKNOWN
        # A term the vocabulary has and pypath has no record for: state it,
        # store no level for it.
        return ResourceLicense(
            license_name=accession,
            license_full_name=getattr(term, 'definition', None),
            license_url=getattr(term, 'source', None) or None,
            origin='cv',
            matched_on=accession,
        )

    db = _license_db()
    record = db[name] if db is not None else None
    return _from_record(name, record, origin='cv', matched_on=accession)


def resolve_license(
    candidates: list[str | None],
    accession: str | None = None,
) -> tuple[ResourceLicense, tuple[str, str] | None]:
    """The terms to store for one resource, and any disagreement behind them.

    The curated ``resources.json`` entry wins wherever it has one; the CV term
    fills the gap where it has none. The second element is ``(curated name, CV
    name)`` when both name a concrete license and the two differ — the caller
    reports it, and the curated entry is what was stored.
    """

    curated = resolve_curated_license(candidates)
    from_cv = resolve_cv_license(accession)

    # `UNSPECIFIED` is the absence of terms, so it disagrees with nothing. Any
    # other CV term names something concrete, including the ones pypath has no
    # record for: `PUBLIC` against a curated `CC0 1.0` is a real disagreement
    # about which instrument applies, and is worth a person's attention.
    cv_concrete = bool(from_cv.license_name) and accession != UNSPECIFIED_TERM

    conflict = None
    if (
        curated.license_name
        and cv_concrete
        and curated.license_name != from_cv.license_name
    ):
        conflict = (curated.license_name, from_cv.license_name)

    if curated.is_known:
        return curated, conflict
    if from_cv.is_known:
        # No curation for this resource, or curation the level model cannot
        # place: the resource's own declared term carries it.
        return from_cv, conflict
    # Neither source reaches usable levels. Keep whichever name exists, for
    # attribution only — `is_known` is false either way.
    return (curated if curated.license_name else from_cv), conflict


def _source_candidates(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> list[tuple[int, str, list[str], str | None]]:
    """Every ``data_source`` row with the names it answers to and its CV term.

    Both come from the ``resources`` table, which ``sync_resources_table``
    fills from ``ResourceConfig``: the names through the three-name model,
    the license accession straight off ``ResourceConfig.license``. Reading them
    here keeps one resolution of each in the build rather than two, and is why
    this module needs no ``inputs_v2`` import — the declared term is already in
    the database, for every resource the build loaded. A database without that
    table (a scratch schema, an older build) falls back to the slug alone and
    to the curated entry.
    """

    cur.execute(
        sql.SQL('SELECT to_regclass({})').format(
            sql.Literal(f'{schema}.resources')
        )
    )
    has_resources = cur.fetchone()[0] is not None

    if has_resources:
        cur.execute(
            sql.SQL(
                """
                SELECT ds.source_id,
                       ds.name,
                       r.resource_short,
                       r.resource_full,
                       r.synonyms,
                       r.license
                FROM {}.data_source AS ds
                LEFT JOIN {}.resources AS r ON r.resource_id = ds.name
                ORDER BY ds.name
                """
            ).format(sql.Identifier(schema), sql.Identifier(schema))
        )
        return [
            (
                source_id,
                name,
                [name, short, full, *(synonyms or ())],
                accession,
            )
            for source_id, name, short, full, synonyms, accession
            in cur.fetchall()
        ]

    cur.execute(
        sql.SQL('SELECT source_id, name FROM {}.data_source ORDER BY name').format(
            sql.Identifier(schema)
        )
    )
    return [(source_id, name, [name], None) for source_id, name in cur.fetchall()]


def sync_data_source_licenses(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> LicenseTableStats:
    """Load one ``data_source_license`` row per resource.

    Every source gets a row, the unmapped ones included: a row saying
    ``is_known = false`` states the exclusion, where a missing row leaves it to
    be inferred. The load is an upsert, so a resource whose terms changed is
    corrected in place and one that lost its mapping loses its levels with it.
    """

    with conn.cursor() as cur:
        sources = _source_candidates(cur, schema)
        rows = []
        unmapped: list[str] = []
        conflicts: list[tuple[str, str, str]] = []
        origins = {'curation': 0, 'cv': 0}
        for source_id, name, candidates, accession in sources:
            resolved, conflict = resolve_license(candidates, accession)
            if conflict is not None:
                curated_name, cv_name = conflict
                # Stored the curated entry, and said so: which of the two is
                # right is a curation question, not one this load can settle.
                logger.warning(
                    'license of %r disagrees between the curated entry (%s) '
                    'and the declared CV term (%s); stored the curated one',
                    name,
                    curated_name,
                    cv_name,
                )
                conflicts.append((name, curated_name, cv_name))
            if not resolved.is_known:
                unmapped.append(name)
            elif resolved.origin in origins:
                origins[resolved.origin] += 1
            rows.append(
                (
                    source_id,
                    resolved.license_name,
                    resolved.license_full_name,
                    resolved.license_url,
                    resolved.purpose_level,
                    resolved.sharing_level,
                    resolved.attrib_level,
                    resolved.is_known,
                )
            )

        cur.executemany(
            sql.SQL(
                """
                INSERT INTO {}.data_source_license (
                  source_id,
                  license_name,
                  license_full_name,
                  license_url,
                  purpose_level,
                  sharing_level,
                  attrib_level,
                  is_known
                )
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (source_id) DO UPDATE SET
                  license_name = EXCLUDED.license_name,
                  license_full_name = EXCLUDED.license_full_name,
                  license_url = EXCLUDED.license_url,
                  purpose_level = EXCLUDED.purpose_level,
                  sharing_level = EXCLUDED.sharing_level,
                  attrib_level = EXCLUDED.attrib_level,
                  is_known = EXCLUDED.is_known
                """
            )
            .format(sql.Identifier(schema))
            .as_string(cur.connection),
            rows,
        )
    conn.commit()

    if unmapped:
        # Named rather than counted: each of these is dropped from every
        # license-filtered result, and silence about which ones would leave
        # that exclusion unauditable.
        logger.warning(
            'no license mapped for %d of %d sources: %s',
            len(unmapped),
            len(rows),
            ', '.join(sorted(unmapped)),
        )

    return LicenseTableStats(
        sources=len(rows),
        known=len(rows) - len(unmapped),
        unmapped=tuple(sorted(unmapped)),
        from_curation=origins['curation'],
        from_cv=origins['cv'],
        conflicts=tuple(sorted(conflicts)),
    )
