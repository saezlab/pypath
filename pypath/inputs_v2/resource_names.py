"""Resource 3-name model: slug / short / full + forgiving filter resolution (Milestone M).

Every resource carries three names:

- **slug**  — all-lowercase canonical identifier & filter key; no ``_``, no spaces.
- **short** — the resource's own exact spelling (display name); no ``_``, no spaces.
- **full**  — the long name; spaces allowed; no ``_``.

``_`` is reserved for ``primary_secondary`` composite labels (e.g.
``PhosphoSite_SIGNOR``) and must not appear in any single resource's names.

Names are sourced from the **authoritative resource metadata**, not hand-curated:
the resource registry (:func:`resource_registry`) is loaded from
``pypath/resources/data/resources.json`` — the same file
:class:`pypath.resources.controller.ResourceController` reads — where each entry
is keyed by the resource's **self-spelled short name** and carries ``full_name``
/ ``label`` / ``synonyms``. ``inputs_v2`` resources may additionally pin explicit
``slug``/``short``/``full`` on their :class:`ResourceConfig`. Resources not found
derive a slug from their ``name`` and fall back to it for short/full; the
build-time validator flags any rule violation. (The legacy
``pypath.resources.descriptions`` module is NOT used.)
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from functools import lru_cache

__all__ = [
    'ResourceNames',
    'resource_registry',
    'slugify',
    'validate_resource_name',
    'resolve_names',
    'resolve_filter',
    'build_filter_index',
]

_SLUG_RE = re.compile(r'^[a-z0-9]+$')
_NO_UNDERSCORE_SPACE_RE = re.compile(r'^[^_\s]+$')


@dataclass(frozen=True)
class ResourceNames:
    slug: str
    short: str
    full: str
    synonyms: tuple[str, ...] = ()


def slugify(name: str) -> str:
    """Derive a canonical slug: lowercase, keep only ``[a-z0-9]``."""
    return re.sub(r'[^a-z0-9]', '', (name or '').lower())


def validate_resource_name(slug: str, short: str, full: str) -> list[str]:
    """Return a list of rule violations (empty == valid)."""
    errors: list[str] = []
    if not slug or not _SLUG_RE.match(slug):
        errors.append(f'slug {slug!r} must be lowercase alphanumeric (no _, no spaces)')
    if not short or not _NO_UNDERSCORE_SPACE_RE.match(short):
        errors.append(f'short {short!r} must contain no _ and no spaces')
    if not full or '_' in full:
        errors.append(f'full {full!r} must contain no _')
    return errors


@lru_cache(maxsize=1)
def _resources_json() -> dict:
    """Load the authoritative resource metadata (same file the controller reads)."""
    try:
        from importlib import resources as importlib_resources

        path = importlib_resources.files('pypath.resources.data') / 'resources.json'
        with path.open('r', encoding='utf-8') as handle:
            return json.load(handle)
    except Exception:  # pragma: no cover - missing/unreadable metadata
        return {}


@lru_cache(maxsize=1)
def resource_registry() -> dict[str, ResourceNames]:
    """``slug → ResourceNames`` built from ``resources.json`` (authoritative).

    The JSON key is the resource's self-spelled **short**; ``full_name`` (else
    ``label``, else the key) is the **full** name; ``synonyms`` carry over.
    """
    registry: dict[str, ResourceNames] = {}
    for key, entry in _resources_json().items():
        slug = slugify(key)
        if not slug or slug in registry:
            continue
        entry = entry if isinstance(entry, dict) else {}
        full = entry.get('full_name') or entry.get('label') or key
        synonyms = tuple(s for s in (entry.get('synonyms') or ()) if s)
        registry[slug] = ResourceNames(
            slug=slug, short=key, full=full, synonyms=synonyms
        )
    return registry


def resolve_names(
    *,
    name: str,
    slug: str | None = None,
    short: str | None = None,
    full: str | None = None,
    synonyms: tuple[str, ...] = (),
) -> ResourceNames:
    """Resolve the 3 names: explicit fields → ``resources.json`` registry → derived.

    ``slug`` (when given) is the lookup key — pass the canonical slug (e.g. from a
    ``ResourceCv`` member name) so an inconsistent ``name`` does not mislead the
    lookup; otherwise the slug is derived from ``short``/``name``.
    """
    lookup_slug = slug or slugify(short or name)
    entry = resource_registry().get(lookup_slug)
    resolved_short = short or (entry.short if entry else (name or lookup_slug))
    resolved_full = full or (entry.full if entry else resolved_short)
    merged_synonyms = tuple(
        dict.fromkeys((*synonyms, *(entry.synonyms if entry else ())))
    )
    return ResourceNames(
        slug=lookup_slug,
        short=resolved_short,
        full=resolved_full,
        synonyms=merged_synonyms,
    )


def build_filter_index(
    resources: list[ResourceNames] | None = None,
) -> dict[str, str]:
    """Map every {slug, slugified short, slugified synonym} → canonical slug."""
    index: dict[str, str] = {}
    for names in (
        resources if resources is not None else resource_registry().values()
    ):
        index[names.slug] = names.slug
        index.setdefault(slugify(names.short), names.slug)
        for synonym in names.synonyms:
            index.setdefault(slugify(synonym), names.slug)
    return index


def resolve_filter(
    query: str,
    index: dict[str, str] | None = None,
) -> str | None:
    """Resolve a resource filter argument (slug/short/synonym, case-insensitive) → slug."""
    if not query:
        return None
    lookup = index if index is not None else build_filter_index()
    return lookup.get(slugify(query))
