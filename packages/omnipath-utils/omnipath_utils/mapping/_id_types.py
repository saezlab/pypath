"""ID type registry -- loads and queries id_types.yaml."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any
import logging
import re

from pkg_infra.data import load

_log = logging.getLogger(__name__)


class IdTypeRegistry:
    """Registry of biological identifier types.

    Loads id_types.yaml and provides lookup by name, alias,
    backend column name, or entity type.
    """

    _instance: IdTypeRegistry | None = None

    def __init__(self) -> None:
        self._types: dict[str, dict] = load(
            'id_types.yaml',
            module='omnipath_utils',
        )
        self._build_indices()

    @classmethod
    def get(cls) -> IdTypeRegistry:
        """Singleton access."""

        if cls._instance is None:
            cls._instance = cls()

        return cls._instance

    def _build_indices(self) -> None:
        """Build lookup indices for fast queries."""

        self._by_alias: dict[str, str] = {}
        self._by_backend: dict[str, dict[str, str]] = {}

        for name, info in self._types.items():
            # Index aliases
            for alias in info.get('aliases', []):
                self._by_alias[alias] = name

            self._by_alias[name] = name

            # Also index with underscores/dashes normalized
            self._by_alias[name.replace('-', '_')] = name
            self._by_alias[name.replace('_', '-')] = name

            # Index backend columns -> our name
            for backend, col in info.get('backends', {}).items():
                if col:
                    self._by_backend.setdefault(backend, {})[col] = name

    def resolve(self, name: str) -> str | None:
        """Resolve an ID type name or alias to canonical name.

        Args:
            name: ID type name or alias.

        Returns:
            Canonical name, or None if not found.
        """

        return self._by_alias.get(name) or self._by_alias.get(name.lower())

    def info(self, name: str) -> dict[str, Any] | None:
        """Get full info dict for an ID type.

        Args:
            name: ID type name or alias.

        Returns:
            Info dict with label, entity_type, curie_prefix, backends,
            aliases; or None if not found.
        """

        canonical = self.resolve(name)

        return self._types.get(canonical) if canonical else None

    def entity_type(self, name: str) -> str | None:
        """Get entity type for an ID type.

        Args:
            name: ID type name or alias.

        Returns:
            Entity type string, or None.
        """

        info = self.info(name)

        return info.get('entity_type') if info else None

    def curie_prefix(self, name: str) -> str | None:
        """Get Bioregistry CURIE prefix.

        Args:
            name: ID type name or alias.

        Returns:
            CURIE prefix string, or None.
        """

        info = self.info(name)

        return info.get('curie_prefix') if info else None

    def url_pattern(self, name: str) -> str | None:
        """Get a URL template that resolves an identifier to its web page.

        An explicit ``url_pattern`` in the registry wins (use it to point at a
        provider's native page). Otherwise, when a Bioregistry CURIE prefix is
        known, the Bioregistry resolver template is derived so every typed
        identifier gets a working link with no per-type curation. ``{$id}`` is
        the placeholder a consumer replaces with the identifier.

        Args:
            name: ID type name or alias.

        Returns:
            URL template string, or None when neither an explicit pattern nor a
            CURIE prefix is available.
        """

        info = self.info(name)

        if not info:
            return None

        explicit = info.get('url_pattern')

        if explicit:
            return explicit

        prefix = info.get('curie_prefix')

        return f'https://bioregistry.io/{prefix}:{{$id}}' if prefix else None

    def id_pattern(self, name: str) -> str | None:
        """Get a regular expression a valid identifier of this type matches.

        Only returned when the registry declares one explicitly.

        Args:
            name: ID type name or alias.

        Returns:
            Regular-expression string, or None.
        """

        info = self.info(name)

        return info.get('id_pattern') if info else None

    def backend_column(self, name: str, backend: str) -> str | None:
        """Get the backend-specific column name for an ID type.

        Args:
            name: ID type name or alias.
            backend: Backend name (e.g. 'uniprot', 'ensembl').

        Returns:
            Backend column name, or None.
        """

        info = self.info(name)

        if info:
            return info.get('backends', {}).get(backend)

        return None

    def from_backend_column(
        self,
        backend: str,
        column: str,
    ) -> str | None:
        """Look up our canonical name from a backend column name.

        Args:
            backend: Backend name.
            column: Backend-specific column name.

        Returns:
            Canonical ID type name, or None.
        """

        return self._by_backend.get(backend, {}).get(column)

    def by_entity_type(self, entity_type: str) -> list[str]:
        """List all ID types for a given entity type.

        Args:
            entity_type: One of protein, gene, transcript,
                small_molecule, mirna, probe.

        Returns:
            List of canonical ID type names.
        """

        return [
            name
            for name, info in self._types.items()
            if info.get('entity_type') == entity_type
        ]

    def by_backend(self, backend: str) -> dict[str, str]:
        """Map of our ID type names to backend column names for a backend.

        Args:
            backend: Backend name.

        Returns:
            Dict mapping canonical name to backend column name.
        """

        return {
            name: info['backends'][backend]
            for name, info in self._types.items()
            if backend in info.get('backends', {}) and info['backends'][backend]
        }

    def all_names(self) -> list[str]:
        """List all canonical ID type names."""

        return list(self._types.keys())

    def __contains__(self, name: str) -> bool:
        return self.resolve(name) is not None

    def __len__(self) -> int:
        return len(self._types)

    def __repr__(self) -> str:
        return f'<IdTypeRegistry [{len(self)} types]>'


# One canonical form per chemical namespace (spec 011 T024-T026). WP1's
# resolver joins the translation database and the build's own evidence on
# this normalized key. The two sides must derive the identical string from
# every raw form a source actually publishes. That covers prefix presence
# and case, and per-namespace padding or its absence. HMDB has 5- and
# 7-digit accessions. KEGG zero-pads its compound number. PubChem's
# canonical CID never carries leading zeros. One implementation, here, is
# what keeps the two sides from silently drifting apart the way the
# un-normalized key already did once.

_VALUE_PATTERNS: dict[str, str] = {
    'chebi': r'CHEBI:\d+',
    'hmdb': r'HMDB\d{7}',
    'swisslipids': r'SLM:\d+',
    'lipidmaps': r'LM[A-Z0-9]+',
    'kegg': r'C\d{5}',
    'pubchem': r'[1-9]\d*|0',
    'pubchem_substance': r'[1-9]\d*|0',
    'chembl': r'CHEMBL\d+',
    'cas': r'\d{2,7}-\d{2}-\d',
    'inchikey': r'[A-Z]{14}-[A-Z]{10}-[A-Z]',
}


def _normalize_chebi(value: str) -> str | None:
    m = re.fullmatch(r'(?:chebi:)?\s*(\d+)', value.strip(), re.IGNORECASE)
    return f'CHEBI:{m.group(1)}' if m else None


def _normalize_hmdb(value: str) -> str | None:
    m = re.fullmatch(r'hmdb(\d+)', value.strip(), re.IGNORECASE)
    return f'HMDB{m.group(1).zfill(7)}' if m else None


def _normalize_swisslipids(value: str) -> str | None:
    m = re.fullmatch(r'(?:slm:)?\s*(\d+)', value.strip(), re.IGNORECASE)
    return f'SLM:{m.group(1)}' if m else None


def _normalize_lipidmaps(value: str) -> str | None:
    v = value.strip()
    return v.upper() if re.fullmatch(r'lm[a-z0-9]+', v, re.IGNORECASE) else None


def _normalize_kegg(value: str) -> str | None:
    m = re.fullmatch(r'c(\d+)', value.strip(), re.IGNORECASE)
    return f'C{m.group(1).zfill(5)}' if m else None


def _normalize_pubchem(value: str) -> str | None:
    m = re.fullmatch(r'0*(\d+)', value.strip())
    return m.group(1) if m else None


def _normalize_chembl(value: str) -> str | None:
    v = value.strip().replace(' ', '')
    m = re.fullmatch(r'chembl(\d+)', v, re.IGNORECASE)
    return f'CHEMBL{m.group(1)}' if m else None


def _normalize_cas(value: str) -> str | None:
    v = value.strip()
    return v if re.fullmatch(r'\d{2,7}-\d{2}-\d', v) else None


def _normalize_inchikey(value: str) -> str | None:
    m = re.fullmatch(
        r'(?:inchikey=)?([a-z]{14}-[a-z]{10}-[a-z])', value.strip(), re.IGNORECASE,
    )
    return m.group(1).upper() if m else None


_NORMALIZERS: dict[str, Callable[[str], str | None]] = {
    'chebi': _normalize_chebi,
    'hmdb': _normalize_hmdb,
    'swisslipids': _normalize_swisslipids,
    'lipidmaps': _normalize_lipidmaps,
    'kegg': _normalize_kegg,
    'pubchem': _normalize_pubchem,
    'pubchem_substance': _normalize_pubchem,
    'chembl': _normalize_chembl,
    'cas': _normalize_cas,
    'inchikey': _normalize_inchikey,
}


def normalize_identifier(id_type: str, value: str | None) -> str | None:
    """The one normalized form of ``value`` for the chemical namespace
    ``id_type``. Every consumer, on both sides of the join, calls this
    instead of writing its own prefix/case/padding rule.

    Returns ``None`` for an unregistered namespace, an empty value, or a
    value that does not match that namespace's raw form at all.
    """

    if not value:
        return None

    canonical = IdTypeRegistry.get().resolve(id_type) or id_type
    normalizer = _NORMALIZERS.get(canonical)

    return normalizer(value) if normalizer else None


def normalize_name(value: str) -> str:
    """The one normalized form of a chemical name or synonym (spec 011 R5).

    Case- and whitespace-folded so the same molecule under different
    capitalization or incidental spacing collapses to one lookup key.
    Every writer into ``id_mapping_long`` and every reader of it (the
    'name' axis in :mod:`omnipath_utils.db._query`) calls this instead of
    writing its own ``.strip().lower()``.
    """

    return str(value).strip().lower()


def value_pattern(id_type: str) -> str | None:
    """The regular expression a namespace's *normalized* form must match.

    Distinct from :meth:`IdTypeRegistry.id_pattern`, which validates a raw,
    pre-normalization value when the registry declares one.
    """

    canonical = IdTypeRegistry.get().resolve(id_type) or id_type

    return _VALUE_PATTERNS.get(canonical)
