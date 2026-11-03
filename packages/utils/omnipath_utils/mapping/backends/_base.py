"""Base class for mapping backends."""

from __future__ import annotations

from abc import ABC, abstractmethod
import logging

from omnipath_utils.mapping._id_types import IdTypeRegistry

_log = logging.getLogger(__name__)


class MappingBackend(ABC):
    """Abstract base for ID mapping data sources.

    Subclasses declare their backend name and YAML key. The base class
    handles:

    - ID type column resolution via IdTypeRegistry
    - Support checking with debug log when columns are not found
    - pypath vs direct HTTP dispatch with ImportError fallback
    - Info logging at start and end of each read
    """

    #: Human-readable name for log messages (e.g. "uniprot", "biomart").
    name: str = ''

    #: Key in the ``backends`` dict inside id_types.yaml
    #: (e.g. "uniprot", "ensembl").
    yaml_key: str = ''

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def read(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        """Read mapping data.

        Resolves backend-specific column names, checks support,
        then dispatches to pypath or direct HTTP.

        Returns:
            Dict mapping source IDs to sets of target IDs.
        """

        reg = IdTypeRegistry.get()
        src_col = reg.backend_column(id_type, self.yaml_key)
        tgt_col = reg.backend_column(target_id_type, self.yaml_key)

        if not src_col or not tgt_col:
            _log.debug(
                '%s backend does not support %s -> %s',
                self.name,
                id_type,
                target_id_type,
            )
            return {}

        _log.info(
            '%s: loading %s -> %s (organism %d)',
            self.name,
            id_type,
            target_id_type,
            ncbi_tax_id,
        )

        try:
            data = self._read_via_pypath(
                id_type,
                target_id_type,
                ncbi_tax_id,
                src_col=src_col,
                tgt_col=tgt_col,
                **kwargs,
            )
        except ImportError:
            _log.debug(
                '%s: pypath not available, using direct HTTP',
                self.name,
            )
            data = self._read_direct(
                id_type,
                target_id_type,
                ncbi_tax_id,
                src_col=src_col,
                tgt_col=tgt_col,
                **kwargs,
            )

        _log.info(
            '%s: loaded %d entries for %s -> %s (organism %d)',
            self.name,
            len(data),
            id_type,
            target_id_type,
            ncbi_tax_id,
        )

        return data

    # ------------------------------------------------------------------
    # Hooks for subclasses
    # ------------------------------------------------------------------

    @abstractmethod
    def _read_via_pypath(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        *,
        src_col: str,
        tgt_col: str,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        """Backend-specific pypath loading.

        Must raise ``ImportError`` if pypath is not available so
        the base class can fall back to :meth:`_read_direct`.
        """

        ...

    @abstractmethod
    def _read_direct(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        *,
        src_col: str,
        tgt_col: str,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        """Backend-specific direct HTTP loading."""

        ...
