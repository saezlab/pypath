"""MappingTable -- wrapper around an ID translation dictionary."""

from __future__ import annotations

import time
from typing import NamedTuple
import logging

_log = logging.getLogger(__name__)


class MappingTableKey(NamedTuple):
    """Unique identifier for a mapping table."""

    id_type: str
    target_id_type: str
    ncbi_tax_id: int


class MappingTable:
    """Wrapper around a dict mapping source IDs to sets of target IDs.

    The core data structure is self.data: dict[str, set[str]].
    Most mappings are unambiguous (one target per source), but one-to-many
    is common (e.g. gene symbol -> multiple UniProt isoforms).

    Tables auto-expire after lifetime seconds of inactivity.
    """

    def __init__(
        self,
        data: dict[str, set[str]],
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        lifetime: int = 300,
    ):
        self.data = data
        self.id_type = id_type
        self.target_id_type = target_id_type
        self.ncbi_tax_id = ncbi_tax_id
        self.lifetime = lifetime
        self._last_used = time.time()

    @property
    def key(self) -> MappingTableKey:
        """The unique key identifying this table."""

        return MappingTableKey(
            self.id_type,
            self.target_id_type,
            self.ncbi_tax_id,
        )

    def __getitem__(self, key: str) -> set[str]:
        self._last_used = time.time()

        return self.data.get(key, set())

    def __contains__(self, key: str) -> bool:
        self._last_used = time.time()

        return key in self.data

    def __len__(self) -> int:
        return len(self.data)

    def __repr__(self) -> str:
        return (
            f'<MappingTable {self.id_type} -> {self.target_id_type}, '
            f'organism={self.ncbi_tax_id}, n={len(self)}>'
        )

    @property
    def expired(self) -> bool:
        """Whether this table has exceeded its inactivity lifetime."""

        return time.time() - self._last_used > self.lifetime

    @property
    def items(self):
        """Proxy to dict.items."""

        return self.data.items

    @property
    def keys(self):
        """Proxy to dict.keys."""

        return self.data.keys

    @property
    def values(self):
        """Proxy to dict.values."""

        return self.data.values
