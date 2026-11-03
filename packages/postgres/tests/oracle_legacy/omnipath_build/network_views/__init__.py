"""Declarative network-view framework + the onboarded networks."""

from __future__ import annotations

from omnipath_build.network_views._framework import (
    NetworkDefinition,
    NetworkViewStats,
    apply_all,
    apply_network,
    ensure_network_registry,
    refresh_all,
    refresh_network,
    register_network,
)
from omnipath_build.network_views._definitions import (
    NETWORKS,
    METALINKSDB,
    LIANA,
    REACTIONS,
)

__all__ = [
    'NetworkDefinition',
    'NetworkViewStats',
    'apply_all',
    'apply_network',
    'ensure_network_registry',
    'refresh_all',
    'refresh_network',
    'register_network',
    'NETWORKS',
    'METALINKSDB',
    'LIANA',
    'REACTIONS',
]
