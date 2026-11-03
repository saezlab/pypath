"""Entity resolution against the reference library.

Thin stateful wrapper around :class:`omnipath_resolver.canonical.LibraryMatcher`
that keeps the interface the writer and pipeline use
(``resolve_entity_targets`` / ``resolution_stats``). Identifier lookup caches live in
the matcher; attributed observation metadata is merged downstream.

The library is located, in order, from the explicit ``library_dir`` argument,
``$OMNIPATH_LIBRARY_DIR``, or ``data/reference/library`` next to the output
root the pipeline passes in.  Without a library every gene/protein and
chemical observation is *unmatched* and keeps its own identifier.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .canonical import LibraryMatcher, Match, ResolutionTracker, get_policy
from .canonical.library import pin_library
from .contracts import RawEntityObservation

logger = logging.getLogger(__name__)


@dataclass(slots=True)
class ResolvedEntityInfo:
    """Canonical identity of one observed entity."""

    canonical_namespace: str
    canonical_identifier: str
    label: str
    aliases: dict[str, list[str]] = field(default_factory=dict)
    taxon: str | None = None
    resolved_by: str = "unmatched"
    node_id: str | None = None
    protein_namespace: str | None = None
    protein_identifier: str | None = None
    protein_node_id: str | None = None
    protein_label: str | None = None
    protein_aliases: dict[str, list[str]] = field(default_factory=dict)
    protein_taxon: str | None = None
    protein_gene_candidates: tuple[str, ...] = ()
    gene_mapping_status: str | None = None
    gene_candidates: tuple[str, ...] = ()
    transcript_namespace: str | None = None
    transcript_identifier: str | None = None

    @property
    def matched(self) -> bool:
        return self.node_id is not None


@dataclass(slots=True)
class ResolvedEntityTarget(ResolvedEntityInfo):
    entity_type: str | None = None
    reference_library: str | None = None


def _target(info: ResolvedEntityInfo) -> ResolvedEntityTarget:
    return ResolvedEntityTarget(
        canonical_namespace=info.canonical_namespace,
        canonical_identifier=info.canonical_identifier,
        label=info.label,
        aliases=info.aliases,
        taxon=info.taxon,
        resolved_by=info.resolved_by,
        node_id=info.node_id,
        **_molecular_fields(info),
    )


def locate_library_dir(
    explicit: str | Path | None = None, data_root: str | Path | None = None
) -> Path | None:
    for candidate in (explicit, os.environ.get("OMNIPATH_LIBRARY_DIR")):
        if candidate:
            return pin_library(candidate)
    if data_root is not None:
        return pin_library(Path(data_root) / "reference" / "library")
    return None


class EntityResolver:
    """Resolve raw entity observations to library nodes."""

    def __init__(
        self,
        library_dir: str | Path | None = None,
        *,
        data_root: str | Path | None = None,
        defer_aliases: bool = False,
        memory_limit: str | None = None,
    ) -> None:
        self.library_dir = locate_library_dir(library_dir, data_root)
        self.matcher = LibraryMatcher(
            self.library_dir, defer_aliases=defer_aliases, memory_limit=memory_limit
        )
        if self.library_dir is not None and not self.matcher.libraries:
            logger.warning(
                "No reference library found under %s; entities keep their native identifiers",
                self.library_dir,
            )
        self._tracker = ResolutionTracker()

    @property
    def libraries(self) -> list[str]:
        return self.matcher.libraries

    def resolve_entities(
        self,
        entities: dict[str, RawEntityObservation],
        progress: bool = True,
    ) -> dict[str, ResolvedEntityInfo]:
        # Scalar base lookup; serving builds use resolve_entity_targets below.
        matches = self.matcher.match(entities)
        results: dict[str, ResolvedEntityInfo] = {}
        for key, obs in entities.items():
            ent_type = (
                str(obs.entity_type or "unknown")
                .strip()
                .lower()
                .replace(" ", "_")
                .replace("-", "_")
            )
            policy = get_policy(ent_type)
            match = matches[key]
            info = _to_info(match)
            # Candidate/node caching belongs to the matcher. Observation-derived
            # metadata must not accumulate here: that makes results history-dependent.
            results[key] = info
            self._tracker.record(key, ent_type, policy, match)
        return results

    def resolve_entity_targets(
        self,
        entities: dict[str, RawEntityObservation],
        progress: bool = True,
    ) -> dict[str, list[ResolvedEntityTarget]]:
        result = {}
        for key, matches in self.matcher.targets(entities).items():
            targets = []
            for match in matches:
                info = _target(_to_info(match))
                info.entity_type = match.entity_type
                info.reference_library = match.reference_library
                targets.append(info)
            result[key] = targets
            obs = entities[key]
            self._tracker.record(key, obs.entity_type, get_policy(obs.entity_type), matches[0])
        return result

    def resolution_stats(self) -> dict[str, Any]:
        stats = self._tracker.summary()
        stats["library_dir"] = str(self.library_dir) if self.library_dir else None
        stats["libraries"] = self.libraries
        stats["lookup_metrics"] = dict(self.matcher.metrics)
        return stats

    def export_resolution_keys(self, path: str | Path) -> Path:
        return self._tracker.export_keys(path)

    def close(self) -> None:
        self._tracker.close()
        self.matcher.close()


def _to_info(match: Match) -> ResolvedEntityInfo:
    return ResolvedEntityInfo(
        canonical_namespace=match.canonical_namespace,
        canonical_identifier=match.canonical_identifier,
        label=match.label,
        aliases={ns: list(vals) for ns, vals in match.aliases.items()},
        taxon=match.taxon,
        resolved_by=match.resolved_by,
        node_id=match.node_id,
        **_molecular_fields(match),
    )


def _molecular_fields(info):
    return {
        name: getattr(info, name)
        for name in (
            "protein_namespace",
            "protein_identifier",
            "protein_node_id",
            "protein_label",
            "protein_aliases",
            "protein_taxon",
            "protein_gene_candidates",
            "gene_mapping_status",
            "gene_candidates",
            "transcript_namespace",
            "transcript_identifier",
        )
    }
