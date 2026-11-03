"""Build-refusing checks for the chemical key-space join (spec 011 T042/T044).

Two guards WP1 exists to add:

* :func:`check_key_space_overlap` -- a namespace present in both the
  ingested evidence and the resolver lookup, yet sharing no key, is exactly
  the shape of this cycle's defect (spec.md US1 acceptance scenario 6). Run
  it every build. WP1's own normalization landing is no guard against a
  *future* namespace drifting the same way.
* :func:`check_absorption_threshold` -- one canonical identity absorbing
  an unreasonable share of source records (spec.md US1 acceptance scenario
  5). This repeats the placeholder-merge failure from the lipid hierarchy
  work, now in the chemical key space.

Both raise ``RuntimeError`` naming every offending namespace/identifier at
once, so the build fails loudly instead of silently losing resolution.
"""

from __future__ import annotations


def check_key_space_overlap(
    evidence_keys: dict[str, set[str]],
    resolver_keys: dict[str, set[str]],
) -> None:
    """Fail when a namespace present in ``evidence_keys`` shares no key
    with ``resolver_keys``.

    A namespace absent from ``evidence_keys`` (nothing ingested for it)
    is not a failure -- there is nothing to resolve, which is a different
    situation from a namespace present on both sides that fails to match.
    """

    offending = sorted(
        namespace
        for namespace, keys in evidence_keys.items()
        if keys and not (keys & resolver_keys.get(namespace, set()))
    )
    if offending:
        raise RuntimeError(
            'resolution preflight: no key overlap between the ingested '
            'evidence and the resolver lookup for namespace(s): '
            + ', '.join(offending)
        )


def check_absorption_threshold(
    counts: dict[str, int],
    threshold: int,
) -> None:
    """Fail when one canonical identity's source-record count exceeds
    ``threshold``."""

    offending = sorted(
        identifier
        for identifier, count in counts.items()
        if count > threshold
    )
    if offending:
        raise RuntimeError(
            f'resolution absorption threshold ({threshold}) exceeded by '
            'identifier(s): ' + ', '.join(offending)
        )
