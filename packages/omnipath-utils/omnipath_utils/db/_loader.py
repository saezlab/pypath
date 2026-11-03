"""On-demand background loading of mapping tables."""

from __future__ import annotations

import logging
import threading

_log = logging.getLogger(__name__)

_pending: set[tuple] = set()  # tables currently being loaded
_lock = threading.Lock()


def request_table(
    source_type: str,
    target_type: str,
    ncbi_tax_id: int,
    db_url: str | None = None,
):
    """Request a mapping table to be loaded in the background.

    If the table is already being loaded or queued, this is a no-op.
    Returns True if a new load was queued, False if already pending.
    """
    key = (source_type, target_type, ncbi_tax_id)

    with _lock:
        if key in _pending:
            return False
        _pending.add(key)

    _log.info(
        'Queuing background load: %s -> %s (organism %d)',
        source_type,
        target_type,
        ncbi_tax_id,
    )

    thread = threading.Thread(
        target=_load_table,
        args=(source_type, target_type, ncbi_tax_id, db_url),
        daemon=True,
    )
    thread.start()
    return True


def is_pending(source_type: str, target_type: str, ncbi_tax_id: int) -> bool:
    """Check if a table is currently being loaded."""
    return (source_type, target_type, ncbi_tax_id) in _pending


def _load_table(
    source_type: str,
    target_type: str,
    ncbi_tax_id: int,
    db_url: str | None,
):
    """Background worker: load a single mapping table into the DB."""
    key = (source_type, target_type, ncbi_tax_id)

    try:
        from omnipath_utils.db._build import DatabaseBuilder

        builder = DatabaseBuilder(db_url=db_url)

        # Try to find a backend for this pair
        from omnipath_utils.mapping._mapper import Mapper

        mapper = Mapper()
        backends = mapper._find_backends(source_type, target_type)

        if not backends:
            _log.warning(
                'No backend found for %s -> %s', source_type, target_type
            )
            return

        # Try each backend
        for backend_name in backends:
            try:
                builder.populate_mapping(
                    source_type,
                    target_type,
                    ncbi_tax_id,
                    backend_name,
                )
                _log.info(
                    'Background load complete: %s -> %s (org %d) via %s',
                    source_type,
                    target_type,
                    ncbi_tax_id,
                    backend_name,
                )
                return
            except Exception as e:
                _log.debug('Backend %s failed: %s', backend_name, e)

        _log.warning(
            'All backends failed for %s -> %s (org %d)',
            source_type,
            target_type,
            ncbi_tax_id,
        )

    except Exception as e:
        _log.error(
            'Background load failed: %s -> %s (org %d): %s',
            source_type,
            target_type,
            ncbi_tax_id,
            e,
        )

    finally:
        with _lock:
            _pending.discard(key)
