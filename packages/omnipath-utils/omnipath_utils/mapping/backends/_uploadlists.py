"""UniProt ID Mapping (uploadlists) backend.

This backend uses the UniProt ID Mapping batch service for targeted
translations of specific source IDs. It does not use the standard
pypath/direct dispatch because its workflow is fundamentally different
(submit job, poll, collect).
"""

from __future__ import annotations

import time
import logging
from collections import defaultdict

import requests

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping._id_types import IdTypeRegistry
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)

IDMAPPING_RUN = 'https://rest.uniprot.org/idmapping/run'
IDMAPPING_STATUS = 'https://rest.uniprot.org/idmapping/status/%s'
IDMAPPING_STREAM = 'https://rest.uniprot.org/idmapping/stream/%s'

CHUNK_SIZE = 5000  # Max IDs per request
POLL_INTERVAL = 3  # Seconds between status checks
MAX_POLL = 60  # Max poll attempts


class UploadListsBackend(MappingBackend):
    """UniProt ID Mapping batch service for targeted translations.

    Overrides ``read()`` because this backend has its own column
    resolution logic (``uniprot_from``/``uniprot_to`` keys) and
    requires ``source_ids`` rather than downloading a full table.
    """

    name = 'uploadlists'
    yaml_key = 'uploadlists'

    def read(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        """Translate specific source IDs via the UniProt batch service."""

        reg = IdTypeRegistry.get()

        from_db = reg.backend_column(
            id_type, 'uniprot_from'
        ) or reg.backend_column(id_type, 'uploadlists')
        to_db = reg.backend_column(
            target_id_type, 'uniprot_to'
        ) or reg.backend_column(target_id_type, 'uploadlists')

        if not from_db or not to_db:
            _log.debug(
                'Uploadlists does not support %s -> %s',
                id_type,
                target_id_type,
            )
            return {}

        # For targeted translation: requires source_ids.
        source_ids = kwargs.get('source_ids')

        if not source_ids:
            _log.debug(
                'Uploadlists requires source_ids for targeted translation',
            )
            return {}

        _log.info(
            '%s: translating %d IDs (%s -> %s)',
            self.name,
            len(source_ids),
            id_type,
            target_id_type,
        )

        data = self._translate_batch(source_ids, from_db, to_db)

        _log.info(
            '%s: loaded %d entries for %s -> %s',
            self.name,
            len(data),
            id_type,
            target_id_type,
        )

        return data

    def _read_via_pypath(
        self,
        id_type,
        target_id_type,
        ncbi_tax_id,
        *,
        src_col,
        tgt_col,
        **kwargs,
    ):
        """Not used -- read() is overridden."""
        raise NotImplementedError  # pragma: no cover

    def _read_direct(
        self,
        id_type,
        target_id_type,
        ncbi_tax_id,
        *,
        src_col,
        tgt_col,
        **kwargs,
    ):
        """Not used -- read() is overridden."""
        raise NotImplementedError  # pragma: no cover

    def _translate_batch(
        self,
        ids: list[str],
        from_db: str,
        to_db: str,
    ) -> dict[str, set[str]]:
        """Submit IDs in chunks, collect results."""

        data: dict[str, set[str]] = defaultdict(set)

        for i in range(0, len(ids), CHUNK_SIZE):
            chunk = ids[i : i + CHUNK_SIZE]
            chunk_result = self._submit_and_collect(chunk, from_db, to_db)

            for src, targets in chunk_result.items():
                data[src].update(targets)

        return dict(data)

    def _submit_and_collect(
        self,
        ids: list[str],
        from_db: str,
        to_db: str,
    ) -> dict[str, set[str]]:
        """Submit one chunk and poll for results."""

        _log.info(
            'Submitting %d IDs to UniProt ID Mapping (%s -> %s)',
            len(ids),
            from_db,
            to_db,
        )

        # Submit job
        resp = requests.post(
            IDMAPPING_RUN,
            data={'from': from_db, 'to': to_db, 'ids': ','.join(ids)},
            timeout=30,
        )
        resp.raise_for_status()
        job_id = resp.json().get('jobId')

        if not job_id:
            _log.warning('No jobId in response: %s', resp.text)
            return {}

        # Poll for completion
        for _attempt in range(MAX_POLL):
            time.sleep(POLL_INTERVAL)

            status_resp = requests.get(
                IDMAPPING_STATUS % job_id,
                timeout=30,
            )
            status_data = status_resp.json()

            if (
                'results' in status_data
                or status_data.get('jobStatus') == 'FINISHED'
            ):
                break

            if status_resp.status_code == 303:  # Redirect to results
                break

            if status_data.get('jobStatus') == 'ERROR':
                _log.warning('ID Mapping job failed: %s', status_data)
                return {}
        else:
            _log.warning(
                'ID Mapping job timed out after %d polls',
                MAX_POLL,
            )
            return {}

        # Fetch results as TSV
        result_resp = requests.get(
            IDMAPPING_STREAM % job_id,
            params={'format': 'tsv'},
            timeout=120,
        )
        result_resp.raise_for_status()

        data: dict[str, set[str]] = defaultdict(set)
        lines = result_resp.text.strip().split('\n')

        for line in lines[1:]:  # skip header
            parts = line.split('\t')

            if len(parts) >= 2:
                src = parts[0].strip()
                tgt = parts[1].strip()

                if src and tgt:
                    data[src].add(tgt)

        _log.info(
            'ID Mapping: got %d mappings for %d source IDs',
            len(data),
            len(ids),
        )

        return dict(data)


register('uploadlists', UploadListsBackend)
