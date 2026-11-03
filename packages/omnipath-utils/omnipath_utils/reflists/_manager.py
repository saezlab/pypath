"""Reference list manager — provides sets of identifiers."""

from __future__ import annotations

import logging

from omnipath_utils._constants import DEFAULT_ORGANISM

_log = logging.getLogger(__name__)


class ReferenceListManager:
    """Manages reference lists (complete sets of IDs for an organism).

    Reference lists are derived from the mapping infrastructure.
    For example, "all human SwissProt IDs" is obtained by querying
    UniProt for all reviewed proteins of organism 9606.
    """

    _instance: ReferenceListManager | None = None

    _loader_methods: dict[str, str] = {
        'swissprot': '_load_swissprot',
        'trembl': '_load_trembl',
        'uniprot': '_load_uniprot_all',
    }

    def __init__(self):
        self._cache: dict[tuple, set[str]] = {}

    @classmethod
    def get(cls) -> ReferenceListManager:
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    def get_reflist(
        self,
        list_name: str,
        ncbi_tax_id: int = DEFAULT_ORGANISM,
    ) -> set[str]:
        """Get a reference list by name.

        Supported list names:
        - 'swissprot': all reviewed UniProt ACs
        - 'trembl': all unreviewed UniProt ACs
        - 'uniprot': all UniProt ACs (SwissProt + TrEMBL)

        Results are cached after first load.
        """
        key = (list_name, ncbi_tax_id)

        if key not in self._cache:
            method_name = self._loader_methods.get(list_name)
            if not method_name:
                _log.warning('Unknown reference list: %s', list_name)
                return set()
            loader = getattr(self, method_name)
            self._cache[key] = loader(ncbi_tax_id)

        return self._cache[key]

    def _load_swissprot(self, ncbi_tax_id: int) -> set[str]:
        """Load all reviewed UniProt ACs for an organism."""
        _log.info(
            'Loading SwissProt reference list for organism %d', ncbi_tax_id
        )

        try:
            from pypath.inputs.uniprot import UniprotQuery

            q = UniprotQuery(
                organism=ncbi_tax_id, fields='accession', reviewed=True
            )
            data = q.perform()
            result = set(data.keys()) if data else set()
        except ImportError:
            # Fallback: use direct HTTP
            result = self._load_uniprot_direct(ncbi_tax_id, reviewed=True)

        _log.info(
            'SwissProt list: %d IDs for organism %d', len(result), ncbi_tax_id
        )
        return result

    def _load_trembl(self, ncbi_tax_id: int) -> set[str]:
        """Load all unreviewed UniProt ACs for an organism."""
        _log.info('Loading TrEMBL reference list for organism %d', ncbi_tax_id)

        try:
            from pypath.inputs.uniprot import UniprotQuery

            q = UniprotQuery(
                organism=ncbi_tax_id, fields='accession', reviewed=False
            )
            data = q.perform()
            result = set(data.keys()) if data else set()
        except ImportError:
            result = self._load_uniprot_direct(ncbi_tax_id, reviewed=False)

        _log.info(
            'TrEMBL list: %d IDs for organism %d', len(result), ncbi_tax_id
        )
        return result

    def _load_uniprot_all(self, ncbi_tax_id: int) -> set[str]:
        """Load all UniProt ACs (SwissProt + TrEMBL)."""
        return self._load_swissprot(ncbi_tax_id) | self._load_trembl(
            ncbi_tax_id
        )

    def load_swissprot_global(self) -> set[str]:
        """Load the complete reviewed (SwissProt) AC set, all organisms.

        Reviewed status is organism-agnostic; one ``reviewed:true`` query (no
        organism filter) returns the whole reviewed set (~570k ACs). Cached
        under the sentinel key ``('swissprot', 0)``.
        """
        key = ('swissprot', 0)
        if key not in self._cache:
            try:
                from pypath.inputs.uniprot import UniprotQuery

                q = UniprotQuery(fields='accession', reviewed=True)
                data = q.perform()
                result = set(data.keys()) if data else set()
            except Exception:
                result = self._load_uniprot_global_direct(reviewed=True)
            self._cache[key] = result
        return self._cache[key]

    def _load_uniprot_global_direct(self, reviewed: bool = True) -> set[str]:
        """Direct HTTP fallback: all (un)reviewed ACs, no organism.

        Uses the *compressed* stream (~5 MB gzip) which is far more robust than
        the plain stream for the full ~570k-row pull (the latter routinely drops
        mid-response). Retries the whole request a few times on connection drop.
        """
        import gzip
        import io

        import requests

        url = 'https://rest.uniprot.org/uniprotkb/stream'
        params = {
            'query': f'reviewed:{str(reviewed).lower()}',
            'fields': 'accession',
            'format': 'list',
            'compressed': 'true',
        }
        last_err: Exception | None = None
        for _attempt in range(4):
            try:
                resp = requests.get(
                    url, params=params, timeout=900,
                    headers={'Accept-Encoding': 'gzip'},
                )
                resp.raise_for_status()
                raw = gzip.GzipFile(fileobj=io.BytesIO(resp.content)).read()
                return {
                    line.strip()
                    for line in raw.decode().splitlines()
                    if line.strip()
                }
            except Exception as err:  # noqa: BLE001 - retry transient drops
                last_err = err
                _log.warning('global reviewed-list fetch failed: %s', err)
        raise RuntimeError(
            f'could not fetch global reviewed list: {last_err}'
        )

    def _load_uniprot_direct(
        self, ncbi_tax_id: int, reviewed: bool | None = None
    ) -> set[str]:
        """Direct HTTP fallback for loading UniProt AC lists."""
        import requests

        url = 'https://rest.uniprot.org/uniprotkb/stream'
        params = {
            'query': f'(organism_id:{ncbi_tax_id})',
            'fields': 'accession',
            'format': 'tsv',
        }
        if reviewed is not None:
            params['query'] += f' AND (reviewed:{str(reviewed).lower()})'

        resp = requests.get(url, params=params, timeout=300)
        resp.raise_for_status()

        lines = resp.text.strip().split('\n')
        return {line.strip() for line in lines[1:] if line.strip()}

    def is_swissprot(
        self, uniprot_ac: str, ncbi_tax_id: int = DEFAULT_ORGANISM
    ) -> bool:
        """Check if a UniProt AC is in the SwissProt (reviewed) set."""
        return uniprot_ac in self.get_reflist('swissprot', ncbi_tax_id)

    def is_trembl(
        self, uniprot_ac: str, ncbi_tax_id: int = DEFAULT_ORGANISM
    ) -> bool:
        """Check if a UniProt AC is in the TrEMBL (unreviewed) set."""
        return uniprot_ac in self.get_reflist('trembl', ncbi_tax_id)
