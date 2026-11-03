"""PubChem structure-bearing mapping backend (Milestone I).

Streams ``pypath.inputs.pubchem.pubchem_mapping(target, source='cid')``, which
yields ``(cid, value)`` rows from the compact PubChem FTP ``Extras`` tables
(``Compound/Extras/CID-<TYPE>.gz``) -- CID -> standard InChIKey / InChI /
SMILES. Because the rows are streamed, ``limit`` (``--pubchem-max-records``)
is honoured by stopping early via :func:`itertools.islice`. No RDKit.
"""

from __future__ import annotations

import logging
from itertools import islice

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)

_CID_COLUMN = 'cid'


class PubchemBackend(MappingBackend):
    name = 'pubchem'
    yaml_key = 'pubchem'

    def _read_via_pypath(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        *,
        src_col: str,
        tgt_col: str,
        limit: int | None = None,
        **kwargs: object,
    ) -> dict[str, set[str]]:
        from pypath.inputs.pubchem import pubchem_mapping

        # One side is always the Compound CID; stream CID -> <structure> and
        # orient the result accordingly.
        if src_col == _CID_COLUMN:
            target, invert = tgt_col, False
        elif tgt_col == _CID_COLUMN:
            target, invert = src_col, True
        else:
            _log.debug(
                'pubchem backend only bridges via CID (%s -> %s)',
                src_col,
                tgt_col,
            )
            return {}

        rows = pubchem_mapping(target, source=_CID_COLUMN)
        if limit is not None:
            rows = islice(rows, limit)

        result: dict[str, set[str]] = {}
        for cid, value in rows:
            if invert:
                result.setdefault(value, set()).add(cid)
            else:
                result.setdefault(cid, set()).add(value)

        return result

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
        raise ImportError('pubchem requires pypath')


register('pubchem', PubchemBackend)
