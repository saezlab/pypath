"""UniChem mapping backend for small molecule identifiers."""

from __future__ import annotations

import logging

from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)

# Map our canonical names to UniChem source labels used by pypath.
# UniChem labels are case-sensitive; pypath uses the label field
# returned by the UniChem REST API.
_UNICHEM_LABELS: dict[str, str] = {
    'chembl': 'ChEMBL',
    'chebi': 'ChEBI',
    'drugbank': 'DrugBank',
    'pubchem': 'PubChem',
    'hmdb': 'HMDB',
    'kegg': 'KEGG Ligand',
    'lipidmaps': 'LIPID MAPS®',
    'surechembl': 'SureChEMBL',
    'swisslipids': 'SwissLipids',
    'brenda': 'Brenda',
    'rhea': 'Rhea',
    'pdb': 'PDBe',
}


def _to_unichem_label(name: str) -> str:
    """Convert our canonical name to a UniChem label for pypath."""
    return _UNICHEM_LABELS.get(name, name)


class UniChemBackend(MappingBackend):
    name = 'unichem'
    yaml_key = 'unichem'

    def read(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        """UniChem mappings are not organism-specific."""
        try:
            return self._read_via_pypath(id_type, target_id_type, ncbi_tax_id)
        except ImportError:
            _log.debug('pypath not available for unichem backend')
            return {}

    def _read_via_pypath(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        from pypath.inputs.unichem import unichem_mapping

        label_a = _to_unichem_label(id_type)
        label_b = _to_unichem_label(target_id_type)

        _log.info('UniChem: %s (%s) -> %s (%s)', id_type, label_a, target_id_type, label_b)

        try:
            data = unichem_mapping(label_a, label_b)
        except (ValueError, KeyError) as e:
            _log.debug('UniChem: %s', e)
            return {}

        if not data:
            return {}

        # unichem_mapping returns {id: set(ids)} already
        return {k: (v if isinstance(v, set) else {v}) for k, v in data.items()}

    def _read_direct(self, id_type, target_id_type, ncbi_tax_id, **kwargs):
        raise ImportError('unichem requires pypath')


register('unichem', UniChemBackend)
