"""Ensembl BioMart mapping backend.

Division-aware: vertebrates and the classic model organisms (fly, worm, yeast)
are served by the main Ensembl BioMart (``www.ensembl.org``, ``*_gene_ensembl``
datasets); the Ensembl Genomes divisions (plants, fungi, metazoa, protists) are
served by their own hosts and ``*_eg_gene`` datasets. Bacteria use per-genome
*collection* datasets whose names cannot be derived from the species name alone,
so they are not routed here (their Entrez anchoring comes from NCBI
gene_info/gene2accession/gene2ensembl instead).

Prefers ``pypath.inputs.biomart`` for the main host; falls back to (and, for
division organisms, always uses) direct HTTP against the appropriate martservice.
"""

from __future__ import annotations

import logging
import textwrap
from collections import defaultdict

from omnipath_utils.taxonomy import ensure_ensembl_name
from omnipath_utils.mapping.backends import register
from omnipath_utils.mapping.backends._base import MappingBackend

_log = logging.getLogger(__name__)

BIOMART_URL = 'https://www.ensembl.org/biomart/martservice'

# Ensembl Genomes divisions: martservice host + virtualSchemaName (mart name).
# Datasets follow the ``{ensembl_name}_eg_gene`` convention. Bacteria is
# deliberately absent — its collection-based dataset names are not derivable
# from the species name.
ENSEMBL_GENOMES_DIVISIONS = {
    'plants': ('https://plants.ensembl.org/biomart/martservice', 'plants_mart'),
    'fungi': ('https://fungi.ensembl.org/biomart/martservice', 'fungi_mart'),
    'metazoa': (
        'https://metazoa.ensembl.org/biomart/martservice', 'metazoa_mart',
    ),
    'protists': (
        'https://protists.ensembl.org/biomart/martservice', 'protists_mart',
    ),
}

# NCBI taxon -> Ensembl Genomes division, ONLY for organisms that are not on the
# main (vertebrate) Ensembl BioMart. Fly (7227), worm (6239) and yeast (4932)
# are on the main host as classic model organisms, so they are omitted here and
# use the vertebrate path. Extend as further non-vertebrate organisms enter scope.
ORGANISM_DIVISION = {
    3702: 'plants',      # Arabidopsis thaliana
    4577: 'plants',      # Zea mays
    39947: 'plants',     # Oryza sativa japonica
    5476: 'fungi',       # Candida albicans
    5061: 'fungi',       # Aspergillus niger
    36329: 'protists',   # Plasmodium falciparum 3D7
}


XML_TEMPLATE = textwrap.dedent("""\
    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query virtualSchemaName="{schema}" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
      <Dataset name="{dataset}" interface="default">
        {attributes}
      </Dataset>
    </Query>""")

ATTR_TEMPLATE = '<Attribute name="{name}" />'


class _Target:
    """Resolved BioMart endpoint for an organism."""

    __slots__ = ('host', 'schema', 'dataset', 'division')

    def __init__(self, host, schema, dataset, division):
        self.host = host
        self.schema = schema
        self.dataset = dataset
        self.division = division  # None for the main (vertebrate) host


class BioMartBackend(MappingBackend):
    """Fetch ID mappings from Ensembl / Ensembl Genomes BioMart."""

    name = 'biomart'
    yaml_key = 'ensembl'

    # ------------------------------------------------------------------
    # Endpoint resolution
    # ------------------------------------------------------------------

    @staticmethod
    def _target(ncbi_tax_id: int) -> _Target | None:
        """Resolve the host / virtualSchema / dataset for an organism.

        Returns ``None`` when the organism has no Ensembl name, or belongs to
        the bacteria division (collection datasets, not handled).
        """

        ensembl_name = ensure_ensembl_name(ncbi_tax_id)

        if not ensembl_name:
            _log.warning('No Ensembl name for organism %d', ncbi_tax_id)
            return None

        division = ORGANISM_DIVISION.get(ncbi_tax_id)

        if division is None:

            return _Target(
                host=BIOMART_URL,
                schema='default',
                dataset=f'{ensembl_name}_gene_ensembl',
                division=None,
            )

        if division not in ENSEMBL_GENOMES_DIVISIONS:
            _log.warning(
                'Organism %d division %r not routable via BioMart; skipping',
                ncbi_tax_id, division,
            )
            return None

        host, schema = ENSEMBL_GENOMES_DIVISIONS[division]

        return _Target(
            host=host,
            schema=schema,
            dataset=f'{ensembl_name}_eg_gene',
            division=division,
        )

    # ------------------------------------------------------------------
    # pypath.inputs path (main/vertebrate host only)
    # ------------------------------------------------------------------

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
        """Read via ``pypath.inputs.biomart`` — main host only.

        pypath's ``biomart_query`` hardwires the vertebrate martservice URL, so
        for Ensembl Genomes divisions we raise ``ImportError`` to make the base
        class fall through to the division-aware :meth:`_read_direct`.
        """

        target = self._target(ncbi_tax_id)

        if target is None:
            return {}

        if target.division is not None:
            # pypath cannot reach the genomes host — force the direct path.
            raise ImportError('division organism: use direct BioMart')

        from pypath.inputs.biomart import biomart_query

        attrs = [src_col, tgt_col] if src_col != tgt_col else [src_col]

        _log.info(
            'BioMart query (pypath): attrs=%s, dataset=%s',
            attrs, target.dataset,
        )

        data: dict[str, set[str]] = defaultdict(set)

        for rec in biomart_query(attrs=attrs, dataset=target.dataset):
            id_a = getattr(rec, src_col, '')
            id_b = getattr(rec, tgt_col, '')

            if id_a and id_b:
                data[id_a].add(id_b)

        return dict(data)

    # ------------------------------------------------------------------
    # Direct HTTP (main host + Ensembl Genomes divisions)
    # ------------------------------------------------------------------

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
        """Query the appropriate martservice directly (host keyed on division)."""

        import requests

        target = self._target(ncbi_tax_id)

        if target is None:
            return {}

        attrs = [src_col, tgt_col] if src_col != tgt_col else [src_col]
        attr_xml = '\n        '.join(
            ATTR_TEMPLATE.format(name=a) for a in attrs
        )

        xml = XML_TEMPLATE.format(
            schema=target.schema, dataset=target.dataset, attributes=attr_xml,
        )
        xml = xml.replace('\n', '').replace('  ', '')  # compact

        _log.info(
            'Querying BioMart (direct HTTP): %s -> %s (host %s, dataset %s)',
            id_type, target_id_type, target.host, target.dataset,
        )

        resp = requests.get(target.host, params={'query': xml}, timeout=120)
        resp.raise_for_status()

        data: dict[str, set[str]] = defaultdict(set)
        lines = resp.text.strip().split('\n')

        for line in lines[1:]:  # skip header
            parts = line.split('\t')

            if len(parts) >= 2:
                src = parts[0].strip()
                tgt = parts[1].strip()

                if src and tgt:
                    data[src].add(tgt)

        return dict(data)


register('biomart', BioMartBackend)
