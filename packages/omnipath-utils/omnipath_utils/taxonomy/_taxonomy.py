#!/usr/bin/env python

#
# This file is part of the `omnipath_utils` Python module
#
# Copyright 2026
# Heidelberg University Hospital
#
# File author(s): Denes Turei (turei.denes@gmail.com)
#
# Distributed under the GPL-3.0-or-later license
# See the file `LICENSE` or read a copy at
# https://www.gnu.org/licenses/gpl-3.0.txt
#

"""Taxonomy manager -- organism name resolution and conversion."""

from __future__ import annotations

import logging

from pkg_infra.data import load

_log = logging.getLogger(__name__)

# Type alias for anything that can identify an organism
TaxonInput = int | str


class TaxonomyManager:
    """Resolve organism identifiers across naming systems.

    Loads organisms.yaml as seed data. Additional organisms can be
    loaded lazily from UniProt/Ensembl/KEGG APIs on demand.
    """

    _instance: TaxonomyManager | None = None

    def __init__(self) -> None:
        raw = load('organisms.yaml', module='omnipath_utils')
        self._by_taxid: dict[int, dict] = {int(k): v for k, v in raw.items()}
        self._build_indices()

    @classmethod
    def get(cls) -> TaxonomyManager:
        """Return the singleton instance, creating it if needed."""

        if cls._instance is None:
            cls._instance = cls()

        return cls._instance

    def _build_indices(self) -> None:
        """Build reverse lookup indices from all name forms."""

        self._to_taxid: dict[str, int] = {}

        for taxid, info in self._by_taxid.items():
            for field in (
                'common_name',
                'latin_name',
                'short_latin',
                'ensembl_name',
                'kegg_code',
                'mirbase_code',
                'oma_code',
                'uniprot_code',
                'dbptm_code',
            ):
                val = info.get(field, '')

                if val:
                    self._to_taxid[val.lower()] = taxid

            # Index synonyms
            for syn in info.get('synonyms', []):
                if syn:
                    self._to_taxid[syn.lower()] = taxid

            # Also index the taxid as string
            self._to_taxid[str(taxid)] = taxid

    def ensure_ncbi_tax_id(self, taxon: TaxonInput) -> int | None:
        """Convert any organism representation to NCBI Taxonomy ID.

        Args:
            taxon:
                Organism name or NCBI Taxonomy ID in any supported
                form: common name, latin name, Ensembl name, KEGG
                code, miRBase code, OMA code, UniProt code, dbPTM
                code, synonym, or integer/string taxonomy ID.

        Returns:
            NCBI Taxonomy ID, or None if the organism is not found.
        """

        if isinstance(taxon, int):
            return taxon if taxon in self._by_taxid else None

        if isinstance(taxon, str):
            taxon = taxon.strip()

            # Try as integer string
            try:
                tid = int(taxon)
                return tid if tid in self._by_taxid else None
            except ValueError:
                pass

            # Direct lookup (case-insensitive)
            result = self._to_taxid.get(taxon.lower())
            if result:
                return result

            # Parenthetical: "Caenorhabditis elegans (PRJNA13758)"
            if '(' in taxon:
                part0 = taxon.split('(', 1)[0].strip()
                part1 = taxon.split('(', 1)[1].rstrip(')').strip()
                return (
                    self._to_taxid.get(part0.lower()) or
                    self._to_taxid.get(part1.lower())
                )

            return None

        return None

    def _get_field(
        self,
        taxon: TaxonInput,
        field: str,
    ) -> str | None:
        """Get a specific field for an organism.

        Args:
            taxon:
                Organism identifier in any supported form.
            field:
                Field name from organisms.yaml.

        Returns:
            The field value, or None if not found.
        """

        taxid = self.ensure_ncbi_tax_id(taxon)

        if taxid is None:
            return None

        info = self._by_taxid.get(taxid)

        return info.get(field) if info else None

    def ensure_common_name(self, taxon: TaxonInput) -> str | None:
        """Common English name for an organism.

        Args:
            taxon:
                Organism identifier in any supported form.

        Returns:
            Common name string, or None if not found.
        """

        return self._get_field(taxon, 'common_name')

    def ensure_latin_name(self, taxon: TaxonInput) -> str | None:
        """Full latin (scientific) name for an organism.

        Args:
            taxon:
                Organism identifier in any supported form.

        Returns:
            Latin name string, or None if not found.
        """

        return self._get_field(taxon, 'latin_name')

    def ensure_ensembl_name(self, taxon: TaxonInput) -> str | None:
        """Ensembl organism name.

        Args:
            taxon:
                Organism identifier in any supported form.

        Returns:
            Ensembl name string, or None if not found.
        """

        return self._get_field(taxon, 'ensembl_name')

    def ensure_kegg_code(self, taxon: TaxonInput) -> str | None:
        """KEGG three-letter organism code.

        Args:
            taxon:
                Organism identifier in any supported form.

        Returns:
            KEGG code string, or None if not found.
        """

        return self._get_field(taxon, 'kegg_code')

    def ensure_mirbase_name(self, taxon: TaxonInput) -> str | None:
        """MiRBase three-letter organism code.

        Args:
            taxon:
                Organism identifier in any supported form.

        Returns:
            miRBase code string, or None if not found.
        """

        return self._get_field(taxon, 'mirbase_code')

    def ensure_oma_code(self, taxon: TaxonInput) -> str | None:
        """OMA five-letter organism code.

        Args:
            taxon:
                Organism identifier in any supported form.

        Returns:
            OMA code string, or None if not found.
        """

        return self._get_field(taxon, 'oma_code')

    def all_organisms(self) -> dict[int, dict]:
        """Return all known organisms.

        Returns:
            Dict mapping NCBI Taxonomy IDs to organism info dicts.
        """

        return dict(self._by_taxid)

    def load_from_ensembl(self) -> None:
        """Load organisms from Ensembl via pypath.inputs.ensembl."""
        try:
            from pypath.inputs.ensembl import ensembl_organisms
            count = 0
            for org in ensembl_organisms():
                taxid = int(org.taxon_id)
                if taxid not in self._by_taxid:
                    self._by_taxid[taxid] = {
                        'common_name': (org.common_name or '').lower(),
                        'latin_name': org.scientific_name or '',
                        'short_latin': '',
                        'ensembl_name': org.ensembl_name or '',
                        'kegg_code': '',
                        'mirbase_code': '',
                        'oma_code': '',
                        'uniprot_code': '',
                        'dbptm_code': '',
                        'synonyms': [],
                    }
                    count += 1
                else:
                    # Enrich
                    info = self._by_taxid[taxid]
                    if not info.get('ensembl_name') and org.ensembl_name:
                        info['ensembl_name'] = org.ensembl_name
                    if not info.get('latin_name') and org.scientific_name:
                        info['latin_name'] = org.scientific_name
                    if not info.get('common_name') and org.common_name:
                        info['common_name'] = org.common_name.lower()

            self._build_indices()
            _log.info(
                'Loaded %d new organisms from Ensembl (total: %d)',
                count,
                len(self._by_taxid),
            )
        except ImportError:
            _log.debug('pypath not available for Ensembl organisms')
        except Exception as e:
            _log.warning('Failed to load Ensembl organisms: %s', e)

    def load_from_uniprot(self) -> None:
        """Load organisms from UniProt speclist via pypath.inputs.uniprot."""
        try:
            from pypath.inputs.uniprot import uniprot_ncbi_taxids_2
            ncbi_data = uniprot_ncbi_taxids_2()
            count = 0
            for taxid, taxon in ncbi_data.items():
                taxid = int(taxid)
                if not taxid:
                    continue
                if taxid not in self._by_taxid:
                    self._by_taxid[taxid] = {
                        'common_name': (taxon.english or '').lower(),
                        'latin_name': taxon.latin or '',
                        'short_latin': '',
                        'ensembl_name': '',
                        'kegg_code': '',
                        'mirbase_code': '',
                        'oma_code': '',
                        'uniprot_code': '',
                        'dbptm_code': '',
                        'synonyms': [],
                    }
                    count += 1
                else:
                    info = self._by_taxid[taxid]
                    if not info.get('latin_name') and taxon.latin:
                        info['latin_name'] = taxon.latin
                    if not info.get('common_name') and taxon.english:
                        info['common_name'] = taxon.english.lower()

            self._build_indices()
            _log.info(
                'Loaded %d new organisms from UniProt (total: %d)',
                count,
                len(self._by_taxid),
            )
        except ImportError:
            _log.debug('pypath not available for UniProt organisms')
        except Exception as e:
            _log.warning('Failed to load UniProt organisms: %s', e)

    def load_from_kegg(self) -> None:
        """Load KEGG codes via pypath.inputs.kegg_organisms."""
        try:
            from pypath.inputs.kegg_organisms import kegg_organisms
            count = 0
            for org in kegg_organisms():
                # Try to match by latin name
                for taxid, info in self._by_taxid.items():
                    if (
                        info.get('latin_name', '').lower() == org.name.lower()
                        and not info.get('kegg_code')
                    ):
                        info['kegg_code'] = org.code
                        count += 1
                        break

            self._build_indices()
            _log.info('Enriched %d organisms with KEGG codes', count)
        except ImportError:
            _log.debug('pypath not available for KEGG organisms')
        except Exception as e:
            _log.warning('Failed to load KEGG organisms: %s', e)

    def load_from_mirbase(self) -> None:
        """Load miRBase codes via pypath.inputs.mirbase."""
        try:
            from pypath.inputs.mirbase import mirbase_organisms
            mirbase_to_ncbi = mirbase_organisms('mirbase', 'ncbi')
            count = 0
            for mirbase_code, ncbi_id in mirbase_to_ncbi.items():
                taxid = int(ncbi_id)
                if taxid in self._by_taxid:
                    if not self._by_taxid[taxid].get('mirbase_code'):
                        self._by_taxid[taxid]['mirbase_code'] = mirbase_code
                        count += 1

            self._build_indices()
            _log.info('Enriched %d organisms with miRBase codes', count)
        except ImportError:
            _log.debug('pypath not available for miRBase organisms')
        except Exception as e:
            _log.warning('Failed to load miRBase organisms: %s', e)

    def load_all(self) -> None:
        """Load organisms from all available sources."""
        self.load_from_ensembl()
        self.load_from_uniprot()
        self.load_from_kegg()
        self.load_from_mirbase()
        _log.info(
            'Taxonomy fully loaded: %d organisms',
            len(self._by_taxid),
        )

    def __contains__(self, taxon: TaxonInput) -> bool:
        return self.ensure_ncbi_tax_id(taxon) is not None

    def __len__(self) -> int:
        return len(self._by_taxid)
