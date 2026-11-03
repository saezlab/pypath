"""Orthology manager -- loads and queries ortholog tables."""

from __future__ import annotations

import logging
import threading
from collections import defaultdict

_log = logging.getLogger(__name__)


class OrthologyTable:
    """Cached orthology translation dict for one organism pair + resource."""

    def __init__(self, data, source, target, id_type, resource, metadata=None):
        self.data = data  # {source_id: set(target_ids)}
        self.source = source
        self.target = target
        self.id_type = id_type
        self.resource = resource
        self.metadata = (
            metadata or {}
        )  # {source_id: {target_id: {rel_type, score, ...}}}

    def __getitem__(self, key):
        return self.data.get(key, set())

    def __len__(self):
        return len(self.data)


class OrthologyManager:
    """Central orthology translation manager."""

    _instance = None
    _lock = threading.Lock()

    def __init__(self):
        self._tables = {}  # (source, target, id_type, resource) -> OrthologyTable

    @classmethod
    def get(cls) -> OrthologyManager:
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance

    def translate(
        self,
        identifiers,
        source=9606,
        target=10090,
        id_type='genesymbol',
        only_swissprot=True,
        resource=None,
        min_sources=1,
        rel_type=None,
        raw=False,
    ):
        table = self._get_table(source, target, id_type, resource, min_sources)

        if not table:
            return {name: set() for name in identifiers}

        result = {}
        for name in identifiers:
            orthologs = table[name]

            # Filter by rel_type if specified
            if rel_type and table.metadata:
                filtered = set()
                for orth in orthologs:
                    meta = table.metadata.get(name, {}).get(orth, {})
                    if meta.get('rel_type') in rel_type:
                        filtered.add(orth)
                orthologs = filtered if filtered else orthologs

            result[name] = orthologs

        return result

    def get_table(self, source, target, id_type, resource=None, min_sources=1):
        table = self._get_table(source, target, id_type, resource, min_sources)
        return dict(table.data) if table else {}

    def _get_table(self, source, target, id_type, resource=None, min_sources=1):
        resources = (
            [resource]
            if resource
            else ['hcop', 'ensembl', 'oma', 'orthodb', 'alliance', 'homologene']
        )

        for res in resources:
            key = (source, target, id_type, res)
            if key in self._tables:
                return self._tables[key]

            table = self._load_table(source, target, id_type, res, min_sources)
            if table:
                self._tables[key] = table
                return table

        return None

    def _load_table(self, source, target, id_type, resource, min_sources=1):
        _log.info(
            'Loading orthology: %d -> %d, id_type=%s, resource=%s',
            source,
            target,
            id_type,
            resource,
        )

        try:
            if resource == 'hcop':
                return self._load_hcop(source, target, id_type, min_sources)
            elif resource == 'oma':
                return self._load_oma(source, target, id_type)
            elif resource == 'ensembl':
                return self._load_ensembl(source, target, id_type)
            elif resource == 'homologene':
                return self._load_homologene(source, target, id_type)
            elif resource == 'orthodb':
                return self._load_orthodb(source, target, id_type)
            elif resource == 'alliance':
                return self._load_alliance(source, target, id_type)
        except ImportError:
            _log.debug(
                'pypath not available for orthology resource %s', resource
            )
        except Exception as e:
            _log.warning('Failed to load orthology from %s: %s', resource, e)

        return None

    def _load_hcop(self, source, target, id_type, min_sources):
        from pypath.inputs.hcop import hcop_orthologs

        # HCOP is human-centric: source must be 9606
        if source != 9606:
            _log.debug('HCOP only supports human as source')
            return None

        records = hcop_orthologs(target, min_sources=min_sources)

        if not records:
            return None

        _ID_FIELD_MAP = {
            'genesymbol': ('human_symbol', 'ortholog_symbol'),
            'entrez': ('human_entrez', 'ortholog_entrez'),
            'ensembl': ('human_ensembl', 'ortholog_ensembl'),
            'ensg': ('human_ensembl', 'ortholog_ensembl'),
        }

        fields = _ID_FIELD_MAP.get(id_type, ('human_symbol', 'ortholog_symbol'))

        data = defaultdict(set)
        metadata = defaultdict(lambda: defaultdict(dict))

        for rec in records:
            src_id = getattr(rec, fields[0])
            tgt_id = getattr(rec, fields[1])
            if src_id and tgt_id:
                data[src_id].add(tgt_id)
                metadata[src_id][tgt_id] = {
                    'support': rec.support,
                    'n_sources': rec.n_sources,
                }

        return OrthologyTable(
            data=dict(data),
            source=source,
            target=target,
            id_type=id_type,
            resource='hcop',
            metadata=dict(metadata),
        )

    def _load_oma(self, source, target, id_type):
        from pypath.inputs.oma import oma_orthologs

        records = oma_orthologs(
            source, target, id_type=id_type, return_df=False
        )

        if not records:
            return None

        data = defaultdict(set)
        metadata = defaultdict(lambda: defaultdict(dict))

        for rec in records:
            src_id = str(rec.id_a) if hasattr(rec, 'id_a') else str(rec[0])
            tgt_id = str(rec.id_b) if hasattr(rec, 'id_b') else str(rec[1])

            if src_id and tgt_id:
                data[src_id].add(tgt_id)

                meta = {}
                if hasattr(rec, 'rel_type'):
                    meta['rel_type'] = rec.rel_type
                if hasattr(rec, 'score'):
                    meta['score'] = rec.score
                if meta:
                    metadata[src_id][tgt_id] = meta

        return OrthologyTable(
            data=dict(data),
            source=source,
            target=target,
            id_type=id_type,
            resource='oma',
            metadata=dict(metadata),
        )

    def _load_ensembl(self, source, target, id_type):
        from pypath.inputs.biomart import biomart_homology

        records = biomart_homology(source, target)

        if not records:
            return None

        # BioMart returns Ensembl gene/peptide IDs
        data = defaultdict(set)
        metadata = defaultdict(lambda: defaultdict(dict))

        for rec in records:
            _src_id = (
                rec.ensembl_peptide_id
                if hasattr(rec, 'ensembl_peptide_id')
                else ''
            )
            _tgt_id = ''
            if hasattr(rec, '_fields'):
                for f in rec._fields:
                    if 'homolog_ensembl_peptide' in f:
                        _tgt_id = getattr(rec, f, '')
                        break

            # Also get gene-level
            src_gene = (
                rec.ensembl_gene_id if hasattr(rec, 'ensembl_gene_id') else ''
            )

            orth_type = ''
            confidence = ''
            for field in getattr(rec, '_fields', ()):
                if 'orthology_type' in field:
                    orth_type = getattr(rec, field, '')
                if 'orthology_confidence' in field:
                    confidence = getattr(rec, field, '')

            # Try to get gene symbol from the record
            tgt_symbol = ''
            for field in getattr(rec, '_fields', ()):
                if 'associated_gene_name' in field:
                    tgt_symbol = getattr(rec, field, '')

            src_symbol = getattr(rec, 'external_gene_name', '')

            if id_type == 'genesymbol' and src_symbol and tgt_symbol:
                data[src_symbol].add(tgt_symbol)
                metadata[src_symbol][tgt_symbol] = {
                    'orth_type': orth_type,
                    'confidence': confidence,
                }
            elif id_type in ('ensg', 'ensembl') and src_gene:
                tgt_gene = ''
                for field in getattr(rec, '_fields', ()):
                    if 'homolog_ensembl_gene' in field:
                        tgt_gene = getattr(rec, field, '')
                if tgt_gene:
                    data[src_gene].add(tgt_gene)
                    metadata[src_gene][tgt_gene] = {
                        'orth_type': orth_type,
                        'confidence': confidence,
                    }

        return (
            OrthologyTable(
                data=dict(data),
                source=source,
                target=target,
                id_type=id_type,
                resource='ensembl',
                metadata=dict(metadata),
            )
            if data
            else None
        )

    def _load_homologene(self, source, target, id_type):
        from pypath.inputs.homologene import homologene_dict

        try:
            data = homologene_dict(source, target, id_type)
        except Exception as e:
            _log.debug('HomoloGene failed: %s', e)
            return None

        if not data:
            return None

        # homologene_dict returns {source_id: set(target_ids)}
        return OrthologyTable(
            data=data,
            source=source,
            target=target,
            id_type=id_type,
            resource='homologene',
        )

    def _load_orthodb(self, source, target, id_type):
        from pypath.inputs.orthodb import orthodb_orthologs

        data = orthodb_orthologs(
            source=source,
            target=target,
            id_type=id_type,
        )

        if not data:
            return None

        return OrthologyTable(
            data=data,
            source=source,
            target=target,
            id_type=id_type,
            resource='orthodb',
        )

    def _load_alliance(self, source, target, id_type):
        from pypath.inputs.alliance import alliance_dict

        data = alliance_dict(
            source=source,
            target=target,
            id_type=id_type,
        )

        if not data:
            return None

        return OrthologyTable(
            data=data,
            source=source,
            target=target,
            id_type=id_type,
            resource='alliance',
        )
