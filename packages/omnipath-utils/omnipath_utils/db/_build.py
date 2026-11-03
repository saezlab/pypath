"""Database build orchestrator."""

from __future__ import annotations

import os
import time
import logging

from sqlalchemy import text
from sqlalchemy.orm import Session

from omnipath_utils.db._schema import (
    Base, IdType, Backend, Organism, BuildInfo, NamespaceExemption, LipidName,
)
from omnipath_utils.db._connection import SCHEMA, get_engine, ensure_schema
from omnipath_utils.mapping._id_types import IdTypeRegistry, normalize_identifier
from omnipath_utils.taxonomy._taxonomy import TaxonomyManager

_log = logging.getLogger(__name__)


class DatabaseBuilder:
    """Orchestrates building the omnipath_utils database."""

    def __init__(
        self,
        db_url: str | None = None,
        max_records: int | None = None,
        pubchem_max_records: int | None = None,
    ):
        self._db_url = db_url
        self._max_records = max_records
        self._pubchem_max_records = pubchem_max_records
        self.engine = get_engine(db_url)
        ensure_schema(self.engine)
        if max_records is not None:
            _log.warning(
                'max_records=%d: mapping tables are CAPPED for testing; '
                'this is NOT a complete build',
                max_records,
            )
        if pubchem_max_records is not None:
            _log.warning(
                'pubchem_max_records=%d: only the PubChem table is capped; '
                'all other resources load in full',
                pubchem_max_records,
            )

    def _effective_limit(self, backend_name: str) -> int | None:
        """Row cap for a backend: PubChem honours its own cap when set,
        every other backend uses the global ``max_records`` (if any)."""
        if backend_name == 'pubchem' and self._pubchem_max_records is not None:
            return self._pubchem_max_records
        return self._max_records

    def _run_mappings_parallel(self, pairs, organisms):
        """Run independent ``populate_mapping`` calls across a thread pool.

        Each ``(source, target, organism, backend)`` combination writes a
        disjoint slice of ``id_mapping`` and opens its own connection, so the
        jobs are safe to run concurrently. Worker count is controlled by
        ``OMNIPATH_BUILD_MAPPING_WORKERS`` (default 8); see
        ``omnipath-build/docs/build-tuning.md``.
        """
        from concurrent.futures import ThreadPoolExecutor

        jobs = [
            (src, tgt, org, backend)
            for src, tgt, backend in pairs
            for org in organisms
        ]
        if not jobs:
            return

        workers = int(os.environ.get('OMNIPATH_BUILD_MAPPING_WORKERS', '8'))
        workers = max(1, min(workers, len(jobs)))
        _log.info(
            'Building %d mapping tables with %d parallel workers',
            len(jobs),
            workers,
        )

        def _one(job):
            src, tgt, org, backend = job
            try:
                self.populate_mapping(src, tgt, org, backend)
            except Exception as e:
                _log.error(
                    'Failed: %s -> %s (org %d): %s', src, tgt, org, e,
                )

        with ThreadPoolExecutor(max_workers=workers) as ex:
            list(ex.map(_one, jobs))

    def create_tables(self):
        """Create all tables."""
        Base.metadata.create_all(self.engine)
        _log.info('All tables created')

    def populate_id_types(self):
        """Populate id_type table from id_types.yaml."""
        registry = IdTypeRegistry.get()

        with Session(self.engine) as session:
            for name in registry.all_names():
                info = registry.info(name)
                # Also write the CURIE + URL/id-pattern metadata, and update it
                # on an existing row so a cumulative build picks up registry
                # changes (the table is otherwise insert-only for id types).
                fields = {
                    'label': info.get('label'),
                    'entity_type': info.get('entity_type'),
                    'curie_prefix': info.get('curie_prefix'),
                    'url_pattern': registry.url_pattern(name),
                    'id_pattern': registry.id_pattern(name),
                }
                existing = session.query(IdType).filter_by(name=name).first()
                if existing:
                    for column, value in fields.items():
                        setattr(existing, column, value)
                else:
                    session.add(IdType(name=name, **fields))
            session.commit()

        _log.info('Populated %d ID types', len(registry))

    def populate_backends(self):
        """Populate backend table."""
        backends = [
            'uniprot',
            'uploadlists',
            'biomart',
            'pro',
            'unichem',
            'ramp',
            'hmdb',
            'array',
            'mirbase',
            'file',
            'uniprot_ftp',
            'gene2ensembl',
            'gene_info',
            'gene2accession',
            'gene_history',
            'ensembl_history',
            'kegg_gene',
            'metanetx',
            'bigg',
            'refmet',
            'reactome',
            'recon3d',
            # Structure-bearing backends (inputs_v2 adapter + PubChem)
            'chebi',
            'chembl',
            'lipidmaps',
            'swisslipids',
            'pubchem',
            'kegg_compound',
        ]

        with Session(self.engine) as session:
            for name in backends:
                existing = session.query(Backend).filter_by(name=name).first()
                if not existing:
                    session.add(Backend(name=name))
            session.commit()

        _log.info('Populated %d backends', len(backends))

    #: Namespaces the build declares (id_types.yaml) but deliberately does
    #: not load, reviewed 2026-09-09 (spec 011 R7/T093 -- "the exemption is
    #: data, not silence").
    _NAMESPACE_EXEMPTIONS: dict[str, str] = {
        'kegg_glycan': (
            'id_types.yaml declares this reachable via ramp, but RaMP\'s '
            'auto-discovered compound ID types (_populate_ramp) do not '
            'currently include a kegg_glycan entry.'
        ),
        'lipidbank': (
            'id_types.yaml declares this reachable via ramp, but RaMP\'s '
            'auto-discovered compound ID types do not currently include a '
            'lipidbank entry.'
        ),
        'plantfa': (
            'id_types.yaml declares this reachable via ramp, but RaMP\'s '
            'auto-discovered compound ID types do not currently include a '
            'plantfa entry.'
        ),
        'polymer': (
            'id_types.yaml declares this reachable via ramp, but RaMP\'s '
            'auto-discovered compound ID types do not currently include a '
            'polymer entry.'
        ),
        'phenol_explorer': (
            'No pypath loader exists; id_types.yaml declares it reachable '
            'via hmdb, but HMDB\'s raw metabolite records carry no '
            'phenol_explorer cross-reference field.'
        ),
        'wikipedia': (
            'No pypath loader exists; id_types.yaml declares it reachable '
            'via hmdb, but HMDB\'s raw metabolite records carry no '
            'wikipedia cross-reference field.'
        ),
        'biocyc': (
            'No pypath loader exists; id_types.yaml declares it reachable '
            'via hmdb, but HMDB\'s raw metabolite records carry no biocyc '
            'cross-reference field.'
        ),
        'fbonto': (
            'No pypath loader exists; id_types.yaml declares it reachable '
            'via hmdb, but HMDB\'s raw metabolite records carry no fbonto '
            'cross-reference field.'
        ),
        'ehmn': (
            'No pypath loader exists for chalmers_gem (the Human-GEM-style '
            'genome-scale metabolic model id_types.yaml declares this '
            'reachable through) at all -- no inputs_v2 or legacy inputs '
            'module.'
        ),
        'hepatonet': (
            'No pypath loader exists for chalmers_gem (see ehmn -- same '
            'root cause, same missing backend).'
        ),
        'hmr2': (
            'No pypath loader exists for chalmers_gem (see ehmn -- same '
            'root cause, same missing backend).'
        ),
        'metabolicatlas': (
            'No pypath loader exists for chalmers_gem (see ehmn -- same '
            'root cause, same missing backend).'
        ),
        'metabolicatlas_comp': (
            'No pypath loader exists for chalmers_gem (see ehmn -- same '
            'root cause, same missing backend).'
        ),
        'knapsack': (
            'No pypath loader exists (no inputs_v2 or legacy inputs module); '
            'id_types.yaml declares it reachable via hmdb, but HMDB\'s raw '
            'metabolite records carry no knapsack cross-reference field.'
        ),
        'chemspider': (
            'No pypath loader exists; id_types.yaml declares it reachable '
            'via hmdb, but HMDB\'s raw metabolite records carry no '
            'chemspider cross-reference field.'
        ),
        'metlin': (
            'No pypath loader exists; id_types.yaml declares it reachable '
            'via hmdb, but HMDB\'s raw metabolite records carry no metlin '
            'cross-reference field.'
        ),
        'vmh': (
            'No dedicated pypath loader; Recon3D (hosted on VMH\'s own '
            'site) does not carry a distinct VMH accession field in its '
            'own schema, only bigg_metabolite_id/hmdb/chebi/kegg_compound/'
            'metanetx.'
        ),
        'wikidata': (
            'No pypath loader exists; id_types.yaml declares it reachable '
            'via ramp, but RaMP\'s auto-discovered compound ID types do '
            'not currently include a wikidata entry.'
        ),
        'recon3d': (
            'pypath.inputs_v2.recon3d has a working schema (name/'
            'bigg_metabolite_id/hmdb/chebi/kegg_compound/metanetx), but '
            'its custom JSON raw_parser (data_type dispatch over the BiGG '
            'Models JSON export) returns zero rows through the generic '
            'inputs_v2 adapter raw_rows() this build uses for every other '
            'CSV-based source -- needs a dedicated reader, not a data gap.'
        ),
        'rhea': (
            'UniChem\'s Rhea source and RaMP\'s rhea-comp type both return '
            'the reaction participant\'s ChEBI id as the "Rhea id" -- Rhea '
            'mints no compound identifier of its own, participants are '
            'cited by ChEBI id directly. Excluded from both auto-discovery '
            'loops (T091) rather than stored as a wrong value.'
        ),
        'pdbe': (
            'UniChem\'s PDBe source returns a display label plus conformer '
            'type (e.g. "GR2 - Ideal conformer") as the "PDBe id", not the '
            'ligand code. Excluded from UniChem auto-discovery (T091) '
            'rather than stored as a wrong value; PDBe\'s own ligand '
            'dictionary would be the correct source for a real loader.'
        ),
    }

    def populate_namespace_exemptions(self):
        """Record every reviewed exemption (spec 011 T093), upserting the
        reason so a re-run picks up a changed one.
        """
        with Session(self.engine) as session:
            for namespace, reason in self._NAMESPACE_EXEMPTIONS.items():
                existing = (
                    session.query(NamespaceExemption)
                    .filter_by(namespace=namespace)
                    .first()
                )
                if existing:
                    existing.reason = reason
                else:
                    session.add(
                        NamespaceExemption(namespace=namespace, reason=reason)
                    )
            session.commit()
        _log.info(
            'Recorded %d namespace exemptions', len(self._NAMESPACE_EXEMPTIONS)
        )

    def _clear_wrong_content_rows(self):
        """Delete the id_mapping rows T091 found were storing the wrong
        thing (rhea holding ChEBI values, pdbe holding descriptive text) --
        the code fix (``_UNICHEM_WRONG_CONTENT``) stops new wrong rows from
        being written; this removes what a prior run already wrote.
        """
        with Session(self.engine) as session:
            for namespace in ('rhea', 'pdbe'):
                id_type = (
                    session.query(IdType).filter_by(name=namespace).first()
                )
                if not id_type:
                    continue
                deleted = session.execute(
                    text(
                        f'DELETE FROM {SCHEMA}.id_mapping'
                        ' WHERE source_type_id = :t OR target_type_id = :t'
                    ),
                    {'t': id_type.id},
                ).rowcount
                _log.info(
                    'Cleared %d wrong-content id_mapping rows for %s',
                    deleted, namespace,
                )
            session.commit()

    def populate_organisms(self):
        """Populate organism table from all available sources."""
        tm = TaxonomyManager.get()
        tm.load_all()
        orgs = tm.all_organisms()

        # Fields that exist in the Organism table
        _db_fields = {
            'common_name', 'latin_name', 'short_latin', 'ensembl_name',
            'kegg_code', 'mirbase_code', 'oma_code', 'uniprot_code',
            'dbptm_code',
        }

        with Session(self.engine) as session:
            for taxid, info in orgs.items():
                existing = (
                    session.query(Organism).filter_by(ncbi_tax_id=taxid).first()
                )
                if not existing:
                    session.add(
                        Organism(
                            ncbi_tax_id=taxid,
                            **{
                                k: v for k, v in info.items()
                                if v and k in _db_fields
                            },
                        )
                    )
            session.commit()

        _log.info('Populated %d organisms', len(orgs))

    def populate_mapping(
        self,
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int,
        backend_name: str,
    ):
        """Populate id_mapping table for one (source, target, organism, backend) combination.

        Uses the mapping infrastructure to load data, then bulk inserts via COPY.
        """
        from omnipath_utils.mapping._reader import MapReader

        _log.info(
            'Building mapping: %s -> %s, organism %d, backend %s',
            id_type,
            target_id_type,
            ncbi_tax_id,
            backend_name,
        )

        start = time.time()

        # Get type and backend IDs
        with Session(self.engine) as session:
            src_type = session.query(IdType).filter_by(name=id_type).first()
            tgt_type = (
                session.query(IdType).filter_by(name=target_id_type).first()
            )
            backend = (
                session.query(Backend).filter_by(name=backend_name).first()
            )

            if not src_type or not tgt_type or not backend:
                _log.error('Missing ID type or backend in database')
                return

            src_type_id = src_type.id
            tgt_type_id = tgt_type.id
            backend_id = backend.id

        # Load mapping data via the existing reader. For inputs_v2 backends
        # the limit is honoured at load time (via islice); other backends
        # ignore the unknown kwarg and are capped at COPY time below.
        limit = self._effective_limit(backend_name)
        reader_params = {}
        if limit is not None:
            reader_params['limit'] = limit
        reader = MapReader(
            id_type=id_type,
            target_id_type=target_id_type,
            ncbi_tax_id=ncbi_tax_id,
            backend=backend_name,
            **reader_params,
        )
        data = reader.load()

        if not data:
            _log.warning('No data loaded for %s -> %s', id_type, target_id_type)
            return

        # Delete existing rows for this combination
        with Session(self.engine) as session:
            session.execute(
                text(
                    f'DELETE FROM {SCHEMA}.id_mapping'
                    ' WHERE source_type_id = :src AND target_type_id = :tgt'
                    ' AND ncbi_tax_id = :tax AND backend_id = :bk'
                ),
                {
                    'src': src_type_id,
                    'tgt': tgt_type_id,
                    'tax': ncbi_tax_id,
                    'bk': backend_id,
                },
            )
            session.commit()

        # Bulk insert via COPY using psycopg3
        from omnipath_utils.db._connection import get_connection

        conn = get_connection(self._db_url)

        row_count = 0
        with conn.cursor() as cur:
            with cur.copy(
                f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
            ) as copy:
                for source_id, target_ids in data.items():
                    stored_source_id = (
                        normalize_identifier(id_type, source_id) or source_id
                    )
                    for target_id in target_ids:
                        stored_target_id = (
                            normalize_identifier(target_id_type, target_id)
                            or target_id
                        )
                        copy.write_row(
                            (
                                src_type_id,
                                tgt_type_id,
                                ncbi_tax_id,
                                stored_source_id[:64],
                                stored_target_id[:64],
                                backend_id,
                            )
                        )
                        row_count += 1
                    if limit is not None and row_count >= limit:
                        break

        conn.commit()
        conn.close()

        duration = time.time() - start

        # Record build info
        with Session(self.engine) as session:
            session.add(
                BuildInfo(
                    table_name='id_mapping',
                    source_type=id_type,
                    target_type=target_id_type,
                    ncbi_tax_id=ncbi_tax_id,
                    backend=backend_name,
                    row_count=row_count,
                    duration_secs=duration,
                    status='done',
                )
            )
            session.commit()

        _log.info(
            'Built %s -> %s: %d rows in %.1fs',
            id_type,
            target_id_type,
            row_count,
            duration,
        )

    def populate_reflists(self, ncbi_tax_id: int):
        """Populate reference lists for an organism.

        Loads SwissProt and TrEMBL ID sets via the reflists module and
        bulk-inserts them into the reflist table using COPY.
        """
        from omnipath_utils.reflists import all_trembls, all_swissprots
        from omnipath_utils.db._connection import get_connection

        start = time.time()

        with Session(self.engine) as session:
            uniprot_type = (
                session.query(IdType).filter_by(name='uniprot').first()
            )
            if not uniprot_type:
                _log.error('ID type "uniprot" not found, skipping reflists')
                return
            type_id = uniprot_type.id

        conn = get_connection(self._db_url)
        cur = conn.cursor()

        # Clear existing for this organism
        cur.execute(
            f'DELETE FROM {SCHEMA}.reflist WHERE ncbi_tax_id = %s',
            (ncbi_tax_id,),
        )
        conn.commit()

        total_rows = 0

        for list_name, loader in [
            ('swissprot', all_swissprots),
            ('trembl', all_trembls),
        ]:
            ids = loader(ncbi_tax_id)
            if self._max_records is not None:
                ids = set(list(ids)[: self._max_records])
            _log.info(
                'Loading %d %s IDs for organism %d',
                len(ids),
                list_name,
                ncbi_tax_id,
            )

            with cur.copy(
                f'COPY {SCHEMA}.reflist'
                ' (identifier, id_type_id, ncbi_tax_id, list_name)'
                ' FROM STDIN'
            ) as copy:
                for ac in ids:
                    copy.write_row((ac, type_id, ncbi_tax_id, list_name))

            conn.commit()
            total_rows += len(ids)

        duration = time.time() - start

        # Record build info
        with Session(self.engine) as session:
            session.add(
                BuildInfo(
                    table_name='reflist',
                    source_type='uniprot',
                    target_type='reflist (swissprot + trembl)',
                    ncbi_tax_id=ncbi_tax_id,
                    backend='reflists',
                    row_count=total_rows,
                    duration_secs=duration,
                    status='done',
                )
            )
            session.commit()

        conn.close()
        _log.info(
            'Reflists for organism %d: %d rows in %.1fs',
            ncbi_tax_id,
            total_rows,
            duration,
        )

    def populate_reflists_global_swissprot(self):
        """Load the global reviewed (SwissProt) AC set under ncbi_tax_id 0.

        SwissProt membership is organism-agnostic — an AC is reviewed or not,
        independent of taxon. One ``reviewed:true`` query (no organism) yields
        the complete reviewed set (~570k ACs), which drives SwissProt-preference
        for *every* organism in the resolver projection (FR-003), not just the
        handful with per-organism reflists.
        """
        from omnipath_utils.reflists import all_swissprots_global
        from omnipath_utils.db._connection import get_connection

        start = time.time()
        with Session(self.engine) as session:
            uniprot_type = (
                session.query(IdType).filter_by(name='uniprot').first()
            )
            if not uniprot_type:
                _log.error('ID type "uniprot" not found, skipping global reflist')
                return
            type_id = uniprot_type.id

        ids = all_swissprots_global()
        if self._max_records is not None:
            ids = set(list(ids)[: self._max_records])

        conn = get_connection(self._db_url)
        cur = conn.cursor()
        cur.execute(
            f"DELETE FROM {SCHEMA}.reflist WHERE ncbi_tax_id = 0 "
            f"AND list_name = 'swissprot'"
        )
        with cur.copy(
            f'COPY {SCHEMA}.reflist'
            ' (identifier, id_type_id, ncbi_tax_id, list_name)'
            ' FROM STDIN'
        ) as copy:
            for ac in ids:
                copy.write_row((ac, type_id, 0, 'swissprot'))
        conn.commit()
        conn.close()
        _log.info(
            'Global SwissProt reflist: %d ACs in %.1fs',
            len(ids),
            time.time() - start,
        )

    def load_uniprot_sec_ac(self):
        """Load the secondary -> primary UniProt AC map, organism-agnostic (tax 0).

        ``sec_ac.txt`` is a single global UniProt file (no organism), so it is
        loaded once at ``ncbi_tax_id = 0`` rather than per-organism (adding it to
        the per-organism ``PROTEIN_CORE`` fan-out would run it once per taxon and
        proteome-filter each, which is wrong). Every UniProt-bearing build needs
        it: the resolver views normalise every UniProt to its primary accession
        via this slice so the delivered mapping tables are primary-only (ADR 0006
        -- Utils owns normalization; the "uniprot cleanup" is always part of
        dealing with UniProts). Idempotent (DELETE + COPY the tax-0 slice).
        """
        from pypath.inputs.uniprot import get_uniprot_sec
        from omnipath_utils.db._connection import get_connection

        with Session(self.engine) as session:
            src = session.query(IdType).filter_by(name='uniprot-sec').first()
            tgt = session.query(IdType).filter_by(name='uniprot-pri').first()
            backend = session.query(Backend).filter_by(name='uniprot').first()
            if not src or not tgt or not backend:
                _log.error(
                    'load_uniprot_sec_ac: uniprot-sec/uniprot-pri id_type or '
                    'uniprot backend missing; skipping'
                )
                return
            src_id, tgt_id, backend_id = src.id, tgt.id, backend.id

        start = time.time()
        # organism=None -> every secondary/primary pair, no proteome filter.
        pairs = list(get_uniprot_sec(organism=None))
        if self._max_records is not None:
            pairs = pairs[: self._max_records]

        with Session(self.engine) as session:
            session.execute(
                text(
                    f'DELETE FROM {SCHEMA}.id_mapping WHERE source_type_id = :s '
                    'AND target_type_id = :t AND ncbi_tax_id = 0 '
                    'AND backend_id = :b'
                ),
                {'s': src_id, 't': tgt_id, 'b': backend_id},
            )
            session.commit()

        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for sec, pri in pairs:
                        copy.write_row(
                            (src_id, tgt_id, 0, sec[:64], pri[:64], backend_id)
                        )
                        n += 1
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d secondary->primary UniProt AC pairs (tax 0) in %.1fs',
            n, time.time() - start,
        )

    def load_uniprot_deleted(self):
        """Load the deleted-UniProt-AC marker list (007 US2), organism-agnostic.

        UniProt's deleted-accession lists (SwissProt always; TrEMBL only when
        ``OMNIPATH_UTILS_DELETED_TREMBL`` is set — it needs >5 GB) name accessions
        that were removed with **no successor**. They are stored self-referentially
        (``source_id = target_id``) under ``uniprot-deleted`` at tax 0 so the
        recovery stage can report such an AC as explicitly *deleted* rather than
        returning an unexplained empty result (FR-008). Idempotent (DELETE + COPY).
        """
        from pypath.inputs.uniprot import swissprot_deleted, trembl_deleted
        from omnipath_utils.db._connection import get_connection

        with Session(self.engine) as session:
            t = session.query(IdType).filter_by(name='uniprot-deleted').first()
            backend = session.query(Backend).filter_by(name='uniprot').first()
            if not t or not backend:
                _log.error(
                    'load_uniprot_deleted: uniprot-deleted id_type or uniprot '
                    'backend missing; skipping'
                )
                return
            type_id, backend_id = t.id, backend.id

        start = time.time()
        deleted = set(swissprot_deleted())
        if os.environ.get('OMNIPATH_UTILS_DELETED_TREMBL'):
            deleted |= set(trembl_deleted(confirm=False))
        deleted = sorted(a for a in deleted if a)
        if self._max_records is not None:
            deleted = deleted[: self._max_records]

        with Session(self.engine) as session:
            session.execute(
                text(
                    f'DELETE FROM {SCHEMA}.id_mapping WHERE source_type_id = :t '
                    'AND target_type_id = :t AND ncbi_tax_id = 0 '
                    'AND backend_id = :b'
                ),
                {'t': type_id, 'b': backend_id},
            )
            session.commit()

        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for ac in deleted:
                        copy.write_row(
                            (type_id, type_id, 0, ac[:64], ac[:64], backend_id)
                        )
                        n += 1
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d deleted UniProt AC markers (tax 0) in %.1fs',
            n, time.time() - start,
        )

    def load_gene_history(self, organisms: list[int] | None = None):
        """Load NCBI ``gene_history`` — retired Entrez -> current Entrez (007 US2).

        ``gene_history.gz`` gives, per taxon, a discontinued GeneID and the current
        GeneID it was merged into (``None`` when the gene was removed with no
        successor). Stored under ``entrez-history``:

        * merged: ``entrez-history`` ``discontinued`` -> ``entrez`` ``current``
          (the recovery target).
        * removed with no successor: self-referential ``entrez-history`` ->
          ``entrez-history`` so the recovery stage can flag it explicitly deleted.

        ``organisms`` restricts to a scope's taxa (``None`` = all taxa, complete
        scope). Idempotent (DELETE the whole ``entrez-history`` source slice + COPY).
        """
        from pypath.inputs.ncbi_gene import gene_history
        from omnipath_utils.db._connection import get_connection

        with Session(self.engine) as session:
            hist = session.query(IdType).filter_by(name='entrez-history').first()
            entrez = session.query(IdType).filter_by(name='entrez').first()
            backend = session.query(Backend).filter_by(name='gene_history').first()
            if not backend:
                backend = Backend(name='gene_history')
                session.add(backend)
                session.commit()
            if not hist or not entrez or not backend:
                _log.error(
                    'load_gene_history: entrez-history/entrez id_type or backend '
                    'missing; skipping'
                )
                return
            hist_id, entrez_id, backend_id = hist.id, entrez.id, backend.id

        taxa = set(organisms) if organisms else None
        start = time.time()

        with Session(self.engine) as session:
            session.execute(
                text(
                    f'DELETE FROM {SCHEMA}.id_mapping WHERE source_type_id = :h '
                    'AND backend_id = :b'
                ),
                {'h': hist_id, 'b': backend_id},
            )
            session.commit()

        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for rec in gene_history():
                        tax = rec.ncbi_tax_id or 0
                        if taxa is not None and tax not in taxa:
                            continue
                        old = rec.discontinued_entrez
                        if not old:
                            continue
                        if rec.entrez:  # merged -> current
                            copy.write_row((
                                hist_id, entrez_id, tax,
                                str(old)[:64], str(rec.entrez)[:64], backend_id,
                            ))
                        else:  # removed with no successor -> deleted marker
                            copy.write_row((
                                hist_id, hist_id, tax,
                                str(old)[:64], str(old)[:64], backend_id,
                            ))
                        n += 1
                        if self._max_records is not None and n >= self._max_records:
                            break
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d gene_history (retired->current Entrez) rows in %.1fs',
            n, time.time() - start,
        )

    #: Model organisms whose Ensembl gene history is loaded when no scope taxa
    #: are given (Ensembl main site species). Division species — plants, fungi,
    #: metazoa — live on a separate FTP and are skipped where absent.
    _ENSEMBL_HISTORY_DEFAULT_TAXA = (
        9606, 10090, 10116, 7955, 9913, 9823, 9031, 9615, 9544, 8364,
    )

    def load_ensembl_history(self, organisms: list[int] | None = None):
        """Load Ensembl ``stable_id_event`` — retired/merged Ensembl gene ids ->
        current gene id (007 US2). One species at a time from the Ensembl MySQL
        FTP:

        * merged / renamed: ``ensembl-history`` old -> ``ensg`` current (the
          recovery target).
        * removed with no successor: self-referential ``ensembl-history`` ->
          ``ensembl-history`` so the recovery stage can flag it deleted.

        ``organisms`` restricts to a scope's taxa; ``None`` uses a default set of
        model organisms Ensembl's main site covers (division species are on a
        separate FTP and are skipped). Idempotent (DELETE the whole
        ``ensembl-history`` source slice + COPY). Honours ``MAX_RECORDS``.
        """
        from pypath.inputs.ensembl import ensembl_id_history
        from omnipath_utils.db._connection import get_connection

        with Session(self.engine) as session:
            hist = (
                session.query(IdType).filter_by(name='ensembl-history').first()
            )
            ensg = session.query(IdType).filter_by(name='ensg').first()
            backend = (
                session.query(Backend).filter_by(name='ensembl_history').first()
            )
            if not backend:
                backend = Backend(name='ensembl_history')
                session.add(backend)
                session.commit()
            if not hist or not ensg or not backend:
                _log.error(
                    'load_ensembl_history: ensembl-history/ensg id_type or '
                    'backend missing; skipping'
                )
                return
            hist_id, ensg_id, backend_id = hist.id, ensg.id, backend.id

        taxa = list(organisms) if organisms else list(
            self._ENSEMBL_HISTORY_DEFAULT_TAXA
        )
        start = time.time()

        with Session(self.engine) as session:
            session.execute(
                text(
                    f'DELETE FROM {SCHEMA}.id_mapping WHERE source_type_id = :h '
                    'AND backend_id = :b'
                ),
                {'h': hist_id, 'b': backend_id},
            )
            session.commit()

        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for tax in taxa:
                        try:
                            events = ensembl_id_history(tax)
                        except Exception as e:
                            _log.error(
                                'load_ensembl_history: taxon %d failed: %s',
                                tax, e,
                            )
                            continue
                        for rec in events:
                            old = rec.old
                            if not old:
                                continue
                            if rec.new:  # merged / renamed -> current
                                copy.write_row((
                                    hist_id, ensg_id, tax,
                                    str(old)[:64], str(rec.new)[:64], backend_id,
                                ))
                            else:  # removed with no successor -> deleted marker
                                copy.write_row((
                                    hist_id, hist_id, tax,
                                    str(old)[:64], str(old)[:64], backend_id,
                                ))
                            n += 1
                            if (
                                self._max_records is not None
                                and n >= self._max_records
                            ):
                                break
                        if (
                            self._max_records is not None
                            and n >= self._max_records
                        ):
                            break
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d ensembl_history (retired->current Ensembl gene) rows '
            'in %.1fs',
            n, time.time() - start,
        )

    def load_gene2ensembl(self):
        """Load NCBI ``gene2ensembl`` — authoritative ensp/ensg -> entrez, all taxa.

        NCBI is the authority for Entrez Gene; this single all-organism file gives,
        per taxon, ENSP->Entrez and ENSG->Entrez DIRECT and for **every** transcript
        (unlike the UniProt idmapping, which cross-references only one canonical ENSP
        per protein and so misses the other-transcript ENSPs resources supply, e.g.
        STITCH). ``resolver_gene`` anchors Ensembl/gene ids to Entrez through THIS
        (gene space), not through UniProt. Ensembl ids are stored versionless.
        Idempotent (DELETE + COPY the ``gene2ensembl`` backend slice).
        """
        from pypath.inputs.ncbi_gene import gene2ensembl
        from omnipath_utils.db._connection import get_connection

        with Session(self.engine) as session:
            ensg = session.query(IdType).filter_by(name='ensg').first()
            ensp = session.query(IdType).filter_by(name='ensp').first()
            entrez = session.query(IdType).filter_by(name='entrez').first()
            backend = (
                session.query(Backend).filter_by(name='gene2ensembl').first()
            )
            if not backend:
                backend = Backend(name='gene2ensembl')
                session.add(backend)
                session.commit()
            if not ensg or not ensp or not entrez:
                _log.error(
                    'load_gene2ensembl: ensg/ensp/entrez id_type missing; skipping'
                )
                return
            ensg_id, ensp_id, entrez_id, backend_id = (
                ensg.id, ensp.id, entrez.id, backend.id,
            )

        start = time.time()
        with Session(self.engine) as session:
            session.execute(
                text(
                    f'DELETE FROM {SCHEMA}.id_mapping WHERE backend_id = :b'
                ),
                {'b': backend_id},
            )
            session.commit()

        # Deduplicate on the wire (the file has one row per transcript, so ENSG and
        # ENSP<->Entrez pairs repeat heavily).
        seen: set[tuple] = set()
        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for rec in gene2ensembl():
                        if rec.entrez is None:
                            continue
                        for src_id, sid in (
                            (ensp_id, rec.ensembl_protein),
                            (ensg_id, rec.ensembl_gene),
                        ):
                            if not sid:
                                continue
                            key = (src_id, rec.ncbi_tax_id, sid, rec.entrez)
                            if key in seen:
                                continue
                            seen.add(key)
                            copy.write_row((
                                src_id, entrez_id, rec.ncbi_tax_id,
                                sid[:64], rec.entrez[:64], backend_id,
                            ))
                            n += 1
                            if self._max_records and n >= self._max_records:
                                break
                        if self._max_records and n >= self._max_records:
                            break
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d gene2ensembl ensp/ensg->entrez rows in %.1fs',
            n, time.time() - start,
        )

    def load_gene_info(self, organisms=None):
        """Load NCBI ``gene_info`` — authoritative genesymbol/synonym -> entrez.

        The true gene-space symbol<->entrez map for **all** organisms (covers
        genes with no UniProt), closing the ``resolver_protein.sql`` follow-up
        gap. ``organisms`` limits the loaded taxa (``None`` = all, for the
        ``complete`` scope). Idempotent (DELETE + COPY the ``gene_info`` slice).
        """
        from pypath.inputs.ncbi_gene import gene_info
        from omnipath_utils.db._connection import get_connection

        taxa = set(organisms) if organisms else None

        with Session(self.engine) as session:
            symbol = session.query(IdType).filter_by(name='genesymbol').first()
            syn = (
                session.query(IdType).filter_by(name='genesymbol-syn').first()
            )
            entrez = session.query(IdType).filter_by(name='entrez').first()
            backend = (
                session.query(Backend).filter_by(name='gene_info').first()
            )
            if not backend:
                backend = Backend(name='gene_info')
                session.add(backend)
                session.commit()
            if not symbol or not entrez:
                _log.error(
                    'load_gene_info: genesymbol/entrez id_type missing; skipping'
                )
                return
            symbol_id, entrez_id, backend_id = (
                symbol.id, entrez.id, backend.id,
            )
            syn_id = syn.id if syn else None

        start = time.time()
        with Session(self.engine) as session:
            session.execute(
                text(f'DELETE FROM {SCHEMA}.id_mapping WHERE backend_id = :b'),
                {'b': backend_id},
            )
            session.commit()

        seen: set[tuple] = set()
        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for rec in gene_info():
                        if rec.entrez is None:
                            continue
                        if taxa is not None and rec.ncbi_tax_id not in taxa:
                            continue
                        pairs = []
                        if rec.symbol:
                            pairs.append((symbol_id, rec.symbol))
                        if syn_id is not None and rec.synonyms:
                            pairs.extend(
                                (syn_id, s) for s in rec.synonyms if s
                            )
                        for src_id, sval in pairs:
                            key = (src_id, rec.ncbi_tax_id, sval, rec.entrez)
                            if key in seen:
                                continue
                            seen.add(key)
                            copy.write_row((
                                src_id, entrez_id, rec.ncbi_tax_id,
                                sval[:64], rec.entrez[:64], backend_id,
                            ))
                            n += 1
                        if self._max_records and n >= self._max_records:
                            break
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d gene_info genesymbol/syn->entrez rows in %.1fs',
            n, time.time() - start,
        )

    def load_gene2accession(self, organisms=None):
        """Load NCBI ``gene2accession`` — RefSeq RNA/protein -> entrez.

        Gene-space RefSeq (``refseqn`` RNA, ``refseqp`` protein) <-> Entrez for
        the scope's organisms (``None`` = all). Accessions are versionless; the
        file repeats a gene per assembly, so pairs are de-duplicated on the wire.
        Idempotent (DELETE + COPY the ``gene2accession`` slice).
        """
        from pypath.inputs.ncbi_gene import gene2accession
        from omnipath_utils.db._connection import get_connection

        taxa = set(organisms) if organisms else None

        with Session(self.engine) as session:
            refseqn = session.query(IdType).filter_by(name='refseqn').first()
            refseqp = session.query(IdType).filter_by(name='refseqp').first()
            entrez = session.query(IdType).filter_by(name='entrez').first()
            backend = (
                session.query(Backend).filter_by(name='gene2accession').first()
            )
            if not backend:
                backend = Backend(name='gene2accession')
                session.add(backend)
                session.commit()
            if not entrez or (not refseqn and not refseqp):
                _log.error(
                    'load_gene2accession: entrez/refseq id_type missing; skipping'
                )
                return
            entrez_id, backend_id = entrez.id, backend.id
            refseqn_id = refseqn.id if refseqn else None
            refseqp_id = refseqp.id if refseqp else None

        start = time.time()
        with Session(self.engine) as session:
            session.execute(
                text(f'DELETE FROM {SCHEMA}.id_mapping WHERE backend_id = :b'),
                {'b': backend_id},
            )
            session.commit()

        seen: set[tuple] = set()
        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for rec in gene2accession():
                        if rec.entrez is None:
                            continue
                        if taxa is not None and rec.ncbi_tax_id not in taxa:
                            continue
                        for src_id, sval in (
                            (refseqn_id, rec.rna_refseq),
                            (refseqp_id, rec.protein_refseq),
                        ):
                            if src_id is None or not sval:
                                continue
                            key = (src_id, rec.ncbi_tax_id, sval, rec.entrez)
                            if key in seen:
                                continue
                            seen.add(key)
                            copy.write_row((
                                src_id, entrez_id, rec.ncbi_tax_id,
                                sval[:64], rec.entrez[:64], backend_id,
                            ))
                            n += 1
                        if self._max_records and n >= self._max_records:
                            break
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d gene2accession refseq->entrez rows in %.1fs',
            n, time.time() - start,
        )

    def load_kegg_genes(self, organisms=None):
        """Load KEGG gene ids -> entrez / uniprot for the scope's organisms.

        KEGG gene relations are organism-specific (one REST ``conv`` call per
        organism), so this iterates the scope organisms, resolving each to its
        KEGG organism code. KEGG gene ids are stored in canonical ``<org>:<n>``
        form. Idempotent (DELETE the ``kegg_gene`` backend slice then COPY).
        """
        from pypath.inputs.kegg import (
            kegg_gene_to_ncbi_gene,
            kegg_gene_to_uniprot,
        )
        from omnipath_utils.db._connection import get_connection

        organisms = organisms or [9606]
        tm = TaxonomyManager.get()

        with Session(self.engine) as session:
            kegg_gene = (
                session.query(IdType).filter_by(name='kegg_gene').first()
            )
            entrez = session.query(IdType).filter_by(name='entrez').first()
            uniprot = session.query(IdType).filter_by(name='uniprot').first()
            backend = (
                session.query(Backend).filter_by(name='kegg_gene').first()
            )
            if not backend:
                backend = Backend(name='kegg_gene')
                session.add(backend)
                session.commit()
            if not kegg_gene or not entrez:
                _log.error(
                    'load_kegg_genes: kegg_gene/entrez id_type missing; skipping'
                )
                return
            kegg_gene_id, entrez_id, backend_id = (
                kegg_gene.id, entrez.id, backend.id,
            )
            uniprot_id = uniprot.id if uniprot else None

        start = time.time()
        with Session(self.engine) as session:
            session.execute(
                text(f'DELETE FROM {SCHEMA}.id_mapping WHERE backend_id = :b'),
                {'b': backend_id},
            )
            session.commit()

        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for org in organisms:
                        code = tm.ensure_kegg_code(org)
                        if not code:
                            _log.warning(
                                'load_kegg_genes: no KEGG code for organism %d',
                                org,
                            )
                            continue
                        targets = [(entrez_id, kegg_gene_to_ncbi_gene(code))]
                        if uniprot_id is not None:
                            targets.append(
                                (uniprot_id, kegg_gene_to_uniprot(code)),
                            )
                        for tgt_id, data in targets:
                            for src_val, tgt_vals in data.items():
                                for tval in tgt_vals:
                                    copy.write_row((
                                        kegg_gene_id, tgt_id, org,
                                        src_val[:64], tval[:64], backend_id,
                                    ))
                                    n += 1
                        if self._max_records and n >= self._max_records:
                            break
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d kegg_gene->entrez/uniprot rows in %.1fs',
            n, time.time() - start,
        )

    def load_mane(self, organisms=None):
        """Load the MANE matched RefSeq<->Ensembl transcript set (human).

        MANE gives, per human gene, one agreed-upon transcript with matching
        RefSeq and Ensembl RNA/protein accessions plus the Entrez/Ensembl gene
        and symbol. Two kinds of rows are stored (007 US3):

        * gene anchors — refseqn/refseqp/enst/ensp/ensg/genesymbol -> ``entrez``
          (the ``resolver_gene`` picks these up like any other gene-space pair);
        * transcript cross-links — ``enst`` <-> ``refseqn`` and ``ensp`` <->
          ``refseqp`` (both directions) so a RefSeq transcript translates to its
          Ensembl counterpart and back.

        MANE is human-only; skipped when 9606 is not in ``organisms`` (``None``
        or the ``complete`` scope always loads it). Idempotent (DELETE the
        ``mane`` backend slice then COPY). Honours ``MAX_RECORDS``.
        """
        if organisms is not None and 9606 not in organisms:
            _log.info('load_mane: human (9606) not in scope; skipping')
            return

        from pypath.inputs.mane import mane_summary
        from omnipath_utils.db._connection import get_connection

        with Session(self.engine) as session:
            ids = {
                name: (t.id if t else None)
                for name, t in (
                    (n, session.query(IdType).filter_by(name=n).first())
                    for n in (
                        'entrez', 'refseqn', 'refseqp', 'enst', 'ensp', 'ensg',
                        'genesymbol',
                    )
                )
            }
            backend = session.query(Backend).filter_by(name='mane').first()
            if not backend:
                backend = Backend(name='mane')
                session.add(backend)
                session.commit()
            if not ids['entrez']:
                _log.error('load_mane: entrez id_type missing; skipping')
                return
            backend_id = backend.id

        entrez_id = ids['entrez']
        start = time.time()
        with Session(self.engine) as session:
            session.execute(
                text(f'DELETE FROM {SCHEMA}.id_mapping WHERE backend_id = :b'),
                {'b': backend_id},
            )
            session.commit()

        seen: set[tuple] = set()
        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for rec in mane_summary():
                        if rec.entrez is None:
                            continue
                        vals = {
                            'refseqn': rec.refseq_rna,
                            'refseqp': rec.refseq_protein,
                            'enst': rec.ensembl_transcript,
                            'ensp': rec.ensembl_protein,
                            'ensg': rec.ensembl_gene,
                            'genesymbol': rec.symbol,
                        }
                        rows = []
                        # gene anchors: every present id -> entrez
                        for name, val in vals.items():
                            if ids[name] is not None and val:
                                rows.append((ids[name], entrez_id, val, rec.entrez))
                        # transcript cross-links (both directions)
                        for a, b in (('enst', 'refseqn'), ('ensp', 'refseqp')):
                            if (
                                ids[a] is not None and ids[b] is not None
                                and vals[a] and vals[b]
                            ):
                                rows.append((ids[a], ids[b], vals[a], vals[b]))
                                rows.append((ids[b], ids[a], vals[b], vals[a]))
                        for src_id, tgt_id, sval, tval in rows:
                            key = (src_id, tgt_id, sval, tval)
                            if key in seen:
                                continue
                            seen.add(key)
                            copy.write_row((
                                src_id, tgt_id, rec.ncbi_tax_id,
                                sval[:64], tval[:64], backend_id,
                            ))
                            n += 1
                        if self._max_records and n >= self._max_records:
                            break
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d MANE RefSeq<->Ensembl rows in %.1fs',
            n, time.time() - start,
        )

    def load_taxonomy_species(self):
        """Load the NCBI Taxonomy strain -> species map into ``taxon_species``.

        Walks ``nodes.dmp`` up the taxonomic tree so any ``tax_id`` (strain,
        sub-species, or the species itself) resolves to its species-level
        ancestor. A species maps to itself; a genus-or-above taxid with no
        species-level ancestor keeps its own taxid. Retired taxids (``merged.dmp``)
        are canonicalised to the species of their current taxid too, so lookups
        on an old id still hit. Organism-agnostic (no scope). Idempotent (DELETE
        the whole derived table + COPY). Honours ``MAX_RECORDS``.
        """
        from pypath.inputs.ncbi_taxonomy import taxonomy_nodes, taxonomy_merged
        from omnipath_utils.db._connection import get_connection

        start = time.time()

        # retired taxid -> current taxid
        merged = {rec.old_tax_id: rec.new_tax_id for rec in taxonomy_merged()}

        # the taxonomic tree (parent) and each taxon's rank
        parent: dict[int, int] = {}
        ranks: dict[int, str] = {}
        for rec in taxonomy_nodes():
            parent[rec.tax_id] = rec.parent_tax_id
            ranks[rec.tax_id] = rec.rank

        def species_of(t):
            t = merged.get(t, t)
            seen_t: set[int] = set()
            cur = t
            while cur is not None and cur not in seen_t:
                seen_t.add(cur)
                if ranks.get(cur) == 'species':
                    return cur
                p = parent.get(cur)
                if p is None or p == cur:  # reached root (taxid 1 is its own parent)
                    break
                cur = p
            return t  # no species-level ancestor -> keep original

        with Session(self.engine) as session:
            session.execute(text(f'DELETE FROM {SCHEMA}.taxon_species'))
            session.commit()

        seen: set[int] = set()
        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.taxon_species (tax_id, species_tax_id) '
                    'FROM STDIN'
                ) as copy:
                    # every known taxon -> its species (a species maps to itself)
                    for tax_id in ranks:
                        if tax_id in seen:
                            continue
                        seen.add(tax_id)
                        copy.write_row((tax_id, species_of(tax_id)))
                        n += 1
                        if self._max_records and n >= self._max_records:
                            break
                    else:
                        # retired taxid -> species of its current taxid
                        for old_tax_id, new_tax_id in merged.items():
                            if old_tax_id in seen:
                                continue
                            seen.add(old_tax_id)
                            copy.write_row((old_tax_id, species_of(new_tax_id)))
                            n += 1
                            if self._max_records and n >= self._max_records:
                                break
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d taxon->species rows in %.1fs',
            n, time.time() - start,
        )

    def load_rnacentral_mirbase_gene(self, organisms=None):
        """Load an RNAcentral miRBase -> gene bridge into ``id_mapping``.

        Cross-references miRBase accessions (``MI####``/``MIMAT####``) to Ensembl
        gene ids (``ensg``) and gene symbols (``genesymbol``) via their shared
        RNAcentral ``URS`` accession: ``ensembl.tsv`` and ``hgnc.tsv`` give
        ``URS -> ENSG``/``URS -> symbol``, ``mirbase.tsv`` gives ``URS ->
        miRBase``. ``organisms`` restricts to a scope's taxa (``None`` = all).
        Stored under the ``rnacentral`` backend. Idempotent (DELETE the backend
        slice + COPY). Honours ``MAX_RECORDS``.
        """
        from pypath.inputs.rnacentral import rnacentral_xref
        from omnipath_utils.db._connection import get_connection

        taxa = set(organisms) if organisms else None

        with Session(self.engine) as session:
            mirbase = session.query(IdType).filter_by(name='mirbase').first()
            ensg = session.query(IdType).filter_by(name='ensg').first()
            symbol = session.query(IdType).filter_by(name='genesymbol').first()
            backend = (
                session.query(Backend).filter_by(name='rnacentral').first()
            )
            if not backend:
                backend = Backend(name='rnacentral')
                session.add(backend)
                session.commit()
            if not mirbase:
                _log.error(
                    'load_rnacentral_mirbase_gene: mirbase id_type missing; '
                    'skipping'
                )
                return
            mirbase_id, backend_id = mirbase.id, backend.id
            ensg_id = ensg.id if ensg else None
            symbol_id = symbol.id if symbol else None

        start = time.time()

        # URS -> ENSG / symbol, scoped to the requested taxa
        urs_to_ensg: dict[str, set[str]] = {}
        for rec in rnacentral_xref('ensembl'):
            if taxa is not None and rec.ncbi_tax_id not in taxa:
                continue
            if not rec.rnacentral or not rec.gene_name:
                continue
            ensg_val = rec.gene_name.split('.', 1)[0]
            if ensg_val:
                urs_to_ensg.setdefault(rec.rnacentral, set()).add(ensg_val)

        urs_to_symbol: dict[str, set[str]] = {}
        for rec in rnacentral_xref('hgnc'):
            if taxa is not None and rec.ncbi_tax_id not in taxa:
                continue
            if not rec.rnacentral or not rec.gene_name:
                continue
            urs_to_symbol.setdefault(rec.rnacentral, set()).add(rec.gene_name)

        with Session(self.engine) as session:
            session.execute(
                text(f'DELETE FROM {SCHEMA}.id_mapping WHERE backend_id = :b'),
                {'b': backend_id},
            )
            session.commit()

        seen: set[tuple] = set()
        conn = get_connection(self._db_url)
        n = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping (source_type_id, target_type_id, '
                    'ncbi_tax_id, source_id, target_id, backend_id) FROM STDIN'
                ) as copy:
                    for rec in rnacentral_xref('mirbase'):
                        if taxa is not None and rec.ncbi_tax_id not in taxa:
                            continue
                        urs = rec.rnacentral
                        mi = rec.external_id
                        if not urs or not mi:
                            continue
                        rows = []
                        if ensg_id is not None:
                            for ensg_val in urs_to_ensg.get(urs, ()):
                                rows.append((ensg_id, ensg_val))
                        if symbol_id is not None:
                            for sym in urs_to_symbol.get(urs, ()):
                                rows.append((symbol_id, sym))
                        for tgt_id, tval in rows:
                            key = (tgt_id, rec.ncbi_tax_id, mi, tval)
                            if key in seen:
                                continue
                            seen.add(key)
                            copy.write_row((
                                mirbase_id, tgt_id, rec.ncbi_tax_id,
                                mi[:64], tval[:64], backend_id,
                            ))
                            n += 1
                        if self._max_records and n >= self._max_records:
                            break
            conn.commit()
        finally:
            conn.close()

        _log.info(
            'Loaded %d RNAcentral miRBase->gene rows in %.1fs',
            n, time.time() - start,
        )

    #: The canonical resolver projections, each in its own SQL file so
    #: :meth:`create_resolver_views` can rebuild them independently (007 R10 /
    #: Phase 3P). Order matters — later stages read earlier ones.
    RESOLVER_NAMES = ('gene', 'global', 'protein', 'combined', 'chemical')

    def create_resolver_views(
        self,
        force_ftp_core: bool = False,
        resolvers: 'set[str] | list[str] | None' = None,
    ):
        """(Re)create the canonical resolver projection views (SQL DDL).

        Read by omnipath-build via DuckDB ATTACH (spec 002/003; idempotent):

        * ``resolver_gene`` — any in-scope source id -> its **NCBI Gene (Entrez)
          anchor** per taxon (US7). A VIEW = ``resolver_gene_ftp`` UNION ALL
          ``resolver_gene_curated`` (spec 007 R10 / Phase 3P).
        * ``resolver_protein`` — per-taxon ``source_id -> UniProt`` (primary
          SwissProt where available); the SQL replacement for the Python uniprot
          cleanup in DB-backed mode.
        * ``resolver_gene_protein_global`` — taxon-agnostic UniProt/Entrez ->
          Entrez (the full global slice, T069/R25).
        * ``resolver_uniprot_taxon`` — UniProt accession -> organism, so the build
          can recover the organism of a protein mention that cited an accession
          but no species. An FTP core alongside ``resolver_gene_ftp``.
        * ``resolver_gene_protein_combined`` — one canonical target per key (only
          the pre-built parquet/duckdb resolver mode reads it, not the live-PG
          keyed path).
        * ``resolver_chemical`` — chemical ``source_id -> InChIKey`` (spec 003 R7).

        Applied at the end of every build mode so an additive/incremental load
        never leaves the views stale or missing.

        Selective + gated rebuild (007 R10 / Phase 3P, T065/T067). The projections
        are split across per-resolver SQL files:

        * FTP cores — ``resolver_gene_ftp``, ``resolver_gene_protein_global`` and
          ``resolver_uniprot_taxon`` — derive from the ~46 GB ``id_mapping_ftp``
          (the ~3 h cost). They change
          ONLY on an FTP reload: a full FTP load DROPs id_mapping_ftp CASCADE,
          dropping them, so their mere absence is the rebuild trigger. Rebuilt only
          when missing or ``force_ftp_core`` (which :meth:`populate_from_ftp`
          passes). ``resolver_gene`` (a VIEW) and ``resolver_gene_curated`` /
          ``resolver_gene_protein_combined`` fall to that CASCADE too, so they are
          always re-derived below.
        * Curated deltas — ``resolver_gene_curated`` (+ the ``resolver_gene`` view),
          ``resolver_protein``, ``resolver_gene_protein_combined``,
          ``resolver_chemical`` — rebuild from the curated ``id_mapping``.

        ``resolvers`` selects which projections to (re)build: a subset of
        :attr:`RESOLVER_NAMES`. ``None`` -> the ``OMNIPATH_UTILS_RESOLVERS`` env
        (comma-separated) if set, else all. For fast P2 iteration on a live utils
        DB, ``resolvers={'gene'}`` rebuilds only the gene FTP core (if the FTP was
        reloaded) + the ~166 s curated delta, skipping the 29-34 GB
        protein/combined/global rebuilds that a gene-space load does not affect.
        The caller owns cross-resolver dependencies: ``combined`` reads the
        ``resolver_gene`` view, so build ``gene`` too (or ensure it already exists).
        """
        import importlib.resources as ir

        from omnipath_utils.db._connection import get_connection

        if resolvers is None:
            env = os.environ.get('OMNIPATH_UTILS_RESOLVERS')
            resolvers = (
                {r.strip() for r in env.split(',') if r.strip()}
                if env
                else set(self.RESOLVER_NAMES)
            )
        resolvers = set(resolvers)
        unknown = resolvers - set(self.RESOLVER_NAMES)
        if unknown:
            raise ValueError(
                f'unknown resolver(s) {sorted(unknown)}; '
                f'valid: {self.RESOLVER_NAMES}'
            )

        conn = get_connection(self._db_url)
        cur = conn.cursor()

        # Session tuning for the resolver (re)builds. The FTP-core file also
        # SET LOCALs, but a delta-only run skips it, so the curated rebuilds would
        # otherwise fall back to the default 64 MB work_mem and spill. SET LOCAL:
        # scoped to this transaction, reset at commit.
        for guc, val in (
            ('work_mem', os.environ.get('OMNIPATH_UTILS_RESOLVER_WORK_MEM', '4GB')),
            ('maintenance_work_mem', '4GB'),
            ('max_parallel_workers_per_gather', '8'),
            ('max_parallel_workers', '8'),
            ('max_parallel_maintenance_workers', '8'),
        ):
            cur.execute(f"SET LOCAL {guc} = '{val}'")

        def present(rel: str) -> bool:
            cur.execute(f"SELECT to_regclass('{SCHEMA}.{rel}') IS NOT NULL")
            return bool(cur.fetchone()[0])

        # (label, sql file) in dependency order; FTP cores gated on absence.
        plan: list[tuple[str, str]] = []
        if 'gene' in resolvers:
            if force_ftp_core or not present('resolver_gene_ftp'):
                plan.append(('gene_ftp', 'sql/resolver_gene_ftp.sql'))
            else:
                _log.info(
                    'resolver_gene_ftp present, no FTP reload — skipping the FTP '
                    'core; rebuilding the curated gene delta only',
                )
            # accession -> organism, an FTP core alongside resolver_gene_ftp:
            # both derive from id_mapping_ftp and change only on an FTP reload.
            if force_ftp_core or not present('resolver_uniprot_taxon'):
                plan.append(('uniprot_taxon', 'sql/resolver_uniprot_taxon.sql'))
            else:
                _log.info(
                    'resolver_uniprot_taxon present, no FTP reload — skipping',
                )
            plan.append(('gene_curated', 'sql/resolver_gene_curated.sql'))
        if 'global' in resolvers:
            if force_ftp_core or not present('resolver_gene_protein_global'):
                plan.append(('global', 'sql/resolver_global.sql'))
            else:
                _log.info(
                    'resolver_gene_protein_global present, no FTP reload — skipping',
                )
        if 'protein' in resolvers:
            plan.append(('protein', 'sql/resolver_protein.sql'))
        if 'combined' in resolvers:
            plan.append(('combined', 'sql/resolver_combined.sql'))
        if 'chemical' in resolvers:
            plan.append(('chemical', 'sql/resolver_chemical.sql'))

        # A gene rebuild DROP ... CASCADEs resolver_gene_protein_combined (it reads
        # the resolver_gene view). If the caller did not also ask to rebuild it,
        # warn: it is left absent until the next full run (harmless for the live
        # keyed path, which does not read it).
        if (
            'gene' in resolvers
            and 'combined' not in resolvers
            and present('resolver_gene_protein_combined')
        ):
            _log.warning(
                'rebuilding gene will drop resolver_gene_protein_combined (it '
                'reads the resolver_gene view) and it is NOT in this rebuild set — '
                'it will be absent until the next full build',
            )

        for label, sql_file in plan:
            _log.info('resolver rebuild: %s', label)
            sql = (
                ir.files('omnipath_utils.db')
                .joinpath(sql_file)
                .read_text(encoding='utf-8')
            )
            cur.execute(sql)
        conn.commit()
        conn.close()
        _log.info(
            'Resolver projections rebuilt: %s',
            ', '.join(label for label, _ in plan) or '(none)',
        )

    def build_reference_tables(self):
        """Build all reference tables (id_types, backends, organisms)."""
        self.create_tables()
        self.populate_id_types()
        self.populate_backends()
        self.populate_organisms()

    def _record_build_scope(self, scope, organisms):
        """Record the 007 build scope in ``build_info`` (idempotent).

        Stores one ``table_name='build_scope'`` row: the scope label in
        ``source_type`` and the organism count in ``row_count`` (0 for the
        all-taxa ``complete`` scope). Lets a coverage report / consumer state
        which scope produced the DB.
        """
        label = scope or 'custom'
        taxa = organisms if organisms is not None else []
        with Session(self.engine) as session:
            session.execute(
                text(
                    f"DELETE FROM {SCHEMA}.build_info"
                    " WHERE table_name = 'build_scope'"
                )
            )
            session.add(
                BuildInfo(
                    table_name='build_scope',
                    source_type=str(label)[:64],
                    target_type=(
                        'complete' if scope == 'complete' else 'organisms'
                    ),
                    ncbi_tax_id=0,
                    backend='cli',
                    row_count=len(taxa),
                    status='done',
                )
            )
            session.commit()
        _log.info(
            'build scope = %s (%s)',
            label,
            'complete/all-taxa'
            if scope == 'complete'
            else f'{len(taxa)} organism(s)',
        )

    def build_all(
        self,
        organisms: list[int] | None = None,
        scope: str | None = None,
    ):
        """Full build: reference tables + mappings + reflists.

        Uses preset infrastructure internally. Equivalent to a custom
        preset with the given organisms, protein core + Ensembl mappings,
        and reference lists.

        ``scope`` is the 007 build-scope label (``only-human`` / ``core-model``
        / ``extended-model`` / ``model-organisms`` / ``complete``, or a custom
        organism list) recorded in ``build_info``. For ``complete`` scope the
        per-organism REST protein loop is skipped — protein mappings come from
        the all-organism FTP idmapping (``populate_from_ftp`` / already present
        in a cumulative base) plus the organism-agnostic gene-space anchors
        loaded below.
        """
        from omnipath_utils.db._presets import (
            ENSEMBL,
            ENSEMBL_GENOMES,
            PROTEIN_CORE,
        )
        from omnipath_utils.mapping.backends._biomart import ORGANISM_DIVISION

        self.build_reference_tables()

        complete = scope == 'complete'
        if not complete:
            organisms = organisms or [9606]  # default: human only

        self._record_build_scope(scope, organisms if not complete else None)

        if complete:
            _log.info(
                'build_all: complete scope — protein mappings from the FTP '
                'idmapping base + all-taxa gene-space anchors; skipping the '
                'per-organism REST mapping/reflist loop',
            )
            organisms = []

        if organisms:
            self._run_mappings_parallel(PROTEIN_CORE + ENSEMBL, organisms)
            # Ensembl Genomes divisions (plants/fungi/metazoa/protists): the
            # eg-mart datasets only exist for these organisms, so run the
            # ensgg/ensgp/ensgt pairs on the division subset of the scope.
            division_organisms = [
                o for o in organisms if o in ORGANISM_DIVISION
            ]
            if division_organisms:
                self._run_mappings_parallel(
                    ENSEMBL_GENOMES, division_organisms,
                )

        # Reference lists
        for org in organisms:
            try:
                self.populate_reflists(org)
            except Exception as e:
                _log.error('Failed reflists for org %d: %s', org, e)

        # Global reviewed (SwissProt) set — drives organism-agnostic
        # SwissProt-preference in the resolver projection (FR-003).
        try:
            self.populate_reflists_global_swissprot()
        except Exception as e:
            _log.error('Failed global SwissProt reflist: %s', e)

        # Secondary->primary UniProt normalization, always part of handling
        # proteins (ADR 0006). Organism-agnostic, loaded once at tax 0.
        try:
            self.load_uniprot_sec_ac()
        except Exception as e:
            _log.error('Failed to load sec_ac (secondary->primary UniProt): %s', e)

        # Deprecated-ID recovery tables (007 US2): deleted UniProt AC markers
        # (organism-agnostic) + retired->current Entrez (gene_history, scope taxa).
        try:
            self.load_uniprot_deleted()
        except Exception as e:
            _log.error('Failed to load uniprot_deleted: %s', e)
        try:
            self.load_gene_history(organisms=None if complete else organisms)
        except Exception as e:
            _log.error('Failed to load gene_history: %s', e)
        try:
            self.load_ensembl_history(
                organisms=None if complete else organisms,
            )
        except Exception as e:
            _log.error('Failed to load ensembl_history: %s', e)

        # Authoritative ensp/ensg->entrez (NCBI gene2ensembl) — gene-space anchor,
        # all organisms, every transcript (the ENSP-coverage fix).
        try:
            self.load_gene2ensembl()
        except Exception as e:
            _log.error('Failed to load gene2ensembl: %s', e)

        # Complete gene-space symbol/synonym + RefSeq -> entrez (007 US1).
        # None = all taxa (complete scope); otherwise the scope's organisms.
        gene_space_taxa = None if complete else organisms
        try:
            self.load_gene_info(organisms=gene_space_taxa)
        except Exception as e:
            _log.error('Failed to load gene_info: %s', e)
        try:
            self.load_gene2accession(organisms=gene_space_taxa)
        except Exception as e:
            _log.error('Failed to load gene2accession: %s', e)

        # MANE matched RefSeq<->Ensembl transcript set (human; 007 US3). Loaded
        # whenever human is in scope (complete always).
        try:
            self.load_mane(organisms=None if complete else organisms)
        except Exception as e:
            _log.error('Failed to load MANE: %s', e)

        # NCBI Taxonomy strain->species map (organism-agnostic).
        try:
            self.load_taxonomy_species()
        except Exception as e:
            _log.error('Failed to load taxonomy_species: %s', e)

        # RNAcentral miRBase->gene (ensg/genesymbol) bridge.
        try:
            self.load_rnacentral_mirbase_gene(
                organisms=None if complete else organisms
            )
        except Exception as e:
            _log.error('Failed to load rnacentral_mirbase_gene: %s', e)

        # KEGG gene ids -> entrez/uniprot. Per-organism REST calls, so only for
        # an explicit organism list (skipped for complete: KEGG has no bulk
        # all-organism gene conv, and iterating every NCBI taxon is infeasible).
        if organisms:
            try:
                self.load_kegg_genes(organisms=organisms)
            except Exception as e:
                _log.error('Failed to load kegg_genes: %s', e)

    # FTP idmapping labels NOT loaded by default (annotation/clustering, not ID
    # translation; they 10x the table). Opt in via OMNIPATH_BUILD_FTP_HEAVY_TYPES.
    _FTP_HEAVY_LABELS = {'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'PDB', 'STRING'}

    def _pg_set(self, cur):
        """Apply session-level (not global) tuning for the in-DB transform."""
        gucs = {
            'work_mem': os.environ.get('OMNIPATH_BUILD_WORK_MEM', '1GB'),
            'maintenance_work_mem': os.environ.get(
                'OMNIPATH_BUILD_MAINTENANCE_WORK_MEM', '4GB'
            ),
            'max_parallel_workers_per_gather': os.environ.get(
                'OMNIPATH_BUILD_MAX_PARALLEL_WORKERS_PER_GATHER', '8'
            ),
            'max_parallel_maintenance_workers': os.environ.get(
                'OMNIPATH_BUILD_MAX_PARALLEL_MAINTENANCE_WORKERS', '8'
            ),
            # Per-worker private hashes instead of one shared hash table in
            # /dev/shm: keeps the billion-row hash join within a modest
            # container shm regardless of deployment.
            'enable_parallel_hash': os.environ.get(
                'OMNIPATH_BUILD_ENABLE_PARALLEL_HASH', 'off'
            ),
            'synchronous_commit': 'off',
        }
        for k, v in gucs.items():
            # SET does not accept bound parameters; values are operator-controlled
            # env/defaults. Quote as a string literal (PG coerces ints/sizes).
            safe = str(v).replace("'", "")
            cur.execute(f"SET {k} = '{safe}'")

    def populate_from_ftp(self, id_types: set[str] | None = None):
        """Load the full UniProt FTP idmapping into a swappable ``id_mapping_ftp``.

        The complete ``idmapping.dat.gz`` (~18 GB, all organisms incl. the long
        tail with no per-organism file) is loaded by **streaming it into a staging
        table and transforming in-database** — no per-row Python, no awk/split:

        1. ``pypath.inputs.uniprot_ftp.stream_full_idmapping`` downloads (dlmachine,
           cached) and streams the **decompressed** bytes; the build client pipes
           them straight into a client-side ``COPY … FROM STDIN`` of the raw
           ``(ac, id_type_label, id_value)`` rows into an UNLOGGED staging table
           (portable: the DB server needs no file access).
        2. In-DB, set-based: build ``ac → taxid`` from the ``NCBI_TaxID`` rows and a
           ``id_type_label → id_type_id`` map (synonyms included), then build the
           new table resolving ``id_type`` and ``ncbi_tax_id`` as **foreign keys**.
        3. Build its indexes **offline**, then **atomically swap** it in
           (``id_mapping_ftp``) behind the ``id_mapping_all`` view — the curated
           ``id_mapping`` keeps serving with full indexes throughout.

        Heavy non-translation labels (GO/UniRef/PDB/STRING) are excluded unless
        ``OMNIPATH_BUILD_FTP_HEAVY_TYPES=1``. Session tuning via ``OMNIPATH_BUILD_*``
        (see ``omnipath-build/docs/build-tuning.md``).
        """
        from omnipath_utils.db._connection import SCHEMA, get_connection
        from pypath.inputs.uniprot_ftp import IDTYPE_MAP, stream_full_idmapping

        _log.info('Starting full FTP idmapping build (staging + in-DB transform)')
        start = time.time()

        # Backend id + uniprot source type id + label -> id_type_id map.
        load_heavy = os.environ.get('OMNIPATH_BUILD_FTP_HEAVY_TYPES') == '1'
        with Session(self.engine) as session:
            backend = (
                session.query(Backend).filter_by(name='uniprot_ftp').first()
            )
            if not backend:
                _log.error('Backend uniprot_ftp not found in database')
                return
            backend_id = backend.id
            uniprot_type = (
                session.query(IdType).filter_by(name='uniprot').first()
            )
            if not uniprot_type:
                _log.error('ID type "uniprot" not found')
                return
            uniprot_type_id = uniprot_type.id

            label_to_id = {}
            for ftp_label, canonical_name in IDTYPE_MAP.items():
                if canonical_name.startswith('_'):  # _taxid handled separately
                    continue
                if ftp_label in self._FTP_HEAVY_LABELS and not load_heavy:
                    continue
                id_type = (
                    session.query(IdType).filter_by(name=canonical_name).first()
                )
                if id_type:
                    label_to_id[ftp_label] = id_type.id
        _log.info(
            'FTP label map: %d ID-type labels (heavy types %s)',
            len(label_to_id),
            'included' if load_heavy else 'excluded',
        )

        stg = 'idmapping_staging'
        tax = 'idmapping_ac_taxid'
        lmap = 'idmapping_label_map'
        new = f'{SCHEMA}.id_mapping_ftp_new'

        conn = get_connection(self._db_url)
        conn.autocommit = False
        try:
            with conn.cursor() as cur:
                self._pg_set(cur)

                # 1. Stream decompressed file -> raw COPY into UNLOGGED staging.
                cur.execute(f'DROP TABLE IF EXISTS {stg}')
                cur.execute(
                    f'CREATE UNLOGGED TABLE {stg} '
                    '(ac text, id_type_label text, id_value text)'
                )
                _log.info('Streaming idmapping into staging via COPY...')
                # OMNIPATH_BUILD_FTP_FILE: stream a pre-staged .gz directly
                # (skips the download-manager lookup).
                ftp_file = os.environ.get('OMNIPATH_BUILD_FTP_FILE')
                with cur.copy(
                    f'COPY {stg} (ac, id_type_label, id_value) FROM STDIN'
                ) as copy:
                    for block in stream_full_idmapping(path=ftp_file):
                        copy.write(block)
                conn.commit()
                cur.execute(f'SELECT count(*) FROM {stg}')
                staged = cur.fetchone()[0]
                _log.info('Staged %d raw idmapping rows', staged)

                # 2a. label -> id_type_id map table (for the FK join).
                cur.execute(f'DROP TABLE IF EXISTS {lmap}')
                cur.execute(
                    f'CREATE UNLOGGED TABLE {lmap} '
                    '(id_type_label text PRIMARY KEY, id_type_id smallint)'
                )
                with cur.copy(
                    f'COPY {lmap} (id_type_label, id_type_id) FROM STDIN'
                ) as copy:
                    for label, tid in label_to_id.items():
                        copy.write_row((label, tid))

                # 2b. ac -> taxid from the NCBI_TaxID rows.
                cur.execute(f'DROP TABLE IF EXISTS {tax}')
                cur.execute(
                    f'CREATE UNLOGGED TABLE {tax} AS '
                    f"SELECT ac, id_value::int AS taxid FROM {stg} "
                    "WHERE id_type_label = 'NCBI_TaxID' AND id_value ~ '^[0-9]+$'"
                )
                cur.execute(f'CREATE INDEX ON {tax} (ac)')
                cur.execute(f'ANALYZE {tax}')
                conn.commit()

                # 2c. Build the new table with id_type + taxon resolved as FKs.
                limit = (
                    f' LIMIT {int(self._max_records)}'
                    if self._max_records else ''
                )
                cur.execute(f'DROP TABLE IF EXISTS {new}')
                cur.execute(
                    # UNLOGGED: skip WAL on this large rebuildable build output.
                    f'CREATE UNLOGGED TABLE {new} AS '
                    f'SELECT {uniprot_type_id}::smallint AS source_type_id, '
                    'lm.id_type_id::smallint AS target_type_id, '
                    'COALESCE(t.taxid, 0)::integer AS ncbi_tax_id, '
                    'left(s.ac, 64)::varchar(64) AS source_id, '
                    # Normalise Ensembl IDs to versionless (ADR 0006): the FTP
                    # file carries `.N` version suffixes (ENSG/ENST/ENSP) but
                    # resources supply versionless IDs, so a versioned target
                    # would never join. Strip the trailing version here so the
                    # delivered mapping table is already canonical.
                    "left(CASE WHEN s.id_type_label IN "
                    "('Ensembl', 'Ensembl_TRS', 'Ensembl_PRO') "
                    "THEN regexp_replace(s.id_value, '\\.[0-9]+$', '') "
                    "ELSE s.id_value END, 64)::varchar(64) AS target_id, "
                    f'{backend_id}::smallint AS backend_id '
                    f'FROM {stg} s '
                    f'JOIN {lmap} lm USING (id_type_label) '
                    f'LEFT JOIN {tax} t USING (ac)'
                    f'{limit}'
                )
                conn.commit()
                cur.execute(f'SELECT count(*) FROM {new}')
                mapping_count = cur.fetchone()[0]
                _log.info('Transformed %d FTP mapping rows', mapping_count)

                # 3a. Build covering indexes offline (mirror id_mapping's
                # secondaries, + INCLUDE the payload so lookups are index-only
                # scans on this billion-row table). Parallel (unlogged -> no WAL).
                cur.execute(
                    f'CREATE INDEX ON {new} '
                    '(source_type_id, target_type_id, ncbi_tax_id, source_id) '
                    'INCLUDE (target_id)'
                )
                cur.execute(
                    f'CREATE INDEX ON {new} '
                    '(target_type_id, source_type_id, ncbi_tax_id, target_id) '
                    'INCLUDE (source_id)'
                )
                cur.execute(f'ANALYZE {new}')
                conn.commit()

                # 3b. Atomic swap behind the id_mapping_all view. CASCADE drops
                # the resolver views (resolver_gene / resolver_protein /
                # resolver_gene_protein_global) that depend on id_mapping_ftp;
                # create_resolver_views() below rebuilds them all from
                # resolver_protein.sql.
                cur.execute(f'DROP VIEW IF EXISTS {SCHEMA}.id_mapping_all CASCADE')
                cur.execute(f'DROP TABLE IF EXISTS {SCHEMA}.id_mapping_ftp CASCADE')
                cur.execute(
                    f'ALTER TABLE {new} RENAME TO id_mapping_ftp'
                )
                cols = (
                    'source_type_id, target_type_id, ncbi_tax_id, '
                    'source_id, target_id, backend_id'
                )
                cur.execute(
                    f'CREATE VIEW {SCHEMA}.id_mapping_all AS '
                    f'SELECT {cols} FROM {SCHEMA}.id_mapping '
                    f'UNION ALL SELECT {cols} FROM {SCHEMA}.id_mapping_ftp'
                )

                # 4. Cleanup staging.
                for t in (stg, tax, lmap):
                    cur.execute(f'DROP TABLE IF EXISTS {t}')
                conn.commit()
                # 5. Crash-safety: a crash TRUNCATEs an UNLOGGED table on
                # recovery — the 744M map was lost to an OOM reboot once. The
                # load + index builds above stay unlogged (fast, no per-row WAL);
                # this one-time SET LOGGED (own txn, after the swap commits)
                # writes the table to WAL so it survives crashes. Reads are
                # unaffected. Opt out for a throwaway build: set
                # OMNIPATH_BUILD_FTP_UNLOGGED=1.
                if os.environ.get('OMNIPATH_BUILD_FTP_UNLOGGED', '') not in (
                    '1', 'true', 'True', 'yes',
                ):
                    _log.info('SET LOGGED on id_mapping_ftp (crash-safety)')
                    cur.execute(
                        f'ALTER TABLE {SCHEMA}.id_mapping_ftp SET LOGGED'
                    )
                    conn.commit()
        except Exception:
            conn.rollback()
            conn.close()
            raise
        conn.close()

        # Secondary->primary UniProt normalization (ADR 0006) must exist before
        # the resolver views, which read it to canonicalise every UniProt.
        try:
            self.load_uniprot_sec_ac()
        except Exception as e:
            _log.error('Failed to load sec_ac (secondary->primary UniProt): %s', e)

        # Authoritative ensp/ensg->entrez (NCBI gene2ensembl) — resolver_gene
        # anchors Ensembl/gene ids through this gene-space map, not UniProt.
        try:
            self.load_gene2ensembl()
        except Exception as e:
            _log.error('Failed to load gene2ensembl: %s', e)

        # The resolver projection reads id_mapping_ftp, so (re)create it now.
        # id_mapping_ftp was just reloaded (its CASCADE drop already removed the
        # FTP core), so force the expensive resolver_gene_ftp rebuild (T065).
        self.create_resolver_views(force_ftp_core=True)

        duration = time.time() - start
        _log.info(
            'Full FTP build complete: %d mappings into id_mapping_ftp, %.1f min',
            mapping_count,
            duration / 60,
        )
        with Session(self.engine) as session:
            session.add(
                BuildInfo(
                    table_name='id_mapping_ftp',
                    source_type='uniprot (all)',
                    target_type='all FTP types',
                    ncbi_tax_id=0,
                    backend='uniprot_ftp',
                    row_count=mapping_count,
                    duration_secs=duration,
                    status='done',
                )
            )
            session.commit()

    def export_parquet(self, tables: list[tuple], output_dir: str):
        """Export mapping tables as Parquet files for fast API delivery."""
        import pyarrow as pa
        import pyarrow.parquet as pq

        os.makedirs(output_dir, exist_ok=True)

        from omnipath_utils.db._query import get_full_table

        for src_type, tgt_type, ncbi_tax_id in tables:
            with Session(self.engine) as session:
                data = get_full_table(
                    session, src_type, tgt_type, ncbi_tax_id
                )

            if not data:
                _log.warning(
                    'No data for %s -> %s (org %d), skipping Parquet',
                    src_type,
                    tgt_type,
                    ncbi_tax_id,
                )
                continue

            # Convert to flat rows for Parquet
            rows = [
                (k, v) for k, targets in data.items() for v in targets
            ]
            if not rows:
                continue

            table = pa.table(
                {
                    src_type: [r[0] for r in rows],
                    tgt_type: [r[1] for r in rows],
                }
            )

            fname = f'{src_type}__{tgt_type}__{ncbi_tax_id}.parquet'
            path = os.path.join(output_dir, fname)
            pq.write_table(table, path, compression='snappy')
            _log.info('Exported %s: %d rows', fname, len(rows))

    def populate_metabolites(self):
        """Populate metabolite ID mappings from UniChem, RaMP, MetaNetX and BiGG.

        Auto-discovers all available ID types from all sources and builds
        all available pairwise mappings.
        """
        self._populate_unichem()
        self._populate_ramp()
        self._populate_metanetx()
        self._populate_bigg()
        self._populate_reactome()
        self._clear_wrong_content_rows()
        self._populate_structures()
        self._populate_chemical_long()
        self._populate_lipid_names()
        self.populate_namespace_exemptions()
        self.record_structure_key_capability()
        self.record_lipid_nomenclature_capability()

    # ------------------------------------------------------------------
    # Long-value chemical mappings: names + structures (id_mapping_long)
    # ------------------------------------------------------------------

    #: Values longer than this (bytes) are skipped, not truncated, to stay
    #: within the B-tree index-key bound and never forge a wrong key
    #: (FR-017/FR-019 edge case). Full names/IUPAC and typical InChI/SMILES are
    #: far below this; only pathological structure strings exceed it.
    _MAX_LONG_KEY = 2000

    @staticmethod
    def _long_rows(
        data: dict,
        src_type_id: int,
        tgt_type_id: int,
        backend_id: int,
        is_name: bool,
        max_key: int = _MAX_LONG_KEY,
        limit: int | None = None,
    ) -> tuple[list[tuple], int]:
        """Pure projection of a ``{source: {targets}}`` mapping into
        ``id_mapping_long`` COPY rows ``(src_type_id, tgt_type_id, 0, source_id,
        source_label, target_id, backend_id)``.

        Name source keys are lowercased + whitespace-stripped with the original
        kept as ``source_label`` (FR-002); structure/db-id keys are verbatim
        (``source_label`` is ``None``; FR-019). Empty and over-``max_key`` values
        are skipped (counted), never truncated. ``limit`` caps the row count
        (FR-012). Returns ``(rows, skipped)``.
        """
        rows: list[tuple] = []
        skipped = 0
        for source_id, target_ids in data.items():
            raw = str(source_id).strip()
            if is_name:
                from omnipath_utils.mapping._id_types import normalize_name

                key, label = normalize_name(raw), raw
            else:
                key, label = raw, None
            if not key or len(key.encode('utf-8')) > max_key:
                skipped += 1
                continue
            for target_id in target_ids:
                tgt = str(target_id).strip()
                if not tgt or len(tgt.encode('utf-8')) > max_key:
                    skipped += 1
                    continue
                rows.append(
                    (src_type_id, tgt_type_id, 0, key, label, tgt, backend_id)
                )
                if limit is not None and len(rows) >= limit:
                    return rows, skipped
        return rows, skipped

    def _long_type_class(self, type_name: str) -> str:
        """Classify an id_type for ``id_mapping_long`` case handling."""
        from omnipath_utils.db._query import NAME_TYPES, STRUCTURE_TYPES

        if type_name in NAME_TYPES:
            return 'name'
        if type_name in STRUCTURE_TYPES:
            return 'structure'
        return 'id'

    @staticmethod
    def _canonicalize_smiles_in_data(
        data: dict[str, set[str]], source_type: str, target_type: str,
    ) -> dict[str, set[str]]:
        """Replace a raw SMILES key/value with its canonical form (T079),
        so the stored lookup key matches what ``_lookup_key`` computes for
        an incoming query (T080) -- two spellings of one molecule become
        one stored row instead of two unmatched ones.

        Runs *before* ``_long_rows``, which stays byte-verbatim for
        structures (its own deliberate, tested contract -- naive folding
        would corrupt SMILES case-sensitivity). Chemistry-aware
        canonicalization is a different transform at a different layer,
        not a relaxation of that contract. A value that fails to
        canonicalize (garbage, or the toolkit unavailable) keeps its raw
        form -- degrades to a byte comparison rather than dropping the row
        (research R14).
        """
        if source_type != 'smiles' and target_type != 'smiles':
            return data

        from omnipath_utils.mapping._chemistry import (
            canonicalize_smiles,
            chemistry_available,
        )

        if not chemistry_available():
            return data

        canon_src = source_type == 'smiles'
        canon_tgt = target_type == 'smiles'
        out: dict[str, set[str]] = {}
        for src, targets in data.items():
            new_src = (canonicalize_smiles(src) or src) if canon_src else src
            new_targets = {
                (canonicalize_smiles(t) or t) if canon_tgt else t
                for t in targets
            }
            out.setdefault(new_src, set()).update(new_targets)
        return out

    def _populate_long_slice(
        self,
        data: dict,
        source_type: str,
        target_type: str,
        backend_name: str,
    ) -> int:
        """Load one ``(source_type, target_type, backend)`` slice into
        ``id_mapping_long`` -- the long-value (name/structure) sibling table.

        Name source keys are lowercased + whitespace-stripped (case-insensitive
        matching, FR-002) with the original case kept in ``source_label``;
        structure source keys (``inchi``/``smiles``) are stored verbatim
        (FR-019). Idempotent: the slice is DELETEd then COPYd (R6). Honours the
        backend's effective ``--max-records`` cap (FR-012). Over-long values are
        skipped and counted, never truncated. Returns the row count written.
        """
        if not data:
            return 0

        with Session(self.engine) as session:
            src_type = session.query(IdType).filter_by(name=source_type).first()
            tgt_type = session.query(IdType).filter_by(name=target_type).first()
            backend = session.query(Backend).filter_by(name=backend_name).first()
            if not src_type or not tgt_type or not backend:
                _log.error(
                    'id_mapping_long: missing id_type/backend for %s -> %s (%s)',
                    source_type, target_type, backend_name,
                )
                return 0
            src_type_id, tgt_type_id, backend_id = (
                src_type.id, tgt_type.id, backend.id,
            )

        is_name = self._long_type_class(source_type) == 'name'
        limit = self._effective_limit(backend_name)
        start = time.time()

        # Idempotent: replace this slice (organism-agnostic, tax 0).
        with Session(self.engine) as session:
            session.execute(
                text(
                    f'DELETE FROM {SCHEMA}.id_mapping_long'
                    ' WHERE source_type_id = :src AND target_type_id = :tgt'
                    ' AND ncbi_tax_id = 0 AND backend_id = :bk'
                ),
                {'src': src_type_id, 'tgt': tgt_type_id, 'bk': backend_id},
            )
            session.commit()

        data = self._canonicalize_smiles_in_data(data, source_type, target_type)
        rows, skipped = self._long_rows(
            data, src_type_id, tgt_type_id, backend_id, is_name,
            self._MAX_LONG_KEY, limit,
        )

        from omnipath_utils.db._connection import get_connection

        conn = get_connection(self._db_url)
        row_count = 0
        try:
            with conn.cursor() as cur:
                with cur.copy(
                    f'COPY {SCHEMA}.id_mapping_long'
                    ' (source_type_id, target_type_id, ncbi_tax_id,'
                    ' source_id, source_label, target_id, backend_id)'
                    ' FROM STDIN'
                ) as copy:
                    for row in rows:
                        copy.write_row(row)
                        row_count += 1
            conn.commit()
        finally:
            conn.close()

        duration = time.time() - start
        with Session(self.engine) as session:
            session.add(
                BuildInfo(
                    table_name='id_mapping_long',
                    source_type=source_type,
                    target_type=target_type,
                    ncbi_tax_id=0,
                    backend=backend_name,
                    row_count=row_count,
                    duration_secs=duration,
                    status='done',
                )
            )
            session.commit()

        _log.info(
            'id_mapping_long %s -> %s (%s): %d rows%s in %.1fs',
            source_type, target_type, backend_name, row_count,
            f' ({skipped} skipped)' if skipped else '', duration,
        )
        return row_count

    #: Builders whose failure fails the whole build, not just a warning.
    #: HMDB carries the richest per-metabolite cross-reference block in the
    #: field (spec 011 T011). Silently continuing without it would repeat
    #: the gap cycle 005 already left unmeasured for a full cycle.
    _REQUIRED_LONG_BUILDERS = {'hmdb'}

    def _populate_chemical_long(self):
        """Build the long-value chemical layer in ``id_mapping_long``:
        names (name/synonym/iupac/traditional_iupac) and structures
        (inchi/smiles) for every targeted resource, with the full pairwise set
        across the three identifier classes (US1/US2).

        Per-resource builders are dispatched here; each writes its slices via
        :meth:`_populate_long_slice` and its counts to ``build_info``. The
        dispatcher records per-class roll-ups at the end (FR-011).
        """
        _log.info('Building long-value chemical mappings (names + structures)...')

        builders = [
            ('chebi', self._long_chebi),
            ('hmdb', self._long_hmdb),
            ('chembl', self._long_chembl),
            ('kegg_compound', self._long_kegg),
            ('ramp', self._long_ramp),
            ('refmet', self._long_refmet),
        ]
        for label, fn in builders:
            try:
                fn()
            except Exception as e:
                if label in self._REQUIRED_LONG_BUILDERS:
                    raise
                _log.warning('chemical_long: %s builder failed: %s', label, e)

        self._record_long_rollups()

    def _emit_long_relations(
        self,
        rows,
        backend: str,
        name_cols: dict,
        id_cols: dict,
        struct_cols: dict | None = None,
    ) -> int:
        """Materialise the full cross-class relation set for one resource (R3).

        From each raw row's name column(s) (``name_cols``: long-type -> column),
        database-ID column(s) (``id_cols``; the first is the resource's own
        primary id) and structure column(s) (``struct_cols``: inchi/smiles ->
        column), emit every derivable pair across the three identifier classes
        into ``id_mapping_long`` -- names<->ids, names<->structures,
        id<->structures, structure<->structure -- not hub-only (FR-004/FR-019).
        Database-ID<->database-ID pairs are intentionally excluded (they belong
        in ``id_mapping`` via UniChem). Loads each slice via the idempotent
        long-value COPY helper and returns the total row count.
        """
        from collections import defaultdict

        from omnipath_utils.mapping.backends._inputs_v2_adapter import _as_values

        struct_cols = struct_cols or {}
        primary = next(iter(id_cols), None)
        struct_types = list(struct_cols)
        pairs: dict = defaultdict(lambda: defaultdict(set))

        for row in rows:
            names = {nt: _as_values(row.get(c)) for nt, c in name_cols.items()}
            ids = {it: _as_values(row.get(c)) for it, c in id_cols.items()}
            structs = {st: _as_values(row.get(c)) for st, c in struct_cols.items()}

            # names -> every database id and structure (forward)
            for nt, nvals in names.items():
                for nv in nvals:
                    for it, ivals in ids.items():
                        if ivals:
                            pairs[(nt, it)][nv].update(ivals)
                    for st, svals in structs.items():
                        if svals:
                            pairs[(nt, st)][nv].update(svals)

            # reverse from the resource's primary id -> names and structures
            for pv in ids.get(primary, []):
                for nt, nvals in names.items():
                    if nvals:
                        pairs[(primary, nt)][pv].update(nvals)
                for st, svals in structs.items():
                    if svals:
                        pairs[(primary, st)][pv].update(svals)
                        for sv in svals:
                            pairs[(st, primary)][sv].add(pv)

            # structure <-> structure (inchi <-> smiles)
            if len(struct_types) == 2:
                a, b = struct_types
                for av in structs[a]:
                    if structs[b]:
                        pairs[(a, b)][av].update(structs[b])
                for bv in structs[b]:
                    if structs[a]:
                        pairs[(b, a)][bv].update(structs[a])

        total = 0
        for (st, tt), data in pairs.items():
            data = {k: v for k, v in data.items() if v}
            if data:
                total += self._populate_long_slice(data, st, tt, backend)
        return total

    # Per-resource long-value builders (US1 T015/T016, US2 T024/T025/T026/T027a).
    def _long_chebi(self):
        """ChEBI names/synonyms <-> chebi + xref ids + structures. US1 (T015)."""
        from omnipath_utils.mapping.backends._inputs_v2_adapter import raw_rows

        rows = raw_rows('chebi', 'molecules', self._effective_limit('chebi'))
        n = self._emit_long_relations(
            rows, 'chebi',
            name_cols={'name': 'name', 'synonym': 'synonyms'},
            id_cols={
                'chebi': 'chebi_id', 'hmdb': 'hmdb', 'kegg': 'kegg_compound',
                'pubchem': 'pubchem_compound', 'lipidmaps': 'lipidmaps',
            },
            struct_cols={'inchi': 'inchi', 'smiles': 'smiles'},
        )
        _log.info('chemical_long chebi: %d rows', n)

    def _long_hmdb(self):
        """HMDB names/synonyms <-> hmdb + xref ids + structures. US1 (T016).

        HMDB carries the richest per-metabolite cross-reference block in the
        field: ChEBI, KEGG, PubChem, DrugBank, FooDB and CAS, each on the
        same record, plus its own IUPAC name and traditional IUPAC name. The
        project mirrors the archive on its own host
        (``pypath.inputs_v2.hmdb``). This reads it directly, unlike the
        legacy ``pypath.inputs.hmdb`` path, which reached only a narrow
        synonym-to-ChEBI slice and needed HMDB's site to stay up.

        T087 (spec 011): ``id_types.yaml`` also declares HMDB a backend for
        knapsack/chemspider/metlin/vmh/wikidata/biocyc/fbonto/
        phenol_explorer/molport/wikipedia -- verified against a live raw
        row that none of those fields actually exist in HMDB's parsed
        schema (confirmed columns: accession, cas_registry_number,
        chebi_id, drugbank_id, foodb_id, inchi, inchikey, iupac_name,
        kegg_id, name, pubchem_compound_id, smiles, synonyms,
        traditional_iupac -- nothing else namespace-shaped). Those ten
        namespaces are declared but HMDB does not, in fact, reach them --
        recorded as exemptions (T093), not silently dropped.
        """
        from omnipath_utils.mapping.backends._inputs_v2_adapter import raw_rows

        rows = raw_rows('hmdb', 'metabolites', self._effective_limit('hmdb'))
        for row in rows:
            synonyms = row.get('synonyms')
            if isinstance(synonyms, str):
                row['synonyms'] = [s for s in synonyms.split(';') if s]
        n = self._emit_long_relations(
            rows, 'hmdb',
            name_cols={
                'name': 'name', 'synonym': 'synonyms',
                'iupac': 'iupac_name', 'traditional_iupac': 'traditional_iupac',
            },
            id_cols={
                'hmdb': 'accession', 'chebi': 'chebi_id', 'kegg': 'kegg_id',
                'pubchem': 'pubchem_compound_id', 'drugbank': 'drugbank_id',
                'foodb': 'foodb_id', 'cas': 'cas_registry_number',
            },
            struct_cols={'inchi': 'inchi', 'smiles': 'smiles'},
        )
        _log.info('chemical_long hmdb: %d rows', n)

    def _long_chembl(self):
        """ChEMBL name/synonym <-> chembl/chebi + structures. US2 (T024)."""
        from omnipath_utils.mapping.backends._inputs_v2_adapter import raw_rows

        rows = raw_rows('chembl', 'molecules', self._effective_limit('chembl'))
        name_cols = {'name': 'pref_name'}
        if rows and 'synonyms' in rows[0]:
            name_cols['synonym'] = 'synonyms'
        n = self._emit_long_relations(
            rows, 'chembl',
            name_cols=name_cols,
            id_cols={'chembl': 'chembl_id', 'chebi': 'chebi_id'},
            struct_cols={'inchi': 'standard_inchi', 'smiles': 'canonical_smiles'},
        )
        _log.info('chemical_long chembl: %d rows', n)

    def _long_refmet(self):
        """RefMet names <-> refmet/chebi/hmdb/kegg/lipidmaps/pubchem ids.

        US5 (T088). 169,859 chemical entities canonicalize on a RefMet name
        with no RefMet row in the translation database at all before this
        (research.md R7). RefMet carries no full InChI/SMILES of its own
        (only an InChIKey, which belongs in id_mapping's short-key table via
        a structure backend, not here -- out of this task's scope).
        """
        from omnipath_utils.mapping.backends._inputs_v2_adapter import raw_rows

        rows = raw_rows('refmet', 'metabolites', self._effective_limit('refmet'))
        n = self._emit_long_relations(
            rows, 'refmet',
            name_cols={'name': 'refmet_name'},
            id_cols={
                # The raw RefMet CSV's first header carries a leading
                # space (" refmet_id", not "refmet_id") -- verified
                # against a live fetch, not a typo. Every other column
                # name is clean.
                'refmet': ' refmet_id', 'chebi': 'chebi_id',
                'hmdb': 'hmdb_id', 'kegg': 'kegg_id',
                'lipidmaps': 'lipidmaps_id', 'pubchem': 'pubchem_cid',
            },
        )
        _log.info('chemical_long refmet: %d rows', n)

    def _long_kegg(self):
        """KEGG compound names <-> kegg/chebi. US2 (T025)."""
        from pypath.inputs.kegg import kegg_compound_chebi, kegg_compound_names

        names = kegg_compound_names()      # {kegg_id: [names]}
        kc = kegg_compound_chebi()         # {kegg_id: 'CHEBI:nnn'}
        name_kegg: dict = {}
        name_chebi: dict = {}
        kegg_name: dict = {}
        for kid, nlist in names.items():
            nlist = nlist if isinstance(nlist, (list, set, tuple)) else [nlist]
            kegg_name.setdefault(kid, set()).update(str(n) for n in nlist if n)
            for nm in nlist:
                if not nm:
                    continue
                name_kegg.setdefault(str(nm), set()).add(kid)
                if kid in kc:
                    name_chebi.setdefault(str(nm), set()).add(kc[kid])
        total = 0
        total += self._populate_long_slice(name_kegg, 'name', 'kegg', 'kegg_compound')
        total += self._populate_long_slice(kegg_name, 'kegg', 'name', 'kegg_compound')
        total += self._populate_long_slice(name_chebi, 'name', 'chebi', 'kegg_compound')
        _log.info('chemical_long kegg_compound: %d rows', total)

    @staticmethod
    def _norm_chem_id(id_type: str, value: str) -> str:
        """Normalise a RaMP curie to the DB's id format (e.g. CHEBI:15377)."""
        v = str(value)
        if ':' in v:
            v = v.split(':', 1)[1]
        if id_type == 'chebi':
            return f'CHEBI:{v}'
        return v

    def _long_ramp(self):
        """RaMP synonyms <-> chebi/hmdb (analytesynonym). US2 (T026).

        ``ramp_synonym_mapping(hub)`` returns ``{hub_id: {synonyms}}``; emit both
        directions (synonym->hub and hub->synonym), normalising the hub id to the
        DB format.
        """
        from pypath.inputs.ramp import ramp_synonym_mapping

        total = 0
        for hub in ('chebi', 'hmdb'):
            try:
                data = ramp_synonym_mapping(hub, curies=True)  # {hub_id:{syn}}
            except Exception as e:
                _log.debug('RaMP synonym<->%s failed: %s', hub, e)
                continue
            syn_hub: dict = {}
            hub_syn: dict = {}
            for hid, syns in data.items():
                hid_n = self._norm_chem_id(hub, hid)
                clean = {str(s) for s in syns if s}
                hub_syn.setdefault(hid_n, set()).update(clean)
                for s in clean:
                    syn_hub.setdefault(s, set()).add(hid_n)
            total += self._populate_long_slice(syn_hub, 'synonym', hub, 'ramp')
            total += self._populate_long_slice(hub_syn, hub, 'synonym', 'ramp')
        _log.info('chemical_long ramp: %d synonym rows', total)

    def _record_long_rollups(self):
        """Write per-class roll-up counts for the long-value layer to
        ``build_info`` (``names.total`` / ``structures.total``) from the slices
        already recorded this build."""
        from omnipath_utils.db._query import NAME_TYPES, STRUCTURE_TYPES

        with Session(self.engine) as session:
            rows = session.execute(
                text(
                    f"SELECT source_type, target_type, row_count"
                    f" FROM {SCHEMA}.build_info"
                    " WHERE table_name = 'id_mapping_long'"
                    " AND status = 'done'"
                )
            ).fetchall()
            names_total = sum(
                r[2] or 0 for r in rows
                if r[0] in NAME_TYPES or r[1] in NAME_TYPES
            )
            struct_total = sum(
                r[2] or 0 for r in rows
                if r[0] in STRUCTURE_TYPES or r[1] in STRUCTURE_TYPES
            )
            for kind, total in (
                ('names', names_total), ('structures', struct_total),
            ):
                session.add(
                    BuildInfo(
                        table_name='id_mapping_long',
                        source_type=kind,
                        target_type='total',
                        ncbi_tax_id=0,
                        backend='rollup',
                        row_count=total,
                        status='done',
                    )
                )
            session.commit()
        _log.info(
            'id_mapping_long roll-up: %d name rows, %d structure rows',
            names_total, struct_total,
        )

    def record_ftp_types(self):
        """Precompute and cache the set of id_types present in
        ``id_mapping_ftp`` into ``build_info`` (``table_name = 'ftp_types'``).

        Read at query time by the FR-018 fallback gate so it never has to scan
        the 744 M-row table. Safe/idempotent: replaces any prior row. A no-op if
        the FTP table is absent.
        """
        with Session(self.engine) as session:
            present = session.execute(
                text(f"SELECT to_regclass('{SCHEMA}.id_mapping_ftp')")
            ).scalar()
            if present is None:
                _log.info('record_ftp_types: id_mapping_ftp absent, skipping')
                return
            names = session.execute(
                text(
                    f'SELECT name FROM {SCHEMA}.id_type WHERE id IN ('
                    f' SELECT DISTINCT source_type_id FROM {SCHEMA}.id_mapping_ftp'
                    f' UNION'
                    f' SELECT DISTINCT target_type_id FROM {SCHEMA}.id_mapping_ftp)'
                    ' ORDER BY name'
                )
            ).scalars().all()
            session.execute(
                text(
                    f"DELETE FROM {SCHEMA}.build_info WHERE table_name='ftp_types'"
                )
            )
            # One row per type name (build_info.source_type is varchar(64)).
            for name in names:
                session.add(
                    BuildInfo(
                        table_name='ftp_types',
                        source_type=name,
                        target_type=None,
                        ncbi_tax_id=0,
                        backend='ftp_types',
                        row_count=len(names),
                        status='done',
                    )
                )
            session.commit()
        _log.info('record_ftp_types: cached %d FTP id_types', len(names))

    #: Name-axis id_type names whose values may carry a lipid shorthand.
    _LIPID_NAME_SOURCE_TYPES = ('name', 'synonym', 'iupac', 'traditional_iupac')

    def _populate_lipid_names(self):
        """Build ``lipid_name`` (spec 011 T113): every distinct name-axis
        value that looks lipid-shaped (has a ``C:D`` chain token), parsed by
        the nomenclature grammar into its standardized identity fields.

        Parallel across cores -- pygoslin parsing is pure-Python and
        CPU-bound (unlike rdkit's C extension, it does not release the GIL),
        so threads would not help; a process pool does, mirroring
        omnipath-metabo's ``_lipid_layer.py._parse_many`` (the same
        candidate-selection + parallel-parse shape, independently applied
        here because this is a different database with its own name
        corpus). Idempotent: replaces the whole table each run.
        """
        from omnipath_utils.mapping._lipid import (
            lipid_nomenclature_available,
            parse_lipid_name,
        )

        if not lipid_nomenclature_available():
            _log.info('lipid_name: pygoslin unavailable, skipping (T114)')
            return 0

        with self.engine.connect() as conn:
            names = [
                row[0]
                for row in conn.execute(
                    text(
                        f"""
                        SELECT DISTINCT coalesce(m.source_label, m.source_id)
                        FROM {SCHEMA}.id_mapping_long m
                        JOIN {SCHEMA}.id_type t ON t.id = m.source_type_id
                        WHERE t.name = ANY(:types)
                          AND coalesce(m.source_label, m.source_id) ~ '\\d+:\\d+'
                        """
                    ),
                    {'types': list(self._LIPID_NAME_SOURCE_TYPES)},
                )
            ]

        _log.info('lipid_name: %d candidate names to parse', len(names))
        if not names:
            return 0

        workers = min(16, max(1, (os.cpu_count() or 2) - 2))
        if workers > 1 and len(names) >= 2000:
            from multiprocessing import Pool

            with Pool(processes=workers) as pool:
                parsed = list(pool.imap(parse_lipid_name, names, chunksize=1000))
        else:
            parsed = [parse_lipid_name(n) for n in names]

        rows = [
            {
                'raw_name': raw,
                'lipid_name': p['lipid_name'],
                'lipid_level': p['lipid_level'],
                'chains_possible': p['chains_possible'],
                'chains_listed': p['chains_listed'],
                'lipid_category': p['lipid_category'],
                'lipid_class': p['lipid_class'],
                'total_carbon': p['total_carbon'],
                'total_db': p['total_db'],
                'sum_formula': p['sum_formula'],
                'parser_version': p['parser_version'],
            }
            for raw, p in zip(names, parsed)
            if p is not None
        ]
        _log.info(
            'lipid_name: %d/%d names parsed successfully', len(rows), len(names)
        )

        with Session(self.engine) as session:
            session.execute(text(f'TRUNCATE {SCHEMA}.lipid_name'))
            if rows:
                session.execute(LipidName.__table__.insert(), rows)
            session.commit()
        return len(rows)

    def record_lipid_nomenclature_capability(self):
        """Record whether this build ran with the lipid nomenclature grammar
        (spec 011 T114) into ``build_info`` -- mirrors
        :meth:`record_structure_key_capability` exactly, one row,
        ``table_name='capability'``, a different ``source_type``.
        """
        from omnipath_utils.mapping._lipid import lipid_nomenclature_available

        available = lipid_nomenclature_available()
        with Session(self.engine) as session:
            session.execute(
                text(
                    f"DELETE FROM {SCHEMA}.build_info"
                    " WHERE table_name='capability'"
                    " AND source_type='lipid_nomenclature'"
                )
            )
            session.add(
                BuildInfo(
                    table_name='capability',
                    source_type='lipid_nomenclature',
                    target_type=None,
                    ncbi_tax_id=0,
                    backend='pygoslin',
                    row_count=None,
                    status='available' if available else 'unavailable',
                )
            )
            session.commit()
        _log.info(
            'lipid_nomenclature capability: %s',
            'available' if available else 'unavailable',
        )

    def record_structure_key_capability(self):
        """Record whether this build ran with the chemistry toolkit
        (spec 011 T082/R14) into ``build_info``
        (``table_name = 'capability'``). Read by omnipath-build's own
        ``build_capability`` detection (a query across databases -- this
        build's own capability, recorded once, is the real signal that
        finally replaces its previous hardcoded 'unavailable').

        Idempotent: replaces any prior row.
        """
        from omnipath_utils.mapping._chemistry import chemistry_available

        available = chemistry_available()
        with Session(self.engine) as session:
            session.execute(
                text(
                    f"DELETE FROM {SCHEMA}.build_info"
                    " WHERE table_name='capability'"
                    " AND source_type='structure_key_computation'"
                )
            )
            session.add(
                BuildInfo(
                    table_name='capability',
                    source_type='structure_key_computation',
                    target_type=None,
                    ncbi_tax_id=0,
                    backend='chemistry',
                    row_count=None,
                    status='available' if available else 'unavailable',
                )
            )
            session.commit()
        _log.info(
            'structure_key_computation capability: %s',
            'available' if available else 'unavailable',
        )

    # ------------------------------------------------------------------
    # UniChem auto-discovery
    # ------------------------------------------------------------------

    # Mapping from UniChem normalised labels to our canonical id_type
    # names where they differ.
    _UNICHEM_NAME_MAP: dict[str, str] = {
        'lipid_maps': 'lipidmaps',
        'probes&drugs': 'probes_drugs',
    }

    # Namespaces UniChem exposes but whose "source_id" is not that
    # namespace's own accession (spec 011 T091): UniChem's Rhea entries
    # return the participant's ChEBI id (Rhea reactions cite ChEBI
    # compounds directly and mint no compound id of their own), and its
    # PDBe entries return a display label plus conformer type (e.g. "GR2 -
    # Ideal conformer") instead of the ligand code. Excluded from
    # auto-discovery rather than stored wrong; see deferred-items.md for
    # what a real loader for each would need.
    _UNICHEM_WRONG_CONTENT = frozenset({'rhea', 'pdbe'})

    def _unichem_canonical(self, label: str) -> str | None:
        """Normalise a UniChem source label to a canonical id_type name."""
        import re

        raw = label.lower().replace(' ', '_')
        raw = re.sub(r'[^a-z0-9_&]', '', raw).strip('_')

        if not raw:
            return None

        canonical = self._UNICHEM_NAME_MAP.get(raw, raw)
        return None if canonical in self._UNICHEM_WRONG_CONTENT else canonical

    def _populate_unichem(self):
        """Auto-discover and build all UniChem pairwise mappings."""
        try:
            from pypath.inputs.unichem import (
                unichem_mapping,
                unichem_sources,
            )
        except ImportError:
            _log.warning('pypath not available for UniChem')
            return

        sources = unichem_sources()
        _log.info('UniChem: %d sources available', len(sources))

        # Register any newly discovered sources as ID types
        self._register_unichem_types(sources)

        # Map UniChem source numbers to our canonical id_type names
        num_to_name: dict[int, str] = {}

        for num, label in sources.items():
            canonical = self._unichem_canonical(label)

            if canonical:
                num_to_name[num] = canonical

        # Get backend ID
        with Session(self.engine) as session:
            backend = (
                session.query(Backend).filter_by(name='unichem').first()
            )

            if not backend:
                _log.error('Backend unichem not found')
                return

            backend_id = backend.id

        from omnipath_utils.db._connection import get_connection

        built = 0
        failed = 0
        source_nums = sorted(num_to_name.keys())

        for i, src_num in enumerate(source_nums):
            for tgt_num in source_nums[i + 1:]:
                src_name = num_to_name[src_num]
                tgt_name = num_to_name[tgt_num]

                try:
                    data = unichem_mapping(src_num, tgt_num)

                    if not data:
                        continue

                    with Session(self.engine) as session:
                        src_type = (
                            session.query(IdType)
                            .filter_by(name=src_name)
                            .first()
                        )
                        tgt_type = (
                            session.query(IdType)
                            .filter_by(name=tgt_name)
                            .first()
                        )

                        if not src_type or not tgt_type:
                            _log.debug(
                                'UniChem: id_type not found'
                                ' for %s or %s',
                                src_name,
                                tgt_name,
                            )
                            continue

                        src_type_id = src_type.id
                        tgt_type_id = tgt_type.id

                        # Idempotent re-run: clear this (src, tgt, tax=0,
                        # backend) slice before COPY so a cumulative rebuild
                        # does not append duplicate rows.
                        session.execute(
                            text(
                                f'DELETE FROM {SCHEMA}.id_mapping'
                                ' WHERE source_type_id = :src'
                                ' AND target_type_id = :tgt'
                                ' AND ncbi_tax_id = 0 AND backend_id = :bk'
                            ),
                            {
                                'src': src_type_id,
                                'tgt': tgt_type_id,
                                'bk': backend_id,
                            },
                        )
                        session.commit()

                    conn = get_connection(self._db_url)
                    row_count = 0

                    with conn.cursor() as cur:
                        with cur.copy(
                            f'COPY {SCHEMA}.id_mapping'
                            ' (source_type_id, target_type_id,'
                            ' ncbi_tax_id, source_id, target_id,'
                            ' backend_id) FROM STDIN'
                        ) as copy:
                            for src_id, tgt_ids in data.items():
                                for tgt_id in tgt_ids:
                                    copy.write_row((
                                        src_type_id,
                                        tgt_type_id,
                                        0,
                                        str(src_id)[:64],
                                        str(tgt_id)[:64],
                                        backend_id,
                                    ))
                                    row_count += 1
                                if (
                                    self._max_records is not None
                                    and row_count >= self._max_records
                                ):
                                    break

                    conn.commit()
                    conn.close()

                    if row_count:
                        built += 1
                        _log.info(
                            'UniChem %s -> %s: %d rows',
                            src_name,
                            tgt_name,
                            row_count,
                        )

                except Exception as e:
                    failed += 1
                    _log.debug(
                        'UniChem %s -> %s: %s',
                        src_name,
                        tgt_name,
                        e,
                    )

        _log.info(
            'UniChem: built %d tables, %d failed/empty',
            built,
            failed,
        )

    def _register_unichem_types(self, sources: dict[int, str]):
        """Register UniChem source IDs as id_types if not present."""
        with Session(self.engine) as session:
            for _num, label in sources.items():
                canonical = self._unichem_canonical(label)

                if not canonical:
                    continue

                existing = (
                    session.query(IdType)
                    .filter_by(name=canonical)
                    .first()
                )

                if not existing:
                    session.add(
                        IdType(
                            name=canonical,
                            label=label,
                            entity_type='small_molecule',
                        )
                    )
                    _log.info(
                        'Registered new ID type from UniChem:'
                        ' %s (%s)',
                        canonical,
                        label,
                    )

            session.commit()

    # ------------------------------------------------------------------
    # RaMP auto-discovery
    # ------------------------------------------------------------------

    # Mapping from RaMP IDtype strings to our canonical id_type names
    # where they differ.
    _RAMP_NAME_MAP: dict[str, str] = {
        'CAS': 'cas',
        'LIPIDMAPS': 'lipidmaps',
        'rhea-comp': 'rhea',
    }

    def _ramp_canonical(self, ramp_type: str) -> str | None:
        """Normalise a RaMP IDtype string to a canonical id_type name.

        ``None`` for a namespace RaMP's own ``rhea-comp`` type does not
        actually carry (spec 011 T091) -- it is the ChEBI id of the
        reaction participant, not a Rhea-native compound id, the same
        wrong-content issue as UniChem's Rhea source
        (``_UNICHEM_WRONG_CONTENT``).
        """
        canonical = self._RAMP_NAME_MAP.get(ramp_type, ramp_type.lower())
        return None if canonical in self._UNICHEM_WRONG_CONTENT else canonical

    def _populate_ramp(self):
        """Auto-discover and build RaMP pairwise mappings."""
        try:
            from pypath.inputs.ramp import ramp_mapping
            from pypath.inputs.ramp._sqlite import (
                id_types as ramp_id_types_fn,
            )
        except ImportError:
            _log.warning('pypath not available for RaMP')
            return

        # Discover available compound ID types
        try:
            ramp_types = sorted(
                ramp_id_types_fn(entity_type='compound')
            )
        except Exception as e:
            _log.error('RaMP discovery failed: %s', e)
            return

        _log.info(
            'RaMP: %d compound ID types available', len(ramp_types)
        )

        # Register any new types
        self._register_ramp_types(ramp_types)

        # Get backend ID
        with Session(self.engine) as session:
            backend = (
                session.query(Backend).filter_by(name='ramp').first()
            )

            if not backend:
                _log.error('Backend ramp not found')
                return

            backend_id = backend.id

        from omnipath_utils.db._connection import get_connection

        # Build important pairs: each type to chebi and hmdb (hub types)
        hubs = ['chebi', 'hmdb']
        built = 0

        for src_ramp_type in ramp_types:
            src_canonical = self._ramp_canonical(src_ramp_type)

            for hub in hubs:
                if src_canonical == hub:
                    continue

                try:
                    data = ramp_mapping(src_ramp_type, hub)

                    if not data:
                        continue

                    with Session(self.engine) as session:
                        src_t = (
                            session.query(IdType)
                            .filter_by(name=src_canonical)
                            .first()
                        )
                        tgt_t = (
                            session.query(IdType)
                            .filter_by(name=hub)
                            .first()
                        )

                        if not src_t or not tgt_t:
                            _log.debug(
                                'RaMP: id_type not found'
                                ' for %s or %s',
                                src_canonical,
                                hub,
                            )
                            continue

                        src_type_id = src_t.id
                        tgt_type_id = tgt_t.id

                        # Idempotent re-run: clear this (src, tgt, tax=0,
                        # backend) slice before COPY so a cumulative rebuild
                        # does not append duplicate rows.
                        session.execute(
                            text(
                                f'DELETE FROM {SCHEMA}.id_mapping'
                                ' WHERE source_type_id = :src'
                                ' AND target_type_id = :tgt'
                                ' AND ncbi_tax_id = 0 AND backend_id = :bk'
                            ),
                            {
                                'src': src_type_id,
                                'tgt': tgt_type_id,
                                'bk': backend_id,
                            },
                        )
                        session.commit()

                    conn = get_connection(self._db_url)
                    row_count = 0

                    with conn.cursor() as cur:
                        with cur.copy(
                            f'COPY {SCHEMA}.id_mapping'
                            ' (source_type_id, target_type_id,'
                            ' ncbi_tax_id, source_id, target_id,'
                            ' backend_id) FROM STDIN'
                        ) as copy:
                            for src_id, tgt_ids in data.items():
                                for tgt_id in tgt_ids:
                                    copy.write_row((
                                        src_type_id,
                                        tgt_type_id,
                                        0,
                                        str(src_id)[:64],
                                        str(tgt_id)[:64],
                                        backend_id,
                                    ))
                                    row_count += 1
                                if (
                                    self._max_records is not None
                                    and row_count >= self._max_records
                                ):
                                    break

                    conn.commit()
                    conn.close()

                    if row_count:
                        built += 1
                        _log.info(
                            'RaMP %s -> %s: %d rows',
                            src_canonical,
                            hub,
                            row_count,
                        )

                except Exception as e:
                    _log.debug(
                        'RaMP %s -> %s: %s',
                        src_canonical,
                        hub,
                        e,
                    )

        _log.info('RaMP: built %d tables', built)

    def _register_ramp_types(self, ramp_types: list[str]):
        """Register RaMP ID types if not already present."""
        with Session(self.engine) as session:
            for rtype in ramp_types:
                canonical = self._ramp_canonical(rtype)
                if canonical is None:
                    continue

                existing = (
                    session.query(IdType)
                    .filter_by(name=canonical)
                    .first()
                )

                if not existing:
                    session.add(
                        IdType(
                            name=canonical,
                            label=rtype,
                            entity_type='small_molecule',
                        )
                    )
                    _log.info(
                        'Registered new ID type from RaMP: %s',
                        canonical,
                    )

            session.commit()

    def populate_mirna(self, organisms: list[int]):
        """Populate miRNA ID mappings."""
        for org in organisms:
            for src, tgt in [
                ('mir-pre', 'mirbase'),
                ('mir-mat-name', 'mirbase'),
            ]:
                try:
                    self.populate_mapping(src, tgt, org, 'mirbase')
                except Exception as e:
                    _log.error(
                        'miRBase %s -> %s (org %d) failed: %s',
                        src,
                        tgt,
                        org,
                        e,
                    )

    def populate_orthology(self, organisms: list[int]):
        """Populate orthology tables for organism pairs."""
        source = 9606  # human as source

        from omnipath_utils.db._connection import get_connection
        from omnipath_utils.orthology._manager import OrthologyManager

        for target in organisms:
            if target == source:
                continue

            _log.info('Building orthology: %d -> %d', source, target)

            mgr = OrthologyManager()
            try:
                table = mgr._get_table(
                    source,
                    target,
                    'genesymbol',
                    resource='hcop',
                    min_sources=1,
                )
                if not table:
                    _log.warning(
                        'No orthology data for %d -> %d', source, target
                    )
                    continue

                conn = get_connection(self._db_url)
                cur = conn.cursor()

                # Delete existing
                cur.execute(
                    f'DELETE FROM {SCHEMA}.orthology'
                    ' WHERE source_tax_id = %s AND target_tax_id = %s',
                    (source, target),
                )
                conn.commit()

                row_count = 0
                with cur.copy(
                    f'COPY {SCHEMA}.orthology'
                    ' (source_id, target_id, source_tax_id, target_tax_id,'
                    ' id_type, resource, n_sources, support)'
                    ' FROM STDIN'
                ) as copy:
                    for src_id, target_ids in table.data.items():
                        for tgt_id in target_ids:
                            meta = table.metadata.get(src_id, {}).get(
                                tgt_id, {}
                            )
                            copy.write_row(
                                (
                                    src_id[:64],
                                    tgt_id[:64],
                                    source,
                                    target,
                                    'genesymbol',
                                    'hcop',
                                    meta.get('n_sources'),
                                    meta.get('support', '')[:500],
                                )
                            )
                            row_count += 1

                conn.commit()
                conn.close()
                _log.info(
                    'Orthology %d -> %d: %d pairs',
                    source,
                    target,
                    row_count,
                )

            except Exception as e:
                _log.error(
                    'Orthology %d -> %d failed: %s', source, target, e
                )

    def build_preset(self, preset: str, parquet_dir: str | None = None):
        """Build using a named preset."""
        from omnipath_utils.db._presets import PRESETS, PARQUET_TABLES

        if preset not in PRESETS:
            raise ValueError(
                f'Unknown preset: {preset}. '
                f'Available: {list(PRESETS.keys())}'
            )

        config = PRESETS[preset]
        _log.info(
            'Building preset "%s": %s', preset, config['description']
        )

        self.build_reference_tables()

        organisms = config['organisms'] or [9606]

        if config.get('ftp'):
            self.populate_from_ftp()  # loads sec_ac + resolver views internally
        else:
            self._run_mappings_parallel(config['mappings'], organisms)
            # Secondary->primary UniProt normalization, always part of handling
            # proteins (ADR 0006); organism-agnostic (tax 0). Skipped for the
            # chemical-only preset (no protein mappings).
            if config['mappings']:
                try:
                    self.load_uniprot_sec_ac()
                except Exception as e:
                    _log.error('Failed to load sec_ac: %s', e)

        if config.get('metabolite'):
            self.populate_metabolites()

        if config.get('mirna'):
            self.populate_mirna(organisms)

        if config.get('orthology'):
            self.populate_orthology(organisms)

        if config.get('reflists'):
            for org in organisms:
                try:
                    self.populate_reflists(org)
                except Exception as e:
                    _log.error('Failed reflists org %d: %s', org, e)

        # Export Parquet files
        if parquet_dir:
            tables = PARQUET_TABLES.get(
                preset, PARQUET_TABLES.get('minimal', [])
            )
            self.export_parquet(tables, parquet_dir)

        _log.info('Preset "%s" build complete', preset)

    # ------------------------------------------------------------------
    # MetaNetX cross-references
    # ------------------------------------------------------------------


    _BIGG_PAIRS = [
        ("bigg", "chebi"),
        ("bigg", "hmdb"),
        ("bigg", "kegg"),
        ("bigg", "metanetx"),
    ]

    _METANETX_PAIRS = [
        ("bigg", "chebi"),
        ("bigg", "hmdb"),
        ("bigg", "kegg"),
        ("kegg", "chebi"),
        ("hmdb", "chebi"),
        ("lipidmaps", "chebi"),
        ("swisslipids", "chebi"),
        ("metanetx", "chebi"),
        ("metanetx", "hmdb"),
        ("metanetx", "kegg"),
        ("metanetx", "bigg"),
    ]

    def _populate_metanetx(self):
        """Build MetaNetX pairwise metabolite mappings."""

        _log.info("Building MetaNetX cross-reference mappings...")

        for src, tgt in self._METANETX_PAIRS:
            try:
                self.populate_mapping(src, tgt, 0, "metanetx")
            except Exception as e:
                _log.warning(
                    "MetaNetX %s -> %s failed: %s", src, tgt, e,
                )

    def _populate_bigg(self):
        """Build BiGG metabolite ID mappings from BiGG Models TSV."""

        _log.info("Building BiGG metabolite cross-reference mappings...")

        for src, tgt in self._BIGG_PAIRS:
            try:
                self.populate_mapping(src, tgt, 0, "bigg")
            except Exception as e:
                _log.warning(
                    "BiGG %s -> %s failed: %s", src, tgt, e,
                )

    def _populate_reactome(self):
        """ChEBI <-> Reactome stable identifiers, from Reactome's own
        ChEBI2Reactome cross-reference file (US5, T089).

        Before this, no ``reactome`` id_type existed at all, which is why
        521 Reactome-only metabolites (a mention citing a bare Reactome
        stable id as its only identifier) could never resolve to anything
        else (research.md R7). Each row is one (ChEBI compound, pathway it
        participates in) pair -- the "pathway_id" field is genuinely the
        Reactome stable id of that participation, despite the generic
        field name shared across ``reactome_old``'s four datasets. Human
        rows only (pathways are organism-scoped; other organisms would
        just add orthologous-pathway noise for the same compound).
        """
        try:
            from pypath.inputs.reactome_old import reactome_chebis
        except ImportError:
            _log.warning('pypath not available for Reactome')
            return

        with Session(self.engine) as session:
            chebi_type = session.query(IdType).filter_by(name='chebi').first()
            reactome_type = (
                session.query(IdType).filter_by(name='reactome').first()
            )
            backend = (
                session.query(Backend).filter_by(name='reactome').first()
            )
            if not chebi_type or not reactome_type or not backend:
                _log.error(
                    'Reactome: missing chebi/reactome id_type or backend'
                )
                return
            chebi_type_id = chebi_type.id
            reactome_type_id = reactome_type.id
            backend_id = backend.id

        chebi_to_reactome: dict[str, set[str]] = {}
        reactome_to_chebi: dict[str, set[str]] = {}
        for row in reactome_chebis():
            if row.organism != 'Homo sapiens' or not row.chebi_id:
                continue
            # id_mapping's own chebi convention is "CHEBI:NNNN" (verified
            # against direct chebi->inchikey and bigg->chebi rows) -- the
            # raw reactome_old row's chebi_id is a bare number.
            chebi_id = f'CHEBI:{row.chebi_id}'
            chebi_to_reactome.setdefault(chebi_id, set()).add(row.pathway_id)
            reactome_to_chebi.setdefault(row.pathway_id, set()).add(chebi_id)

        from omnipath_utils.db._connection import get_connection

        limit = self._effective_limit('reactome')
        total = 0
        for src_type_id, tgt_type_id, data in (
            (chebi_type_id, reactome_type_id, chebi_to_reactome),
            (reactome_type_id, chebi_type_id, reactome_to_chebi),
        ):
            if not data:
                continue
            with Session(self.engine) as session:
                session.execute(
                    text(
                        f'DELETE FROM {SCHEMA}.id_mapping'
                        ' WHERE source_type_id = :src'
                        ' AND target_type_id = :tgt'
                        ' AND ncbi_tax_id = 0 AND backend_id = :bk'
                    ),
                    {'src': src_type_id, 'tgt': tgt_type_id, 'bk': backend_id},
                )
                session.commit()

            conn = get_connection(self._db_url)
            row_count = 0
            try:
                with conn.cursor() as cur:
                    with cur.copy(
                        f'COPY {SCHEMA}.id_mapping'
                        ' (source_type_id, target_type_id, ncbi_tax_id,'
                        ' source_id, target_id, backend_id) FROM STDIN'
                    ) as copy:
                        for src_id, tgt_ids in data.items():
                            for tgt_id in tgt_ids:
                                copy.write_row((
                                    src_type_id, tgt_type_id, 0,
                                    str(src_id)[:64], str(tgt_id)[:64],
                                    backend_id,
                                ))
                                row_count += 1
                                if limit is not None and row_count >= limit:
                                    break
                            if limit is not None and row_count >= limit:
                                break
                conn.commit()
            finally:
                conn.close()
            total += row_count

        _log.info('Reactome: %d chebi<->reactome rows', total)

    # ------------------------------------------------------------------
    # Structure-bearing namespaces (inputs_v2 adapter backends + PubChem)
    # ------------------------------------------------------------------

    # Each (source_namespace, backend) maps the namespace's own ID to the
    # Standard InChIKey -- the chemical resolver's canonical target (see
    # resolver_policy.yaml). InChIKey is the only structure id that fits the
    # 64-char id_mapping column (InChI/SMILES are kept in the metabo layer,
    # not here). Organism-agnostic (ncbi_tax_id=0).
    _STRUCTURE_BACKENDS = [
        "chebi",
        "chembl",
        "lipidmaps",
        "swisslipids",
        "pubchem",
        # T081: hmdb has its own precomputed inchikey column
        # (id_types.yaml's inchikey: backends), same shape as the five
        # above -- it was simply missing from this list.
        "hmdb",
    ]

    def _populate_structures(self):
        """Build namespace -> Standard InChIKey mappings via the structure
        backends (ChEBI, ChEMBL, LIPID MAPS, SwissLipids, PubChem).

        These feed the chemical resolver's canonical InChIKey projection.
        PubChem is the largest by far and honours ``--pubchem-max-records``.
        """

        _log.info("Building structure-bearing (-> InChIKey) mappings...")

        for backend in self._STRUCTURE_BACKENDS:
            try:
                if backend == "pubchem":
                    # PubChem cid->inchikey is ~119M rows -- the single largest
                    # chemical mapping. Use the parallel, streamed loader rather
                    # than the serial materialise-dict-then-one-COPY path.
                    self._populate_pubchem_inchikey()
                else:
                    self.populate_mapping(backend, "inchikey", 0, backend)
            except Exception as e:
                _log.warning(
                    "%s -> inchikey failed: %s", backend, e,
                )

    def _copy_pairs_parallel(
        self,
        pairs,
        src_type_id: int,
        tgt_type_id: int,
        ncbi_tax_id: int,
        backend_id: int,
        n_workers: int | None = None,
        batch_size: int = 20000,
    ) -> int:
        """Fan a stream of ``(source_id, target_id)`` pairs into N concurrent COPYs.

        One feeder (this call) pulls from ``pairs`` -- typically a generator that
        streams + decompresses straight off disk -- batches the rows, and
        dispatches them over a bounded queue to ``n_workers`` worker threads,
        each holding its own connection and a single open ``COPY ... FROM
        STDIN``. This parallelises the write (the bottleneck for a heavily
        indexed target like ``id_mapping``: 3 indexes => ~45k rows/s on one
        COPY) and keeps memory bounded (no full ``{source: {targets}}`` dict).

        Each worker commits its own portion, so an interrupted load can leave a
        partial slice; callers DELETE the slice first and a re-run is idempotent.
        Returns the total row count written.
        """
        import queue
        import threading
        from omnipath_utils.db._connection import get_connection

        if n_workers is None:
            n_workers = min(
                int(os.environ.get("OMNIPATH_BUILD_MAPPING_WORKERS", "8")),
                16,
            )
        n_workers = max(1, n_workers)

        q: queue.Queue = queue.Queue(maxsize=n_workers * 4)
        _SENTINEL = object()
        counts = [0] * n_workers
        errors: list[Exception] = []
        copy_sql = (
            f"COPY {SCHEMA}.id_mapping "
            "(source_type_id, target_type_id, ncbi_tax_id, "
            "source_id, target_id, backend_id) FROM STDIN"
        )

        def worker(idx: int):
            try:
                conn = get_connection(self._db_url)
                with conn.cursor() as cur:
                    cur.execute("SET synchronous_commit = off")
                    with cur.copy(copy_sql) as cp:
                        while True:
                            batch = q.get()
                            if batch is _SENTINEL:
                                q.task_done()
                                break
                            for src_id, tgt_id in batch:
                                cp.write_row((
                                    src_type_id, tgt_type_id, ncbi_tax_id,
                                    src_id, tgt_id, backend_id,
                                ))
                            counts[idx] += len(batch)
                            q.task_done()
                conn.commit()
                conn.close()
            except Exception as e:  # noqa: BLE001
                errors.append(e)
                _log.error("COPY worker %d failed: %s", idx, e)

        threads = [
            threading.Thread(target=worker, args=(i,), daemon=True)
            for i in range(n_workers)
        ]
        for t in threads:
            t.start()

        batch: list = []
        for pair in pairs:
            batch.append(pair)
            if len(batch) >= batch_size:
                q.put(batch)
                batch = []
                if errors:
                    break
        if batch and not errors:
            q.put(batch)
        for _ in range(n_workers):
            q.put(_SENTINEL)
        for t in threads:
            t.join()

        if errors:
            raise errors[0]

        return sum(counts)

    def _populate_pubchem_inchikey(self, n_workers: int | None = None):
        """Parallel, streamed load of the full PubChem CID -> InChIKey mapping.

        The single largest chemical namespace (~119M rows). Rather than
        materialising the whole ``{cid: {inchikey}}`` dict in RAM (~11 GB) and
        draining it through one serial COPY into the 3-index ``id_mapping``
        table, this streams ``pubchem_mapping`` straight into N concurrent COPY
        workers -- bounded memory, ~N x the write throughput. Honours
        ``--pubchem-max-records`` (the per-backend effective limit). (T020/T021)
        """
        from itertools import islice
        from pypath.inputs.pubchem import pubchem_mapping

        with Session(self.engine) as session:
            src = session.query(IdType).filter_by(name="pubchem").first()
            tgt = session.query(IdType).filter_by(name="inchikey").first()
            bk = session.query(Backend).filter_by(name="pubchem").first()
            if not src or not tgt or not bk:
                _log.error("pubchem/inchikey id_type or backend missing")
                return
            src_id, tgt_id, bk_id = src.id, tgt.id, bk.id

        # Replace the existing pubchem -> inchikey slice: delete (committed),
        # then the parallel COPY appends. A re-run is idempotent.
        with Session(self.engine) as session:
            session.execute(
                text(
                    f"DELETE FROM {SCHEMA}.id_mapping WHERE source_type_id = :s "
                    "AND target_type_id = :t AND ncbi_tax_id = 0 "
                    "AND backend_id = :b"
                ),
                {"s": src_id, "t": tgt_id, "b": bk_id},
            )
            session.commit()

        limit = self._effective_limit("pubchem")
        rows = pubchem_mapping("inchikey", source="cid")  # yields (cid, inchikey)
        if limit is not None:
            rows = islice(rows, limit)

        start = time.time()
        n = self._copy_pairs_parallel(rows, src_id, tgt_id, 0, bk_id, n_workers)
        _log.info(
            "PubChem cid->inchikey: %d rows in %.1fs (parallel streamed COPY)",
            n, time.time() - start,
        )
