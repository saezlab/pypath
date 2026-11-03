"""Main downstream order over published, already-resolved relational tables.

The loader owns connections, schema locks, committed checkpoints and retries.
Imported main steps retain their commits until the loader supplies its commit
facade. No input parser, resource build or resolver entrypoint is called here.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, is_dataclass
import time
from types import SimpleNamespace

from psycopg2 import sql
from psycopg2.extras import Json

from omnipath_core.resource_metadata import resource_metadata
from .chemical_resolution_level import (
    rebuild_chemical_resolution_levels,
    rebuild_chemical_ambiguous_name_candidates,
)
from .classify import (
    classify_chemical_class,
    classify_metabolic_domain,
    classify_interaction_class,
)
from .db import (
    rebuild_derived_tables,
    rebuild_interaction_tables,
    rebuild_bitmap_tables,
    rebuild_resource_overlap_summary,
    sync_resources_table,
    sync_data_source_licenses,
    emit_build_manifest,
)
from .labels import populate_entity_labels, populate_entity_name, populate_chemical_labels
from .authority import populate_identifier_authority
from .resource_declarations import IDENTIFIER_MINTS, SOURCE_COMMIT as RESOURCE_DECLARATION_COMMIT
from .db.resources import _populate_build_capability_table, INTERACTION_RECORD_STEP
from .db.derived_tables import InteractionDeriveStats


def _search_path(conn, schema):
    with conn.cursor() as cur:
        cur.execute(sql.SQL("SET search_path = {}, public").format(sql.Identifier(schema)))


def _reported(value):
    return asdict(value) if is_dataclass(value) else value


def _metadata_inputs(published_metadata):
    records = published_metadata or {}
    if isinstance(records, Mapping) and isinstance(records.get("resources"), Mapping):
        records = records["resources"]
    discovered = {}
    for source, manifest in records.items():
        if not isinstance(manifest, Mapping):
            continue
        metadata = resource_metadata(source, dict(manifest))
        config = SimpleNamespace(
            name=metadata.get("name") or source,
            description=metadata.get("description"),
            url=metadata.get("website") or metadata.get("url"),
            license=metadata.get("license"),
            pubmed=metadata.get("pubmed"),
            annotation_ontologies=metadata.get("annotation_ontologies", []),
            mints=metadata.get("mints") or IDENTIFIER_MINTS.get(source, ()),
        )
        # These are data-only resource descriptors; there is no callable input.
        discovered[source] = [
            SimpleNamespace(
                function_name="resource", qualified_module=None, call=SimpleNamespace(config=config)
            )
        ]
    return discovered


def _availability(conn, schema):
    """Declare source information absent from unchanged serving artifacts."""
    limitations = [
        {
            "capability": "original_entity_occurrences",
            "available": False,
            "provider": "published Parquet",
            "reason": "Serving files aggregate entities; original occurrence identity is absent",
        },
        {
            "capability": "original_resolver_status",
            "available": False,
            "provider": "published Parquet",
            "reason": "Published identities have status published; no original matched flag exists",
        },
        {
            "capability": "gene_anchored_identity",
            "available": False,
            "provider": "published Parquet",
            "reason": "Source-typed identities remain unchanged; supported gene references are stored separately in entity_reference_context",
        },
        {
            "capability": "metabolic_domain_source_priority",
            "available": False,
            "provider": "published Parquet",
            "reason": "Original compound superclass/class and subsystem priority attributes absent",
        },
        {
            "capability": "participant_role_precision",
            "available": False,
            "provider": "published Biolink statements",
            "reason": "Some ligand/receptor sources publish generic verbs without original occurrence-role flags",
        },
        {
            "capability": "mechanism_origin_for_activity_claims",
            "available": False,
            "provider": "published Biolink statements",
            "reason": "Mechanisms dataset is explicit; activity claims may omit their original mechanism identifier",
        },
        {
            "capability": "intact_confidence_original_value_shape",
            "available": False,
            "provider": "published Biolink statements",
            "reason": "IntAct numeric confidence extraction does not preserve the original prefixed scalar spelling",
        },
        {
            "capability": "pharmacological_class_precision",
            "available": False,
            "provider": "published Biolink statements",
            "reason": "Direction qualifiers alone do not distinguish orthosteric/allosteric roles",
        },
        {
            "capability": "transport_event_type",
            "available": True,
            "provider": "published Biolink statements",
            "reason": "A molecular_activity states its SBO reaction type (Biolink type slot)",
        },
        {
            "capability": "ontology_definition_comment_distinction",
            "available": False,
            "provider": "published Parquet",
            "reason": "Ontology definitions and comments share description; source distinction absent",
        },
    ]
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL("SELECT capabilities FROM {}.build_manifest").format(sql.Identifier(schema))
        )
        row = cur.fetchone()
        current = row[0] or [] if row else []
        cur.execute(
            sql.SQL("UPDATE {}.build_manifest SET capabilities=%s").format(sql.Identifier(schema)),
            [Json([*current, *limitations])],
        )
        _populate_build_capability_table(cur, schema, [*current, *limitations])
    conn.commit()
    return limitations


def run_main_derivations(conn, *, schema="public", progress=False, published_metadata=None):
    """Run main's shared derivations in main CLI order; return real phase timings."""
    results = {}
    phase_seconds = {}

    def run(name, call):
        _search_path(conn, schema)
        started = time.monotonic()
        if progress:
            print(f"[main-align] phase={name} event=start", flush=True)
        result = call()
        phase_seconds[name] = time.monotonic() - started
        results[name] = _reported(result)
        if progress:
            print(
                f"[main-align] phase={name} event=done seconds={phase_seconds[name]:.6f}",
                flush=True,
            )
        return result

    # Resource descriptors are supplied from pinned metadata. Discovery/importing
    # parsers is an earlier-build concern and is deliberately not repeated.
    discovered = _metadata_inputs(published_metadata)
    if not discovered:
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL("SELECT name FROM {}.data_source ORDER BY name").format(
                    sql.Identifier(schema)
                )
            )
            discovered = _metadata_inputs({name: {} for (name,) in cur.fetchall()})
    configs = {source: functions[0].call.config for source, functions in discovered.items()}

    def authority():
        with conn.cursor() as cur:
            rows = populate_identifier_authority(cur, schema, configs, progress=progress)
        conn.commit()
        return {
            "rows": rows,
            "declarations_available": any(config.mints for config in configs.values()),
            "fallback_declaration_commit": RESOURCE_DECLARATION_COMMIT,
        }

    run("identifier_authority", authority)
    run(
        "derived_tables",
        lambda: rebuild_derived_tables(
            conn, schema=schema, progress=progress, interactions=False, use_external_mappings=False
        ),
    )
    run(
        "chemical_resolution_levels",
        lambda: rebuild_chemical_resolution_levels(conn, schema=schema, progress=progress),
    )
    run(
        "chemical_ambiguous_name_candidates",
        lambda: rebuild_chemical_ambiguous_name_candidates(conn, schema=schema, progress=progress),
    )
    run("chemical_class", lambda: classify_chemical_class(conn, schema=schema))
    run("interaction_class", lambda: classify_interaction_class(conn, schema=schema))
    interaction = run(
        "interactions", lambda: rebuild_interaction_tables(conn, schema=schema, progress=progress)
    )
    run("entity_labels", lambda: populate_entity_labels(conn, schema=schema))
    run("entity_name", lambda: populate_entity_name(conn, schema=schema))
    run("chemical_labels", lambda: populate_chemical_labels(conn, schema=schema))
    # After the labels: the domain is read from the labels of the classes and
    # pathways a chemical is linked to.
    run("metabolic_domain", lambda: classify_metabolic_domain(conn, schema=schema))
    run("bitmaps", lambda: rebuild_bitmap_tables(conn, schema=schema, progress=progress))
    run(
        "resource_overlap",
        lambda: rebuild_resource_overlap_summary(conn, schema=schema, progress=progress),
    )
    run(
        "resources",
        lambda: sync_resources_table(conn, discovered, schema=schema, prefer_bitmaps=True),
    )
    run("licenses", lambda: sync_data_source_licenses(conn, schema=schema))
    run(
        "build_manifest",
        lambda: emit_build_manifest(
            conn,
            schema=schema,
            derive_cost=_interaction_derive_cost(interaction),
            deferral_cost=_interaction_deferral_cost(interaction),
            utils_db_url=None,
        ),
    )
    results["compatibility_limitations"] = _availability(conn, schema)
    return {"phase_seconds": phase_seconds, "results": results}


def _interaction_derive_cost(
    stats: InteractionDeriveStats | None,
) -> dict[str, dict[str, object]] | None:
    """What the interaction projection cost this run, for the build manifest.

    Read off the ``InteractionDeriveStats`` the derive step already returned
    rather than re-derived. **Seconds and rows are reported per step**, each
    taken from ``step_seconds``, rather than the projection's whole wall clock
    being attributed to one table — these are the numbers the projection's
    share of the build-cost budget is argued against, and a step total cannot
    say which step overran.

    The projection writes ``interaction_fact_resource`` and stops, so the
    ``interaction_fact_combined`` line is gone with the table it priced. The
    ``scope_cost`` its neighbour used to hand the manifest went with it — that
    neighbour reported one scope, the all-resources collapse, and the fold is a
    query now and materialises nothing, so this build has no scope to report
    rather than a scope reported as free.

    Returns ``None`` when the projection did not run, so ``--no-interactions``
    records no cost instead of claiming zeros.
    """
    if stats is None:
        return None
    step_seconds = stats.step_seconds or {}
    return {
        INTERACTION_RECORD_STEP: {
            "seconds": step_seconds.get(INTERACTION_RECORD_STEP),
            "rows": stats.records,
        },
        "interaction_party": {"rows": stats.parties},
        "interaction_header": {
            "seconds": step_seconds.get("interaction_header"),
            "rows": stats.interactions,
        },
        # The two staging steps carry no rows of their own and are reported
        # anyway, so the step seconds add up to the projection's wall clock
        # rather than leaving an unattributed remainder. A ceiling argument
        # made against a partial accounting is an argument about the wrong
        # number.
        "interaction_class_evidence": {
            "seconds": step_seconds.get("interaction_class_evidence"),
        },
        "interaction_evidence_fold": {
            "seconds": step_seconds.get("interaction_evidence_fold"),
        },
    }


def _interaction_deferral_cost(
    stats: InteractionDeriveStats | None,
) -> dict[str, object] | None:
    """What deferring the constraints over the load bought, for the manifest.

    Read off the record the derive step already returned, and handed to
    ``emit_build_manifest`` unchanged: the normalising is that function's job,
    including the rule that an unmeasured field stays ``null`` rather than
    becoming a zero, so that a build which ran without the deferral and a
    deferral that saved nothing do not read alike.

    Returns ``None`` when the projection did not run at all, so
    ``--no-interactions`` records no deferral rather than claiming one that
    saved nothing.
    """
    if stats is None:
        return None
    return dict(stats.deferral) or None
