"""Build chemical identifier resolver mappings from supported sources.

Chemical resolver rows normalize source-specific identifiers such as ChEBI,
ChEMBL, HMDB, LipidMaps, RaMP, RefMet, SwissLipids, and PubChem. Structural
chemicals canonicalize to standard InChIKey; remaining ChEBI and RefMet
identifiers can canonicalize through conservative non-structural cross-reference
clusters so residues, fragments, and classes still resolve as chemical-domain
entities.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from dataclasses import dataclass
from collections.abc import Callable, Iterable

import polars as pl

from pypath.inputs_v2.hmdb import resource as hmdb_resource
from pypath.inputs_v2.chebi import resource as chebi_resource
from pypath.inputs_v2.chembl import resource as chembl_resource
from pypath.inputs_v2.refmet import resource as refmet_resource
from pypath.internals.cv_terms import (
    IdentifierNamespaceCv,
    cv_term_label_accession,
)
from pypath.inputs_v2.lipidmaps import resource as lipidmaps_resource
from pypath.inputs_v2.swisslipids import resource as swisslipids_resource
from omnipath_build.resolver.paths import (
    ensure_chemicals_data_dir,
    activate_raw_download_data_dir,
)
from omnipath_build.resolver.parquet import write_parquet_from_dict_rows
from omnipath_build.resolver.identifier_types import (
    IDENTIFIER_TYPE_SCHEMA,
    identifier_type_id,
    identifier_type_rows,
)

CHEMICAL_IDENTIFIER_LOOKUP_SCHEMA: dict[str, pl.DataType] = {
    'key_identifier_type_id': pl.UInt32,
    'key_value': pl.Utf8,
    'canonical_identifier_type_id': pl.UInt32,
    'canonical_identifier': pl.Utf8,
}

CHEBI_SOURCE = 'chebi'
HMDB_SOURCE = 'hmdb'
LIPIDMAPS_SOURCE = 'lipidmaps'
SWISSLIPIDS_SOURCE = 'swisslipids'
CHEMBL_SOURCE = 'chembl'
REFMET_SOURCE = 'refmet'
RAMP_SOURCE = 'ramp'
PUBCHEM_SOURCE = 'pubchem'
CHEMICAL_SOURCES: tuple[str, ...] = (
    CHEBI_SOURCE,
    HMDB_SOURCE,
    LIPIDMAPS_SOURCE,
    SWISSLIPIDS_SOURCE,
    CHEMBL_SOURCE,
    REFMET_SOURCE,
    RAMP_SOURCE,
    PUBCHEM_SOURCE,
)
NONSTRUCTURAL_CROSS_REFERENCE_SOURCE_NAMES = frozenset(
    {CHEBI_SOURCE, REFMET_SOURCE}
)
IDENTIFIER_TYPE_OUTPUT_FILENAME = 'identifier_type.parquet'
CHEMICAL_IDENTIFIER_LOOKUP_PARTITION_DIRNAME = 'lookup'
CHEMICAL_IDENTIFIER_LOOKUP_AMBIGUOUS_PARTITION_DIRNAME = 'ambiguous'
CHEBI_TYPE = cv_term_label_accession(IdentifierNamespaceCv.CHEBI)
CHEMBL_COMPOUND_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.CHEMBL_COMPOUND
)
HMDB_TYPE = cv_term_label_accession(IdentifierNamespaceCv.HMDB)
LIPIDMAPS_TYPE = cv_term_label_accession(IdentifierNamespaceCv.LIPIDMAPS)
RAMP_ID_TYPE = cv_term_label_accession(IdentifierNamespaceCv.RAMP_ID)
REFMET_TYPE = cv_term_label_accession(IdentifierNamespaceCv.REFMET)
SWISSLIPIDS_TYPE = cv_term_label_accession(IdentifierNamespaceCv.SWISSLIPIDS)
STANDARD_INCHI_KEY_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.STANDARD_INCHI_KEY
)
PUBCHEM_COMPOUND_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.PUBCHEM_COMPOUND
)
KEGG_COMPOUND_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.KEGG_COMPOUND
)
CAS_TYPE = cv_term_label_accession(IdentifierNamespaceCv.CAS)
CHEMICAL_SOURCE_IDENTIFIER_TYPES = {
    CHEBI_SOURCE: CHEBI_TYPE,
    CHEMBL_SOURCE: CHEMBL_COMPOUND_TYPE,
    HMDB_SOURCE: HMDB_TYPE,
    LIPIDMAPS_SOURCE: LIPIDMAPS_TYPE,
    PUBCHEM_SOURCE: PUBCHEM_COMPOUND_TYPE,
    RAMP_SOURCE: RAMP_ID_TYPE,
    REFMET_SOURCE: REFMET_TYPE,
    SWISSLIPIDS_SOURCE: SWISSLIPIDS_TYPE,
}


@dataclass(frozen=True)
class _NonstructuralCrossReferencePolicy:
    raw_rows: Callable[[], Iterable[dict]]
    rows: Callable[[Iterable[dict], set[str]], Iterable[dict]]
    structural_hit_sources: frozenset[str] = frozenset()


def _clean(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _clean_inchikey(value: object) -> str | None:
    text = _clean(value)
    if text is None or text.lower() in {'none', 'inchikey=none'}:
        return None
    return text.removeprefix('InChIKey=')


def _clean_inchi(value: object) -> str | None:
    text = _clean(value)
    if text is None or text.lower() in {'none', 'inchi=none'}:
        return None
    return text


def _clean_sequence(value: object, *, upper: bool = False) -> tuple[str, ...]:
    if value is None:
        return ()
    values = value if isinstance(value, list | tuple | set) else (value,)
    result = []
    for item in values:
        text = _clean(item)
        if text is None:
            continue
        result.append(text.upper() if upper else text)
    return tuple(dict.fromkeys(result))


def _chebi_identifier_values(row: dict) -> tuple[str, ...]:
    values = []
    for raw_id in (
        row.get('chebi_id') or row.get('id'),
        *(row.get('alt_ids') or []),
    ):
        match = re.fullmatch(
            r'(?:CHEBI:)?(\d+)',
            str(raw_id or '').strip(),
        )
        if match:
            values.append(match.group(1))
    return tuple(dict.fromkeys(values))


def _structural_hits_by_node(
    rows: Iterable[tuple[set[tuple[str, str]], str | None]],
) -> dict[tuple[str, str], set[str]]:
    hits: dict[tuple[str, str], set[str]] = {}
    for nodes, standard_inchi_key in rows:
        if not standard_inchi_key:
            continue
        for node in nodes:
            hits.setdefault(node, set()).add(standard_inchi_key)
    return hits


def _cluster_lookup_rows(
    candidate_rows: Iterable[tuple[int, dict, set[tuple[str, str]]]],
    structural_hits: dict[tuple[str, str], set[str]],
    *,
    emit_key_types: set[str] | None = None,
) -> tuple[list[dict], set[int], set[tuple[str, str]]]:
    row_sets = list(candidate_rows)
    if not row_sets:
        return [], set(), set()

    node_rows: dict[tuple[str, str], set[int]] = {}
    row_nodes: dict[int, set[tuple[str, str]]] = {}
    for index, _row, nodes in row_sets:
        row_nodes[index] = nodes
        for node in nodes:
            node_rows.setdefault(node, set()).add(index)

    emitted_rows: list[dict] = []
    clustered_row_indices: set[int] = set()
    clustered_nodes: set[tuple[str, str]] = set()
    seen_nodes: set[tuple[str, str]] = set()
    for start_node in sorted(node_rows):
        if start_node in seen_nodes:
            continue
        component_nodes: set[tuple[str, str]] = set()
        component_rows: set[int] = set()
        stack = [start_node]
        seen_nodes.add(start_node)
        while stack:
            node = stack.pop()
            component_nodes.add(node)
            for row_index in node_rows.get(node, ()):
                component_rows.add(row_index)
                for next_node in row_nodes.get(row_index, ()):
                    if next_node in seen_nodes:
                        continue
                    seen_nodes.add(next_node)
                    stack.append(next_node)

        reachable_structures = {
            structure
            for node in component_nodes
            for structure in structural_hits.get(node, ())
        }
        if len(reachable_structures) == 1:
            canonical_type = STANDARD_INCHI_KEY_TYPE
            canonical_identifier = next(iter(reachable_structures))
        else:
            canonical_type, canonical_identifier = (
                _cluster_canonical_identifier(component_nodes)
            )
        emitted_rows.extend(
            {
                'key_type': key_type,
                'key_value': key_value,
                'canonical_type': canonical_type,
                'canonical_identifier': canonical_identifier,
            }
            for key_type, key_value in sorted(component_nodes)
            if emit_key_types is None or key_type in emit_key_types
        )
        clustered_row_indices.update(component_rows)
        clustered_nodes.update(component_nodes)

    return emitted_rows, clustered_row_indices, clustered_nodes


def _cluster_canonical_identifier(
    nodes: set[tuple[str, str]],
) -> tuple[str, str]:
    for preferred_type in (CHEBI_TYPE, REFMET_TYPE):
        values = [
            value for node_type, value in nodes if node_type == preferred_type
        ]
        if values:
            return preferred_type, min(values, key=_identifier_sort_key)
    return min(nodes)


def _identifier_sort_key(value: str) -> tuple[int, int | str]:
    return (0, int(value)) if value.isdigit() else (1, value)


def _chebi_row(row: dict) -> dict | None:
    match = re.fullmatch(
        r'(?:CHEBI:)?(\d+)',
        str(row.get('chebi_id') or row.get('id') or '').strip(),
    )
    key_value = match.group(1) if match else None
    standard_inchi = _clean_inchi(row.get('inchi'))
    standard_inchi_key = _clean_inchikey(row.get('inchikey'))
    if not key_value:
        return None
    if standard_inchi_key:
        return {
            'key_type': CHEBI_TYPE,
            'key_value': key_value,
            'canonical_type': STANDARD_INCHI_KEY_TYPE,
            'canonical_identifier': standard_inchi_key,
            'standard_inchi_key': standard_inchi_key,
            'standard_inchi': standard_inchi,
        }
    return {
        'key_type': CHEBI_TYPE,
        'key_value': key_value,
        'canonical_type': CHEBI_TYPE,
        'canonical_identifier': key_value,
    }


def _chebi_rows(row: dict) -> Iterable[dict]:
    base = _chebi_row(row)
    if base is None:
        return

    seen = set()
    for raw_id in (
        row.get('chebi_id') or row.get('id'),
        *(row.get('alt_ids') or []),
    ):
        match = re.fullmatch(
            r'(?:CHEBI:)?(\d+)',
            str(raw_id or '').strip(),
        )
        if not match:
            continue
        key_value = match.group(1)
        if key_value in seen:
            continue
        seen.add(key_value)
        yield {**base, 'key_value': key_value}


def _chebi_cross_reference_nodes(row: dict) -> set[tuple[str, str]]:
    nodes = {
        (CHEBI_TYPE, key_value)
        for key_value in _chebi_identifier_values(row)
    }
    nodes.update(
        (KEGG_COMPOUND_TYPE, value)
        for value in _clean_sequence(row.get('kegg_compound'), upper=True)
    )
    nodes.update(
        (CAS_TYPE, value)
        for value in _clean_sequence(row.get('cas'))
    )
    return nodes


def _chebi_has_structural_namespace_xref(row: dict) -> bool:
    return any(
        _clean_sequence(row.get(field), upper=field != 'pubchem_compound')
        for field in ('pubchem_compound', 'hmdb', 'lipidmaps')
    )


def _chebi_rows_with_nonstructural_clusters(
    raw_rows: Iterable[dict],
    _structural_hit_sources: set[str] | None = None,
) -> Iterable[dict]:
    rows = list(raw_rows)
    nonstructural_rows: list[tuple[int, dict, set[tuple[str, str]]]] = []
    structural_hits = _structural_hits_by_node(
        (
            _chebi_cross_reference_nodes(row),
            _clean_inchikey(row.get('inchikey')),
        )
        for row in rows
    )

    for index, row in enumerate(rows):
        if _clean_inchikey(row.get('inchikey')):
            yield from _chebi_rows(row)
            continue

        nodes = _chebi_cross_reference_nodes(row)
        if (
            any(node_type != CHEBI_TYPE for node_type, _ in nodes)
            and not _chebi_has_structural_namespace_xref(row)
        ):
            nonstructural_rows.append((index, row, nodes))

    (
        cluster_rows,
        clustered_row_indices,
        clustered_nodes,
    ) = _cluster_lookup_rows(
        nonstructural_rows,
        structural_hits,
    )
    yield from cluster_rows
    for index, row in enumerate(rows):
        if _clean_inchikey(row.get('inchikey')):
            continue
        if (
            index not in clustered_row_indices
            and not (_chebi_cross_reference_nodes(row) & clustered_nodes)
        ):
            yield from _chebi_rows(row)


def _chembl_row(row: dict) -> dict | None:
    key_value = _clean(row.get('chembl_id'))
    standard_inchi_key = _clean_inchikey(row.get('standard_inchi_key'))
    if not key_value or not standard_inchi_key:
        return None
    return {
        'key_type': CHEMBL_COMPOUND_TYPE,
        'key_value': key_value,
        'standard_inchi_key': standard_inchi_key,
    }


def _hmdb_row(row: dict) -> dict | None:
    key_value = _clean(row.get('accession'))
    standard_inchi = _clean_inchi(row.get('inchi'))
    standard_inchi_key = _clean_inchikey(row.get('inchikey'))
    if not key_value or not standard_inchi or not standard_inchi_key:
        return None
    return {
        'key_type': HMDB_TYPE,
        'key_value': key_value,
        'standard_inchi_key': standard_inchi_key,
        'standard_inchi': standard_inchi,
    }


def _lipidmaps_row(row: dict) -> dict | None:
    key_value = _clean(row.get('LM_ID'))
    standard_inchi = _clean_inchi(row.get('INCHI'))
    standard_inchi_key = _clean_inchikey(row.get('INCHI_KEY'))
    if not key_value or not standard_inchi or not standard_inchi_key:
        return None
    return {
        'key_type': LIPIDMAPS_TYPE,
        'key_value': key_value,
        'standard_inchi_key': standard_inchi_key,
        'standard_inchi': standard_inchi,
    }


def _ramp_row(row: dict) -> dict | None:
    key_value = _clean(row.get('ramp_id'))
    standard_inchi_key = _clean_inchikey(row.get('inchi_key'))
    if not key_value or not standard_inchi_key:
        return None
    return {
        'key_type': RAMP_ID_TYPE,
        'key_value': key_value,
        'standard_inchi_key': standard_inchi_key,
    }


def _refmet_row(row: dict) -> dict | None:
    key_value = _clean(row.get('refmet_id') or row.get(' refmet_id'))
    standard_inchi_key = _clean_inchikey(row.get('inchi_key'))
    if not key_value or not standard_inchi_key:
        return None
    return {
        'key_type': REFMET_TYPE,
        'key_value': key_value,
        'standard_inchi_key': standard_inchi_key,
    }


def _refmet_cross_reference_nodes(row: dict) -> set[tuple[str, str]]:
    nodes: set[tuple[str, str]] = set()
    key_value = _clean(row.get('refmet_id') or row.get(' refmet_id'))
    if key_value:
        nodes.add((REFMET_TYPE, key_value))
    nodes.update(
        (CHEBI_TYPE, key_value)
        for key_value in _clean_sequence(row.get('chebi_id'))
    )
    return nodes


def _refmet_rows_with_nonstructural_clusters(
    raw_rows: Iterable[dict],
    structural_hit_sources: set[str] | None = None,
) -> Iterable[dict]:
    rows = list(raw_rows)
    nonstructural_rows: list[tuple[int, dict, set[tuple[str, str]]]] = []
    structural_hits = _structural_hits_by_node(
        (
            _refmet_cross_reference_nodes(row),
            _clean_inchikey(row.get('inchi_key')),
        )
        for row in rows
    )
    _merge_structural_hits(
        structural_hits,
        _structural_hits_for_sources(structural_hit_sources or set()),
    )

    for index, row in enumerate(rows):
        if _clean_inchikey(row.get('inchi_key')):
            mapped = _refmet_row(row)
            if mapped is not None:
                yield mapped
            continue

        nodes = _refmet_cross_reference_nodes(row)
        refmet_id = _clean(row.get('refmet_id') or row.get(' refmet_id'))
        if (REFMET_TYPE, refmet_id) in nodes:
            nonstructural_rows.append((index, row, nodes))

    cluster_rows, _clustered_row_indices, _clustered_nodes = (
        _cluster_lookup_rows(
            nonstructural_rows,
            structural_hits,
            emit_key_types={REFMET_TYPE},
        )
    )
    yield from cluster_rows


def _chebi_structural_hits_by_node() -> dict[tuple[str, str], set[str]]:
    return _structural_hits_by_node(
        (
            _chebi_cross_reference_nodes(row),
            _clean_inchikey(row.get('inchikey')),
        )
        for row in chebi_resource.molecules.raw()
    )


STRUCTURAL_HIT_SOURCE_BUILDERS: dict[
    str,
    Callable[[], dict[tuple[str, str], set[str]]],
] = {
    CHEBI_SOURCE: _chebi_structural_hits_by_node,
}


def _structural_hits_for_sources(
    sources: Iterable[str],
) -> dict[tuple[str, str], set[str]]:
    structural_hits: dict[tuple[str, str], set[str]] = {}
    for source in sources:
        builder = STRUCTURAL_HIT_SOURCE_BUILDERS[source]
        _merge_structural_hits(structural_hits, builder())
    return structural_hits


def _merge_structural_hits(
    target: dict[tuple[str, str], set[str]],
    source: dict[tuple[str, str], set[str]],
) -> None:
    for node, structures in source.items():
        target.setdefault(node, set()).update(structures)


def _swisslipids_row(row: dict) -> dict | None:
    key_value = _clean(row.get('Lipid ID'))
    standard_inchi = _clean_inchi(row.get('InChI (pH7.3)'))
    standard_inchi_key = _clean_inchikey(row.get('InChI key (pH7.3)'))
    if not key_value or not standard_inchi or not standard_inchi_key:
        return None
    return {
        'key_type': SWISSLIPIDS_TYPE,
        'key_value': key_value,
        'standard_inchi_key': standard_inchi_key,
        'standard_inchi': standard_inchi,
    }


NONSTRUCTURAL_CROSS_REFERENCE_POLICIES: dict[
    str,
    _NonstructuralCrossReferencePolicy,
] = {
    CHEBI_SOURCE: _NonstructuralCrossReferencePolicy(
        raw_rows=chebi_resource.molecules.raw,
        rows=_chebi_rows_with_nonstructural_clusters,
    ),
    REFMET_SOURCE: _NonstructuralCrossReferencePolicy(
        raw_rows=refmet_resource.metabolites.raw,
        rows=_refmet_rows_with_nonstructural_clusters,
        structural_hit_sources=frozenset({CHEBI_SOURCE}),
    ),
}


_CHEMICAL_DATASETS: dict[str, tuple[object, Callable[[dict], dict | None]]] = {
    CHEBI_SOURCE: (chebi_resource.molecules, _chebi_row),
    CHEMBL_SOURCE: (chembl_resource.molecules, _chembl_row),
    HMDB_SOURCE: (hmdb_resource.metabolites, _hmdb_row),
    LIPIDMAPS_SOURCE: (lipidmaps_resource.lipids, _lipidmaps_row),
    REFMET_SOURCE: (refmet_resource.metabolites, _refmet_row),
    SWISSLIPIDS_SOURCE: (swisslipids_resource.lipids, _swisslipids_row),
}


def _validate_chemical_sources(sources: Iterable[str]) -> tuple[str, ...]:
    selected = tuple(sources)
    unsupported = sorted(set(selected) - set(CHEMICAL_SOURCES))
    if unsupported:
        raise ValueError(f'Unsupported chemical source(s): {unsupported}')
    return selected


def _order_chemical_sources(sources: Iterable[str]) -> tuple[str, ...]:
    """Return sources in dependency order, independent lookups before filters."""

    selected = set(_validate_chemical_sources(sources))
    return tuple(source for source in CHEMICAL_SOURCES if source in selected)


def _chemical_identifier_rows(
    sources: Iterable[str],
    max_records: int | None = None,
    pubchem_url: str | Path | None = None,
    pubchem_shards: int | None = None,
    chemical_lookup_sources: Iterable[str] | None = None,
) -> Iterable[dict]:
    selected_sources = _validate_chemical_sources(sources)
    completed_sources = set(chemical_lookup_sources or ())
    available_cluster_sources = (
        set(selected_sources) | completed_sources
    ) & NONSTRUCTURAL_CROSS_REFERENCE_SOURCE_NAMES
    for source in selected_sources:
        cluster_policy = NONSTRUCTURAL_CROSS_REFERENCE_POLICIES.get(source)
        if cluster_policy is not None:
            emitted = 0
            structural_hit_sources = (
                available_cluster_sources & cluster_policy.structural_hit_sources
            )
            for row in cluster_policy.rows(
                cluster_policy.raw_rows(),
                structural_hit_sources,
            ):
                yield row
                emitted += 1
                if max_records is not None and emitted >= max_records:
                    break
            continue
        if source == PUBCHEM_SOURCE:
            from omnipath_build.resolver.sources.pubchem import (
                iter_pubchem_compound_rows,
            )

            rows = iter_pubchem_compound_rows(
                pubchem_url,
                shard_count=pubchem_shards,
            )
            emitted = 0
            for row in rows:
                yield row
                emitted += 1
                if max_records is not None and emitted >= max_records:
                    break
            continue
        if source == RAMP_SOURCE:
            from pypath.inputs_v2.rampdb import resource as ramp_resource

            emitted = 0
            for raw_row in ramp_resource.chem_props.raw():
                row = _ramp_row(raw_row)
                if row is None:
                    continue
                yield row
                emitted += 1
                if max_records is not None and emitted >= max_records:
                    break
            continue

        dataset, mapper = _CHEMICAL_DATASETS[source]
        emitted = 0
        for raw_row in dataset.raw():
            rows = (row for row in (mapper(raw_row),) if row is not None)
            for row in rows:
                yield row
                emitted += 1
                if max_records is not None and emitted >= max_records:
                    break
            if max_records is not None and emitted >= max_records:
                break


def build_chemical_identifier_lookup(
    sources: Iterable[str] = CHEMICAL_SOURCES,
    max_records: int | None = None,
    pubchem_url: str | Path | None = None,
    pubchem_shards: int | None = None,
) -> pl.DataFrame:
    """Return normalized chemical identifier lookup rows as a dataframe."""

    activate_raw_download_data_dir()
    rows = list(
        _chemical_identifier_rows(
            sources,
            max_records=max_records,
            pubchem_url=pubchem_url,
            pubchem_shards=pubchem_shards,
        )
    )
    if not rows:
        return pl.DataFrame(schema=CHEMICAL_IDENTIFIER_LOOKUP_SCHEMA)
    return _split_chemical_identifier_lookup(rows)[0]


def materialize_chemical_sources(
    sources: Iterable[str],
    output_dir: str | Path | None = None,
    max_records: int | None = None,
    pubchem_url: str | Path | None = None,
    pubchem_shards: int | None = None,
    jobs: int = 1,
    skip_existing: bool = True,
    continue_on_error: bool = False,
) -> dict[str, int]:
    """Write chemical resolver parquet files and return output row counts."""

    selected = _order_chemical_sources(sources)
    output_dir = (
        Path(output_dir)
        if output_dir is not None
        else ensure_chemicals_data_dir()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    activate_raw_download_data_dir()
    identifier_type_path = output_dir / IDENTIFIER_TYPE_OUTPUT_FILENAME
    existing_identifier_types = _read_existing_identifier_types(
        identifier_type_path
    )
    existing_sources = _loaded_partitioned_chemical_sources(output_dir)
    sources_to_skip = existing_sources if skip_existing else set()

    rows: list[dict] = []
    completed_sources: list[str] = [
        source for source in CHEMICAL_SOURCES if source in existing_sources
    ]
    written_lookup_rows = 0
    written_ambiguous_rows = 0

    for source in selected:
        if source in sources_to_skip:
            print(
                f'[resolver] skip source={source} existing_dir={output_dir}',
                flush=True,
            )
            if existing_identifier_types is None:
                existing_identifier_types = _merge_identifier_types(
                    existing_identifier_types,
                    _identifier_types_for_source(source),
                )
            continue
        lookup_sources = tuple(completed_sources)

        if source == 'pubchem':
            try:
                streamed_counts = _write_streaming_pubchem_lookup_files(
                    output_dir,
                    source=pubchem_url,
                    max_records=max_records,
                    pubchem_shards=pubchem_shards,
                    jobs=jobs,
                )
            except Exception as exc:
                if not continue_on_error:
                    raise
                print(
                    '[warning] '
                    f'[resolver.{source}] materialize failed; continuing: '
                    f'{exc.__class__.__name__}: {exc}',
                    file=sys.stderr,
                    flush=True,
                )
                continue
            else:
                written_lookup_rows += streamed_counts[
                    'chemical_identifier_lookup_rows'
                ]
                written_ambiguous_rows += streamed_counts[
                    'chemical_identifier_lookup_ambiguous_rows'
                ]
                existing_identifier_types = _merge_identifier_types(
                    existing_identifier_types,
                    _identifier_types_for_source(source),
                )
                completed_sources.append(source)
            continue

        try:
            rows = list(
                _chemical_identifier_rows(
                    (source,),
                    max_records=max_records,
                    pubchem_url=pubchem_url,
                    pubchem_shards=pubchem_shards,
                    chemical_lookup_sources=lookup_sources,
                )
            )
        except Exception as exc:
            if not continue_on_error:
                raise
            print(
                '[warning] '
                f'[resolver.{source}] materialize failed; continuing: '
                f'{exc.__class__.__name__}: {exc}',
                file=sys.stderr,
                flush=True,
            )
            continue
        else:
            lookup, ambiguous, identifier_types = (
                _write_chemical_source_partition_files(
                    source,
                    rows,
                    output_dir,
                )
            )
            written_lookup_rows += lookup.height
            written_ambiguous_rows += ambiguous.height
            existing_identifier_types = _merge_identifier_types(
                existing_identifier_types,
                identifier_types,
            )
            completed_sources.append(source)
            rows = []

    if existing_identifier_types is None:
        existing_identifier_types = pl.DataFrame(
            identifier_type_rows({STANDARD_INCHI_KEY_TYPE}),
            schema=IDENTIFIER_TYPE_SCHEMA,
        )
    _write_parquet_atomic(existing_identifier_types, identifier_type_path)

    return {
        'chemical_identifier_lookup_rows': written_lookup_rows,
        'chemical_identifier_lookup_ambiguous_rows': written_ambiguous_rows,
        'identifier_type_rows': existing_identifier_types.height,
    }


def _write_chemical_source_partition_files(
    source: str,
    rows: Iterable[dict],
    output_dir: Path,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Write one source partition, replacing any earlier partition atomically."""

    lookup, ambiguous, identifier_types = _split_chemical_identifier_lookup(
        rows
    )

    _write_parquet_atomic(
        lookup,
        _chemical_lookup_partition_path(output_dir, source),
    )
    _write_parquet_atomic(
        ambiguous,
        _chemical_ambiguous_partition_path(output_dir, source),
    )

    return lookup, ambiguous, identifier_types


def _write_streaming_pubchem_lookup_files(
    output_dir: Path,
    *,
    source: str | Path | None,
    max_records: int | None,
    pubchem_shards: int | None,
    jobs: int,
) -> dict[str, int]:
    from omnipath_build.resolver.sources.pubchem import (
        iter_pubchem_compound_rows,
        iter_pubchem_lookup_parquet_rows,
        materialize_pubchem_compound_shards,
    )

    lookup_path = _chemical_lookup_partition_path(output_dir, PUBCHEM_SOURCE)
    ambiguous_path = _chemical_ambiguous_partition_path(
        output_dir,
        PUBCHEM_SOURCE,
    )

    if max_records is not None or jobs <= 1:
        pubchem_rows = iter_pubchem_compound_rows(
            source,
            shard_count=pubchem_shards,
        )
    else:
        shard_paths = materialize_pubchem_compound_shards(
            output_dir,
            source=source,
            pubchem_shards=pubchem_shards,
            jobs=jobs,
        )
        pubchem_rows = iter_pubchem_lookup_parquet_rows(
            path for path, _ in shard_paths
        )

    if max_records is not None:
        pubchem_rows = _take(pubchem_rows, max_records)

    tmp_lookup_path = _temporary_parquet_path(lookup_path)
    lookup_row_count = write_parquet_from_dict_rows(
        (
            pubchem_rows
            if max_records is None and jobs > 1
            else _normalized_chemical_lookup_rows(pubchem_rows)
        ),
        CHEMICAL_IDENTIFIER_LOOKUP_SCHEMA,
        tmp_lookup_path,
    )
    tmp_lookup_path.replace(lookup_path)

    ambiguous = _empty_chemical_lookup()
    _write_parquet_atomic(ambiguous, ambiguous_path)

    return {
        'chemical_identifier_lookup_rows': lookup_row_count,
        'chemical_identifier_lookup_ambiguous_rows': ambiguous.height,
        'identifier_type_rows': _identifier_types_for_source(
            PUBCHEM_SOURCE
        ).height,
    }


def _chemical_lookup_partition_path(output_dir: Path, source: str) -> Path:
    return (
        output_dir
        / CHEMICAL_IDENTIFIER_LOOKUP_PARTITION_DIRNAME
        / f'{source}.parquet'
    )


def _chemical_ambiguous_partition_path(output_dir: Path, source: str) -> Path:
    return (
        output_dir
        / CHEMICAL_IDENTIFIER_LOOKUP_AMBIGUOUS_PARTITION_DIRNAME
        / f'{source}.parquet'
    )


def _temporary_parquet_path(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.with_name(f'.{path.name}.tmp')


def _write_parquet_atomic(frame: pl.DataFrame, path: Path) -> None:
    tmp_path = _temporary_parquet_path(path)
    frame.write_parquet(tmp_path)
    tmp_path.replace(path)


def _read_existing_identifier_types(path: Path) -> pl.DataFrame | None:
    return pl.read_parquet(path) if path.exists() else None


def _identifier_types_for_source(source: str) -> pl.DataFrame:
    type_names = {
        STANDARD_INCHI_KEY_TYPE,
        CHEMICAL_SOURCE_IDENTIFIER_TYPES[source],
    }
    return pl.DataFrame(
        identifier_type_rows(type_names),
        schema=IDENTIFIER_TYPE_SCHEMA,
    )


def _merge_identifier_types(
    existing: pl.DataFrame | None,
    update: pl.DataFrame,
) -> pl.DataFrame:
    frames = [frame for frame in (existing, update) if frame is not None]
    return (
        pl.concat(frames, how='vertical_relaxed')
        .unique(subset=['identifier_type_id'])
        .sort('identifier_type_id')
    )


def _loaded_partitioned_chemical_sources(output_dir: Path) -> set[str]:
    lookup_dir = output_dir / CHEMICAL_IDENTIFIER_LOOKUP_PARTITION_DIRNAME
    if not lookup_dir.exists():
        return set()
    return {
        path.stem
        for path in lookup_dir.glob('*.parquet')
        if path.stem in CHEMICAL_SOURCES
    }


def _take(rows: Iterable[dict], max_records: int) -> Iterable[dict]:
    emitted = 0
    for row in rows:
        if emitted >= max_records:
            break
        yield row
        emitted += 1


def _normalized_chemical_lookup_rows(rows: Iterable[dict]) -> Iterable[dict]:
    for row in rows:
        key_type = row.get('key_type')
        if key_type is None:
            continue
        key_type = str(key_type)
        canonical_type = str(
            row.get('canonical_type')
            or row.get('canonical_identifier_type')
            or STANDARD_INCHI_KEY_TYPE
        )
        canonical_identifier = (
            row.get('canonical_identifier')
            or row.get('standard_inchi_key')
        )
        yield {
            'key_identifier_type_id': identifier_type_id(key_type),
            'key_value': row.get('key_value'),
            'canonical_identifier_type_id': identifier_type_id(canonical_type),
            'canonical_identifier': canonical_identifier,
        }


def _split_chemical_identifier_lookup(
    rows: Iterable[dict],
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    normalized_rows: list[dict[str, object]] = []
    type_names = {STANDARD_INCHI_KEY_TYPE}
    for row in rows:
        key_type = row.get('key_type')
        if key_type is None:
            continue
        key_type = str(key_type)
        type_names.add(key_type)
        type_names.add(
            str(
                row.get('canonical_type')
                or row.get('canonical_identifier_type')
                or STANDARD_INCHI_KEY_TYPE
            )
        )
        normalized_rows.extend(_normalized_chemical_lookup_rows((row,)))

    identifier_types = pl.DataFrame(
        identifier_type_rows(type_names),
        schema=IDENTIFIER_TYPE_SCHEMA,
    )
    if not normalized_rows:
        return (
            _empty_chemical_lookup(),
            _empty_chemical_lookup(),
            identifier_types,
        )

    lookup = (
        pl.DataFrame(normalized_rows, schema=CHEMICAL_IDENTIFIER_LOOKUP_SCHEMA)
        .filter(
            pl.col('key_value').is_not_null()
            & (pl.col('key_value') != '')
            & pl.col('canonical_identifier').is_not_null()
            & (pl.col('canonical_identifier') != '')
        )
        .unique()
    )
    unambiguous, ambiguous = _split_chemical_lookup_frame(lookup)
    return unambiguous, ambiguous, identifier_types


def _empty_chemical_lookup() -> pl.DataFrame:
    return pl.DataFrame(schema=CHEMICAL_IDENTIFIER_LOOKUP_SCHEMA)


def _split_chemical_lookup_frame(
    lookup: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    if lookup.is_empty():
        empty = _empty_chemical_lookup()
        return empty, empty

    join_keys = [
        'key_identifier_type_id',
        'key_value',
        'canonical_identifier_type_id',
    ]
    ambiguous_keys = (
        lookup.group_by(join_keys)
        .agg(
            pl.col('canonical_identifier')
            .n_unique()
            .alias('canonical_identifier_count')
        )
        .filter(pl.col('canonical_identifier_count') > 1)
        .select(join_keys)
    )
    if ambiguous_keys.is_empty():
        return lookup, _empty_chemical_lookup()

    ambiguous = lookup.join(ambiguous_keys, on=join_keys, how='semi')
    unambiguous = lookup.join(ambiguous_keys, on=join_keys, how='anti')
    return unambiguous, ambiguous
