"""Build miRNA identifier resolver mappings from miRBase.

Interaction resources refer to microRNAs by name (e.g. ``hsa-mir-21`` for the
precursor, ``hsa-miR-21-5p`` for a mature product). This source resolves those
names to the correct miRBase accession **without collapsing the maturation
stages**: precursor names resolve to ``MI#`` (MIRBASE_PRECURSOR) and mature
names to ``MIMAT#`` (MIRBASE_MATURE), each canonicalised onto its own
namespace. Accessions map to themselves (identity), so already-resolved keys
pass through.

Data come from the legacy ``pypath.inputs.mirbase`` tables (on the dlmachine
download stack); all organisms are emitted (``organism=None``), which also
avoids the taxonomy lookup. miRNA names are organism-prefixed and unique, so no
taxonomy column or ambiguity split is needed.
"""

from __future__ import annotations

from pathlib import Path
from collections.abc import Iterable

import polars as pl

from pypath.inputs.mirbase import mirbase_mirna, mirbase_mirna_mature
from pypath.internals.cv_terms import (
    IdentifierNamespaceCv,
    cv_term_label_accession,
)
from omnipath_build.resolver.paths import (
    ensure_mirna_data_dir,
    activate_raw_download_data_dir,
)
from omnipath_build.resolver.identifier_types import (
    IDENTIFIER_TYPE_SCHEMA,
    identifier_type_id,
    identifier_type_rows,
)

MIRNA_IDENTIFIER_LOOKUP_SCHEMA: dict[str, pl.DataType] = {
    'key_identifier_type_id': pl.UInt32,
    'key_value': pl.Utf8,
    'canonical_identifier_type_id': pl.UInt32,
    'canonical_identifier': pl.Utf8,
}

NAME_TYPE = cv_term_label_accession(IdentifierNamespaceCv.NAME)
MIRBASE_PRECURSOR_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.MIRBASE_PRECURSOR
)
MIRBASE_MATURE_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.MIRBASE_MATURE
)

MIRNA_IDENTIFIER_LOOKUP_OUTPUT_FILENAME = 'mirna_identifier_lookup.parquet'
IDENTIFIER_TYPE_OUTPUT_FILENAME = 'identifier_type.parquet'


def _mirna_identifier_rows() -> Iterable[dict]:
    """Yield name/accession -> canonical miRBase accession resolver rows."""

    # Precursors (MI#): primary name, synonym, and accession identity.
    for row in mirbase_mirna(None):
        mi_accession = row[1] if len(row) > 1 else None
        if not mi_accession:
            continue
        name = row[2] if len(row) > 2 else None
        synonym = row[3] if len(row) > 3 else None
        if name:
            yield _row(NAME_TYPE, name, MIRBASE_PRECURSOR_TYPE, mi_accession)
        # Alternative precursor name emitted under NAME too: miRNA names are
        # organism-prefixed and unique, so a synonym can be resolved by the
        # same NAME lookup without a dedicated synonym key type.
        if synonym and synonym != name:
            yield _row(NAME_TYPE, synonym, MIRBASE_PRECURSOR_TYPE, mi_accession)
        yield _row(
            MIRBASE_PRECURSOR_TYPE,
            mi_accession,
            MIRBASE_PRECURSOR_TYPE,
            mi_accession,
        )

    # Mature products (MIMAT#): own name and accession identity. The parent
    # precursor name is deliberately not emitted, so it resolves to MI# only.
    for row in mirbase_mirna_mature(None):
        mimat_accession = row[3] if len(row) > 3 else None
        if not mimat_accession:
            continue
        name = row[1] if len(row) > 1 else None
        if name:
            yield _row(NAME_TYPE, name, MIRBASE_MATURE_TYPE, mimat_accession)
        yield _row(
            MIRBASE_MATURE_TYPE,
            mimat_accession,
            MIRBASE_MATURE_TYPE,
            mimat_accession,
        )


def _row(
    key_type: str,
    key_value: str,
    canonical_type: str,
    canonical_identifier: str,
) -> dict[str, object]:
    return {
        'key_type': key_type,
        'key_value': key_value,
        'canonical_type': canonical_type,
        'canonical_identifier': canonical_identifier,
    }


def _build_mirna_lookup() -> tuple[pl.DataFrame, pl.DataFrame]:
    normalized_rows: list[dict[str, object]] = []
    type_names = {
        NAME_TYPE,
        MIRBASE_PRECURSOR_TYPE,
        MIRBASE_MATURE_TYPE,
    }
    for row in _mirna_identifier_rows():
        normalized_rows.append(
            {
                'key_identifier_type_id': identifier_type_id(row['key_type']),
                'key_value': row['key_value'],
                'canonical_identifier_type_id': identifier_type_id(
                    row['canonical_type']
                ),
                'canonical_identifier': row['canonical_identifier'],
            }
        )

    identifier_types = pl.DataFrame(
        identifier_type_rows(type_names),
        schema=IDENTIFIER_TYPE_SCHEMA,
    )
    if not normalized_rows:
        empty = pl.DataFrame(schema=MIRNA_IDENTIFIER_LOOKUP_SCHEMA)
        return empty, identifier_types

    lookup = (
        pl.DataFrame(normalized_rows, schema=MIRNA_IDENTIFIER_LOOKUP_SCHEMA)
        .filter(
            pl.col('key_value').is_not_null()
            & (pl.col('key_value') != '')
            & pl.col('canonical_identifier').is_not_null()
            & (pl.col('canonical_identifier') != '')
        )
        .unique()
    )
    return lookup, identifier_types


def build_mirna_identifier_lookup() -> pl.DataFrame:
    """Return the miRNA name/accession -> accession lookup as a dataframe."""

    activate_raw_download_data_dir()
    return _build_mirna_lookup()[0]


def materialize_mirna(
    output_dir: str | Path | None = None,
    skip_existing: bool = True,
) -> dict[str, int]:
    """Write miRNA resolver parquet files and return output row counts."""

    output_dir = (
        Path(output_dir)
        if output_dir is not None
        else ensure_mirna_data_dir()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    lookup_path = output_dir / MIRNA_IDENTIFIER_LOOKUP_OUTPUT_FILENAME
    identifier_type_path = output_dir / IDENTIFIER_TYPE_OUTPUT_FILENAME
    if skip_existing and lookup_path.exists() and identifier_type_path.exists():
        print(
            f'[resolver] skip source=mirbase existing_dir={output_dir}',
            flush=True,
        )
        return {
            'mirna_identifier_lookup_rows': _parquet_row_count(lookup_path),
            'identifier_type_rows': _parquet_row_count(identifier_type_path),
        }

    activate_raw_download_data_dir()
    lookup, identifier_types = _build_mirna_lookup()
    lookup.write_parquet(lookup_path)
    identifier_types.write_parquet(identifier_type_path)

    return {
        'mirna_identifier_lookup_rows': lookup.height,
        'identifier_type_rows': identifier_types.height,
    }


def _parquet_row_count(path: Path) -> int:
    return pl.scan_parquet(path).select(pl.len()).collect().item()
