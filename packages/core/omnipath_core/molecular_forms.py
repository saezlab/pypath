"""Occurrence-scoped molecular descriptions shared by inputs and serving outputs.

Missing fields (including empty feature lists) mean unspecified, never canonical
isoform, unmodified or wild type. Product entity keys are assigned by resolution
only. Positions retain their source coordinate reference; unknown references
must not be silently replaced with a selected reference protein.
"""

from __future__ import annotations

from collections.abc import Mapping
import re
from typing import Any, NamedTuple

import pyarrow as pa


class SequenceIdentifier(NamedTuple):
    ns: str
    id: str


class CoordinateReference(NamedTuple):
    identifier: SequenceIdentifier | None = None
    coordinate_system: str = "unknown"
    position_base: int | None = None


class MolecularModification(NamedTuple):
    term: str | None = None
    residue: str | None = None
    position: int | None = None
    end_position: int | None = None
    coordinate_reference: CoordinateReference | None = None
    description: str | None = None


class MolecularVariant(NamedTuple):
    identifier: SequenceIdentifier | None = None
    reference: str | None = None
    alternate: str | None = None
    position: int | None = None
    end_position: int | None = None
    coordinate_reference: CoordinateReference | None = None
    description: str | None = None


class MolecularRegion(NamedTuple):
    type: str | None = None
    identifier: SequenceIdentifier | None = None
    position: int | None = None
    end_position: int | None = None
    coordinate_reference: CoordinateReference | None = None
    description: str | None = None


class MolecularForm(NamedTuple):
    protein_entity_key: str | None = None
    transcript_entity_key: str | None = None
    isoform_identifier: SequenceIdentifier | None = None
    sequence_identifiers: list[SequenceIdentifier] | None = None
    modifications: list[MolecularModification] | None = None
    variants: list[MolecularVariant] | None = None
    regions: list[MolecularRegion] | None = None


SEQUENCE_IDENTIFIER_STRUCT = pa.struct([("ns", pa.string()), ("id", pa.string())])
COORDINATE_REFERENCE_STRUCT = pa.struct(
    [
        ("identifier", SEQUENCE_IDENTIFIER_STRUCT),
        ("coordinate_system", pa.string()),
        ("position_base", pa.int32()),
    ]
)
MODIFICATION_STRUCT = pa.struct(
    [
        ("term", pa.string()),
        ("residue", pa.string()),
        ("position", pa.int64()),
        ("end_position", pa.int64()),
        ("coordinate_reference", COORDINATE_REFERENCE_STRUCT),
        ("description", pa.string()),
    ]
)
VARIANT_STRUCT = pa.struct(
    [
        ("identifier", SEQUENCE_IDENTIFIER_STRUCT),
        ("reference", pa.string()),
        ("alternate", pa.string()),
        ("position", pa.int64()),
        ("end_position", pa.int64()),
        ("coordinate_reference", COORDINATE_REFERENCE_STRUCT),
        ("description", pa.string()),
    ]
)
# A source-reported sequence segment: the extent of a processed chain, fragment
# or mature product within its precursor, or an annotated segment such as a
# signal peptide. ``type`` keeps the source's region kind; position fields share
# the feature semantics, so an unknown start or end stays unspecified.
REGION_STRUCT = pa.struct(
    [
        ("type", pa.string()),
        ("identifier", SEQUENCE_IDENTIFIER_STRUCT),
        ("position", pa.int64()),
        ("end_position", pa.int64()),
        ("coordinate_reference", COORDINATE_REFERENCE_STRUCT),
        ("description", pa.string()),
    ]
)
MOLECULAR_FORM_STRUCT = pa.struct(
    [
        ("protein_entity_key", pa.string()),
        ("transcript_entity_key", pa.string()),
        ("isoform_identifier", SEQUENCE_IDENTIFIER_STRUCT),
        ("sequence_identifiers", pa.list_(SEQUENCE_IDENTIFIER_STRUCT)),
        ("modifications", pa.list_(MODIFICATION_STRUCT)),
        ("variants", pa.list_(VARIANT_STRUCT)),
        ("regions", pa.list_(REGION_STRUCT)),
    ]
)


def _mapping(value: Any) -> Mapping:
    if hasattr(value, "_asdict"):
        return value._asdict()
    if isinstance(value, Mapping):
        return value
    raise TypeError(f"Molecular form values must be mappings or named tuples: {type(value)}")


def _identifier(value: Any) -> dict | None:
    if value is None:
        return None
    data = _mapping(value)
    ns = data.get("ns", data.get("type"))
    identifier = data.get("id", data.get("value"))
    if not ns or identifier is None or not str(identifier).strip():
        raise ValueError("A molecular identifier needs both a namespace and an identifier")
    return {"ns": str(ns).strip(), "id": str(identifier).strip()}


def _coordinate(value: Any) -> dict | None:
    if value is None:
        return None
    data = _mapping(value)
    system = data.get("coordinate_system") or "unknown"
    if system not in {"protein", "transcript", "genomic", "unknown"}:
        raise ValueError(f"Unsupported coordinate system: {system!r}")
    base = data.get("position_base")
    if base is not None and (type(base) is not int or base not in {0, 1}):
        raise ValueError("position_base must be 0, 1 or unspecified")
    return {
        "identifier": _identifier(data.get("identifier")),
        "coordinate_system": system,
        "position_base": base,
    }


def _features(values: Any, schema: pa.StructType) -> list[dict] | None:
    result = []
    for value in values or []:
        data = _mapping(value)
        unknown = set(data) - {field.name for field in schema}
        if unknown:
            raise ValueError(f"Unknown molecular feature fields: {sorted(unknown)}")
        item = {field.name: data.get(field.name) for field in schema}
        for name in ("position", "end_position"):
            position = item[name]
            if position is not None and (type(position) is not int or position < 0):
                raise ValueError(f"{name} must be a non-negative integer or unspecified")
        if (
            item["position"] is not None
            and item["end_position"] is not None
            and item["end_position"] < item["position"]
        ):
            raise ValueError("end_position must not precede position")
        if "identifier" in item:
            item["identifier"] = _identifier(item["identifier"])
        item["coordinate_reference"] = _coordinate(item["coordinate_reference"])
        if item["coordinate_reference"] and item["coordinate_reference"]["position_base"] == 1:
            if item["position"] == 0 or item["end_position"] == 0:
                raise ValueError("One-based coordinates must be positive")
        # A position with no supplied reference explicitly keeps unknown context.
        if item["coordinate_reference"] is None and (
            item["position"] is not None or item["end_position"] is not None
        ):
            item["coordinate_reference"] = _coordinate({})
        result.append(item)
    return result or None


def normalize_molecular_form(value: Any, *, allow_resolved: bool = True) -> dict | None:
    """Return an Arrow/JSON-ready copy without inventing molecular specificity.

    Input builders use ``allow_resolved=False``: a source cannot assign product
    entity keys. Feature order and pairings are retained. Unknown fields fail
    loudly rather than losing scientifically meaningful source information.
    """
    if value is None:
        return None
    data = _mapping(value)
    unknown = set(data) - {field.name for field in MOLECULAR_FORM_STRUCT}
    if unknown:
        raise ValueError(f"Unknown molecular form fields: {sorted(unknown)}")
    result = {field.name: data.get(field.name) for field in MOLECULAR_FORM_STRUCT}
    if not allow_resolved and (result["protein_entity_key"] or result["transcript_entity_key"]):
        raise ValueError("Product entity keys may only be assigned by resolution")
    result["isoform_identifier"] = _identifier(result["isoform_identifier"])
    result["sequence_identifiers"] = [
        _identifier(item) for item in result["sequence_identifiers"] or []
    ] or None
    result["modifications"] = _features(result["modifications"], MODIFICATION_STRUCT)
    result["variants"] = _features(result["variants"], VARIANT_STRUCT)
    result["regions"] = _features(result["regions"], REGION_STRUCT)
    return result if any(item is not None for item in result.values()) else None


def molecular_form_from_identifiers(identifiers: Any) -> dict | None:
    """Retain asserted product/sequence identifiers before normalization.

    The caller supplies participant identifiers, not an arbitrary catalogue
    cross-reference collection. Multiple isoforms are preserved as sequences;
    the singular isoform field is populated only for a unique assertion.
    """
    sequences = []
    isoforms = []
    for value in identifiers or []:
        identifier = _identifier(value)
        ns, accession = identifier["ns"], identifier["id"]
        if ns in {"uniprot", "uniprotkb"}:
            if not re.fullmatch(
                r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
                r"(?:-\d+)?(?:-PRO_\d+)?",
                accession,
            ):
                continue
            identifier["ns"] = "uniprot"
            isoform = re.fullmatch(r"(.+-\d+)(?:-PRO_\d+)?", accession)
            if isoform:
                isoforms.append({"ns": "uniprot", "id": isoform[1]})
            elif "-PRO_" not in accession:
                # An entry accession does not assert its canonical sequence.
                continue
        elif ns.startswith("refseq"):
            if not re.fullmatch(r"(?:AP|NP|XP|YP|WP|ZP|NM|XM|NR|XR)_[0-9]+(?:\.\d+)?", accession):
                continue
        elif ns in {"ensembl", "ensp", "enst", "ensembl_protein", "ensembl_transcript"}:
            if not re.fullmatch(r"ENS[A-Z]*[PT]\d+(?:\.\d+)?", accession):
                continue
        else:
            continue
        if identifier not in sequences:
            sequences.append(identifier)
    unique_isoforms = [item for i, item in enumerate(isoforms) if item not in isoforms[:i]]
    return normalize_molecular_form(
        {
            "isoform_identifier": unique_isoforms[0] if len(unique_isoforms) == 1 else None,
            "sequence_identifiers": sequences,
        },
        allow_resolved=False,
    )
