"""The filter and cursor contracts of the queries."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

# Unknown fields pass through, as in the API's models.
_CONFIG = ConfigDict(extra="allow", populate_by_name=True)


class EntitySearchCursor(BaseModel):
    """Where a search continues: its match phase and the entities already returned."""

    model_config = _CONFIG

    phase: Literal["all", "key", "prefix", "contains"]
    offset: int = Field(ge=0)
    queryKey: str


class SearchFilters(BaseModel):
    """One internal filter contract, accepting historical wire aliases once."""

    model_config = _CONFIG

    reference_entity_keys: list[str] = Field(default_factory=list)
    protein_entity_keys: list[str] = Field(default_factory=list)
    transcript_entity_keys: list[str] = Field(default_factory=list)
    isoform_identifiers: list[str] = Field(
        default_factory=list, description="Exact namespace:id; unspecified forms do not match."
    )
    molecular_endpoint_mode: Literal["any", "both", "source", "target"] = "any"
    entity_ids: list[str] = Field(default_factory=list)
    entity_pks: list[str] = Field(default_factory=list)
    scope_entity_ids: list[str] = Field(default_factory=list)
    entity_types: list[str] = Field(default_factory=list)
    sources: list[str] = Field(default_factory=list)
    taxonomy_ids: list[str] = Field(default_factory=list)
    predicates: list[str] = Field(default_factory=list)
    relation_categories: list[str] = Field(default_factory=list)
    object_aspect_qualifier: list[str] = Field(default_factory=list)
    object_direction_qualifier: list[str] = Field(default_factory=list)
    causal_mechanism_qualifier: list[str] = Field(default_factory=list)
    interaction_types: list[str] = Field(default_factory=list)
    # With interaction_types: the type at the other end of a relation, so that
    # interaction_types=[protein] and partner_types=[chemical] is protein-chemical.
    partner_types: list[str] = Field(default_factory=list)
    annotation_term_ids: list[str] = Field(default_factory=list)
    gene_mode: Literal["any", "both", "source", "source_only", "target", "target_only"] = "any"
    scope_endpoint_mode: Literal[
        "any", "both", "source", "source_only", "target", "target_only"
    ] = "any"
    is_directed: bool | None = None
    sign: int | None = None
    signs: list[int] = Field(default_factory=list)
    min_evidence_count: int | None = Field(default=None, ge=0)

    @model_validator(mode="before")
    @classmethod
    def normalize_aliases(cls, value: Any) -> Any:
        if not isinstance(value, dict):
            return value
        data = dict(value)
        aliases = {
            "entity_ids": ("entities", "genes", "entityIds"),
            "entity_pks": ("entityPks",),
            "entity_types": ("entityTypes",),
            "sources": ("resources",),
            "taxonomy_ids": ("taxonomyIds", "ncbi_tax_id"),
            "relation_categories": ("categories",),
            "interaction_types": ("interactionTypes",),
            "partner_types": ("partnerTypes",),
            "annotation_term_ids": ("annotationTermIds",),
        }
        for canonical, legacy in aliases.items():
            merged = []
            for key in (canonical, *legacy):
                values = data.pop(key, [])
                if values is None:
                    values = []
                if not isinstance(values, list):
                    raise ValueError(f"{key} must be a list")
                for item in values:
                    if not isinstance(item, (str, int)) or isinstance(item, bool):
                        raise ValueError(f"{key} must contain strings or integers")
                    text = str(item).strip()
                    if text and text not in merged:
                        merged.append(text)
            data[canonical] = merged
        if "endpointMode" in data:
            mode = data.pop("endpointMode")
            data.setdefault("gene_mode", mode)
            data.setdefault("scope_endpoint_mode", mode)
        elif "scope_endpoint_mode" in data:
            data.setdefault("gene_mode", data["scope_endpoint_mode"])
        return data


def normalize_filters(value: SearchFilters | dict[str, Any] | None) -> dict[str, Any]:
    model = value if isinstance(value, SearchFilters) else SearchFilters.model_validate(value or {})
    return model.model_dump()
