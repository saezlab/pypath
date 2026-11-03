# Duplicated Identifier Values Within Entities

## Summary

`data/combined/entity.parquet` contains no exact duplicate identifier records for the same entity when comparing `(identifier_type, identifier)` pairs.

However, the same identifier string is often attached to the same entity multiple times under different identifier types.

## Observed counts

- Total entities: `409,568`
- Entities with exact duplicate `(identifier_type, identifier)` pairs: `0`
- Entities where the same identifier value appears under multiple types: `36,619`
- Distinct repeated identifier values within entities: `50,695`
- Extra repeated occurrences: `56,111`

## Common patterns

Most duplicated values are caused by inconsistent assignment of semantic identifier/name types, for example:

```text
OM:0202:Name + OM:0203:Synonym
MI:1097:Uniprot + OM:0007:Signor
OM:0201:Gene Name Synonym + OM:0203:Synonym
OM:0200:Gene Name Primary + OM:0202:Name
OM:0200:Gene Name Primary + OM:0202:Name + OM:0203:Synonym
```

Example:

```text
entity_pk: A0A087WXS9
identifier value: TBC1D3I
assigned as:
  - OM:0200:Gene Name Primary
  - OM:0202:Name
```

## Likely cause

This should probably be fixed in the input modules where identifier types are assigned. Different inputs appear to classify the same value inconsistently, e.g. as a primary gene name, display name, synonym, source-specific identifier, etc.

## Proposed direction

Make identifier type assignment more consistent before entity merging/combination, ideally by defining clearer rules for:

- canonical source identifiers vs source cross-references
- primary names vs display names
- gene name synonyms vs generic synonyms
- source-specific identifiers such as SIGNOR vs UniProt accessions

The combined output should avoid adding the same identifier string multiple times to one entity unless the distinction between identifier types is intentional and useful.
