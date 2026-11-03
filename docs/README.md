# Documentation

Start with the [core documentation](../core_documentation/README.md): schema,
entity resolution, PostgreSQL, subsets, decisions and glossary. Then see the
[repository README](../README.md) and [architecture](architecture.md).

## Current guides

- [Build pipeline](../packages/build/README.md)
- [Identity layer](../packages/build/REFERENCE.md)
- [Native resolution](../packages/resolver/README.md)
- [Identity layer spec](identity-layer-spec.md)
- [API](../packages/api/README.md)
- [Python client](../packages/query/python/README.md)
- [Web explorer](../packages/web/README.md)
- [PostgreSQL](../packages/postgres/README.md)
- [Parquet-to-PostgreSQL column map](postgres-parquet-column-map.md)
- [Subset products](../packages/subsets/README.md)
- [Versioning](../VERSIONING.md)
- [Serving deployment](../deploy/README.md)

## Molecular forms

- [Entity resolution and molecular forms: conversation decisions](entity-resolution-and-molecular-forms.md)
- [Migration-branch integration and validation](reports/molecular-parquet-migration-integration-20261005.md)
- [Earlier isolated pilot and measurements](reports/molecular-form-implementation-20261004.md)
- [Resource input coverage](reports/molecular-form-input-coverage-20261004.md)

Dated reports record specific experiments and deployments. Their paths, producer
revisions and numbers describe those runs. Earlier exploratory artifacts are
listed in [history](history/README.md); the numbered [pipeline notes](pipeline/README.md)
describe the previous implementation, and [reference methods](reference-resolver.md)
and the [compiler semantics brief](reference-compiler-semantics.md) describe the
previous compact reference that the identity layer replaced; its code has been removed.
