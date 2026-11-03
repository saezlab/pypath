# Resource browsing metadata

Everything the resources page shows about a resource is declared by its pypath
`inputs_v2` module, in its `ResourceConfig`: the display name (`name`), a one-line
`summary`, the tags (`topics`, `content`, `knowledge_origin`), and the
license, whose `LicenseCV` term says whether academic and commercial use are allowed.
The `inputs_v2` README has the rules.

`resource_tags.yaml` is the controlled vocabulary of tags: stable IDs, labels, dimensions
(`topic`, `content`, `knowledge_origin`) and descriptions.
`omnipath_core.resource_metadata` has one enum of IDs per dimension (`Topic`, `Content`,
`KnowledgeOrigin`), and the config has one field per dimension. Tags describe what
the module imports rather than every dataset available on a provider's website.
Knowledge-origin tags mean a resource contains that kind of knowledge, not that every
imported record has that evidence. Snapshot-specific dataset subsets remain separately
visible in build metadata; these resource-level tags do not promise that every subset
includes all tagged content.

The page's Entities filter is not a tag: it lists the entity types each build recorded
in `resolution_stats.json` (`by_entity_type`), so it describes the data of that build.

`resource_metadata.yaml` is generated from the configs by
`scripts/generate_resource_metadata.py` and must not be edited by hand.
`make check-generated` fails while it is stale. The API reads it, so it needs no pypath.
Both files use JSON syntax (a YAML subset).

The loader first selects build metadata (or the generated file for legacy builds), then
adds the current display name, description and tags. Exact license, download time and
other acquisition fields remain from that selected metadata. To change browsing
metadata, edit the config, regenerate, deploy the updated core package and restart the
API. No resource rebuild is required.

License use is a property of the license, declared once per `LicenseCV` term. These are
browsing classifications, subject to the exact license conditions. `allowed` means the
recorded data license permits that category of use; `not_allowed` means it excludes it;
`requires_permission` means a separate authorization is specified; `unknown` means the
available metadata does not settle the question. A noncommercial license is not a
guarantee that every use by an academic institution is noncommercial. Software licenses
alone are left unknown for imported data. Vague commercial restrictions also stay
unknown. The loader classifies a build by its own recorded license; a missing license,
or one no current input declares, yields unknown. Review the source terms when changing
a classification.
