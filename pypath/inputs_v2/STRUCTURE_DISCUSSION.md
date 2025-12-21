# inputs_v2 Structure Discussion

## Snapshot of current state
- Each resource module exposes a `resource()` metadata generator and one or more data generators (e.g., `uniprot_proteins`, `guidetopharma_targets`).
- The core mapping logic is declarative via `EntityBuilder` + `IdentifiersBuilder` + `AnnotationsBuilder` in each function.
- Downloading is inlined per generator via `download_and_open`, with ad hoc flags per source.
- There is no explicit registry of resources or their outputs; discovery is mostly import-by-name via `inputs_v2.get_method`.

## Option A: Keep function-based modules, add lightweight conventions
**Idea**: Preserve the current module layout but formalize a minimal contract and helper utilities.

**What changes**
- Define a standard naming scheme and simple module-level constants, e.g. `RESOURCE_ID`, `RESOURCE_NAME`, `RESOURCE_URL`, `DEFAULT_SUBFOLDER`.
- Add helper functions in a shared module, e.g. `inputs_v2/utils.py`:
  - `resource_metadata(...) -> Entity`
  - `download_tabular(...) -> opener`
  - `iter_dict_rows(opener, delimiter='\t')`
- Encourage a common signature for generators, e.g. `(force_refresh=False, **kwargs)`.

**Pros**
- Minimal refactor; high compatibility with `get_method`.
- Very low cognitive overhead for new resource authors.
- Keeps full flexibility for odd or complex sources.

**Cons**
- More duplication remains (download flags, metadata shape).
- No explicit object to introspect; registry still implicit.
- No built-in separation between raw rows and mapped entities.

## Option B: Introduce a Resource class (your proposal)
**Idea**: Each resource is represented by a class that owns metadata, download, raw row iterators, and output generators.

**Possible shape**
```python
class Resource:
    id: ResourceCv
    name: str
    url: str
    license: LicenseCV
    update_category: UpdateCategoryCV
    pubmed: str | None
    description: str

    def metadata(self) -> Generator[Entity]:
        ...

    def download(self, *, force_refresh: bool = False):
        ...

    def raw(self, *, force_refresh: bool = False):
        """Yield raw rows (dicts or lists)."""
        ...

    def outputs(self) -> dict[str, Callable[..., Generator[Entity]]]:
        """Return named generators (e.g., targets, ligands, interactions)."""
        ...
```

**How it integrates with existing `get_method`**
- Keep module-level functions as thin wrappers:
  - `def resource(): return GuideToPharma().metadata()`
  - `def guidetopharma_targets(...): return GuideToPharma().targets(...)`
- Alternatively, add a registry and update `get_method` to resolve names to `Resource.outputs()`.

**Pros**
- Clear separation: metadata / download / raw / output mappings.
- Easier to cache or reuse shared downloads across outputs.
- A natural place to house per-resource configuration and state (file paths, columns, etc.).

**Cons**
- More code around simple sources; class boilerplate can obscure the mapping.
- You will need to decide on lifecycle (singleton vs. per-call instance).
- If multiple datasets per resource need different downloads, the class can get large.

## Option C: ResourceSpec dataclass + registry (hybrid)
**Idea**: Avoid a full class hierarchy; keep functions but centralize metadata and download config.

**What changes**
- Create a `ResourceSpec` dataclass:
  - metadata fields (resource id, license, url, description)
  - download config (url, filename, subfolder, large, ext, encoding)
- Each module exposes a `SPEC` and a set of generator functions.

**Pros**
- Minimal boilerplate; readable at top-of-file.
- Explicit metadata and download configuration can be reused by multiple generators.
- Easy to build a registry for documentation or automated checks.

**Cons**
- Not as strong a “single place” as a class with methods.
- Still relies on conventions for raw vs. mapped outputs.

## Option D: Config-driven resource declarations
**Idea**: Move most of the resource metadata and builder schema into a config file (YAML/JSON), with thin Python loaders.

**Pros**
- Non-code updates for simple mappings.
- Can enable automated validation or generation of documentation.

**Cons**
- Many resources require custom parsing (regex, special rules).
- Would create a split between declarative config and imperative Python that might be harder to debug.

## Suggested direction (pragmatic)
A hybrid of **Option B** and **Option C** tends to work well:
- Introduce a `ResourceSpec` (Option C) now to reduce duplication and standardize metadata/download config.
- For complex sources (Guide to Pharmacology, BindingDB), consider a `Resource` class that uses `SPEC` internally and exposes multiple output generators.
- Keep module-level functions as a stable API so `inputs_v2.get_method` continues to work.

This provides structure without forcing a full refactor of every module. Over time, the most complex resources can migrate to classes, while simple ones remain function-based.

## Concrete next steps if you want to prototype the Resource class
1. Create `pypath/pypath/inputs_v2/base.py` with `Resource` + `resource_metadata(...)` helper.
2. Pick a multi-output resource (Guide to Pharmacology) and wrap it in a class, leaving wrappers at module scope.
3. Add a tiny registry to map resource names to class instances for introspection and documentation.
4. Keep `get_method` unchanged initially; update it later if registry proves useful.

If you want, I can sketch a minimal `Resource` base and a `GuideToPharmaResource` skeleton next.
