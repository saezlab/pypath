"""Resource 3-name model: slug / short / full + forgiving filter resolution (Milestone M).

Every resource carries three names:

- **slug**  — all-lowercase canonical identifier & filter key; no ``_``, no spaces.
- **short** — the resource's own exact spelling (display name); no ``_``, no spaces.
- **full**  — the long name; spaces allowed; no ``_``.

``_`` is reserved for ``primary_secondary`` composite labels (e.g.
``PhosphoSite_SIGNOR``) and must not appear in any single resource's names.

The curated registry below (:data:`RESOURCE_NAMES`) is the reviewed audit
artifact (FR-049) — short = the resource's self-spelling, full reconciled from
the resource documentation / the legacy ``pypath.resources.descriptions``
(``label`` / ``full_name``). Resources not in the registry derive a slug from
their ``name`` and fall back to ``name`` for short/full; the build-time validator
flags any rule violation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

__all__ = [
    'ResourceNames',
    'RESOURCE_NAMES',
    'slugify',
    'validate_resource_name',
    'resolve_names',
    'resolve_filter',
    'build_filter_index',
]

_SLUG_RE = re.compile(r'^[a-z0-9]+$')
_NO_UNDERSCORE_SPACE_RE = re.compile(r'^[^_\s]+$')


@dataclass(frozen=True)
class ResourceNames:
    slug: str
    short: str
    full: str
    synonyms: tuple[str, ...] = ()


def slugify(name: str) -> str:
    """Derive a canonical slug: lowercase, keep only ``[a-z0-9]``."""
    return re.sub(r'[^a-z0-9]', '', (name or '').lower())


def validate_resource_name(slug: str, short: str, full: str) -> list[str]:
    """Return a list of rule violations (empty == valid)."""
    errors: list[str] = []
    if not slug or not _SLUG_RE.match(slug):
        errors.append(f'slug {slug!r} must be lowercase alphanumeric (no _, no spaces)')
    if not short or not _NO_UNDERSCORE_SPACE_RE.match(short):
        errors.append(f'short {short!r} must contain no _ and no spaces')
    if not full or '_' in full:
        errors.append(f'full {full!r} must contain no _')
    return errors


# Curated audit registry (slug → names). The reviewed source of truth; extend as
# the audit reconciles every resource. short = self-spelling; full = long name.
RESOURCE_NAMES: dict[str, ResourceNames] = {
    n.slug: n
    for n in (
        ResourceNames('ramp', 'RaMP', 'Relational Database for Metabolomic Pathways',
                      ('RaMP-DB', 'rampdb')),
        ResourceNames('hmdb', 'HMDB', 'Human Metabolome Database', ()),
        ResourceNames('chebi', 'ChEBI', 'Chemical Entities of Biological Interest', ()),
        ResourceNames('chembl', 'ChEMBL', 'ChEMBL', ()),
        ResourceNames('lipidmaps', 'LipidMaps', 'LIPID Metabolites And Pathways Strategy',
                      ('LIPID MAPS', 'lipid-maps')),
        ResourceNames('swisslipids', 'SwissLipids', 'SwissLipids', ()),
        ResourceNames('pubchem', 'PubChem', 'PubChem', ()),
        ResourceNames('uniprot', 'UniProt', 'Universal Protein Resource', ()),
        ResourceNames('signor', 'SIGNOR',
                      'Signaling Network Open Resource', ()),
        ResourceNames('intact', 'IntAct', 'IntAct Molecular Interaction Database', ()),
        ResourceNames('reactome', 'Reactome', 'Reactome', ()),
        ResourceNames('bindingdb', 'BindingDB', 'BindingDB', ()),
        ResourceNames('cellphonedb', 'CellPhoneDB', 'CellPhoneDB', ()),
        ResourceNames('wikipathways', 'WikiPathways', 'WikiPathways', ()),
        ResourceNames('stitch', 'STITCH',
                      'Search Tool for Interactions of Chemicals', ()),
        ResourceNames('mrclinksdb', 'MRClinksDB', 'MRClinksDB', ('mrclinkdb',)),
        ResourceNames('guidetopharma', 'GuideToPharmacology',
                      'IUPHAR Guide to Pharmacology', ('guide2pharma', 'gtopdb', 'guidetopharmacology')),
        ResourceNames('cellinker', 'Cellinker', 'Cellinker', ()),
        ResourceNames('tcdb', 'TCDB', 'Transporter Classification Database', ()),
        ResourceNames('recon3d', 'Recon3D', 'Recon3D', ()),
        ResourceNames('kegg', 'KEGG', 'Kyoto Encyclopedia of Genes and Genomes', ()),
        ResourceNames('ptfi', 'PTFI', 'Periodic Table of Food Initiative', ('ptfidiscover',)),
        ResourceNames('metatlas', 'MetAtlas', 'Metabolic Atlas', ()),
        ResourceNames('foodb', 'FooDB', 'FooDB', ()),
        ResourceNames('refmet', 'RefMet', 'Reference list of Metabolite names', ()),
        ResourceNames('psimi', 'PSI-MI', 'PSI-MI Controlled Vocabulary',
                      ('psi_mi', 'psi-mi')),
        ResourceNames('phenolexplorer', 'Phenol-Explorer',
                      'Phenol-Explorer Database', ('phenol_explorer',)),
        ResourceNames('omnipathontology', 'OmniPath', 'OmniPath Ontology',
                      ('omnipath_ontology',)),
    )
}


def resolve_names(
    *,
    name: str,
    slug: str | None = None,
    short: str | None = None,
    full: str | None = None,
    synonyms: tuple[str, ...] = (),
) -> ResourceNames:
    """Resolve the 3 names: explicit fields → curated registry → derived from ``name``."""
    derived_slug = slug or slugify(name)
    curated = RESOURCE_NAMES.get(derived_slug)
    resolved_slug = derived_slug
    resolved_short = short or (curated.short if curated else name)
    resolved_full = full or (curated.full if curated else name)
    merged_synonyms = tuple(
        dict.fromkeys(
            (*synonyms, *(curated.synonyms if curated else ()))
        )
    )
    return ResourceNames(
        slug=resolved_slug,
        short=resolved_short,
        full=resolved_full,
        synonyms=merged_synonyms,
    )


def build_filter_index(
    resources: list[ResourceNames] | None = None,
) -> dict[str, str]:
    """Map every {slug, slugified short, slugified synonym} → canonical slug."""
    index: dict[str, str] = {}
    for names in (resources if resources is not None else RESOURCE_NAMES.values()):
        index[names.slug] = names.slug
        index[slugify(names.short)] = names.slug
        for synonym in names.synonyms:
            index[slugify(synonym)] = names.slug
    return index


def resolve_filter(
    query: str,
    index: dict[str, str] | None = None,
) -> str | None:
    """Resolve a resource filter argument (slug/short/synonym, case-insensitive) → slug."""
    if not query:
        return None
    lookup = index if index is not None else build_filter_index()
    return lookup.get(slugify(query))
