"""Database build presets for common deployment scenarios."""

from __future__ import annotations

# Organism groups
HUMAN = [9606]
CORE = [9606, 10090]  # human + mouse

# --- Build-scope tiers (007-resolution-coverage) ---
# Nested: only-human ⊂ core-model ⊂ extended-model ⊂ model-organisms ⊂ complete.
# Each tier is a strict superset of the previous one.
CORE_MODEL = [9606, 10090, 10116]  # human, mouse, rat
EXTENDED_MODEL = CORE_MODEL + [
    7955,   # zebrafish (Danio rerio)
    7227,   # fruit fly (Drosophila melanogaster)
    6239,   # C. elegans
]
# Full common model-organism set (007 clarification). Superset of EXTENDED_MODEL.
MODEL_ORGANISMS = EXTENDED_MODEL + [
    9823,     # pig (Sus scrofa)
    9544,     # rhesus macaque (Macaca mulatta)
    8364,     # Xenopus tropicalis
    4932,     # budding yeast (S. cerevisiae)
    511145,   # E. coli str. K-12 substr. MG1655 (NCBI Gene reference)
    9615,     # dog (Canis lupus familiaris)
]

# Named build scopes → organism taxid list (None = complete: all organisms).
SCOPES = {
    'only-human': HUMAN,
    'core-model': CORE_MODEL,
    'extended-model': EXTENDED_MODEL,
    'model-organisms': MODEL_ORGANISMS,
    'complete': None,
}

# Common-name → NCBI taxid, so an explicit --scope list may use names
# (e.g. "human,pig,chimpanzee") as well as numeric taxids.
COMMON_NAMES = {
    'human': 9606, 'mouse': 10090, 'rat': 10116, 'pig': 9823,
    'rhesus': 9544, 'macaque': 9544, 'xenopus': 8364, 'frog': 8364,
    'zebrafish': 7955, 'fruitfly': 7227, 'fly': 7227, 'drosophila': 7227,
    'worm': 6239, 'celegans': 6239, 'yeast': 4932, 'ecoli': 511145,
    'dog': 9615, 'chimpanzee': 9598, 'chimp': 9598, 'cow': 9913,
    'chicken': 9031,
}

# Mapping pair definitions
# Each tuple: (source_type, target_type, backend)

PROTEIN_CORE = [
    ('genesymbol', 'uniprot', 'uniprot'),
    ('entrez', 'uniprot', 'uniprot'),
    ('hgnc', 'uniprot', 'uniprot'),
    ('refseqp', 'uniprot', 'uniprot'),
    # SwissProt/TrEMBL for cleanup pipeline
    ('genesymbol', 'swissprot', 'uniprot'),
    ('trembl', 'genesymbol', 'uniprot'),
    # Synonyms
    ('genesymbol-syn', 'uniprot', 'uniprot'),
]

ENSEMBL = [
    ('ensg', 'genesymbol', 'biomart'),
    ('ensp', 'ensg', 'biomart'),
    ('enst', 'ensg', 'biomart'),
    ('ensp', 'uniprot', 'biomart'),
    ('ensg', 'uniprot', 'biomart'),
]

# Ensembl Genomes (non-vertebrate divisions: plants/fungi/metazoa/protists).
# Same biomart backend, but the division-typed ensgg/ensgp/ensgt id-types so the
# Ensembl Genomes stable IDs stay distinguishable from the vertebrate ensg/ensp.
# Run only for organisms the BioMart backend routes to a genomes division
# (see ORGANISM_DIVISION in mapping/backends/_biomart.py).
ENSEMBL_GENOMES = [
    ('ensgg', 'genesymbol', 'biomart'),
    ('ensgp', 'ensgg', 'biomart'),
    ('ensgt', 'ensgg', 'biomart'),
    ('ensgp', 'uniprot', 'biomart'),
    ('ensgg', 'uniprot', 'biomart'),
]

METABOLITE = [
    # Metabolite mappings are organism-agnostic (ncbi_tax_id=0) and handled
    # by populate_metabolites() which auto-discovers all pairs from
    # UniChem, RaMP, MetaNetX, and BiGG backends.
]

MIRNA = [
    # miRBase -- organism-specific
    # Handled by mirbase backend, not standard mapping pairs
]


# COSMOS PKN-specific protein mappings (cleanup pipeline)
COSMOS_CLEANUP = [
    # Secondary -> primary UniProt AC resolution is loaded organism-agnostically
    # (tax 0) by DatabaseBuilder.load_uniprot_sec_ac() for every protein build
    # (ADR 0006), not as a per-organism mapping pair here.
    # UniProt -> gene symbol (for TrEMBL -> SwissProt resolution)
    ('uniprot', 'genesymbol', 'uniprot'),
    # SwissProt lookup (already in PROTEIN_CORE via genesymbol -> swissprot)
]

# Preset definitions
PRESETS = {
    'minimal': {
        'description': 'Basic human protein ID mappings',
        'organisms': HUMAN,
        'mappings': PROTEIN_CORE + ENSEMBL,
        'metabolite': False,
        'mirna': False,
        'orthology': False,
        'reflists': True,
    },

    'chemical': {
        'description': (
            'Complete chemical stack: database IDs + names + structures '
            '(InChI/SMILES), no protein/gene/miRNA'
        ),
        'organisms': HUMAN,   # chemicals are organism-agnostic (rows at tax 0)
        'mappings': [],       # no protein/gene mappings
        'metabolite': True,   # UniChem/RaMP/MetaNetX/BiGG/structures + names
        'mirna': False,
        'orthology': False,
        'reflists': False,
    },

    'cosmos': {
        'description': 'All mappings needed by COSMOS PKN build (human + mouse)',
        'organisms': CORE,
        'mappings': PROTEIN_CORE + ENSEMBL + COSMOS_CLEANUP,
        'metabolite': True,
        'mirna': False,
        'orthology': True,
        'reflists': True,
    },
    'standard': {
        'description': 'Human + mouse protein mappings with metabolites',
        'organisms': CORE,
        'mappings': PROTEIN_CORE + ENSEMBL,
        'metabolite': True,
        'mirna': False,
        'orthology': True,
        'reflists': True,
    },
    'model': {
        'description': 'All model organisms + metabolites + miRNA + orthology',
        'organisms': MODEL_ORGANISMS,
        'mappings': PROTEIN_CORE + ENSEMBL,
        'metabolite': True,
        'mirna': True,
        'orthology': True,
        'reflists': True,
    },
    'full': {
        'description': 'All organisms from UniProt FTP + everything',
        'organisms': None,  # None means use FTP for all
        'mappings': PROTEIN_CORE + ENSEMBL,
        'metabolite': True,
        'mirna': True,
        'orthology': True,
        'reflists': True,
        'ftp': True,
    },
}

# Common Parquet export tables (source_type, target_type, ncbi_tax_id)
# These are the most frequently requested via the API
PARQUET_TABLES = {
    'minimal': [
        ('genesymbol', 'uniprot', 9606),
        ('entrez', 'uniprot', 9606),
        ('uniprot', 'genesymbol', 9606),
    ],
    'standard': [
        ('genesymbol', 'uniprot', 9606),
        ('entrez', 'uniprot', 9606),
        ('uniprot', 'genesymbol', 9606),
        ('genesymbol', 'uniprot', 10090),
        ('entrez', 'uniprot', 10090),
        ('ensg', 'genesymbol', 9606),
        ('ensg', 'genesymbol', 10090),
    ],
    'model': [
        ('genesymbol', 'uniprot', org)
        for org in MODEL_ORGANISMS
    ] + [
        ('entrez', 'uniprot', org)
        for org in MODEL_ORGANISMS
    ] + [
        ('ensg', 'genesymbol', org)
        for org in MODEL_ORGANISMS
    ],
}
PARQUET_TABLES['full'] = PARQUET_TABLES['model']


def resolve_scope(value):
    """Resolve a ``--scope`` value to an organism taxid list.

    Accepts either a named tier (see :data:`SCOPES` —
    ``only-human`` / ``core-model`` / ``extended-model`` /
    ``model-organisms`` / ``complete``) or an explicit comma/semicolon
    separated list mixing NCBI taxonomy ids and common organism names,
    e.g. ``"9606,9823,9598"`` or ``"human,pig,chimpanzee"``.

    Returns a list of taxids, or ``None`` for ``complete`` (all organisms /
    FTP-all). Raises ``ValueError`` on an unrecognised token.
    """
    if value is None:
        return None

    key = value.strip().lower()

    if key in SCOPES:
        return SCOPES[key]

    taxa: list[int] = []

    for token in value.replace(';', ',').split(','):
        token = token.strip()

        if not token:
            continue

        if token.isdigit():
            taxa.append(int(token))
            continue

        name = ''.join(ch for ch in token.lower() if ch.isalnum())

        if name in COMMON_NAMES:
            taxa.append(COMMON_NAMES[name])
            continue

        raise ValueError(
            f'Invalid --scope {value!r}: expected a preset '
            f'({", ".join(SCOPES)}) or a comma-separated list of NCBI '
            f'taxonomy ids / known organism names; got unrecognised '
            f'token {token!r}.'
        )

    if not taxa:
        raise ValueError(f'Invalid --scope {value!r}: no organisms resolved.')

    # De-duplicate while preserving order.
    seen: set[int] = set()
    return [t for t in taxa if not (t in seen or seen.add(t))]
