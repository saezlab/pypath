"""
One small release for the metabolite class product, written straight to Parquet.

Each entity stands for one case of the rules: cortisol in HMDB and its 17S stereoisomer in
KEGG, lactate in HMDB and its charged form in a Rhea reaction, a ChEMBL compound that no
source claims and one that only cross-references HMDB, a native fallback and a lipid
without a structure, a molecule only a reaction knows, and one that only the ChEBI
metabolite role claims. Entity keys are readable names.
"""

from __future__ import annotations

from pathlib import Path

from parquet_release import entity_row, relation_row, write_resources

CORTISOL_KEY = 'JYGXADMDTFJGBT-VWUMJDOOSA-N'
CORTISOL_17S_KEY = 'JYGXADMDTFJGBT-CZFMHFDVSA-N'
LACTATE_KEY = 'JVTAAEKCZFNVCJ-UHFFFAOYSA-N'
LACTATE_ION_KEY = 'JVTAAEKCZFNVCJ-UHFFFAOYSA-M'
PYRUVATE_KEY = 'LCTONWCANYUPML-UHFFFAOYSA-N'
MENTHOL_KEY = 'NOOLISFMXDJSKH-UHFFFAOYSA-N'
COMPOUND_KEY = 'AAAAAAAAAAAAAA-UHFFFAOYSA-N'
XREF_KEY = 'BBBBBBBBBBBBBB-UHFFFAOYSA-N'

CHEMICAL = 'chemical_entity'


def _molecule(key, inchikey, label, identifiers = ()):
    """
    A chemical identified by its InChIKey, in the connectivity group of that key.
    """

    return entity_row(
        key,
        CHEMICAL,
        'inchikey',
        inchikey,
        label,
        [
            ('inchikey', inchikey, True),
            *((ns, value, False) for ns, value in identifiers),
        ],
        inchikey,
    )


def _term(key, namespace, identifier, label, entity_type = CHEMICAL):
    """
    An entity without a structure key.
    """

    return entity_row(
        key,
        entity_type,
        namespace,
        identifier,
        label,
        [(namespace, identifier, True)],
    )


def _relation(number, source, subject, predicate, obj, *, kind = 'relation', **types):
    """
    A relation of ``source``, keyed by the source and a number.
    """

    return relation_row(
        f'{source}:{number}',
        subject,
        predicate,
        obj,
        kind = kind,
        source = source,
        **types
    )


def resources() -> dict:
    """
    The entities and relations of each resource of the release.
    """

    cortisol_ids = [
        ('hmdb', 'HMDB00063'),
        ('hmdb', 'HMDB0000063'),
        ('chebi', 'CHEBI:17650'),
    ]
    reaction = {'subject_type': 'molecular_activity', 'object_type': CHEMICAL}
    ontology = {'kind': 'ontology', 'subject_type': CHEMICAL, 'object_type': CHEMICAL}
    protein = {'object_type': 'protein'}

    return {
        'hmdb': (
            [
                _molecule('cortisol', CORTISOL_KEY, 'Cortisol', cortisol_ids),
                _molecule(
                    'lactate',
                    LACTATE_KEY,
                    'L-Lactic acid',
                    [('hmdb', 'HMDB0000190')],
                ),
                _term('lipid', 'goslin', 'PC 34:1', 'PC(34:1)'),
                # HMDB's own disease terms are not chemicals and are not claimed.
                _term('hmdb_disease', 'hmdb_disease', 'D1', 'A disease', 'ontology_class'),
            ],
            [],
        ),
        'chebi': (
            [
                _molecule(
                    'cortisol',
                    CORTISOL_KEY,
                    '11beta,17alpha,21-trihydroxypregn-4-ene-3,20-dione',
                    cortisol_ids,
                ),
                _molecule('menthol', MENTHOL_KEY, 'menthol'),
                _term('role_metabolite', 'chebi', 'CHEBI:25212', 'metabolite'),
                _term('role_human_metabolite', 'chebi', 'CHEBI:77746', 'human metabolite'),
            ],
            [
                _relation(
                    1,
                    'chebi',
                    'role_human_metabolite',
                    'subclass_of',
                    'role_metabolite',
                    **ontology
                ),
                _relation(
                    2,
                    'chebi',
                    'menthol',
                    'has_chemical_role',
                    'role_human_metabolite',
                    **ontology
                ),
                _relation(
                    3,
                    'chebi',
                    'cortisol',
                    'has_chemical_role',
                    'role_metabolite',
                    **ontology
                ),
            ],
        ),
        'kegg': (
            [_molecule('cortisol_17s', CORTISOL_17S_KEY, 'Cortisol (17S)')],
            [],
        ),
        'rhea': (
            [
                _molecule('lactate_ion', LACTATE_ION_KEY, '(S)-lactate'),
                _molecule('pyruvate', PYRUVATE_KEY, 'pyruvate'),
                _term(
                    'rhea_reaction',
                    'rhea',
                    'RHEA:1',
                    'lactate dehydrogenase',
                    'molecular_activity',
                ),
            ],
            [
                _relation(1, 'rhea', 'rhea_reaction', 'has_input', 'pyruvate', **reaction),
                _relation(2, 'rhea', 'rhea_reaction', 'has_output', 'lactate_ion', **reaction),
            ],
        ),
        'metatlas': (
            [
                # Resolution abstained: a native fallback with no structure key.
                _term('fallback', 'chebi', 'CHEBI:99999', 'unresolved metabolite'),
                _term(
                    'metatlas_reaction',
                    'metatlas',
                    'MAR1',
                    'a reaction',
                    'molecular_activity',
                ),
            ],
            [
                _relation(
                    1,
                    'metatlas',
                    'metatlas_reaction',
                    'has_input',
                    'fallback',
                    **reaction
                ),
            ],
        ),
        'chembl': (
            [
                _molecule('screening_compound', COMPOUND_KEY, 'CHEMBL1'),
                # Another resource's cross-reference to HMDB does not make an HMDB molecule.
                _molecule('xref_compound', XREF_KEY, 'CHEMBL2', [('hmdb', 'HMDB0099999')]),
                _term('receptor', 'uniprot', 'P04150', 'NR3C1', 'protein'),
            ],
            [
                _relation(1, 'chembl', 'cortisol', 'interacts_with', 'receptor', **protein),
                _relation(
                    2,
                    'chembl',
                    'screening_compound',
                    'interacts_with',
                    'receptor',
                    **protein
                ),
            ],
        ),
        # A resource that sorts before ChEBI gives the release label of menthol.
        'bindingdb': (
            [_molecule('menthol', MENTHOL_KEY, 'MENTHOL (BindingDB)')],
            [],
        ),
    }


ALL_SIGNALS = """\
classes:
  - {name: hmdb,              precedence: 1}
  - {name: chebi_role,        precedence: 2}
  - {name: metabolic_network, precedence: 3}
  - {name: none,              precedence: 4}
default: none
sources:
  - {resource: hmdb, priority: 1, class: hmdb}
  - {chebi_role_root: 'CHEBI:25212', priority: 2, class: chebi_role}
  - {predicates: [has_input, has_output], priority: 3, class: metabolic_network,
     resources: [kegg, rhea, recon3d, metatlas]}
"""


def write_release(root: Path, version: str = '2026.11') -> str:
    """
    Write the release under ``root``.

    Args:
        root: The data root.
        version: The version of the release.

    Returns:
        The version of the release.
    """

    return write_resources(root, resources(), version)


def write_rules(directory: Path, text: str = ALL_SIGNALS) -> Path:
    """
    Write a rule file into ``directory``; returns its path.
    """

    path = Path(directory) / 'metabolite_class.yaml'
    path.write_text(text)

    return path
