"""The chemistry toolkit binding (spec 011 R6/R14) -- an optional extra
(``omnipath-utils[chem]``, rdkit). Every function here degrades to ``None``
rather than raising when the toolkit is unavailable: "structure-to-anything
translation routes report themselves unavailable and say so; identifier-to-
identifier translation is unaffected" (research R14). No caller needs to
check :func:`chemistry_available` itself -- these functions already do, and
a caller that does want to report the capability explicitly (a route, a
capability-detection query) can call it directly.
"""

from __future__ import annotations

from functools import cache


@cache
def chemistry_available() -> bool:
    """Whether ``rdkit`` is importable in this process. Cached -- it's a
    static fact about the running environment, not something that changes
    mid-process.
    """

    try:
        import rdkit  # noqa: F401
    except ImportError:
        return False
    return True


def compute_inchikey(
    *, inchi: str | None = None, smiles: str | None = None,
) -> str | None:
    """The Standard InChIKey for a structure, computed from its InChI when
    given (preferred -- InChI is already canonical by construction) or
    else its SMILES. Returns ``None`` when the toolkit is unavailable,
    both inputs are empty, or the input does not parse.
    """

    if not chemistry_available():
        return None
    if not inchi and not smiles:
        return None

    from rdkit import Chem
    from rdkit.Chem import inchi as rdkit_inchi

    try:
        if inchi:
            key = rdkit_inchi.InchiToInchiKey(inchi)
            return key or None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return rdkit_inchi.MolToInchiKey(mol) or None
    except Exception:  # noqa: BLE001 -- rdkit raises a mix of exception types
        return None


def canonicalize_smiles(smiles: str) -> str | None:
    """The canonical form of a SMILES string, so two spellings of the same
    molecule compare equal (C7). Returns ``None`` when the toolkit is
    unavailable or the string does not parse.
    """

    if not chemistry_available() or not smiles:
        return None

    from rdkit import Chem

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True) or None
    except Exception:  # noqa: BLE001
        return None
