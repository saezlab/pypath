"""Name ranking for choosing one display name among several labels (a chemical
structure group's members). Lexical ties are deterministic."""

import re

NAME_NAMESPACES = ("name", "om:0202", "entry_name", "entry name")
NAME_RULES = (
    (r"(?i)^https?://", 100),
    (r"(?i)^(MLS|SMR|cid_|ZINC|SID_|CID_|InChI=)", 100),
    (r"^[A-Z0-9._-]{1,3}$", 10),
    (r"^[A-Z][0-9A-Z]{3,}$", 8),
    (r"^[A-Z][a-z]", -3),
)
CHEMICAL_TYPES = {
    "chemicalentity",
    "chemical",
    "smallmolecule",
    "compound",
    "metabolite",
    "drug",
    "lipid",
}


def name_rank(name):
    name = name.strip()
    score = sum(penalty for pattern, penalty in NAME_RULES if re.search(pattern, name))
    score += 25 * (len(name) > 80) + 8 * (len(name) > 40) + 4 * (len(name) < 4)
    return score, len(name), name


def preferred_name(names):
    names = {name.strip() for name in names if name and name.strip()}
    return min(names, key=name_rank) if names else None
