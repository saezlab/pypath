"""Recognize sequence-archive accessions in nominal UniProt protein columns."""
import re
from omnipath_core.naming import Namespace


def accession_namespace(value):
    value = str(value or '').strip()
    if re.fullmatch(r'(?:AP|NP|XP|YP|WP|ZP)_[0-9]+(?:\.[0-9]+)?', value):
        return Namespace.REFSEQ_PROTEIN
    if re.fullmatch(r'[A-Z]{2}_[0-9]+(?:\.[0-9]+)?', value):
        return Namespace.REFSEQ
    if re.fullmatch(r'[A-Z]{3}[0-9]{5,}(?:\.[0-9]+)?', value):
        return Namespace.GENBANK
    return Namespace.UNIPROT
