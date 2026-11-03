#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Export OmniPath CV terms as an extension of PSI-MI OBO.

This script generates a combined OBO file containing:
1. The full PSI-MI ontology (fetched from GitHub)
2. OmniPath-specific terms (OM:* accessions) appended with proper is_a relationships

OmniPath terms that extend PSI-MI categories use is_a relationships to their
MI parent terms, enabling proper hierarchy traversal across both ontologies.

Usage:
    python scripts/export_omnipath_obo.py [output_path]
    
    output_path: Optional path for the OBO file (default: omnipath_mi.obo)
"""

from __future__ import annotations

import inspect
import re
import sys
import urllib.request
from datetime import datetime
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from pypath.internals.cv_terms import CvEnum
from pypath.internals import cv_terms


# PSI-MI OBO source
PSI_MI_URL = "https://raw.githubusercontent.com/HUPO-PSI/psi-mi-CV/master/psi-mi.obo"


def fetch_psi_mi() -> str:
    """Fetch PSI-MI OBO content from GitHub.
    
    Returns:
        The PSI-MI OBO file content as a string.
    """
    print("Fetching PSI-MI OBO...")
    with urllib.request.urlopen(PSI_MI_URL) as response:
        content = response.read().decode('utf-8')
    print(f"  Fetched {len(content)} bytes")
    return content


def fix_malformed_dates(content: str) -> str:
    """Fix malformed date formats in OBO content.
    
    PSI-MI has dates in non-standard format like "date: 15:04:2021 22:57"
    (DD:MM:YYYY instead of ISO format) which causes pronto/fastobo to fail.
    
    Args:
        content: The OBO file content.
        
    Returns:
        Content with malformed dates removed.
    """
    # Remove malformed header date lines (DD:MM:YYYY format)
    content = re.sub(r'^date: \d{2}:\d{2}:\d{4}.*\n', '', content, flags=re.MULTILINE)
    # Remove all creation_date fields (parser gets corrupted by malformed header)
    content = re.sub(r'^creation_date:.*\n', '', content, flags=re.MULTILINE)
    return content


def format_name(name: str) -> str:
    """Convert ENUM_NAME to readable 'enum name' format."""
    return name.lower().replace("_", " ")


def escape_obo_string(s: str) -> str:
    """Escape special characters for OBO format."""
    if not s:
        return ""
    # Escape backslashes first, then quotes
    return s.replace("\\", "\\\\").replace('"', '\\"')


def extract_om_terms() -> list[dict]:
    """Extract all OM:* terms from cv_terms module.
    
    Returns:
        List of dicts with keys: accession, name, definition, is_a
        Note: We now use is_a for ALL parent relationships (both OM and MI).
    """
    terms = {}  # accession -> term_data
    parent_terms = {}  # Track parent terms defined as tuples
    
    # Get all CvEnum subclasses from the cv_terms module
    for _, obj in inspect.getmembers(cv_terms):
        if not (inspect.isclass(obj) and issubclass(obj, CvEnum) and obj is not CvEnum):
            continue
        
        # Get the parent term if it exists
        parent_accession = None
        if hasattr(obj, "parent_cv_term") and obj.parent_cv_term:
            parent_term = obj.parent_cv_term
            if isinstance(parent_term, tuple):
                parent_accession = parent_term[0]
                # If it's an OM parent term, store it for later
                if parent_accession.startswith("OM:") and parent_accession not in parent_terms:
                    parent_name = parent_term[1] if len(parent_term) > 1 else ""
                    parent_def = parent_term[2] if len(parent_term) > 2 else ""
                    parent_terms[parent_accession] = {
                        "accession": parent_accession,
                        "name": format_name(parent_name) if parent_name else parent_accession,
                        "definition": parent_def,
                        "is_a": "MI:0000",  # Link OM roots to MI root
                    }
            elif isinstance(parent_term, str):
                parent_accession = parent_term
        
        # Iterate through all enum members
        for member in obj:
            accession = member.value
            
            # Only include OmniPath-specific terms (OM accessions)
            if not accession.startswith("OM:"):
                continue
            
            definition = getattr(member, "definition", None) or ""
            
            # Use is_a for ALL parent relationships (both OM and MI)
            # This allows proper hierarchy traversal in the combined ontology
            # Use is_a for ALL parent relationships (both OM and MI)
            # Default to MI:0000 (molecular interaction) if no parent is specified
            is_a = parent_accession if parent_accession else "MI:0000"
            
            terms[accession] = {
                "accession": accession,
                "name": format_name(member.name),
                "definition": definition,
                "is_a": is_a,
            }
    
    # Add parent terms that aren't already included as regular terms
    for parent_acc, parent_data in parent_terms.items():
        if parent_acc not in terms:
            terms[parent_acc] = parent_data
    
    return list(terms.values())


def format_obo(terms: list[dict]) -> str:
    """Format OM terms as OBO content to append to PSI-MI.
    
    Args:
        terms: List of term dicts with accession, name, definition, is_a
        
    Returns:
        OBO format string (just the [Term] entries, no header)
    """
    lines = []
    
    # Sort terms by accession for consistent output
    for term in sorted(terms, key=lambda t: t["accession"]):
        lines.append("[Term]")
        lines.append(f"id: {term['accession']}")
        lines.append(f"name: {term['name']}")
        
        if term["definition"]:
            escaped_def = escape_obo_string(term["definition"])
            lines.append(f'def: "{escaped_def}" []')
        
        if term["is_a"]:
            lines.append(f"is_a: {term['is_a']}")
        
        lines.append("")
    
    return "\n".join(lines)


def main(output_path: str = "omnipath_mi.obo") -> None:
    """Generate combined PSI-MI + OmniPath OBO file.
    
    Args:
        output_path: Path for the output OBO file
    """
    # Fetch and fix PSI-MI
    psi_mi = fetch_psi_mi()
    psi_mi = fix_malformed_dates(psi_mi)
    print("  Fixed malformed dates")
    
    # Extract OM terms
    print("Extracting OM terms from cv_terms...")
    terms = extract_om_terms()
    print(f"Found {len(terms)} OM terms")
    
    # Format OM terms
    print("Generating OBO format...")
    om_obo = format_obo(terms)
    
    # Combine: PSI-MI + OmniPath terms
    combined = psi_mi + "\n" + om_obo
    
    output = Path(output_path)
    output.write_text(combined)
    print(f"Wrote {output.resolve()}")
    
    # Print summary
    root_om_terms = [t for t in terms if t["is_a"] is None]
    mi_child_terms = [t for t in terms if t["is_a"] and t["is_a"].startswith("MI:")]
    om_child_terms = [t for t in terms if t["is_a"] and t["is_a"].startswith("OM:")]
    
    print(f"\nSummary:")
    print(f"  OM terms: {len(terms)}")
    print(f"  OM root terms (no parent): {len(root_om_terms)}")
    print(f"  OM terms extending MI: {len(mi_child_terms)}")
    print(f"  OM terms under OM parents: {len(om_child_terms)}")
    
    if root_om_terms:
        print(f"\n  OmniPath root terms:")
        for rt in sorted(root_om_terms, key=lambda t: t["accession"])[:10]:
            print(f"    - {rt['accession']}: {rt['name']}")
        if len(root_om_terms) > 10:
            print(f"    ... and {len(root_om_terms) - 10} more")


if __name__ == "__main__":
    output = sys.argv[1] if len(sys.argv) > 1 else "omnipath_mi.obo"
    main(output)
