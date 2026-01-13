#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Export OmniPath CV terms to OBO format.

This script extracts all OmniPath-specific controlled vocabulary terms (OM:* accessions)
from pypath.internals.cv_terms and generates an OBO format ontology file.

Usage:
    python scripts/export_omnipath_obo.py [output_path]
    
    output_path: Optional path for the OBO file (default: omnipath.obo)
"""

from __future__ import annotations

import inspect
import re
import sys
from datetime import datetime
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from pypath.internals.cv_terms import CvEnum
from pypath.internals import cv_terms


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
        List of dicts with keys: accession, name, definition, parent, xref
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
                        "parent": None,  # Root terms have no parent
                        "xref": None,
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
            
            # Determine is_a (OM:* only) and xref (MI:* references)
            is_a = None
            xref = None
            if parent_accession:
                if parent_accession.startswith("OM:"):
                    is_a = parent_accession
                elif parent_accession.startswith("MI:"):
                    xref = parent_accession
            
            terms[accession] = {
                "accession": accession,
                "name": format_name(member.name),
                "definition": definition,
                "parent": is_a,
                "xref": xref,
            }
    
    # Add parent terms that aren't already included as regular terms
    for parent_acc, parent_data in parent_terms.items():
        if parent_acc not in terms:
            terms[parent_acc] = parent_data
    
    return list(terms.values())


def format_obo(terms: list[dict]) -> str:
    """Format terms as OBO file content.
    
    Args:
        terms: List of term dicts with accession, name, definition, parent
        
    Returns:
        OBO format string
    """
    lines = [
        "format-version: 1.2",
        f"date: {datetime.now().strftime('%d:%m:%Y %H:%M')}",
        "ontology: omnipath",
        "default-namespace: OMNIPATH",
        "",
    ]
    
    # Sort terms by accession for consistent output
    for term in sorted(terms, key=lambda t: t["accession"]):
        lines.append("[Term]")
        lines.append(f"id: {term['accession']}")
        lines.append(f"name: {term['name']}")
        
        if term["definition"]:
            escaped_def = escape_obo_string(term["definition"])
            lines.append(f'def: "{escaped_def}" []')
        
        if term["parent"]:
            lines.append(f"is_a: {term['parent']}")
        
        if term.get("xref"):
            lines.append(f"xref: {term['xref']}")
        
        lines.append("")
    
    return "\n".join(lines)


def main(output_path: str = "omnipath.obo") -> None:
    """Generate OmniPath OBO file from cv_terms.
    
    Args:
        output_path: Path for the output OBO file
    """
    print("Extracting OM terms from cv_terms...")
    terms = extract_om_terms()
    print(f"Found {len(terms)} OM terms")
    
    print("Generating OBO format...")
    obo_content = format_obo(terms)
    
    output = Path(output_path)
    output.write_text(obo_content)
    print(f"Wrote {output.resolve()}")
    
    # Print summary
    root_terms = [t for t in terms if t["parent"] is None]
    print(f"\nSummary:")
    print(f"  Total terms: {len(terms)}")
    print(f"  Root terms: {len(root_terms)}")
    for rt in sorted(root_terms, key=lambda t: t["accession"]):
        print(f"    - {rt['accession']}: {rt['name']}")


if __name__ == "__main__":
    output = sys.argv[1] if len(sys.argv) > 1 else "omnipath.obo"
    main(output)
