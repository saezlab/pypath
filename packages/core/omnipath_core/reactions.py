"""Vocabulary of reaction records: their type, their compartments, their catalysts.

A reaction is a Biolink ``molecular_activity``: participants hang off it through
``has_input`` / ``has_output``, catalysts through ``enabled_by`` (or the inverse
``catalyzes`` a resource states separately), and each participant link says
where the participant is. Input modules emit reactions in this shape and nothing
else; binary statements such as transport are projected from it by the build
(:mod:`omnipath_build.projections`).

The type of a reaction is an SBO term in the Biolink ``type`` slot (rdf:type):
Biolink has no subclasses of ``molecular_activity``. Recon3D publishes SBO on
every reaction, and the BioPAX conversion classes map onto it one to one.
"""

from __future__ import annotations

# The slot a reaction's SBO terms are stated under, as an entity annotation.
REACTION_TYPE = "type"
# Where a participant is (on the participant link) and where a reaction takes
# place (on the reaction entity): GO cellular-component CURIEs. A Biolink
# qualifier, so on a participant link it is part of the statement's identity.
COMPARTMENT = "anatomical_context_qualifier"

BIOCHEMICAL_OR_TRANSPORT_REACTION = "SBO:0000167"
BIOCHEMICAL_REACTION = "SBO:0000176"
# "translocation reaction". SBO also has "transport reaction" (SBO:0000655) as
# its sibling; Recon3D, like every COBRA model, types transport with 0000185,
# so that is the term the other resources are mapped to.
TRANSPORT = "SBO:0000185"
NON_COVALENT_BINDING = "SBO:0000177"
DEGRADATION = "SBO:0000179"
DISSOCIATION = "SBO:0000180"
TRANSCRIPTION = "SBO:0000183"
TRANSLATION = "SBO:0000184"
EXCHANGE = "SBO:0000627"
DEMAND = "SBO:0000628"
BIOMASS_PRODUCTION = "SBO:0000629"
PSEUDOREACTION = "SBO:0000631"
SINK = "SBO:0000632"

TRANSPORT_TYPES = frozenset({TRANSPORT, "SBO:0000587", "SBO:0000588", "SBO:0000655"})
# Model devices rather than biology: they open the model's boundary or drain
# its biomass. Kept, typed, so a consumer can leave them out.
PSEUDOREACTION_TYPES = frozenset({PSEUDOREACTION, EXCHANGE, DEMAND, SINK, BIOMASS_PRODUCTION})

SBO_LABELS = {
    BIOCHEMICAL_OR_TRANSPORT_REACTION: "biochemical or transport reaction",
    BIOCHEMICAL_REACTION: "biochemical reaction",
    TRANSPORT: "translocation reaction",
    "SBO:0000655": "transport reaction",
    NON_COVALENT_BINDING: "non-covalent binding",
    DEGRADATION: "degradation",
    DISSOCIATION: "dissociation",
    TRANSCRIPTION: "transcription",
    TRANSLATION: "translation",
    EXCHANGE: "exchange reaction",
    DEMAND: "demand reaction",
    BIOMASS_PRODUCTION: "biomass production",
    PSEUDOREACTION: "pseudoreaction",
    SINK: "sink reaction",
}

# The compartments of the genome-scale models (Recon3D, Human-GEM), by the GO
# name their parsers already translate the compartment codes to.
# "extracellular space" (GO:0005615) is obsolete, replaced by
# "extracellular region".
GO_COMPARTMENTS = {
    "cytosol": "GO:0005829",
    "extracellular space": "GO:0005576",
    "extracellular region": "GO:0005576",
    "mitochondrion": "GO:0005739",
    "mitochondrial intermembrane space": "GO:0005758",
    "endoplasmic reticulum": "GO:0005783",
    "peroxisome": "GO:0005777",
    "Golgi apparatus": "GO:0005794",
    "lysosome": "GO:0005764",
    "nucleus": "GO:0005634",
}


def go_compartment(name: str | None) -> str | None:
    """The GO cellular component of a model compartment name, or None."""
    return GO_COMPARTMENTS.get((name or "").strip())


def compartment_context(compartments) -> str:
    """The identity a reaction's compartments add to its entity key.

    The same chemistry in different compartments, or the same cargo across
    different membranes, is a different reaction. The set is unordered: a
    resource's uptake and efflux of one cargo across one membrane are one
    reaction, as Rhea's direction-free master reaction is.
    """
    return "|".join(sorted({c.strip() for c in compartments or () if c and c.strip()}))
