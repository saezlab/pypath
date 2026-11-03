"""Ordinary source attributes retained without changing Biolink statement identity.

BioPAX supplies reaction context terms; source-specific attributes retain their
own names where Biolink has no equivalent. These are not edge qualifiers.
"""

CONVERSION_DIRECTION = "biopax:conversionDirection"
CELLULAR_LOCATION = "biopax:cellularLocation"
# Side of the membrane a transported participant is on in a Rhea equation:
# "in" or "out". A side, not a compartment.
TRANSPORT_SIDE = "rhea:transport_side"
TRAIT_TYPE = "macdb:trait_type"
PARTICIPANT_ROLE = "connectomedb:participant_role"
# The source's own name for a reaction's type (a BioPAX class, a Human-GEM
# subsystem), kept beside the SBO term it is mapped to.
SOURCE_REACTION_TYPE = "omnipath:source_reaction_type"
# The source's identifier of one reaction. A resource that states a reaction's
# catalysts apart from its participants (Reactome's controls) carries it on
# both, so the two pair up even where several of its reactions resolve to one
# entity.
SOURCE_REACTION = "omnipath:source_reaction"
# Where a transported cargo comes from and goes to: GO cellular components, or
# Rhea's in/out where the source names no compartment.
TRANSPORT_FROM = "omnipath:transport_from"
TRANSPORT_TO = "omnipath:transport_to"
# The reaction entity a projected relation was read from.
PROJECTED_FROM = "omnipath:projected_from"
