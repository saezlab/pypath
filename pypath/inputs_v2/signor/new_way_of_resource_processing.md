New: 
Entity(
    source="my_database",
    entity_type=EntityType(EntityTypeCv.PROTEIN_COMPLEX),
    identifiers=Identifiers(
        Column(0, 
               cv=IdentifierNamespaceCv.NAME),
        
        Column(1,
               delimiter=";",
               processing={
                   "extract_prefix": r"^([a-z]+):",   # capture prefix before colon
                   "extract_value": r"^[a-z]+:(.*)",  # capture value after prefix
               },
               cv={
                   "chebi": IdentifierNamespaceCv.CHEBI,
                   "go": IdentifierNamespaceCv.GO,
                   "ec": IdentifierNamespaceCv.ENZYME,
                   "default": IdentifierNamespaceCv.SYNONYM,
               }),
    ),
    annotations=Annotations(
        Column(3,
               delimiter=";",
               processing={
                   "extract_value": r"([\d.]+)",
                   "extract_unit": r"[ \t]*([a-zA-Zμ]+)$",
               },
               term_cv=AnnotationTypeCv.MOLECULAR_WEIGHT,
               unit_cv=IdentifierNamespaceCv.UNIT),
    ),
    members=Members(
        Entities(  # Plural - creates multiple Entity objects, one per split value
            entity_type=EntityTypeCv.PROTEIN,
            entity_source="member",
            identifiers=Identifiers(
                Column(2,                # "P12345;Q67890;P11111"
                       delimiter=";",    # Splits into 3 member entities
                       cv=IdentifierNamespaceCv.UNIPROT)
            ),
            annotations=Annotations(
                Column(4,                # "2;1;4" (aligned with column 2 splits)
                       delimiter=";",
                       processing={
                           "extract_value": r"(\d+)"
                       },
                       term_cv=AnnotationTypeCv.STOICHIOMETRY),
                Column(5,                # "catalytic;regulatory;structural"
                       delimiter=";",
                       term_cv=AnnotationTypeCv.ROLE)
            )
        )
    )
)

Old:

# If unit is in separate column, first join it to the value column
COLUMN_SCHEMA = {
    0: (Identifier, {
        "cv": IdentifierNamespaceCv.NAME,
    }),
    1: (Identifier, {
        "delimiter": ";",
        "processing": {
            "extract_prefix": r"^([a-z]+):",   # capture prefix before colon
            "extract_value":  r"^[a-z]+:(.*)", # capture value after prefix
        },
        "cv": {
            "chebi": IdentifierNamespaceCv.CHEBI,
            "go": IdentifierNamespaceCv.GO,
            "ec": IdentifierNamespaceCv.ENZYME,
            "default": IdentifierNamespaceCv.SYNONYM,
        },
    }),
    2: (Membership.member.identifiers[0], {
        "delimiter": ";",
        "cv": IdentifierNamespaceCv.UNIPROT,
    }),
    3: (Annotation, {
        "delimiter": ";",
        "processing": {
            "extract_value": r"([\d.]+)",
            "extract_unit":  r"[ \t]*([a-zA-Zμ]+)$",
        },
        "cv": {
            'value': AnnotationTypeCv.MOLECULAR_WEIGHT,
            'unit': IdentifierNamespaceCv.UNIT,
        },
    }),
}
