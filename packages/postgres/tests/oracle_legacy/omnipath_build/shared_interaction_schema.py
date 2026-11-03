"""Shared CV term sets for interpreting interaction direction and roles."""

from __future__ import annotations

from pypath.internals.cv_terms import (
    BiologicalEffectCv,
    BiologicalRoleCv,
    CausalStatementCv,
    cv_term_label_accession,
    InteractionParameterCv,
    LigandTypeCv,
    ParticipantMetadataCv,
    PharmacologicalActionCv,
)


def _term_forms(term: object) -> tuple[str, ...]:
    label_accession = cv_term_label_accession(term)
    raw = str(term)
    return (raw,) if label_accession == raw else (raw, label_accession)


def _optional_accession(enum_cls: type, name: str) -> tuple[str, ...]:
    member = getattr(enum_cls, name, None)
    return () if member is None else _term_forms(member)


POSITIVE_SIGN_ACCESSIONS = frozenset(
    form
    for term in (
        CausalStatementCv.UP_REGULATES,
        CausalStatementCv.UP_REGULATES_ACTIVITY,
        CausalStatementCv.UP_REGULATES_QUANTITY,
        CausalStatementCv.UP_REGULATES_QUANTITY_BY_EXPRESSION,
        CausalStatementCv.UP_REGULATES_QUANTITY_BY_STABILIZATION,
        BiologicalEffectCv.UP_REGULATES_ACTIVITY,
        BiologicalEffectCv.UP_REGULATES_QUANTITY,
        BiologicalRoleCv.STIMULATOR,
        PharmacologicalActionCv.AGONIST,
        PharmacologicalActionCv.FULL_AGONIST,
        PharmacologicalActionCv.PARTIAL_AGONIST,
        PharmacologicalActionCv.BIASED_AGONIST,
        PharmacologicalActionCv.ACTIVATION,
        PharmacologicalActionCv.POSITIVE,
        PharmacologicalActionCv.POTENTIATION,
        LigandTypeCv.ACTIVATOR,
        LigandTypeCv.AGONIST,
    )
    for form in _term_forms(term)
)

NEGATIVE_SIGN_ACCESSIONS = frozenset(
    form
    for term in (
        CausalStatementCv.DOWN_REGULATES,
        CausalStatementCv.DOWN_REGULATES_ACTIVITY,
        CausalStatementCv.DOWN_REGULATES_QUANTITY,
        CausalStatementCv.DOWN_REGULATES_QUANTITY_BY_DESTABLIZATION,
        CausalStatementCv.DOWN_REGULATES_QUANTITY_BY_REPRESSION,
        BiologicalEffectCv.DOWN_REGULATES_ACTIVITY,
        BiologicalEffectCv.DOWN_REGULATES_QUANTITY,
        BiologicalRoleCv.INHIBITOR,
        PharmacologicalActionCv.ANTAGONIST,
        PharmacologicalActionCv.INVERSE_AGONIST,
        PharmacologicalActionCv.INHIBITION,
        PharmacologicalActionCv.IRREVERSIBLE_INHIBITION,
        PharmacologicalActionCv.FEEDBACK_INHIBITION,
        PharmacologicalActionCv.NEGATIVE,
        PharmacologicalActionCv.PORE_BLOCKER,
        LigandTypeCv.INHIBITOR,
        LigandTypeCv.ANTAGONIST,
        LigandTypeCv.CHANNEL_BLOCKER,
        LigandTypeCv.GATING_INHIBITOR,
    )
    for form in _term_forms(term)
)

SOURCE_ROLE_ACCESSIONS = frozenset(
    [
        form
        for term in (
            ParticipantMetadataCv.SOURCE,
            BiologicalRoleCv.CONTROLLER,
            BiologicalRoleCv.REACTANT,
            BiologicalRoleCv.TEMPLATE,
            BiologicalRoleCv.ENZYME,
            BiologicalRoleCv.CATALYST,
            BiologicalRoleCv.REGULATOR,
            BiologicalRoleCv.INHIBITOR,
            BiologicalRoleCv.STIMULATOR,
            BiologicalRoleCv.ALLOSTERIC_EFFECTOR,
            LigandTypeCv.INHIBITOR,
            LigandTypeCv.ACTIVATOR,
            LigandTypeCv.AGONIST,
            LigandTypeCv.ANTAGONIST,
        )
        for form in _term_forms(term)
    ]
    + list(_optional_accession(ParticipantMetadataCv, 'LIGAND'))
)

TARGET_ROLE_ACCESSIONS = frozenset(
    [
        form
        for term in (
            ParticipantMetadataCv.TARGET,
            BiologicalRoleCv.CONTROLLED,
            BiologicalRoleCv.PRODUCT,
            BiologicalRoleCv.SUBSTRATE,
            BiologicalRoleCv.REGULATOR_TARGET,
        )
        for form in _term_forms(term)
    ]
    + list(_optional_accession(ParticipantMetadataCv, 'RECEPTOR'))
)

ACTIVATORY_PARAMETER_ACCESSIONS = frozenset({
    InteractionParameterCv.EC50.value,
})

INHIBITORY_PARAMETER_ACCESSIONS = frozenset({
    InteractionParameterCv.KI.value,
    InteractionParameterCv.IC50.value,
})
