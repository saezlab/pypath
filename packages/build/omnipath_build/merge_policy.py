"""Associative serving projections; attributed assertions remain separate."""


def preferred_label(identifier: str, *labels: str | None) -> str:
    values = {str(label) for label in labels if label}
    if not values:
        return identifier
    return min(
        values,
        key=lambda value: (value == identifier or (value.isascii() and value.isdigit()), value),
    )
