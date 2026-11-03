"""Exact Goslin lipid-name normalization shared by construction and runtime."""

from __future__ import annotations

_PARSER = None


def normalize(name):
    global _PARSER
    from pygoslin.parser.Parser import LipidParser
    from pygoslin.domain.LipidLevel import LipidLevel
    from pygoslin.domain.LipidExceptions import LipidException

    if _PARSER is None:
        _PARSER = LipidParser()
    result = dict(name=name, goslin=None, level=None, status="unparsed", specificity=0)
    try:
        lipid = _PARSER.parse(name)
    except LipidException:
        return result
    if lipid is None:
        return result
    level = lipid.lipid.info.level
    if level.value < LipidLevel.SPECIES.value:
        result["status"] = "class_only"
        return result
    if lipid.adduct is not None:
        result["status"] = "adduct_or_isotope"
        return result
    lipid.sort_fatty_acyl_chains()  # only unordered chains; preserves known sn positions
    label = lipid.get_lipid_string()
    result.update(
        goslin=level.name.lower() + ":" + label,
        level=level.name.lower(),
        specificity=level.value,
        status="parsed",
    )
    return result
