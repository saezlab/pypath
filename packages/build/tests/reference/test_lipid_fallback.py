from omnipath_resolver import observations as replay
from omnipath_resolver.canonical.match import Vote


def test_lipid_fallback_is_only_for_identifier_poor_chemicals(monkeypatch):
    def result(name):
        return dict(
            status="parsed", goslin=name, specificity={"species": 4, "sn_position": 8}[name]
        )

    monkeypatch.setattr(replay, "lipid_name_result", result)
    votes = replay.lipid_fallback_votes(
        [], {"name": ["species"], "synonym": ["sn_position"]}, "chemical"
    )
    assert [(v.ns, v.id) for v in votes] == [("goslin", "sn_position")]
    marker = Vote("inchikey", "AAAAAAAAAAAAAA-BBBBBBBBSA-N", None, "id", False)
    assert replay.lipid_fallback_votes([marker], {"name": ["must not parse"]}, "chemical") == [
        marker
    ]
    assert replay.lipid_fallback_votes([], {"name": ["must not parse"]}, "gene_protein") == []
    native = Vote("refmet", "RM1", None, "id", False)
    assert [
        (v.ns, v.id)
        for v in replay.lipid_fallback_votes([native], {"name": ["species"]}, "chemical")
    ] == [("refmet", "RM1"), ("goslin", "species")]
