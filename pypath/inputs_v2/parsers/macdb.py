"""Parsers and row-enrichment helpers for MACdb inputs_v2."""

from __future__ import annotations

from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.base import Download


def publication_pmc_map(
    publication_download: Download,
    *,
    force_refresh: bool = False,
    **kwargs: Any,
) -> dict[str, str]:
    """Map PubMed IDs to PubMed Central IDs from the publication table."""
    result = {}
    opener = publication_download.open(force_refresh=force_refresh, **kwargs)

    for row in iter_tsv(opener):
        pmid = row.get('PMID')
        pmcid = row.get('PMC_ID')

        if pmid and pmcid:
            result[pmid] = pmcid

    return result


def study_map(
    study_download: Download,
    publication_download: Download,
    *,
    force_refresh: bool = False,
    **kwargs: Any,
) -> dict[str, dict[str, str]]:
    """Map MACdb cohort IDs to study metadata and publication references."""
    pubmed_to_pmc = publication_pmc_map(
        publication_download,
        force_refresh=force_refresh,
        **kwargs,
    )
    result = {}
    opener = study_download.open(force_refresh=force_refresh, **kwargs)

    for row in iter_tsv(opener):
        cohort_id = row.get('Cohort_id')
        pmid = row.get('pubmed_id') or row.get('pubmed') or ''

        if cohort_id:
            row['pmc_id'] = pubmed_to_pmc.get(pmid, '')
            result[cohort_id] = row

    return result


def trait_name_map(
    trait_download: Download,
    *,
    force_refresh: bool = False,
    **kwargs: Any,
) -> dict[str, str]:
    """Map MACdb trait ontology IDs to trait labels."""
    result = {}
    opener = trait_download.open(force_refresh=force_refresh, **kwargs)

    for row in iter_tsv(opener):
        trait_id = row.get('Trait_Ontology_ID')
        trait_name = row.get('Trait_Ontology')

        if trait_id and trait_name:
            result[trait_id] = trait_name

    return result


def iter_associations(
    opener,
    *,
    study_download: Download,
    publication_download: Download,
    trait_download: Download,
    force_refresh: bool = False,
    **kwargs: Any,
):
    """Yield metabolite-trait association rows enriched with study evidence."""
    cohort_to_study = study_map(
        study_download,
        publication_download,
        force_refresh=force_refresh,
        **kwargs,
    )
    trait_names = trait_name_map(
        trait_download,
        force_refresh=force_refresh,
        **kwargs,
    )

    for row in iter_tsv(opener):
        study = cohort_to_study.get(row.get('Cohort_id'), {})
        row.update({f'study_{key}': value for key, value in study.items()})
        trait_id = study.get('Trait_onto_ID', '')
        row['trait_id'] = trait_id
        row['trait_name'] = trait_names.get(trait_id, '')
        row['pubmed_id'] = study.get('pubmed_id') or study.get('pubmed') or ''
        row['pmc_id'] = study.get('pmc_id', '')
        yield row
