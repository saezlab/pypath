"""Parser for PHI-base 5 pathogen-host interaction phenotype annotations. Filters for human and mouse hosts."""

from __future__ import annotations

import json

from collections.abc import Generator
from typing import Any


def iter_phi_phenotypes(
        opener,
        **_kwargs: Any,
    ) -> Generator[dict[str, Any], None, None]:
    """Yield flat pathogen protein-phenotype observations from PHI-base JSON."""

    #identify the json file in the zip
    members = getattr(opener, 'result', {})
    handle = None

    for name, member in members.items():

        if name.endswith('.json') and not name.endswith('.schema.json'):
            handle = member
            break

    if handle is None:
        raise ValueError('PHI-base: JSON data member not found in archive')

    sessions = json.load(handle).get('curation_sessions', {}) # session is a group of related rows from one paper

    for session_id, session in sessions.items():
        annotations = session.get('annotations', [])
        metagenotypes = session.get('metagenotypes', {})
        genotypes = session.get('genotypes', {})
        alleles = session.get('alleles', {})
        genes = session.get('genes', {})
        organisms = session.get('organisms', {})

        for annotation in annotations:
            if annotation.get('type') != 'pathogen_host_interaction_phenotype':
                continue

            metagenotype_id = annotation.get('metagenotype')
            metagenotype = metagenotypes.get(metagenotype_id, {})
            pathogen_genotype = genotypes.get(
                metagenotype.get('pathogen_genotype'),
                {},
            )
            host_genotype_id = metagenotype.get('host_genotype')
            host_genotype = genotypes.get(host_genotype_id, {})
            pathogen_taxon_id = pathogen_genotype.get('organism_taxonid')
            host_taxon_id = host_genotype.get('organism_taxonid')

            if host_taxon_id not in {9606, 10090}: # take only human and mouse
                continue

            for locus in pathogen_genotype.get('loci', []):
                for allele_ref in locus:
                    allele = alleles.get(allele_ref.get('id'), {})
                    gene = genes.get(allele.get('gene'), {})
                    uniprot_data = gene.get('uniprot_data', {})
                    uniprot_id = uniprot_data.get('uniprot_id')
                    term_id = annotation.get('term')

                    if not uniprot_id or not term_id or ':' not in term_id:
                        continue

                    record = {
                        'pathogen_uniprot': uniprot_id,
                        'pathogen_phig_id': gene.get('phig_id'),
                        'pathogen_gene_name': uniprot_data.get('name'),
                        'pathogen_product': uniprot_data.get('product'),
                        'pathogen_taxon_id': pathogen_taxon_id,
                        'pathogen_organism': organisms.get(
                            str(pathogen_taxon_id), {}
                        ).get('full_name'),
                        'pathogen_strain': pathogen_genotype.get(
                            'organism_strain'
                        ),
                        'host_taxon_id': host_taxon_id,
                        'host_organism': organisms.get(
                            str(host_taxon_id), {}
                        ).get('full_name'),
                        'host_strain': host_genotype.get('organism_strain'),
                        'host_genotype_id': host_genotype_id,
                        'metagenotype_id': metagenotype_id,
                        'phenotype_term_id': term_id,
                        'annotation_extensions': [
                            ' '.join(
                                str(value)
                                for value in (
                                    extension.get('relation'),
                                    extension.get('rangeDisplayName'),
                                    extension.get('rangeValue'),
                                )
                                if value
                            )
                            for extension in annotation.get('extension', [])
                        ],
                        'conditions': annotation.get('conditions', []),
                        'evidence_code': annotation.get('evidence_code'),
                        'figure': annotation.get('figure'),
                        'high_level_terms': annotation.get('high_level_terms', []),
                        'phi4_ids': annotation.get('phi4_id', []),
                        'publication': annotation.get('publication'),
                        'curation_session_id': session_id,
                        'curation_status': annotation.get('status'),
                        'checked': annotation.get('checked'),
                    }

                    for high_level_term in record.pop('high_level_terms'):
                        yield {**record, 'high_level_term': high_level_term}
