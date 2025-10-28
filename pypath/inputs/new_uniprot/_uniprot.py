#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

from __future__ import annotations

import csv
from collections.abc import Generator

from pypath.share.downloads import download_and_open
from ._records import UniProtRecord

__all__ = [
    'uniprot_data',
]

# UniProt REST API URL for comprehensive protein data
UNIPROT_DATA_URL = (
    'https://rest.uniprot.org/uniprotkb/stream?'
    'compressed=true&'
    'fields=accession%2Cid%2Cprotein_name%2Clength%2Cmass%2Csequence%2C'
    'gene_primary%2Cgene_synonym%2Corganism_id%2Ccc_disease%2Cft_mutagen%2C'
    'cc_subcellular_location%2Ccc_ptm%2Clit_pubmed_id%2Ccc_function%2C'
    'xref_ensembl%2Cxref_kegg%2Ccc_pathway%2Ccc_activity_regulation%2C'
    'keyword%2Cec%2Cgo%2Cft_transmem%2Cprotein_families%2Cxref_refseq%2C'
    'xref_alphafolddb%2Cxref_pdb%2Cxref_chembl%2Cxref_phosphositeplus%2C'
    'xref_signor%2Cxref_pathwaycommons%2Cxref_intact%2Cxref_biogrid%2C'
    'xref_complexportal&'
    'format=tsv&'
    'query=%28%28taxonomy_id%3A9606%29+OR+%28taxonomy_id%3A10090%29+'
    'OR+%28taxonomy_id%3A10116%29%29+AND+%28reviewed%3Atrue%29'
)


def uniprot_data(organism: int | list[int] | None = None) -> Generator[UniProtRecord]:
    """
    Download UniProt data in TSV format with comprehensive fields.

    Args:
        organism: NCBI taxonomy ID(s). Can be a single int, list of ints, or None.
                  If None, defaults to [9606, 10090, 10116] (human, mouse, rat).
                  Common values: 9606 for human, 10090 for mouse, 10116 for rat.
                  Note: Currently the URL is hardcoded for human, mouse, and rat.

    Yields:
        UniProtRecord named tuples with protein information
    """
    # Handle organism parameter
    if organism is None:
        organisms = [9606, 10090, 10116]  # Default: human, mouse, rat
    elif isinstance(organism, int):
        organisms = [organism]
    else:
        organisms = list(organism)

    url = UNIPROT_DATA_URL

    # Use download_and_open which handles gzip decompression automatically
    opener = download_and_open(
        url,
        filename=f'uniprot_data_{"_".join(map(str, organisms))}.tsv.gz',
        subfolder='uniprot',
        large=True,
        encoding='utf-8',
        default_mode='r',
        ext='gz',
    )

    if opener and opener.result:
        # opener.result is a file handle for gz files
        reader = csv.DictReader(opener.result, delimiter='\t')

        for row in reader:
            # Helper function to safely get and strip values
            def get_field(field_name: str, converter=None):
                value = row.get(field_name, '').strip()
                if not value or value == '':
                    return None
                if converter:
                    try:
                        return converter(value)
                    except (ValueError, TypeError):
                        return None
                return value

            yield UniProtRecord(
                accession=get_field('Entry'),
                entry_name=get_field('Entry Name'),
                protein_name=get_field('Protein names'),
                length=get_field('Length', int),
                mass=get_field('Mass', int),
                sequence=get_field('Sequence'),
                gene_primary=get_field('Gene Names (primary)'),
                gene_synonym=get_field('Gene Names (synonym)'),
                organism_id=get_field('Organism (ID)', int),
                cc_disease=get_field('Involvement in disease'),
                ft_mutagen=get_field('Mutagenesis'),
                cc_subcellular_location=get_field('Subcellular location [CC]'),
                cc_ptm=get_field('Post-translational modification'),
                lit_pubmed_id=get_field('PubMed ID'),
                cc_function=get_field('Function [CC]'),
                xref_ensembl=get_field('Ensembl'),
                xref_kegg=get_field('KEGG'),
                cc_pathway=get_field('Pathway'),
                cc_activity_regulation=get_field('Activity regulation'),
                keyword=get_field('Keywords'),
                ec=get_field('EC number'),
                go=get_field('Gene Ontology IDs'),
                ft_transmem=get_field('Transmembrane'),
                protein_families=get_field('Protein families'),
                xref_refseq=get_field('RefSeq'),
                xref_alphafolddb=get_field('AlphaFoldDB'),
                xref_pdb=get_field('PDB'),
                xref_chembl=get_field('ChEMBL'),
                xref_phosphositeplus=get_field('PhosphoSitePlus'),
                xref_signor=get_field('SIGNOR'),
                xref_pathwaycommons=get_field('PathwayCommons'),
                xref_intact=get_field('IntAct'),
                xref_biogrid=get_field('BioGRID'),
                xref_complexportal=get_field('ComplexPortal'),
            )
