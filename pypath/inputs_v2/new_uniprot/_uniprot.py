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
UNIPROT_DATA_URL = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession,id,protein_name,length,mass,sequence,gene_primary,gene_synonym,organism_id,cc_disease,ft_mutagen,cc_subcellular_location,cc_ptm,lit_pubmed_id,cc_function,xref_ensembl,xref_kegg,cc_pathway,cc_activity_regulation,keywordid,ec,go_id,ft_transmem,protein_families,xref_refseq,xref_alphafolddb,xref_pdb,xref_chembl,xref_phosphositeplus,xref_signor,xref_pathwaycommons,xref_intact,xref_biogrid,xref_complexportal&format=tsv&query=(taxonomy_id:9606 OR taxonomy_id:10090 OR taxonomy_id:10116) AND reviewed:true"


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

            # Helper function to split semicolon-separated values
            def get_list_field(field_name: str):
                value = row.get(field_name, '').strip()
                if not value or value == '':
                    return None
                return [item.strip() for item in value.split(';') if item.strip()]

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
                ec=get_field('EC number'),
                keyword=get_list_field('Keywords IDs'),
                go=get_list_field('Gene Ontology IDs'),
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
