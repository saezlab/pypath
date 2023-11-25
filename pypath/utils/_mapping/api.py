#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
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

from pypath.utils._mapping.app import get_mapper


def map_name(
        name,
        id_type,
        target_id_type,
        ncbi_tax_id = None,
        strict = False,
        expand_complexes = True,
        uniprot_cleanup = True,
    ):
    """
    Translates one instance of one ID type to a different one.
    Returns set of the target ID type.

    This function should be used to convert individual IDs.
    It takes care about everything and ideally you don't need to
    think on the details.

    How does it work: looks up dictionaries between the original
    and target ID type, if doesn't find, attempts to load from the
    predefined inputs.
    If the original name is genesymbol, first it looks up among the
    preferred gene names from UniProt, if not found, it takes an
    attempt with the alternative gene names. If the gene symbol
    still couldn't be found, and strict = False, the last attempt
    only the first 5 characters of the gene symbol matched. If the
    target name type is uniprot, then it converts all the ACs to
    primary. Then, for the Trembl IDs it looks up the preferred gene
    names, and find Swissprot IDs with the same preferred gene name.

    Args
        name (str): The original name to be converted.
        id_type (str): The type of the name. Available by default:
            - genesymbol (gene name)
            - entrez (Entrez Gene ID \[#\])
            - refseqp (NCBI RefSeq Protein ID \[NP\_\*|XP\_\*\])
            - ensp (Ensembl protein ID \[ENSP\*\])
            - enst (Ensembl transcript ID \[ENST\*\])
            - ensg (Ensembl genomic DNA ID \[ENSG\*\])
            - hgnc (HGNC ID \[HGNC:#\])
            - gi (GI number \[#\])
            - embl (DDBJ/EMBL/GeneBank CDS accession)
            - embl_id (DDBJ/EMBL/GeneBank accession)
            And many more, see the code of
            ``pypath.internals.input_formats``
        target_id_type (str): The name type to translate to, more or
            less the same values are available as for ``id_type``.
        ncbi_tax_id (int): NCBI Taxonomy ID of the organism.
        strict (bool): In case a Gene Symbol can not be translated,
            try to add number "1" to the end, or try to match only
            its first five characters. This option is rarely used,
            but it makes possible to translate some non-standard
            gene names typically found in old, unmaintained resources.
        expand_complexes (bool): When encountering complexes,
            translated the IDs of its components and return a set
            of IDs. The alternative behaviour is to return the
            `Complex` objects.
        uniprot_cleanup (bool): When the `target_id_type` is UniProt
            ID, call the `uniprot_cleanup` function at the end.
    """

    mapper = get_mapper()

    return mapper.map_name(
        name = name,
        id_type = id_type,
        target_id_type = target_id_type,
        ncbi_tax_id = ncbi_tax_id,
        strict = strict,
        expand_complexes = expand_complexes,
        uniprot_cleanup = uniprot_cleanup,
    )


def map_name0(
        name,
        id_type,
        target_id_type,
        ncbi_tax_id = None,
        strict = False,
        expand_complexes = True,
        uniprot_cleanup = True,
    ):
    """
    Translates the name and returns only one of the resulted IDs. It
    means in case of ambiguous ID translation, a random one of them
    will be picked and returned. Recommended to use only if the
    translation between the given ID types is mostly unambigous and
    the loss of information can be ignored. See more details at
    `map_name`.
    """

    mapper = get_mapper()

    return mapper.map_name0(
        name = name,
        id_type = id_type,
        target_id_type = target_id_type,
        ncbi_tax_id = ncbi_tax_id,
        strict = strict,
        expand_complexes = expand_complexes,
        uniprot_cleanup = uniprot_cleanup,
    )


def map_names(
        names,
        id_type = None,
        target_id_type = None,
        ncbi_tax_id = None,
        strict = False,
        expand_complexes = True,
        uniprot_cleanup = True,
    ):
    """
    Same as ``map_name`` but translates multiple IDs at once. These two
    functions could be seamlessly implemented as one, still I created
    separate functions to always make it explicit if a set of translated
    IDs come from multiple original IDs.

    Args
        name (str): The original name to be converted.
        id_type (str): The type of the name. Available by default:
            - genesymbol (gene name)
            - entrez (Entrez Gene ID \[#\])
            - refseqp (NCBI RefSeq Protein ID \[NP\_\*|XP\_\*\])
            - ensp (Ensembl protein ID \[ENSP\*\])
            - enst (Ensembl transcript ID \[ENST\*\])
            - ensg (Ensembl genomic DNA ID \[ENSG\*\])
            - hgnc (HGNC ID \[HGNC:#\])
            - gi (GI number \[#\])
            - embl (DDBJ/EMBL/GeneBank CDS accession)
            - embl_id (DDBJ/EMBL/GeneBank accession)
            And many more, see the code of
            ``pypath.internals.input_formats``
        target_id_type (str): The name type to translate to, more or
            less the same values are available as for ``id_type``.
        ncbi_tax_id (int): NCBI Taxonomy ID of the organism.
        strict (bool): In case a Gene Symbol can not be translated,
            try to add number "1" to the end, or try to match only
            its first five characters. This option is rarely used,
            but it makes possible to translate some non-standard
            gene names typically found in old, unmaintained resources.
        expand_complexes (bool): When encountering complexes,
            translated the IDs of its components and return a set
            of IDs. The alternative behaviour is to return the
            `Complex` objects.
        uniprot_cleanup (bool): When the `target_id_type` is UniProt
            ID, call the `Mapper.uniprot_cleanup` function at the end.
    """

    mapper = get_mapper()

    return mapper.map_names(
        names = names,
        id_type = id_type,
        target_id_type = target_id_type,
        ncbi_tax_id = ncbi_tax_id,
        strict = strict,
        expand_complexes = expand_complexes,
        uniprot_cleanup = uniprot_cleanup,
    )


def label(name, id_type = None, entity_type = None, ncbi_tax_id = 9606):
    """
    For any kind of entity, either protein, miRNA or protein complex,
    returns the preferred human readable label. For proteins this means
    Gene Symbols, for miRNAs miRNA names, for complexes a series of
    Gene Symbols.
    """

    mapper = get_mapper()

    return mapper.label(
        name = name,
        id_type = id_type,
        entity_type = entity_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def guess_type(name, entity_type = None):
    """
    From a string, tries to guess the ID type and optionally the entity
    type. Returns a tuple of strings: ID type and entity type.
    """

    mapper = get_mapper()

    return mapper.guess_type(name = name, entity_type = entity_type)


def id_from_label(label, label_id_type = 'genesymbol', ncbi_tax_id = None):
    """
    For a label (e.g. Gene Symbol) returns the corresponding IDs (e.g.
    UniProt IDs).
    """

    mapper = get_mapper()

    return mapper.id_from_label(
        label = label,
        label_id_type = label_id_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def id_from_label0(label, label_id_type = 'genesymbol', ncbi_tax_id = None):
    """
    For a label (e.g. Gene Symbol) returns a single ID (e.g. UniProt IDs).
    """

    mapper = get_mapper()

    return mapper.id_from_label0(
        label = label,
        label_id_type = label_id_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def translation_dict(
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int | None = None,
    ) -> MappingTable | None:
    """
    Identifier translation table as a dict of sets.
    """

    mapper = get_mapper()

    return mapper.translation_dict(
        id_type = id_type,
        target_id_type = target_id_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def translation_df(
        id_type: str,
        target_id_type: str,
        ncbi_tax_id: int | None = None,
    ) -> MappingTable | None:
    """
    Identifier translation table as a `pandas.DataFrame`.
    """

    mapper = get_mapper()

    return mapper.translation_df(
        id_type = id_type,
        target_id_type = target_id_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def mapping_tables() -> list[MappingTableDefinition]:
    """
    A list of built-in mapping tables.

    If `id_type_b` is `None`, that means translation to all other ID types
    provided by the same resource is possible.
    """

    return get_mapper().mapping_tables()


def id_types() -> list[IdType]:
    """
    Identifier types with their labels.
    """

    return get_mapper().id_types()
