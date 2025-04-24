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

from future.utils import iteritems

import time
import datetime
import itertools

import timeloop

import pypath.share.common as common
import pypath.share.session as session
import pypath.share.settings as settings
import pypath_common._constants as _const
import pypath.inputs.uniprot as uniprot_input
import pypath.inputs.ensembl as ensembl_input

_logger = session.Logger(name = 'taxonomy')
_log = _logger._log

db = {}
_cleanup_period = settings.get('mapper_cleanup_interval')
_lifetime = 300
_last_used = {}
NOT_ORGANISM_SPECIFIC = _const.NOT_ORGANISM_SPECIFIC

# XXX: Shouldn't we keep all functions and variables separated
#      (together among them)?
taxids = {
    9606: 'human',
    10090: 'mouse',
    10116: 'rat',
    9031: 'chicken',
    9913: 'cow',
    9986: 'rabbit',
    9940: 'sheep',
    10141: 'guinea pig',
    10036: 'hamster',
    7227: 'fruit fly',
    9615: 'dog',
    9823: 'pig',
    8355: 'frog',
    9091: 'quail',
    9796: 'horse',
    9925: 'goat',
    89462: 'water buffalo',
    9598: 'monkey',
    9103: 'turkey',
    9685: 'cat',
    7604: 'starfish',
    7609: 'spiny starfish',
    1213717: 'torpedo',
    9669: 'ferret',
    8839: 'duck',
    9593: 'gorilla',
    7460: 'honeybee',
    8407: 'european common frog',
    9544: 'rhesus macaque',
}

taxids2 = dict(
    (
        t.taxon_id,
        t.common_name.lower()
    )
    for t in ensembl_input.ensembl_organisms()
)

taxa = common.swap_dict_simple(taxids)
taxa2 = common.swap_dict_simple(taxids2)


taxa_synonyms = {
    'bovine': 'cow',
    'western gorilla': 'gorilla',
}


phosphoelm_taxids = {
    9606: 'Homo sapiens',
    10090: 'Mus musculus',
    9913: 'Bos taurus',
    9986: 'Oryctolagus cuniculus',
    9615: 'Canis familiaris',
    10029: 'Cricetulus griseus',
    9267: 'Didelphis virginiana',
    9031: 'Gallus gallus',
    10036: 'Mesocricetus auratus',
    9940: 'Ovis aries',
    10116: 'Rattus norvegicus',
    9823: 'Sus scrofa',
    8355: 'Xenopus laevis',
    10141: 'Cavia porcellus',
    9796: 'Equus caballus',
    7227: 'Drosophila melanogaster',
    487: 'Neisseria meningitidis',
    562: 'Escherichia coli',
    5207: 'Cryptococcus neoformans',
    470: 'Acinetobacter baumannii',
    1280: 'Staphylococcus aureus',
    4932: 'Saccharomyces cerevisiae',
    34: 'Myxococcus xanthus',
    1392: 'Bacillus anthracis',
    210: 'Helicobacter pylori',
    6239: 'Caenorhabditis elegans',
}


phosphoelm_taxids.update(
    [
        (
            t.taxon_id,
            t.scientific_name
        )
        for t in ensembl_input.ensembl_organisms()
    ]
)


dbptm_taxids = {
    9606: 'HUMAN',
    10090: 'MOUSE',
    7227: 'DROME',
    10116: 'RAT',
    559292: 'YEAST',
    284812: 'SCHPO',
    4081: 'SOLLC',
    3702: 'ARATH',
    9940: 'SHEEP',
    9913: 'BOVIN',
    9925: 'CAPHI',
    44689: 'DICDI',
    4577: 'MAIZE',
    9823: 'PIG',
    9615: 'CANLF',
    6239: 'CAEEL',
    8455: 'XENLA',
    83333: 'ECOLI',
    1891767: 'SV40',
}


mirbase_taxids = {
    9606: 'hsa',
    10090: 'mmu',
    10116: 'rno',
    7227: 'dme',
}


ensembl_taxids = dict(
    (
        t.taxon_id,
        t.ensembl_name
    )
    for t in ensembl_input.ensembl_organisms()
)


nonstandard_taxids = {
    'drosophila': 7227,
    'c.elegans': 6239,
    'xenopus': 8355,
    'Synechocystis_sp.': 1142,
}


def shorten_latin_name(name: str, dot: bool = True) -> str:
    """
    For a complete latin name, returns its shortened version.

    In short latin names the genus name is marked only by its initial.
    """

    if name:

        name = name.split()

        return f'{name[0][0].upper()}{"." if dot else ""} {"".join(name[1:])}'


def short_latin_names(long_names: dict[str, int]) -> dict[str, int]:
    """
    For a dict of long latin names returns a dict with all names shortened.
    """

    return {
        shorten_latin_name(k, dot = dot): v
        for k, v in long_names.items()
        for dot in (True, False)
    }


def ensure_common_name(taxon_id: str | int, lower: bool = False) -> str | None:
    """
    Common English name of an organism.

    Args:
        taxon_id:
            Organism name or NCBI Taxonomy ID.
        lower:
            Return lowercase name. Default is capitalized.
    """

    common_name = (
        # priority for these common names
        taxids.get(taxon_id, None) or
        _ensure_name(taxon_id, 'common')
    )

    if common_name:

        method = 'lower' if lower else 'capitalize'
        common_name = getattr(common_name, method)()

    return common_name


def ensure_latin_name(taxon_id):

    return _ensure_name(taxon_id, 'latin')


def ensure_ensembl_name(taxon_id):

    return _ensure_name(taxon_id, 'ensembl')


def _ensure_name(taxon_id, name_type):

    ncbi_tax_id = ensure_ncbi_tax_id(taxon_id)

    ncbi_to_name = get_db('ncbi_to_%s' % name_type)

    if ncbi_tax_id in ncbi_to_name:

        return ncbi_to_name[ncbi_tax_id]

    _log(
        'Could not find %s taxon name for `%s`.' % (
            name_type,
            str(taxon_id),
        )
    )


def taxid_from_common_name(taxon_name):

    if common.is_str(taxon_name):

        taxon_name = taxon_name.strip()
        taxon_name_l = taxon_name.lower()
        taxon_name_c = taxon_name.capitalize()

    if (
        taxon_name is None or
        not taxon_name_l or
        taxon_name in {'none', 'unknown'}
    ):

        return None

    if taxon_name_l in taxa_synonyms:

        return taxid_from_common_name(taxa_synonyms[taxon_name_l])

    if taxon_name_l in taxa:

        return taxa[taxon_name_l]

    if taxon_name_l in taxa2:

        return taxa2[taxon_name_l]

    common_to_ncbi = get_db('common')

    if taxon_name in common_to_ncbi:

        return common_to_ncbi[taxon_name]

    if taxon_name_c in common_to_ncbi:

        return common_to_ncbi[taxon_name_c]


def taxid_from_latin_name(taxon_name):

    if taxon_name in latin_name_to_ncbi_tax_id:

        return latin_name_to_ncbi_tax_id[taxon_name]

    if taxon_name in short_latin_name_to_ncbi_tax_id:

        return short_latin_name_to_ncbi_tax_id[taxon_name]

    latin_to_ncbi = get_db('latin')

    if taxon_name in latin_to_ncbi:

        return latin_to_ncbi[taxon_name]


def taxid_from_dbptm_taxon_name(taxon_name):

    if taxon_name in dbptm_to_ncbi_tax_id:

        return dbptm_to_ncbi_tax_id[taxon_name]


def taxid_from_nonstandard(taxon_name):

    if taxon_name in nonstandard_taxids:

        return nonstandard_taxids[taxon_name]


def taxid_from_ensembl_name(taxon_name):

    if taxon_name in ensembl_name_to_ncbi_tax_id:

        return ensembl_name_to_ncbi_tax_id[taxon_name]


def ensure_ncbi_tax_id(taxon_id):
    """
    For taxon names of various formats returns NCBI Taxonomy ID if possible.
    Handles English names, scientific names and other common language
    synonyms and database specific codenames.
    """

    if isinstance(taxon_id, int):

        return taxon_id

    else:

        if hasattr(taxon_id, 'strip'):

            taxon_id = taxon_id.strip()

        if common.is_str(taxon_id) and '(' in taxon_id:

            part0, part1 = taxon_id.split('(', maxsplit = 1)

            ncbi_tax_id = (
                ensure_ncbi_tax_id(part0) or
                ensure_ncbi_tax_id(part1.split(')', maxsplit = 1)[0])
            )

        elif hasattr(taxon_id, 'isdigit') and taxon_id.isdigit():

            ncbi_tax_id = int(taxon_id)

        else:

            ncbi_tax_id = (
                taxid_from_dbptm_taxon_name(taxon_id) or
                taxid_from_nonstandard(taxon_id) or
                taxid_from_common_name(taxon_id) or
                taxid_from_latin_name(taxon_id) or
                taxid_from_ensembl_name(taxon_id)
            )

        if not ncbi_tax_id:

            _log('Could not map to NCBI Taxonomy ID: `%s`.' % str(taxon_id))

        return ncbi_tax_id


def uniprot_taxid(uniprot):
    """
    For a UniProt ID returns its NCBI Taxonomy ID.
    """

    uniprot_to_taxid = get_db('swissprot')

    if uniprot in uniprot_to_taxid:

        return uniprot_to_taxid[uniprot]


dbptm_to_ncbi_tax_id = common.swap_dict_simple(dbptm_taxids)
latin_name_to_ncbi_tax_id = common.swap_dict_simple(phosphoelm_taxids)
short_latin_name_to_ncbi_tax_id = short_latin_names(latin_name_to_ncbi_tax_id)
ensembl_name_to_ncbi_tax_id = common.swap_dict_simple(ensembl_taxids)

_cleanup_timeloop = timeloop.Timeloop()
_cleanup_timeloop.logger.setLevel(9999)

@_cleanup_timeloop.job(
    interval = datetime.timedelta(
        seconds = _cleanup_period
    )
)
def _cleanup():

    keys = list(globals()['db'].keys())

    for key in keys:

        if time.time() - globals()['_last_used'][key] > _lifetime:

            _remove(key)

_cleanup_timeloop.start(block = False)


def _remove(key):

    if key in globals()['db']:

        _logger._log(
            'Removing taxonomy data `%s`.' % key
        )
        del globals()['db'][key]

    if key in globals()['_last_used']:

        del globals()['_last_used'][key]


def get_db(key):

    if key not in globals()['db']:

        init_db(key)

    if key in globals()['db']:

        globals()['_last_used'][key] = time.time()

        return globals()['db'][key]

    else:

        return {}


def init_db(key):

    ncbi_data = uniprot_input.uniprot_ncbi_taxids_2()
    this_db = None
    swap = False
    _key = key

    if key.startswith('ncbi_to_'):

        swap = True
        _key = key.rsplit('_', maxsplit = 1)[-1]

    if _key == 'latin':

        this_db = dict(
            (
                taxon.latin,
                taxon.ncbi_id,
            )
            for taxon in ncbi_data.values()
        )

        if not swap:

            this_db.update(short_latin_names(this_db))

    elif _key == 'common':

        this_db = (
            dict(
                (
                    k.capitalize(), v
                )
                for k, v in itertools.chain(
                    iteritems(taxa),
                    iteritems(taxa2)
                )
            )
        )

        this_db.update(
            dict(
                (
                    taxon.english,
                    taxon.ncbi_id,
                )
                for taxon in ncbi_data.values()
                if taxon.english
            )
        )

    elif _key == 'swissprot':

        uniprot_data = uniprot_input.uniprot_taxonomy()
        latin_to_ncbi = get_db('latin')

        this_db = dict(
            (
                swissprot,
                latin_to_ncbi[name],
            )
            for swissprot, names in iteritems(uniprot_data)
            for name in names
            if name in latin_to_ncbi
        )

    elif _key == 'ensembl':

        this_db = ensembl_name_to_ncbi_tax_id

    if swap:

        this_db = common.swap_dict(this_db, force_sets = True)
        this_db = {k: min(v, key = len) for k, v in this_db.items()}

    else:

        this_db.update({k.lower(): v for k, v in this_db.items()})

    if this_db:

        globals()['db'][key] = this_db
        globals()['_last_used'][key] = time.time()
