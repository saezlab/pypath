#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#           Tennur Kılıç
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from typing import Optional

import re
import csv
import collections
import base64

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.share.settings as settings
import pypath.inputs.credentials as credentials

_logger = session.Logger(name = 'drugbank_input')
_log = _logger._log


def _drugbank_credentials(
        user: Optional[str] = None,
        passwd: Optional[str] = None,
        credentials_fname: Optional[str] = None,
    ) -> tuple[str, str]:

    return credentials.credentials(
        user = user,
        passwd = passwd,
        resource = 'DrugBank',
        from_file = credentials_fname,
    )


def _drugbank_download(
        *args,
        user: Optional[str] = None,
        passwd: Optional[str] = None,
        credentials_fname: Optional[str] = None,
        **kwargs
    ) -> Optional[curl.Curl]:

    try:

        cred = _drugbank_credentials(
            user = user,
            passwd = passwd,
            credentials_fname = credentials_fname,
        )

    except RuntimeError:

        _log('No credentials available for the DrugBank website.')

        return None

    defaults = {
        'large': True,
        'silent': False,
        'compr': 'zip',
    }

    defaults.update(kwargs)

    auth_str = base64.b64encode(f"{cred['user']}:{cred['passwd']}".encode())

    defaults['req_headers'] = [
        f'Authorization: Basic {auth_str.decode()}',
        settings.get('user_agent'),
    ]

    return curl.Curl(*args, **defaults)


def drugbank_raw_interactions(
        user: Optional[str] = None,
        passwd: Optional[str] = None,
        credentials_fname: Optional[str] = None,
        pharma_active: bool = False,
    ) -> list[tuple] :
    """
    Retrieves protein identifiers from Drugbank.

    Args:
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.
        pharma_active:
            Only pharmacologically active relations.

    Returns:
        List of drug-protein relations.
    """

    csv_name = 'pharmacologically_active.csv' if pharma_active else 'all.csv'

    fields = (
        'drugbank_id',
        'uniprot_id',
        'relation',
    )

    DrugbankRawInteraction = collections.namedtuple(
        'DrugbankRawInteraction',
        fields,
        defaults = (None,) * len(fields),
    )

    result = []

    for rel in ('carrier', 'enzyme', 'target', 'transporter'):

        url = urls.urls['drugbank'][f'drug_{rel}_identifiers']

        c = _drugbank_download(
            url = url,
            user = user,
            passwd = passwd,
            credentials_fname = credentials_fname,
            files_needed = (csv_name,),
        )

        if not c: continue

        _ = next(c.result[csv_name])

        for l in c.result[csv_name]:

            drugs, uniprot = l.strip().split(',')[-1], l.strip().split(',')[5]

            drugs = drugs.strip().split(';')

            result.extend(
                DrugbankRawInteraction(
                    drugbank_id = drug.strip(),
                    uniprot_id = uniprot,
                    relation = rel,
                )
                for drug in drugs
            )

    return result


def drugbank_interactions(
        user: Optional[str] = None,
        passwd: Optional[str] = None,
        credentials_fname: Optional[str] = None,
        pharma_active: bool = False,
    ) -> list[tuple] :
    """
    Drug-protein and protein-drug interactions from Drugbank.

    Args:
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.
        pharma_active:
            Only pharmacologically active interactions.

    Returns:
        List of drug-protein and protein-drug interactions.
    """

    raw = drugbank_raw_interactions(
        user = user,
        passwd = passwd,
        pharma_active = pharma_active,
        credentials_fname = credentials_fname,
    )

    drugs = dict(
        (d.drugbank, d)
        for d in drugbank_drugs(user = user, passwd = passwd)
    )

    DrugbankInteraction = collections.namedtuple(
        'DrugbankInteraction',
        (
            'source',
            'target',
            'source_entity_type',
            'target_entity_type',
            'interaction_type',
        )
    )

    result = []

    for r in raw:

        drug = drugs.get(r.drugbank_id, None)

        # TODO: later engage the mapping module here
        if drug and drug.pubchem_cid:

            src_tgt = reversed if r.relation == 'target' else lambda x: x

            result.append(
                DrugbankInteraction(
                    *src_tgt((r.uniprot_id, drug.pubchem_cid)),
                    *src_tgt(('protein', 'drug')),
                    interaction_type = r.relation,
                )
            )

    return result


def drugbank_drugs(
        user: Optional[str] = None,
        passwd: Optional[str] = None,
        credentials_fname: Optional[str] = None,
    ) -> list[tuple]:
    """
    Retrieves drug identifiers from Drugbank.

    Each drug is annotated by its various database cross-references.

    Args:
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.

    Returns:
        List of named tuples, each field corresponding to various identifiers.
    """

    fields = (
        'drugbank',
        'name',
        'type',
        'groups',
        'cas',
        'inchikey',
        'inchi',
        'smiles',
        'formula',
        'kegg_compound',
        'kegg_drug',
        'pubchem_cid',
        'pubchem_sid',
        'chebi',
        'chembl',
        'pharmgkb',
        'het',
    )

    raw = {}

    for table in ('drug', 'structure'):

        csv_ = f'{table} links.csv'

        c = _drugbank_download(
            url = urls.urls['drugbank'][f'all_{table}s'],
            user = user,
            passwd = passwd,
            credentials_fname = credentials_fname,
            files_needed = (csv_,),
        )

        if not c: continue

        raw[table] = dict(
            (rec['DrugBank ID'], rec)
            for rec in csv.DictReader(c.result[csv_], delimiter = ',')
        )

    DrugbankDrug = collections.namedtuple(
        'DrugbankDrug',
        fields,
        defaults = (None,) * len(fields),
    )

    result = []

    for dbid, struct in raw['structure'].items():

        drug = raw['drug'].get(dbid, {})

        result.append(
            DrugbankDrug(
                drugbank = dbid,
                name = struct['Name'],
                type = drug.get('Drug Type', None),
                groups = struct['Drug Groups'],
                cas = struct['CAS Number'],
                inchikey = struct['InChIKey'],
                inchi = struct['InChI'],
                smiles = struct['SMILES'],
                formula = struct['Formula'],
                kegg_compound = struct['KEGG Compound ID'],
                kegg_drug = struct['KEGG Drug ID'],
                pubchem_cid = struct['PubChem Compound ID'],
                pubchem_sid = struct['PubChem Substance ID'],
                chebi = struct['ChEBI ID'],
                chembl = struct['ChEMBL ID'],
                pharmgkb = drug.get('PharmGKB ID', None),
                het = drug.get('HET ID', None),
            )
        )

    return result


def drugbank_annotations(
        user: Optional[str] = None,
        passwd: Optional[str] = None,
        credentials_fname: Optional[str] = None,
    ) -> dict[str, set[tuple]]:
    """
    Drug annotations from Drugbank.

    The annotations are restricted to the drug molecule type and drug status.

    Args:
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.

    Returns:
        List of drug annotations.
    """

    drugs = drugbank_drugs(
        user = user,
        passwd = passwd,
        credentials_fname = credentials_fname,
    )

    DrugbankAnnotation = collections.namedtuple(
        'DrugbankAnnotation',
        (
            'type',
            'status',
        )
    )

    result = collections.defaultdict(set)

    for d in drugs:

        if d.pubchem_cid:

            result[d.pubchem_cid].add(
                DrugbankAnnotation(
                    type = d.type,
                    status = re.sub(',\s*', ';', d.groups),
                )
            )

    return dict(result)


def drugbank_mapping(
        id_type: str,
        target_id_type: str,
        user: Optional[str] = None,
        passwd: Optional[str] = None,
        credentials_fname: Optional[str] = None,
    ) -> dict[str, set[str]]:
    """
    Identifier translation table from DrugBank.

    Available ID types: drugbank, name, type, groups, cas, inchikey,
    inchi, smiles, formula, kegg_compound, kegg_drug, pubchem_compound,
    pubchem_substance, chebi, chembl, pharmgkb, het.

    Args:
        id_type:
            The identifier type to be used as keys.
        target_id_type:
            The identifier type that will be collected into the values.
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.
        credentials_fname:
            File name or path to a file with DrugBank login credentials.

    Returns:
        An identifier translation table.
    """

    synonyms = {
        'pubchem_compound': 'pubchem_cid',
        'pubchem_substance': 'pubchem_sid',
    }


    def id_type_proc(_id_type):

        _id_type = re.sub('[^cs]id$', '', _id_type.lower()).replace(' ', '_')

        return synonyms.get(_id_type, _id_type)


    drugs = drugbank_drugs(
        user = user,
        passwd = passwd,
        credentials_fname = credentials_fname,
    )

    result = collections.defaultdict(set)

    id_type = id_type_proc(id_type)
    target_id_type = id_type_proc(target_id_type)

    for d in drugs:

        the_id = getattr(d, id_type)
        target_id = getattr(d, target_id_type)

        if the_id and target_id:

            result[the_id].add(target_id)

    return dict(result)
