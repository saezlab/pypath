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

from typing import Optional

import os
import re
import csv
import collections
import base64
from lxml import etree
from zipfile import ZipFile

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
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

    Args
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.
        pharma_active:
            Only pharmacologically active relations.

    Returns
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

    Args
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.
        pharma_active:
            Only pharmacologically active interactions.

    Returns
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

    Args
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.

    Returns
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

    Args
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.

    Returns
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

    Args
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

    Returns
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


class DrugbankFull:
    """
    This is a wrapper around the Drugbank full database XML file.
    Provides access to the full Drugbank database.
    The class provides two methods: drugbank_drugs_full and drugbank_targets_full.
    The first method returns a list of namedtuples, each of which represents a drug.
    The second method returns a list of namedtuples, each of which represents a drug's target.

    Args
        user:
            E-mail address with registered DrugBank account.
        passwd:
            Password for the DrugBank account.
    """

    def __init__(
        self, 
        user: Optional[str] = None,
        passwd: Optional[str] = None,
        credentials_fname: Optional[str] = None,
    ):

        path = _drugbank_download(
            url = urls.urls['drugbank']['full_database'],
            user = user,
            passwd = passwd,
            credentials_fname = credentials_fname,
            ).fileobj.name

        with ZipFile(path, 'r') as zip_ref:
            zip_ref.extractall(os.path.dirname(path))

        file = os.path.join(os.path.dirname(path), 'full database.xml')

        self.tree = etree.ElementTree(file = file)
        self.ns = self.tree.getroot().nsmap
        self.ns['db'] = self.ns[None]
        del self.ns[None]

        self.drugs = self.tree.xpath('db:drug', namespaces=self.ns)


    def drugbank_drugs_full(
            self,
            fields: str | list[str] | None = None,
        ) -> list[tuple]:
        """
        Returns a list of namedtuples containing detailed information about drugs.

        Args
        fields:
            The fields to return. If None, all XML fields are returned.
            Default: None

        Returns
            A list of namedtuples containing information about drugs.
        """

        basic_fields = [
            'drugbank_id', 'type', 'name', 'description', 'cas_number', 'unii', 'average_mass',
            'monoisotopic_mass', 'state', 'synthesis_reference', 'indication', 'pharmacodynamics', 
            'mechanism_of_action', 'toxicity', 'metabolism', 'absorption', 'half_life',
            'protein_binding', 'route_of_elimination', 'volume_of_distribution', 'clearance',
            'fda_label', 'msds',
        ]

        fields_w_subfields = {
            'groups':  {'path': '/db:group'},
            'general_references': {'path': '/db:articles/db:article/db:pubmed-id'}, 
            'classification': {'path': '/db:class'},
            'synonyms': {'path': '/db:synonym'},
            'products': {'path': '/db:product/db:name'},
            'international_brands': {'path': '/db:international-brand/db:name'},
            'mixtures': {'path': '/db:mixture/db:name'},
            'packagers': {'path': '/db:packager/db:name'},
            'manufacturers': {'path': '/db:manufacturer/db:name'},
            'categories': {'path': '/db:category/db:mesh-id'},
            'affected_organisms': {'path': '/db:affected-organism'},
            'atc_codes': {'path': '/db:atc-code', 'key': 'code'},
            'ahfs_codes': {'path': '/db:ahfs-code', 'key': 'code'},
            'pdb_entries': {'path': '/db:pdb-entry'},
            'patents': {'path': '/db:patent/db:number'},
            'food_interactions': {'path': '/db:food-interaction'},
            'drug_interactions': {'path': '/db:drug-interaction/db:drugbank-id'},
            'pathways': {'path': '/db:pathway/db:smpdb-id'},
        }

        # TODO: later process and engage fields below
        # future_fields: 'salts', 'prices', 'dosages', 'sequences',
        #               'experimental_properties', 'external_links',
        #               'reactions', 'snp_effects', 'snp_adverse_drug_reactions'
            
        fields = fields or basic_fields + list(fields_w_subfields.keys())
        fields = common.to_list(fields)
        if 'drugbank_id' not in fields:
            fields.insert(0, 'drugbank_id')

        result = []

        record = collections.namedtuple('DrugbankDrug', fields)

        for drug in self.drugs:

            field_dict = {}

            for field in fields:
                
                if field == 'drugbank_id':
                    field_dict[field] = [i for i in drug.xpath('db:drugbank-id', namespaces=self.ns) if i.attrib.get('primary') == 'true'][0].text

                elif field == 'type':
                    field_dict[field] = drug.get('type')

                else:
                    if field in fields_w_subfields:
                        path_to_field = f"db:{field.replace('_', '-')}{fields_w_subfields[field]['path']}"

                        if 'key' in fields_w_subfields[field]:                    
                            field_dict[field] = {f.get(fields_w_subfields[field]['key']) for f in drug.xpath(path_to_field, namespaces=self.ns)}
                        
                        else:
                            field_dict[field] = {f.text for f in drug.xpath(path_to_field, namespaces=self.ns)}

                    else:
                        path_to_field = f"db:{field.replace('_', '-')}"
                        field_dict[field] = {f.text for f in drug.xpath(path_to_field, namespaces=self.ns)}
                
                
            for k, v in field_dict.items():

                if v and type(v) != str:
                    field_dict[k] = [elem.replace('\r\n', ' ') for elem in v if elem]

                    if len(field_dict[k]) == 1:
                        field_dict[k] = field_dict[k][0]

                if not field_dict[k]:
                    field_dict[k] = None

            result.append(record(**field_dict))
        
        return result


    def drugbank_targets_full(
            self,
            fields: str | list[str] | None = None,
        ) -> list[tuple]:
        """
        Returns a list of namedtuples containing detailed information about drug-target interactions.

        Args
            fields:
                The fields to return.
                Default: None

        Returns
            A list of namedtuples containing information about the target of drugs.
        """

        result = []

        all_fields = [
            'drugbank_id',
            'id',
            'name',
            'organism',
            'actions',
            'references',
            'known_action',
            'polypeptide',
        ]
        
        fields = fields or all_fields
        fields = common.to_list(fields)
        if 'drugbank_id' not in fields:
            fields.insert(0, 'drugbank_id')

        record = collections.namedtuple('DrugbankTarget', fields)

        
        for drug in self.drugs:

            db_id = [i for i in drug.xpath('db:drugbank-id', namespaces=self.ns) if i.attrib.get('primary') == 'true'][0].text

            for target in drug.xpath('db:targets/db:target', namespaces=self.ns):

                target_dict = {}
                target_dict['drugbank_id'] = db_id
                for field in fields:

                    if field in ['id', 'name', 'organism', 'known_action']:
                        target_dict[field] = [f.text for f in target.xpath(f"db:{field.replace('_', '-')}", namespaces=self.ns)]

                    elif field == 'actions':
                        target_dict[field] =  [f.text for f in target.xpath('db:actions/db:action', namespaces=self.ns)]

                    elif field == 'references':
                        target_dict[field] =  [f.text for f in target.xpath('db:references/db:articles/db:article/db:pubmed-id', namespaces=self.ns)]
                    
                    elif field == 'polypeptide':
                        target_dict[field] = [(f.get('id'), f.get('source')) for f in target.xpath('db:polypeptide', namespaces=self.ns)]

                for k, v in target_dict.items():

                    if v and len(v) == 1:
                        target_dict[k] = v[0]
                    
                    if not v:
                        target_dict[k] = None

                result.append(record(**target_dict))
 
        return result


    def drugbank_external_ids_full(
            self,
        ) -> dict[str, dict]:
        """
            Returns a dictionary containing all external identifiers of drugs.
        """

        result = {}

        for drug in self.drugs:

            db_id = [i for i in drug.xpath('db:drugbank-id', namespaces=self.ns) if i.attrib.get('primary') == 'true'][0].text

            for ext_id in drug.xpath('db:external-identifiers/db:external-identifier', namespaces=self.ns):
                source = ext_id.xpath('db:resource', namespaces=self.ns)[0].text
                identifier = ext_id.xpath('db:identifier', namespaces=self.ns)[0].text

                if db_id not in result:
                    result[db_id] = {}

                result[db_id][source] = identifier

        return result


    def drugbank_properties_full(
            self,
        ) -> dict[str, dict]:
        """
            Returns a dictionary containing calculated properties of drugs.
        """

        result = {}

        for drug in self.drugs:

            db_id = [i for i in drug.xpath('db:drugbank-id', namespaces=self.ns) if i.attrib.get('primary') == 'true'][0].text

            for prop in drug.xpath('db:calculated-properties/db:property', namespaces=self.ns):
                kind = prop.xpath('db:kind', namespaces=self.ns)[0].text
                identifier = prop.xpath('db:value', namespaces=self.ns)[0].text

                if db_id not in result:
                    result[db_id] = {}

                result[db_id][kind] = identifier

        return result