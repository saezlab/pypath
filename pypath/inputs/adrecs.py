import gzip
import collections
import pandas as pd

import pypath.share.curl as curl
import pypath.resources.urls as urls


def adrecs_drug_information() -> list[tuple]:
    '''
    Retrieves drug name and synonyms. Drug ids matched in different resources.
    Returns:
        list of namedtuples with all attributes 
        (drug_id, drug_name, drugbank_id, pubchem_id, mesh_id, kegg_id, ttd_id)
    '''

    url = urls.urls['adrecs']['url_drug_information']
    path = curl.Curl(
        url,
        silent = False,
        large = True
    ).fileobj.name

    with gzip.open(path, 'rb') as f_in:

        df = pd.read_excel(f_in)
        df.replace({'Not Available': None}, inplace=True)
        df['DRUG_SYNONYMS'] =  df['DRUG_SYNONYMS'].str.split('|')
        df['DRUG_SYNONYMS'] =  df['DRUG_SYNONYMS'].apply(lambda x: tuple([element.strip() for element in x]) if x is not None else None)
        df.columns= df.columns.str.lower()
        
        return list(df.itertuples(name='AdrecsDrugInfo', index=False))


def adrecs_terminology() -> list[tuple]:
    '''
    Retrieves ADR terminology and hierachy.
    Returns:
        list of namedtuples with all attributes 
        (adrecs_id, adr_id, adr_term, adr_synonyms, meddra_code)
    '''

    url = urls.urls['adrecs']['url_terminology']
    path = curl.Curl(
        url,
        silent = False,
        large = True
    ).fileobj.name

    with gzip.open(path, 'rb') as f_in:

        df = pd.read_excel(f_in)
        df.replace({'Not Available': None}, inplace=True)
        df['ADR_SYNONYMS'] = df['ADR_SYNONYMS'].str.split('|')
        df['ADR_SYNONYMS'] = df['ADR_SYNONYMS'].apply(lambda x: tuple([element.strip() for element in x]) if x is not None else None)
        df.columns= df.columns.str.lower()
        
        return list(df.itertuples(name='AdrecsOntology', index=False))


def adrecs_drugs() -> list[tuple]:
    '''
    Retrieves Drug-ADR pairs
    Returns:
        list of namedtuples with drug_id and adr_id
    '''
    
    url = urls.urls['adrecs']['url_adrecs_drugs']
    c = curl.Curl(url, large = True, silent = False)
    drugs = c.result
    
    fields = [
        'drug_id',
        'adr_id'
    ]

    result = set()
    record = collections.namedtuple('AdrecsDrug', fields)

    # requisite fields from all_fields
    indices = [0, 1]

    for line in drugs:

        if not line.strip():
            continue
        
        line = line.strip().split('\t')
        line = dict(zip(fields, (line[i] for i in indices)))

        result.add(
            record(**dict(zip(fields, (line.get(f, None) for f in fields))))
        )

    return list(result)