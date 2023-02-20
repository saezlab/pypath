from collections import namedtuple
from tqdm.notebook import tqdm

import json
import bs4
import requests
import json

import pypath.share.curl as curl
import pypath.resources.urls as urls


def get_record_number():
    '''
    Retrieves number of records in ddinter datasource url
    '''
    cookies = {
        'csrftoken': '975RT1K1FXpLX6wZdcyTMIRdOAMwNNkx1MNFAbeJrlrRBFI5DbNmaQPz7tj2UNFZ',
    }

    data = {
        'csrfmiddlewaretoken': 'ToKztkgeCHgBfAbezQL7GULrQABtRmHgL3snauKWo5iHT9nkZP0A42JN9t8ZYm2I',
        'draw': '1',
        'columns[0][data]': '',
        'columns[0][name]': '',
        'columns[0][searchable]': 'true',
        'columns[0][orderable]': 'false',
        'columns[0][search][value]': '',
        'columns[0][search][regex]': 'false',
        'start': '0',
        'length': '24',
        'search[value]': '',
        'search[regex]': 'false',
        'custom-length': '25',
    }

    url = urls.urls['ddinter']['url_source']
    response = requests.post(url, cookies=cookies, data=data, verify=False)

    return int(response.json()["recordsTotal"])


def get_mappings(drug):
    '''
    Retrieves mappings of input drug in drugbank, chembl, pubchem
    Returns:
        list of mapping namedtuple (single)
    '''
    url = urls.urls['ddinter']['url_mapping'] % drug

    c = curl.Curl(
        url, 
        silent = False, 
        large = True, 
        )

    soup = bs4.BeautifulSoup(c.fileobj, 'html.parser')
    refs = soup.find_all('a')

    links = [link.get('href', '') for link in refs]
    result = set()
    
    Mapping = namedtuple(
        'Mapping',
        [
            'drugbank',
            'chembl',
            'pubchem'
        ],
        defaults=None
    )
    mapping_targets = ["drugbank", "chembl", "pubchem"]
    mapping_dict = {}

    for link in links:
        for target in mapping_targets:
            if target in link:
                mapping_dict[target] = link.split("/")[-1]
    
    result.add(Mapping(**mapping_dict))
                
    return list(result)


def check_hashable(data):

    if isinstance(data, (dict, list, set)):

        return tuple(data)

    return data

def get_interactions(drug):
    '''
    Retrieves interactions of input drug
    Returns:
        list of interaction namedtuple with drug ids, drug names, interaction level and actions
    '''
    url = urls.urls['ddinter']['url_interaction'] % drug
    c = curl.Curl(url)
    data = json.loads(c.result)

    result = set() # first decleration is set to prevent recurrency. But at the end it will return as a list 
    Interaction = namedtuple(
            'DDInterInteraction',
            [
                'drug1_id',
                'drug1_name',
                'drug2_id',
                'drug2_name',
                'level',
                'actions'
            ],
            defaults = None
        )

    drug1_id = data['info']['id']
    drug1_name = data['info']['Name']
    
    for interaction_fe in data['interactions']:

        interaction = Interaction(
            drug1_id = drug1_id,
            drug1_name = drug1_name,
            drug2_id = check_hashable(interaction_fe['id']),
            drug2_name =  check_hashable(interaction_fe['name']),
            level = check_hashable(interaction_fe['level']),
            actions = check_hashable(interaction_fe['actions'])
        )

        result.add(interaction)
            
    return list(result)


def get_all_interactions():
    '''
    Retrieves all interactions
    Returns:
        list of interaction namedtuple with drug ids, drug names, interaction level and actions
    '''
    result = set() 
    
    Interaction = namedtuple(
            'DDInterInteraction',
            [
                'drug1_id',
                'drug1_name',
                'drug2_id',
                'drug2_name',
                'level',
                'actions'
            ],
            defaults = None
            )
    
    for index in tqdm(range(1, get_record_number()+1)):

        ddinter_drug = 'DDInter'+str(index)

        url = urls.urls['ddinter']['url_interaction'] % ddinter_drug

        c = curl.Curl(url)
        data = json.loads(c.result)

        drug1_id = data['info']['id']
        drug1_name = data['info']['Name']
    
        for interaction_fe in data['interactions']:

            interaction = Interaction(
                drug1_id = drug1_id,
                drug1_name = drug1_name,
                drug2_id = check_hashable(interaction_fe['id']),
                drug2_name =  check_hashable(interaction_fe['name']),
                level = check_hashable(interaction_fe['level']),
                actions = check_hashable(interaction_fe['actions'])
            )

            result.add(interaction)
            
    return list(result)


