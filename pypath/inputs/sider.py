import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common


def sider_drug() -> list[tuple]:
    """
    Retrieves drug information from the Sider database.
    Returns:
        Drug cid, name and atc information as a list of named tuples.
    """

    Drug = collections.namedtuple(
        'Drug',
        [   
            'cid',
            'drug_name',
            'drug_atc',
        ],
        defaults=None
    )

    result = set()
    dct = {}

    url = urls.urls['sider']['url_drug_names']
    c = curl.Curl(url, large = True, silent = False)

    drug_names = c.result
    
    for drug in drug_names:

        q = drug.strip('\n').split('\t')
        cid = q[0]
        name = q[1]

        if cid not in dct:
            dct[cid] = {'name':name, 'atc':None}

    url_drug_atc = urls.urls['sider']['url_drug_atc']
    c2 = curl.Curl(url_drug_atc, large = True, silent = False)

    drug_atc = c2.result

    for drug in drug_atc:

        q = drug.strip('\n').split('\t')
        cid = q[0]
        atc = q[1]
        
        if cid in dct:
            if dct[cid]['atc'] is not None:
                dct[cid]['atc'].append(atc)
            else:
                dct[cid]['atc'] = [atc]

    for key in dct:

        drug = Drug(
            cid = key,
            drug_name = dct[key]['name'],
            drug_atc = tuple(dct[key]['atc']) if dct[key]['atc'] is not None else None
        )

        result.add(drug)

    return list(result)


def sider_meddra() -> list[tuple]:
    """
    Retrieves drug information from the Sider database.
    Returns:
        Drug cid, umls concept ids both for label and MedDra and side effect name.
    """

    fields = (
            'cid',
            'umls_concept_id_on_label',
            'umls_concept_id_on_MedDRA',
            'side_effect_name',
    )
       
    fields = common.to_list(fields)

    url_meddra_all = urls.urls['sider']['url_meddra_all']

    c = curl.Curl(
        url_meddra_all, 
        large = True, 
        silent = False
    )

    result = set()
    record = collections.namedtuple('Drug', fields)

    # essential features' indices
    indices = [0, 2, 4, 5]

    for line in c.result:

        if not line.strip():
            continue
        
        line = line.strip().split('\t')
        line = dict(zip(fields, (line[i] for i in indices)))

        result.add(
            record(**dict(zip(fields, (line.get(f, None) for f in fields))))
        )

    return list(result)


def sider_meddra_with_freq() -> list[tuple]:
    """
    Retrieves drug information from the Sider database.
    Returns:
        Drug cid, umls concept ids both for label and MedDra,
        frequency information and side effect name.
    Attention! -> sider_meddra_all() function has about 20k row bigger 
    than this dataset without frequency information
    """

    fields = (
            'cid',
            'umls_concept_id_on_label',
            'frequency',
            'umls_concept_id_on_MedDRA',
            'side_effect_name',
    )

    fields = common.to_list(fields)

    url_meddra_with_freq = urls.urls['sider']['url_meddra_freq']

    c = curl.Curl(
        url_meddra_with_freq, 
        large = True, 
        silent = False
    )

    result = set()
    record = collections.namedtuple('Drug', fields)

    # essential features' indices
    indices = [0, 2, 4, 8, 9]

    for line in c.result:

        if not line.strip():
            continue
        
        line = line.strip().split('\t')
        line = dict(zip(fields, (line[i] for i in indices)))

        result.add(
            record(**dict(zip(fields, (line.get(f, None) for f in fields))))
        )

    return list(result)