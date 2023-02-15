import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls


def twosides() -> list[tuple]:
    """
    Retrieves drug-drug interaction safety signals mined from the FDA's Adverse Event Reporting System
    Returns:
        interaction information as list of tuples
    """

    url = urls.urls['twosides']['url']
    c = curl.Curl(url, large = True, silent = False)
    drugs = c.result

    fields = [
        'drug_1_rxnorn_id', 
        'drug_1_concept_name',
        'drug_2_rxnorm_id',
        'drug_2_concept_name',
        'condition_meddra_id', 
        'condition_concept_name', 
        'PRR', 
        'PRR_error', 
        'mean_reporting_frequency'
    ]

    result = set()
    record = collections.namedtuple('TwoSides', fields)

    # requisite fields from all_fields
    indices = [0, 1, 2, 3, 4, 5, 10, 11, 12]

    for line in drugs:

        if not line.strip():
            continue
        
        line = line.strip().split(',')
        line = dict(zip(fields, (line[i] for i in indices)))

        result.add(
            record(**dict(zip(fields, (line.get(f, None) for f in fields))))
        )
    
    return list(result)
