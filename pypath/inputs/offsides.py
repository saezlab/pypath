import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls


def offsides() -> list[tuple]:
    """
    Retrieves individual drug side effect signals mined from the FDA's Adverse Event Reporting System
    Returns:
        side effect information as list of tuples
    """

    url = urls.urls['offsides']['url']
    c = curl.Curl(url, large = True, silent = False)
    drug_names = c.result

    fields = [
        'drug_rxnorn_id', 
        'drug_concept_name', 
        'condition_meddra_id', 
        'condition_concept_name', 
        'PRR', 
        'PRR_error', 
        'mean_reporting_frequency'
    ]

    result = set()
    record = collections.namedtuple('Offside', fields)

    # requisite fields from all_fields
    indices = [0, 1, 2, 3, 8, 9, 10]

    for line in drug_names:
        if not line.strip():
            continue
        
        line = line.strip().split(',')
        line = dict(zip(fields, (line[i] for i in indices)))

        result.add(
            record(**dict(zip(fields, (line.get(f, None) for f in fields))))
        )
    
    return list(result)
