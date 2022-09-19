from typing import Optional, Union, Literal

from collections import namedtuple
import urllib
from math import ceil

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.common as inputs_common


def oma_pairs(
        genome_id_a: Union[str, int],
        genome_id_b: Union[str, int],
        rel_type: Optional[Literal['1:1', '1:n', 'm:1', 'm:n']] = None,
        score: int = None
    ) -> list[tuple]:
    """
    Downloads Oma data -> pairwise relations among two genomes

    Args:
        genome_id_a: NCBI taxon id of the first genome
        genome_id_b: NCBI taxon id of the second genome
        rel_type: limit relations to a certain type of relations (1:1, 1:n, m:1 or m:n)
        score: resemblance metric

    Returns:
        A list with tuples. tuples are pairing genomes
    """

    result = set() # first decleration is set to prevent recurrency. But at the end it will return as a list 
    OmaPairs = namedtuple(
            'OmaPairs',
            [
                'id_a',
                'id_b',
                'taxon_id_a',
                'taxon_id_b',
                'rel_type',
                'score'
            ],
            defaults = None
        )
    
    url = urls.urls['oma']['url']
    page_url = '%s%s/%s/' % (url,genome_id_a, genome_id_b) # first page url for finding 
    page_size = ceil(float(urllib.request.urlopen(page_url).headers['X-total-count']) / 100) # sample size / 100(each page has)

    for i in range(1, page_size + 1):

        page_url = '%s%s/%s/?page=%u' % (url,genome_id_a, genome_id_b, i)
        c = curl.Curl(page_url, silent = False)
        data = inputs_common.json_read(c.result)

        if score != None:

            data = filter(lambda x: x['score'] >= score != 0, data)

        if rel_type != None:

            data = filter(lambda x: x['rel_type'] == rel_type, data)
            
        for i in data:

            oma_pair = OmaPairs(
                id_a = i['entry_1']['canonicalid'],
                id_b = i['entry_2']['canonicalid'],
                taxon_id_a = i['entry_1']['species']['taxon_id'],
                taxon_id_b = i['entry_2']['species']['taxon_id'],
                rel_type = i['rel_type'],
                score = round(i['score'],2)
            )

            result.add(oma_pair)

    return list(result)