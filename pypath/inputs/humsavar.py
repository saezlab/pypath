import collections 

import pypath.resources.urls as urls
import pypath.share.curl as curl


def parse(description):
        name = None
        symbol = None
        omim_id = None

        if ('(' in description) and ('[' in description):
            name = description[:description.find("(") - 1]
            symbol = description[description.find("(")+1:description.find(")")]
            omim_id = description[description.find("[")+1:description.find("]")]
        
        elif ('(' in description) and ('[' not in description):
            name = description[:description.find("(") - 1]
            symbol = description[description.find("(")+1:description.find(")")]

        elif ('(' not in description) and ('[' in description):
            name = description[:description.find("[") - 1]
            omim_id = description[description.find("[")+1:description.find("]")]

        return (name, symbol, omim_id)


def humsavar() -> list[tuple]:
    '''
    Retrieves all missense variants annotated in UniProtKB/Swiss-Prot human entries.

    Returns:
        Drug attributes in below as a list of named tuples.
    '''
    Variant = collections.namedtuple(
            'Variant',
            [   
                'main_gene_name',
                'swiss_prot_ac',
                'ftid',
                'AA_change',
                'variant_category',
                'dbSNP',
                'disease_name',
                'disease_symbol',
                'disease_omim_id'
            ],
            defaults=None
        )

    url = urls.urls['humsavar']['url']

    data = curl.Curl(url, large=True, silent=False)

    data = data.result

    result = set()

    # skipping data description information in txt file
    for r in data:
        if r.startswith('_'):
            break

    for row in data:
        description = ''
        info = row.split()

        # make sure not only contains first item '-'
        if len(info) > 1:

            # make sure swiss_prot_ac does not mix the main_gene_name because
            # some of them concat with main gene name unsensibly
            if len(info[1]) == 6:

                description = ' '.join(info[6:])
            
                name, symbol, omim_id = parse(description)

                variant = Variant(
                        main_gene_name = info[0],
                        swiss_prot_ac = info[1],
                        ftid = info[2],
                        AA_change = info[3],
                        variant_category=info[4],
                        dbSNP=info[5],
                        disease_name= name,
                        disease_symbol=symbol,
                        disease_omim_id=omim_id
                )

                result.add(variant)

    return list(result)
