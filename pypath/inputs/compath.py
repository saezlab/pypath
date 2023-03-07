import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl

def kegg_control(dct):

    if dct['source_db'] == 'kegg':
        dct['id_1'] = dct['id_1'][5:]

    elif dct['target_db'] == 'kegg':
        dct['id_2'] = dct['id_2'][5:]

    return dct


def get_all_mappings(source_db = None, target_db = None) -> list[tuple]:
    """
    Retrieves proposed and accepted mappings by the users/curators between a pair of pathways.
    Parameters:
        source_db & target_db are database names. It specifies direction 
    Returns:
        interaction information as list of tuples
    """

    url = urls.urls['compath']['url']
    c = curl.Curl(url, large=True).result

    allowed_types = ['kegg', 'wikipathways', 'reactome', None]

    if (source_db not in allowed_types) or (target_db not in allowed_types):
        raise ValueError('Invalid Input')

    result = set()
    fields = [
            'pathway1',
            'id_1',
            'source_db',
            'relation',
            'pathway2',
            'id_2',
            'target_db',     
        ]
    record = collections.namedtuple('Mapping', fields)


    for mapping in c:

        elements = mapping.split('\t')
        elements[-1] = elements[-1].strip('\n')

        zipped = dict(zip(fields, elements))
        zipped = kegg_control(zipped)

        if (source_db is None) and (target_db is None):

            result.add(
                record(**dict(zip(fields, (zipped.get(f, None) for f in fields))))
            )

        
        elif (source_db is None) and (target_db is not None):

            if zipped['target_db'] == target_db:
                
                result.add(
                    record(**dict(zip(fields, (zipped.get(f, None) for f in fields))))
                )
        
        elif (source_db is not None) and (target_db is None):

            if zipped['source_db'] == source_db:

                result.add(
                    record(**dict(zip(fields, (zipped.get(f, None) for f in fields))))
                )

        else:
            if (zipped['source_db'] == source_db) and (zipped['target_db'] == target_db):

                result.add(
                    record(**dict(zip(fields, (zipped.get(f, None) for f in fields))))
                )

    return list(result)

