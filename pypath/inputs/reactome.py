import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls


class Reactome:

    def __init__(self, source):

        if source == 'uniprot2all':

            self.url = urls.urls['reactome']['url_uniprot2all']
            self.fields = [
                'uniprot_id',
                'pathway_id',
                'url',
                'pathway_name',
                'evidence_code',
                'organism'
            ]
        
        elif source == 'chebi2all':

            self.url = urls.urls['reactome']['url_chebi2all']
            self.fields = [
                'chebi_id',
                'pathway_id',
                'url',
                'pathway_name',
                'evidence_code',
                'organism'
            ] 
        
        elif source == 'pathways':

            self.url = urls.urls['reactome']['url_pathways']
            self.fields = [
                'pathway_id',
                'pathway_name',
                'organism'
            ]

        elif source == 'pathway_relations':
            
            self.url = urls.urls['reactome']['url_pathway_relations']
            self.fields = [
                'parent',
                'child'
            ]

        else:
            raise ValueError('Invalid Input')


def reactome_data(source) -> list[tuple]:
    '''
    Retrieves pathway information by given input
    Args:
        source:
        - uniprot2all (fields: uniprot_id, pathway_id, url, pathway_name, evidence_code, organism)
        - chebi2all (fields: chebi_id, pathway_id, url, pathway_name, evidence_code, organism)
        - pathways (fields: pathway_id, pathway_name, organism)
        - pathway_relations (fields: parent, child)

    Returns:
        list of namedtuples according to source and its fields.
    '''

    source = Reactome(source = source)
    data = curl.Curl(source.url, large = True).result

    fields = source.fields

    result = set()
    record = collections.namedtuple('Pathway', fields)

    for line in data:

        elements = line.strip('\n').split('\t')
        zipped = dict(zip(fields, elements))

        result.add(
            record(**dict(zip(fields, (zipped.get(f, None) for f in fields))))
        )
    
    return list(result)