import collections

from SPARQLWrapper import SPARQLWrapper, JSON

import pypath.resources.urls as urls


def disease_pathogen_interactions():
    """
        Downloads PathoPhenoDb data -> disease pathogen pairs
        Returns:
        A list with tuples. Tuples are pairing disease and pathogen
    """

    result = set()

    DiseasePathogen = collections.namedtuple(
        'DiseasePathogen',
        [
            'disease_id',
            'disease_name',
            'pathogen_tax_id',
            'pathogen_name',
            'evidence_type'
        ],
        defaults=None
    )

    url = urls.urls['pathophenodb']["url"]
    sparql = SPARQLWrapper(url)
    sparql.setReturnFormat(JSON)
    sparql.setQuery("""
            #EX3:List all diseases which caused by pathogens
            PREFIX SIO: <http://semanticscience.org/resource/SIO_>
            PREFIX RO: <http://purl.obolibrary.org/obo/RO_>
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            SELECT distinct ?Disease_ID ?Disease ?Pathogen_ID ?Pathogen ?evidence_Code    
            FROM <http://patho.phenomebrowser.net>
            WHERE 
            {
                ?Disease_ID SIO:000255 ?o .  
                ?o RO:0002558 ?o1 .
                ?o RO:0002556  ?Pathogen_ID .
                ?Disease_ID rdfs:label ?Disease .
                ?Pathogen_ID rdfs:label ?Pathogen .
                ?o1 rdfs:label ?evidence_Code .
            } 
        """
    )

    ret = sparql.queryAndConvert()

    for r in ret["results"]["bindings"]:
        pair = DiseasePathogen(
            disease_id = r['Disease_ID']['value'].split("/")[-1].replace('_',':'),
            disease_name = r['Disease']['value'],
            pathogen_tax_id = r['Pathogen_ID']['value'].split("/")[-1].split('_')[1],
            pathogen_name = r['Pathogen']['value'],
            evidence_type = r['evidence_Code']['value']
        )
        result.add(pair)
    
    return list(result)
    