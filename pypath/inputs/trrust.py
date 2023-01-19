
from pypath.share import curl
from pypath.resources.urls import urls
import re

def trrust_scraping(org):

    print(f'Fetching records for {org}')
    choices = {'human', 'mouse'}
    
    if org not in choices:
        print('Not yet supported')
        exit(0)

    url = urls['trrust']['scraping_url'] % org

    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding="utf-8",
        default_mode="r",
    )

    html_generator = c.result
    records = list()

    regex = r'[^<]*<td[^>]*>([^<]*)<\/td>'
    matcher = re.compile(regex)

    for line in html_generator:

        attributes = matcher.findall(line)

        if not attributes:
            continue

        attributes = [e if e != '-' else None for e in attributes]

        try:
            aliases = attributes[2].split(',')
            aliases = [e.strip() for e in aliases]
            aliases = tuple(aliases)
        except AttributeError:
            aliases = None
        
        record = {
            'gene_symbol': attributes[0],
            'entrez_id': attributes[1],
            'aliases': aliases,
            'full_name': attributes[3]
        }

        records.append(record)

    return records

def trrust_general(org):
    
    choices = {'human', 'mouse'}
    
    if org not in choices:
        print('Not yet supported')
        exit(0)
    
    url = urls['trrust']['tsv_url'] % org

    c = curl.Curl(
        url,
        silent=False,
        large=True,
        encoding="utf-8",
        default_mode="r",
    )

    interactions = dict()

    for line in c.result:

        line = line.strip('\n ')
        data = line.split("\t")

        key_gene = data[0]

        interaction = {
            'gene_symbol': data[1],
            'effect': data[2]
        }

        try:
            interactions[key_gene].append(interaction)
        except KeyError:
            interactions[key_gene] = [interaction]

    return interactions

def trrust_human():
    return trrust_general('human')

def trrust_mouse():
    return trrust_general('mouse')

def scrape_human():
    return trrust_scraping('human')

def scrape_mouse():
    return trrust_scraping('mouse')
