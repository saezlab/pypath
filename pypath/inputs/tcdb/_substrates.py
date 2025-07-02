from collections.abc import Generator
import re
import pathlib as pl
import collections
import pandas as pd
import pyparsing
from pypath_common import _misc as _common
from pypath.share import curl
from pypath.utils import mapping, taxonomy
from pypath.resources import urls


TcsbSubstrate = collections.namedtuple(
    'TcsbSubstrate',
    ['transporter', 'substrate'],
)


def tcdb_substrate() -> Generator[tuple[str, str], None, None]:

    url_substrates = urls.urls['tcdb']['url_substrates']
    c = curl.Curl(url_substrates, large = True)
    TC2Sub = c.result

    url_acc2tc = urls.urls['tcdb']['url_acc2tc']
    c = curl.Curl(url_acc2tc, large = True)
    TC2Uni = c.result

    #ID conversion
    dic_TCDB_Uni = collections.defaultdict(list)
    for line in TC2Uni:

        uniprot, tcdb = line.split('\t')
        dic_TCDB_Uni[tcdb].append(uniprot)

    for line in TC2Sub:

        tcdb, substrate = line.split('\t')

        for uniprot in dic_TCDB_Uni[tcdb]:

            yield TcsbSubstrate(uniprot, substrate)


def protein_location():
    # protein location
    import requests
    url = "https://rest.uniprot.org/uniprotkb/search?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Cft_intramem%2Ccc_subcellular_location%2Cft_transmem%2Cft_topo_dom&format=tsv&query=%28%28taxonomy_id%3A10090%29%29&size=500"
    from requests.adapters import HTTPAdapter, Retry
    import re
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)

    interactions = []
    for batch, total in get_batch(url):
        for line in batch.text.splitlines()[1:]:
            interactions.append(line)
        print(f'{len(interactions)} / {total}')

    fileObject = open('data/uniprot_subLocation.txt', 'w')
    for i in interactions:
        fileObject.write(i)
        fileObject.write("\n")
    fileObject.close()



