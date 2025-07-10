from collections.abc import Generator
import collections
from pypath.share import curl
from pypath.resources import urls


TcdbSubstrate = collections.namedtuple(
    'TcdbSubstrate',
    ['transporter_uniprot', 'substrate_id', 'substrate_name'],
)


def tcdb_substrate() -> Generator[TcdbSubstrate, None, None]:

    url_substrates = urls.urls['tcdb']['url_substrates']
    c = curl.Curl(url_substrates, large = True)
    TC2Sub = c.result

    url_acc2tc = urls.urls['tcdb']['url_acc2tc']
    c = curl.Curl(url_acc2tc, large = True)
    TC2Uni = c.result

    #ID conversion
    dic_TCDB_Uni = collections.defaultdict(list)
    for line in TC2Uni:

        uniprot, tcdb = line.strip('\n').split('\t')

        if uniprot:
            dic_TCDB_Uni[tcdb].append(uniprot)

    for line in TC2Sub:

        tcdb, substrates = line.strip('\n').split('\t')
        for substrate in substrates.split('|'):
            substrate_id, substrate_name = substrate.split(';')
            for transporter_uniprot in dic_TCDB_Uni[tcdb]:
                yield TcdbSubstrate(transporter_uniprot, substrate_id,substrate_name)


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



