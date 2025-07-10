from collections.abc import Generator
import collections
from pypath.share import curl
from pypath.resources import urls
from pypath.inputs import uniprot
import pandas as pd

TcdbSubstrate = collections.namedtuple(
    'TcdbSubstrate',
    ['transporter_uniprot', 'substrate_id', 'substrate_name','location','features'],
)

SLCSubstrate = collections.namedtuple(
    'SLCSubstrate',
    ['transporter_uniprot', 'substrate_id','substrate_name','location','features'],
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
        uniprot_symbol, tcdb = line.strip('\n').split('\t')
        if uniprot_symbol:
            dic_TCDB_Uni[tcdb].append(uniprot_symbol)
    for line in TC2Sub:
        tcdb, substrates = line.strip('\n').split('\t')
        for substrate in substrates.split('|'):
            substrate_id, substrate_name = substrate.split(';')
            for transporter_uniprot in dic_TCDB_Uni[tcdb]:
                try:
                    uploc = uniprot.uniprot_locations(accession=transporter_uniprot.upper(), organism=None)
                    if any(uploc):
                        loc = list(uploc[transporter_uniprot.upper()])[0]
                        location = loc.location
                        features = loc.features
                    else:
                        location = None
                        features = None
                except ValueError:
                    location = None
                    features = None
                yield TcdbSubstrate(transporter_uniprot, substrate_id, substrate_name,location,features)

def slc_table() -> Generator[SLCSubstrate, None, None]:
    url = "https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209652&file=msb209652-sup-0003-TableEV2.xlsx"
    df = pd.read_excel(url)
    substrate_annotation = pd.read_excel("/Users/priscillabai/Downloads/msb209652-sup-0002-tableev1.xlsx", sheet_name=2)
    slc_annotation = pd.read_excel("/Users/priscillabai/Downloads/msb209652-sup-0002-tableev1.xlsx", sheet_name=1)
    slc_annotation = slc_annotation[slc_annotation['Substrate_class'] != 'Orphan']
    for transporter_id, substrate_name,location in zip(
            slc_annotation["Ensembl_ID"],
            slc_annotation["Substrates"],
            slc_annotation["Subcellular_localization"]):
        #ID conversion
        ID_conversion(df_new, "Act", "Activator", Name_list[0])


        yield SLCSubstrate(transporter_id, substrate_name, location,)






