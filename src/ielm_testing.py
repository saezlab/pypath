from pypath import dataio
from pypath import data_formats
import bs4

lst = []
with open('list', 'r') as f:
    for l in f:
        l = l.strip().split('\t')
        lst.append(l)

dmi = dataio.get_ielm(lst)

mydomains = 'HMMS'
id_type = 'UniProtKB AC'

data = dataio.get_ielm(lst[:1])

url = data_formats.urls['proteomic_ielm']['url']
network = ''
for pp in lst:
    network += '%s %s\r\n' % (pp[0], pp[1])



post = {'network': network, 'databases': id_type, 'mydomains': mydomains}
data = dataio.curl(url, post = post, silent = False, cache = False, headers = headers)




post = {'session_ID': sessid,'database': 'UniProtKB_ID', 'number': '', 'domains':' HMMS'}

data2 = dataio.curl('http://i.elm.eu.org/wait_2/', post = post, silent = False)

soup = bs4.BeautifulSoup(data)

wait = 0
while True:
    data3 = dataio.curl('http://i.elm.eu.org/proteomic_results/%s'%sessid, silent = False)
    soup3 = bs4.BeautifulSoup(data3)
    if len(soup3.find_all('table')) > 0:
        return soup3
    if wait > 600:
        return None
    time.sleep(3)
    

headers = ['Host: i.elm.eu.org', 
           'User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:30.0) Gecko/20110304 Firefox/30.0',
           'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8', 'Accept-Language: en-US,en;q=0.5',
           'DNT: 1',
           'Referer: http://i.elm.eu.org/wait_2/',
           'Connection: keep-alive']



curl 'http://i.elm.eu.org/test_submit/' -H 'Host: i.elm.eu.org' -H 'User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:30.0) Gecko/20110304 Firefox/30.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' -H 'Accept-Language: en-US,en;q=0.5' --compressed -H 'DNT: 1' -H 'Referer: http://i.elm.eu.org/search/' -H 'Connection: keep-alive' --data 'mydomains=HMMS&network=P31749+P42858%0D%0AP31749+P10070%0D%0AP31749+Q16512%0D%0AP31749+P31947%0D%0AP31749+P31946%0D%0AP31749+Q16584%0D%0AP31749+Q12778%0D%0AP31749+O15530%0D%0AP31749+Q6P1N0%0D%0AP31749+Q9UNE7%0D%0AP31749+Q92793%0D%0AP31749+P07900%0D%0AP31749+P19838%0D%0AP31749+P27348%0D%0AP31749+P46937%0D%0AP31749+P51617%0D%0AP31749+Q04917%0D%0AP31749+Q05513%0D%0AP31749+Q04912%0D%0AP31749+P38936%0D%0AP31749+P01375%0D%0AP31749+Q02078%0D%0AP31749+O15519%0D%0AP31749+Q06413%0D%0AP31749+P98170%0D%0AP31749+Q9UHD2%0D%0AP31749+P12755%0D%0AP31749+O43521%0D%0AP31749+O43524%0D%0AP31749+Q12824%0D%0AP31749+O14492%0D%0AP31749+P42336%0D%0AP31749+O14746%0D%0AP31749+P42345%0D%0AP31749+P41279%0D%0AP31749+Q13322%0D%0AP31749+P32121%0D%0AP31749+P21860%0D%0AP31749+O95999%0D%0AP31749+P08238%0D%0AP31749+P14778%0D%0AP31749+P62258%0D%0AP31749+Q8TAQ2%0D%0AP31749+Q9BZQ8%0D%0AP31749+P30291%0D%0AP31749+Q8NCD3%0D%0AP31749+Q53GL0%0D%0AP31749+P68400%0D%0AP31749+Q14653%0D%0AP31749+Q13362%0D%0AP31749+Q12888%0D%0AP31749+Q92547%0D%0AP31749+Q92574%0D%0AP31749+P67809%0D%0AP31749+O15360%0D%0AP31749+P53365%0D%0AP31749+Q16236%0D%0AP31749+Q9BVI0%0D%0AP31749+Q86V81%0D%0AP31749+Q15672%0D%0AP31749+P29317%0D%0AP31749+P07550%0D%0AP31749+Q9Y3C5%0D%0AP31749+P26358%0D%0AP31749+P40818%0D%0AP31749+P54274%0D%0AP31749+Q9H6Z4%0D%0AP31749+O15151%0D%0AP31749+Q15910%0D%0AP31749+P49815%0D%0AP31749+P19634%0D%0AP31749+Q7L5N1%0D%0AP31749+P14672%0D%0AP31749+O94875%0D%0AP31749+P25445%0D%0AP31749+Q15942%0D%0AP31749+O43464%0D%0AP31749+Q9UQF2%0D%0AP31749+Q99418%0D%0AP31749+P23396%0D%0AP31749+P17542%0D%0AP31749+O95793%0D%0AP31749+P09601%0D%0AP31749+Q16656%0D%0AP31749+Q7Z6J0%0D%0AP31749+O14558%0D%0AP31749+O60343%0D%0AP31749+P21453%0D%0AP31749+Q05195%0D%0AP31749+O60825%0D%0AP31749+P49918%0D%0AP31749+Q9UKV3%0D%0AP31749+P13631%0D%0AP31749+Q9NPC1%0D%0AP31749+Q07352%0D%0AP31749+Q9H4X1%0D%0AP31749+Q14814%0D%0AP31749+Q9BWT1%0D%0AP31749+P49840%0D%0AP31749+Q15027%0D%0A&databases=UniProtKB+AC' > curl1.out
'network=GRB2_HUMAN++++++++SOS1_HUMAN&databases=UniProtKB_ID&mydomains=HMMS' > curl1.out

curl 'http://i.elm.eu.org/wait_2/' -H 'Host: i.elm.eu.org' -H 'User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:30.0) Gecko/20110304 Firefox/30.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' -H 'Accept-Language: en-US,en;q=0.5' --compressed -H 'DNT: 1' -H 'Referer: http://i.elm.eu.org/test_submit/' -H 'Connection: keep-alive' --data 'session_IDa7e3642a8352a248dfe6ae35dc6d70ecc7f87468=&database=UniProtKB_ID&number=&domains=HMMS' > curl2.out

a7e3642a8352a248dfe6ae35dc6d70ecc7f87468

curl 'http://i.elm.eu.org/proteomic_results/a7e3642a8352a248dfe6ae35dc6d70ecc7f87468' -H 'Host: i.elm.eu.org' -H 'User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:30.0) Gecko/20110304 Firefox/30.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' -H 'Accept-Language: en-US,en;q=0.5' --compressed -H 'DNT: 1' -H 'Referer: http://i.elm.eu.org/test_submit/' -H 'Connection: keep-alive' > curl3.out