import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls

def guide2pharma_table():
    """"""
    url = urls.urls['gtp']['url']

    c = curl.Curl(url, silent = False, large = True, encoding = 'utf-8')

    return c

    # data = csv.DictReader(c.result)    

