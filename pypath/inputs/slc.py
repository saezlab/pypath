
#SLC table
import requests
url = "https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209652&file=msb209652-sup-0003-TableEV2.xlsx"
c = curl.Curl(url)
SLC_table = c.result

url = "https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209652&file=msb209652-sup-0002-TableEV1.xlsx"
