
from future.utils import iteritems
from __future__ import print_function

import pypath.curl

tgzurl = 'http://www.ebi.ac.uk/~denes/54b510889336eb2591d8beff/test.tar.gz'
zipurl = 'http://www.ebi.ac.uk/~denes/54b510889336eb2591d8beff/test.zip'
txturl = 'http://www.ebi.ac.uk/~denes/54b510889336eb2591d8beff/testfile.txt'
gzurl = 'http://www.ebi.ac.uk/~denes/54b510889336eb2591d8beff/testfile.txt.gz'
isourl = 'http://www.ebi.ac.uk/~denes/54b510889336eb2591d8beff/testfile-iso.txt'

test_urls = {
    'tgz': tgzurl,
    'zip': zipurl,
    'txt': txturl,
    'gz': gzurl,
    'iso': isourl
}

ls = [False, True]
cs = [False, True]

for typ, url in iteritems(test_urls):
    for large in ls:
        for cache in cs:
            print('')
            c = pypath.curl.Curl(url, large = large, cache = cache, silent = False)
            print('cache %s' % cache)
            print('large %s' % large)
            print(c.type)
            print(type(c.result))
            print(c.result)
            if type(c.result) is dict:
                for fname, content in iteritems(c.result):
                    print('\t%s' % fname)
                    print('\t', content)
