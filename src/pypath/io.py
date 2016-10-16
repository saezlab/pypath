#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pycurl
from StringIO import StringIO


def curl(url):
    result = StringIO()
    header = StringIO()
    c = pycurl.Curl()
    c.setopt(c.URL, url)
    c.setopt(c.FOLLOWLOCATION, True)
    c.setopt(c.WRITEFUNCTION, result.write)
    c.setopt(c.HEADERFUNCTION, header.write)
    for i in range(3):
        try:
            c.perform()
            status = c.getinfo(pycurl.HTTP_CODE)
            break
        except:
            status = 500
    c.close()
    if status == 200:
        return result.getvalue()
    else:
        return None
