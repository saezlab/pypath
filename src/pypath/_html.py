#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
# this module makes possible a
# dynamic data integration, download 
# files from various resources, in standard
# or non-standard text based and xml formats,
# process them, sometimes parse html
#

import bs4
import os

_fonts = open('fonts.css', 'r').read().decode('utf-8') \
    if os.path.exists('fonts.css') else ''

_ebi_logo = u'''data:image/png;base64,
    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'''.replace('\n', '')\
    .replace(' ', '')

_valid_html5 = u'''data:image/png;base64,
    iVBORw0KGgoAAAANSUhEUgAAAFAAAAAPCAIAAAD8q9/YAAAA6ElEQVRIieVXUQqEIBAdl74D57cD
    BB6mY3iojuFhgg7Qrx+eYD/eIjJF2rYVrI9B5KXjPGdUUlprqgkNEXnvnw7jJjBzg94y9OfddW4O
    1pz3cxHacSKi19Nh3I3qBXduXvc7N8PESMEUAqUlOjDBbw4WUwSfxXaGhZJl6GHpFgjmDHDy0/MP
    GftigjWwQ2s1heMgbH23fXfbZdMSrGnHCe2+k6sEQxgSm/JrpgSIMpvAkq/Yl/KlpWBUaWxBxrqN
    TCzmn7xnhUhTGjfraIaV1tp7X8k7zMzVP0t/j09JPx3GTWBmVdvf0hvVzH8tNHhSIAAAAABJRU5E
    rkJggg=='''.replace('\n', '').replace(' ', '')

_valid_css3 = u'''data:image/png;base64,
    iVBORw0KGgoAAAANSUhEUgAAAFAAAAAPCAIAAAD8q9/YAAAA7ElEQVRIieWX0Q2EIAyGC2EBuoAD
    MIcrOIejOIcrOIcDuAAjmHvopWlaLycS9YHvwdQ/hPSnrUQXY4SWCACQc347jYdAxECRH5f67fap
    T/NWv89NrEMHAP7tNJ6mOcNBve9Tz+3N8T71pMjOJ7FmFqjH5BSQIkWl2AWlaMOEtA3ClTwCFZSy
    Dh17UIFdw8iDuOb52LDlVz3rv3acd5o3VXOrwFFfFHF2hv24+HHh3masUgp3KQCkeWOT55UidIXJ
    FT9JtDPM5i9X2Fbv7wAfKqW4GGPOuZF7GBGbu5aaM/xt6bfTeAhEdK39LX0ALc+PofdKE1EAAAAA
    SUVORK5CYII='''.replace('\n', '').replace(' ', '')

_favicon = u'''data:image/png;base64,
    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'''\
    .replace('\n', '').replace(' ', '')

_default_css = u'''
                %s
                .code {
                    font-family: monospace;
                }
                .box {
                    border-style: solid;
                    border-width: 1px 1px 1px 3px;
                    margin: 0.75em 1em;
                    padding: 0px 0.75em;
                    color: #444444;
                }
                .quotebox {
                    border-color: #6FA6A9;
                    background-color: #E4ECEC;
                }
                .codebox {
                    border-color: #A9C98B;
                    background-color: #EFF4E9;
                }
                .small {
                    font-size: small;
                }
                body {
                    margin: 6ex 23.61%% 18ex 14.59%%;
                    width: 61.8%%;
                    color: #646567;
                    font-family: "HelveticaNeueLT Pro", "Helvetica Neue", Helvetica, Arial, sans-serif;
                }
                a:link, a:active, a:visited {
                    color: #646567;
                    text-decoration: underline;
                    text-decoration-style: dotted;
                }
                h1 {
                    color: #6FA6A9;
                }
                h2.omnipath, a.omnipath:link, a.omnipath:active, a.omnipath:visited {
                    color: #6EA945;
                    text-decoration: none;
                }
                h2.base, a.base:link, a.base:active, a.base:visited {
                    color: #646567;
                    text-decoration: none;
                }''' % _fonts

_header = u'''<!DOCTYPE html>
    <html lang="en">
        <head>
            <meta http-equiv="X-UA-Compatible" content="IE=edge">
            <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
            <meta name="dc.language" content="en">
            <meta name="viewport" content="width=device-width, 
            initial-scale=1.0">
            <link rel="icon" type="image/png" href="%s" />
            <style type="text/css">%s
            </style>
            <title>%s</title>
        </head>
        <body>
        '''

_footer = u'''
        <br>
        <br>
        <p> <a href="http://www.ebi.ac.uk/~denes">Denes Turei, 2016.</a> 
            Feedback: denes@ebi.ac.uk </p>
        <p>
            <a href="https://validator.w3.org/check/referer">
                <img 
                    src="%s"
                     alt="Valid HTML5" />
            </a>
            <a href="https://jigsaw.w3.org/css-validator/check/referer">
                <img 
                    src="%s"
                     alt="Valid CSS3" />
            </a>
        </p>
        \t</body>\n
    </html>\n''' % (_valid_html5, _valid_css3)

def get_header(title = ''):
    return _header % (_favicon, _default_css, 'OmniPath :: %s'%title)

def default_template(content, page_title, title = ''):
    page = u'%s\n\t<img src="%s" alt="EMBL-EBI"/>\n\t'\
            '<h1>%s</h1>\n\t%s\n%s' % \
            (get_header(title), _ebi_logo, page_title, content, _footer)
    return bs4.BeautifulSoup(page, 'lxml', from_encoding = 'utf-8')\
        .prettify().encode('utf-8')

def main_page():
    doc = u'<p>Welcome to the home of OmniPath, a comprehensive collection of '\
        'literature curated human signaling pathways. And pypath, the powerful '\
        'Python module for molecular networks and pathways analysis.</p>\n'\
        '<h2><a class="omnipath" href="/info">Metainformation about signaling pathway resources'\
        '</a></h2>\n'\
        '<h2><a class="omnipath" href="http://github.com/saezlab/pypath" target="_blank">pypath'\
        ' code</a></h2>\n'\
        '<h2><a class="omnipath" href="http://pypath.omnipathdb.org/" target = "_blank">pypath'\
        ' documentation</a></h2>\n'\
        '<h2><a class="omnipath" href="" target="_blank">The article</a></h2>\n'\
        '<p>D Turei, T Korcsmaros and J Saez-Rodriguez: Benchmark of literature'\
        ' curated signaling pathway resources (submitted February 2016)</p>\n'\
        '<h2><a class="omnipath" href="https://github.com/saezlab/pypath/blob/master/webservice.rst" '\
        'target="_blank">How to access the data?</a></h2>\n'
    return default_template(doc, 
        'OmniPath: literature curated human signaling pathways', 
        'literature curated human signaling pathways')
