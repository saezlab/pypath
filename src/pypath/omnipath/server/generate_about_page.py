#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module.
#  Contains descriptions for all resources: e.g. their URLs
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import sys
from future.utils import iteritems

import pypath._html as _html
import pypath.urls as urls

import pypath.share.session as session_mod
_logger = session_mod.Logger(name='generate_about_page')

__all__ = ['generate_about_html', 'write_html']

if 'long' not in __builtins__:
    long = int

if 'unicode' not in __builtins__:
    unicode = str

pypath_methods = {
    'data': 'Data source (URLs and files)',
    'format': 'Data format definition',
    'intr': 'Interactions',
    'input': 'Data input methods',
    'ptm': 'Enzyme-substrate relationships and PTMs',
    'dmi': 'Domain-motif interactions',
    'ddi': 'Domain-domain interactions',
    'misc': 'Miscellaneous'
}


def generate_about_html(descriptions, format='b'):
    """
    Generates a HTML page from the `descriptions` array.
    This HTML is provided by the webservice under `/info`,
    or can be saved locally with `write_html()`.
    """
    # Header
    title = 'Metadata about signaling pathway resources'
    doc = (
        '<div class="yellowbox box">\n'
        '<p>\n'
            '<em>\n'
            'Information on this page has been last revised in Nov 2016.\n'
            'As of Oct 2019 we are working on updating and extending this\n'
            'page and will publish the new version soon.\n'
            'About updates of the OmniPath database content please refer to\n'
            '<a href="http://archive.omnipathdb.org/README.txt">\n'
                'our archive.\n'
            '</a>\n'
            '</em>\n'
        '</p>\n'
        '</div>\n'
        '<p>This collection was created during the construction '
        'of OmniPath when we considered more than 50 resources and '
        'selected the ones containing literature curation effort. '
        'OmniPath is a network '
        'of signaling pathways intending to '
        'combine all high quality, manually curated efforts. The '
        'descriptions here cite the relevant sentences '
        'about the curation protocols from the original articles and webpages. '
        'URLs pointing to the articles and the webpages, and some '
        'additional metadata are provided where available. '
        'The resources with green title are included by default in '
        'OmniPath. <span class="code">pypath</span> methods are listed '
        ' where available, to know more please look at <a '
        'target="_blank" href="http://pypath.omnipathdb.org/">'
        'pypath documentation.</a> This list is only about network '
        'resources. <span class="code">pypath</span> is able to '
        'process and integrate many other resources, please see '
        'the paper and the documentation to know more.</p>'
        '<p class="small">We searched for license information '
        'in the main, About, Download and FAQ sections of the webpages, '
        'and run Google searches for the database name and license. '
        'Where we could not find anything about licensing, we assumed '
        'no license. Unfortunately due to todays restrictive copyright '
        'legislations, users don\'t have the freedom to use, modify and '
        'redistribute the data without a license explicitely granting '
        'these to them. Despite the clear intention from the authors to '
        'make their data public, and statements on the webpage like '
        '"free to use" or "available for download".</p>\n'
    )
    doc += '\t<h2>Contents</h2>\n'
    doc += '\t<ul>\n'

    # Table of Content
    for k, v in sorted(descriptions.items(), key=lambda x: x[0].lower()):
        doc += '\t\t\t<li><a href="#%s" class="%s">%s</a></li>\n' % \
            (k, 'omnipath' if 'omnipath' in v and v['omnipath'] else 'base',
             v['label'] if 'label' in v else k)
    doc += '\t</ul>\n'

    # Sections
    for k, v in sorted(descriptions.items(), key=lambda x: x[0].lower()):

        doc += u'\t\t<br>\n\t\t<h2 id="%s" class="%s">%s%s</h2>\n' % \
            (k, 'omnipath' if 'omnipath' in v and v['omnipath'] else 'base',
             v['label'] if 'label' in v else k,
             (u' – %s' % (v['full_name'],)) if 'full_name' in v else '')
        doc += '\t\t\t<p><b>Category || Subcategory &gt;&gt;&gt;</b> %s || %s</p>\n' % \
            (v['type'].capitalize() if 'type' in v else 'Undefined',
                v['subtype'].capitalize() if 'subtype' in v else 'Undefined')
        if 'year' in v:
            doc += '\t\t\t<h3>Last released: %u<\h3>\n' % v['year']
        if 'releases' in v:
            doc += '\t\t\t<p><b>Released in years: </b>%s</p>\n' % \
                ', '.join(['%u' % y for y in v['releases']])
        if 'authors' in v and v['authors'] is not None:
            doc += '\t\t\t<p><b>Created by </b>%s</p>\n' % ', '.join(v[
                'authors'])
        if 'emails' in v and v['emails'] is not None:
            doc += '\t\t\t<p><b>Contact: </b></p>\n\n\t\t\t\t<ul>\n%s\n' % \
                ''.join(['\t\t\t\t<li><a href="mailto:%s">%s &lt;%s&gt;</li>\n' %
                         (email, contact_name, email) for email, contact_name in
                         zip(v['emails'][::2], v['emails'][1::2])])
            doc += '\t\t\t\t</ul>\n'

            emails = (
                ", ".join("%s <%s>" % (contact_name, email)
                          for email, contact_name in
                          zip(v['emails'][::2], v['emails'][1::2]))  # zip makes (0,1),(2,3) ...
                if 'emails' in v else ''
            )
        if 'license' in v:
            try:
                doc += '\t\t\t<p><b>License:</b> %s%s%s\n' % (
                    ('<a href="%s" target="_blank">' % v['license']['url']) if
                    'url' in v['license'] else '', v['license']['name'], '</a>'
                    if 'url' in v['license'] else '')
                doc += '\t\t\t%s\n</p>' % "".join(v['license']['comment'])\
                        if 'comment' in v['license'] else ''

            except KeyError:
                sys.stdout.write('Wrong license format or incomplete information for %s\n' % k)
                sys.stdout.flush()
        if 'urls' in v:
            for uk, uv in iteritems(v['urls']):
                if len(uv) > 0 and uk != 'omictools':
                    try:
                        doc += '\t\t\t<h3>%s</h3>\n' % (uk.capitalize())
                        doc += '\t\t\t<ul>\n'
                        for a in uv:
                            doc += (
                                    '\t\t\t\t<li><a href="%s" '
                                    'target="_blank">%s</a></li>\n' % (
                                        a, a
                                    )
                            )
                        doc += '\t\t\t</ul>\n'
                    except UnicodeDecodeError:
                        sys.stdout.write('UnicodeDecodeError at %s\n' % k)
                        sys.stdout.flush()
        if 'pubmeds' in v:
            doc += '\t\t\t<h3>PubMed</h3>\n'
            doc += '\t\t\t<ul>\n'
            for pmid in v['pubmeds']:
                doc += '\t\t\t\t<li><a href="%s" '\
                    'target="_blank">%s</a></li>\n' % (
                        'http://www.ncbi.nlm.nih.gov/pubmed/%u' % pmid,
                        'http://www.ncbi.nlm.nih.gov/pubmed/%u' % pmid
                    )
            doc += '\t\t\t</ul>\n'
        if 'urls' in v and 'omictools' in v['urls'] or 'pathguide' in v:
            doc += '\t\t\t<h3>Collections</h3>\n\t\t\t<ul>'
            if 'omictools' in v['urls']:
                doc += '\t\t\t<li><a href="%s" target="_blank">OmicTools</a></li>\n' % \
                    v['urls']['omictools'][0]
            if 'pathguide' in v:
                doc += '\t\t\t<li><a href="%s" target="_blank">PathGuide</a></li>\n' % \
                    (urls.urls['pathguide']['url'] % v['pathguide'])
            doc += '\t\t\t</ul>\n'
        if 'taxons' in v:
            doc += '<p><b>Taxons: </b><em>%s</em></p>' % \
                ', '.join(['%s%s' % (t[0].upper(), t[1:])
                           for t in v['taxons']])
        if 'size' in v and type(v['size']) is dict and \
                v['size']['nodes'] is not None and v['size']['edges'] is not None:
            doc += '<p><b>Nodes: </b>%s, <b>Edges:</b>%s</p>' % (
                v['size']['nodes'], v['size']['edges'])
        if 'data_import' in v:
            doc += '\t\t\t<p><b>Direct data import from: </b>%s</p>\n' % \
                ', '.join(v['data_import'])
        if 'includes' in v:
            doc += '\t\t\t<p><b>Includes data from: </b>%s</p>\n' % \
                ', '.join(v['includes'])
        if 'descriptions' in v or 'notes' in v:
            doc += '\t\t\t<h3>Quotes</h3>\n'
            if 'descriptions' in v:
                doc += '\t\t\t\t<div class="quotebox box">\n'
                description_full = " ".join(v['descriptions'])
                pars = description_full.split('\n')
                if len(pars) > 0:
                    doc += '\t\t\t\t<p>%s</p>\n' % pars[0]
                # for p in pars:
                #     p = p.strip()
                #     if len(p) > 0:
                #         doc += '\t\t\t\t<p>%s</p>\n' % p
                doc += '\t\t\t\t</div>\n'
            if 'notes' in v:
                doc += '\t\t\t<h3>Notes</h3>\n'
                doc += '\t\t\t\t<div class="quotebox box">\n'
                notes_full = " ".join(v['notes'])
                pars = notes_full.split('\n')
                for p in pars:
                    p = p.strip()
                    if len(p) > 0:
                        doc += '\t\t\t\t<p>%s</p>\n' % p
                doc += '\t\t\t\t</div>\n'
        if 'data_integration' in v:
            doc += '\t\t\t<p><b>Data integration in '\
                '<span class="code">pypath:</span></b> %s</p>' % \
                v['data_integration']
        if 'pypath' in v:
            doc += '\t\t\t<h3>Methods in <span class="code">pypath'\
                '</span></h3>\n'
            doc += '\t\t\t\t<div class="codebox box">\n'
            for cat in sorted(pypath_methods.keys()):
                name = pypath_methods[cat]
                if cat in v['pypath']:
                    doc += '\t\t\t\t\t<p>%s</p>\n\t\t\t\t\t<ul>\n' % name
                    for met in v['pypath'][cat]:
                        doc += '\t\t\t\t\t\t<li><span class="code">%s'\
                            '</span></li>\n' % met
                    doc += '\t\t\t\t\t</ul>\n'
            doc += '\t\t\t\t</div>\n'

    if format == 'b':
        return _html.default_template(doc, title, title)
    return _html.default_template(doc, title, title, format='str')


def write_html(descriptions, filename='resources.html'):
    """
    Saves the HTML descriptions to custom local file.
    """
    if not descriptions:
        _logger._console("Empty description was provided. Cannot write html.")
        return None

    html = generate_about_html(descriptions, format="str")
    #with codecs.open(filename, encoding='utf-8', mode='w') as f:
    with open(filename, 'w') as f:
        f.write(html)


def resource_list_latex(descriptions,
                        filename='resource-list.tex',
                        latex_hdr=True,
                        fontsize=8,
                        font='HelveticaNeueLTStd-LtCn'):
    """
    Generates Supplementary Table 3 (The list of the 52 resources considered) for the article.
    """
    if not descriptions:
        _logger._console("Empty description was provided. Cannot generate table.")
        return None

    _latex_hdr = r'''\documentclass[a4paper,%upt]{extarticle}
        \usepackage{fontspec}
        \usepackage{xunicode}
        \usepackage{polyglossia}
        \setdefaultlanguage{english}
        \usepackage{xltxtra}
        \usepackage{microtype}
        \usepackage[margin=5pt,portrait,paperwidth=15cm,paperheight=18cm]{geometry}
        \usepackage{amsmath}
        \usepackage{amssymb}
        \usepackage{textcomp}
        \usepackage[usenames,dvipsnames,svgnames,table]{xcolor}
        \usepackage{color}
        \usepackage{booktabs}
        \usepackage{tabularx}
        \setmainfont{%s}
        \definecolor{grey875}{gray}{0.125}
        \begin{document}
        \color{grey875}
        \thispagestyle{empty}
        \vfill
    ''' % (fontsize, font) if latex_hdr else ''
    _latex_end = r'''
            \end{document}
        ''' if latex_hdr else ''
    tex = r'''\begin{tabularx}{0.94\textwidth}{>{\raggedright\scriptsize\arraybackslash\hsize=.15\hsize}X>{\raggedright\scriptsize\arraybackslash\hsize=.35\hsize}X>{\raggedright\scriptsize\arraybackslash\hsize=.15\hsize}X>{\raggedright\scriptsize\arraybackslash\hsize=.15\hsize}X}
    \toprule
    Resource name & Class, subclass & Resource name & Class, subclass \\
    \midrule
    '''
    res = sorted(
        [(v['label'] if 'label' in v else k,
          '%s, %s' % (v['type'].capitalize(), v['subtype'].capitalize())
          if 'type' in v and 'subtype' in v else '')
         for k, v in iteritems(descriptions)],
        key=lambda x: x[0].lower())
    if len(res) % 2 != 0:
        res.append('')
    res2 = zip(res[:int(len(res) / 2)], res[int(len(res) / 2):])
    for r in res2:
        tex += r'%s & %s & %s & %s \\' % (
            r[0][0].replace('&', '\&'),
            r[0][1].replace('&', '\&'),
            (r[1][0] if len(r[1]) else '').replace('&', '\&'),
            (r[1][1] if len(r[1]) else '').replace('&', '\&'),
        ) + '\n'
    tex += r'\bottomrule' + '\n'
    tex += r'\end{tabularx}' + '\n'
    with open(filename, 'w') as f:
        f.write('%s%s%s' % (_latex_hdr if latex_hdr else '', tex, _latex_end
                            if latex_hdr else ''))


def export_licenses(descriptions, outfile='licenses.tsv'):

    if not descriptions:
        _logger._console("Empty description was provided. Cannot export licenses.")
        return None

    hdr = [
        'Name',
        'License',
        'License URL',
        'Contact',
    ]
    rows = []

    for k, v in iteritems(descriptions):
        name = v['label'] if 'label' in v else k
        license_name = v['license']['name'] if 'license' in v and 'name' in v['license'] else ''
        license_comment = (
            "".join(v['license']['comment'])
            if 'license' in v and 'comment' in v['license'] else ''
        )
        license_url = (
            v['license']['url']
                if 'license' in v and 'url' in v['license'] else
            ''
        )
        emails = (
            ", ".join("%s <%s>" % (contact_name, email)
                for email, contact_name in zip(v['emails'][::2], v['emails'][1::2]))  # zip makes (0,1),(2,3) ...
            if 'emails' in v else ''
        )

        rows.append([
            name,
            license_name,
            license_url,
            emails,
            license_comment
        ])

    with open(outfile, 'w') as fp:

        _ = fp.write('\t'.join(hdr) + '\n')
        [print(row) for row in rows]
        _ = fp.write('\n'.join(
            '\t'.join(row) for row in rows
        ))


