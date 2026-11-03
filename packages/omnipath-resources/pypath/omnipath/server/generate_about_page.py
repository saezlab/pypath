#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

from future.utils import iteritems

import pypath.omnipath.server._html as _html
import pypath.resources.urls as urls

import pypath.share.common as common
import pypath.share.session as session_mod
_logger = session_mod.Logger(name = 'generate_about_page')
_log = _logger._log

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
    'misc': 'Miscellaneous',
}


def generate_about_html(descriptions, format = 'b'):
    """
    Generates a HTML page from the `descriptions` array.
    This HTML is provided by the webservice under `/info`,
    or can be saved locally with `write_html()`.
    """
    # Header
    title = 'Metadata about resources in OmniPath and pypath'
    doc = (
        '<div class="yellowbox box">\n'
        '<p>\n'
            '<em>\n'
            'We created this webpage in Nov 2016.\n'
            'Since the end of 2019 we have been gradually updating and \n'
            'extending the information. As of July 2020 we updated the URLs\n'
            'and licensing terms of most of the resources.\n'
            'Most of the "pypath methods" and the\n'
            'years of releases are out of date.\n'
            'We will keep updating these, if you find any wrong information\n'
            'please notify us at omnipathdb@gmail.com.\n'
            'About updates of the OmniPath database content please refer to\n'
            '<a href="http://archive.omnipathdb.org/README.txt">\n'
                'our archive.\n'
            '</a>\n'
            '</em>\n'
        '</p>\n'
        '</div>\n'
        '<p>This collection is a byproduct of the development of\n'
        'OmniPath, a database built from above 100 resources.\n'
        'Initially OmniPath focused on the literature curatied activity\n'
        'flow networks. Today it covers a much broader range of molecular\n'
        'interaction data, and besides its network database OmniPath\n'
        'has four other databases: enzyme-PTM relationships, protein\n'
        'complexes, molecular annotations (function, localization,\n'
        'structure, etc) and intercellular communication roles.\n'
        'The "omnipath" dataset of the network database follows the\n'
        'principles of the initial release of OmniPath, focusing on \n'
        'high quality, manually curated signaling pathways. The '
        'descriptions here cite the relevant sentences '
        'about the curation protocols from the original articles and\n'
        'webpages. URLs pointing to the articles and the webpages, and some\n'
        'additional metadata are provided where available.\n'
        'The resources with green title are included by default in\n'
        'OmniPath. <span class="code">pypath</span> methods are listed\n'
        ' where available, to know more please look at\n<a '
        'target="_blank" href="https://pypath.omnipathdb.org/">'
        'pypath documentation.</a>\n'
        '<p class="small">\n'
        '<strong>How we collected the license information?</strong>\n'
        'We searched for license information\n'
        'in the main, About, Download and FAQ sections of the webpages,\n'
        'and run Google searches for the database name and license.\n'
        'Where we could not find anything about licensing, we assumed\n'
        'no license. Unfortunately due to todays restrictive copyright\n'
        'legislations, users don\'t have the freedom to use, modify and\n'
        'redistribute the data without a license explicitely granting\n'
        'these to them. Despite the clear intention from the authors to\n'
        'make their data public, and statements on the webpage like\n'
        '"free to use" or "available for download". In these cases\n'
        'we contacted the authors for permission to redistribute their\n'
        'data.</p>\n'
    )
    doc += '\t<h2>Contents</h2>\n'
    doc += '\t<ul>\n'

    # Table of Content
    for k, v in sorted(descriptions.items(), key=lambda x: x[0].lower()):
        doc += '\t\t\t<li><a href="#%s" class="%s">%s%s</a></li>\n' % \
            (
                k,
                'omnipath'
                    if 'omnipath' in v and v['omnipath'] else 'base',
                v['label'] if 'label' in v else k,
                (' – %s' % v['full_name']) if 'full_name' in v else ''
            )
    doc += '\t</ul>\n'

    # Sections
    for k, v in sorted(descriptions.items(), key=lambda x: x[0].lower()):

        resource_type = (
            common.to_list(v['type']) if 'type' in v else ['Undefined']
        )
        resource_subtype = (
            common.to_list(v['subtype']) if 'subtype' in v else ['Undefined']
        )
        resource_type = ', '.join(t.capitalize() for t in resource_type)
        resource_subtype = ', '.join(t.capitalize() for t in resource_subtype)

        doc += u'\t\t<br>\n\t\t<h2 id="%s" class="%s">%s%s</h2>\n' % \
            (k, 'omnipath' if 'omnipath' in v and v['omnipath'] else 'base',
             v['label'] if 'label' in v else k,
             (u' – %s' % (v['full_name'],)) if 'full_name' in v else '')
        doc += '\t\t\t<p><b>Category || Subcategory &gt;&gt;&gt;</b> %s || %s</p>\n' % \
            (resource_type, resource_subtype)
        if 'year' in v:
            doc += '\t\t\t<h3>Last released: %u<\\h3>\n' % v['year']
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
                license = v['license']

                doc += (
                    '\t\t\t<p><b>License:</b> '
                    '<a href="%s" target="_blank">%s%s</a>%s</p>\n' % (
                        license.url,
                        license.full_name,
                        (' (%s)' % license.name) if license.name != k else '',
                        ('</p>\n\t\t\t<p><em>"%s"</em>' % license.note)
                            if hasattr(license, 'note') else
                        ''
                    )
                )

            except KeyError:

                _log(
                    'Wrong license format or incomplete '
                    'information for `%s`.' % k
                )

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
                        _log('UnicodeDecodeError at %s' % k)

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
        # turned off temporarily, until we update this information
        #if 'pypath' in v:
        if False:
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
            r[0][0].replace('&', r'\&'),
            r[0][1].replace('&', r'\&'),
            (r[1][0] if len(r[1]) else '').replace('&', r'\&'),
            (r[1][1] if len(r[1]) else '').replace('&', r'\&'),
        ) + '\n'
    tex += r'\bottomrule' + '\n'
    tex += r'\end{tabularx}' + '\n'
    with open(filename, 'w') as f:
        f.write('%s%s%s' % (_latex_hdr if latex_hdr else '', tex, _latex_end
                            if latex_hdr else ''))


def export_licenses(descriptions, outfile='licenses.tsv'):

    if not descriptions:
        _logger._console(
            'Empty description was provided. Cannot export licenses.'
        )
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
        license_name = (
            v['license']['name']
                if 'license' in v and 'name' in v['license'] else
            ''
        )
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
                for email, contact_name
                # zip makes (0,1),(2,3) ...
                in zip(v['emails'][::2], v['emails'][1::2]))
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
