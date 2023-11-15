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
from past.builtins import xrange, range, reduce

from pypath.share import session as session_mod
_logger = session_mod.log()

import re
import sys
import os
import itertools
import collections
import importlib as imp
import copy
import subprocess
import warnings
from datetime import date

import math
import numpy as np
from numpy.random import randn
import matplotlib as mpl
import scipy.cluster.hierarchy as hc

try:
    import igraph

except ModuleNotFoundError:
    sys.stdout.write('Module `igraph` is not available.'
                     '\nSome plotting functionalities won\'t be accessible.\n')

try:
    import hcluster as hc2
except:
    _logger.msg('Module `hcluster` not available.')
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
from matplotlib import ticker
from scipy import stats

import pypath.share.common as common
import pypath.share.session as session
import pypath_common._constants as _const
from pypath.visual.igraph_drawing import DefaultGraphDrawerFFsupport
import pypath.resources.descriptions
import pypath.share.progress as progress
import pypath.internals.refs as _refs


def is_opentype_cff_font(filename):
    """
    This is necessary to fix a bug in matplotlib:
    https://github.com/matplotlib/matplotlib/pull/6714
    Returns True if the given font is a Postscript Compact Font Format
    Font embedded in an OpenType wrapper.  Used by the PostScript and
    PDF backends that can not subset these fonts.
    """
    if os.path.splitext(filename)[1].lower() == '.otf':
        result = _is_opentype_cff_font_cache.get(filename)
        if result is None:
            with open(filename, 'rb') as fd:
                tag = fd.read(4)
            result = (tag == b'OTTO')
            _is_opentype_cff_font_cache[filename] = result
        return result
    return False


mpl.font_manager.is_opentype_cff_font = is_opentype_cff_font

if not 'next' in __builtins__:

    def next(gen):
        return gen.next()

# helper functions


def rotate_labels(angles=(0, -90, -135, -180, -270, -315)):
    i = 0
    while True:
        yield angles[i % len(angles)]
        i += 1


def move_labels(dist=(0, 10, 20, 30, 40, 50, 60, 70)):
    i = 0
    while True:
        yield dist[i % len(dist)] if i < 20 else 20
        i += 1


def overlap(bbox1, bbox2):
    return (bbox1._points[0][0] > bbox2._points[0][0] and
            bbox1._points[0][0] < bbox2._points[1][0] or
            bbox1._points[1][0] > bbox2._points[0][0] and
            bbox1._points[1][0] < bbox2._points[1][0] or
            bbox2._points[0][0] > bbox1._points[0][0] and
            bbox2._points[0][0] < bbox1._points[1][0] or
            bbox2._points[1][0] > bbox1._points[0][0] and
            bbox2._points[1][0] < bbox1._points[1][0]) and \
        (bbox1._points[0][1] > bbox2._points[0][1] and
         bbox1._points[0][1] < bbox2._points[1][1] or
         bbox1._points[1][1] > bbox2._points[0][1] and
         bbox1._points[1][1] < bbox2._points[1][1] or
         bbox2._points[0][1] > bbox1._points[0][1] and
         bbox2._points[0][1] < bbox1._points[1][1] or
         bbox2._points[1][1] > bbox1._points[0][1] and
         bbox2._points[1][1] < bbox1._points[1][1])


def get_moves(bbox1, bbox2):
    xmove = 0
    ymove = 0
    if bbox1._points[0][0] > bbox2._points[0][0] and \
            bbox1._points[0][0] < bbox2._points[1][0] or \
            bbox1._points[1][0] > bbox2._points[0][0] and \
            bbox1._points[1][0] < bbox2._points[1][0] or \
            bbox2._points[0][0] > bbox1._points[0][0] and \
            bbox2._points[0][0] < bbox1._points[1][0] or \
            bbox2._points[1][0] > bbox1._points[0][0] and \
            bbox2._points[1][0] < bbox1._points[1][0]:
        if (bbox1._points[0][0] + bbox1._points[1][0]) / 2.0 < \
                (bbox2._points[0][0] + bbox2._points[1][0]) / 2.0:
            xmove = bbox1._points[1][0] - bbox2._points[0][0]
        else:
            xmove = bbox1._points[0][0] - bbox2._points[1][0]
    if bbox1._points[0][1] > bbox2._points[0][1] and \
            bbox1._points[0][1] < bbox2._points[1][1] or \
            bbox1._points[1][1] > bbox2._points[0][1] and \
            bbox1._points[1][1] < bbox2._points[1][1] or \
            bbox2._points[0][1] > bbox1._points[0][1] and \
            bbox2._points[0][1] < bbox1._points[1][1] or \
            bbox2._points[1][1] > bbox1._points[0][1] and \
            bbox2._points[1][1] < bbox1._points[1][1]:
        if (bbox1._points[0][1] + bbox1._points[1][1]) / 2.0 < \
                (bbox2._points[0][1] + bbox2._points[1][1]) / 2.0:
            ymove = bbox1._points[1][1] - bbox2._points[0][1]
        else:
            ymove = bbox1._points[0][1] - bbox2._points[1][1]
    return (xmove, ymove)


class Plot(object):

    def __init__(self,
                 fname=None,
                 font_family='Helvetica Neue LT Std',
                 font_style='normal',
                 font_weight='normal',
                 font_variant='normal',
                 font_stretch='normal',
                 palette=None,
                 context='poster',
                 lab_size=(9, 9),
                 axis_lab_size=10.0,
                 rc={}):
        for k, v in iteritems(locals()):
            if not hasattr(self, k) or getattr(self, k) is None:
                setattr(self, k, v)
        if type(self.lab_size) is not tuple:
            self.lab_size = (self.lab_size, ) * 2
        if 'axes.labelsize' not in self.rc:
            self.rc['axes.labelsize'] = self.axis_lab_size
        if 'ytick.labelsize' not in self.rc:
            self.rc['ytick.labelsize'] = self.lab_size[0]
        if 'ytick.labelsize' not in self.rc:
            self.rc['ytick.labelsize'] = self.lab_size[1]
        self.rc['font.family'] = self.font_family
        self.rc['font.style'] = self.font_style
        self.rc['font.variant'] = self.font_variant
        self.rc['font.weight'] = self.font_weight
        self.rc['font.stretch'] = self.font_stretch
        self.palette = palette or self.embl_palette()
        self.fp = mpl.font_manager.FontProperties(
            family=self.font_family,
            style=self.font_style,
            variant=self.font_variant,
            weight=self.font_weight,
            stretch=self.font_stretch)

    def embl_palette(self, inFile='embl_colors'):
        cols = []
        inFile = os.path.join(session.session().module_root, 'data', inFile)
        with open(inFile, 'r') as f:
            series = []
            for i, l in enumerate(f):
                l = [x.strip() for x in l.split(',')]
                series.append(
                    self.rgb2hex(tuple([256 * float(x) for x in l[0:3]])))
                if len(series) == 7:
                    cols.append(series)
                    series = []
        return cols


    @staticmethod
    def rgb2hex(rgb):

        rgb = tuple(int(i) for i in rgb)
        return '#%02x%02x%02x' % rgb


    def finish(self):
        '''
        Saves and closes a figure.
        '''

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        self.fig.savefig(self.fname)
        plt.close(self.fig)


class MultiBarplot(Plot):

    def __init__(self,
                 x,
                 y,
                 categories=None,
                 cat_names=None,
                 cat_ordr=None,
                 fname=None,
                 figsize=(12, 4),
                 xlab='',
                 ylab='',
                 title='',
                 lab_angle=90,
                 color='#007b7f',
                 order=False,
                 desc=True,
                 ylog=False,
                 legend=None,
                 fin=True,
                 rc={},
                 axis_lab_font={},
                 bar_args={},
                 ticklabel_font={},
                 legend_font={},
                 title_font={},
                 title_halign='center',
                 title_valign='top',
                 y2=None,
                 color2=None,
                 ylim=None,
                 grouped=False,
                 group_labels=[],
                 summary=False,
                 summary_name='',
                 uniform_ylim=False,
                 do=True,
                 legloc=0,
                 maketitle=True,
                 **kwargs):
        """
        Plots multiple barplots side-by-side.
        Not all options are compatible with each other.


        """

        for k, v in iteritems(locals()):
            setattr(self, k, v)

        for k, v in iteritems(kwargs):
            setattr(self, k, v)

        super(MultiBarplot, self).__init__()

        self.axes = []

        self.bar_args_default = {
            'width': 0.8,
            'edgecolor': 'none',
            'linewidth': 0.0,
            'align': 'center'
        }
        self.axis_lab_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'x-large'
        }
        self.ticklabel_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'large'
        }
        self.legend_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'small'
        }
        self.title_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'xx-large'
        }

        self.bar_args = common.merge_dicts(bar_args, self.bar_args_default)
        self.axis_lab_font = common.merge_dicts(axis_lab_font,
                                                self.axis_lab_font_default)
        self.ticklabel_font = common.merge_dicts(ticklabel_font,
                                                 self.ticklabel_font_default)
        self.legend_font = common.merge_dicts(legend_font,
                                              self.legend_font_default)
        self.title_font = common.merge_dicts(title_font,
                                             self.title_font_default)
        self.fp_axis_lab = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.axis_lab_font))
        self.fp_ticklabel = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.ticklabel_font))
        self.fp_legend = \
            mpl.font_manager.FontProperties(**copy.deepcopy(self.legend_font))
        self.fp_title = \
            mpl.font_manager.FontProperties(**copy.deepcopy(self.title_font))
        if hasattr(self, 'xticklabel_font'):
            self.xticklabel_font = \
                common.merge_dicts(self.xticklabel_font, self.ticklabel_font)
            self.fp_xticklabel = \
                mpl.font_manager.FontProperties(
                    **copy.deepcopy(self.xticklabel_font))
        else:
            self.fp_xticklabel = self.fp_ticklabel

        self.x = np.array(self.x, dtype=np.object)

        if self.grouped:
            self.grouped_y = self.y
            self.y = self.y[0]
            for i, gy in enumerate(self.grouped_y):
                self.grouped_y[i] = np.array(gy)
                setattr(self, 'y_g%u' % i, self.grouped_y[i])
            self.grouped_colors = self.color
            self.color = self.color[0]

        self.y = np.array(self.y)

        if hasattr(self, 'y2') and self.y2 is not None:
            self.y2 = np.array(self.y2)

        if self.do:
            self.plot()

    def reload(self):
        """
        Reloads the module and updates the class instance.
        """
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        """
        The total workflow of this class.
        Calls all methods in the correct order.
        """
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        """
        Executes all necessary tasks before plotting in the correct order.
        """
        self.set_categories()
        self.plots_order()
        self.set_colors()
        self.sort()
        self.by_plot()

    def do_plot(self):
        """
        Calls the plotting methods in the correct order.
        """
        self.set_figsize()
        self.init_fig()
        self.set_grid()
        self.make_plots()
        self.set_ylims()
        self.set_title()
        self.groups_legend()

    def post_plot(self):
        """
        Saves the plot into file, and closes the figure.
        """
        self.finish()

    def set_categories(self):
        """
        Sets a list with category indices (integers) of length equal of x,
        and sets dicts to translate between category names and indices.
        """

        # setting category names
        self.cnames = None

        if type(self.categories) is dict:
            # self.cnames: name -> number dict
            # each category has a number
            if (
                # we have the summary panel of all categories
                self.summary and
                any(
                    map(
                        lambda c: type(c) is tuple,
                        self.categories.keys()
                    )
                )
            ):
                self.cnames = dict(
                    map(
                        reversed,
                        enumerate(
                            [self.summary_name] +
                            sorted(list(set(self.categories.values())))
                        )
                    )
                )
            else:
                # only the categories, no summary panel
                self.cnames = dict(
                    map(
                        reversed,
                        enumerate(sorted(list(set(
                            self.categories.values()
                        ))))
                    )
                )

            # which element belongs to which category
            # list of integers
            # each category has a number
            self.cats = list(
                map(
                    lambda name: self.cnames[self.categories[name]],
                    self.x
                )
            )

        elif type(self.categories) is list:
            if type(self.categories[0]) is int:
                self.cats = self.categories
                self.cnames = dict(map(reversed, enumerate(self.cat_names)))
            else:
                self.cnames = dict(
                    map(reversed,
                        enumerate(sorted(list(set(self.categories))))))
                self.cats = list(
                    map(lambda name: self.cnames[name], self.categories))
        elif type(self.x[0]) is list:
            self.cats = []
            _x = []
            for i, c in enumerate(self.x):
                self.cats.extend([i] * len(c))
                _x.extend(c)
            self.x = _x
        else:
            self.cats = [0] * len(self.x)

        # this is the number of subplots including the summary plot
        self.numof_cats = len(set(self.cats)) + (1 if self.summary else 0)

        # if category names provided directly
        if self.cnames is None:
            if self.cat_names is not None:
                self.cnames = dict(zip(self.cat_names, list(set(self.cats))))
            else:
                self.cnames = dict(
                    map(lambda c: ('#%u' % c, c),
                        sorted(common.unique_list(self.cats))))
                if len(self.cnames) == 1:
                    self.cnames = {self.xlab: 0}

        # dict to translate category number to category name
        self.cnums = dict(map(reversed, iteritems(self.cnames)))

        self.cats = np.array(self.cats)

    def set_colors(self, colseries=''):
        """
        Compiles an array of colors equal length of x.
        """

        # calls for each group at grouped barplot
        if self.grouped and colseries == '':
            for i, gcol in enumerate(self.grouped_colors):
                setattr(self, 'color_g%u' % i, gcol)
                self.set_colors(colseries='_g%u' % i)
            self.col = self.col_g0
            return None

        colorattr = 'color%s' % colseries
        colattr = 'col%s' % colseries
        color = getattr(self, colorattr)
        ccol = None

        if type(color) is str:
            # only one color is given, all bars are the same color
            ccol = dict(map(lambda name: (name, color), self.cat_ordr))
        elif len(color) == len(self.cnames):
            # one color given for each category
            if type(color[0]) is str:
                # one color for one bar
                ccol = dict(
                    map(lambda c: (c[1], color[c[0]]),
                        enumerate(self.cat_ordr)))
            elif type(color[0]) is list:
                # 2 different colored bars overlayed
                setattr(self, colattr, [])
                for ccols in self.colors:
                    getattr(self, colattr).extend(ccols)
                setattr(self, colattr, np.array(getattr(self, colattr)))
        if type(color) not in _const.SIMPLE_TYPES and \
                len(color) == len(self.x):
            # color is a list of colors for each bar
            setattr(self, colattr, np.array(color))
        elif ccol is not None:
            # ccol is dict of colors for each category
            setattr(
                self, colattr,
                np.array(
                    list(map(lambda cnum: ccol[self.cnums[cnum]], self.cats))))

        # if got a second data series:
        if colseries == '' \
                and hasattr(self, 'color2') \
                and self.color2 is not None:
            self.set_colors(colseries='2')

    def plots_order(self):
        """
        Defines the order of the subplots.
        """

        if self.cat_ordr is None and self.cat_names is not None:

            self.cat_ordr = common.uniq_ord_list(sorted(self.cat_names))

        elif self.cat_ordr is None:

            self.cat_ordr = common.unique_list(sorted(self.cnames.keys()))

    def by_plot(self):
        """
        Sets list of lists with x and y values and colors by category.
        """

        # print(self.cnames)
        # print(self.cnums)
        # print(self.cats)
        # print(self.cat_ordr)
        # print(list(sorted(filter(lambda s: type(s[1]) is tuple,
        # zip(self.cats, self.x, self.y)),
        # key = lambda s: self.cat_ordr.index(self.cnums[s[0]]))))

        # variables to set
        attrs = ['x', 'y', 'col', 'cats']

        # if we have second series of data for shaded bars
        if hasattr(self, 'y2') and self.y2 is not None:
            attrs.extend(['y2', 'col2'])

        # for grouped barplot
        if self.grouped:

            for i in xrange(len(self.grouped_y)):

                attrs.append('y_g%u' % i)
                attrs.append('col_g%u' % i)

        for dim in attrs:

            attr = 'cat_%s' % dim
            # this only if we do summary plot
            # (the one on the left with summary of categories)
            if self.summary:

                setattr(
                    self,
                    attr,
                    [
                        list(
                            map(
                                lambda s: s[2],
                                sorted(
                                    filter(
                                        lambda s: isinstance(s[1], tuple),
                                           zip(
                                               self.cats,
                                               self.x,
                                               getattr(self, dim)
                                            )
                                        ),
                                    key = lambda s:
                                        self.cat_ordr.index(self.cnums[s[0]])
                                )
                            )
                        )
                    ]
                )

            else:
                setattr(self, attr, [])

        if self.summary:
            # only here we set the elements of the summary subplot
            # to the dummy category of summary
            self.cats = list(
                map(
                    lambda name:
                        self.cnames[self.summary_name]
                            if type(name) is tuple else
                        self.cnames[self.categories[name]],
                    self.x
                )
            )

        for dim in attrs:
            attr = 'cat_%s' % dim
            # this is what is always necessary:
            setattr(self, attr,
                    getattr(self, attr) +
                    list(map(
                        lambda name:
                        list(map(
                            lambda n_lab:
                                n_lab[1],
                            filter(
                                lambda n_lab:
                                    n_lab[0] == self.cnames[name],
                                zip(self.cats, getattr(self, dim))
                            )
                        )),
                        self.cat_ordr
                    ))
                    )
        self.xlabs = ([self.summary_name]
                      if self.summary else []) + self.cat_ordr


    def sort(self):
        """
        Finds the defined or default order, and
        sorts the arrays x, y and col accordingly.
        """
        if type(self.order) is str:
            if self.order == 'x':
                self.ordr = np.array(self.x.argsort())
            elif self.order == 'y':
                self.ordr = np.array(self.y.argsort())
        elif hasattr(self.order, '__iter__') and \
                len(set(self.order) & set(self.x)) == len(self.x):
            self.ordr = np.array(list(map(
                lambda i:
                # this is ugly, but needed a quick
                # solution when introducing tuples...
                list(self.x).index(i),
                self.order
            )))
        else:
            self.ordr = np.array(xrange(len(self.x)))
        # descending order
        if self.desc:
            self.ordr = self.ordr[::-1]
        # ordering all variables
        self.x = self.x[self.ordr]
        self.y = self.y[self.ordr]
        self.col = self.col[self.ordr]
        self.cats = self.cats[self.ordr]
        # if we have second data series (shaded parts of bars)
        # order also those variables
        if hasattr(self, 'y2') and self.y2 is not None:
            self.y2 = self.y2[self.ordr]
        if hasattr(self, 'col2'):
            self.col2 = self.col2[self.ordr]
        # for grouped barplot
        if self.grouped:
            for g in xrange(len(self.grouped_y)):
                yattr = 'y_g%u' % g
                colattr = 'col_g%u' % g
                setattr(self, yattr, getattr(self, yattr)[self.ordr])
                setattr(self, colattr, getattr(self, colattr)[self.ordr])


    def set_figsize(self):
        """
        Converts width and height to a tuple which can be used as figsize.
        """

        if hasattr(self, 'width') and hasattr(self, 'height'):

            self.figsize = (self.width, self.height)

    def init_fig(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """

        self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)

    def set_grid(self):
        """
        Sets up a grid according to the number of subplots,
        with proportions according to the number of elements
        in each subplot.
        """

        self.gs = mpl.gridspec.GridSpec(
            2,
            self.numof_cats,
            height_ratios=[1, 0],
            width_ratios=list(map(len, self.cat_x)))
        self.axes = [[None] * self.numof_cats, [None] * self.numof_cats]

    def get_subplot(self, i, j=0):
        if self.axes[j][i] is None:
            self.axes[j][i] = self.fig.add_subplot(self.gs[j, i])
        self.ax = self.axes[j][i]

    def make_plots(self):
        """
        Does the actual plotting.
        """
        original_width = self.bar_args['width']
        if self.grouped:
            self.bar_args['width'] = self.bar_args['width'] / float(
                len(self.grouped_y))

        width = self.bar_args['width']
        correction = len(self.grouped_y) * width / 2.0 if self.grouped else 0.0

        for i, x in enumerate(self.cat_x):
            self.get_subplot(i)
            xcoo = np.arange(len(x)) - self.bar_args['width'] / 2.0
            xtlabs = \
                np.array(
                    list(
                        map(
                            lambda l:
                                l if type(l) is not tuple
                            else self.categories[l]
                            if self.summary
                            else l[0],
                            x
                        )
                    )
                )
            self.ax.bar(x = xcoo - correction,
                        height = self.cat_y[i],
                        color = self.cat_col[i],
                        tick_label = xtlabs,
                        **copy.deepcopy(self.bar_args))
            if self.grouped:
                for j in xrange(1, len(self.grouped_y)):
                    self.ax.bar(x = xcoo + width * j - correction,
                                height = getattr(self, 'cat_y_g%u' % j)[i],
                                color = getattr(self, 'cat_col_g%u' % j)[i],
                                **copy.deepcopy(self.bar_args))
            if hasattr(self, 'y2') and self.y2 is not None:
                self.ax.bar(x = xcoo - correction,
                            height = self.cat_y2[i],
                            color = self.cat_col2[i],
                            **copy.deepcopy(self.bar_args))
            self.labels()
            self.ax.xaxis.grid(False)
            self.ax.set_xlim([-1, max(xcoo) + 0.5])
            if self.ylim is not None:
                self.ax.set_ylim(self.ylim)
            if self.grouped:
                self.ax.set_xticks(self.ax.get_xticks() + correction)
            if self.ylog:
                self.ax.set_yscale('symlog')
                self.ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

            self.ax.yaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
            self.ax.xaxis.grid(False)
            self.ax.set_axisbelow(True)
            self.ax.set_facecolor('#EAEAF2')
            list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
            self.ax.tick_params(which='both', length=0)

            self.get_subplot(i, 1)
            self.ax.xaxis.set_ticklabels([])
            self.ax.yaxis.set_ticklabels([])
            self.ax.set_xlabel(self.xlabs[i], fontproperties=self.fp_axis_lab)
            self.ax.xaxis.label.set_verticalalignment('bottom')
            list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
            self.ax.tick_params(which='both', length=0)

        self.bar_args['width'] = original_width

    def set_ylims(self):
        if self.uniform_ylim:
            maxy = max(
                map(lambda ax: ax.get_ylim()[1], self.axes[0][(
                    1 if self.summary else 0):]))
            _ = list(
                map(lambda ax: ax.set_ylim([0, maxy]), self.axes[0][(
                    1 if self.summary else 0):]))

    def labels(self):
        """
        Sets properties of axis labels and ticklabels.
        """
        list(
            map(lambda tick: tick.label.set_fontproperties(self.fp_xticklabel) or (self.lab_angle == 0 or self.lab_angle == 90) and (tick.label.set_rotation(self.lab_angle) or tick.label.set_horizontalalignment('center')),
                self.ax.xaxis.get_major_ticks()))
        list(
            map(lambda tick: tick.label.set_fontproperties(self.fp_ticklabel),
                self.ax.yaxis.get_major_ticks()))
        self.ax.set_ylabel(self.ylab, fontproperties=self.fp_axis_lab)
        # self.ax.yaxis.label.set_fontproperties(self)

    def set_title(self):
        """
        Sets the main title.
        """
        if self.maketitle:
            self.title_text = self.fig.suptitle(self.title)
            self.title_text.set_fontproperties(self.fp_title)
            self.title_text.set_horizontalalignment(self.title_halign)
            self.title_text.set_verticalalignment(self.title_valign)

    def groups_legend(self):
        if self.grouped:
            lhandles = \
                list(
                    map(
                        lambda g:
                            mpl.patches.Patch(color=g[1], label=g[0]),
                        zip(self.group_labels, self.grouped_colors)
                    )
                )
            broadest_ax = max(self.axes[0],
                              key=lambda ax: len(ax.get_xticks()))
            broadest_ax.legend(
                handles=lhandles,
                prop=self.fp_legend,
                frameon=False,
                loc=self.legloc)

    def align_x_labels(self):
        self.lowest_ax = min(self.axes[0],
                             key=lambda ax: ax.xaxis.label.get_position()[1])
        self.minxlabcoo = self.lowest_ax.xaxis.label.get_position()[1]
        self.lowest_xlab_dcoo = self.lowest_ax.transData.transform(
            self.lowest_ax.xaxis.label.get_position())
        list(
            map(
                lambda ax:
                ax.xaxis.set_label_coords(
                    self.fig.transFigure.inverted().transform(
                        ax.transAxes.transform((0.5, 0.5)))[0],
                    self.fig.transFigure.inverted().transform(
                        self.lowest_xlab_dcoo)[1],
                    transform=self.fig.transFigure
                ),
                self.axes[0]
            )
        )
        for ax in self.axes[0]:
            ax.xaxis._autolabelpos = False

    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        if self.maketitle:
            self.fig.subplots_adjust(top=0.85)

        self.cvs.draw()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.cvs.print_figure(self.pdf)

        self.pdf.close()
        self.fig.clf()

# ## ## ##
# ## ## ##


class Barplot(Plot):
    def __init__(self,
                 x,
                 y,
                 data=None,
                 fname=None,
                 font_family='Helvetica Neue LT Std',
                 font_style='normal',
                 font_weight='normal',
                 font_variant='normal',
                 font_stretch='normal',
                 xlab='',
                 ylab='',
                 axis_lab_size=10.0,
                 lab_angle=90,
                 lab_size=(9, 9),
                 color='#007b7f',
                 order=False,
                 desc=True,
                 legend=None,
                 fin=True,
                 y_break=None,
                 rc={},
                 palette=None,
                 context='poster',
                 do_plot=True,
                 **kwargs):
        '''
        y_break : tuple
        If not None, the y-axis will have a break. 2 floats in the tuple, < 1.0,
        mean the lower and upper proportion of the plot shown. The part between
        them will be hidden. E.g. y_break = (0.3, 0.1) shows the lower 30% and
        upper 10%, but 60% in the middle will be cut out.
        '''
        for k, v in iteritems(locals()):
            setattr(self, k, v)
        self.sns = sns
        self.rc = self.rc or {
            'lines.linewidth': 1.0,
            'patch.linewidth': 0.0,
            'grid.linewidth': 1.0
        }
        super(Barplot, self).__init__()
        self.color = self.color or self.palette[0][0]
        if type(self.color) is list:
            self.palette = sns.color_palette(self.color)
        elif self.color is not None:
            self.palette = sns.color_palette([self.color] * len(self.x))
        self.color = None
        if self.do_plot:
            self.plot(**kwargs)

    def plot(self, x=None, y=None, **kwargs):
        if x is not None:
            self.x = x
        if y is not None:
            self.y = y
        if type(self.x) is list or type(self.x) is tuple:
            self.x = np.array(self.x)
        if type(self.y) is list or type(self.y) is tuple:
            self.y = np.array(self.y)
        self.seaborn_style()
        self.fig, self.ax = plt.subplots()
        self.sort()
        if self.y_break:
            self._break_y_gs()
        self.ax = sns.barplot(
            self.x,
            y=self.y,
            data=None,
            color=self.color,
            order=self.ordr,
            ax=self.ax,
            palette=self.palette,
            #fontproperties = self.fp,
            **kwargs)
        if self.y_break:
            self._break_y_axis(**kwargs)
        self.labels()
        if self.fin:
            self.finish()

    def sort(self):
        colcyc = itertools.cycle(list(self.palette))
        palcyc = [next(colcyc) for _ in xrange(len(self.x))]
        if self.order == 'x':
            self.ordr = np.array([self.x[i] for i in self.x.argsort()])
            self.palette = sns.color_palette(
                [palcyc[i] for i in self.x.argsort()])
        elif self.order == 'y':
            self.ordr = np.array([self.x[i] for i in self.y.argsort()])
            self.palette = sns.color_palette(
                [palcyc[i] for i in self.y.argsort()])
        elif len(set(self.order) & set(self.x)) == len(self.x):
            self.ordr = np.array(self.order)
            xl = list(self.x)
            self.palette = sns.color_palette(
                [palcyc[xl.index(i)] for i in self.ordr])
        else:
            self.ordr = self.x
        if self.desc:
            self.ordr = self.ordr[::-1]
            self.palette = sns.color_palette(list(self.palette)[::-1])

    def _break_y_gs(self):
        self.gs = gridspec.GridSpec(
            2,
            1,
            height_ratios=[
                self.y_break[1] / sum(self.y_break),
                self.y_break[0] / sum(self.y_break)
            ])
        self.fig = plt.figure()
        self.ax2 = self.fig.add_subplot(self.gs[0])
        self.ax = self.fig.add_subplot(self.gs[1])

    def __break_y_axis(self, **kwargs):
        self.ax2 = self.sns.barplot(
            self.x,
            y=self.y,
            data=None,
            color=self.color,
            order=self.ordr,
            ax=self.ax2,
            palette=self.palette,
            **kwargs)
        self.ax2.yaxis.set_major_locator(
            ticker.MaxNLocator(
                nbins=int(9 / sum(self.y_break) + 1), steps=[1, 2, 5, 10]))
        self._originalYticks = self.ax2.get_yticks()
        ymin, ymax = self.ax.get_ylim()
        ymax = min(ytick for ytick in self.ax.get_yticks()
                   if ytick > max(self.y))
        self.ax.set_ylim((ymin, ymax * self.y_break[0]))
        self.ax2.set_ylim((ymax - ymax * self.y_break[1], ymax))
        self.lower_y_min, self.lower_y_max = self.ax.get_ylim()
        self.upper_y_min, self.upper_y_max = self.ax2.get_ylim()
        plt.subplots_adjust(hspace=0.08)
        self.ax2.spines['bottom'].set_visible(False)
        plt.setp(self.ax2.xaxis.get_majorticklabels(), visible=False)
        self.ax.spines['top'].set_visible(False)
        self.ax.set_yticks([
            yt for yt in self._originalYticks
            if yt >= self.lower_y_min and yt <= self.lower_y_max
        ])
        self.ax2.set_yticks([
            yt for yt in self._originalYticks
            if yt >= self.upper_y_min and yt <= self.upper_y_max
        ])

    def _break_y_axis(self, **kwargs):
        self.ax2 = self.sns.barplot(
            self.x,
            y=self.y,
            data=None,
            color=self.color,
            order=self.ordr,
            ax=self.ax2,
            palette=self.palette,
            **kwargs)
        self.ax2.yaxis.set_major_locator(
            ticker.MaxNLocator(
                nbins=int(9 / sum(self.y_break) + 1), steps=[1, 2, 5, 10]))
        self._originalYticks = self.ax2.get_yticks()
        ymin, ymax = self.ax.get_ylim()
        yticks = [
            ytick for ytick in self.ax.get_yticks() if ytick > max(self.y)
        ]
        if len(yticks) > 0:
            ymax = min(yticks)
        else:
            ymax = max(self.ax.get_yticks())
        self.ax.set_ylim((ymin, ymax * self.y_break[0]))
        self.ax2.set_ylim((ymax - ymax * self.y_break[1], ymax))
        self.lower_y_min, self.lower_y_max = self.ax.get_ylim()
        self.upper_y_min, self.upper_y_max = self.ax2.get_ylim()
        plt.subplots_adjust(hspace=0.08)
        self.ax2.spines['bottom'].set_visible(False)
        plt.setp(self.ax2.xaxis.get_majorticklabels(), visible=False)
        self.ax.spines['top'].set_visible(False)
        self.ax.set_yticks([
            yt for yt in self._originalYticks
            if yt >= self.lower_y_min and yt <= self.lower_y_max
        ])
        self.ax2.set_yticks([
            yt for yt in self._originalYticks
            if yt >= self.upper_y_min and yt <= self.upper_y_max
        ])
        # further adjusting of upper ylims:
        yticks = [yt for yt in self.ax2.get_yticks() if yt > max(self.y)]
        if len(yticks) > 0:
            ymax = min(yticks)
            self.ax.set_ylim((ymin, ymax * self.y_break[0]))
            self.ax2.set_ylim((ymax - ymax * self.y_break[1], ymax))
            self.lower_y_min, self.lower_y_max = self.ax.get_ylim()
            self.upper_y_min, self.upper_y_max = self.ax2.get_ylim()
            self.ax2.set_yticks([
                yt for yt in self._originalYticks
                if yt >= self.upper_y_min and yt <= self.upper_y_max
            ])

    def seaborn_style(self, context=None, rc=None):
        self.sns.set(font=self.font_family, rc=rc or self.rc)
        self.sns.set_context(context or self.context, rc=rc or self.rc)

    def labels(self):
        for tick in self.ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(self.lab_size[0])
        for tick in self.ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(self.lab_size[1])
        if self.y_break:
            for tick in self.ax2.yaxis.get_major_ticks():
                tick.label.set_fontsize(self.lab_size[1])
        self.ax.set_ylabel(self.ylab)
        self.ax.yaxis.get_label().set_fontproperties(self.fp)
        self.ax.yaxis.get_label().set_fontsize(self.axis_lab_size)
        self.ax.set_xlabel(self.xlab)
        self.ax.xaxis.get_label().set_fontproperties(self.fp)
        self.ax.xaxis.get_label().set_fontsize(self.axis_lab_size)
        plt.setp(self.ax.xaxis.get_majorticklabels(), rotation=self.lab_angle)
        if type(self.legend) is dict:
            legend_patches = [
                mpatches.Patch(
                    color=col, label=lab)
                for lab, col in iteritems(self.legend)
            ]
            self.ax.legend(handles=legend_patches)


def boxplot(data,
            labels,
            xlab,
            ylab,
            fname,
            fontfamily='Helvetica Neue LT Std',
            textcol='black',
            violin=False):
    fig, ax = plt.subplots()
    sns.set(font=fontfamily)
    if violin:
        ax = sns.violinplot(
            data,
            names=labels,
            color=embl_colors,
            linewidth=0.1,
            saturation=0.66)
    else:
        ax = sns.boxplot(
            data,
            names=labels,
            color=embl_colors,
            linewidth=0.1,
            saturation=0.66)
    ax.set_xlabel(
        xlab,
        weight='light',
        fontsize=12,
        variant='normal',
        color=textcol,
        stretch='normal')
    ax.set_ylabel(
        ylab,
        weight='light',
        fontsize=12,
        variant='normal',
        color=textcol,
        stretch='normal')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
        tick.label.set_color(textcol)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(11)
        tick.label.set_color(textcol)
    fig.savefig(fname)


class StackedBarplot(object):
    def __init__(self,
                 x,
                 y,
                 fname,
                 names,
                 colors,
                 xlab='',
                 ylab='',
                 title='',
                 title_halign='center',
                 title_valign='top',
                 bar_args={},
                 axis_lab_font={},
                 ticklabel_font={},
                 title_font={},
                 legend_font={},
                 lab_angle=90,
                 figsize=(9, 6),
                 legend=True,
                 order=False,
                 desc=True):

        for k, v in iteritems(locals()):
            setattr(self, k, v)

        self.bar_args_default = {
            'width': 0.8,
            'edgecolor': 'none',
            'linewidth': 0.0,
            'align': 'center'
        }
        self.axis_lab_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'x-large'
        }
        self.ticklabel_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'large'
        }
        self.legend_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'small'
        }
        self.title_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'xx-large'
        }

        self.bar_args = common.merge_dicts(bar_args, self.bar_args_default)
        self.axis_lab_font = common.merge_dicts(axis_lab_font,
                                                self.axis_lab_font_default)
        self.ticklabel_font = common.merge_dicts(ticklabel_font,
                                                 self.ticklabel_font_default)
        self.legend_font = common.merge_dicts(legend_font,
                                              self.legend_font_default)
        self.title_font = common.merge_dicts(title_font,
                                             self.title_font_default)

        self.plot()

    def reload(self):
        """
        Reloads the module and updates the class instance.
        """
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        """
        The total workflow of this class.
        Calls all methods in the correct order.
        """
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        self.set_fontproperties()
        self.sort()
        self.set_figsize()
        self.init_fig()

    def do_plot(self):
        self.set_background()
        self.set_gridlines()
        self.make_plot()
        self.make_legend()
        self.set_ticklabels()
        self.set_axis_labels()
        self.set_title()

    def post_plot(self):
        self.finish()

    def set_fontproperties(self):
        self.fp_axis_lab = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.axis_lab_font))
        self.fp_ticklabel = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.ticklabel_font))
        self.fp_legend = \
            mpl.font_manager.FontProperties(**copy.deepcopy(self.legend_font))
        self.fp_title = \
            mpl.font_manager.FontProperties(**copy.deepcopy(self.title_font))

    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        if hasattr(self, 'width') and hasattr(self, 'height'):
            self.figsize = (self.width, self.height)

    def init_fig(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """
        self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)
        self.ax = self.fig.add_subplot(1, 1, 1)

    def set_background(self):
        self.ax.yaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
        self.ax.xaxis.grid(False)
        # self.ax.yaxis.grid(True, color = '#FFFFFF', linewidth = 2)
        self.ax.set_axisbelow(True)

    def set_gridlines(self):
        self.ax.set_facecolor('#EAEAF2')
        list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
        self.ax.tick_params(which='both', length=0)

    def sort(self):
        self.x = np.array(self.x)
        self.y = list(map(np.array, self.y))
        self.total = reduce(lambda l1, l2: l1.__add__(l2), self.y)

        if type(self.order) is str:
            if self.order == 'x':
                self.ordr = self.x.argsort()
            elif self.order == 'y':
                self.ordr = self.total.argsort()
        elif type(self.order) is int:
            self.ordr = self.y[i].argsort()
        elif len(set(self.order) & set(self.x)) == len(self.x):
            self.ordr = np.array(
                list(map(lambda l: np.where(self.x == l)[0][0], self.order)))
        else:
            self.ordr = np.arange(len(self.x))
        if self.desc:
            self.ordr = self.ordr[::-1]
        self.x = self.x[self.ordr]
        self.y = list(map(lambda iy: iy[self.ordr], self.y))

    def make_plot(self):

        self.xcoo = np.arange(len(self.x))

        for j in xrange(len(self.y), 0, -1):
            this_level = reduce(lambda l1, l2: l1.__add__(l2), self.y[:j])
            self.ax.bar(x=self.xcoo,
                        height=this_level,
                        tick_label=self.x,
                        color=self.colors[j - 1],
                        label=self.names[j - 1],
                        **copy.deepcopy(self.bar_args))

    def set_ticklabels(self):
        list(
            map(lambda l: l.set_fontproperties(self.fp_ticklabel) or l.set_rotation(self.lab_angle),
                self.ax.xaxis.get_majorticklabels()))
        list(
            map(lambda l: l.set_fontproperties(self.fp_ticklabel),
                self.ax.yaxis.get_majorticklabels()))

    def set_axis_labels(self):
        self.ax.set_ylabel(self.ylab, fontproperties=self.fp_axis_lab)
        self.ax.set_xlabel(self.xlab, fontproperties=self.fp_axis_lab)

    def set_xlim(self):
        self.ax.set_xlim([-1, max(self.xcoo) + 0.5])

    def set_title(self):
        """
        Sets the main title.
        """
        self.title_text = self.fig.suptitle(self.title)
        self.title_text.set_fontproperties(self.fp_title)
        self.title_text.set_horizontalalignment(self.title_halign)
        self.title_text.set_verticalalignment(self.title_valign)

    def make_legend(self):
        if self.legend:
            self.lhandles = \
                list(
                    map(
                        lambda i:
                            mpl.patches.Patch(color=self.colors[
                                              i], label=self.names[i]),
                        xrange(len(self.y))
                    )
                )
            self.leg = self.ax.legend(
                handles=self.lhandles, prop=self.fp_legend, frameon=False)
            self.leg.get_title().set_fontproperties(self.fp_axis_lab)

    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        self.fig.subplots_adjust(top=0.85)
        self.cvs.draw()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.cvs.print_figure(self.pdf)

        self.pdf.close()
        self.fig.clf()

## ## ##


class ScatterPlus(object):
    def __init__(self,
                 x,
                 y,
                 size=None,
                 color='#114477',
                 labels=None,
                 xlog=False,
                 ylog=False,
                 xlim=None,
                 ylim=None,
                 xtickscale=None,
                 ytickscale=None,
                 legscale=None,
                 fname=None,
                 confi=True,
                 title_font={},
                 ticklabel_font={},
                 legend_font={},
                 axis_lab_font={},
                 annot_font={},
                 xlab='',
                 ylab='',
                 axis_lab_size=10.0,
                 min_size=5.0,
                 max_size=30.0,
                 log_size=False,
                 alpha=0.5,
                 size_scaling=0.8,
                 lab_angle=90,
                 order=False,
                 desc=True,
                 legend=True,
                 legtitle='',
                 legstrip=(None, None),
                 color_labels=[],
                 legloc=4,
                 size_to_value=lambda x: x,
                 value_to_size=lambda x: x,
                 figsize=(10.0, 7.5),
                 title='',
                 title_halign='center',
                 title_valign='top',
                 fin=True,
                 rc={},
                 **kwargs):

        for k, v in iteritems(locals()):
            setattr(self, k, v)
        self.rc_default = {
            'lines.linewidth': 1.0,
            'patch.linewidth': 0.0,
            'grid.linewidth': 1.0
        }
        self.rc = common.merge_dicts(self.rc, self.rc_default)

        self.axis_lab_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'xx-large'
        }
        self.ticklabel_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'x-large'
        }
        self.legend_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'medium'
        }
        self.title_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'xx-large'
        }
        self.annot_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'large'
        }

        self.axis_lab_font = common.merge_dicts(axis_lab_font,
                                                self.axis_lab_font_default)
        self.ticklabel_font = common.merge_dicts(ticklabel_font,
                                                 self.ticklabel_font_default)
        self.legend_font = common.merge_dicts(legend_font,
                                              self.legend_font_default)
        self.title_font = common.merge_dicts(title_font,
                                             self.title_font_default)
        self.annot_font = common.merge_dicts(annot_font,
                                             self.annot_font_default)

        self.x = np.array(x)
        self.y = np.array(y)
        self.labels = np.array(labels)

        self.plot()

    def reload(self):
        """
        Reloads the module and updates the class instance.
        """
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        self.set_colors()
        self.set_fontproperties()
        self.set_size()

    def do_plot(self):
        self.set_figsize()
        self.init_fig()
        self.set_log()
        self.make_plot()
        self.confidence_interval()
        self.axes_limits()
        self.set_ticklocs()
        self.set_background()
        self.set_gridlines()
        self.annotations()
        self.axes_labels()
        self.axes_limits()
        self.axes_ticklabels()
        self.set_title()
        self.axes_limits()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        self.axes_limits()
        self.remove_annotation_overlaps()
        self.make_legend()

    def post_plot(self):
        self.finish()

    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        if hasattr(self, 'width') and hasattr(self, 'height'):
            self.figsize = (self.width, self.height)

    def init_fig(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """
        self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)
        self.ax = self.fig.add_subplot(1, 1, 1)

    def set_colors(self):
        if type(self.color) is str:
            self.color = [self.color] * len(self.x)
        if type(self.color) is dict:
            self.color = \
                list(
                    map(
                        lambda c:
                            self.color[c],
                        self.categories
                    )
                )

    def set_fontproperties(self):
        self.fp_axis_lab = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.axis_lab_font))
        self.fp_ticklabel = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.ticklabel_font))
        self.fp_legend = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.legend_font))
        self.fp_title = \
            mpl.font_manager.FontProperties(**copy.deepcopy(self.title_font))
        self.fp_annot = \
            mpl.font_manager.FontProperties(**copy.deepcopy(self.annot_font))

    def set_title(self):
        """
        Sets the main title.
        """
        self.title_text = self.fig.suptitle(self.title)
        self.title_text.set_fontproperties(self.fp_title)
        self.title_text.set_horizontalalignment(self.title_halign)
        self.title_text.set_verticalalignment(self.title_valign)

    def scale(self, scale=None, q=1):
        scale = [1, 2, 5] if scale is None else scale

        def _scaler(scale, q):
            while True:
                for i in scale:
                    yield i * q
                q *= 10.0

        _scale = _scaler(scale, q)
        return _scale

    def set_size(self):
        if self.size is None:
            self.size = self.min_size

        if type(self.size) in common.numeric_types:
            self.size = [self.size] * len(self.x)

        self.size_values = np.array(self.size)

        self.size = self.values_to_sizes(self.size)

    def set_log(self):
        if self.ylog:
            self.ax.set_yscale('symlog' if self.ylog == 'symlog' else 'log')
        if self.xlog:
            self.ax.set_xscale('symlog' if self.xlog == 'symlog' else 'log')

    def make_plot(self):

        self.scatter = self.ax.scatter(
            self.x,
            self.y,
            s=self.size,
            c=self.color,
            alpha=self.alpha,
            edgecolors='none')

    def annotations(self):
        if self.labels is not None:
            self.annots = []
            dists = move_labels()
            for label, xx, yy, yf in \
                    zip(self.labels, self.x, self.y, self.y_fit):
                dst = next(dists)
                d = 1.0 if yy > 10**yf else -1.0
                coo = ((-7 - dst) * d / 3.0, (21 + dst) * d)
                self.annots.append(
                    self.ax.annotate(
                        label,
                        xy=(xx, yy),
                        xytext=coo,
                        xycoords='data',
                        textcoords='offset points',
                        ha='center',
                        va='bottom',
                        color='#114477',
                        arrowprops=dict(
                            arrowstyle='-',
                            connectionstyle='arc,rad=.0',
                            color='#114477',
                            edgecolor='#114477',
                            alpha=1.0,
                            visible=True,
                            linewidth=0.2), ))
            for ann in self.annots:
                ann.set_fontproperties(self.fp_annot)

    def set_ticklocs(self):
        if self.xlog:
            xscaler = self.scale(self.xtickscale)
            self.xtickloc = []
            while True:
                tickloc = next(xscaler)
                if tickloc >= self._xlim[0]:
                    self.xtickloc.append(tickloc)
                if tickloc > self._xlim[1]:
                    break
            self._xtickloc = self.ax.set_xticks(self.xtickloc)
        if self.ylog:
            yscaler = self.scale(self.ytickscale)
            self.ytickloc = []
            while True:
                tickloc = next(yscaler)
                if tickloc >= self._ylim[0]:
                    self.ytickloc.append(tickloc)
                if tickloc > self._ylim[1]:
                    break
            self._ytickloc = self.ax.set_yticks(self.ytickloc)

    def axes_limits(self, xlim=None, ylim=None):
        xlim = xlim if xlim is not None else self.xlim \
            if self.xlim is not None else self.ax.get_xlim()
        ylim = ylim if ylim is not None else self.ylim \
            if self.ylim is not None else self.ax.get_ylim()
        if xlim is not None:
            self._xlim = self.ax.set_xlim(xlim)
        if ylim is not None:
            self._ylim = self.ax.set_ylim(ylim)

    def set_background(self):
        self.ax.grid(True, color='#FFFFFF', lw=1, ls='solid')
        # self.ax.yaxis.grid(True, color = '#FFFFFF', linewidth = 2)
        self.ax.set_axisbelow(True)

    def set_gridlines(self):
        self.ax.set_facecolor('#EAEAF2')
        list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
        self.ax.tick_params(which='both', length=0)

    def axes_ticklabels(self):
        if self.xlog:
            self.xticklabs = []
            for i, t in enumerate(
                    list(self.ax.xaxis.get_major_locator().locs)):
                tlab = str(int(t)) if t - int(t) == 0 or t >= 10.0 else str(t)
                self.xticklabs.append(tlab[:-2]
                                      if tlab.endswith('.0') else tlab)
            self._xticklabels = self.ax.set_xticklabels(self.xticklabs)
        if self.ylog:
            self.yticklabs = []
            for i, t in enumerate(
                    list(self.ax.yaxis.get_major_locator().locs)):
                tlab = str(int(t)) if t - int(t) == 0 or t >= 10.0 else str(t)
                self.yticklabs.append(tlab[:-2]
                                      if tlab.endswith('.0') else tlab)
            self._yticklabels = self.ax.set_yticklabels(self.yticklabs)
        for t in self.ax.xaxis.get_major_ticks():
            t.label.set_fontproperties(self.fp_ticklabel)
        for t in self.ax.yaxis.get_major_ticks():
            t.label.set_fontproperties(self.fp_ticklabel)

    def confidence_interval(self):
        # the points:
        def remove_inf(a, log=False):
            return \
                np.array(
                    list(
                        map(
                            lambda i:
                                0.0 if np.isinf(i) or i == 0.0
                            else np.log10(i) if log
                            else i,
                            a
                        )
                    )
                )

        self._x = remove_inf(self.x, self.xlog)
        self._y = remove_inf(self.y, self.ylog)
        # (log)linear fit with confidence and prediction interval:
        (self.m, self.b), self.V = np.polyfit(self._x, self._y, 1, cov=True)
        self.n = self._x.size
        self.y_fit = np.polyval((self.m, self.b), self._x)
        self.df = self.n - 2
        self.t = stats.t.ppf(0.95, self.df)
        self.resid = self._y - self.y_fit
        self.chi2 = np.sum((self.resid / self.y_fit)**2)
        self.chi2_red = self.chi2 / self.df
        self.s_err = np.sqrt(np.sum(self.resid**2) / self.df)
        self.x2 = np.linspace(np.min(self._x), np.max(self._x), 100)
        self.y2 = np.linspace(np.min(self.y_fit), np.max(self.y_fit), 100)
        # confidence interval
        self.ci = self.t * self.s_err * np.sqrt(1 / self.n + (
            self.x2 - np.mean(self._x))**2 / np.sum((self._x - np.mean(self._x)
                                                     )**2))
        # prediction interval
        self.pi = self.t * self.s_err * np.sqrt(1 + 1 / self.n + (
            self.x2 - np.mean(self._x))**2 / np.sum((self._x - np.mean(self._x)
                                                     )**2))
        # regression line
        self.rline_x = [10**xx for xx in self._x] if self.xlog else self._x
        self.rline_y = [10**yy
                        for yy in self.y_fit] if self.ylog else self.y_fit
        self.rline = self.ax.plot(
            self.rline_x, self.rline_y, '-', color='#B6B7B9', alpha=0.5)
        # confidence interval
        self.ci_rfill_x = [10**xx for xx in self.x2] if self.xlog else self.x2
        self.ci_rfill_y_upper = [10**yy for yy in (self.y2 + self.ci)] \
            if self.ylog else self.y2 + self.ci
        self.ci_rfill_y_lower = [10**yy for yy in (self.y2 - self.ci)] \
            if self.ylog else self.y2 - self.ci
        self.rfill = self.ax.fill_between(
            self.ci_rfill_x,
            self.ci_rfill_y_upper,
            self.ci_rfill_y_lower,
            color='#B6B7B9',
            edgecolor='',
            alpha=0.2)
        # prediction intreval
        self.pi_y_upper = [10**yy for yy in (self.y2 + self.pi)] \
            if self.ylog else self.y2 + self.pi
        self.pi_y_lower = [10**yy for yy in (self.y2 - self.pi)] \
            if self.ylog else self.y2 - self.pi
        self.pilowerline = self.ax.plot(
            self.ci_rfill_x,
            self.pi_y_lower,
            '--',
            color='#B6B7B9',
            linewidth=0.5)
        self.piupperline = self.ax.plot(
            self.ci_rfill_x,
            self.pi_y_upper,
            '--',
            color='#B6B7B9',
            linewidth=0.5)

    def values_to_sizes(self, values):
        """
        Transformation converts from size values in data dimension
        to dimension of the size graphical parameter.
        """
        values = np.array(values)
        svals = self.size_values
        mins = self.min_size
        maxs = self.max_size

        if self.log_size:
            values = np.log2(values)
            svals = np.log2(svals)
            mins = np.log2(mins)
            maxs = np.log2(maxs)

        return \
            (np.array(values) - min(svals)) / \
            float(max(svals) - min(svals)) * \
            (self.max_size - self.min_size) + self.min_size

    def make_legend(self):
        if self.size is not None and self.legend:

            self.lhandles1 = []
            self.llabels1 = []
            self.lhandles2 = []
            self.llabels2 = []
            sizemin = min(self.size_values)
            self.leglower = 10**int(np.floor(np.log10(sizemin)))
            self.leglower = self.leglower \
                if abs(sizemin - self.leglower) < abs(sizemin - self.leglower * 10) \
                else self.leglower * 10
            self.legscaler = self.scale(self.legscale, q=self.leglower)
            self.legsizes = []

            while True:
                legvalue = next(self.legscaler)
                self.legsizes.append(legvalue)
                if legvalue >= max(self.size_values):
                    break

            self.legsizes = self.legsizes[self.legstrip[0]:-self.legstrip[
                1] if self.legstrip[1] is not None else None]

            if len(self.legsizes) > 1 and \
                    abs(max(self.size_values) - self.legsizes[-1]) > \
                    abs(max(self.size_values) - self.legsizes[-2]):
                self.legsizes = self.legsizes[:-1]

            self.real_legsizes = np.array(self.legsizes)

            self.legsizes = self.values_to_sizes(self.legsizes)

            self.legsizes = np.sqrt(self.legsizes / np.pi)

            for lab, col in self.color_labels:
                self.lhandles1.append(
                    mpl.patches.Patch(
                        color=col, label=lab, alpha=.5, linewidth=0.0))
                self.llabels1.append(lab)

            for i, s in enumerate(self.legsizes):
                rs = self.real_legsizes[i]
                self.lhandles2.append(
                    mpl.legend.Line2D(
                        range(1),
                        range(1),
                        color='none',
                        marker='o',
                        markersize=s,
                        markerfacecolor='#114477',
                        markeredgecolor='none',
                        alpha=.5,
                        label=str(int(rs))
                        if rs - int(rs) == 0 or rs >= 10.0 else str(rs)))
                self.llabels2.append(
                    str(int(rs))
                    if rs - int(rs) == 0 or rs >= 10.0 else str(rs))

            self.leg2 = self.ax.legend(
                self.lhandles2 + self.lhandles1,
                self.llabels2 + self.llabels1,
                title=self.legtitle,
                labelspacing=.9,
                borderaxespad=.9,
                loc=self.legloc,
                prop=self.fp_legend,
                markerscale=1.0,
                frameon=False,
                numpoints=1)

            # self.ax.add_artist(self.leg2)

            #bbleg2 = self.leg2.legendPatch.get_bbox().inverse_transformed(self.fig.transFigure)

            #upperright = bbleg2.corners()[3]

            # print(upperright)

            # self.leg1 = self.ax.legend(
            # self.lhandles1,
            # self.llabels1,
            #labelspacing = .9,
            #borderaxespad = .9,
            #prop = self.fp_legend,
            #markerscale = 1.0,
            #frameon = False,
            #numpoints = 1,
            #loc = 'upper left',
            #bbox_to_anchor =list(upperright),
            #bbox_transform = self.fig.transFigure
            #)

            self.leg2.get_title().set_fontproperties(self.fp_ticklabel)

    def axes_labels(self):
        if self.xlab is not None:
            self._xlab = self.ax.set_xlabel(self.xlab)
            self.ax.xaxis.label.set_fontproperties(self.fp_axis_lab)
        if self.ylab is not None:
            self._ylab = self.ax.set_ylabel(self.ylab)
            self.ax.yaxis.label.set_fontproperties(self.fp_axis_lab)

    def remove_annotation_overlaps(self):

        return None

        if self.labels is not None:

            self.fig.savefig(self.fname)
            self.ax.figure.canvas.draw()
            steps = [0] * len(self.annots)
            for i, a2 in enumerate(self.annots):
                overlaps = False
                for z in xrange(100):
                    for a1 in self.annots[:i]:
                        print(a1, a2)
                        if a1 is None or a2 is None:

                            continue

                        if overlap(a1.get_window_extent(),
                                   a2.get_window_extent()):
                            #print('Overlapping labels: %s and %s' % (a1._text, a2._text))
                            mv = get_moves(a1.get_window_extent(),
                                           a2.get_window_extent())
                            if steps[i] % 2 == 0:
                                a2.xyann = (a2.xyann[0] + mv[0] * 1.1 *
                                            (z / 2 + 1), a2.xyann[1])
                            else:
                                a2.xyann = (a2.xyann[0],
                                            a2.xyann[1] + mv[1] * 1.1 *
                                            (z / 2 + 1))
                            steps[i] += 1
                        else:
                            #print('OK, these do not overlap: %s and %s' % (a1._text, a2._text))
                            pass
                    if not overlaps:
                        #print('No more overlaps')
                        break

    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        self.fig.subplots_adjust(top=0.92)
        self.cvs.draw()
        self.cvs.print_figure(self.pdf)
        self.pdf.close()
        self.fig.clf()


class Histogram(Plot):
    def __init__(self,
                 data,
                 labels,
                 fname,
                 font_family='Helvetica Neue LT Std',
                 font_style='normal',
                 font_weight='normal',
                 font_variant='normal',
                 font_stretch='normal',
                 xlab='',
                 ylab='',
                 title='',
                 axis_lab_size=10.0,
                 lab_angle=90,
                 lab_size=(9, 9),
                 color=None,
                 palette=None,
                 rc={},
                 context='poster',
                 figsize=(5.0, 3.0),
                 bins=None,
                 nbins=None,
                 x_log=False,
                 y_log=False,
                 tone=2,
                 alpha=0.5,
                 legend_size=6,
                 xlim=None,
                 kde_base=0.2,
                 kde_perc=12.0,
                 **kwargs):
        self.data = data
        if type(self.data[0]) in common.numeric_types:
            self.data = [data]
        for i, d in enumerate(self.data):
            if type(d) is list:
                self.data[i] = np.array(d)
        self.labels = labels
        if type(self.labels) in _const.CHAR_TYPES:
            self.labels = [labels]
        for k, v in iteritems(locals()):
            setattr(self, k, v)
        self.sns = sns
        self.rc = self.rc or {
            'lines.linewidth': 1.0,
            'patch.linewidth': 0.0,
            'grid.linewidth': 1.0
        }
        super(Histogram, self).__init__(
            fname=fname,
            font_family=font_family,
            font_style=font_style,
            font_weight=font_weight,
            font_variant=font_variant,
            font_stretch=font_stretch,
            palette=palette,
            context=context,
            lab_size=self.lab_size,
            axis_lab_size=self.axis_lab_size,
            rc=self.rc)
        if self.color is None:
            self.set_palette()
        self.data_range()
        self.set_bins(bins)
        self.plot_args = kwargs
        self.plot(**kwargs)

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def set_bins(self, bins=None, nbins=None):
        self.bins = bins
        if self.bins is None:
            self.bins_default()
        self.nbins = nbins
        if self.nbins is None:
            self.nbins_default()

    def data_range(self):
        self.lowest = min(
            map(lambda d: np.nanmin(d[np.where(d != 0.0)]), self.data))
        self.highest = max(map(max, self.data))

    def nbins_default(self):
        self.nbins = min(
            map(
                lambda d:
                    len(d) / (150.0 + 3.0 * np.log10(len(d) / np.log10(len(d)))),
                self.data
            )
        )

    def bins_default(self):
        if self.x_log:
            self.bin_limits = sorted([
                np.log10(self.lowest) if self.lowest > 0.0 else
                np.log10(0.000001), np.log10(self.highest)
            ])
        else:
            self.bin_limits = [self.lowest, self.highest]
        self.bins = \
            np.logspace(self.bin_limits[0], self.bin_limits[1], self.nbins) \
            if self.x_log else \
            np.linspace(self.lowest, self.highest, self.nbins)

    def set_palette(self, palette=None):
        self.palette = self.palette if palette is None else palette
        self.colors = \
            map(
                lambda x:
                    x if len(x) < 100 else '%s%s' % (
                        x, '%02x' % (self.alpha * 255.0)),
                map(
                    lambda i:
                        self.palette[i % len(self.palette)]
                        [min(len(self.palette[i]) - 1, self.tone)],
                    xrange(len(self.data))
                )
            )

    def remove_borders(self):
        for patches in self.hist[2]:
            map(lambda p: p.set_linewidth(0.0), patches)

    def set_labels(self):
        self.ax.set_ylabel(self.ylab)
        self.ax.set_xlabel(self.xlab)
        self.ax.set_title(self.title)

    def add_density_lines(self, **kwargs):
        for i, d in enumerate(self.data):
            self.kde_bandwidth = self.kde_base / d.std(ddof=1)
            density = stats.gaussian_kde(d, bw_method=self.kde_bandwidth)
            x = np.arange(self.lowest, self.highest,
                          self.highest / len(self.hist[0][i]))
            y = np.array(density(x))
            limit = np.percentile(x, self.kde_perc)
            y = y[np.where(x < limit)]
            x = x[np.where(x < limit)]
            #y2 = mpl.mlab.normpdf(x, np.mean(d), np.std(d))
            ylim = self.ax.get_ylim()
            xlim = self.ax.get_xlim()
            self.ax.plot(
                x,
                y,
                ls='--',
                lw=.5,
                c=self.palette[i][0],
                label='%s, density' % self.labels[i])
            #self.ax.plot(x, y2, ls = ':', lw = .5, c = self.palette[i][0])
            self.ax.set_ylim(ylim)
            self.ax.set_xlim(xlim)

    def set_log(self):
        if self.y_log:
            self.ax.set_yscale('log')
        if self.x_log:
            self.ax.set_xscale('log')

    def set_ticklabels(self):
        # self.ax.yaxis.set_ticklabels(
        # map(
        # lambda x:
        #'%.01f%%' % x if x >= 0.1 else '',
        # self.ax.get_yticks()
        #)
        #)
        self.ax.xaxis.set_ticklabels(
            map(lambda x: '{:,g}'.format(x), self.ax.get_xticks()))

    def set_xlim(self):
        if self.xlim is not None:
            self.ax.set_xlim(self.xlim)

    def add_legend(self):
        self.ax.legend(prop={'size': self.legend_size})

    def plot(self, **kwargs):
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.ax = self.fig.add_subplot(111)
        # for i, d in enumerate(self.data):
        #    sns.distplot(d, ax = self.ax,
        #                 axlabel = False, color = self.colors[i])
        self.hist = self.ax.hist(
            self.data,
            bins=self.bins,
            label=self.labels,
            color=self.colors,
            log=self.y_log,
            alpha=self.alpha,
            **kwargs)
        self.remove_borders()
        self.set_labels()
        self.set_log()
        self.add_density_lines()
        self.set_xlim()
        self.set_ticklabels()
        self.add_legend()


class SimilarityGraph(object):

    def __init__(self,
                 pp,
                 fname,
                 similarity,
                 size,
                 layout_method='fruchterman_reingold',
                 layout_param={},
                 width=1024,
                 height=1024,
                 margin=124,
                 **kwargs):

        for k, v in iteritems(locals()):
            setattr(self, k, v)

        self.graph = self.pp.graph
        self.size_param_defaults = {
            'vertex': (1.12, 0.55, 0.040),
            'edge': (1.25, 0.48, 0.065),
            # 'edge': (1.25, 0.48, 0.065),
            'curation': (0.125, 0.44, 0.08)
        }
        self.scale_defaults = {
            'vertex': {
                'vscale': [50, 100, 500, 1000, 5000],
                'escale': [0.05, 0.1, 0.2, 0.5],
                'ew': lambda x: (x * 10.0)**1.8
            },
            'edge': {
                #'vscale': [10, 50, 100],
                'vscale': [50, 100, 500, 1000, 2000],
                'escale': [0.05, 0.1, 0.2, 0.5],
                'ew': lambda x: (x * 10.0)**1.3 # this was 1.8
            },
            'curation': {
                'vscale': [100, 1000, 5000, 10000, 20000],
                'escale': [5.0, 10.0, 15.0, 30.0],
                'ew': lambda x: (x * 0.35)**1.3
            }
        }

        self.plot()

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        self.get_similarity()
        self.init_sgraph()
        self.get_size()
        self.make_layout()
        self.layout_limits()
        self.build_legend()
        self.legend_coordinates()
        self.init_pdf()

    def do_plot(self):
        self.make_plot()

    def post_plot(self):
        self.finish()

    def get_similarity(self):
        if type(self.similarity) is str and hasattr(self, '%s_sim' %
                                                    self.similarity):
            self.size_param = self.size_param_defaults[self.similarity]
            self.vscale = self.scale_defaults[self.similarity]['vscale']
            self.escale = self.scale_defaults[self.similarity]['escale']
            getattr(self, '%s_sim' % self.similarity)()
        elif hasattr(self.similarity, '__call__'):
            self.edges = self.similarity(self.pp)

    def get_size(self):
        if type(self.size) is str and hasattr(self, 'sizes_%s' % self.size):
            getattr(self, 'sizes_%s' % self.size)()
        elif hasattr(self.similarity, '__call__'):
            self.sgraph.vs['size'] = self.size(self.pp)

    def vertex_sim(self):
        self.sim = self.pp.databases_similarity()
        self.edges = [
            e
            for e in [(it[0], iit[0], iit[1], self.sim['nodes'][it[0]][iit[0]])
                      for it in self.sim['edges'].items()
                      for iit in it[1].items()]
            if (e[2] > 0.15 or (e[0] == 'MatrixDB' and e[2] > 0.02
                                ) or (e[0] == 'NRF2ome' and e[2] > 0.07
                                      ) or (e[0] == 'ACSN' and e[2] > 0.10)
                ) and e[0] != e[1]
        ]

    def edge_sim(self):
        self.edges = [(s1, s2, pypath.common.simpson_index(
            list(
                map(lambda e: e.index,
                    filter(lambda e: s1 in e['sources'], self.graph.es))),
            list(
                map(lambda e: e.index,
                    filter(lambda e: s2 in e['sources'], self.graph.es))), ))
                      for s1 in self.pp.sources for s2 in self.pp.sources]

        self.edges = [e for e in self.edges if e[2] > 0.0545 and e[0] != e[1]]

    def refs_sim(self):
        self.edges = [(s1, s2, pypath.common.simpson_index([
            r.pmid
            for r in common.unique_list(
                common.flat_list([[] if s1 not in e['refs_by_source'] else e[
                    'refs_by_source'][s1] for e in self.graph.es]))
        ], [
            r.pmid
            for r in common.unique_list(
                common.flat_list([[] if s2 not in e['refs_by_source'] else e[
                    'refs_by_source'][s2] for e in self.graph.es]))
        ])) for s1 in self.pp.sources for s2 in self.pp.sources]

        self.edges = [e for e in self.edges if e[2] > 0.0545 and e[0] != e[1]]

    def curation_sim(self):
        self.edges = [(s1, s2, (sum([
            0.0 if s1 not in e['refs_by_source'] or
            s2 not in e['refs_by_source'] else len(
                set([r1.pmid
                     for r1 in e['refs_by_source'][s1]]).symmetric_difference(
                         set([r2.pmid for r2 in e['refs_by_source'][s2]])))
            for e in self.graph.es
        ]) / float(
            len([
                e for e in self.graph.es
                if s1 in e['sources'] and s2 in e['sources']
            ]) + 0.001))) for s1 in self.pp.sources for s2 in self.pp.sources]

        self.edges = [e for e in self.edges if e[2] >= 4.0 and e[0] != e[1]]

    def init_sgraph(self):
        self.sgraph = igraph.Graph.TupleList(self.edges, edge_attrs=['weight'])
        self.sgraph.simplify(combine_edges='mean')

    def make_layout(self):
        self.layout_param_defaults = {
            'fruchterman_reingold': {
                'weights': 'weight' if self.similarity != 'curation' else None,
                'repulserad': self.sgraph.vcount()**2.8,
                'maxiter': 1000,
                'area': self.sgraph.vcount()**2.3
            }
        }
        if self.layout_method in self.layout_param_defaults:
            self.layout_param = common.merge_dicts(
                self.layout_param_defaults[self.layout_method],
                self.layout_param)
        self.layout = getattr(self.sgraph, 'layout_%s' %
                              self.layout_method)(**self.layout_param)

    def sizes_vertex(self):

        self.sgraph.vs['size'] = \
            list(
                map(
                    lambda v:
                        len(
                            list(
                                filter(
                                    lambda e:
                                        v['name'] in e['sources'],
                                    self.graph.es
                                )
                            )
                        )**self.size_param[1],
                    self.sgraph.vs
                )
        )

    def sizes_edge(self):
        """
        Sets the size according to number of edges for each resource.
        """

        self.sgraph.vs['size'] = \
            list(
                map(
                    lambda v:
                        len(
                            list(
                                filter(
                                    lambda e:
                                        v['name'] in e['sources'],
                                    self.graph.es
                                )
                            )
                        )**0.48,
                    self.sgraph.vs
                )
        )

    def sizes_refs(self):
        self.sgraph.vs['size'] = \
            [len(
                common.unique_list(
                    common.flat_list([
                        e['refs_by_source'][v['name']]
                        for e in self.graph.es
                        if v['name'] in e['refs_by_source']
                    ])
                )
            )**0.48 for v in self.sgraph.vs]

    def sizes_curation(self):
        sizes = []
        for v in self.sgraph.vs:
            allrefs = \
                len(
                    common.unique_list(
                        common.flat_list(
                            [[r.pmid for r in e1['refs_by_source'][v['name']]]
                                for e1 in self.graph.es
                                if v['name'] in e1['refs_by_source']
                             ]
                        )
                    )
                )
            alledges = float(
                len([
                    e2.index for e2 in self.graph.es
                    if v['name'] in e2['sources']
                ]))
            uniqcits = sum([
                len([rr.pmid for rr in e3['refs_by_source'][v['name']]])
                for e3 in self.graph.es if v['name'] in e3['refs_by_source']
            ])

            sizes.append((allrefs / alledges * uniqcits)**self.size_param[1])

        self.sgraph.vs['size'] = sizes

    def build_legend(self):
        self.sgraph.add_vertices([str(i) for i in self.vscale])
        self.sgraph.add_vertices(
            ['%.2f_%u' % (i, a) for i in self.escale for a in [0, 1]])
        self.sgraph.add_edges(
            [('%.2f_%u' % (i, 0), '%.2f_%u' % (i, 1)) for i in self.escale])

    def layout_limits(self):
        self.xmax = max([c[0] for c in self.layout._coords])
        self.ymin = min([c[1] for c in self.layout._coords])
        self.xrng = self.xmax - min([c[0] for c in self.layout._coords])
        self.yrng = max([c[1] for c in self.layout._coords]) - self.ymin
        self.xleg = self.xmax + self.xrng * 0.2

    def legend_coordinates(self):
        for i, s in enumerate(self.vscale):
            v = self.sgraph.vs[self.sgraph.vs['name'].index(str(s))]
            v['size'] = s**self.size_param[1]
            self.layout._coords.append([
                self.xleg,
                # start from ymin
                self.ymin \
                # plus a constant distance at each step
                # scaled by size_ param[0]


                + i * self.size_param[0] \
                # add discances dependent on the marker size
                # this is scaled by size_param[1], just like
                # the graph's vertices



                + sum(self.vscale[:i + 1])**self.size_param[1] \
                # and the ratio between the constant and
                # the increasing component set by
                # size_param[2]:



                * self.size_param[2]]
            )

        self.sgraph.es['label'] = ['' for _ in self.sgraph.es]

        for i, s in enumerate(self.escale):
            v1 = self.sgraph.vs[self.sgraph.vs['name'].index('%.2f_%u' % (s, 0
                                                                          ))]
            v2 = self.sgraph.vs[self.sgraph.vs['name'].index('%.2f_%u' % (s, 1
                                                                          ))]
            e = self.sgraph.es[self.sgraph.get_eid(v1.index, v2.index)]
            e['weight'] = s
            e['label'] = '%.2f' % s
            ycoo = self.ymin + self.yrng * 0.7 + i * 1.8
            self.layout._coords.append([self.xleg - self.xrng * 0.07, ycoo])
            self.layout._coords.append([self.xleg + self.xrng * 0.07, ycoo])
            v1['size'] = 0.0
            v2['size'] = 0.0
            v1['name'] = ''
            v2['name'] = ''

    def init_pdf(self):
        import cairo
        self.surface = cairo.PDFSurface(self.fname, self.width, self.height)
        self.bbox = igraph.drawing.utils.BoundingBox(self.margin, self.margin,
                                                     self.width - self.margin,
                                                     self.height - self.margin)

    def make_plot(self):
        self.fig = igraph.plot(
            self.sgraph,
            vertex_label=self.sgraph.vs['name'],
            layout=self.layout,
            bbox=self.bbox,
            target=self.surface,
            drawer_factory=DefaultGraphDrawerFFsupport,
            vertex_size=self.sgraph.vs['size'],
            vertex_frame_width=0,
            vertex_color='#77AADD',
            vertex_label_color='#000000FF',
            vertex_label_family='Sentinel Book',
            edge_label_color='#000000FF',
            edge_label_family='Sentinel Book',
            vertex_label_size=27,
            vertex_label_dist=1.4,
            edge_label_size=27,
            edge_label=self.sgraph.es['label'],
            edge_width=list(
                map(self.scale_defaults[self.similarity]['ew'], self.sgraph.es[
                    'weight'])),
            edge_color='#77441155',
            edge_curved=False)

    def finish(self):
        self.fig.redraw()
        self.fig.save()


class HistoryTree(object):
    def __init__(self, fname, **kwargs):

        for k, v in iteritems(locals()):
            setattr(self, k, v)

        for k, v in iteritems(kwargs):
            setattr(self, k, v)

        self.tikzfname = fname

        self.get_years()

        self.defaults = {
            'yearbarwidth': 0.4,
            'sepwidth': 0.04,
            'lineheight': [1.0] * len(self.years),
            'labelbg': 'twilightblue',
            'labelfg': 'teal',
            'nodelabbg': 'teal',
            'nodelabfg': 'white',
            'dotcol': 'teal',
            'linecol': 'teal',
            'dataimportcol': 'mantis',
            'dotsize': 3.0,
            'linewidth': 1.0,
            'rowbg': 'twilightblue',
            'width': 20.0,
            'xoffset': 0.5,
            'compile': True,
            'dotlineopacity': 1.0,
            'horizontal':
            True,  # whether the timeline should be the horizontal axis,
            'latex': '/usr/bin/xelatex'
        }

        for attr, val in iteritems(self.defaults):
            if not hasattr(self, attr):
                setattr(self, attr, val)

        self.plot()

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        self.build_elements()
        self.compose_tikz()
        self.write_tex()
        self.run_latex()

    def get_years(self):
        self.d = pypath.descriptions.descriptions
        self.firstyear = min(
            common.flat_list([
                [r['year']] for r in self.d.values() if 'year' in r
            ] + [r['releases'] for r in self.d.values() if 'releases' in r]))
        self.lastyear = self.lastyear if hasattr(
            self, 'lastyear') else date.today().year
        self.years = list(range(self.lastyear - self.firstyear + 1))

    def build_elements(self):
        self.set_styles()
        self.set_latex_header()
        self.sort()
        self.set_grid()
        self.get_timelines()
        self.get_layer_lines()
        self.get_connections()
        self.get_layer_nodes()
        self.get_legend()
        self.get_closing()

    def compose_tikz(self):
        self.tikz = \
            self.latex_header + \
            self.grid + \
            self.layer_nodes + \
            self.timelines + \
            self.legend + \
            self.layer_lines + \
            self.connections + \
            self.closing

    def set_styles(self):
        # TikZ styles
        self.tikzstyles = r'''
            \tikzstyle{omnipath}=[rectangle, anchor = center, inner sep = 2pt, fill = %s,
                rotate = 90, text = %s, draw = %s]
            \tikzstyle{others}=[rectangle, anchor = center, inner sep = 2pt, fill = %s,
                rotate = 90, text = %s, draw = %s]
        ''' % (self.nodelabfg, self.nodelabbg, self.nodelabbg, self.nodelabbg,
               self.nodelabfg, self.nodelabbg)

    def set_latex_header(self):
        # LaTeX preamble for XeLaTeX
        self.latex_header = r'''\documentclass[10pt]{extarticle}
            \usepackage{fontspec}
            \usepackage{xunicode}
            \usepackage{polyglossia}
            \setdefaultlanguage{english}
            \usepackage{xltxtra}
            \usepackage{microtype}
            \usepackage[cm]{fullpage}
            \usepackage{rotating}
            \usepackage{fullpage}
            \usepackage[usenames,dvipsnames,svgnames,table]{xcolor}
            \usepackage[
                voffset=0cm,
                hoffset=0cm,
                left=0cm,
                top=0cm,
                bottom=0cm,
                right=0cm,
                paperwidth=18cm,
                paperheight=23cm,
                footskip=0pt,
                headheight=0pt,
                headsep=0pt
                ]{geometry}
            \usepackage{color}
            \setmainfont{HelveticaNeueLTStd-LtCn}
            \usepackage{tikz}
            \definecolor{zircon}{RGB}{228, 236, 236}
            \definecolor{teal}{RGB}{17, 68, 119}
            \definecolor{twilightblue}{RGB}{227, 238, 248}
            \definecolor{mantis}{RGB}{119, 204, 204}%s
            \begin{document}
            \thispagestyle{empty}
            \pgfdeclarelayer{background}
            \pgfdeclarelayer{nodes}
            \pgfdeclarelayer{lines}
            \pgfsetlayers{background,lines,nodes}%s
            \begin{tikzpicture}
            \begin{pgfonlayer}{background}
        ''' % (self.tikzstyles, r'''
            \begin{turn}{-90}''' if self.horizontal else '')

    def sort(self):
        self.ordr = sorted(
            [(lab, r['releases'] if 'releases' in r else [] + [r['year']]
              if 'year' in r else [], r['label']
              if 'label' in r else lab, r['data_import']
              if 'data_import' in r else [], 'omnipath'
              if 'omnipath' in r and r['omnipath'] else 'others')
             for lab, r in iteritems(self.d)
             if 'year' in r or 'releases' in r],
            key=lambda x: min(x[1]))

    def set_grid(self):
        # the background grid and year labels
        self.grid = ''
        for i in self.years:
            self.grid += r'''        \fill[anchor = south west, fill = %s,
                inner sep = 0pt, outer sep = 0pt] '''\
                r'''(%f, %f) rectangle (%f, %f);
                \node[anchor = north west, rotate = 90, text width = %fcm,
                    fill = %s, inner sep = 0pt, outer sep = 0pt, align = center,
                    minimum size = %fcm] at (0.0, %f) {\small{\color{%s}%u}};
                \fill[fill = red] (%f, %f) circle (0.0pt);
        ''' % (
                self.rowbg,  # background of row
                self.yearbarwidth + self.sepwidth,  # left edge of row
                sum(self.lineheight[:i]),  # top edge of row
                self.width,  # right edge of the row
                sum(self.lineheight[:i + 1]) - \
                self.sepwidth,  # bottom edge of row
                self.lineheight[i] - self.sepwidth,  # height of year label
                self.labelbg,  # background of label
                self.yearbarwidth,  # width of year label
                sum(self.lineheight[:i]),  # top of year label
                self.labelfg,  # text color of label
                self.firstyear + i,  # year
                0.0, sum(self.lineheight[:i])  # red dot
            )

    def get_layer_nodes(self):
        # new layer for nodes
        self.layer_nodes = r'''    \end{pgfonlayer}
            \begin{pgfonlayer}{nodes}
            '''

    def get_timelines(self):
        # horizontal distance between vertical columns
        self.xdist = (
            self.width - self.yearbarwidth - self.sepwidth - self.xoffset
        ) / float(len(self.ordr))
        self.nodelabels = []
        # drawing vertical dots, labels and connecting lines:
        self.timelines = ''
        for i, r in enumerate(self.ordr):
            coox = self.xdist * i + self.yearbarwidth + \
                self.sepwidth + self.xdist / 2.0 + self.xoffset
            ymax = max(r[1])
            ydots = [y for y in r[1] if y != ymax]
            ylaby = ymax - self.firstyear
            cooylab = sum(self.lineheight[:ylaby]) + \
                self.lineheight[ylaby] / 2.0
            ydots = [
                sum(self.lineheight[:y - self.firstyear]) +
                self.lineheight[y - self.firstyear] / 2.0 for y in ydots
            ]
            for j, cooy in enumerate(ydots):
                self.timelines += r'''        \node[circle, fill = %s, minimum size = %f, opacity = %f]
                    (%s) at (%f, %f) {};
                ''' % (
                    self.dotcol,  # fill color for dot
                    self.dotsize,  # size of the dot
                    self.dotlineopacity,  # opacity of dot
                    '%s%u' % (r[0].lower(), j),  # label
                    coox,  # x coordinate
                    cooy  # y coordinate
                )
            self.timelines += r'''        \node[%s]
                    (%s) at (%f, %f)
                    {\footnotesize %s};
            ''' % (
                r[4],  # node style
                r[0].lower(),  # node name
                coox,  # node x coordinate
                cooylab,  # node y coordinate
                r[0]  # label text
            )
            self.nodelabels.append(r[0].lower())
            if len(r[1]) > 1:
                self.timelines += r'''        \draw[draw = %s, line width = %fpt, opacity = %f] (%s%s);
                ''' % (self.linecol, self.linewidth, self.dotlineopacity,
                       '%s) -- (' % r[0].lower(), ') -- ('.join([
                           '%s%u' % (r[0].lower(), j) for j in xrange(
                               len(ydots))
                       ]))

    def get_legend(self):
        # legend
        self.legend = ''
        self.legend += r'''        \node[circle, anchor = south, minimum size = %f,
                opacity = %f, fill = %s] at (%f, %f) {};
            ''' % (self.dotsize, self.dotlineopacity, self.dotcol,
                   self.width - 1.5, 0.5)
        self.legend += r'''        \node[anchor = west, rotate = 90] at (%f, %f) {\color{teal} Release/update year};
            ''' % (self.width - 1.5, 1.2)

        self.legend += r'''
                \draw[-latex, draw = %s, line width = %fpt, opacity = %f]
                (%f, %f) -- (%f, %f);
                \node[anchor = west, rotate = 90] at (%f, %f) {\color{teal} Data transfer};
            ''' % (self.dataimportcol, self.linewidth, self.dotlineopacity,
                   self.width - 0.9, 0.5, self.width - 0.9, 1.0,
                   self.width - 0.9, 1.2)

        self.legend += r'''        \node[others, anchor = west] at (%f, %f) {Other resource};
            ''' % (self.width - 2.1, 0.5)
        self.legend += r'''        \node[omnipath, anchor = west] at (%f, %f) {Resource in OmniPath};
            ''' % (self.width - 2.7, 0.5)

    def get_layer_lines(self):
        # new layer for crossing lines showing data transfers:
        self.layer_lines = r'''\end{pgfonlayer}
            \begin{pgfonlayer}{lines}
            '''

    def get_connections(self):
        # drawing data transfer lines:
        self.connections = ''
        for r in self.ordr:
            for s in r[3]:
                if r[0].lower() in self.nodelabels and s.lower(
                ) in self.nodelabels:
                    self.connections += r'''        \draw[-latex, draw = %s, line width = %fpt, opacity = %f]
                        (%s) -- (%s);
                    ''' % (self.dataimportcol, self.linewidth,
                           self.dotlineopacity, s.lower(), r[0].lower())

    def get_closing(self):
        # closing layer, tikzpicture and LaTeX document:
        self.closing = r'''    \end{pgfonlayer}
            \end{tikzpicture}%s
            \end{document}
        ''' % (r'''
            \end{turn}''' if self.horizontal else '')

    def write_tex(self):
        """
        Writes the TeX markup to file.
        """
        with open(self.tikzfname, 'w') as f:
            f.write(self.tikz)

    def run_latex(self):
        """
        Runs LaTeX to compile the TeX file.
        """
        if self.compile:
            subprocess.call([self.latex, self.tikzfname])
            sys.stdout.write(
                '\t:: TeX has been written to `%s`, PDF to `%s.pdf`.\n' %
                (self.tikzfname, '.'.join(self.tikzfname.split('.')[:-1])))


class HtpCharacteristics(object):
    def __init__(self,
                 pp,
                 fname,
                 upper=200,
                 lower=5,
                 axis_lab_font={},
                 ticklabel_font={},
                 title_font={},
                 title='',
                 htdata={},
                 **kwargs):

        for k, v in iteritems(locals()):
            setattr(self, k, v)

        for k, v in iteritems(kwargs):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.defaults = {
            'figsize': (7, 14),
            'title_halign': 'center',
            'title_valign': 'top'
        }

        for k, v in iteritems(self.defaults):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.plot_param = [('rnum', 'Number of references', '#332288'),
                           ('enum', 'Number of edges', '#44AA99'),
                           ('snum', 'Number of resources', '#117733'),
                           ('lenum', 'LT interactions', '#882255'),
                           ('lvnum', 'LT proteins', '#AA4499')]

        self.graph = pp.graph

        self.axis_lab_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'x-large'
        }
        self.ticklabel_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'large'
        }
        self.title_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'xx-large'
        }

        self.axis_lab_font = common.merge_dicts(axis_lab_font,
                                                self.axis_lab_font_default)
        self.ticklabel_font = common.merge_dicts(ticklabel_font,
                                                 self.ticklabel_font_default)
        self.title_font = common.merge_dicts(title_font,
                                             self.title_font_default)

        self.plot()

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
        if hasattr(self, 'pp'):
            self.pp.reload()

    def plot(self):
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        self.set_fontproperties()
        self.htp_calculations()
        self.set_figsize()

    def do_plot(self):
        self.init_fig()
        self.set_grid()
        self.make_plots()
        self.set_ticklabels()
        self.set_title()

    def post_plot(self):
        self.finish()
        if hasattr(self, 'prg'):
            self.prg.terminate()

    def set_fontproperties(self):
        self.fp_axis_lab = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.axis_lab_font))
        self.fp_ticklabel = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.ticklabel_font))
        self.fp_title = \
            mpl.font_manager.FontProperties(**copy.deepcopy(self.title_font))

    def reset(self):
        self.htdata = {}

    def get_point(self, htlim):

        if htlim not in self.htdata:
            htrefs = set(
                [i[0] for i in self.refc.most_common() if i[1] > htlim])
            htedgs = set([
                e.index for e in self.graph.es
                if len(set([r.pmid for r in e['references']]) - htrefs) == 0
            ])
            htvtcs = [
                v.index for v in self.graph.vs
                if len(
                    set(self.graph.incident(
                        v.index, mode='ALL')) - htedgs) == 0
            ]
            htsrcs = common.unique_list(
                common.flat_list([self.graph.es[e]['sources'] for e in htedgs]))
            self.htdata[htlim] = {
                'rnum': len(htrefs),
                'enum': len(htedgs),
                'snum': len(htsrcs),
                'htsrcs': set(htsrcs),
                'htrefs': htrefs,
                'lenum': self.graph.ecount() - len(htedgs),
                'lvnum': self.graph.vcount() - len(htvtcs)
            }

    def htp_calculations(self):
        if not len(self.htdata):
            self.refc = collections.Counter(
                common.flat_list((r.pmid for r in e['references'])
                                for e in self.pp.graph.es))

            # percentage of high throughput interactions
            htsrcs_prev = set(self.pp.sources)
            self.prg = progress.Progress(self.upper - self.lower,
                                         'Analysing HTP refs/interactions', 1)
            for htlim in reversed(xrange(self.lower, self.upper + 1)):
                self.prg.step()
                self.get_point(htlim)
                htsrcs_new = self.htdata[htlim]['htsrcs']
                diff = htsrcs_new - htsrcs_prev
                htsrcs_prev = htsrcs_new
                if len(diff):
                    sys.stdout.write(
                        '\n\t:: %s: no references with more interaction than %u\n'
                        % (', '.join(list(diff)), htlim - 1))
                    sys.stdout.flush()

    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        if hasattr(self, 'width') and hasattr(self, 'height'):
            self.figsize = (self.width, self.height)

    def init_fig(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """
        self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)
        self.axes = {}

    def set_grid(self):
        """
        Sets up a grid according to the number of subplots,
        with one additional column of zero width on the left
        to have aligned y axis labels.
        """
        self.nrows = len(self.plot_param)
        self.gs = mpl.gridspec.GridSpec(
            self.nrows, 2, height_ratios=[1] * self.nrows,
            width_ratios=[0, 1])
        self.axes = [[None] * self.nrows, [None] * self.nrows]

    def get_subplot(self, i, j=1):
        if self.axes[j][i] is None:
            self.axes[j][i] = self.fig.add_subplot(self.gs[i, j])
        self.ax = self.axes[j][i]

    def make_plots(self):
        if hasattr(self, 'prg'):
            self.prg.set_status('drawing plots')
        for i, (key, ylab, col) in enumerate(self.plot_param):
            self.get_subplot(i)
            self.ax.plot(
                sorted(self.htdata.keys()),
                list(
                    map(lambda h: self.htdata[h][key],
                        sorted(self.htdata.keys()))),
                '-',
                color=col,
                linewidth=2)

            if i == self.nrows - 1:
                # xlabel only for the lowest subplot
                self.ax.set_xlabel(
                    'HT limit [interaction/reference]',
                    fontproperties=self.fp_axis_lab)

            if i == 0:
                self.ax.set_yscale('symlog')
                self.ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

            self.ax.set_xscale('log')
            self.ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
            self.ax.yaxis.grid(True, color='#FFFFFF', lw=2, ls='solid')
            self.ax.xaxis.grid(False)
            self.ax.set_axisbelow(True)
            self.ax.set_facecolor('#EAEAF2')
            list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
            self.ax.tick_params(which='both', length=0)

            self.get_subplot(i, 0)
            self.ax.set_ylabel(ylab, fontproperties=self.fp_axis_lab)
            self.ax.xaxis.set_ticklabels([])
            self.ax.yaxis.set_ticklabels([])
            self.ax.xaxis.label.set_verticalalignment('bottom')
            list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
            self.ax.tick_params(which='both', length=0)

    def set_ticklabels(self):
        for ax in self.axes[1]:
            list(
                map(lambda l: l.set_fontproperties(self.fp_ticklabel),
                    ax.xaxis.get_majorticklabels()))
            list(
                map(lambda l: l.set_fontproperties(self.fp_ticklabel),
                    ax.yaxis.get_majorticklabels()))
        for ax in self.axes[0]:
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.label.set_verticalalignment('bottom')

    def set_title(self):
        """
        Sets the main title.
        """
        self.title_text = self.fig.suptitle(self.title)
        self.title_text.set_fontproperties(self.fp_title)
        self.title_text.set_horizontalalignment(self.title_halign)
        self.title_text.set_verticalalignment(self.title_valign)

    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        self.fig.subplots_adjust(top=0.94)
        self.cvs.draw()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.cvs.print_figure(self.pdf)

        self.pdf.close()
        self.fig.clf()


class RefsComposite(object):

    def __init__(self,
                 pp,
                 fname,
                 axis_lab_font={},
                 ticklabel_font={},
                 title_font={},
                 legend_font={},
                 bar_args={},
                 title='',
                 color='#88CCEE',
                 hcolor='#88CCEE',
                 all_color='#44AA99',
                 htp_threshold=20,
                 figsize=(12.8, 8.8),
                 all_name='All',
                 curation_plot=False,
                 pubmeds=None,
                 earliest=None,
                 **kwargs):

        for k, v in iteritems(locals()):
            if not hasattr(self, k):
                setattr(self, k, v)

        for k, v in iteritems(kwargs):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.defaults = {'title_halign': 'center', 'title_valign': 'top'}

        for k, v in iteritems(self.defaults):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.graph = self.pp.graph

        self.boxplot_whiskerprops = {
            'linewidth': 0.5,
            'color': '#777777',
            'linestyle': '-'
        }
        self.boxplot_medianprops = {'linewidth': 0.5, 'color': '#777777'}
        self.boxplot_capprops = {'linewidth': 0.5, 'color': '#777777'}
        self.boxplot_flierprops = {
            'markersize': 3,
            'marker': 'D',
            'markerfacecolor': '#777777',
            'markeredgecolor': 'none'
        }
        self.bar_args_default = {
            'align': 'center',
            'edgecolor': 'none',
            'linewidth': 0.0,
        }

        self.axis_lab_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'x-large'
        }
        self.ticklabel_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'medium'
        }
        self.title_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'xx-large'
        }
        self.legend_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'small'
        }

        self.axis_lab_font = common.merge_dicts(axis_lab_font,
                                                self.axis_lab_font_default)
        self.ticklabel_font = common.merge_dicts(ticklabel_font,
                                                 self.ticklabel_font_default)
        self.title_font = common.merge_dicts(title_font,
                                             self.title_font_default)
        self.legend_font = common.merge_dicts(legend_font,
                                              self.legend_font_default)

        self.bar_args = common.merge_dicts(self.bar_args,
                                           self.bar_args_default)

        self.plot()

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        self.set_fontproperties()
        self.set_figsize()
        self.get_data()

    def do_plot(self):
        self.init_fig()
        self.set_grid()
        self.plot_npubmeds_db()
        self.plot_npubmeds_y()
        self.plot_newint()
        self.plot_boxplot()
        self.plot_curation()
        self.set_title()

    def post_plot(self):
        self.finish()

    def set_fontproperties(self):
        self.fp_axis_lab = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.axis_lab_font))
        self.fp_ticklabel = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.ticklabel_font))
        self.fp_title = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.title_font))
        self.fp_legend = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.legend_font))

    def get_data(self):
        if self.pubmeds is None or self.earliest is None:
            self.pubmeds, self.earliest = _refs.get_pubmed_data(
                self.pp, htp_threshold=self.htp_threshold)

    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        if hasattr(self, 'width') and hasattr(self, 'height'):
            self.figsize = (self.width, self.height)

    def init_fig(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """
        self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)
        self.axes = {}

    def set_grid(self):
        """
        Sets up a grid according to the number of subplots,
        with one additional column of zero width on the left
        to have aligned y axis labels.
        """
        self.gs = mpl.gridspec.GridSpec(
            2, 3, height_ratios=[2, 8], width_ratios=[1.1, 1.1, 8])
        self.axes = [[None, None, None], [None, None, None]]

    def get_subplot(self, i, j):
        if self.axes[i][j] is None:
            self.axes[i][j] = self.fig.add_subplot(self.gs[i, j])
        self.ax = self.axes[i][j]

    def get_hcolor(self, names):
        for name in names:
            if type(self.hcolor) is dict and name in self.hcolor:
                yield self.hcolor[name]
            elif name == self.all_name and self.all_color is not None:
                yield self.all_color
            else:
                yield self.hcolor

    def plot_npubmeds_db(self):

        self.refc_by_db = dict(self.pubmeds.database.value_counts())

        # all means all refs in the array
        self.refc_by_db[self.all_name] = self.pubmeds.shape[0]
        self.medpubyr = \
            [(s, np.median(self.pubmeds[(self.pubmeds.database == s)]['year']),
              np.mean(self.pubmeds[(self.pubmeds.database == s)]['year']))
                for s in self.pubmeds.database.unique()] \
            + [(self.all_name, np.median(self.pubmeds.year), np.mean(self.pubmeds.year))]

        # order by median of publications years
        self.medpubyr = sorted(self.medpubyr, key=lambda x: (x[1], x[2]))
        self.boxplot_ordr = list(map(lambda x: x[0], self.medpubyr))
        self.refc_db = list(
            map(lambda s: self.refc_by_db[s], self.boxplot_ordr))

        self.get_subplot(0, 2)
        self.ax.bar(np.arange(len(self.refc_db)) + 0.5,
                    self.refc_db,
                    width=0.5,
                    color=list(self.get_hcolor(self.boxplot_ordr)),
                    **copy.deepcopy(self.bar_args))
        self.ax.set_xticks(np.arange(len(self.boxplot_ordr)) + 0.56)
        self.ax.set_xticklabels(
            self.boxplot_ordr, rotation=90, fontproperties=self.fp_axis_lab)
        self.ax.set_yscale('symlog')
        self.ax.set_xlim([0.0, len(self.boxplot_ordr)])
        self.ax.set_xlabel(
            'Number of PubMed IDs', fontproperties=self.fp_axis_lab)
        self.ax.xaxis.set_label_position('top')
        self.ax.xaxis.grid(False)
        self.ax.set_ylim([
            1, (100000 if max(self.refc_db) < 100000 else max(self.refc_db))
        ])
        self.ax.set_yticks([1, 10, 100, 1000, 10000] + ([100000] if max(
            self.refc_db) > 100000 else []))
        self.ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

        _ = list(
            map(lambda tl: tl.set_fontproperties(self.fp_ticklabel) or tl.set_fontsize('small'),
                self.ax.get_yticklabels()))
        _ = list(
            map(lambda tl: tl.set_fontproperties(self.fp_ticklabel),
                self.ax.get_xticklabels()))

        self.ax.yaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
        self.ax.xaxis.grid(False)
        self.ax.set_axisbelow(True)
        self.ax.set_facecolor('#EAEAF2')
        list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
        self.ax.tick_params(which='both', length=0)

    def plot_npubmeds_y(self):

        # the refs by year barplot
        self.get_subplot(1, 1)

        self.refc_by_y = self.pubmeds.year.value_counts()
        self.refc_y = \
            list(
                map(
                    lambda y:
                        self.refc_by_y[y] if y in self.refc_by_y else 0.0,
                    np.arange(min(self.refc_by_y.index),
                              max(self.refc_by_y.index) + 1)
                )
            )

        self.ax.barh(
            np.arange(len(self.refc_y)) + 0.7,
            self.refc_y,
            tick_label=np.arange(
                min(self.refc_by_y.index), max(self.refc_by_y.index) + 1),
            color=[self.color] * len(self.refc_y),
            **copy.deepcopy(self.bar_args))
        self.ax.set_xscale('log')
        self.ax.set_ylim([-0.3, len(self.refc_y) + 0.7])
        self.ax.set_xlim([max(self.refc_y), 1])
        self.ax.set_yticklabels([])
        self.ax.set_ylabel(
            'Number of PubMed IDs', fontproperties=self.fp_axis_lab)
        self.ax.yaxis.set_label_position('left')
        self.ax.set_xticks([1, 100, 1000] + ([10000] if max(self.refc_y) >
                                             10000 else []))
        self.ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        list(
            map(lambda tl: tl.set_fontproperties(self.fp_ticklabel) or tl.set_rotation(90),
                self.ax.get_xticklabels()))

        self.ax.xaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
        self.ax.yaxis.grid(False)
        self.ax.set_axisbelow(True)
        self.ax.set_facecolor('#EAEAF2')
        list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
        self.ax.tick_params(which='both', length=0)

    def plot_newint(self):

        # refs adding new interaction by year barplot [%]
        self.get_subplot(1, 0)

        self.newr_by_y = self.earliest.year.value_counts()
        self.refc_y_n = \
            list(
                map(
                    lambda y:
                        self.newr_by_y[y] / float(self.refc_by_y[y]) * 100.0
                    if y in self.newr_by_y
                    else 100.0 if y > min(self.refc_by_y.index)
                    else 0.0,
                    np.arange(min(self.refc_by_y.index),
                              max(self.refc_by_y.index) + 1)
                )
            )
        self.ax.barh(
            np.arange(len(self.refc_y_n)) + 0.7,
            self.refc_y_n,
            color=[self.color] * len(self.refc_y_n),
            **copy.deepcopy(self.bar_args))
        self.ax.set_ylim([-0.3, len(self.refc_y_n) + 0.7])
        self.ax.set_yticklabels([])
        self.ax.set_ylabel(
            '%s of PubMed IDs\nadding new interactions' % r'%',
            fontproperties=self.fp_axis_lab)
        self.ax.yaxis.set_label_position('left')
        self.ax.set_xlim(self.ax.get_xlim()[::-1])
        self.ax.set_xticks([0, 50, 100])
        list(
            map(lambda tl: tl.set_fontproperties(self.fp_ticklabel) or tl.set_rotation(90),
                self.ax.get_xticklabels()))

        self.ax.xaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
        self.ax.yaxis.grid(False)
        self.ax.set_axisbelow(True)
        self.ax.set_facecolor('#EAEAF2')
        list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
        self.ax.tick_params(which='both', length=0)

    def plot_boxplot(self):

        # here the boxplot
        self.get_subplot(1, 2)
        data = \
            list(
                map(
                    lambda s:
                        list(self.pubmeds.year[self.pubmeds.database == s])
                    if s in self.pubmeds.database.values else
                    list(self.pubmeds.year) if s == self.all_name else
                    [],
                    self.boxplot_ordr
                )
            )
        self.boxplot = self.ax.boxplot(
            data,
            patch_artist=True,
            whiskerprops=self.boxplot_whiskerprops,
            capprops=self.boxplot_capprops,
            medianprops=self.boxplot_medianprops,
            flierprops=self.boxplot_flierprops)

        for patch, col in zip(self.boxplot['boxes'],
                              list(self.get_hcolor(self.boxplot_ordr))):
            patch.set_facecolor(col)
            patch.set_linewidth(0.0)

        self.ax.set_xlabel('Resources', fontproperties=self.fp_axis_lab)
        self.ax.xaxis.set_ticklabels([])
        self.ax.set_ylabel('Years', fontproperties=self.fp_axis_lab)
        self.ax.yaxis.set_label_position('right')
        self.ax.set_ylim(min(self.pubmeds.year), max(self.pubmeds.year))
        self.ax.xaxis.grid(False)
        self.ax.set_ylim(
            [self.ax.get_ylim()[0] - 1, self.ax.get_ylim()[1] + 1])
        list(
            map(lambda tl: tl.set_fontproperties(self.fp_ticklabel),
                self.ax.get_yticklabels()))

        self.ax.yaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
        self.ax.xaxis.grid(False)
        self.ax.set_axisbelow(True)
        self.ax.set_facecolor('#EAEAF2')
        list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
        self.ax.tick_params(which='both', length=0)

    def plot_curation(self):

        if not self.curation_plot:
            return None

        self.numof_curated = dict((y, 0) for y in np.unique(self.pubmeds.year))
        self.all_by_y = dict((y, 0) for y in np.unique(self.pubmeds.year))

        for y in np.unique(self.pubmeds.year):
            ydata = self.pubmeds[self.pubmeds.year == y]
            for e in np.unique(ydata.eid):
                edata = ydata[ydata.eid == e]
                if edata.shape[1] > 1:
                    ndb = len(np.unique(edata.database))
                    self.all_by_y[y] += ((ndb - 1)**2) / 2.0
                    for pm in np.unique(edata.pmid):
                        pmdata = edata[edata.pmid == pm]
                        if pmdata.shape[1] > 1:
                            ndb = len(np.unique(pmdata.database))
                            self.numof_curated[y] += ((ndb - 1)**2) / 2.0

        self.ratio_curated = dict((y, self.numof_curated[y] / total)
                                  if total > 0.0 else (y, 0.0)
                                  for y, total in iteritems(self.all_by_y))

        self.first_y_edges = dict(self.pubmeds.groupby(['eid']).year.min())

        self.ecount_y = dict(
            (y, len([_ for fy in self.first_y_edges.values() if fy <= y]))
            for y in np.unique(self.pubmeds.year))

        self.vcount_y = dict((y, len(
            set([
                p for e, fy in iteritems(self.first_y_edges)
                if fy <= y
                for p in [self.graph.es[e].source, self.graph.es[e].target]
            ]))) for y in np.unique(self.pubmeds.year))

        self.ccount_by_y = dict(
            (y, self.pubmeds[self.pubmeds.year == y].shape[0])
            for y in np.unique(self.pubmeds.year))

        self.ax = self.fig.add_subplot(self.gs[0, :2])
        self.axes[0][0] = self.ax
        self.axes[0][1] = self.ax

        self.ax.plot(
            np.unique(self.pubmeds.year),
            [self.ecount_y[y] for y in np.unique(self.pubmeds.year)],
            marker='.',
            color='#6EA945',
            label='Interactions',
            alpha=0.4)
        self.ax.plot(
            np.unique(self.pubmeds.year),
            [self.vcount_y[y] for y in np.unique(self.pubmeds.year)],
            marker='.',
            color='#007B7F',
            label='Proteins',
            alpha=0.4)
        self.ax.plot(
            np.unique(self.pubmeds.year),
            [self.ccount_by_y[y] for y in np.unique(self.pubmeds.year)],
            marker='.',
            color='#FCCC06',
            label='Curation effort',
            alpha=0.4)

        self.leg = self.ax.legend(loc=2)
        list(
            map(lambda t: t.set_fontproperties(self.fp_legend),
                self.leg.get_texts()))

        self.ax.set_xlabel('Year', fontproperties=self.fp_axis_lab)
        self.ax.set_xlim([min(self.pubmeds.year), 2014])
        self.ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(5000))
        self.ax.xaxis.grid(False)
        list(
            map(lambda tl: tl.set_rotation(90) or tl.set_fontproperties(self.fp_ticklabel),
                self.ax.get_xticklabels()))
        list(
            map(lambda tl: tl.set_fontproperties(self.fp_ticklabel),
                self.ax.get_yticklabels()))

        self.ax.yaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
        self.ax.xaxis.grid(False)
        self.ax.set_axisbelow(True)
        self.ax.set_facecolor('#EAEAF2')
        list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
        self.ax.tick_params(which='both', length=0)

    def set_title(self):
        """
        Sets the main title.
        """
        self.title_text = self.fig.suptitle(self.title)
        self.title_text.set_fontproperties(self.fp_title)
        self.title_text.set_horizontalalignment(self.title_halign)
        self.title_text.set_verticalalignment(self.title_valign)

    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        self.fig.subplots_adjust(top=0.92)
        self.cvs.draw()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.cvs.print_figure(self.pdf)

        self.pdf.close()
        self.fig.clf()


class CurationPlot(object):
    def __init__(self,
                 pp,
                 fname,
                 colors,
                 pubmeds=None,
                 earliest=None,
                 axis_lab_font={},
                 ticklabel_font={},
                 legend_font={},
                 **kwargs):

        for k, v in iteritems(locals()):
            if not hasattr(self, k):
                setattr(self, k, v)

        for k, v in iteritems(kwargs):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.defaults = {
            'title_halign': 'center',
            'title_valign': 'top',
            'figsize': (4, 3)
        }

        for k, v in iteritems(self.defaults):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.graph = self.pp.graph

        self.axis_lab_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'x-large'
        }
        self.ticklabel_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'medium'
        }
        self.legend_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'small'
        }

        self.axis_lab_font = common.merge_dicts(axis_lab_font,
                                                self.axis_lab_font_default)
        self.ticklabel_font = common.merge_dicts(ticklabel_font,
                                                 self.ticklabel_font_default)
        self.legend_font = common.merge_dicts(legend_font,
                                              self.legend_font_default)

        self.plot()

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        self.set_fontproperties()
        self.get_data()
        self.set_figsize()

    def do_plot(self):
        self.init_fig()
        self.make_plot()

    def post_plot(self):
        self.finish()

    def set_fontproperties(self):
        self.fp_axis_lab = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.axis_lab_font))
        self.fp_ticklabel = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.ticklabel_font))
        self.fp_legend = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.legend_font))

    def get_data(self):
        if self.pubmeds is None:
            self.pubmeds, self.earliest = _refs.get_pubmed_data(
                self.pp, htp_threshold=self.htp_threshold)

    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        if hasattr(self, 'width') and hasattr(self, 'height'):
            self.figsize = (self.width, self.height)

    def init_fig(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """
        self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)

    def make_plot(self):

        self.numof_curated = dict((y, 0) for y in np.unique(self.pubmeds.year))
        self.all_by_y = dict((y, 0) for y in np.unique(self.pubmeds.year))

        self.first_y_edges = dict(self.pubmeds.groupby(['eid']).year.min())

        self.ecount_y = dict(
            (y, len([fy for fy in self.first_y_edges.values() if fy <= y]))
            for y in np.unique(self.pubmeds.year))

        self.vcount_y = dict((y, len(
            set([
                p for e, fy in iteritems(self.first_y_edges)
                if fy <= y
                for p in [self.graph.es[e].source, self.graph.es[e].target]
            ]))) for y in np.unique(self.pubmeds.year))

        self.ccount_by_y = dict(
            (y, self.pubmeds[self.pubmeds.year == y].shape[0])
            for y in np.unique(self.pubmeds.year))

        self.ax = self.fig.add_subplot(1, 1, 1)

        self.ax.plot(
            np.unique(self.pubmeds.year),
            [self.ecount_y[y] for y in np.unique(self.pubmeds.year)],
            marker='.',
            color=self.colors[0],
            label='Interactions',
            alpha=1.0,
            markeredgecolor='none',
            markersize=10)
        self.ax.plot(
            np.unique(self.pubmeds.year),
            [self.vcount_y[y] for y in np.unique(self.pubmeds.year)],
            marker='.',
            color=self.colors[1],
            label='Proteins',
            alpha=1.0,
            markeredgecolor='none',
            markersize=10)
        self.ax.plot(
            np.unique(self.pubmeds.year),
            [self.ccount_by_y[y] for y in np.unique(self.pubmeds.year)],
            marker='.',
            color=self.colors[2],
            label='Curated papers',
            alpha=1.0,
            markeredgecolor='none',
            markersize=10)

        self.leg = self.ax.legend(loc=2, frameon=False)
        list(
            map(lambda t: t.set_fontproperties(self.fp_legend),
                self.leg.get_texts()))

        self.ax.set_xlabel('Year', fontproperties=self.fp_axis_lab)
        #self.ax.set_ylabel('Number of proteins\nor interactions', fontproperties = self.fp_axis_lab)
        self.ax.set_xlim([1980, 2019])
        ymax = max(self.ecount_y.values())
        self.ax.yaxis.set_major_locator(
            mpl.ticker.MultipleLocator(20000 if ymax > 80000 else 10000
                                       if ymax > 50000 else 5000))
        self.ax.xaxis.grid(False)
        list(
            map(lambda tl: tl.set_rotation(90) or tl.set_fontproperties(self.fp_ticklabel),
                self.ax.get_xticklabels()))
        list(
            map(lambda tl: tl.set_fontproperties(self.fp_ticklabel),
                self.ax.get_yticklabels()))

        self.ax.yaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
        self.ax.xaxis.grid(False)
        self.ax.set_axisbelow(True)
        self.ax.set_facecolor('#EAEAF2')
        list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
        self.ax.tick_params(which='both', length=0)

    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        self.cvs.draw()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.cvs.print_figure(self.pdf)

        self.pdf.close()
        self.fig.clf()


class BarplotsGrid(object):

    def __init__(self,
                 pp,
                 x,
                 by,
                 fname,
                 ylab,
                 data=None,
                 color='#77AADD',
                 xlim=None,
                 uniform_xlim=True,
                 uniform_ylim=False,
                 full_range_x=True,
                 sort=False,
                 desc=False,
                 ylog=False,
                 hoffset=0.0,
                 woffset=0.0,
                 axis_lab_font={},
                 ticklabel_font={},
                 small_ticklabel_font={},
                 title_font={},
                 bar_args={},
                 htp_threshold=20,
                 xmin=None,
                 **kwargs):

        for k, v in iteritems(locals()):
            if not hasattr(self, k):
                setattr(self, k, v)

        for k, v in iteritems(kwargs):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.defaults = {
            'title_halign': 'center',
            'title_valign': 'top',
            'small_xticklabels': False
        }

        for k, v in iteritems(self.defaults):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.bar_args_default = {
            'align': 'center',
            'edgecolor': 'none',
            'linewidth': 0.0,
        }

        self.axis_lab_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'x-large'
        }
        self.ticklabel_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'large'
        }
        self.small_ticklabel_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'medium'
        }
        self.title_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'xx-large'
        }

        self.bar_args = common.merge_dicts(self.bar_args,
                                           self.bar_args_default)

        self.axis_lab_font = common.merge_dicts(axis_lab_font,
                                                self.axis_lab_font_default)
        self.ticklabel_font = common.merge_dicts(ticklabel_font,
                                                 self.ticklabel_font_default)
        self.small_ticklabel_font = common.merge_dicts(
            small_ticklabel_font, self.small_ticklabel_font_default)
        self.title_font = common.merge_dicts(title_font,
                                             self.title_font_default)

        self.x_by_plots = []
        self.y_by_plots = []
        self.level_by_plots = []

        self.plot()

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        self.set_fontproperties()
        self.get_data()
        self.number_of_plots()
        self.set_figsize()
        self.set_xlim()

    def do_plot(self):
        self.init_fig()
        self.set_grid()
        self.make_plots()
        self.set_title()
        self.set_ylim()

    def post_plot(self):
        self.finish()

    def set_fontproperties(self):
        self.fp_axis_lab = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.axis_lab_font))
        self.fp_ticklabel = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.ticklabel_font))
        self.fp_small_ticklabel = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.small_ticklabel_font))
        self.fp_title = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.title_font))

    def get_data(self):
        if self.data is None:
            self.data, none = _refs.get_pubmed_data(
                self.pp, htp_threshold=self.htp_threshold)
        self.levels = np.array(sorted(getattr(self.data, self.by).unique()))
        self.ilevels = np.array(
            list(
                map(lambda l: getattr(self.data, self.by) == l, self.levels)))

    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        if hasattr(self, 'width') and hasattr(self, 'height'):
            self.figsize = (self.width, self.height)

    def number_of_plots(self):
        self.numof_plots = len(getattr(self.data, self.by).unique())
        self.nrows = 1
        self.ncols = 1
        while True:
            if self.nrows * self.ncols >= self.numof_plots:
                break
            self.nrows += 1
            if self.nrows * self.ncols >= self.numof_plots:
                break
            self.ncols += 1
        self.figsize = (8 + self.woffset + self.ncols * 2,
                        8 + self.hoffset + self.nrows * 2)

    def init_fig(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """
        self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)
        self.axes = {}

    def set_grid(self):
        """
        Sets up a grid according to the number of subplots,
        with one additional column of zero width on the left
        to have aligned y axis labels.
        """
        self.gs = mpl.gridspec.GridSpec(
            self.nrows,
            self.ncols + 1,
            height_ratios=[1.0] * self.nrows,
            width_ratios=[0.0] + [1.0] * self.ncols)
        self.axes = list(
            map(lambda _: [None] * (self.ncols + 1), xrange(self.nrows)))

    def get_subplot(self, i, j):
        if self.axes[i][j] is None:
            self.axes[i][j] = self.fig.add_subplot(self.gs[i, j])
        self.ax = self.axes[i][j]

    def get_color(self, x):
        if type(self.color) is dict:
            return np.array(list(map(lambda l: self.color[l], x)))
        elif type(self.color) is list:
            return self.color[:len(x)]
        else:
            return self.color

    def set_xlim(self):
        if self.xlim is None and self.uniform_xlim:
            allx = getattr(self.data, self.x).unique()
            xmin = allx.min() if self.xmin is None else self.xmin
            xlen = allx.max() - xmin if self.full_range_x \
                else len(list(filter(lambda x: x >= xmin, allx)))
            self.xlim = [-1.0, xlen + 1.0]

    def make_plots(self):

        for row in xrange(self.nrows):
            for col in xrange(self.ncols):
                nplot = row * self.ncols + col
                if nplot >= len(self.levels):
                    continue
                self.get_subplot(row, col + 1)
                level = self.levels[nplot]
                if self.uniform_xlim:
                    x = np.array(sorted(getattr(self.data, self.x).unique()))
                else:
                    x = np.array(
                        sorted(
                            getattr(self.data, self.x)[np.where(self.ilevels[
                                nplot])[0]].unique()))
                if self.full_range_x and x.dtype == np.int:
                    xmin = self.xmin if self.xmin is not None else min(x)
                    x = np.arange(xmin, max(x) + 1)
                y = dict(
                    getattr(self.data, self.x)[self.ilevels[nplot]]
                    .value_counts())
                y = np.array(list(map(lambda l: y[l] if l in y else 0.0, x)))
                if self.sort:
                    ordr = y.argsort()
                    if self.desc:
                        ordr = ordr[::-1]
                    y = y[ordr]
                    x = x[ordr]
                self.ax.bar(x = np.arange(len(x)),
                            height = y,
                            tick_label = x,
                            color=self.get_color(x),
                            **copy.deepcopy(self.bar_args))
                self.level_by_plots.append(level)
                self.x_by_plots.append(x)
                self.y_by_plots.append(y)
                if self.uniform_xlim or self.xlim is not None:
                    self.ax.set_xlim(self.xlim)
                if self.full_range_x:
                    step = 10 if max(x) - min(x) > 50 else 5
                    xticks = list(
                        filter(lambda i: i[0] % step == 0,
                               zip(x, np.arange(len(x)))))
                    self.ax.set_xticks(list(map(lambda i: i[1], xticks)))
                    self.ax.set_xticklabels(list(map(lambda i: i[0], xticks)))
                if self.xmin is not None:
                    self.ax.set_xlim(
                        [list(x).index(self.xmin), self.ax.get_xlim()[1]])
                list(
                    map(lambda tl:
                        tl.set_fontproperties(
                            self.fp_small_ticklabel
                                if self.small_xticklabels else
                            self.fp_ticklabel
                        ) or tl.set_rotation(90),
                        self.ax.get_xticklabels()
                    )
                )
                list(
                    map(lambda tl: tl.set_fontproperties(self.fp_ticklabel),
                        self.ax.get_yticklabels()))

                self.ax.xaxis.grid(True)
                self.ax.yaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
                self.ax.xaxis.grid(False)
                self.ax.set_axisbelow(True)
                self.ax.set_facecolor('#EAEAF2')
                list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
                self.ax.tick_params(which='both', length=0)
                self.ax.set_title(level, fontproperties=self.fp_axis_lab)

                # this adds huge whitespace because matplotlib's
                # tight_layout does not explicitely consider ticklabels
                # this for the aligned xlabels
                #self.get_subplot(row * 2 + 1, col + 1)
                # self.ax.xaxis.set_ticklabels([])
                # self.ax.yaxis.set_ticklabels([])
                #self.ax.set_xlabel(level, fontproperties = self.fp_axis_lab)
                # self.ax.xaxis.label.set_verticalalignment('bottom')
                #list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
                #self.ax.tick_params(which = 'both', length = 0)
                # self.ax.set_facecolor('#CCCCCC')

                if col == 0:
                    self.get_subplot(row, 0)
                    self.ax.yaxis.set_ticklabels([])
                    self.ax.xaxis.set_ticklabels([])
                    self.ax.set_ylabel(
                        self.ylab, fontproperties=self.fp_axis_lab)
                    self.ax.yaxis.label.set_verticalalignment('top')
                    list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
                    self.ax.tick_params(which='both', length=0)

    def set_title(self):
        """
        Sets the main title.
        """
        self.title_text = self.fig.suptitle(self.title)
        self.title_text.set_fontproperties(self.fp_title)
        self.title_text.set_horizontalalignment(self.title_halign)
        self.title_text.set_verticalalignment(self.title_valign)

    def set_ylim(self):
        if self.uniform_ylim:
            ymin = min(map(lambda row:
                           min(map(lambda ax: ax.get_ylim()[0],
                                   filter(lambda ax: ax is not None, row[1:]))),
                           self.axes))
            ymax = max(map(lambda row:
                           max(map(lambda ax: ax.get_ylim()[1],
                                   filter(lambda ax: ax is not None, row[1:]))),
                           self.axes[1:]))
            for aa in self.axes:
                for ax in aa[1:]:
                    if ax is not None:
                        ax.set_ylim([ymin, ymax])

    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """
        #self.gs.update(wspace=0.1, hspace=0.1)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        self.fig.subplots_adjust(top=0.94)
        self.cvs.draw()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.cvs.print_figure(self.pdf)

        self.pdf.close()
        self.fig.clf()


class Dendrogram(object):
    def __init__(self,
                 fname,
                 data,
                 color='#4477AA',
                 axis_lab_font={},
                 ticklabel_font={},
                 **kwargs):

        for k, v in iteritems(locals()):
            if not hasattr(self, k):
                setattr(self, k, v)

        for k, v in iteritems(kwargs):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.defaults = {
            'title_halign': 'center',
            'title_valign': 'top',
            'figsize': (6, 4)
        }

        for k, v in iteritems(self.defaults):
            if not hasattr(self, k):
                setattr(self, k, v)

        self.axis_lab_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'bold',
            'variant': 'normal',
            'size': 'x-large'
        }
        self.ticklabel_font_default = {
            'family': ['Helvetica Neue LT Std'],
            'style': 'normal',
            'stretch': 'condensed',
            'weight': 'roman',
            'variant': 'normal',
            'size': 'medium'
        }

        self.axis_lab_font = common.merge_dicts(axis_lab_font,
                                                self.axis_lab_font_default)
        self.ticklabel_font = common.merge_dicts(ticklabel_font,
                                                 self.ticklabel_font_default)

        self.plot()

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def plot(self):
        self.pre_plot()
        self.do_plot()
        self.post_plot()

    def pre_plot(self):
        self.set_fontproperties()
        self.set_figsize()

    def do_plot(self):
        self.init_fig()
        self.make_plot()

    def post_plot(self):
        self.finish()

    def set_fontproperties(self):
        self.fp_axis_lab = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.axis_lab_font))
        self.fp_ticklabel = \
            mpl.font_manager.FontProperties(
                **copy.deepcopy(self.ticklabel_font))

    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        if hasattr(self, 'width') and hasattr(self, 'height'):
            self.figsize = (self.width, self.height)

    def init_fig(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """
        self.pdf = mpl.backends.backend_pdf.PdfPages(self.fname)
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)

    def make_plot(self):

        self.z = hc.linkage(self.data, method='average')

        self.ax = self.fig.add_subplot(1, 1, 1)

        self.dendro = \
            hc.dendrogram(self.z,
                          labels=self.data.columns,
                          color_threshold=0,
                          orientation='left',
                          ax=self.ax,
                          link_color_func=lambda x: self.color)

        _ = [
            tl.set_fontproperties(self.fp_ticklabel)
            for tl in self.ax.get_yticklabels()
        ]
        _ = [
            tl.set_fontproperties(self.fp_ticklabel)
            for tl in self.ax.get_xticklabels()
        ]

        self.ax.xaxis.grid(True, color='#FFFFFF', lw=1, ls='solid')
        self.ax.yaxis.grid(False)
        self.ax.set_axisbelow(True)
        self.ax.set_facecolor('#EAEAF2')
        list(map(lambda s: s.set_lw(0), self.ax.spines.values()))
        self.ax.tick_params(which='both', length=0)

    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """
        #self.gs.update(wspace=0.1, hspace=0.1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.fig.tight_layout()

        #self.fig.subplots_adjust(top = 0.94)
        self.cvs.draw()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.cvs.print_figure(self.pdf)

        self.pdf.close()
        self.fig.clf()
