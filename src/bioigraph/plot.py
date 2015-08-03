#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import re
import sys
import os
import itertools

import numpy as np
from numpy.random import randn
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as hc
import hcluster as hc2
import matplotlib.gridspec as gridspec
from matplotlib import ticker

from progress import Progress
from common import uniqList, ROOT

# color converting functions
def rgb2hex(rgb):
    return '#%02x%02x%02x' % rgb

def hex2rgb(rgbhex):
    rgbhex = rgbhex.lstrip('#')
    lv = len(rgbhex)
    return tuple(int(rgbhex[i:i + 2], 16) for i in range(0, lv, 2))

def rgb1(rgb256):
    return rgb256 if not any([i > 1 for i in rgb256]) \
        else tuple([x / float(255) for x in rgb256])

def rgb256(rgb1):
    return rgb1 if any([i > 1.0 for i in rgb1]) \
        else tuple([x * 255.0 for x in rgb1])

class Plot(object):
    
    def __init__(self, fname = None, font_family = 'Helvetica Neue LT Std', 
        palette = None, context = 'poster', rc = None):
        for k, v in locals().iteritems():
            if not hasattr(self, k) or getattr(self, k) is not None:
                setattr(self, k, v)
        self.palette = palette or self.embl_palette()
    
    def embl_palette(self, inFile = 'embl_colors'):
        cols = []
        inFile = os.path.join(ROOT, 'data', inFile)
        with open(inFile, 'r') as f:
            series = []
            for i, l in enumerate(f):
                l = [x.strip() for x in l.split(',')]
                series.append(rgb2hex(tuple([256 * float(x) for x in l[0:3]])))
                if len(series) == 7:
                    cols.append(series)
                    series = []
        return cols
    
    def finish(self):
        '''
        Saves and closes a figure.
        '''
        self.fig.tight_layout()
        self.fig.savefig(self.fname)
        plt.close(self.fig)


class Barplot(Plot):
    
    def __init__(self, x, y, data = None, fname = None, font_family = 'Helvetica Neue LT Std', 
        xlab = '', ylab = '', lab_angle = 90, lab_size = 9, color = '#007b7f', 
        order = False, desc = True, legend = None, fin = True, 
        y_break = None, rc = None, palette = None, context = 'poster', **kwargs):
        '''
        y_break : tuple
        If not None, the y-axis will have a break. 2 floats in the tuple, < 1.0, 
        mean the lower and upper proportion of the plot shown. The part between
        them will be hidden. E.g. y_break = (0.3, 0.1) shows the lower 30% and 
        upper 10%, but 60% in the middle will be cut out.
        '''
        for k, v in locals().iteritems():
            setattr(self, k, v)
        self.sns = sns
        self.rc = self.rc or {'lines.linewidth': 1.0, 'patch.linewidth': 0.0, 
            'grid.linewidth': 1.0}
        super(Barplot, self).__init__(fname = fname, font_family = font_family, 
            palette = palette, context = context, rc = self.rc)
        self.color = self.color or self.palette[0][0]
        if type(self.color) is list:
            self.palette = sns.color_palette(self.color)
        elif self.color is not None:
            self.palette = sns.color_palette([self.color] * len(self.x))
        self.color = None
        self.plot(**kwargs)
    
    def plot(self, x = None, y = None, **kwargs):
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
        self.ax = sns.barplot(self.x, y = self.y, data = None, color = self.color, 
            order = self.ordr, ax = self.ax, palette = self.palette, **kwargs)
        if self.y_break:
            self._break_y_axis(**kwargs)
        self.labels()
        if self.fin:
            self.finish()
    
    def sort(self):
        colcyc = itertools.cycle(list(self.palette))
        palcyc = [colcyc.next() for _ in xrange(len(self.x))]
        if self.order == 'x':
            self.ordr = np.array([self.x[i] for i in self.x.argsort()])
            self.palette = sns.color_palette([palcyc[i] for i in self.x.argsort()])
        elif self.order == 'y':
            self.ordr = np.array([self.x[i] for i in self.y.argsort()])
            #print 'palette 1 %s' % str(self.palette)
            #print 'palcyc %s' % str(palcyc)
            self.palette = sns.color_palette([palcyc[i] for i in self.y.argsort()])
            #print 'argsort %s' % str(list(self.y.argsort()))
            #print 'palette 2 %s' % str(self.palette)
        else:
            self.ordr = self.x
        if self.desc:
            self.ordr = self.ordr[::-1]
            self.palette = sns.color_palette(list(self.palette)[::-1])
    
    def _break_y_gs(self):
        self.gs = gridspec.GridSpec(2, 1, 
            height_ratios = [self.y_break[1] / sum(self.y_break),
                self.y_break[0] / sum(self.y_break)])
        self.fig = plt.figure()
        self.ax2 = self.fig.add_subplot(self.gs[0])
        self.ax = self.fig.add_subplot(self.gs[1])
    
    def _break_y_axis(self, **kwargs):
        self.ax2 = self.sns.barplot(self.x, y = self.y, data = None, 
            color = self.color, order = self.ordr, ax = self.ax2, 
            palette = self.palette, **kwargs)
        self.ax2.yaxis.set_major_locator(
            ticker.MaxNLocator(nbins = int(9/sum(self.y_break) + 1), steps = [1, 2, 5, 10]))
        self._originalYticks = self.ax2.get_yticks()
        ymin, ymax = self.ax.get_ylim()
        ymax = min(ytick for ytick in self.ax.get_yticks() if ytick > max(self.y))
        self.ax.set_ylim((ymin, ymax * self.y_break[0]))
        self.ax2.set_ylim((ymax - ymax * self.y_break[1], ymax))
        self.lower_y_min, self.lower_y_max = self.ax.get_ylim()
        self.upper_y_min, self.upper_y_max = self.ax2.get_ylim()
        plt.subplots_adjust(hspace = 0.08)
        self.ax2.spines['bottom'].set_visible(False)
        plt.setp(self.ax2.xaxis.get_majorticklabels(), visible = False)
        self.ax.spines['top'].set_visible(False)
        self.ax.set_yticks([yt for yt in self._originalYticks \
            if yt >= self.lower_y_min and yt <= self.lower_y_max])
        self.ax2.set_yticks([yt for yt in self._originalYticks \
            if yt >= self.upper_y_min and yt <= self.upper_y_max])
    
    def seaborn_style(self, context = None, rc = None):
        self.sns.set(font = self.font_family)
        self.sns.set_context(context or self.context, rc = rc or self.rc)
    
    def labels(self):
        for tick in self.ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(self.lab_size)
        self.ax.set_ylabel(self.ylab)
        self.ax.set_xlabel(self.xlab)
        plt.setp(self.ax.xaxis.get_majorticklabels(), rotation = self.lab_angle)
        if type(self.legend) is dict:
            legend_patches = [mpatches.Patch(color = col, label = lab) \
                for lab, col in self.legend.iteritems()]
            self.ax.legend(handles = legend_patches)

def boxplot(data, labels, xlab, ylab, fname, fontfamily = 'Helvetica Neue LT Std',
    textcol = 'black', violin = False):
    fig, ax = plt.subplots()
    sns.set(font = fontfamily)
    if violin:
        ax = sns.violinplot(data, names = labels, 
            color = embl_colors, linewidth = 0.1, saturation = 0.66)
    else:
        ax = sns.boxplot(data, names = labels, 
            color = embl_colors, linewidth = 0.1, saturation = 0.66)
    ax.set_xlabel(xlab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    ax.set_ylabel(ylab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
        tick.label.set_color(textcol)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(11)
        tick.label.set_color(textcol)
    fig.savefig(fname)

def stacked_barplot(x, y, data, fname, names, font_family = 'Helvetica Neue LT Std', 
    xlab = '', ylab = '', lab_angle = 90, lab_size = 9, legend = True, 
    colors = ['#007b7f', '#6ea945', '#fccc06', '#818284', '#da0025'], 
    order = False, desc = True):
    if type(x) is list or type(x) is tuple:
        x = np.array(x)
    for i, yi in enumerate(y):
        if type(yi) is list or type(yi) is tuple:
            y[i] = np.array(yi)
    total = np.array([sum([y[j][i] for j in xrange(len(y))]) for i in xrange(len(y[0]))])
    if order == 'x':
        ordr = np.array([x[i] for i in x.argsort()])
    elif order == 'y':
        ordr = np.array([x[i] for i in total.argsort()])
    elif type(order) is int:
        ordr = np.array([x[i] for i in y[order].argsort()])
    else:
        ordr = x
    if desc:
        ordr = ordr[::-1]
    fig, ax = plt.subplots()
    sns.set(font = font_family)
    sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
    for j in xrange(len(y), 0, -1):
        this_level = np.array([sum([y[jj][i] for jj in xrange(j)]) \
            for i in xrange(len(y[0]))])
        ax = sns.barplot(x, y = this_level, data = data, 
            color = colors[j-1], order = ordr)
    #ax = sns.barplot(x, y = y, data = data, color = color, x_order = ordr)
    #plt.bar(range(len(ordr)), [y[i] for i in ordr], align = 'center')
    #plt.xticks(list(ordr), [x[i] for i in ordr])
    sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(lab_size)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = lab_angle)
    self.finish(fig, fname)