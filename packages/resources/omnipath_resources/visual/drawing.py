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

import os
import sys

try:
    import cairo
except:
    sys.stdout.write('Module `cairo` is not available.'
                     '\nSome plotting functionalities won\'t be accessible.\n')

import math
import time

from pypath.share.common import *
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.share.session as session_mod

__all__ = ['Plot', 'InterSet']

try:
    import igraph

except ModuleNotFoundError:
    sys.stdout.write('Module `igraph` is not available.'
                     '\nSome plotting functionalities won\'t be accessible.\n')


class Plot(session_mod.Logger):

    def __init__(self,
                 graph=None,
                 filename=None,
                 graphix_dir="pdf",
                 graphix_format="pdf",
                 name=None,
                 title_text=None,
                 title_font_family=None,
                 title_font_size=None,
                 title_color='#646567',
                 size=None,
                 layout="fruchterman_reingold",
                 layout_param=None,
                 vertex_label=None,
                 vertex_size=None,
                 vertex_label_size='degree_label_size',
                 edge_width=None,
                 vertex_color='#6EA945',
                 vertex_label_color='#007B7F',
                 vertex_alpha='AA',
                 vertex_frame_color='#FFFFFF00',
                 vertex_frame_width=0,
                 edge_label=None,
                 edge_label_size=None,
                 edge_label_color='#007B7F',
                 edge_curved=None,
                 edge_color='#818284',
                 edge_alpha='AA',
                 autocurve=None,
                 vertex_label_font="sans-serif",
                 edge_label_font="sans-serif",
                 edge_arrow_size=1.0,
                 edge_arrow_width=1.0,
                 palettes={},
                 bbox=None,
                 margin=10,
                 small=None,
                 dimensions=(1280, 1280),
                 grouping=None,
                 **kwargs):
        # setting parameters:
        for key, val in iteritems(locals()):
            setattr(self, key, val)
        self.default_alpha = {
            'vertex_color': 'AA',
            'edge_color': 'AA',
            'vertex_label_color': 'FF',
            'vertex_frame_color': '00',
            'edge_label_color': 'FF'
        }
        self.default_vertex_label_size = 6.0
        self.plots = []
        self.session = common.random_string()
        self.name = self.name if self.name is not None else self.session
        self.label_sizes = {
            'small': (15.0, 13.7),
            'medium': (13.0, 10.0),
            'large': (9.0, 6.0)
        }
        self.palettes = {
            'vertex': ['#6EA945', '#007B7F', '#FCCC06', '#DA0025', '#000000'],
            'edge': ['#007B7F', '#6EA945', '#DA0025'],
            'vertex_label': ['#454447'],
            'edge_label': ['#454447']
        }
        self.small_param = {
            'vertex_size': 21,
            'edge_width': 0.051,
            'autocurve': True,
            'vertex_label_dist': 1.5
        }
        self.medium_param = {
            'vertex_size': 7,
            'edge_width': 0.051,
            'autocurve': True,
            'vertex_label_dist': 1.33,
            'edge_label_size': 1.0
        }
        self.large_param = {
            'vertex_size': 2,
            'edge_width': 0.051,
            'vertex_label_dist': 1.0,
            'edge_label_size': 1.0
        }
        self.layout_defaults = {
            'fruchterman_reingold': {
                'repulserad': self.graph.vcount()**2.8,
                'maxiter': 1000,
                'area': self.graph.vcount()**2.3
            }
        }
        self.update_page()
        self.update_graph(graph)

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def _has_graph(self):
        return type(self.graph) is igraph.Graph

    def update_page(self, size=None, margin=None):
        self.dimensions = self.dimensions if size is None else size
        self.margin = self.margin if margin is None else margin
        if type(self.margin) is int:
            self.margin = tuple([self.margin] * 4)
        self.bbox = igraph.drawing.utils.BoundingBox(
            self.margin[0], self.margin[1],
            self.dimensions[0] - self.margin[2],
            self.dimensions[1] - self.margin[3])

    def title(self, title, family=None, size=None, color=None):
        self.title_text = title
        self.title_font_family = family if family is not None else self.vertex_label_font
        self.title_font_size = size if size is not None else \
            math.ceil(self.bbox.width / 50.0)
        self.title_color = color if color is not None else self.title_color

    def update_graph(self, graph, title=None):
        self.graph = graph
        if title is not None:
            self.title_text = title
        if self._has_graph():
            if type(self.vertex_label) is str and \
                    self.vertex_label in self.graph.vs.attributes():
                self.vertex_label = self.graph.vs[self.vertex_label]
            self.graph.vertex_label_font = self.vertex_label_font
            self.graph.edge_label_font = self.edge_label_font
            self.size = 'small' if self.graph.vcount() <= 100 \
                else 'medium' if self.graph.vcount() <= 500 \
                else 'large'
            getattr(self, '%s_defaults' % self.size)()
            self.update_layout()
            self.colorize('vertex', self.vertex_color, self.vertex_alpha)
            self.colorize('edge', self.edge_color, self.edge_alpha)
            self.colorize('vertex', self.vertex_label_color, attr='label')
            self.colorize('edge', self.edge_label_color, attr='label')
            self.autocurve = True if self.small else False
        if self.title_text is not None:
            self.title(self.title_text, self.title_font_family,
                       self.title_font_size, self.title_color)

    def update_layout(self, layout=None, param=None):
        g = self.graph
        if self.layout_update_needed(layout, param):
            start_time = time.time()
            self.layout_param = self.layout_param if self.layout_param is not None \
                else self.layout_defaults[self.layout] if self.layout in self.layout_defaults \
                else {}
            msg = """Calculating %s layout... (numof nodes/edges: %u/%u)""" % \
                (self.layout, g.vcount(), g.ecount())
            sys.stdout.write('\t::%s' % msg)
            sys.stdout.flush()
            self._log(msg)
            if self.grouping is not None:
                if self.layout in [
                        "intergroup", "modular_fr", "modular_circle"
                ]:
                    f = getattr(gr_plot, "layout_%s" % self.layout)
                    f(g, self.grouping, self.layout_param)
                else:
                    if self.layout not in [
                            "fruchterman_reingold", "fr", "circle"
                    ]:
                        self.layout = "fr"
                    g['layout_data'] = layout_intergroup(g, self.grouping,
                                                         **self.layout_param)
                    g['layout_type'] = "layout_intergroup"
                if self.vertex_color == "groups":
                    self.colorize(what='vertex', coldef=self.grouping)
            self.layout = self.layout if layout is None else layout
            self.layout_param = self.layout_param if param is None else param
            self.layout_data = g.layout(self.layout, **self.layout_param)
            self.layout_data = igraph.Layout(self.layout_data)
            time_elapsed = time.time() - start_time
            m, s = divmod(time_elapsed, 60)
            time_elapsed_str = "%02d min %02d sec" % (m, s)
            sys.stdout.write(' Done in %s. \n' % time_elapsed_str)
            sys.stdout.flush()

    def update_filename(self):
        if self.name is None:
            self.name = 'plot'
        seq = str(len(self.plots) + 1)
        self.nextfile = self.filename if self.filename is not None \
            else 'network-%s-%u.%s' % \
            (self.name, len(self.plots) + 1, self.graphix_format)
        if os.path.sep not in self.nextfile:
            self.nextfile = os.path.join(self.graphix_dir, self.nextfile)

    def layout_update_needed(self, layout, param):
        return not hasattr(self, "layout") or \
            (layout is not None and self.layout != layout) or \
            not hasattr(self, "layout_data") or \
            len(self.layout_data) != self.graph.vcount() or \
            not hasattr(self, "layout_param") or \
            (param is not None and self.layout_param != param)

    def colorize(self,
                 what='vertex',
                 coldef=None,
                 alpha=None,
                 attr=None,
                 palette=None):
        attr = '%scolor' % ('' if attr is None else '%s_' % attr)
        seq = self.graph.vs if what == 'vertex' else self.graph.es
        coldef = coldef if coldef is not None else self.palette[0]
        alpha = alpha if alpha is not None else self.default_alpha['%s_%s' % (
            what, attr)]
        pal = self.palettes[what] if palette is None else palette
        if type(coldef) in _const.CHAR_TYPES and coldef in seq.attributes():
            lev = list(set(seq[coldef]))
            seq[attr] = [pal[lev.index(i[coldef])] for i in seq]
        elif type(coldef) in _const.CHAR_TYPES and len(coldef) <= 9:
            seq[attr] = [coldef for _ in seq]
        elif type(coldef) is list and len(coldef) == len(seq):
            seq[attr] = coldef
        elif type(coldef) is dict and '__attr__' in coldef:
            seq[attr] = [coldef[i[coldef['__attr__']]] for i in seq]
        elif hasattr(coldef, '__call__'):
            seq[attr] = [coldef(i) for i in seq]
        if min([len(i[attr]) for i in seq]) == 7:
            self.set_alpha(seq, alpha, attr)

    def set_param(self, param, value):
        if not hasattr(self, param) or getattr(self, param) is None:
            setattr(self, param, value)

    def small_defaults(self):
        self.set_defaults('small_param')

    def medium_defaults(self):
        self.set_defaults('medium_param')

    def large_defaults(self):
        self.set_defaults('large_param')

    def set_defaults(self, preset):
        if hasattr(self, preset):
            for k, v in iteritems(getattr(self, preset)):
                self.set_param(k, v)

    def set_alpha(self, seq, alpha, attr):
        seq[attr] = ['%s%s' % (c[0:7], alpha) for c in seq[attr]]

    def hex2rgb(self, rgbhex):
        rgbhex = rgbhex.lstrip('#')
        lv = len(rgbhex)
        return tuple(int(rgbhex[i:i + 2], 16) for i in range(0, lv, 2))

    def rgb1(self, rgb256):
        return rgb256 if not any([i > 1 for i in rgb256]) \
            else tuple([x / float(255) for x in rgb256])

    def make_title(self):
        ctx = cairo.Context(self.plots[-1].surface)
        ctx.set_font_size(self.title_font_size)
        ctx.select_font_face(self.title_font_family, cairo.FONT_SLANT_NORMAL,
                             cairo.FONT_WEIGHT_NORMAL)
        ctx.set_source_rgba(*self.rgb1(self.hex2rgb(self.title_color)))
        title_drawer = igraph.drawing.text.TextDrawer(
            ctx, self.title_text, halign=igraph.drawing.text.TextDrawer.CENTER)
        title_drawer.draw_at(0, 40, width=self.bbox.width)

    def draw(self, return_data=False, **kwargs):
        if not self._has_graph():
            return None
        self.update_filename()
        g = self.graph
        msg = """Plotting %s to file %s...""" % \
            (self.graphix_format, self.nextfile)
        sys.stdout.write('\t::%s' % msg)
        sys.stdout.flush()
        self._log(msg)
        if self.graphix_format == "pdf":
            sf = cairo.PDFSurface(self.nextfile, self.dimensions[0],
                                  self.dimensions[1])
        else:
            # currently doing only pdf
            sf = self.nextfile
        if self.vertex_label_size == "degree_label_size":
            # TODO
            dgr = g.vs.degree()
            maxDgr = float(max(dgr))
            g.vs["label_size"] = [None]
            g.vs["label_size"] = [
                math.log(float(v.degree()) / maxDgr + 1.0) *
                self.label_sizes[self.size][0] + self.label_sizes[self.size][1]
                for v in g.vs
            ]
            self.vertex_label_size = g.vs["label_size"]
        elif type(self.vertex_label_size) is not int:
            self.vertex_label_size = self.default_vertex_label_size
        if type(self.edge_curved) is float:
            self.kwargs['edge_curved'] = self.edge_curved
        else:
            self.kwargs['autocurve'] = self.autocurve

        if igraph.__version__ < '0.9.4':
            import pypath.visual.igraph_drawing as ig_drawing
            self.kwargs['drawer_factory'] = (
                ig_drawing.DefaultGraphDrawerFFsupport
            )

        self.plots.append(
            igraph.plot(
                g,
                layout=self.layout_data,
                target=sf,
                bbox=self.bbox,
                vertex_size=self.vertex_size,
                vertex_frame_width=self.vertex_frame_width,
                vertex_label=self.vertex_label,
                vertex_label_size=self.vertex_label_size,
                edge_label=self.edge_label,
                edge_width=self.edge_width,
                edge_arrow_size=self.edge_arrow_size,
                edge_arrow_width=self.edge_arrow_width,
                vertex_label_dist=self.vertex_label_dist,
                **self.kwargs))
        self.plots[-1].redraw()
        if self.title_text is not None:
            self.make_title()
        self.plots[-1].save()
        sys.stdout.write('Ready.\n')
        sys.stdout.flush()
        self._log("Plot saved to %s" % self.nextfile)
        if return_data:
            return (self.plots[-1], g, self.layout_data, sf, self.bbox,
                    self.vertex_size,
                    self.vertex_frame_width, self.vertex_label,
                    self.vertex_label_size, self.edge_width, self.edge_curved,
                    self.edge_arrow_size, self.edge_arrow_width, self.kwargs)


#TODO this class may be dropped?
class InterSet(object):
    def __init__(self,
                 xsizes,
                 intersects,
                 outf='cairotest.pdf',
                 width=1024,
                 height=1024,
                 bgcol='embl_gray125',
                 cols=None,
                 interscols=None,
                 ysizes=None,
                 ycols=None,
                 skip=3.5,
                 margin=24,
                 mincircle=5,
                 cellpadding=4):
        for key, val in iteritems(locals()):
            setattr(self, key, val)
        self.colors = {
            'embl_green': (115, 179, 96, 255),
            'embl_blue': (0, 102, 102, 255),
            'embl_yellow': (250, 183, 0, 255),
            'embl_red': (227, 62, 62, 255),
            'embl_black': (0, 0, 0, 255),
            'embl_gray875': (32, 32, 32, 255),
            'embl_gray75': (64, 64, 64, 255),
            'embl_gray625': (96, 96, 96, 255),
            'embl_gray50': (128, 128, 128, 255),
            'embl_gray25': (192, 192, 192, 255),
            'embl_gray125': (224, 224, 224, 255),
            'white': (255, 255, 255, 255)
        }
        # positions of circle labels:
        self.clabels = [(0.5, math.sqrt(3) / -2.0), (math.sqrt(3) / 2.0, -0.5),
                        (math.sqrt(2) / 2.0, math.sqrt(2) / -2.0)]
        self.palette = [
            self.colors[x]
            for x in ['embl_green', 'embl_blue', 'embl_yellow']
        ]
        self.fontpal = [
            self.colors[x] for x in ['embl_gray875', 'white', 'embl_gray875']
        ]
        self.font = 'HelveticaNeueLT Std Lt'
        self.bgcol = self.bgcol if type(self.bgcol) is tuple else self.colors[
            self.bgcol]
        # set parameters for x:
        # list of column labels (first elements of tuples in x list)
        self.xlabs = [x[0] for x in self.xsizes]
        # colors of sets in column headers
        self.xcols = self.get_colors(self.xsizes)
        # set sizes:
        self.xsizes = self.get_sizes(self.xsizes)
        # same for y:
        self.ylabs = [y[0] for y in self.ysizes
                      ] if self.ysizes is not None else self.xlabs
        self.ycols = self.get_colors(self.ysizes) \
            if self.ysizes is not None else self.xcols
        self.ysizes = self.get_sizes(self.ysizes) \
            if self.ysizes is not None else self.xsizes
        # margin:
        # margin is either a single integer, or a tuple of 4 integers:
        self.margin = self.margin if type(self.margin) is tuple else (
            self.margin, ) * 4
        # table:
        # proportions of cell sizes:
        self.xcellsize = [3, 1] + [3] * len(self.xsizes)
        self.ycellsize = [3, 1] + [3] * len(self.ysizes)
        # sizes of table cells:
        self.xcoo = self.cells(self.xcellsize, self.margin[0],
                               self.width - self.margin[1])
        self.ycoo = self.cells(self.ycellsize, self.margin[2],
                               self.height - self.margin[3])
        # width and height of diagram cells:
        self.cellw = self.xcoo[0]
        self.cellh = self.ycoo[0]
        # largest circle fit in the cells:
        self.maxcircle = min(self.cellw, self.cellh) / 2 - 2 * self.cellpadding
        self.maxarea = pow(self.maxcircle, 2) * math.pi
        self.minarea = pow(self.mincircle, 2) * math.pi
        # scaling circle sizes between min and max circle size
        self.xcircsize = self.scale_sizes(self.xsizes)
        self.ycircsize = self.scale_sizes(self.ysizes)
        ssize = self.scale_sizes([x['size'] for x in self.intersects.values()])
        for i, k in enumerate(self.intersects.keys()):
            self.intersects[k]['ssize'] = ssize[i]
            if 'color' not in self.intersects[k]:
                self.intersects[k]['color'] = self.palette[0:len(ssize[i])]

    def draw(self):
        '''
        main function of this class
        '''
        self.srf = cairo.PDFSurface(self.outf, self.width, self.height)
        self.ctx = cairo.Context(self.srf)
        self.draw_table()
        self.colnames()
        self.rownames()
        self.draw_circles()
        self.srf.finish()
        self.srf.flush()

    def max_text(self, labels, width):
        pts = []
        for lab in labels:
            pts.append(self.fit_text(lab, width))
        return min(pts)

    def fit_text(self, txt, width, pt=24, padding=2):
        overf = 1
        while overf > 0:
            self.ctx.select_font_face(self.font, cairo.FONT_SLANT_NORMAL,
                                      cairo.FONT_WEIGHT_NORMAL)
            self.ctx.set_font_size(pt)
            overf = self.ctx.text_extents(txt)[2] - width + padding
            pt *= 0.95
        return pt

    def get_colors(self, sizes, palette=None):
        palette = palette if palette is not None else self.palette
        return [[
            xx[1] if xx[1] is not None else palette[i]
            for i, xx in enumerate(x[1:])
        ] for x in sizes]

    def get_sizes(self, sizes):
        return [[xx[0] for i, xx in enumerate(x[1:])] for x in sizes]

    def scale_sizes(self, sizes):
        scaled = [
            math.sqrt(x / math.pi)
            for x in self.scale([i for ii in sizes
                                 for i in ii], self.minarea, self.maxarea)
        ]
        lens = [len(s) for s in sizes]
        return [
            scaled[sum(lens[:i]):sum(lens[:i]) + l] for i, l in enumerate(lens)
        ]

    def rgb1(self, rgb256):
        return rgb256 if not any([i > 1 for i in rgb256]) \
            else tuple([x / float(255) for x in rgb256])

    def set_alpha(self, col, alpha):
        return tuple(list(col[0:3]) + [alpha])

    def draw_circles(self):
        self.clabpt = self.fit_text('00000', self.cellw / 7, padding=0)
        y = self.margin[2] + self.cellh / float(2)
        for xi in range(2, len(self.xcoo)):
            x = self.margin[0] + sum(self.xcoo[:xi]) + self.cellw / float(2)
            for i, sz in enumerate(self.xcircsize[xi - 2]):
                self.circle(x, y, sz, self.rgb1(self.xcols[xi - 2][i]))
                self.label(
                    str(self.xsizes[xi - 2][i]),
                    *tuple([
                        a + b * sz
                        for a, b in zip([x, y], list(self.clabels[i]))
                    ]),
                    pt=self.clabpt,
                    c=self.colors['embl_gray875'],
                    center=True,
                    vcenter=True,
                    rot=-45)
        x = self.margin[0] + self.cellw / float(2)
        for yi in range(2, len(self.ycoo)):
            y = self.margin[2] + sum(self.ycoo[:yi]) + self.cellh / float(2)
            for i, sz in enumerate(self.ycircsize[yi - 2]):
                self.circle(x, y, sz, self.rgb1(self.ycols[yi - 2][i]))
                self.label(
                    str(self.ysizes[yi - 2][i]),
                    *tuple([
                        a + b * sz
                        for a, b in zip([x, y], list(self.clabels[i]))
                    ]),
                    pt=self.clabpt,
                    c=self.colors['embl_gray875'],
                    center=True,
                    vcenter=True,
                    rot=-45)
        # intersections
        for yi in range(2, len(self.ycoo)):
            y = self.margin[2] + sum(self.ycoo[:yi]) + self.cellh / float(2)
            for xi in range(2, len(self.xcoo)):
                x = self.margin[0] + \
                    sum(self.xcoo[:xi]) + self.cellw / float(2)
                for i, sz in enumerate(self.intersects[(self.xlabs[
                        xi - 2], self.ylabs[yi - 2])]['ssize']):
                    self.circle(x, y, sz,
                                self.rgb1(self.intersects[(self.xlabs[
                                    xi - 2], self.ylabs[yi - 2])]['color'][i]))
                    self.label(
                        str(self.intersects[(self.xlabs[xi - 2], self.ylabs[
                            yi - 2])]['size'][i]),
                        *tuple([
                            a + b * sz
                            for a, b in zip([x, y], list(self.clabels[i]))
                        ]),
                        pt=self.clabpt,
                        c=self.colors['embl_gray875'],
                        center=True,
                        vcenter=True,
                        rot=-45)
                    '''
                    self.label(str(self.intersects[(self.xlabs[xi-2],
                        self.ylabs[yi-2])]['size']), x, y,
                        self.intersects[(self.xlabs[xi-2],
                            self.ylabs[yi-2])]['ssize'][i] * 0.5 + 3, self.colors['white'])
                    '''

    def circle(self, x, y, r, c):
        c = self.set_alpha(c, 0.5)
        self.ctx.set_source_rgba(*c)
        self.ctx.arc(x, y, r, 0, 2 * math.pi)
        self.ctx.fill()

    def label(self, txt, x, y, pt, c, center=True, rot=0.0, vcenter=False):
        c = self.rgb1(c)
        self.ctx.select_font_face(self.font, cairo.FONT_SLANT_NORMAL,
                                  cairo.FONT_WEIGHT_NORMAL)
        self.ctx.set_font_size(pt)
        self.ctx.save()
        self.ctx.translate(x, y)
        self.ctx.rotate(math.radians(rot))
        if center:
            ext = self.ctx.text_extents(txt)
            self.ctx.translate(-ext[2] / 2.0, 0.0)
        if vcenter:
            ext = self.ctx.text_extents(txt)
            self.ctx.translate(0.0, -ext[3] / 2.0)
        self.ctx.move_to(0, 0)
        self.ctx.set_source_rgba(*c)
        self.ctx.show_text(txt)
        self.ctx.restore()

    def colnames(self):
        # font size for labels:
        self.xlabpt = self.max_text(self.xlabs, self.cellw)
        y = self.margin[2] + self.ycoo[0] + self.ycoo[1] / 2.0
        for xi in range(2, len(self.xcoo)):
            x = self.margin[0] + sum(self.xcoo[:xi]) + self.cellw / 2.0
            lab = self.xlabs[xi - 2]
            self.label(lab, x, y, self.xlabpt, self.colors['embl_gray875'])

    def rownames(self):
        # font size for labels:
        self.xlabpt = self.max_text(self.xlabs, self.cellw)
        x = self.margin[0] + self.xcoo[0] + self.xcoo[1] / 2.0
        for yi in range(2, len(self.ycoo)):
            y = self.margin[2] + sum(self.ycoo[:yi]) + self.cellh / 2.0
            lab = self.ylabs[yi - 2]
            self.label(
                lab, x, y, self.xlabpt, self.colors['embl_gray875'], rot=-90)

    def draw_table(self):
        bg = self.rgb1(self.bgcol)
        self.ctx.set_source_rgba(*bg)
        for xi in range(0, len(self.xcoo)):
            for yi in range(0, len(self.ycoo)):
                ulx = self.margin[0] + sum(self.xcoo[:xi])
                uly = self.margin[2] + sum(self.ycoo[:yi])
                # print 'Drawing rectangle at (%f, %f), size (%f, %f)' % \
                #    (ulx + self.skip, uly + self.skip, self.xcoo[xi] - self.skip,
                #    self.ycoo[yi] - self.skip)
                self.ctx.rectangle(ulx + self.skip, uly + self.skip,
                                   self.xcoo[xi] - self.skip,
                                   self.ycoo[yi] - self.skip)
                self.ctx.stroke_preserve()
                self.ctx.fill()

    def cells(self, props, mi, ma):
        rng = ma - mi
        uni = rng / float(sum(props))
        return [x * uni for x in props]

    def scale(self, lst, lower, upper):
        return [((x - min(set(lst) - set([0]))
                  ) / float(max(lst) - min(set(lst) - set([0]))) *
                 (upper - lower) + lower) if x != 0 else 0 for x in lst]
