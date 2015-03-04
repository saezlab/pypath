#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
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

import cairo
import math

class InterSet(object):
    
    def __init__(self, xsizes, intersects, outf = 'cairotest.pdf', 
            width = 1024, height = 1024, bgcol = 'embl_gray125', 
            cols = None, interscols = None, ysizes = None, 
            ycols = None, skip = 3.5, margin = 24, 
            mincircle = 5, cellpadding = 4):
        for key, val in locals().iteritems():
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
        self.clabels = [
            (0.5, math.sqrt(3) / -2.0), 
            (math.sqrt(3) / 2.0, -0.5),
            (math.sqrt(2) / 2.0, math.sqrt(2) / -2.0)
        ]
        self.palette = [self.colors[x] for x in \
            ['embl_green', 'embl_blue', 'embl_yellow']]
        self.fontpal = [self.colors[x] for x in \
            ['embl_gray875', 'white', 'embl_gray875']]
        self.font = 'HelveticaNeueLT Std Lt'
        self.bgcol = self.bgcol if type(self.bgcol) is tuple else self.colors[self.bgcol]
        ### set parameters for x:
        # list of column labels (first elements of tuples in x list)
        self.xlabs = [x[0] for x in self.xsizes]
        # colors of sets in column headers
        self.xcols = self.get_colors(self.xsizes)
        # set sizes:
        self.xsizes = self.get_sizes(self.xsizes)
        ### same for y:
        self.ylabs = [y[0] for y in self.ysizes] if self.ysizes is not None else self.xlabs
        self.ycols = self.get_colors(self.ysizes) \
            if self.ysizes is not None else self.xcols
        self.ysizes = self.get_sizes(self.ysizes) \
            if self.ysizes is not None else self.xsizes
        ### margin:
        # margin is either a single integer, or a tuple of 4 integers:
        self.margin = self.margin if type(self.margin) is tuple else (self.margin,) * 4
        ### table:
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
    
    def fit_text(self, txt, width, pt = 24, padding = 2):
        overf = 1
        while overf > 0:
            self.ctx.select_font_face(self.font, cairo.FONT_SLANT_NORMAL, \
                cairo.FONT_WEIGHT_NORMAL)
            self.ctx.set_font_size(pt)
            overf = self.ctx.text_extents(txt)[2] - width + padding
            pt *= 0.95
        return pt
    
    def get_colors(self, sizes, palette = None):
        palette = palette if palette is not None else self.palette
        return [[xx[1] if xx[1] is not None else palette[i] 
            for i, xx in enumerate(x[1:])] for x in sizes]
    
    def get_sizes(self, sizes):
        return [[xx[0] for i, xx in enumerate(x[1:])] for x in sizes]
    
    def scale_sizes(self, sizes):
        print sizes
        scaled = [math.sqrt(x / math.pi) for x in 
            self.scale([i for ii in sizes for i in ii], 
                self.minarea, 
                self.maxarea)]
        lens = [len(s) for s in sizes]
        return [scaled[sum(lens[:i]):sum(lens[:i]) + l] for i, l in enumerate(lens)]
    
    def rgb1(self, rgb256):
        return rgb256 if not any([i > 1 for i in rgb256]) \
            else tuple([x/float(255) for x in rgb256])
    
    def set_alpha(self, col, alpha):
        return tuple(list(col[0:3]) + [alpha])
    
    def draw_circles(self):
        self.clabpt = self.fit_text('00000', self.cellw / 7, padding = 0)
        y = self.margin[2] + self.cellh / float(2)
        for xi in range(2, len(self.xcoo)):
            x = self.margin[0] + sum(self.xcoo[:xi]) + self.cellw / float(2)
            for i, sz in enumerate(self.xcircsize[xi-2]):
                self.circle(x, y, sz, self.rgb1(self.xcols[xi-2][i]))
                self.label(str(self.xsizes[xi-2][i]), *tuple([a + b * sz 
                    for a, b in zip([x,y], list(self.clabels[i]))]), 
                    pt = self.clabpt, c = self.colors['embl_gray875'], 
                    center = True, vcenter = True, rot = -45)
        x = self.margin[0] + self.cellw / float(2)
        for yi in range(2, len(self.ycoo)):
            y = self.margin[2] + sum(self.ycoo[:yi]) + self.cellh / float(2)
            for i, sz in enumerate(self.ycircsize[yi-2]):
                self.circle(x, y, sz, self.rgb1(self.ycols[yi-2][i]))
                self.label(str(self.ysizes[yi-2][i]), *tuple([a + b * sz 
                    for a, b in zip([x,y], list(self.clabels[i]))]), 
                    pt = self.clabpt, c = self.colors['embl_gray875'], 
                    center = True, vcenter = True, rot = -45)
        # intersections
        for yi in range(2, len(self.ycoo)):
            y = self.margin[2] + sum(self.ycoo[:yi]) + self.cellh / float(2)
            for xi in range(2, len(self.xcoo)):
                x = self.margin[0] + sum(self.xcoo[:xi]) + self.cellw / float(2)
                for i, sz in enumerate(self.intersects[(self.xlabs[xi-2], 
                    self.ylabs[yi-2])]['ssize']):
                    self.circle(x, y, sz, 
                        self.rgb1(self.intersects[(self.xlabs[xi-2], 
                        self.ylabs[yi-2])]['color'][i]))
                    self.label(str(self.intersects[(self.xlabs[xi-2], 
                        self.ylabs[yi-2])]['size'][i]), *tuple([a + b * sz 
                        for a, b in zip([x,y], list(self.clabels[i]))]), 
                        pt = self.clabpt, c = self.colors['embl_gray875'], 
                        center = True, vcenter = True, rot = -45)
                    '''
                    self.label(str(self.intersects[(self.xlabs[xi-2], 
                        self.ylabs[yi-2])]['size']), x, y, 
                        self.intersects[(self.xlabs[xi-2], 
                            self.ylabs[yi-2])]['ssize'][i] * 0.5 + 3, self.colors['white'])
                    '''
    
    def circle(self, x, y, r, c):
        c = self.set_alpha(c, 0.5)
        self.ctx.set_source_rgba(*c)
        print 'Drawing circle at (%f, %f) of radius %f, color %s' % (x, y, r, str(c))
        self.ctx.arc(x, y, r, 0, 2 * math.pi)
        self.ctx.fill()
    
    def label(self, txt, x, y, pt, c, center = True, rot = 0.0, vcenter = False):
        c = self.rgb1(c)
        self.ctx.select_font_face(self.font, cairo.FONT_SLANT_NORMAL, \
            cairo.FONT_WEIGHT_NORMAL)
        self.ctx.set_font_size(pt)
        print 'Font color set to', c
        print 'Drawing label at %s'%str((x, y))
        self.ctx.save()
        self.ctx.translate(x, y)
        self.ctx.rotate(math.radians(rot))
        if center:
            ext = self.ctx.text_extents(txt)
            self.ctx.translate(- ext[2] / 2.0, 0.0)
        if vcenter:
            ext = self.ctx.text_extents(txt)
            self.ctx.translate(0.0, - ext[3] / 2.0)
        self.ctx.move_to(0, 0)
        self.ctx.set_source_rgba(*c)
        self.ctx.show_text(txt)
        self.ctx.restore()
    
    def colnames(self):
        # font size for labels:
        self.xlabpt = self.max_text(self.xlabs, self.cellw)
        y = self.margin[2] + self.ycoo[0] + self.ycoo[1]/2.0
        for xi in range(2, len(self.xcoo)):
            x = self.margin[0] + sum(self.xcoo[:xi]) + self.cellw / 2.0
            lab = self.xlabs[xi - 2]
            self.label(lab, x, y, self.xlabpt, self.colors['embl_gray875'])
    
    def rownames(self):
        # font size for labels:
        self.xlabpt = self.max_text(self.xlabs, self.cellw)
        x = self.margin[0] + self.xcoo[0] + self.xcoo[1]/2.0
        for yi in range(2, len(self.ycoo)):
            y = self.margin[2] + sum(self.ycoo[:yi]) + self.cellh / 2.0
            lab = self.ylabs[yi - 2]
            self.label(lab, x, y, self.xlabpt, self.colors['embl_gray875'], rot = -90)
    
    def draw_table(self):
        bg = self.rgb1(self.bgcol)
        self.ctx.set_source_rgba(*bg)
        for xi in range(0, len(self.xcoo)):
            for yi in range(0, len(self.ycoo)):
                ulx = self.margin[0] + sum(self.xcoo[:xi])
                uly = self.margin[2] + sum(self.ycoo[:yi])
                print 'Drawing rectangle at (%f, %f), size (%f, %f)' % \
                    (ulx + self.skip, uly + self.skip, self.xcoo[xi] - self.skip, 
                    self.ycoo[yi] - self.skip)
                self.ctx.rectangle(ulx + self.skip, uly + self.skip, 
                    self.xcoo[xi] - self.skip, self.ycoo[yi] - self.skip)
                self.ctx.stroke_preserve()
                self.ctx.fill()
    
    def cells(self, props, mi, ma):
        rng = ma - mi
        uni = rng / float(sum(props))
        return [x * uni for x in props]

    def scale(self, lst, lower, upper):
        print lst
        return [((x - min(set(lst) - set([0]))) / \
            float(max(lst) - min(set(lst) - set([0]))) * (upper - lower) + lower) \
            if x != 0 else 0 for x in lst]
    