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

"""
Drawing routines to draw graphs.

This module contains routines to draw graphs on:

  - Cairo surfaces (L{DefaultGraphDrawer})
  - UbiGraph displays (L{UbiGraphDrawer}, see U{http://ubietylab.net/ubigraph})

It also contains routines to send an igraph graph directly to
(U{Cytoscape<http://www.cytoscape.org>}) using the
(U{CytoscapeRPC plugin<http://gforge.nbic.nl/projects/cytoscaperpc/>}), see
L{CytoscapeGraphDrawer}. L{CytoscapeGraphDrawer} can also fetch the current
network from Cytoscape and convert it to igraph format.
"""

import sys
from collections import defaultdict
try:
    from itertools import izip
except ImportError:
    izip = zip

from math import atan2, cos, pi, sin, tan
from warnings import warn

try:
    from igraph._igraph import convex_hull, VertexSeq
    try:
        from igraph.compat import property
    except ModuleNotFoundError:
        pass
    from igraph.configuration import Configuration
    from igraph.drawing.baseclasses import AbstractDrawer, AbstractCairoDrawer, \
        AbstractXMLRPCDrawer
    from igraph.drawing.colors import color_to_html_format, color_name_to_rgb
    from pypath.visual.igraph_drawing.edge import ArrowEdgeDrawer
    from igraph.drawing.text import TextAlignment, TextDrawer
    from igraph.drawing.metamagic import AttributeCollectorBase
    from igraph.drawing.shapes import PolygonDrawer
    from igraph.drawing.utils import Point
    from pypath.visual.igraph_drawing.vertex import DefaultVertexDrawer
    from igraph.layout import Layout
    from igraph.drawing.graph import AbstractCairoGraphDrawer
    try:
        import cairo
    except ModuleNotFoundError:
        # No cairo support is installed. Create a fake module
        # pylint: disable-msg=C0103
        from igraph.drawing.utils import FakeModule
        cairo = FakeModule("igraph module could not be imported")

except ModuleNotFoundError:
    sys.stdout.write('Module `igraph` is not available.'
                     '\nSome plotting functionalities won\'t be accessible.\n')
    class AbstractCairoGraphDrawer:
        pass
    DefaultVertexDrawer = ArrowEdgeDrawer = TextDrawer = AbstractCairoGraphDrawer

__all__ = ["DefaultGraphDrawerFFsupport"]
__license__ = "GPL"



#####################################################################


class DefaultGraphDrawerFFsupport(AbstractCairoGraphDrawer):
    """Class implementing the default visualisation of a graph.

    The default visualisation of a graph draws the nodes on a 2D plane
    according to a given L{Layout}, then draws a straight or curved
    edge between nodes connected by edges. This is the visualisation
    used when one invokes the L{plot()} function on a L{Graph} object.

    See L{Graph.__plot__()} for the keyword arguments understood by
    this drawer."""

    def __init__(self,
                 context,
                 bbox,
                 vertex_drawer_factory=DefaultVertexDrawer,
                 edge_drawer_factory=ArrowEdgeDrawer,
                 label_drawer_factory=TextDrawer):
        """Constructs the graph drawer and associates it to the given
        Cairo context and the given L{BoundingBox}.

        @param context: the context on which we will draw
        @param bbox:    the bounding box within which we will draw.
                        Can be anything accepted by the constructor
                        of L{BoundingBox} (i.e., a 2-tuple, a 4-tuple
                        or a L{BoundingBox} object).
        @param vertex_drawer_factory: a factory method that returns an
                        L{AbstractCairoVertexDrawer} instance bound to a
                        given Cairo context. The factory method must take
                        three parameters: the Cairo context, the bounding
                        box of the drawing area and the palette to be
                        used for drawing colored vertices. The default
                        vertex drawer is L{DefaultVertexDrawer}.
        @param edge_drawer_factory: a factory method that returns an
                        L{AbstractEdgeDrawer} instance bound to a
                        given Cairo context. The factory method must take
                        two parameters: the Cairo context and the palette
                        to be used for drawing colored edges. You can use
                        any of the actual L{AbstractEdgeDrawer}
                        implementations here to control the style of
                        edges drawn by igraph. The default edge drawer is
                        L{ArrowEdgeDrawer}.
        @param label_drawer_factory: a factory method that returns a
                        L{TextDrawer} instance bound to a given Cairo
                        context. The method must take one parameter: the
                        Cairo context. The default label drawer is
                        L{TextDrawer}.
        """
        AbstractCairoGraphDrawer.__init__(self, context, bbox)
        self.vertex_drawer_factory = vertex_drawer_factory
        self.edge_drawer_factory = edge_drawer_factory
        self.label_drawer_factory = label_drawer_factory

    def _determine_edge_order(self, graph, kwds):
        """Returns the order in which the edge of the given graph have to be
        drawn, assuming that the relevant keyword arguments (C{edge_order} and
        C{edge_order_by}) are given in C{kwds} as a dictionary. If neither
        C{edge_order} nor C{edge_order_by} is present in C{kwds}, this
        function returns C{None} to indicate that the graph drawer is free to
        choose the most convenient edge ordering."""
        if "edge_order" in kwds:
            # Edge order specified explicitly
            return kwds["edge_order"]

        if kwds.get("edge_order_by") is None:
            # No edge order specified
            return None

        # Order edges by the value of some attribute
        edge_order_by = kwds["edge_order_by"]
        reverse = False
        if isinstance(edge_order_by, tuple):
            edge_order_by, reverse = edge_order_by
            if isinstance(reverse, basestring):
                reverse = reverse.lower().startswith("desc")
        attrs = graph.es[edge_order_by]
        edge_order = sorted(
            range(len(attrs)), key=attrs.__getitem__, reverse=bool(reverse))

        return edge_order

    def _determine_vertex_order(self, graph, kwds):
        """Returns the order in which the vertices of the given graph have to be
        drawn, assuming that the relevant keyword arguments (C{vertex_order} and
        C{vertex_order_by}) are given in C{kwds} as a dictionary. If neither
        C{vertex_order} nor C{vertex_order_by} is present in C{kwds}, this
        function returns C{None} to indicate that the graph drawer is free to
        choose the most convenient vertex ordering."""
        if "vertex_order" in kwds:
            # Vertex order specified explicitly
            return kwds["vertex_order"]

        if kwds.get("vertex_order_by") is None:
            # No vertex order specified
            return None

        # Order vertices by the value of some attribute
        vertex_order_by = kwds["vertex_order_by"]
        reverse = False
        if isinstance(vertex_order_by, tuple):
            vertex_order_by, reverse = vertex_order_by
            if isinstance(reverse, basestring):
                reverse = reverse.lower().startswith("desc")
        attrs = graph.vs[vertex_order_by]
        vertex_order = sorted(
            range(len(attrs)), key=attrs.__getitem__, reverse=bool(reverse))

        return vertex_order

    # pylint: disable-msg=W0142,W0221,E1101
    # W0142: Used * or ** magic
    # W0221: argument number differs from overridden method
    # E1101: Module 'cairo' has no 'foo' member - of course it does :)
    def draw(self, graph, palette, *args, **kwds):
        # Some abbreviations for sake of simplicity
        directed = graph.is_directed()
        context = self.context

        # Calculate/get the layout of the graph
        layout = self.ensure_layout(kwds.get("layout", None), graph)

        # Determine the size of the margin on each side
        margin = kwds.get("margin", 0)
        try:
            margin = list(margin)
        except TypeError:
            margin = [margin]
        while len(margin) < 4:
            margin.extend(margin)

        # Contract the drawing area by the margin and fit the layout
        bbox = self.bbox.contract(margin)
        layout.fit_into(
            bbox, keep_aspect_ratio=kwds.get("keep_aspect_ratio", False))

        # Decide whether we need to calculate the curvature of edges
        # automatically -- and calculate them if needed.
        autocurve = kwds.get("autocurve", None)
        if autocurve or (autocurve is None and "edge_curved" not in kwds and
                         "curved" not in graph.edge_attributes() and
                         graph.ecount() < 10000):
            from igraph import autocurve
            default = kwds.get("edge_curved", 0)
            if default is True:
                default = 0.5
            default = float(default)
            kwds["edge_curved"] = autocurve(
                graph, attribute=None, default=default)

        # Construct the vertex, edge and label drawers
        vertex_drawer = self.vertex_drawer_factory(context, bbox, palette,
                                                   layout)
        edge_drawer = self.edge_drawer_factory(context, palette)
        label_drawer = self.label_drawer_factory(context)

        # Construct the visual vertex/edge builders based on the specifications
        # provided by the vertex_drawer and the edge_drawer
        vertex_builder = vertex_drawer.VisualVertexBuilder(graph.vs, kwds)
        edge_builder = edge_drawer.VisualEdgeBuilder(graph.es, kwds)

        # Determine the order in which we will draw the vertices and edges
        vertex_order = self._determine_vertex_order(graph, kwds)
        edge_order = self._determine_edge_order(graph, kwds)

        # Draw the highlighted groups (if any)
        if "mark_groups" in kwds:
            mark_groups = kwds["mark_groups"]

            # Figure out what to do with mark_groups in order to be able to
            # iterate over it and get memberlist-color pairs
            if isinstance(mark_groups, dict):
                group_iter = mark_groups.iteritems()
            elif hasattr(mark_groups, "__iter__"):
                # Lists, tuples, iterators etc
                group_iter = iter(mark_groups)
            else:
                # False
                group_iter = {}.iteritems()

            # We will need a polygon drawer to draw the convex hulls
            polygon_drawer = PolygonDrawer(context, bbox)

            # Iterate over color-memberlist pairs
            for group, color_id in group_iter:
                if not group or color_id is None:
                    continue

                color = palette.get(color_id)

                if isinstance(group, VertexSeq):
                    group = [vertex.index for vertex in group]
                if not hasattr(group, "__iter__"):
                    raise TypeError("group membership list must be iterable")

                # Get the vertex indices that constitute the convex hull
                hull = [
                    group[i]
                    for i in convex_hull([layout[idx] for idx in group])
                ]

                # Calculate the preferred rounding radius for the corners
                corner_radius = 1.25 * \
                    max(vertex_builder[idx].size for idx in hull)

                # Construct the polygon
                polygon = [layout[idx] for idx in hull]

                if len(polygon) == 2:
                    # Expand the polygon (which is a flat line otherwise)
                    a, b = Point(*polygon[0]), Point(*polygon[1])
                    c = corner_radius * (a - b).normalized()
                    n = Point(-c[1], c[0])
                    polygon = [a + n, b + n, b - c, b - n, a - n, a + c]
                else:
                    # Expand the polygon around its center of mass
                    center = Point(*[
                        sum(coords) / float(len(coords))
                        for coords in zip(*polygon)
                    ])
                    polygon = [
                        Point(*point).towards(center, -corner_radius)
                        for point in polygon
                    ]

                # Draw the hull
                context.set_source_rgba(color[0], color[1], color[2],
                                        color[3] * 0.25)
                polygon_drawer.draw_path(polygon, corner_radius=corner_radius)
                context.fill_preserve()
                context.set_source_rgba(*color)
                context.stroke()

        # Construct the iterator that we will use to draw the edges
        es = graph.es
        if edge_order is None:
            # Default edge order
            edge_coord_iter = izip(es, edge_builder)
        else:
            # Specified edge order
            edge_coord_iter = ((es[i], edge_builder[i]) for i in edge_order)

        # Draw the edges
        if directed:
            drawer_method = edge_drawer.draw_directed_edge
        else:
            drawer_method = edge_drawer.draw_undirected_edge
        for edge, visual_edge in edge_coord_iter:
            src, dest = edge.tuple
            src_vertex, dest_vertex = vertex_builder[src], vertex_builder[dest]
            drawer_method(visual_edge, src_vertex, dest_vertex)

        # Construct the iterator that we will use to draw the vertices
        vs = graph.vs
        if vertex_order is None:
            # Default vertex order
            vertex_coord_iter = izip(vs, vertex_builder, layout)
        else:
            # Specified vertex order
            vertex_coord_iter = ((vs[i], vertex_builder[i], layout[i])
                                 for i in vertex_order)

        # Draw the vertices
        drawer_method = vertex_drawer.draw
        context.set_line_width(1)
        for vertex, visual_vertex, coords in vertex_coord_iter:
            drawer_method(visual_vertex, vertex, coords)

        # Set the font we will use to draw the labels
        vertex_label_family = 'sans-serif' if not hasattr(graph, "vertex_label_family") \
            else graph.vertex_label_family

        # Decide whether the labels have to be wrapped
        wrap = kwds.get("wrap_labels")
        if wrap is None:
            wrap = Configuration.instance()["plotting.wrap_labels"]
        wrap = bool(wrap)

        # Construct the iterator that we will use to draw the vertex labels
        if vertex_order is None:
            # Default vertex order
            vertex_coord_iter = izip(vertex_builder, layout)
        else:
            # Specified vertex order
            vertex_coord_iter = ((vertex_builder[i], layout[i])
                                 for i in vertex_order)

        # Draw the vertex labels
        for vertex, coords in vertex_coord_iter:
            if vertex.label is None:
                continue

            if hasattr(vertex, 'label_family'):
                context.select_font_face(vertex.label_family,
                                         cairo.FONT_SLANT_NORMAL,
                                         cairo.FONT_WEIGHT_NORMAL)
            context.set_font_size(vertex.label_size)
            context.set_source_rgba(*vertex.label_color)
            label_drawer.text = vertex.label

            if vertex.label_dist:
                # Label is displaced from the center of the vertex.
                _, yb, w, h, _, _ = label_drawer.text_extents()
                w, h = w / 2.0, h / 2.0
                radius = vertex.label_dist * vertex.size / 2.
                # First we find the reference point that is at distance `radius'
                # from the vertex in the direction given by `label_angle'.
                # Then we place the label in a way that the line connecting the
                # center of the bounding box of the label with the center of the
                # vertex goes through the reference point and the reference
                # point lies exactly on the bounding box of the vertex.
                alpha = vertex.label_angle % (2 * pi)
                cx = coords[0] + radius * cos(alpha)
                cy = coords[1] - radius * sin(alpha)
                # Now we have the reference point. We have to decide which side
                # of the label box will intersect with the line that connects
                # the center of the label with the center of the vertex.
                if w > 0:
                    beta = atan2(h, w) % (2 * pi)
                else:
                    beta = pi / 2.
                gamma = pi - beta
                if alpha > 2 * pi - beta or alpha <= beta:
                    # Intersection at left edge of label
                    cx += w
                    cy -= tan(alpha) * w
                elif alpha > beta and alpha <= gamma:
                    # Intersection at bottom edge of label
                    try:
                        cx += h / tan(alpha)
                    except:
                        pass  # tan(alpha) == inf
                    cy -= h
                elif alpha > gamma and alpha <= gamma + 2 * beta:
                    # Intersection at right edge of label
                    cx -= w
                    cy += tan(alpha) * w
                else:
                    # Intersection at top edge of label
                    try:
                        cx -= h / tan(alpha)
                    except:
                        pass  # tan(alpha) == inf
                    cy += h
                # Draw the label
                label_drawer.draw_at(cx - w, cy - h - yb, wrap=wrap)
            else:
                # Label is exactly in the center of the vertex
                cx, cy = coords
                half_size = vertex.size / 2.
                label_drawer.bbox = (cx - half_size, cy - half_size,
                                     cx + half_size, cy + half_size)
                label_drawer.draw(wrap=wrap)

        # Construct the iterator that we will use to draw the edge labels
        es = graph.es
        if edge_order is None:
            # Default edge order
            edge_coord_iter = izip(es, edge_builder)
        else:
            # Specified edge order
            edge_coord_iter = ((es[i], edge_builder[i]) for i in edge_order)

        # Draw the edge labels
        for edge, visual_edge in edge_coord_iter:
            if visual_edge.label is None:
                continue

            # Set the font size, color and text
            if hasattr(visual_edge, 'label_family'):
                context.select_font_face(visual_edge.label_family,
                                         cairo.FONT_SLANT_NORMAL,
                                         cairo.FONT_WEIGHT_NORMAL)
            context.set_font_size(visual_edge.label_size)
            context.set_source_rgba(*visual_edge.label_color)
            label_drawer.text = visual_edge.label

            # Ask the edge drawer to propose an anchor point for the label
            src, dest = edge.tuple
            src_vertex, dest_vertex = vertex_builder[src], vertex_builder[dest]
            (x, y), (halign, valign) = \
                edge_drawer.get_label_position(edge, src_vertex, dest_vertex)

            # Measure the text
            _, yb, w, h, _, _ = label_drawer.text_extents()
            w /= 2.0
            h /= 2.0

            # Place the text relative to the edge
            if halign == TextAlignment.RIGHT:
                x -= w
            elif halign == TextAlignment.LEFT:
                x += w
            if valign == TextAlignment.BOTTOM:
                y -= h - yb / 2.0
            elif valign == TextAlignment.TOP:
                y += h

            # Draw the edge label
            label_drawer.halign = halign
            label_drawer.valign = valign
            label_drawer.bbox = (x - w, y - h, x + w, y + h)
            label_drawer.draw(wrap=wrap)

#####################################################################
