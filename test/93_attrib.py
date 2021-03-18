#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True)
g.add_vertices(1, 2, 3)
g.add_edge(1, 2)
g.add_edge(1, 2)
g.add_edge(2, 3)

name, val, val2 = 'foo', 'bar', 'buzz'
g.set_graph_attribute(name, val)
eq(g.get_graph_attribute(name), val)

g.set_vertex_attribute(1, name, val)
eq(g.get_vertex_attribute(1, name), val)

g.set_edge_attribute_by_id(1, 2, 0, name, val)
eq(g.get_edge_attribute_by_id(1, 2, 0, name), val)

g.set_edge_attribute_by_id(1, 2, 1, name, val2)
eq(g.get_edge_attribute_by_id(1, 2, 1, name), val2)

g.set_edge_weight_by_id(1, 2, 0, val)
eq(g.get_edge_weight_by_id(1, 2, 0), val)

g.set_edge_weight_by_id(1, 2, 1, val2)
eq(g.get_edge_weight_by_id(1, 2, 1), val2)
