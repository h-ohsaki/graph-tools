#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(3, 1)
g.set_vertex_attribute(1, 'foo', 123)
g.set_edge_attribute_by_id(1, 2, 0, 'bar', 456)
T = g.copy_graph()
ok(T.directed())
ok(T.multiedged())
ok(T.has_edge(1, 2))
ok(not T.has_edge(2, 1))
eq(T.get_vertex_attribute(1, 'foo'), 123)
eq(T.get_edge_attribute_by_id(1, 2, 0, 'bar'), 456)

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(3, 1)
g.set_vertex_attribute(1, 'foo', 123)
g.set_edge_attribute_by_id(1, 2, 0, 'bar', 456)
T = g.copy_graph()
ok(T.undirected())
ok(T.multiedged())
ok(T.has_edge(1, 2))
ok(T.has_edge(2, 1))
eq(T.get_vertex_attribute(1, 'foo'), 123)
eq(T.get_edge_attribute_by_id(1, 2, 0, 'bar'), 456)
