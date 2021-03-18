#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
ok(g)
ok(not g.is_directed())
ok(g.is_undirected())

g.add_vertices(1, 2, 3, 4)
ok(g.has_vertex(1))
ok(g.has_vertex(2))
ok(g.has_vertex(3))
ok(g.has_vertex(4))
ok(not g.has_vertex(5))
ok(not g.has_vertex(0))

g.add_edge(1, 2)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(3, 4)
ok(g.has_edge(1, 2))
ok(g.has_edge(2, 3))
ok(g.has_edge(2, 1))
ok(not g.has_edge(1, 3))
ok(not g.has_edge(4, 1))

eq(len(g.vertices()), 4)
eq(len(g.neighbors(1)), 1)
eq(len(g.neighbors(2)), 2)

eq(len(g.edges()), 4)
eq(len(g.unique_edges()), 3)
eq(len(g.edges_at(1)), 2)
eq(len(g.edges_at(2)), 3)

eq(g.degree(1), 2)
eq(g.degree(2), 3)
