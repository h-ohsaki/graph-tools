#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True)
g.add_vertices(1, 2, 3)
g.add_edge(1, 2)
g.add_edge(1, 2)
g.add_edge(2, 3)

eq(len(g.neighbors(1)), 1)
eq(len(g.neighbors(2)), 2)
eq(len(g.neighbors(3)), 1)

g.delete_vertex(2)

eq(len(g.neighbors(1)), 0)
eq(len(g.neighbors(3)), 0)

g = Graph(directed=False)
g.add_vertices(1, 2, 3)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.delete_vertex(1)
eq(len(g.vertices()), 2)
eq(len(g.neighbors(2)), 1)
eq(len(g.neighbors(3)), 1)
