#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
ok(g)
g.add_vertices(1, 2, 3, 4)

h = g.complete_graph()
ok(h.is_connected())
eq(len(h.vertices()), 4)
eq(len(h.edges()), 6)
eq(h.vertex_degree(1), 3)
