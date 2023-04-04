#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
g = g.create_graph('barabasi', 10, 2, 2)
ok(g)
eq(g.nvertices(), 10)
eq(g.nedges(), 2 * (2 - 1) / 2 + (10 - 2) * 2)
ok(g.is_connected())

g = Graph(directed=False)
g = g.create_graph('barabasi', 100, 3, 4)
ok(g)
eq(g.nvertices(), 100)
eq(g.nedges(), 3 * (3 - 1) / 2 + (100 - 3) * 4)
ok(g.is_connected())

g = Graph(directed=False)
g = g.create_graph('ba', 10, 2, 2)
ok(g)
eq(g.nvertices(), 10)
eq(g.nedges(), 2 * (2 - 1) / 2 + (10 - 2) * 2)
ok(g.is_connected())
