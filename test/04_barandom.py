#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False, multiedged=True)
g = g.create_graph('barandom', 10, 20, 2)
ok(g)
eq(len(g.vertices()), 10)
ok(g.is_connected())

g = Graph(directed=False, multiedged=True)
g = g.create_graph('barandom', 100, 200, 2)
ok(g)
eq(len(g.vertices()), 100)
ok(g.is_connected())
