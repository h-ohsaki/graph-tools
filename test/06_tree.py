#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True, multiedged=True)
g = g.create_graph('tree', 10)
ok(g)
eq(len(g.vertices()), 10)
eq(len(g.edges()), 10 - 1)

g = Graph(directed=True, multiedged=True)
g = g.create_graph('tree', 100)
ok(g)
eq(len(g.vertices()), 100)
eq(len(g.edges()), 100 - 1)

g = Graph(directed=False, multiedged=True)
g = g.create_graph('tree', 10)
ok(g)
eq(len(g.vertices()), 10)
eq(len(g.edges()), 10 - 1)
ok(g.is_connected())
