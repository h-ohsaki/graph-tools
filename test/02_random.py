#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False, multiedged=True)
g = g.create_graph('unknown')
ok(not g)

g = Graph(directed=True, multiedged=True)
g = g.create_graph('random', 4, 10)
ok(g)
eq(len(g.vertices()), 4)
eq(len(g.edges()), 10)
ok(g.is_connected())

g = Graph(directed=True, multiedged=True)
g = g.create_graph('random', 40, 100)
ok(g)
eq(len(g.vertices()), 40)
eq(len(g.edges()), 100)
ok(g.is_connected())

g = Graph(directed=True, multiedged=True)
g = g.create_graph('random', 4, 10)
ok(g)
eq(len(g.vertices()), 4)
eq(len(g.edges()), 10)
