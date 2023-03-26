#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True, multiedged=True)
g = g.create_graph('tree', 10)
ok(g)
eq(g.nvertices(), 10)
eq(g.nedges(), 10 - 1)

g = Graph(directed=True, multiedged=True)
g = g.create_graph('tree', 100)
ok(g)
eq(g.nvertices(), 100)
eq(g.nedges(), 100 - 1)

g = Graph(directed=False, multiedged=True)
g = g.create_graph('tree', 10)
ok(g)
eq(g.nvertices(), 10)
eq(g.nedges(), 10 - 1)
ok(g.is_connected())
