#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True)
g = g.create_graph('btree', 15)
ok(g)
eq(g.nvertices(), 15)
eq(g.nedges(), 15 - 1)

g = Graph(directed=True)
g = g.create_graph('btree', 100)
ok(g)
eq(g.nvertices(), 100)
eq(g.nedges(), 100 - 1)

g = Graph(directed=False)
g = g.create_graph('btree', 10)
ok(g)
eq(g.nvertices(), 10)
eq(g.nedges(), 10 - 1)
ok(g.is_connected())
