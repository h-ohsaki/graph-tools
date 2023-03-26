#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True, multiedged=True)
g = g.create_graph('ring', 10, 1)
ok(g)
eq(g.nvertices(), 10)
eq(g.nedges(), 10)
for v in range(1, 11):
    eq(g.in_degree(v), 1)
    eq(g.out_degree(v), 1)

g = Graph(directed=False, multiedged=True)
g = g.create_graph('ring', 10, 1)
ok(g)
eq(g.nvertices(), 10)
eq(g.nedges(), 10)
for v in range(1, 11):
    eq(g.degree(v), 2)
ok(g.is_connected())
