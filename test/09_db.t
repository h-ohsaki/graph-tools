#!/usr/bin/env python3

from test_more import ok, eq
from graphtools import Graph

g = Graph(directed=True)
g.create_degree_bounded_graph(10, 20)
for v in g.vertices():
    ok(g.degree(v) >= 2)

g = Graph(directed=True)
g.create_degree_bounded_graph(100, 200)
for v in g.vertices():
    ok(g.degree(v) >= 2)
