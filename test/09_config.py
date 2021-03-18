#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
g.create_configuration_graph([3, 2, 2, 1, 1, 1, 0])
ok(g.degree(1) == 3)
ok(g.degree(2) == 2)
ok(g.degree(3) == 2)
ok(g.degree(4) == 1)
ok(g.degree(5) == 1)
ok(g.degree(6) == 1)
ok(g.degree(7) == 0)

g = Graph(directed=False)
g.create_configuration_graph([3, 3, 3, 3, 3, 3])
for v in g.vertices():
    ok(g.degree(v) == 3)

g = Graph(directed=False)
g.create_random_regular_graph(10, 3)
for v in g.vertices():
    ok(g.degree(v) == 3)
