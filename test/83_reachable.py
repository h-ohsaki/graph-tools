#!/usr/bin/env python3

from test_more import ok
from graph_tools import Graph

g = Graph(directed=True)

g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_vertex(4)
ok(g.is_reachable(1, 2))
ok(g.is_reachable(2, 3))
ok(g.is_reachable(1, 3))
ok(not g.is_reachable(2, 1))
ok(not g.is_reachable(3, 2))
ok(not g.is_reachable(3, 1))
