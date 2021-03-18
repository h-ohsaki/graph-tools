#!/usr/bin/env python3

from test_more import ok
from graph_tools import Graph

g = Graph(directed=True)

g.add_edge(1, 2)
g.add_edge(2, 3)
g.floyd_warshall()
ok(g.T[1][2] == 1)
ok(g.T[2][3] == 1)
ok(g.T[1][3] == 2)
ok(not g.T[2].get(1, None))
ok(not g.T[3].get(2, None))
ok(not g.T[3].get(1, None))
