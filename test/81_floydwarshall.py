#!/usr/bin/env python3

import math

from test_more import ok
from graph_tools import Graph

g = Graph(directed=True)

g.add_edge(1, 2)
g.add_edge(2, 3)
g.floyd_warshall()
ok(g.T[1][2] == 1)
ok(g.T[2][3] == 1)
ok(g.T[1][3] == 2)
ok(g.T[2][1] == math.inf)
ok(g.T[3][2] == math.inf)
ok(g.T[3][1] == math.inf)

g = Graph(directed=False)

g.add_edge(1, 2)
g.add_edge(2, 3)
g.floyd_warshall()
ok(g.T[1][2] == 1)
ok(g.T[2][3] == 1)
ok(g.T[1][3] == 2)
ok(g.T[2][1] == 1)
ok(g.T[3][2] == 1)
ok(g.T[3][1] == 2)
