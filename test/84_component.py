#!/usr/bin/env python3

from test_more import ok
from graph_tools import Graph

g = Graph(directed=False)

g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_vertex(4)
ok(len(g.components()) == 2)
ok(len(g.maximal_component()) == 3)
