#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
g.add_vertices(1, 2, 3, 4)
g.complete_graph()

eq(g.spectral_radius(), 3)
# print(g.spectral_gap())
# print(g.natural_connectivity())
# print(g.algebraic_connectivity())
eq(g.effective_resistance(), 3)
# print(g.spanning_tree_count())
