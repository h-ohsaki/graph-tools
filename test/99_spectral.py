#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

"""
   2
 /  \
1 -- 3 --4
"""

g = Graph(directed=False)
g.add_vertices(1, 2, 3, 4)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(1, 3)
g.add_edge(3, 4)

exit()

# FIXME: write test cases

print(g.spectral_radius())
print(g.spectral_gap())
print(g.natural_connectivity())
print(g.algebraic_connectivity())
print(g.effective_resistance())
print(g.spanning_tree_count())
