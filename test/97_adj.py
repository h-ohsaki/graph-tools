#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True)
g.add_vertices(1, 2, 3)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(1, 3)

adj = g.adjacency_matrix()
print(adj)

adj.astype('float32').tofile('/tmp/np.dat')
# BinaryReadList["/tmp/np.dat", "Real32"]

print(g.diagonal_matrix())

