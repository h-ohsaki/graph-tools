#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True)
g.add_vertices(1, 2, 3)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(1, 3)

adj = g.adjacency_matrix()
ok(list(adj[0]) == [0, 1, 1])
ok(list(adj[1]) == [0, 0, 1])
ok(list(adj[2]) == [0, 0, 0])

diag = g.diagonal_matrix()

ok(diag[0][0] == 2)
ok(diag[1][1] == 2)
ok(diag[2][2] == 2)

g = Graph(directed=False)
g.add_vertices(1, 2, 3)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(1, 3)

adj = g.adjacency_matrix()
ok(list(adj[0]) == [0, 1, 1])
ok(list(adj[1]) == [1, 0, 1])
ok(list(adj[2]) == [1, 1, 0])

diag = g.diagonal_matrix()

ok(diag[0][0] == 2)
ok(diag[1][1] == 2)
ok(diag[2][2] == 2)
