#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

"""
1 -- 2 -- 3
"""
g = Graph(directed=False)
g.add_vertices(1, 2, 3)
g.add_edge(1, 2)
g.add_edge(2, 3)
n = g.nvertices()
ok(g.betweenness_centrality(1) == 0)
ok(g.betweenness_centrality(2) == 2 / (n - 1) / (n - 2))
ok(g.betweenness_centrality(3) == 0)


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
n = g.nvertices()
ok(g.betweenness_centrality(1) == 0)
ok(g.betweenness_centrality(2) == 0)
ok(g.betweenness_centrality(3) == 4 / (n - 1) / (n - 2))
ok(g.betweenness_centrality(4) == 0)

"""
   2
 // \
1 -- 3 --4
"""

g = Graph(directed=False)
g.add_vertices(1, 2, 3, 4)
g.add_edge(1, 2)
g.add_edge(1, 2) # multiedge
g.add_edge(2, 3)
g.add_edge(1, 3)
g.add_edge(3, 4)
n = g.nvertices()
ok(g.betweenness_centrality(1) == 0)
ok(g.betweenness_centrality(2) == 0)
ok(g.betweenness_centrality(3) == 4 / (n - 1) / (n - 2))
ok(g.betweenness_centrality(4) == 0)
