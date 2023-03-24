#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
g.add_vertices(1, 2, 3)
g.add_edge(1, 2)
g.add_edge(2, 3)
ok(g.betweenness(1) == 0)
ok(g.betweenness(2) == 2)
ok(g.betweenness(3) == 0)


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
ok(g.betweenness(1) == 0)
ok(g.betweenness(2) == 0)
ok(g.betweenness(3) == 4)
ok(g.betweenness(4) == 0)

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
ok(g.betweenness(1) == 0)
ok(g.betweenness(2) == 0)
ok(g.betweenness(3) == 4)
ok(g.betweenness(4) == 0)
