#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

"""
1 -- 2 -- 3
"""

g = Graph(directed=False)
g.add_edge(1, 2)
g.add_edge(2, 3)
eq(g.nvertices(), 3)
eq(g.nedges(), 2)

g.merge_vertices(1, 2)
eq(g.nvertices(), 2)
eq(g.nedges(), 1)
eq(g.get_edge_weight(1, 3), None)

"""
         3
       /
1 -- 2   |
       \ 
         4
"""

g = Graph(directed=False)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(2, 4)
g.add_edge(3, 4)
eq(g.nvertices(), 4)
eq(g.nedges(), 4)

g.merge_vertices(2, 3)
eq(g.nvertices(), 3)
eq(g.nedges(), 2)
eq(g.get_edge_weight(1, 2), None)
eq(g.get_edge_weight(2, 4), 2)

g = Graph(directed=False)
g.add_edge(1, 2, weight=3)
g.add_edge(2, 3, weight=5)
g.add_edge(2, 4, weight=6)
g.add_edge(3, 4, weight=7)
eq(g.nvertices(), 4)
eq(g.nedges(), 4)

g.merge_vertices(2, 3)
eq(g.nvertices(), 3)
eq(g.nedges(), 2)
eq(g.get_edge_weight(1, 2), 3)
eq(g.get_edge_weight(2, 4), 6 + 7)
