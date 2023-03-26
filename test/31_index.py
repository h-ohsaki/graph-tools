#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_vertex(4)
ok(g.vertex_index(1) == 0)
ok(g.vertex_index(2) == 1)
ok(g.vertex_index(3) == 2)
ok(g.vertex_index(4) == 3)

g.delete_vertex(1)
ok(g.vertex_index(1) == None)
ok(g.vertex_index(2) == 0)
ok(g.vertex_index(3) == 1)
ok(g.vertex_index(4) == 2)
