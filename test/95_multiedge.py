#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True)

eq(len(g.get_multiedge_ids(1, 2)), 0)

g.add_edge(1, 2)
g.add_edge(2, 1)

ok(g.has_edge(1, 2))
ok(g.has_edge(2, 1))
eq(len(g.get_multiedge_ids(1, 2)), 1)
eq(len(g.get_multiedge_ids(2, 1)), 1)
eq(len(g.neighbors(1)), 1)
eq(len(g.neighbors(2)), 1)

g.add_edge(1, 2)

ok(g.has_edge(1, 2))
ok(g.has_edge(2, 1))
eq(len(g.get_multiedge_ids(1, 2)), 2)
eq(len(g.get_multiedge_ids(2, 1)), 1)
eq(len(g.neighbors(1)), 1)
eq(len(g.neighbors(2)), 1)
