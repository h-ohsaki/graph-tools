#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
ok(g)
ok(not g.directed())
ok(g.multiedged())

g = Graph(directed=True)
ok(g)
ok(g.directed())
ok(g.multiedged())
