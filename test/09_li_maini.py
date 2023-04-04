#!/usr/bin/env python3

from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
g.create_li_maini_graph()
ok(g.is_connected())
