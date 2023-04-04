#!/usr/bin/env python3

import re
from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False, multiedged=True)
g.add_vertex(1)
astr = g.export_graph('dot')
ok(astr)
astr = re.sub(r'^//.*\n', '', astr, flags=re.M)
eq(astr, """graph export_dot {
  node [color=gray90,style=filled];
  "1";
}
""")

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
astr = g.export_graph('dot')
ok(astr)
astr = re.sub(r'^//.*\n', '', astr, flags=re.M)
eq(astr, """graph export_dot {
  node [color=gray90,style=filled];
  "1";
  "2";
  "1" -- "2";
}
""")

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
astr = g.export_graph('dot')
ok(astr)
astr = re.sub(r'^//.*\n', '', astr, flags=re.M)
eq(astr, """digraph export_dot {
  node [color=gray90,style=filled];
  "1";
  "2";
  "1" -> "2";
}
""")

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
g.add_edge(2, 1)
astr = g.export_graph('dot')
ok(astr)
astr = re.sub(r'^//.*\n', '', astr, flags=re.M)
eq(astr, """digraph export_dot {
  node [color=gray90,style=filled];
  "1";
  "2";
  "1" -> "2";
  "2" -> "1";
}
""")

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
g.set_edge_weight_by_id(1, 2, 0, 123)
astr = g.export_graph('dot')
ok(astr)
astr = re.sub(r'^//.*\n', '', astr, flags=re.M)
eq(astr, """graph export_dot {
  node [color=gray90,style=filled];
  "1";
  "2";
  "1" -- "2" [weight=123];
}
""")

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
g.add_edge(1, 2)
astr = g.export_graph('dot')
ok(astr)
astr = re.sub(r'^//.*\n', '', astr, flags=re.M)
eq(astr, """digraph export_dot {
  node [color=gray90,style=filled];
  "1";
  "2";
  "1" -> "2";
  "1" -> "2";
}
""")
