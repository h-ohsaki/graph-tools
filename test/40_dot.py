#!/usr/bin/env python3

import re
from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=True)
buf = """digraph sample {
  1;
  2;
  1 -> 2;
}
""".splitlines()
g.import_graph('dot', buf)
eq(g.nvertices(), 2)
eq(g.nedges(), 1)
ok(g.is_directed())
ok(g.has_edge(1, 2))

g = Graph(directed=True)
buf = """// comment here
digraph sample {
  1;
/* another comment
   here */
  2;
  4;
  1 -> 2;
  1 -> 4;
}
""".splitlines()
g.import_graph('dot', buf)
eq(g.nvertices(), 3)
eq(g.nedges(), 2)
ok(g.is_directed())
ok(g.has_edge(1, 2))
ok(g.has_edge(1, 4))

g = Graph(directed=True, multiedged=True)
buf = """// comment here
digraph sample {
  1 [color=yellow];
  2;
  1 -> 2 [bw="1.5Mb",delay="10ms",
          type="RED"];
  1 -> 2;
}
""".splitlines()
g.import_graph('dot', buf)
eq(g.nvertices(), 2)
eq(g.nedges(), 2)
eq(g.is_directed(), 1)
ok(g.has_edge(1, 2))
eq(g.get_vertex_attribute(1, 'color'), 'yellow')
eq(g.get_edge_attribute_by_id(1, 2, 0, 'bw'), '1.5Mb')
eq(g.get_edge_attribute_by_id(1, 2, 0, 'delay'), '10ms')
eq(g.get_edge_attribute_by_id(1, 2, 0, 'type'), 'RED')
