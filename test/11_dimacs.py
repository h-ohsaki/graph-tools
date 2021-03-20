use Test::Simple 'tests' : 8

use Graph::Enhanced

g = Graph(directed=False, multiedged=True)
g.add_vertex(1)
str = g.export_graph('dimacs')
ok(defined str)
str =~ snot ^//.*\nnot not mg
ok(str eq <<'EOF')
p dot2dimacs 1 0
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
str = g.export_graph('dimacs')
ok(defined str)
str =~ snot ^//.*\nnot not mg
ok(str eq <<'EOF')
p dot2dimacs 2 1
a 1 2 1
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
g.set_edge_weight_by_id(1, 2, 0, 123)
str = g.export_graph('dimacs')
ok(defined str)
str =~ snot ^//.*\nnot not mg
ok(str eq <<'EOF')
p dot2dimacs 2 1
a 1 2 123
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
g.add_edge(1, 2)
str = g.export_graph('dimacs')
ok(defined str)
str =~ snot ^//.*\nnot not mg
ok(str eq <<'EOF')
p dot2dimacs 2 2
a 1 2 1
a 1 2 1
EOF