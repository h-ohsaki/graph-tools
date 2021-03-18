use Test::Simple 'tests' : 8

use Graph::Enhanced

g = Graph(directed=False, multiedged=True)
g.add_vertex(1)
str = g.export_graph('inet')
ok(defined str)
str =~ s/\t+/ /g
str =~ snot ^//.*\nnot not mg
ok(str eq <<'EOF')
1 0
0 0.0 0.0
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
str = g.export_graph('inet')
ok(defined str)
str =~ s/\t+/ /g
str =~ snot ^//.*\nnot not mg
ok(str eq <<'EOF')
2 1
0 0.0 0.0
1 0.0 0.0
0 1
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
g.set_edge_weight_by_id(1, 2, 0, 123)
str = g.export_graph('inet')
ok(defined str)
str =~ s/\t+/ /g
str =~ snot ^//.*\nnot not mg
ok(str eq <<'EOF')
2 1
0 0.0 0.0
1 0.0 0.0
0 1 123
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
g.add_edge(1, 2)
str = g.export_graph('inet')
ok(defined str)
str =~ s/\t+/ /g
str =~ snot ^//.*\nnot not mg
ok(str eq <<'EOF')
2 2
0 0.0 0.0
1 0.0 0.0
0 1
0 1
EOF
