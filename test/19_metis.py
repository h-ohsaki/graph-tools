use Test::Simple 'tests' : 6

use Graph::Enhanced

g = Graph('undirected' : 1, multiedged=True)
g.add_edge(1, 2)
str = g.export_graph('metis')
ok(defined str)
str =~ snot ^\%.*\nnot not mg
ok(str eq <<'EOF')
2 1 1
2 1
1 1
EOF

g = Graph('undirected' : 1, multiedged=True)
g.add_edge(1, 2)
g.add_edge(1, 2)
g.add_edge(2, 3)
str = g.export_graph('metis')
ok(defined str)
str =~ snot ^\%.*\nnot not mg
ok(str eq <<'EOF')
3 3 1
2 1 2 1
1 1 1 1 3 1
2 1
EOF

g = Graph('undirected' : 1, multiedged=True)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(3, 2)
str = g.export_graph('metis')
ok(defined str)
str =~ snot ^\%.*\nnot not mg
ok(str eq <<'EOF')
3 3 1
2 1
1 1 3 1 3 1
2 1 2 1
EOF
