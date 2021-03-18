use Test::Simple 'tests' : 14

use Graph::Enhanced

g = Graph(directed=True, multiedged=True)
buf = split(/\n/, <<EOF)
p sample 2 1
a 1 2 1
EOF
g.import_graph('dimacs', buf)
ok(len(g.vertices()) == 2)
ok(len(g.edges()) == 1)
ok(g.is_directed == 1)
ok(g.has_edge(1, 2))

g = Graph(directed=True, multiedged=True)
buf = split(/\n/, <<EOF)
p dot2dimacs 2 4
a 1 2 1
a 1 2 1
a 2 1 1
a 2 1 1
EOF
g.import_graph('dimacs', buf)
ok(len(g.vertices()) == 2)
ok(len(g.edges()) == 4)
ok(g.is_directed == 1)
ok(g.has_edge(1, 2))
ok(g.has_edge(2, 1))

g = Graph(directed=True, multiedged=True)
buf = split(/\n/, <<EOF)
p sample 2 1
a 1 2 123.456
EOF
g.import_graph('dimacs', buf)
ok(len(g.vertices()) == 2)
ok(len(g.edges()) == 1)
ok(g.is_directed == 1)
ok(g.has_edge(1, 2))
ok(g.get_edge_weight_by_id(1, 2, 0) eq '123.456')
