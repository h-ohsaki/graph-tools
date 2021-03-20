use Test::Simple 'tests' : 12

use Graph::Enhanced

g = Graph(directed=False, multiedged=True)
g.add_vertex(1)
str = g.export_graph('nsnode')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set node(1) [ns node]
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
str = g.export_graph('nsnode')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set node(1) [ns node]
set node(2) [ns node]
ns duplex-link node(1) node(2) 1.5Mb 10ms DropTail
EOF

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
g.add_edge(2, 3)
str = g.export_graph('nsnode')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set node(1) [ns node]
set node(2) [ns node]
set node(3) [ns node]
ns simplex-link node(1) node(2) 1.5Mb 10ms DropTail
ns simplex-link node(2) node(3) 1.5Mb 10ms DropTail
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.set_edge_attribute_by_id(1, 2, 0, 'bw',    '123Mb')
g.set_edge_attribute_by_id(1, 2, 0, 'delay', '456ms')
g.set_edge_attribute_by_id(2, 3, 0, 'type',  'RED')
str = g.export_graph('nsnode')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set node(1) [ns node]
set node(2) [ns node]
set node(3) [ns node]
ns duplex-link node(1) node(2) 123Mb 456ms DropTail
ns duplex-link node(2) node(3) 1.5Mb 10ms RED
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
g.add_edge(1, 2)
str = g.export_graph('nsnode')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set node(1) [ns node]
set node(2) [ns node]
ns duplex-link node(1) node(2) 1.5Mb 10ms DropTail
ns duplex-link node(1) node(2) 1.5Mb 10ms DropTail
EOF

g = Graph(directed=False, multiedged=True)
g.add_edge(1, 2)
g.add_edge(2, 1)
str = g.export_graph('nsnode')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set node(1) [ns node]
set node(2) [ns node]
ns duplex-link node(1) node(2) 1.5Mb 10ms DropTail
ns duplex-link node(1) node(2) 1.5Mb 10ms DropTail
EOF