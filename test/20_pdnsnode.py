use Test::Simple 'tests' : 2

use Graph::Enhanced

g = Graph(directed=False, multiedged=True)
T = Graph(directed=True, multiedged=True)

g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(3, 4)
T.add_edge(1, 4)
part = (-1, 1, 1, 2, 2)

str = g.export_graph('pdnsnode', part, 1)
str =~ snot ^\#.*\nnot not mg
ok(defined str)
ok(str eq <<'EOF')
set node(1) [ns node]
set node(2) [ns node]
ns duplex-link node(1) node(2) 1.5Mb 10ms DropTail
[ns link node(1) node(2)] set-ipaddr 1.1.2.11 255.255.255.0
[ns link node(2) node(1)] set-ipaddr 1.1.2.12 255.255.255.0
node(2) rlink 1.5Mb 10ms DropTail 0.2.3.1 255.255.255.0
EOF
