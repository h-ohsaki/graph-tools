use Test::Simple 'tests' : 6

use Graph::Enhanced

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
str = g.export_graph('nsagent_edge:2')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set node(3) [ns node]
ns duplex-link node(3) node(1) 100M 0.001ms DropTail
set node(4) [ns node]
ns duplex-link node(4) node(2) 100M 0.001ms DropTail
set tcp(1) [ns create-connection TCP node(3) TCPSink node(4) 1]
ns at [uniform 0 1] "tcp(1) send -1"

EOF

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
g.add_edge(2, 3)
str = g.export_graph('nsagent_edge:2')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set node(3) [ns node]
ns duplex-link node(3) node(1) 100M 0.001ms DropTail
set node(4) [ns node]
ns duplex-link node(4) node(2) 100M 0.001ms DropTail
set tcp(1) [ns create-connection TCP node(3) TCPSink node(4) 1]
ns at [uniform 0 1] "tcp(1) send -1"

set node(5) [ns node]
ns duplex-link node(5) node(2) 100M 0.001ms DropTail
set node(6) [ns node]
ns duplex-link node(6) node(3) 100M 0.001ms DropTail
set tcp(2) [ns create-connection TCP node(5) TCPSink node(6) 2]
ns at [uniform 0 1] "tcp(2) send -1"

EOF

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
g.add_edge(1, 2)
str = g.export_graph('nsagent_edge:2')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set node(3) [ns node]
ns duplex-link node(3) node(1) 100M 0.001ms DropTail
set node(4) [ns node]
ns duplex-link node(4) node(2) 100M 0.001ms DropTail
set tcp(1) [ns create-connection TCP node(3) TCPSink node(4) 1]
ns at [uniform 0 1] "tcp(1) send -1"

set node(5) [ns node]
ns duplex-link node(5) node(1) 100M 0.001ms DropTail
set node(6) [ns node]
ns duplex-link node(6) node(2) 100M 0.001ms DropTail
set tcp(2) [ns create-connection TCP node(5) TCPSink node(6) 2]
ns at [uniform 0 1] "tcp(2) send -1"

EOF
