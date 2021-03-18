use Test::Simple 'tests' : 6

use Graph::Enhanced

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
str = g.export_graph('nsagent')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set tcp(1) [ns create-connection TCP node(1) TCPSink node(2) 1]
ns at [uniform 0 1] "tcp(1) send -1"
EOF

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
g.add_edge(2, 3)
str = g.export_graph('nsagent')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set tcp(1) [ns create-connection TCP node(1) TCPSink node(2) 1]
ns at [uniform 0 1] "tcp(1) send -1"
set tcp(2) [ns create-connection TCP node(2) TCPSink node(3) 2]
ns at [uniform 0 1] "tcp(2) send -1"
EOF

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
g.add_edge(1, 2)
str = g.export_graph('nsagent')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set tcp(1) [ns create-connection TCP node(1) TCPSink node(2) 1]
ns at [uniform 0 1] "tcp(1) send -1"
set tcp(2) [ns create-connection TCP node(1) TCPSink node(2) 2]
ns at [uniform 0 1] "tcp(2) send -1"
EOF
