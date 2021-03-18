use Test::Simple 'tests' : 6

use Graph::Enhanced

g = Graph(directed=True, multiedged=True)
buf = split(/\n/, <<'EOF')
# create TCP connections
set tcp(1) [ns create-connection TCP node(1) TCPSink node(2) 1]
# initiate TCP connections after random delay
for _ in { set i 1 } { i <= 1 } { incr i } :
  ns at [uniform 0 1] "tcp(i) send -1"

EOF
g.import_graph('nsagent', buf)
ok(len(g.vertices()) == 2)
ok(len(g.edges()) == 1)
ok(g.has_edge(1, 2))

g = Graph(directed=True, multiedged=True)
buf = split(/\n/, <<'EOF')
set node(1) [ns node]
set node(2) [ns node]
set agent(1) [new Agent/TCP]
set agent(2) [new Agent/TCPSink]
ns attach-agent node(1) agent(1)
ns attach-agent node(2) agent(2)
ns connect agent(1) agent(2)
EOF
g.import_graph('nsagent', buf)
ok(len(g.vertices()) == 2)
ok(len(g.edges()) == 1)
ok(g.has_edge(1, 2))
