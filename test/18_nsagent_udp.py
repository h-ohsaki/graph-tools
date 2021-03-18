use Test::Simple 'tests' : 4

use Graph::Enhanced

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
str = g.export_graph('nsagent_udp')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set udp(1) [ns create-connection UDP node(1) UDP node(2) 1]
set cbr(1) [new Application/Traffic/CBR]
cbr(1) attach-agent udp(1)
ns at [uniform 0 1] "cbr(1) start"

EOF

g = Graph(directed=True, multiedged=True)
g.add_edge(1, 2)
g.add_edge(1, 2)
str = g.export_graph('nsagent_udp')
ok(defined str)
str =~ snot ^\#.*\nnot not mg
ok(str eq <<'EOF')
set udp(1) [ns create-connection UDP node(1) UDP node(2) 1]
set cbr(1) [new Application/Traffic/CBR]
cbr(1) attach-agent udp(1)
ns at [uniform 0 1] "cbr(1) start"

set udp(2) [ns create-connection UDP node(1) UDP node(2) 2]
set cbr(2) [new Application/Traffic/CBR]
cbr(2) attach-agent udp(2)
ns at [uniform 0 1] "cbr(2) start"

EOF
