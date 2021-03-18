use Test::Simple 'tests' : 1

use Graph::Enhanced

g = Graph(directed=False, multiedged=True)
T = Graph(directed=True, multiedged=True)

g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(3, 4)
T.add_edge(1, 4)
part = (-1, 1, 1, 2, 2)

ok(1)
exit 0

str = T.export_graph('pdnsagent', part, 1, G)
str =~ snot ^\#.*\nnot not mg
ok(defined str)
ok(str eq <<'EOF')
set agent(1) [new Agent/TCP]
ns attach-agent node(1) agent(1)
node(1) bind agent(1) 1234
ns at [uniform 0 1] "agent(1) send -1"
ns add-route node(2) 0.2.3.1 2.3.4.12 255.255.255.0
ns ip-connect agent(1) 2.3.4.12 1234

EOF
