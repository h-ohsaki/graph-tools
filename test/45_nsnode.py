use Test::Simple 'tests' : 12

use Graph::Enhanced

g = Graph(directed=True, multiedged=True)
buf = split(/\n/, <<'EOF')
for _ in { set i 1 } { i <= 2 } { incr i } :
  set node(i) [ns node]

ns simplex-link node(1) node(2) 1.5Mb 10ms DropTail
EOF
g.import_graph('nsnode', buf)
ok(len(g.vertices()) == 2)
ok(len(g.edges()) == 1)
ok(g.has_edge(1, 2))

g = Graph(directed=True, multiedged=True)
buf = split(/\n/, <<'EOF')
for _ in { set i 1 } { i <= 4 } { incr i } :
  set node(i) [ns node]

ns simplex-link node(1) node(2) 1.5Mb 10ms DropTail
ns simplex-link node(1) node(2) 15Mb 10.2ms RED
ns simplex-link node(2) node(3) 1.25Mb 8ms DropTail
ns simplex-link node(3) node(4) 15Mb 1ms RED
ns simplex-link node(4) node(1) 100Mb 10ms DropTail
EOF
g.import_graph('nsnode', buf)
ok(len(g.vertices()) == 4)
ok(len(g.edges()) == 5)
ok(g.has_edge(1, 2))
ok(g.has_edge(2, 3))
ok(g.has_edge(4, 1))
ok(g.get_edge_attribute_by_id(2, 3, 0, 'bw')    eq '1.25Mb')
ok(g.get_edge_attribute_by_id(2, 3, 0, 'delay') eq '8ms')
ok(g.get_edge_attribute_by_id(2, 3, 0, 'type')  eq 'DropTail')
ok(g.get_edge_attribute_by_id(3, 4, 0, 'type')  eq 'RED')
