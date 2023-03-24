#!/usr/bin/env python3

from test_more import ok
from graph_tools import Graph

g = Graph(directed=True)

g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_vertex(4)
dist, prev = g.dijkstra(1)
ok(dist[1] == 0)
ok(dist[2] == 1)
ok(dist[3] == 2)
ok(not dist.get(4, None))
ok(prev[1] == [])
ok(prev[2] == [1])
ok(prev[3] == [2])
ok(prev[4] == [])
ok(list(g.shortest_paths(1, 3)) == [[1, 2, 3]])
ok(list(g.shortest_paths(1, 4)) == [])

g = Graph(directed=True)

g.add_edge(1, 2)
g.add_edge(2, 4)
g.add_edge(1, 3)
g.add_edge(3, 4)
dist, prev = g.dijkstra(1)
ok(dist[1] == 0)
ok(dist[2] == 1)
ok(dist[3] == 1)
ok(dist[4] == 2)
ok(prev[1] == [])
ok(prev[2] == [1])
ok(prev[3] == [1])
ok(set(prev[4]) == set([2, 3]))
ok(list(g.shortest_paths(1, 4)) == [[1, 2, 4], [1, 3, 4]])

g = Graph(directed=True)

g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(1, 3)
dist, prev = g.dijkstra(1)
ok(dist[1] == 0)
ok(dist[2] == 1)
ok(dist[3] == 1)
ok(prev[1] == [])
ok(prev[2] == [1])
ok(prev[3] == [1])
ok(list(g.shortest_paths(1, 3)) == [[1, 3]])

g = Graph(directed=True)

g.add_edge(1, 2)
g.set_edge_weight(1, 2, .5)
g.add_edge(2, 4)
g.set_edge_weight(2, 4, 2)
g.add_edge(1, 3)
g.add_edge(3, 4)
dist, prev = g.dijkstra(1)
ok(dist[1] == 0)
ok(dist[2] == .5)
ok(dist[3] == 1)
ok(dist[4] == 2)
ok(prev[1] == [])
ok(prev[2] == [1])
ok(prev[3] == [1])
ok(prev[4] == [3])
ok(list(g.shortest_paths(1, 4)) == [[1, 3, 4]])
