#!/usr/bin/env python3

import numpy
from test_more import ok, eq
from graph_tools import Graph

g = Graph(directed=False)
g.add_vertices(1, 2, 3, 4)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(1, 3)
g.add_edge(3, 4)

for v in [1, 2, 3, 4]:
    print(g.node2vec(v))

for u in [1, 2, 3, 4]:
    e_u = g.node2vec(u)
    print(f'embedding for node {u}: e_{u} = {e_u}')
    for v in [1, 2, 3, 4]:
        e_v = g.node2vec(v)
        norm = numpy.linalg.norm(e_u - e_v, ord=1)
        print(f'  |e_{u} - e_{v}|_1 = {norm:.2f}')
