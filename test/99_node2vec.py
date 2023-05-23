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

vecs = g.node2vecs(dimension=5)

for ui in range(4):
    emb = vecs[ui]
    print(f'embedding for node {ui+1}: e_{ui+1} = {emb}')
    for vi in range(4):
        norm = numpy.linalg.norm(vecs[ui] - vecs[vi], ord=1)
        print(f'  |e_{ui+1} - e_{vi+1}|_1 = {norm:.2f}')
