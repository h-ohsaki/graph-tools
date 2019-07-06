# grapht-tools Package

graph_tools - tools for graph theory and network science with many generation models

# DESCRIPTION

This manual page documents **graph-tools** module, a Python module that
provides a number of features for handling directed/undirected graphs and
complex networks.  **graph-tools** was initially developed for networking
researchers, who perform experiments in the field of graph theory and network
science.  **graph-tools** provides Graph class, which supports both directed
and undirected graphs with multi-edges, vertex weights, edge weights, and
graph attributes.  A number of graph/network generation models and graph
algorithms are supported.

Major features of **graph-tools* are:

- directed/undirected graph with multi-edges, vertex weights, edge weights,
  and graph attributes
  
- vertex operations (add, delete, degree, neighbors, random vertex, and
  set/get vertex attributes)

- edge operations (add, delete, random edge, and set/get edge attributes)

- graph operations (copy, adjacency matrix, diagonal matrix, Laplacian matrix)

- major graph algorithms (exploration, connectivity, components, maximal
  component, Dijkstra, Floyd-Warshall, betweenness centrality)

- spectral graph theory (spectral radius, spectral gap, natural connectivity,
  algebraic connectivity, effective_resistance, and spanning tree count)

- a number of graph/network generation models (random graph, ER (Erdos Renyi),
  BA (Barabasi Albart), randomized BA, ring, tree, binary tree, BA tree,
  generalized BA, latent, lattice, Voronoi, DB (Degree Bounded), configuration
  model, random regular graph, Li-Miani graph)

- graph import/export in DOT (GraphViz) format

# HISTORY

The development of **graph-tools** started in 2007, which was initially an
extension to Graph module in CPAN (Comprehensive Perl Archvie Network) by
Jarkko Hietaniemi.  Our Perl module has been called **graphtools** for long
time and Perl module names were Graph::Util and Graph::Enhanced.
**graphtools** in Perl has been developed until 2018.  Python version of
**graph-tools** was born in 2018 by porting **graphtools** in Perl to Python.
Hence, the internal structure and the coding style receives significant
influence from Graph module by Jarkko Hietaniemi.

# EXAMPLE

```python
from graph_tools import Graph

# create a graph with four nodes and two edges
g = Graph(directed=True)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_vertex(4)
print(g)

# find the all shortest paths from vertex 1
dist, prev = g.dijkstra(1)
print(dist)

# generate BA graph with 100 vertices
g = Graph(directed=False).create_graph('barabasi', 100)

# check if all vertices are mutually connected
print(g.is_connected())

# compute the betweenness centrality of vertex 1
print(g.betweenness(1))
```

# INSTALLATION

```python
pip3 install graph-tools
```

# AVAILABILITY

The latest version of **graph-tools** module is available at PyPI
(https://pypi.org/project/graph-tools/) .

# SEE ALSO

graphviz - graph visualization software (https://graphviz.org/)

# AUTHOR

Hiroyuki Ohsaki <ohsaki[atmark]lsnl.jp>
