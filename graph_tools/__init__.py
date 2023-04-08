#!/usr/bin/env python3
#
# Tools for graph theory and network science with many generation models.
# Copyright (c) 2018-2023, Hiroyuki Ohsaki.
# All rights reserved.
#

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Contributors
# Ryo Nakamura <r-nakamura[atmark]kwansei.ac.jp>
# Yuichi Yasuda <yuichi[atmark]kwansei.ac.jp>

from collections import deque, defaultdict
import functools
import heapq
import itertools
import math
import random
import re
import statistics
import time

from perlcompat import warn, die
import numpy
import tbdump

VERSION = 1.10

CREATE_SUBP = {
    'random': 'create_random_graph',
    'random_sparse': 'create_random_sparse_graph',
    'erdos_renyi': 'create_erdos_renyi_graph',
    'er': 'create_erdos_renyi_graph',
    'barabasi': 'create_barabasi_graph',
    'ba': 'create_barabasi_graph',
    'barandom': 'create_barabasi_random_graph',
    'ring': 'create_ring_graph',
    'tree': 'create_tree_graph',
    'btree': 'create_btree_graph',
    'general_ba': 'create_generalized_barabasi_graph',
    'latent': 'create_latent_graph',
    'treeba': 'create_treeba_graph',
    'lattice': 'create_lattice_graph',
    'voronoi': 'create_voronoi_graph',
    'db': 'create_degree_bounded_graph',
    'degree_bounded': 'create_degree_bounded_graph',
    'configuration': 'create_configuration_graph',
    'regular': 'create_random_regular_graph',
    'li_maini': 'create_li_maini_graph',
    'preset': 'create_preset_graph',
    'star': 'create_star_graph',
}
CREATE_TYPES = sorted(CREATE_SUBP.keys())

IMPORT_SUBP = {
    'dot': 'import_dot',
    'edges': 'import_edge_list',
    'cell': 'import_cell',
}
IMPORT_FORMATS = sorted(IMPORT_SUBP.keys())

EXPORT_SUBP = {
    'dot': 'export_dot',
    'edges': 'export_edge_list',
    'cell': 'export_cell',
}
EXPORT_FORMATS = sorted(EXPORT_SUBP.keys())

MAX_RETRIES = 100

class Graph:
    # FIXME: Support non-multiedged graph; multiedged option is always
    # regarded as True.
    def __init__(self, directed=True, multiedged=True):
        self.G = {}  # Graph attributes.
        self.V = {}  # Vertices.
        self.EI = {}  # Incoming edges.
        self.EO = {}  # Outgoing edges.
        self._directed = directed
        self.invalidate_cache()

    def __repr__(self):
        return self.export_dot()

    def directed(self):
        return self._directed

    def undirected(self):
        return not self._directed

    def multiedged(self):
        return True

    def expect_directed(self):
        if not self.directed():
            die('directed graph expected')

    def expect_undirected(self):
        if not self.undirected():
            die('undirected graph expected')

    def expect_multiedged(self):
        if not self.multiedged():
            die('multiedged graph expected')

    # graph ----------------------------------------------------------------
    def set_graph_attribute(self, attr, val):
        """Define the graph attribute ATTR as value VAL."""
        self.G[attr] = val

    def get_graph_attribute(self, attr):
        """Extract and return the graph attribute named ATTR.  Return None if
        the attribute does not exist."""
        return self.G.get(attr, None)

    def degrees(self):
        """Return an (unsorted) list of vertex degrees."""
        return [self.degree(v) for v in self.vertices()]

    def average_degree(self):
        """Return the average degree of all vertices in the graph."""
        degrees = self.degrees()
        if degrees:
            return statistics.mean(degrees)
        else:
            return 0

    def ntriads(self, v):
        degree = self.degree(v)
        if degree <= 1:
            return 0
        # Count the number of triads among neighbors.
        ntriads = 0
        for u, w in itertools.combinations(self.neighbors(v), 2):
            if self.has_edge(u, w) or self.has_edge(w, u):
                ntriads += 1
        return ntriads

    def clustering_coefficient(self):
        ratios = []
        for v in self.vertices():
            degree = self.degree(v)
            if degree <= 1:
                ratio = 0
            else:
                ratio = self.ntriads(v) / ((degree * (degree - 1)) / 2)
            ratios.append(ratio)
        return statistics.mean(ratios)

    def average_path_length(self):
        lengths = []
        for u, v in itertools.permutations(self.vertices(), 2):
            lengths.append(self.shortest_path_length(u, v))
        return statistics.mean(lengths)

    def assortativity(self):
        """Calcurate the assortativity of a graph.  See Eq. (4) in
        M. E. J. Newman, `Assortative mixing in networks'."""
        M = self.nedges()
        prod_sum, add_sum, square_sum = 0, 0, 0
        for e in self.edges():
            j, k = [self.degree(v) for v in e]
            prod_sum += j * k
            add_sum += (j + k) / 2
            square_sum += (j**2 + k**2) / 2
        r = prod_sum / M - (add_sum / M)**2
        r /= square_sum / M - (add_sum / M)**2
        return r

    # vertex ----------------------------------------------------------------
    def vertices(self):
        """Return all vertices in the graph as a dictionary view object."""
        return self.V.keys()

    def vertex_index(self, v):
        # FIXME: This code is inefficient.
        for n, u in enumerate(self.vertices()):
            if u == v:
                return n
        return None

    def nvertices(self):
        return len(self.vertices())

    def has_vertex(self, v):
        """Check if vertex V exists in the graph."""
        return v in self.V

    def add_vertex(self, v):
        """Attach vertex V to the graph if it does not exist."""
        if not self.has_vertex(v):
            self.V[v] = {}  # Default vertex attribute.

    def add_vertices(self, *vertices):
        """Attach vertices VERTICES to the graph while avoiding duplicates."""
        for v in vertices:
            self.add_vertex(v)

    def predecessors(self, v, ignore=False):
        """Return the list of predecessor vertices of vetex V.  This method
        fails if the graph is undirected.  Error checking is bypassed if
        IGNORE is true."""
        if not ignore:
            self.expect_directed()
        if v not in self.EI:
            return []
        return self.EI[v].keys()

    def successors(self, u, ignore=False):
        """Return the list of successor vertices of vetex V.  This method
        fails if the graph is undirected.  Error checking is bypassed if
        IGNORE is true."""
        if not ignore:
            self.expect_directed()
        if u not in self.EO:
            return []
        return self.EO[u].keys()

    def neighbors(self, v, k=1):
        """Return all K-hop neighbor nodes of vetex V."""
        if k <= 1:
            in_vertices = set(self.predecessors(v, ignore=True))
            out_vertices = set(self.successors(v, ignore=True))
            return in_vertices.union(out_vertices)
        else:
            vertices = set(self.neighbors(v))
            for u in self.neighbors(v):
                found = self.neighbors(u, k - 1)
                vertices.union(found)
            return vertices

    def set_vertex_attribute(self, v, attr, val):
        """Define the vertex attribute of vetex V named ATTR as value VAL."""
        if not self.has_vertex(v):
            die(f'set_vertex_attribute: no vertex {v}')
        self.V[v][attr] = val

    def get_vertex_attribute(self, v, attr):
        """Return the vertex attribute of vetex V named ATTR."""
        return self.V[v].get(attr, None)

    def set_vertex_attributes(self, v, adict):
        """Set multiple vertex attributes of vertex V.  Attributes are passed
        by dictionary ADICT."""
        for key, val in adict.items():
            self.set_vertex_attribute(v, key, val)

    def get_vertex_attributes(self, v):
        """Return all vertex attributes of vertex V as dictionary."""
        return self.V.get(v, {})

    def set_vertex_weight(self, v, val):
        """Set the weight of vertex V to VAL.  The vertex weight is stored in
        vertex attributes with the name 'weight'."""
        return self.set_vertex_attribute(v, 'weight', val)

    def get_vertex_weight(self, v):
        """Return the vertex weight of vertex V."""
        return self.get_vertex_attribute(v, 'weight')

    def delete_vertex(self, v):
        """Delete vertex V from the graph.  All incoming/outgoing edges are
        deleted."""
        if not self.has_vertex(v):
            return
        for u in self.neighbors(v):
            try:
                del self.EO[u][v]
            except KeyError:
                pass
            try:
                del self.EI[u][v]
            except KeyError:
                pass
        del self.V[v]
        try:
            del self.EO[v]
        except KeyError:
            pass
        try:
            del self.EI[v]
        except KeyError:
            pass
        self.invalidate_cache()

    def delete_vertices(self, alist):
        """Delete all vertices ALST."""
        for v in alist:
            self.delete_vertex(v)

    def delete_all_vertices(self):
        self.delete_vertices(self.vertices())

    def random_vertex(self):
        """Randomly choose a vertex from all vertices."""
        return random.choice(list(self.vertices()))

    # edge ----------------------------------------------------------------
    def edges(self):
        """Return all edges in the graph as a list."""
        for u in self.EO:
            for v in self.EO[u]:
                for id_ in self.get_multiedge_ids(u, v):
                    yield u, v

    def nedges(self):
        return len(list(self.edges()))

    def edge_indices(self):
        for u, v in self.edges():
            yield self.vertex_index(u), self.vertex_index(v)

    def unique_edges(self):
        """Return all unique edges in the graph as a list.  All multi-edges
        are unified into a single edge."""
        found = []
        for u in self.EO:
            for v in self.EO[u]:
                found.append([u, v])
        return found

    def has_edge(self, u, v):
        """Check if the graph has edge (u, v)."""
        if self.undirected() and u > v:
            u, v = v, u
        if u not in self.EO:
            return False
        return v in self.EO[u]

    def get_multiedge_ids(self, u, v):
        """Return the edge identifiers (starting from zero) of edges between
        vertex U and vertex V as a list.  For instance, two vertices connected
        by a single edge yields [0].  Note that the order of edge identifers
        are random."""
        if self.undirected() and u > v:
            u, v = v, u
        if not self.has_edge(u, v):
            return []
        return self.EO[u][v].keys()

    def get_edge_count(self, u, v):
        """Return the number of multi-edges connecting vertices U and V."""
        ids = self.get_multiedge_ids(u, v)
        if ids:
            return len(ids)
        else:
            return 0

    def add_edge(self, u, v, weight=None):
        """Add an edge from vertex U to vertex V."""
        if self.undirected() and u > v:
            u, v = v, u
        count = self.get_edge_count(u, v)
        self.add_vertices(u, v)
        if u not in self.EO:
            self.EO[u] = {}
        if v not in self.EO[u]:
            self.EO[u][v] = {}
        self.EO[u][v][count] = {}  # default edge attributes
        if v not in self.EI:
            self.EI[v] = {}
        if u not in self.EI[v]:
            self.EI[v][u] = {}
        self.EI[v][u][count] = {}  # default edge attributes
        if weight:
            self.set_edge_weight_by_id(u, v, count, weight)
        self.invalidate_cache()
        return count

    def delete_edge(self, u, v):
        """Delete an edge between vertices U and V.  If vertices are connected
        with multiple edges, the one with the largest edge identifier is
        deleted."""
        if self.undirected() and u > v:
            u, v = v, u
        if not self.has_edge(u, v):
            return
        count = self.get_edge_count(u, v) - 1
        del self.EO[u][v][count]
        del self.EI[v][u][count]
        if count == 0:
            del self.EO[u][v]
            del self.EI[v][u]
        self.invalidate_cache()
        return count

    def edges_from(self, u, ignore=False):
        """Return the list of edges leaving from vertex U.  This method fails
        if the graph is undirected.  Error checking is bypassed if
        IGNORE is true."""
        if not ignore:
            self.expect_directed()
        found = []
        for v in self.successors(u, ignore):
            for id_ in self.get_multiedge_ids(u, v):
                found.append([u, v])
        return found

    def edges_to(self, v, ignore=False):
        """Return the list of edges coming to vertex U.  This method fails
        if the graph is undirected.  Error checking is bypassed if
        IGNORE is true."""
        if not ignore:
            self.expect_directed()
        found = []
        for u in self.predecessors(v, ignore):
            for id_ in self.EI[v][u]:
                found.append([u, v])
        return found

    def edges_at(self, v):
        """Return the list of all edges connected to vertex V."""
        found = []
        found.extend(self.edges_from(v, ignore=True))
        found.extend(self.edges_to(v, ignore=True))
        return found

    def out_degree(self, u, ignore=False):
        """Return the number outgoing edges connected to vertex U.  This
        method fails if the graph is undirected.  Error checking is bypassed
        if IGNORE is true."""
        if not ignore:
            self.expect_directed()
        return len(self.edges_from(u))

    def in_degree(self, v, ignore=False):
        """Return the number incoming edges connected to vertex V.  This
        method fails if the graph is undirected.  Error checking is bypassed
        if IGNORE is true."""
        if not ignore:
            self.expect_directed()
        return len(self.edges_to(v))

    def degree(self, v):
        """Return the number of edges connected to vetex V."""
        if self.undirected():
            return len(self.edges_at(v))
        else:
            return self.in_degree(v) + self.out_degree(v)

    vertex_degree = degree

    def random_edge(self):
        """Randomly choose an edge from all edges."""
        return random.choice(list(self.edges()))

    def set_edge_attribute_by_id(self, u, v, n, attr, val):
        """Define the attribute of the N-th edge between vertices U and V
        named ATTR as value VAL."""
        if not attr:
            die('set_edge_attribute_by_id: no attribute specified.')
        if self.undirected() and u > v:
            u, v = v, u
        if not self.has_edge(u, v):
            die(f'set_edge_attribute_by_id: edge ({u}, {v}) not found')
        self.EO[u][v].setdefault(n, {})
        self.EO[u][v][n][attr] = val

    def get_edge_attribute_by_id(self, u, v, n, attr):
        """Return the attribute of the N-th edge between vertices U and V
        named ATTR."""
        if not attr:
            die('get_edge_attribute_by_id: no attribute specified.')
        if self.undirected() and u > v:
            u, v = v, u
        if not self.has_edge(u, v):
            die(f'get_edge_attribute_by_id: edge ({u}, {v}) not found')
        return self.EO[u][v][n].get(attr, None)

    def set_edge_attributes_by_id(self, u, v, n, adict):
        """Define attributes of the N-th edge between vertices U and V from
        dictionary ADICT."""
        for key, val in adict.items():
            self.set_edge_attribute_by_id(u, v, n, key, val)

    def get_edge_attributes_by_id(self, u, v, n):
        """Return all attributes of the N-th edge between vertices U and V as
        dictionary."""
        if self.undirected() and u > v:
            u, v = v, u
        return self.EO[u][v].get(n, {})

    def set_edge_weight_by_id(self, u, v, n, val):
        """Set the edge weight of the N-th edge between vertices U and V to
        value VAL."""
        return self.set_edge_attribute_by_id(u, v, n, 'weight', val)

    def get_edge_weight_by_id(self, u, v, n):
        """Return the edge weight of the N-th edge between vertices U and
        V."""
        return self.get_edge_attribute_by_id(u, v, n, 'weight')

    def set_edge_weight(self, u, v, w):
        """Set the weight of the first edge between vertices U and V to weight
        W."""
        return self.set_edge_weight_by_id(u, v, 0, w)

    def get_edge_weight(self, u, v):
        """Return the weight of the first edge between vertices U and V to
        weight
        W."""
        return self.get_edge_weight_by_id(u, v, 0)

    # algorithm ----------------------------------------------------------------
    def invalidate_cache(self):
        self.T = {}  # Shortest path cache (total distances from vertex).
        self.P = {}  # Shortest path cache (preceeding vertices list).
        self.cache = {'between': {}, 'close': {}, 'eigen': {}, 'eccent': {}}

    # https://stackoverflow.com/questions/22897209/dijkstras-algorithm-in-python
    def _dijkstra(self, s):
        """Compute all shortest paths from source vertex S using Dijkstra's
        algorithm.  Return two dictionaries: DIST and PREV.  Dictionary DIST
        records the distance to every vertex in the graph, and dictionary PREV
        records the *previous* node along the shortest-path tree from source
        vertex S.  For instance, DIST['4'] indicates the number of hops from
        source vertex S to vertex '4'.  PREV['4'] indicates the preceeding
        vertex in the shortest path from source vertex S to vertex 4.  You can
        obtain the shortest path to vertex V by traversing dictionary PREV
        from vertex V back to source node S."""
        # Distance from the source to other vertices.
        dist = {v: math.inf for v in self.vertices()}
        dist[s] = 0
        prev = {v: [] for v in self.vertices()}

        visited = {}
        queue = []
        heapq.heappush(queue, (dist[s], s))
        while queue:
            _, u = heapq.heappop(queue)
            visited[u] = True
            neighbors = self.neighbors(
                u) if self.is_undirected() else self.successors(u)
            for v in neighbors:
                if v in visited:
                    continue
                # FIXME: Must reject multi-edged graph.
                w = self.get_edge_weight_by_id(u, v, 0) or 1
                if dist[u] + w < dist[v]:  # Is vertex V closer from vertex U?
                    dist[v] = dist[u] + w
                    prev[v] = [u]
                elif dist[u] + w == dist[v]:  # Handle equal paths.
                    prev[v].append(u)
                heapq.heappush(queue, (dist[u] + w, v))

        self.T[s] = dist
        self.P[s] = prev
        return dist, prev

    def dijkstra(self, s):
        if s in self.T and s in self.P:
            return self.T[s], self.P[s]
        else:
            return self._dijkstra(s)

    def shortest_paths(self, s, t):
        """Return the all shortest-paths from vertex S to vertex T.  Note that
        the shortest path tree from source vertex S is cached for efficiency.
        So, if the network topology is changed since the previous invocation
        of shortest_paths(), you must explicitly call dijkstra() to renew the
        shorte-path tree cache."""
        def find_path(s, t):
            # P[s] stores the shortest-path tree from vertex S.
            # P{s][t] is a set of previous nodes in the shortest-path tree.
            for prev in self.P[s][t]:
                if prev == s:
                    yield [s, t]
                else:
                    for path in find_path(s, prev):
                        yield path + [t]

        # Build shortest-path tree if not cached yet.
        if not s in self.T:
            self.dijkstra(s)
        return find_path(s, t)

    def shortest_path_length(self, s, t):
        # Build shortest-path tree if not cached yet.
        if s not in self.T:
            self.dijkstra(s)
        if t in self.T[s]:
            return self.T[s][t]
        else:
            return math.inf

    def dijkstra_all_pairs(self):
        """Compute all-pairs shortest paths using Dijkstra's algorithm."""
        for v in self.vertices():
            self.dijkstra(v)

    def floyd_warshall(self):
        """Compute all-pairs shortest paths using Floyd-Warshall algorithm."""
        # Initialize weight matrix.
        dist = {}
        next_ = {}
        for v in self.vertices():
            dist[v] = defaultdict(lambda: math.inf)
            next_[v] = {}

        for u, v in self.edges():
            w = self.get_edge_weight_by_id(u, v, 0) or 1
            dist[u][v] = w
            if self.is_undirected():
                dist[v][u] = w

        # Run Floyd-Warshall algorithm to find all-pairs shortest paths.
        for k in self.vertices():
            for u in self.vertices():
                for v in self.vertices():
                    if dist[u][k] + dist[k][v] < dist[u][v]:
                        dist[u][v] = dist[u][k] + dist[k][v]
                        next_[u][v] = k
        self.T = dist

    def is_reachable(self, u, v):
        """Check if any path exists from vertex U to vertex V."""
        return self.shortest_path_length(u, v) < math.inf

    def explore(self, s):
        """Return the list of all vertices reachable from vertex S."""
        explored = set()
        need_visit = set()
        need_visit.add(s)
        while need_visit:
            u = need_visit.pop()
            explored.add(u)
            for v in self.neighbors(u):
                if v not in explored:
                    need_visit.add(v)
        return explored

    def is_connected(self):
        """Check if all vertices in the graph are mutually connected."""
        v = self.random_vertex()
        explored = self.explore(v)
        return len(explored) == self.nvertices()

    def components(self):
        """Return all components (i.e., connected subgraphs) of the graph.
        Components are returned as list of vertices set."""
        components = []
        # Record unvisisted vertices.
        unvisited = set(self.vertices())
        while unvisited:
            # Start exploration from one of unvisited vertices.
            v = unvisited.pop()
            explored = self.explore(v)
            components.append(explored)
            # Remove all visisted vertices.
            unvisited -= set(explored)
        return components

    def maximal_component(self):
        """Return the largest component (i.e., the component with the largest
        number of vertices)."""
        maximal = max(self.components(), key=lambda x: len(x))
        return maximal

    def degree_centrality(self, v):
        return self.degree(v) / (self.nvertices() - 1)

    def betweenness_centralities(self):
        """Return the betweenness centrality for vertex v.  This program
        implements Algorithm 1 (betweenness centrality in unweighted graphs)
        in U. Brandes, `A Fast Algorithm for Betweeness Centrality,' Journal
        of Mathematical Sociology, 2001."""
        self.expect_undirected()
        # Check if the graph is unweighted.
        for _ in range(10):  # test 10 sample edges
            u, v = self.random_edge()
            w = self.get_edge_weight(u, v)
            if w is not None and w != 1:
                die('only supports unweighted graphs.')

        # Betweenness centrality for vertices.
        betweenness = {v: 0 for v in self.vertices()}

        for s in self.vertices():
            S = []  # Empty stack.
            P = {w: [] for w in self.vertices()}
            sigma = {t: 0 for t in self.vertices()}
            sigma[s] = 1
            d = {t: -1 for t in self.vertices()}
            d[s] = 1
            Q = deque()  # Empty queue.
            Q.append(s)

            while Q:
                v = Q.popleft()
                S.append(v)
                for w in self.neighbors(v):
                    # Found for the first time?
                    if d[w] < 0:
                        Q.append(w)
                        d[w] = d[v] + 1
                    # Shortest path to w via v?
                    if d[w] == d[v] + 1:
                        sigma[w] += sigma[v]
                        P[w].append(v)

            delta = {v: 0 for v in self.vertices()}
            # S returns vertices in order of non-increasing distance from s.
            while S:
                w = S.pop()
                for v in P[w]:
                    delta[v] += sigma[v] / sigma[w] * (1 + delta[w])
                if w != s:
                    betweenness[w] += delta[w]
        return betweenness

    def betweenness_centrality(self, v, normalize=True):
        if not v in self.cache['between']:
            self.cache['between'] = self.betweenness_centralities()
        if normalize:
            n = self.nvertices()
            return self.cache['between'][v] / ((n - 1) * (n - 2))
        else:
            return self.cache['between'][v]

    def closeness_centrality(self, v):
        if not v in self.cache['close']:
            d = sum([self.shortest_path_length(v, u) for u in self.vertices()])
            self.cache['close'][v] = self.nvertices() / d
        return self.cache['close'][v]

    def eigenvector_centrality(self, v):
        if not v in self.cache['eigen']:
            adj = self.adjacency_matrix()
            eigvec = self.maximum_eigvec(adj)
            eigvec = self.negate_if_negative(eigvec)
            vi = self.vertex_index(v)
            self.cache['eigen'][v] = float(eigvec[vi])
        return self.cache['eigen'][v]

    def eccentricity(self, v):
        if not v in self.cache['eccent']:
            lengths = [
                self.shortest_path_length(v, u) for u in self.vertices()
            ]
            self.cache['eccent'][v] = max(lengths)
        return self.cache['eccent'][v]

    def eccentricities(self):
        for v in self.vertices():
            yield self.eccentricity(v)

    # graph ----------------------------------------------------------------
    def copy_graph(self, directed=None):
        """Return a copy of the graph."""
        if directed is None:
            directed = self.directed()
        T = Graph(directed)
        # FIXME: Preserve graph attributes.
        for v in self.vertices():
            T.add_vertex(v)
            T.set_vertex_attributes(v, self.get_vertex_attributes(v))
        directed_from_undirected = directed and not self.directed()
        for u, v in self.edges():
            T.add_edge(u, v)
            if directed_from_undirected:
                T.add_edge(v, u)
            for n in self.get_multiedge_ids(u, v):
                T.set_edge_attributes_by_id(
                    u, v, n, self.get_edge_attributes_by_id(u, v, n))
                if directed_from_undirected:
                    T.set_edge_attributes_by_id(
                        v, u, n, self.get_edge_attributes_by_id(u, v, n))
        return T

    def directed_copy(self):
        """Return a directed copy of the graph."""
        return self.copy_graph(directed=True)

    def complete_graph(self):
        """Add edges to all vertex pairs to make the graph fully-meshed."""
        for u in self.vertices():
            for v in self.vertices():
                # FIXME: Works for directed/undirected graphs?
                if u >= v:
                    continue
                if not self.has_edge(u, v):
                    self.add_edge(u, v)
        return self

    # spectral measures ----------------------------------------------------------------
    def adjacency_matrix(self):
        """Return the adjacency matrix of the graph as NumPy.ndarray
        object."""
        N = self.nvertices()
        m = numpy.zeros((N, N), int)
        for u, v in self.edges():
            w = self.get_edge_weight(u, v) or 1
            ui = self.vertex_index(u)
            vi = self.vertex_index(v)
            m[ui, vi] += w
            if self.undirected():
                m[vi, ui] += w
        return m

    def diagonal_matrix(self):
        """Return the diagonal matrix of the graph as NumPy.ndarray object."""
        N = self.nvertices()
        m = numpy.zeros((N, N), int)
        for v in self.vertices():
            vi = self.vertex_index(v)
            m[vi, vi] = self.degree(v)
        return m

    def laplacian_matrix(self):
        """Return the Laplacian matrix of the graph as NumPy.ndarray
        object."""
        return self.diagonal_matrix() - self.adjacency_matrix()

    def weight_matrix(self):
        """Return the weight matrix of the graph."""
        nvertices = self.nvertices()
        m = numpy.zeros((nvertices, nvertices), int)
        for u, v in self.edges():
            w = self.get_edge_weight(u, v) or 1
            ui = self.vertex_index(u)
            vi = self.vertex_index(v)
            m[ui][vi] = m[vi][ui] = w
        return m

    def eigenvals(self, m):
        """Return the eigenvalues of matrix M.  Eigenvalues are sorted in the
        ascending order."""
        eigvals, eigvecs = numpy.linalg.eig(m)
        return sorted(eigvals)

    def maximum_eigenval(self, m):
        """Return the maximum eigenvalue of matrix M."""
        return self.eigenvals(m)[-1]

    def maximum_eig(self, m):
        """Return the maximum eigenvalue and the corresponding eigenvector of
        matrix M."""
        eigvals, eigvecs = numpy.linalg.eig(m)
        # Find the index of the largest eigenvalue.
        i = max(enumerate(eigvals), key=lambda x: x[1])[0]
        eigval = eigvals[i]
        eigvec = eigvecs[:, i]
        return eigval, eigvec

    def maximum_eigvec(self, m):
        """Return the principal eigenvector of M corresponding the maximum
        eigenvalue."""
        eigvals, eigvec = self.maximum_eig(m)
        return eigvec

    def negate_if_negative(self, vec):
        if any([x < 0 for x in vec]):
            vec = -vec
        return vec

    def adjacency_matrix_eigvals(self):
        """Return eigenvalues of the adjacency matrix."""
        return self.eigenvals(self.adjacency_matrix())

    def laplacian_matrix_eigvals(self):
        """Return eigenvalues of the Laplacian matrix."""
        return self.eigenvals(self.laplacian_matrix())

    def spectral_radius(self):
        """Return the spectral raduis from spectral graph theory."""
        return self.maximum_eigenval(self.adjacency_matrix())

    def spectral_gap(self):
        """Return the spectral gap from spectral graph theory."""
        lmbda = self.adjacency_matrix_eigvals()
        return lmbda[-1] - lmbda[-2]

    def natural_connectivity(self):
        """Return the natural connectivity from spectral graph theory."""
        N = self.nvertices()
        lmbda = self.adjacency_matrix_eigvals()
        return math.log(sum(numpy.exp(lmbda)) / N)

    def algebraic_connectivity(self):
        """Return the argebraic connectivity from spectral graph theory."""
        mu = self.laplacian_matrix_eigvals()
        return mu[1]

    def effective_resistance(self):
        """Return the effective resistance from spectral graph theory."""
        N = self.nvertices()
        mu = self.laplacian_matrix_eigvals()
        return N * sum([1 / (mu + 1e-100) for mu in mu[1:]])

    def spanning_tree_count(self):
        """Return the spanning tree count from spectral graph theory."""
        N = self.nvertices()
        mu = self.laplacian_matrix_eigvals()
        return functools.reduce(lambda x, y: x * y,
                                [1 / (mu + 1e-100) for mu in mu[1:]]) / N

    # coarsening ----------------------------------------------------------------
    def merge_vertices(self, u, v):
        """Delete vertex V after connecting all neighbors (except vertex U) of
        vertex V to vertex U."""
        if u == v:
            return
        for t in self.neighbors(v):
            if t == u:
                continue
            # Note: w might be None (i.e., no weight specified)
            w = self.get_edge_weight(v, t)
            if self.has_edge(u, t):
                new_w = self.get_edge_weight(u, t) or 1
                new_w += w or 1
                self.set_edge_weight(u, t, new_w)
            else:
                self.add_edge(u, t)
                self.set_edge_weight(u, t, w)
            self.delete_edge(v, t)
        self.delete_edge(u, v)
        self.delete_vertex(v)
        super_ = self.get_vertex_attribute(u, 'super') or ''
        super_ += f'+{v}'
        self.set_vertex_attribute(u, 'super', super_)

    def RM_coarsening(self, alpha=.5):
        # FIXME: Support directed graphs.
        self.expect_undirected()
        n = self.nvertices()
        nretries = 0
        while self.nvertices() > n * alpha and nretries < 100:
            u, v = self.random_edge()
            u_merged = self.get_vertex_attribute(u, 'super')
            v_merged = self.get_vertex_attribute(v, 'super')
            if u_merged or v_merged:
                nretries += 1
                continue
            self.merge_vertices(u, v)

    def COARSENET_coarsening(self, alpha=.5):
        def d_lambda(lmbda, u, v, a, b):
            # Weight of edges (a, b) and (b, a), respectively.
            # FIXME: This code assumes vertex IDs are 1...N.
            ai = self.vertex_index(a)
            bi = self.vertex_index(a)
            beta1 = adj[ai][bi]
            beta2 = adj[bi][ai]
            # Calcurate \Delta \lambda using Proposition 5.2.
            term1 = (u[ai] * v[ai] + u[bi] * v[bi])
            ut_co = (1 + beta2) / 2 * (lmbda * u[ai] - beta1 * u[bi])
            ut_co += (1 + beta1) / 2 * (lmbda * u[bi] - beta2 * u[ai])
            term2 = v[ai] * ut_co
            term3 = beta2 * u[ai] * v[bi]
            term4 = beta1 * u[bi] * v[ai]
            denom = numpy.dot(v.transpose(), u)
            delta = (-lmbda * term1 + term2 + term3 + term4) / (denom - term1)
            return float(delta)

        # FIXME: Support directed graphs.
        self.expect_undirected()
        n = self.nvertices()
        adj = self.adjacency_matrix()
        # lmbda: lambda (maximum eigenvalue of adjacency matrix)
        # u: right eigenvector corresponding lambda
        # v: left eigenvector corresponding lambda
        lmbda, u = self.maximum_eig(adj)
        # NOTE: u = v if adjacency matrix is symmetric?
        _, v = self.maximum_eig(adj.transpose())
        # u and v must be positive unit vector.
        u = self.negate_if_negative(u)
        v = self.negate_if_negative(v)
        scores = [(d_lambda(lmbda, u, v, a, b), a, b) for a, b in self.edges()]
        scores = sorted(scores, key=lambda x: abs(x[0]))
        # for s, u, v in scores:
        #     warn(f'coarsenet: score={s:6.3f} edge=({u}, {v})')

        while self.nvertices() > n * alpha and len(scores) > 0:
            _, u, v = scores.pop(0)
            # Vertex might have been merged with another vertex.
            if not self.has_vertex(u) or not self.has_vertex(v):
                continue
            self.merge_vertices(u, v)

    def MGC_coarsening(self, alpha=.5):
        # FIXME: Seems to fail appropriate coarsening.
        def find_vertices_with_minimum_distance():
            adj = self.adjacency_matrix()
            dmin = None
            for v in self.vertices():
                # Collect all 2-hop neighbors from vertex V.
                neighbors = self.neighbors(v, 2)
                dv = self.degree(v)
                iv = self.vertex_index(v)
                lv = adj[iv]
                for u in neighbors:
                    du = self.degree(u)
                    iu = self.vertex_index(u)
                    lu = adj[iu]
                    d = numpy.linalg.norm(lv / dv - lu / du, ord=1)
                    if dmin is None or d < dmin:
                        dmin = d
                        umin, vmin = u, v
            return dmin, umin, vmin

        # FIXME: Support directed graphs.
        self.expect_undirected()
        nvertices = self.nvertices()
        while self.nvertices() > nvertices * alpha:
            d, u, v = find_vertices_with_minimum_distance()
            # warn(f'MGC: d={d} edge=({u}, {v})')
            self.merge_vertices(u, v)

    def add_edges_from_sparse_graph(self, sg):
        edges = sg.get_edge_list()
        for ui, vi, w in zip(*edges):
            u, v = ui + 1, vi + 1
            self.add_edge(u, v)
            self.set_edge_weight(u, v, w)

    def local_variation_coarsening(self,
                                   alpha=.5,
                                   method='variation_neighborhoods'):
        import pygsp
        import graph_coarsening
        # Create a sparse graph from the weight matrix.
        g = pygsp.graphs.Graph(self.weight_matrix())
        r = 1 - alpha  # r=0 means no reduction.
        # C: coarsening matrix
        # Gc: coarsed graph
        # Call: coarsening matrix at all levels
        # Gall: coarsed graph at all levels
        C, Gc, Call, Gall = graph_coarsening.coarsen(g, r=r, method=method)
        # Rewrite all edges according to Gc.
        self.delete_vertices(list(self.vertices()))
        self.add_edges_from_sparse_graph(Gc)

    def LVN_coarsening(self, alpha=.5):
        self.local_variation_coarsening(alpha=alpha,
                                        method='variation_neighborhood')

    def LVE_coarsening(self, alpha=.5):
        self.local_variation_coarsening(alpha=alpha, method='variation_edges')

    def kron_coarsening(self, alpha=.5):
        self.local_variation_coarsening(alpha=alpha, method='kron')

    def HEM_coarsening(self, alpha=.5):
        self.local_variation_coarsening(alpha=alpha, method='heavy_edge')

    def affinity_coarsening(self, alpha=.5):
        self.local_variation_coarsening(alpha=alpha, method='affinity_JC')

    def minimum_degree_coarsening(self, degree=2):
        to_delete = set()
        for v in self.vertices():
            if self.degree(v) <= degree:
                # Mark neighbor vertices as supernode.
                for u in self.neighbors(v):
                    if u not in to_delete:
                        super_ = self.get_vertex_attribute(u, 'super') or ''
                        super_ += f'+{v}'
                        self.set_vertex_attribute(u, 'super', super_)
                to_delete.add(v)

        self.delete_vertices(to_delete)

    # util ----------------------------------------------------------------
    def header_string(self, comment='# '):
        date = time.strftime('%Y/%M/%D %H:%M:%S', time.localtime())
        atype = 'directed' if self.is_directed() else 'undirected'
        vcount = self.nvertices()
        ecount = self.nedges()
        astr = f"{comment}Generated by graph-tools (version {VERSION}) at {date}]\n{comment}{atype}, {vcount} vertices, {ecount} edges\n"
        return astr

    # create ----------------------------------------------------------------
    def create_graph(self, atype, *args, **kwargs):
        name = CREATE_SUBP.get(atype, None)
        if not name:
            warn(f"create_graph: no graph creation support for type '{atype}'")
            return None
        func = getattr(self, name, None)
        if not func:
            warn(f"create_graph: graph creation method '{name}' not found")
            return None
        return func(*args, **kwargs)

    def create_random_graph(self, N=10, E=20, multiedged=False):
        """Create an instance of *connected* random graphs with N vertices and
        E edges.  A generated graph is non-multiedged.  You can specify
        MULTIEDGED to allow a multiedged graph."""
        if E < N:
            die('create_random_graph: too small number of edges')
        self.add_vertices(*range(1, N + 1))

        # Add first (N - 1) edges for making sure connectivity.
        for i in range(1, N):
            u = i + 1
            v = random.randrange(1, u)
            if random.random() <= .5:
                self.add_edge(u, v)
            else:
                self.add_edge(v, u)

        # Randomly add remaining (E - (N - 1)) edges.
        for i in range(1, E - (N - 1) + 1):
            # FIXME: Avoid cycle edges, but this may take log time.
            ntries = 1
            while ntries < MAX_RETRIES:
                u = self.random_vertex()
                v = self.random_vertex()
                if u == v:
                    continue
                if multiedged:
                    break
                else:
                    if not self.has_edge(u, v):
                        break
            self.add_edge(u, v)
        return self

    def create_random_sparse_graph(self, N=10, E=20, no_multiedge=False):
        """Create a random graph with N vertices and E edges.  A generated
        graph might be multiedged."""
        for i in range(1, N + 1):
            self.add_vertex(i)

        # Randomly add remaining Eedges.
        for i in range(1, E + 1):
            # FIXME: Avoid cycle edges, but this may take log time.
            ntries = 1
            while ntries < MAX_RETRIES:
                u = random.randrange(1, N + 1)
                v = random.randrange(1, N + 1)
                if not no_multiedge and u != v:
                    break
                if no_multiedge and u != v and not self.has_edge(u, v):
                    break
            self.add_edge(u, v)
        return self

    def create_erdos_renyi_graph(self, N=10, p=.5):
        """Create a random graph using the Erdos-Renyi model."""
        self.add_vertices(*range(1, N + 1))
        for u, v in itertools.combinations(self.vertices(), 2):
            if random.random() < p:
                self.add_edge(u, v)
        return self

    def create_barabasi_graph(self, N=10, m0=2, m=2):
        """Create a scale-free graph using the BA (Barabasi Albert) model."""
        self.expect_undirected()
        # Create complete graph with m0 vertices.
        for v in range(1, m0 + 1):
            self.add_vertex(v)
        self = self.complete_graph()

        step = N - m0
        for k in range(1, step + 1):
            # Add a new vertex with m edges.
            u = m0 + k
            self.add_vertex(u)

            # Attach to a vertex using preferential attachment.
            edges = list(self.edges())
            for i in range(1, m + 1):
                # NOTE: Degree-preferential attachment is realized by
                # selecting a vertex connected to a randomly-chosen edge.
                edge = random.choice(edges)
                v = random.choice(edge)
                self.add_edge(u, v)
        return self

    def create_barabasi_random_graph(self, N=10, E=20, m0=2):
        """Create a scal-free graph using the randomized BA model.  The
        randomized BA model is a generalization of the original BA model; The
        number edges added at every preferential-attachment phase is given by
        a random number rather than a constant parameter m.  This enables
        generation of BA-model-like scale-free graphs with an arbitrary
        average degree."""
        self.expect_undirected()
        # Create complete graph with m0 vertices.
        for v in range(1, m0 + 1):
            self.add_vertex(v)
        self = self.complete_graph()

        # Calcurate number of edges to be connected per vertex.
        E0 = m0 * (m0 - 1) / 2
        nedges = (E - E0) / (N - m0)

        # Add remaining (N - m0) vertices.
        for u in range(m0 + 1, N + 1):
            self.add_vertex(u)

            # Attach to a vertex using preferential attachment.
            # NOTE: Degree-preferential attachment is realized by
            # selecting a vertex connected to a randomly-chosen edge.
            edges = list(self.edges())
            while True:
                edge = random.choice(edges)
                v = random.choice(edge)
                self.add_edge(u, v)

                # NOTE: Using the fact that the average number of
                # successes of infinite Bernoulli traials with probability
                # p is given by 1/p.
                if random.uniform(0, 1) <= 1 / nedges:
                    break
        return self

    def create_ring_graph(self, N=10, step=1):
        """Create a ring graph with N vertices.  Every edge is connected
        between vertices whose identifiers are STEP-th apart; e.g., when STEP
        is k, vertex i is connected with vertex (i - k) and (i + k)."""
        for v in range(1, N + 1):
            self.add_vertex(v)
        # Add (N - 1) edges for making circular topology.
        for _ in range(0, N):
            u = _ + 1
            v = ((_ + step) % N) + 1
            self.add_edge(u, v)
        return self

    def create_tree_graph(self, N=10, *kargs, **kwargs):
        """Create a random tree graph with N vertices."""
        self.add_vertex(1)
        # Add (N - 1) edges for _ in making tree topology.
        for v in range(2, N + 1):
            u = self.random_vertex()
            self.add_edge(u, v)
        return self

    # Binary tree graph.
    def create_btree_graph(self, N=10):
        """Create a balanced binay tree grph with N vertices."""
        depth = 0
        nedges = 1
        finished = False
        while not finished:
            vleft = 2**depth
            for count in range(1, 2**depth + 1):
                v = vleft + (count - 1)
                parent = int(v / 2)
                if (parent == 0):
                    continue
                self.add_edge(v, parent)
                # Record the depth of the vertex from the root.
                # FIXME: This code can be removed?
                self.set_vertex_attribute(v, 'latent', 1 / depth)
                nedges += 1
                if nedges >= N:
                    finished = True
                    break
            depth += 1
        return self

    # Tree BA graph.
    def create_treeba_graph(self, N=10, alpha=1):
        # self.expect_directed()
        attract = [0] * N
        # Create an initial vertex.
        v = 1
        self.add_vertex(v)
        attract[v - 1] = alpha + len(self.in_edges(v))

        # Create a vertex and attach to another using preferential attachment.
        for u in range(2, N + 1):
            # Randomly choose a vertex with a probability proportional to attract.
            v = random.choices(range(1, N + 1), weights=attract, k=1)[0]
            self.add_edge(u, v)
            attract[u - 1] = alpha + len(self.in_edges(u))
            attract[v - 1] = alpha + len(self.in_edges(v))
        return self

    def create_generalized_barabasi_graph(self, N=10, m=2, gamma=3, m0=2):
        """Create a scale-free graph using the generalized BA model proposed
        in S. N. Dorogovtsev, ``Structure of growing networks with
        preferential linking,'' Phisical Review Letters, vol. 85, no. 21,
        pp. 4633-14636, Nov. 2000."""
        A = m * (gamma - 2)
        # Create complete graph with m0 vertices.
        self.add_vertices(*range(1, m0 + 1))
        self = self.complete_graph()

        step = N - m0
        for i in range(1, step + 1):
            # Add a new vertex with m edges.
            u = m0 + i
            self.add_vertex(u)

            # Attach to a vertex using preferential attachment.
            vcount = self.nvertices() - 1
            ecount = self.nedges()
            for j in range(1, m + 1):
                # NOTE: Preferential-attachement with probability A + in_degree.
                total = A * vcount + ecount
                chosen = random.uniform(0, total)
                accum = 0
                for v in range(1, u):
                    accum += A + self.in_degree(v)
                    if chosen < accum:
                        # Make sure newly added node has at least single link.
                        if j == 1:
                            self.add_edge(u, v)
                        else:
                            self.add_edge(random.randrange(1, u + 1), v)
                        break
        return self

    def create_latent_graph(self,
                            N=10,
                            E=20,
                            error_ratio=0,
                            confer='linear',
                            dist='normal',
                            alpha=10):
        # Add N vertices.
        self.add_vertices(*range(1, N + 1))
        # Assign latent variables to all vertices.
        if dist == 'uniform':
            alist = [random.uniform(0, 1) for v in self.vertices()]
        elif dist == 'normal':
            alist = [
                random.normalvariate(1 / 2, 1 / 6) for v in self.vertices()
            ]
        elif dist == 'exponential':
            alist = [random.expovariate(1 / 3) for v in self.vertices()]
        else:
            die('invalid latent variable distribution {dist}.')
        alist = sorted(alist)
        for v in self.vertices():
            self.set_vertex_attribute(v, 'latent', alist[v - 1])

        nedges = 0
        while nedges < E * (1 - error_ratio):
            u = self.random_vertex()
            v = self.random_vertex()
            if u == v:
                continue
            # Calcurate the connecting probability between vertices U and V.
            lu = self.get_vertex_attribute(u, 'latent')
            lv = self.get_vertex_attribute(v, 'latent')
            prob = 1.
            if confer == 'abs':
                prob = lv
            elif confer == 'binary':
                if lv <= lu:
                    prob = 0
            elif confer == 'linear':
                if lv > lu:
                    prob = lv - lu
                else:
                    prob = 0
            elif confer == 'sigmoid':
                prob = 1 / (1 + math.exp(-alpha * (lv - lu)))
            # Probabilistically connect two vertices.
            if random.random() <= prob:
                self.add_edge(u, v)
                nedges += 1

        # Add disturbance; a fraction of edges are added randomly.
        while nedges < E:
            u = self.random_vertex()
            v = self.random_vertex()
            if u == v:
                continue
            self.add_edge(u, v)
            nedges += 1
        return self

    def _lattice_vertex(self, dim, n, *positions):
        """Return the vertex located at *POSITIONS in DIM-dimensional lattice
        with N vetices per side.  For instance, the top-right corner (1, 5) in
        5x5 lattice is vertex 5, and the bottom-left corner (5, 1) in 5x5
        lattice is vertex 21."""
        v = 0
        for i in positions:
            v *= n
            if i > n:
                i -= n
            if i < 1:
                i += n
            v += i - 1
        return v + 1

    def create_lattice_graph(self, dim=2, n=5, is_torus=False):
        """Create a DIM-dimensional lattice graph with N vertices per side.
        The graph has N^DIM vertices in total.  If IS_TORUS is specified, all
        oppsosite ends are connected."""
        if dim == 1:
            for col in range(1, n + 1):
                u = self._lattice_vertex(dim, n, col)
                v = self._lattice_vertex(dim, n, col + 1)
                # Connect with the next right vertex.
                if u < v or is_torus:
                    self.add_edge(u, v)
        elif dim == 2:
            for row in range(1, n + 1):
                for col in range(1, n + 1):
                    u = self._lattice_vertex(dim, n, row, col)
                    # Connect with the next lower vertex.
                    v = self._lattice_vertex(dim, n, row + 1, col)
                    if u < v or is_torus:
                        self.add_edge(u, v)
                    # Connect with the next right vertex.
                    v = self._lattice_vertex(dim, n, row, col + 1)
                    if u < v or is_torus:
                        self.add_edge(u, v)
        return self

    def create_voronoi_graph(self, npoints=10, width=1, height=1):
        import pytess
        points = [(random.uniform(0, width), random.uniform(0, height))
                  for n in range(npoints)]
        polys = pytess.voronoi(points)
        vmax = 1
        vmap = {}
        for orig_pnt, voronoi_pnts in polys:
            for pnt in voronoi_pnts:
                if pnt not in vmap:
                    vmap[pnt] = vmax
                    vmax += 1
            last_pnt = None
            for pnt in voronoi_pnts:
                self.add_vertex(vmap[pnt])
                x, y = pnt
                # FIXME: Quick hack to pack within WIDTH x HEIGHT field.
                x = max(min(x, width), 0)
                y = max(min(y, height), 0)
                self.set_vertex_attribute(vmap[pnt], 'pos', f"+{x}+{y}")
                if last_pnt:
                    self.add_edge(vmap[last_pnt], vmap[pnt])
                last_pnt = pnt
        return self

    def create_degree_bounded_graph(self, N=10, E=20):
        """Generate a DB (Degree-Bounded) random network with N vertices and E
        edges.  For details of the algorithm, refer to K. Yamashita et al.,
        `Revisiting the Robustness of Complex Networks against Random Node
        Removal,' Journal of Information Processing, 2019."""
        k = 2 * E / N  # Average degree.
        kmin = int(k / 2)  # Minimum degree.
        if k != kmin * 2:
            die(f"Average degree {k} must be multiple of 2.")

        # Initially add N vertices.
        self.add_vertices(*range(1, N + 1))

        for u in self.vertices():
            # Randomly connect with other KMIN vertices to make sure that the
            # minimum degree is no less than KMIN.
            vlist = list(self.vertices())
            vlist.remove(u)
            for _ in range(kmin):
                v = random.choice(vlist)
                self.add_edge(u, v)
                vlist.remove(v)
        return self

    def create_configuration_graph(self, degree_seq=None):
        """Generate a graph using Newman's configuration model.  DEGREE_SEQ is
        a list of degrees for every vertex.  Different from common
        implementations, this code never generates self-loops and multi-edges."""
        def _connect_randomly(N, stubs_):
            self.__init__(directed=False)
            self.add_vertices(*range(1, N + 1))

            # Randomly connect two stubs while prohibiting self-loops and multi-edges.
            stubs = stubs_.copy()
            random.shuffle(stubs)
            while stubs:
                u = stubs.pop()
                ntries = 0
                while True:
                    v = stubs.pop()
                    if u != v and not self.has_edge(u, v):
                        self.add_edge(u, v)
                        break
                    else:
                        # FIXME: Rewrite with deque for better efficiency.
                        stubs = [v] + stubs
                        ntries += 1
                        if ntries > len(stubs):
                            return False
            return True

        if degree_seq is None:
            degree_seq = [6, 5, 4, 3, 3, 3, 2, 2, 1, 1]
        self.expect_undirected()

        # First, allocate stubs (i.e., connectors) for every vertex.  If the degree
        # of vertex v is k, STUB has the number k of v's; e.g., if the degree of vertex
        # 4 is 3, STUB contains three 4's.
        stubs = []
        for i, k in enumerate(degree_seq):
            stubs += [i + 1] * k
        if len(stubs) % 2:
            die('Total degree must be even.')

        N = len(degree_seq)
        # Loop until a realization is obtained.
        # FIXME: This code might loop indetinitely.
        while True:
            if _connect_randomly(N, stubs):
                break
        return self

    def create_random_regular_graph(self, N=10, k=3):
        degree_seq = [k] * N
        return self.create_configuration_graph(degree_seq)

    def create_li_maini_graph(self, T=200, M=4, m0=4, m=1, alpha=.1, n=1):
        """Create a graph with M clusters using the evolutionary network
        generation model: C. Li and P. K. Maini, ``An evolving network model
        with community structure,'' Journal of Physics, 2005.  

        Parameters:
          T: the number of steps
          M: the number of initial communities
          m0: the number of vertices in each initial community
          m: the number inner-edges added per every step
          alpha: the probability of adding inter-edges
          n: the number of inter-edges added per every step"""
        self.expect_undirected()

        vmax = 0
        # Community of a vertex.
        community_of = {}
        # Vertices in a community.
        vertices_in = [[] for _ in range(M)]
        # The number of links of a vertex closed within its community.
        inner_degree = {}
        # The number of links of a vertex connected with other communities.
        inter_degree = {}

        def _create_vertex_in(c):
            """Create a new vertex and initialize its attributes."""
            nonlocal vmax
            vmax += 1
            v = vmax
            self.add_vertex(v)
            inner_degree[v] = 0
            inter_degree[v] = 0
            community_of[v] = c
            vertices_in[c].append(v)
            return v

        def _add_edge(u, v):
            """Add an edge between vetex U and vertex V while updating their
            attributes."""
            self.add_edge(u, v)
            if community_of[u] == community_of[v]:
                inner_degree[u] += 1
                inner_degree[v] += 1
            else:
                inter_degree[u] += 1
                inter_degree[v] += 1

        # 1. Initialization
        # Start from a small number m0 of fully connected nodes in each
        # community.
        for c in range(M):
            for _ in range(m0):
                v = _create_vertex_in(c)
            for u, v in itertools.combinations(vertices_in[c], 2):
                _add_edge(u, v)

        # Use M(M-1)/2 inter-community links to connect each community to the
        # other M-1 communities.
        for c1, c2 in itertools.combinations(range(M), 2):
            # The nodes that the inter-community links connected to are
            # selected randomly in each community.
            u, v = random.choice(vertices_in[c1]), random.choice(
                vertices_in[c2])
            _add_edge(u, v)

        # 2. Growth
        for t in range(T):
            # At each time step, a new node is added to a randomly selected
            # community.
            c = random.randrange(M)
            u = _create_vertex_in(c)

            # 3. Preferential attachments
            # The new node will be connected to m nodes in the same community
            # through m inner-community links.
            degrees = [inner_degree[v] for v in vertices_in[c]]
            total = sum(degrees)
            prob = [d / total for d in degrees]
            for v in numpy.random.choice(vertices_in[c],
                                         size=m,
                                         replace=False,
                                         p=prob):
                _add_edge(u, v)

            # With probability alpha connected to n nodes (none with
            # probability 1 - alpha) in the other M - 1 communities through
            # inter-community links.
            if random.random() <= alpha:
                vertices = []
                degrees = []
                for v in self.vertices():
                    if community_of[v] != c:
                        vertices.append(v)
                        degrees.append(inter_degree[v])
                total = sum(degrees)
                prob = [d / total for d in degrees]
                for v in numpy.random.choice(vertices,
                                             size=m,
                                             replace=False,
                                             p=prob):
                    _add_edge(u, v)

        return self

    def create_preset_graph(self, n=0):
        self.import_edge_list("""\
1 2
2 3
2 4
3 4""".splitlines())
        return self

    def create_star_graph(self, n=10):
        for v in range(2, n + 1):
            self.add_edge(1, v)
        return self

    # import ----------------------------------------------------------------
    def import_graph(self, fmt, *args):
        name = IMPORT_SUBP.get(fmt, None)
        func = getattr(self, name, None)
        if not name or not func:
            die(f"import_graph: no import support for graph format '{fmt}'")
        return func(*args)

    def import_dot(self, lines):
        buf = ''
        for line in lines:
            # Remove C++-style comment.
            pos = line.find('//')
            if pos >= 0:
                line = line[pos:]
            line = line.strip()
            buf += line
        # Remove C-style comment.
        buf = re.sub(r'/\*.*?\*/', '', buf)
        m = re.search(r'graph\s+(\S+)\s*{(.*)}', buf)
        if not m:
            die('import_dot: invalid graph format (missing dot graph header)')
        body = m.group(2)
        return self._import_dot_body(body)

    def _import_dot_body(self, body_str):
        for line in body_str.split(';'):
            line = line.strip()
            if not line:
                continue
            if 'graph' in line or 'node' in line or 'edge' in line:
                continue

            # Must be either of graph, node, or edge definition.
            m = re.match(r'([^\[]*)\s*(\[(.*)\])?', line)
            if not m:
                continue
            val, opts = m.group(1), m.group(3) or ''
            # Parse attributes [name1=val1,name2=val2...].
            attrs = self._decode_attributes(opts)

            # Parse vertex/edge definition.
            # FIXME: This might be problematic...
            if '--' in val or '->' in val:  # vertex -- vertex [-- vertex...]
                vertices = re.split(r'\s*-[->]\s*', val)
                while len(vertices) >= 2:
                    u, v = vertices[0], vertices[1]
                    u, v = self._destringfy(u), self._destringfy(v)
                    i = self.add_edge_get_id(u, v)
                    self.set_edge_attributes_by_id(u, v, i, attrs)
                    vertices.pop(0)
            else:  # vertex
                v = self._destringfy(val)
                self.add_vertex(v)
                self.set_vertex_attributes(v, attrs)

    def import_edge_list(self, lines):
        """Creaete undirected graph from LINES in which every line describes
        an edge (a pair of vertices) with its optional weight.  A prt of the
        line after # are ignored.

        Examples:
        # unweighted K3 (complete graph with three vertices)
        1 2
        2 3
        1 3

        # weighted C4 (cubic)
        1 2 1.5
        2 3 1.1
        3 4     # default to 1
        4 1 2
"""
        for line in lines:
            line = line.rstrip()
            line = re.sub(r'#.*$', '', line)
            u, v, *rest = line.split()
            self.add_edge(u, v)
            if rest:
                w = float(rest[0])
                self.set_edge_weight(u, v, w)

    def import_cell(self, lines):
        for line in lines:
            line = line.rstrip()
            m = re.search(r'define\s+e\d+\s+link\s+v(\d+)\s+v(\d+)', line)
            if m:
                u, v = map(int, m.groups())
                self.add_edge(u, v)

    # ----------------------------------------------------------------
    def export_graph(self, fmt, *args):
        name = EXPORT_SUBP.get(fmt, None)
        func = getattr(self, name, None)
        if not name or not func:
            die(f"export_graph: no export support for graph format '{fmt}'")
        return func(*args)

    def _stringfy(self, obj):
        if obj is None:
            return ''
        elif type(obj) == str:
            return f'"{obj}"'
        elif type(obj) == float:
            return f'{obj:g}'
        else:
            return f'{obj}'

    def _destringfy(self, s):
        if not s:
            return None
        s = s.strip()
        m = re.search(r'^"(.*)"$', s)
        if m:
            return m.group(1)
        if s.startswith('0x'):
            return int(s, 16)
        if re.match(r'[\d+-]+$', s):
            return int(s)
        if re.match(r'[\d.eE+-]+$', s):
            return float(s)
        return s

    def _encode_attributes(self, attrs):
        alist = []
        for key, val in attrs.items():
            if val is None:
                continue
            astr = self._stringfy(val)
            alist.append(f'{key}={astr}')
        if not alist:
            return ''
        return ' [' + ','.join(alist) + ']'

    def _decode_attributes(self, astr):
        # Parse attributes [name1=val1,name2=val2...].
        # Note: ASTR must not contain [ and ].
        attrs = {}
        for pair in astr.split(','):
            if not pair:
                continue
            key, val = pair.split('=', 2)
            attrs[key] = self._destringfy(val)
        return attrs

    def export_dot(self, *args):
        astr = self.header_string('// ')
        head = 'digraph' if self.is_directed() else 'graph'
        astr += head + ' export_dot {\n  node [color=gray90,style=filled];\n'
        for v in sorted(self.vertices()):
            astr += '  ' + self._stringfy(v)
            attrs = self.get_vertex_attributes(v)
            astr += self._encode_attributes(attrs)
            astr += ';\n'

        for edge in sorted(self.unique_edges(), key=lambda e: e[0]):
            u = edge[0]
            v = edge[1]
            if self.undirected() and v < u:
                u, v = v, u
            l = '->' if self.is_directed() else '--'
            for i in self.get_multiedge_ids(u, v):
                astr += '  ' + self._stringfy(u) + f' {l} ' + self._stringfy(v)
                # Dump edge attributes.
                attrs = self.get_edge_attributes_by_id(u, v, i)
                astr += self._encode_attributes(attrs)
                astr += ';\n'
        astr += '}\n'
        return astr

    def export_edge_list(self, *args):
        astr = ''
        for u, v in self.edges():
            w = self.get_edge_weight(u, v)
            if w:
                astr += f'{u} {v} {w}\n'
            else:
                astr += f'{u} {v}\n'
        return astr

    def export_cell(self, *args):
        out = """\
palette tcolor 1 1 1
palette vcolor 1 .8 0 .3
palette ecolor 0 .6 1 .5
"""
        max_degree = max([self.degree(v) for v in self.vertices()])
        for v in sorted(self.vertices()):
            # NOTE: Vertices with the maximum degree are drawn with 10% height/width.
            r = .05 * self.degree(v) / max_degree
            out += f'define v{v} ellipse {r} {r} vcolor\n'
            out += f'define t{v} text {v} 12 tcolor\n'
            out += f'attach t{v} v{v}\n'
        n = 1
        for u, v in sorted(self.edges()):
            w = (self.get_edge_weight(u, v) or 1) * 2
            out += f'define e{n} link v{u} v{v} {w} ecolor\n'
            n += 1
        out += 'spring /^v/\n'
        out += 'display\n'
        out += 'wait\n'
        return out

    # Aliases.
    is_directed = directed
    is_undirected = undirected
    is_multiedged = multiedged
    add_edge_get_id = add_edge
    out_edges = edges_from
    in_edges = edges_to

def main():
    g = Graph()
    g.create_voronoi_graph()
    s = g.export_dot()
    print(s)

if __name__ == "__main__":
    main()
