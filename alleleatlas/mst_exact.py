"""Exact MST utilities for cgMLST distances.

Implements Kruskal's algorithm with a deterministic tie-break key:
sort by (distance, node_label_min, node_label_max) so behavior matches
Grapetree's deterministic selection when equal distances occur.

Exports:
- build_mst_from_distance_matrix(D, labels=None) -> networkx.Graph
- mst_to_newick(G, labels=None) -> str
"""
from __future__ import annotations

import networkx as nx
import numpy as np
from typing import Optional, Sequence, List, Tuple


def _edge_list_from_dist_matrix(D: np.ndarray, labels: Optional[Sequence[str]] = None) -> List[Tuple[float, int, int]]:
    n = D.shape[0]
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            dij = D[i, j]
            if np.isnan(dij):
                continue
            edges.append((float(dij), i, j))
    return edges


def build_mst_from_distance_matrix(D: np.ndarray, labels: Optional[Sequence[str]] = None) -> nx.Graph:
    """Build MST (undirected) from full distance matrix D.

    labels: optional sequence of string labels for nodes; used for
    lexicographic tie-breaking when distances equal. Returns an
    undirected NetworkX Graph with 'weight' attribute on edges.
    """
    n = D.shape[0]
    G = nx.Graph()
    for i in range(n):
        G.add_node(i)

    edges = _edge_list_from_dist_matrix(D, labels)

    # Sorting key: distance, then lexicographic node labels (min,max)
    def sort_key(e):
        dist, i, j = e
        if labels is not None:
            a = labels[i]
            b = labels[j]
        else:
            a = str(i)
            b = str(j)
        if a <= b:
            return (dist, a, b)
        else:
            return (dist, b, a)

    edges_sorted = sorted(edges, key=sort_key)

    # Union-Find (Disjoint Set)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra = find(a)
        rb = find(b)
        parent[rb] = ra

    for dist, i, j in edges_sorted:
        if find(i) != find(j):
            G.add_edge(i, j, weight=dist)
            union(i, j)
        # stop early when we have n-1 edges
        if G.number_of_edges() >= n - 1:
            break

    return G


def mst_to_newick(G: nx.Graph, labels: Optional[Sequence[str]] = None) -> str:
    """Serialize an (unrooted) tree graph to Newick string.

    We choose node 0 as root and produce a rooted Newick. For leaf nodes
    we output the label (or index). Branch lengths are taken from edge
    'weight' attributes when present.
    """
    # Build rooted tree with node 0 as root
    root = 0
    T = nx.DiGraph()

    # BFS to orient edges away from root
    from collections import deque
    q = deque([root])
    visited = {root}
    while q:
        u = q.popleft()
        for v in G.neighbors(u):
            if v in visited:
                continue
            visited.add(v)
            w = G[u][v].get('weight', 0.0)
            T.add_edge(u, v, weight=w)
            q.append(v)

    def node_label(i):
        return labels[i] if labels is not None else str(i)

    def recurse(u):
        children = list(T.successors(u))
        if not children:
            return node_label(u)
        parts = []
        for v in children:
            w = T[u][v].get('weight', 0.0)
            parts.append(f"{recurse(v)}:{w}")
        return '(' + ','.join(parts) + ')' + node_label(u)

    newick = recurse(root) + ';'
    return newick


if __name__ == '__main__':
    # quick self-test using a toy distance matrix
    import numpy as _np
    D = _np.array([[0,1,2],[1,0,3],[2,3,0]], dtype=float)
    G = build_mst_from_distance_matrix(D)
    print(mst_to_newick(G))
