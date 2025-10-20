#!/usr/bin/env python3
"""Cluster k-NN graph, compute medoids, build MST of medoids and draw spring + Grapetree plots.

Usage:
 python 03_cluster_and_plot.py --knn path --meta meta.json --outdir out --winsor 95 --anchor-medoids

"""
import argparse
from pathlib import Path
import json
import numpy as np
import scipy.sparse as sp
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import math

parser = argparse.ArgumentParser()
parser.add_argument('--knn', required=True)
parser.add_argument('--meta', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--winsor', type=float, default=95.0)
parser.add_argument('--anchor-medoids', action='store_true')
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# load knn
knn = sp.load_npz(args.knn)
# load names
meta = json.load(open(args.meta))
names = [l.strip() for l in open(meta['names'])]

# cluster (Leiden if available, else connected components)
try:
    import igraph as ig
    import leidenalg
    use_leiden = True
except Exception:
    use_leiden = False

# Build NetworkX graph from sparse matrix in a way that works across NX versions
try:
    G = nx.from_scipy_sparse_matrix(knn)
except AttributeError:
    try:
        # newer networkx may expose from_scipy_sparse_array
        G = nx.from_scipy_sparse_array(knn)
    except Exception:
        # last resort: convert to dense (small test inputs only)
        G = nx.from_numpy_array(knn.toarray())

if use_leiden:
    import numpy as np
    import igraph as ig
    import leidenalg
    # build igraph from csr: use nonzero edges and weights
    sources, targets = knn.nonzero()
    weights = knn.data
    g = ig.Graph(n=knn.shape[0], edges=list(zip(sources.tolist(), targets.tolist())), directed=False)
    g.es['weight'] = weights.tolist()
    part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights='weight', resolution_parameter=1.0)
    labels = np.array(part.membership)
    print('Used Leiden clustering')
else:
    # connected components
    comp = list(nx.connected_components(G))
    labels = np.full(len(names), -1, dtype=int)
    for i, c in enumerate(comp):
        for v in c:
            labels[v] = i
    print('Used connected-components for clustering')

# medoids: per-cluster choose node with minimal sum distance to cluster members (approx via knn)
clusters = {}
for i,l in enumerate(labels):
    clusters.setdefault(l, []).append(i)
medoids = []
for l, members in clusters.items():
    # use degree heuristic: choose highest degree within knn
    degs = knn[members].getnnz(axis=1)
    med = members[int(np.argmax(degs))]
    medoids.append(med)

# write clusters and medoids
with open(outdir / 'clusters.tsv','w') as fh:
    for i,lab in enumerate(labels):
        fh.write(f"{names[i]}\t{lab}\n")
with open(outdir / 'medoids.tsv','w') as fh:
    for i,m in enumerate(medoids):
        fh.write(f"{i}\t{m}\t{names[m]}\n")

# build sample MST from knn
from scipy.sparse.csgraph import minimum_spanning_tree
mst = minimum_spanning_tree(knn).tocoo()
G_mst = nx.Graph()
G_mst.add_nodes_from(range(len(names)))
for u,v,w in zip(mst.row, mst.col, mst.data):
    G_mst.add_edge(int(u), int(v), weight=float(w))

# connect clusters via medoid MST
k = len(medoids)
M = np.zeros((k,k), dtype=float)
# read condensed dist (if available) else approximate via knn
meta_j = json.load(open(args.meta))
if 'dist0' in meta_j:
    dist0 = meta_j['dist0']
    import os, math
    filesize=os.path.getsize(dist0)
    m = filesize // 4
    n = int((1 + math.isqrt(1 + 8*m))//2)
    cond = np.memmap(dist0, dtype=np.float32, mode='r', shape=(m,))
    def idxc(i,j):
        if i==j: return None
        if i>j: i,j=j,i
        return int(n*i - (i*(i+1))//2 + (j - i - 1))
    for a,i in enumerate(medoids):
        for b,jdx in enumerate(medoids):
            if a==b: continue
            M[a,b] = float(cond[idxc(i,jdx)])
else:
    # fallback: use shortest path length in knn
    for a,i in enumerate(medoids):
        lengths = nx.single_source_shortest_path_length(G, i)
        for b,jdx in enumerate(medoids):
            M[a,b] = lengths.get(jdx, 1e9)

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
mst_med = minimum_spanning_tree(csr_matrix(M)).tocoo()
# add medoid edges to sample MST
for u,v,w in zip(mst_med.row, mst_med.col, mst_med.data):
    i = medoids[int(u)]; j = medoids[int(v)]; w = float(w)
    if G_mst.has_edge(i,j):
        if G_mst[i][j].get('weight', float('inf')) > w:
            G_mst[i][j]['weight']=w; G_mst[i][j]['centroid']=True
    else:
        G_mst.add_edge(i,j, weight=w, centroid=True)

# winsorize edges
weights = np.array([d['weight'] for u,v,d in G_mst.edges(data=True)])
limit = np.percentile(weights, args.winsor)
for u,v,d in G_mst.edges(data=True):
    d['w'] = min(d['weight'], limit)

# optional anchor medoids
pos = None
if args.anchor_medoids:
    # compute Grapetree positions, use medoids as fixed anchors
    try:
        from scripts.draw_mst import StaticGrapeTreeLayout
    except Exception:
        # try alleleatlas.scripts or load by path
        try:
            from alleleatlas.scripts.draw_mst import StaticGrapeTreeLayout
        except Exception:
            import importlib.util, pathlib
            p = pathlib.Path(__file__).resolve().parents[1] / 'scripts' / 'draw_mst.py'
            spec = importlib.util.spec_from_file_location('draw_mst_local', str(p))
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            StaticGrapeTreeLayout = mod.StaticGrapeTreeLayout
    layout = StaticGrapeTreeLayout(G_mst)
    pos = layout.run()
    fixed = [m for m in medoids]
    print('Anchoring medoids:', fixed)
    pos = nx.spring_layout(G_mst, pos=pos, fixed=fixed, iterations=200, seed=42)
else:
    pos = nx.spring_layout(G_mst, iterations=200, seed=42)

# save images
cmap = plt.get_cmap('tab20')
label_to_color = {lab: cmap(i%20) for i,lab in enumerate(sorted(clusters.keys()))}
node_colors = [label_to_color.get(labels[i], (0.6,0.6,0.6,1.0)) for i in range(len(names))]
node_sizes = [max(8, math.log(1+len(clusters[labels[i]]))*25) for i in range(len(names))]

plt.figure(figsize=(12,12))
for u,v,data in G_mst.edges(data=True):
    x1,y1 = pos[u]; x2,y2 = pos[v]
    c = 'red' if data.get('centroid') else 'gray'
    plt.plot([x1,x2],[y1,y2], color=c, linewidth=0.6, alpha=0.6)
xs=[pos[i][0] for i in range(len(names))]
ys=[pos[i][1] for i in range(len(names))]
plt.scatter(xs, ys, s=node_sizes, c=node_colors, linewidths=0)
plt.axis('off')
plt.title('MST - spring layout')
plt.savefig(outdir / 'mst_spring.png', dpi=200, bbox_inches='tight')
plt.close()

# Grapetree
pos_g = layout.run()
try:
    from scripts.draw_mst import StaticGrapeTreeLayout as _SGL
except Exception:
    try:
        from alleleatlas.scripts.draw_mst import StaticGrapeTreeLayout as _SGL
    except Exception:
        import importlib.util, pathlib
        p = pathlib.Path(__file__).resolve().parents[1] / 'scripts' / 'draw_mst.py'
        spec = importlib.util.spec_from_file_location('draw_mst_local', str(p))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _SGL = mod.StaticGrapeTreeLayout

layout = _SGL(G_mst)
pos_g = layout.run()
# blend
alpha = 0.6
pos_blend = {i: (alpha*pos_g[i][0] + (1-alpha)*pos[i][0], alpha*pos_g[i][1] + (1-alpha)*pos[i][1]) for i in range(len(names))}
plt.figure(figsize=(12,12))
for u,v,data in G_mst.edges(data=True):
    x1,y1 = pos_blend[u]; x2,y2 = pos_blend[v]
    c = 'red' if data.get('centroid') else 'gray'
    plt.plot([x1,x2],[y1,y2], color=c, linewidth=0.6, alpha=0.6)
xs=[pos_blend[i][0] for i in range(len(names))]
ys=[pos_blend[i][1] for i in range(len(names))]
plt.scatter(xs, ys, s=node_sizes, c=node_colors, linewidths=0)
plt.axis('off')
plt.title('MST - Grapetree-blend')
plt.savefig(outdir / 'mst_grapetree_blend.png', dpi=200, bbox_inches='tight')
plt.close()

print('Wrote outputs to', outdir)
