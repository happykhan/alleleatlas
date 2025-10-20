#!/usr/bin/env python3
"""Build a k-NN graph either from exact distance memmaps or from approximate ANN.

Usage:
  python 02_build_knn.py --meta meta.json --out knn.npz --k 50 --mode exact
  python 02_build_knn.py --profiles profiles_int.npz --out knn.npz --k 50 --mode approx --backend hnswlib

"""
import argparse
import json
from pathlib import Path
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--meta', help='meta.json pointing to dist memmaps (exact mode)')
parser.add_argument('--profiles', help='profiles_int.npz (approx mode)')
parser.add_argument('--out', required=True)
parser.add_argument('--k', type=int, default=50)
parser.add_argument('--mode', choices=['exact','approx'], default='exact')
parser.add_argument('--backend', choices=['hnswlib','faiss','sklearn'], default='hnswlib')
args = parser.parse_args()

out = Path(args.out)

if args.mode == 'exact':
    if not args.meta:
        raise SystemExit('meta required for exact mode')
    meta = json.load(open(args.meta))
    dist0 = meta.get('dist0')
    if not dist0:
        raise SystemExit('dist0 not found in meta')
    # load condensed memmap and construct kNN by scanning rows (works for moderate N)
    import os
    import math
    filesize = os.path.getsize(dist0)
    m = filesize // 4
    n = int((1 + math.isqrt(1 + 8*m))//2)
    cond = np.memmap(dist0, dtype=np.float32, mode='r', shape=(m,))
    # For each i, find k smallest j>i/j<i
    rows = []
    cols = []
    data = []
    for i in range(n):
        # collect distances to all others
        arr = np.empty(n, dtype=float)
        for j in range(n):
            if i==j: arr[j]=np.inf
            else:
                a=min(i,j); b=max(i,j)
                idx = int(n*a - (a*(a+1))//2 + (b - a - 1))
                arr[j]=float(cond[idx])
        idxs = np.argpartition(arr, args.k)[:args.k]
        for j in idxs:
            rows.append(i); cols.append(j); data.append(arr[j])
    # make symmetric
    rows_sym = np.concatenate([rows, cols])
    cols_sym = np.concatenate([cols, rows])
    data_sym = np.concatenate([data, data])
    import scipy.sparse as sp
    knn = sp.coo_matrix((data_sym, (rows_sym, cols_sym)), shape=(n,n)).tocsr()
    sp.save_npz(out, knn)
    print('Wrote knn npz', out)

else:
    # approx mode via hnswlib
    import numpy as np
    data = np.load(args.profiles, allow_pickle=True)
    profiles = data['profiles']
    n, d = profiles.shape
    print('profiles shape', profiles.shape)
    if args.backend == 'hnswlib':
        try:
            import hnswlib
        except Exception:
            raise SystemExit('hnswlib not installed')
        # convert int profiles to float vectors (simple scaling); better to PCA or MinHash first
        X = profiles.astype(float)
        # normalize or scale as needed
        p = hnswlib.Index(space='l2', dim=d)
        p.init_index(max_elements=n, ef_construction=200, M=16)
        p.add_items(X, np.arange(n))
        p.set_ef(200)
        labels, distances = p.knn_query(X, k=args.k+1)
        # drop self
        labels = labels[:,1:]
        distances = distances[:,1:]
        import scipy.sparse as sp
        rows = np.repeat(np.arange(n), args.k)
        cols = labels.reshape(-1)
        data = distances.reshape(-1)
        M = sp.coo_matrix((data, (rows, cols)), shape=(n,n)).tocsr()
        sp.save_npz(out, M)
        print('Wrote knn npz', out)
    elif args.backend == 'sklearn':
        # fallback: use scikit-learn NearestNeighbors (works for small-to-moderate n)
        try:
            from sklearn.neighbors import NearestNeighbors
        except Exception:
            raise SystemExit('scikit-learn not available')
        X = profiles.astype(float)
        nn = NearestNeighbors(n_neighbors=min(args.k+1, n), algorithm='auto', metric='euclidean')
        nn.fit(X)
        distances, labels = nn.kneighbors(X, n_neighbors=min(args.k+1, n))
        # drop self-match at position 0 if present
        if labels.shape[1] > 0 and (labels[:,0] == np.arange(n)).all():
            labels = labels[:,1:]
            distances = distances[:,1:]
        import scipy.sparse as sp
        k_eff = labels.shape[1]
        rows = np.repeat(np.arange(n), k_eff)
        cols = labels.reshape(-1)
        data = distances.reshape(-1)
        M = sp.coo_matrix((data, (rows, cols)), shape=(n,n)).tocsr()
        sp.save_npz(out, M)
        print('Wrote knn npz (sklearn) to', out)
    else:
        raise SystemExit('faiss backend not implemented in this script')

print('done')
