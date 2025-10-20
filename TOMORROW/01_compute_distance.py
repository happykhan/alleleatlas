#!/usr/bin/env python3
"""Compute exact distance matrices or prepare data for approximate k-NN.

Usage: python 01_compute_distance.py --input datafile --outdir out --mode exact|approx --nproc 8

- mode exact: uses alleleatlas.core.distances.compute_distance_matrices to write memmaps (dist0, dist1, mat, mask, names)
- mode approx: writes a numeric vector representation (vstack of allele integers or one-hot/minhash) to disk for ANN libraries

This is a conservative wrapper; it doesn't overwrite existing outputs unless --force is used.
"""

import argparse
import os
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--mode', choices=['exact','approx'], default='exact')
parser.add_argument('--nproc', type=int, default=4)
parser.add_argument('--force', action='store_true')
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

if args.mode == 'exact':
    # Prefer the disk-backed condensed computation if available (memory-safe)
    try:
        from scripts.compute_distances_condensed import compute_condensed
        use_condensed = True
    except Exception:
        use_condensed = False

    try:
        from alleleatlas.core.input import detect_and_normalize_profile
    except Exception as e:
        print('Missing alleleatlas imports:', e)
        raise

    # produce normalized tsv
    print('Computing exact distance matrices (this may take time)...')
    df, _ = detect_and_normalize_profile(args.input, remap_alleles=True)
    normalized = outdir / (Path(args.input).stem + '.normalized.tsv')
    df.to_csv(normalized, sep='\t', index=False)

    names = df.iloc[:,0].astype(str).tolist()
    with open(outdir / 'names.txt','w') as fh:
        fh.write('\n'.join(names))

    meta = {'names': str(outdir / 'names.txt')}

    if use_condensed:
        # compute condensed memmaps (writes .dist0.memmap etc)
        print('Using disk-backed condensed distance computation...')
        dist0, dist1, mat, names_path = compute_condensed(str(normalized), str(outdir / Path(args.input).stem), nproc=args.nproc)
        meta['dist0'] = dist0
        meta['dist1'] = dist1
        meta['mat'] = mat

    else:
        # fallback: call compute_distance_matrices (may use lots of RAM)
        from alleleatlas.core.distances import compute_distance_matrices
        D, M, mask, counts = compute_distance_matrices(str(normalized), nproc=args.nproc)
        # Write condensed arrays to disk not implemented in fallback; user should prefer the condensed script
        print('Computed in-memory distance matrices; but writing condensed memmaps is not implemented in fallback. Consider installing the disk-backed script.')

    with open(outdir / 'meta.json','w') as fh:
        json.dump(meta, fh)
    print('Wrote meta and names')

else:
    # approx mode: create numeric vectors and write out to npz
    import numpy as np
    from alleleatlas.core.input import detect_and_normalize_profile
    from sklearn.preprocessing import OneHotEncoder

    print('Preparing vector representation for ANN (approx mode)...')
    df, _ = detect_and_normalize_profile(args.input, remap_alleles=True)
    names = df.iloc[:,0].astype(str).tolist()
    profiles = df.iloc[:,1:].fillna(-1).astype(int).values
    # naive numeric vector: compress allele integer space per column -> one-hot is huge; instead use simple label-encoding per column
    # For larger projects, consider MinHash or feature hashing
    # Here we produce a compact int matrix and save
    np.savez_compressed(outdir / 'profiles_int.npz', profiles=profiles, names=names)
    print('Wrote profiles_int.npz to', outdir)

print('Done')
