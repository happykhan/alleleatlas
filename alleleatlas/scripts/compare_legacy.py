"""Compare new mst_workflow distances and MST against legacy old_code/MSTrees.py output.

Usage: run from repository root:
    python -m alleleatlas.scripts.compare_legacy
"""
from __future__ import annotations
import importlib.util
import numpy as np
from pathlib import Path


def load_module_from_path(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_legacy_module(path):
    spec = importlib.util.spec_from_file_location('legacy_mstrees', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_legacy_distance(text: str):
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return None
    # first line is count
    try:
        header = lines[0].strip().split()
        n = int(header[-1])
        mat = np.zeros((n, n), dtype=float)
        for i, ln in enumerate(lines[1: 1 + n]):
            parts = ln.strip().split()
            # first token is name
            values = [float(x) for x in parts[1:1 + n]]
            mat[i, :] = values
        return mat
    except (ValueError, IndexError):
        return None


def main():
    repo_root = Path(__file__).resolve().parents[2]
    profile_files = list((repo_root / 'cgmlst_data').glob('*.profile'))
    if not profile_files:
        profile_files = list((repo_root / 'cgmlst_data').glob('*.profile'))
    if not profile_files:
        print('No .profile file found in cgmlst_data')
        return
    profile = str(profile_files[0])
    print('Using profile:', profile)

    legacy = load_legacy_module(str(repo_root / 'old_code' / 'MSTrees.py'))
    mst_workflow = load_module_from_path('mst_workflow', str(repo_root / 'alleleatlas' / 'mst_workflow.py'))
    # get legacy distance text
    print('Computing legacy distance (may take a while)...')
    # force single-process to avoid multiprocessing pickling issues when loading module dynamically
    legacy_text = legacy.backend(profile=profile, method='distance', n_proc=1)
    legacy_dist = parse_legacy_distance(legacy_text)
    if legacy_dist is None:
        print('Failed to parse legacy distances')
        return

    # compute new distances (asymmetric default to mirror MSTreeV2)
    names, profiles = mst_workflow.read_profile(profile)
    print('Computing new distances...')
    new_dist = mst_workflow.distance_asymmetric(profiles)

    # symmetrize for direct comparison
    legacy_sym = (legacy_dist + legacy_dist.T) / 2.0
    new_sym = (new_dist + new_dist.T) / 2.0

    # compare shapes and statistics
    print('legacy_sym shape', legacy_sym.shape, 'new_sym shape', new_sym.shape)
    m = min(legacy_sym.shape[0], new_sym.shape[0])
    diff = np.abs(legacy_sym[:m, :m] - new_sym[:m, :m])
    print('mean abs diff:', float(np.mean(diff)))
    print('max abs diff:', float(np.max(diff)))

    # build MST from new_sym and save figure
    print('Building MST and plotting...')
    # weights heuristic: use harmonic as default
    embeded = [[n] for n in names]
    weights = mst_workflow.heuristic_weights(new_sym, embeded, method='harmonic')
    edges = mst_workflow.build_mst(new_sym, weights, matrix_type='symmetric')
    # collapse at a few example cutoffs
    node2comp, comp_edges, comps = mst_workflow.collapse_components(edges, cutoff=1.0)
    out = mst_workflow.plot_mst(edges, names, node2comp=node2comp, comps=comps, out_file=str(repo_root / 'mst_new.png'))
    print('Wrote new MST plot to', out)

    # save distance matrices for inspection
    np.save(str(repo_root / 'legacy_dist.npy'), legacy_sym)
    np.save(str(repo_root / 'new_dist.npy'), new_sym)
    print('Saved legacy_dist.npy and new_dist.npy in repo root')


if __name__ == '__main__':
    main()
