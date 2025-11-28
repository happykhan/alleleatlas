import numpy as np
import scipy.sparse as sp
from pathlib import Path


def make_small_knn(path: str, n=20, k=3):
    rng = np.random.default_rng(1)
    X = rng.random((n, 8))
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=min(k+1, n), algorithm='auto').fit(X)
    dists, labels = nn.kneighbors(X, n_neighbors=min(k+1, n))
    if labels.shape[1] > 0 and (labels[:, 0] == np.arange(n)).all():
        labels = labels[:, 1:]
        dists = dists[:, 1:]
    rows = np.repeat(np.arange(n), labels.shape[1])
    cols = labels.reshape(-1)
    data = dists.reshape(-1)
    M = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    sp.save_npz(path, M)


def test_analyze_knn_cli(tmp_path):
    knn_file = str(tmp_path / 'small_knn.npz')
    make_small_knn(knn_file, n=30, k=4)
    out_prefix = str(tmp_path / 'out')
    # run analyze CLI (call Typer function directly)
    from alleleatlas.analyze_knn import run_analysis
    run_analysis(knn_file, out_prefix)
    # verify expected files
    assert Path(out_prefix + '_dist_hist.png').exists()
    assert Path(out_prefix + '_deg_hist.png').exists()
    assert Path(out_prefix + '_components.png').exists()


def test_mst_and_plot_cli(tmp_path):
    knn_file = str(tmp_path / 'small_knn.npz')
    make_small_knn(knn_file, n=25, k=4)
    out_img = str(tmp_path / 'mst.png')
    from alleleatlas.mst_and_plot import cli
    cli(knn_file, out_img)
    assert Path(out_img).exists()
