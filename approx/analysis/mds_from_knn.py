from alleleatlas.mst_and_plot import compute_embedding as compute_embedding  # reuse existing helper if present

def load_knn_npz(path: str):
    import scipy.sparse as sp
    return sp.load_npz(path).tocsr()
