import gzip
from pathlib import Path

import pytest

from alleleatlas.cluster.pHierCC import prepare_mat


def test_numeric_only(tmp_path):
    p = tmp_path / "num.tsv"
    p.write_text("ST\tl1\tl2\n1\t1\t2\n2\t3\t4\n")
    mat, names = prepare_mat(str(p))
    assert mat.dtype == int
    assert mat.shape == (2, 3)
    assert int(mat[0, 0]) == 1
    assert list(names) == [1, 2]


def test_mixed_numeric_and_string(tmp_path):
    p = tmp_path / "mixed.tsv"
    md5 = "ca2c932154669b55f0f4b6b75b8df71c"
    p.write_text("ST\tl1\tl2\nA\t1\t{}\nB\t2\t3\n".format(md5))
    mat, names = prepare_mat(str(p))
    # fallback path preserves original names
    assert list(names) == ["A", "B"]
    # numeric tokens should be preserved (1,2,3)
    numeric_values = {1, 2, 3}
    found = set()
    for v in mat[:, 1:].ravel():
        if v in numeric_values:
            found.add(int(v))
    assert found == numeric_values
    # the MD5 string should be mapped to a positive integer greater than existing numeric max
    max_numeric = 3
    # find mapped value at row for A, column l2
    mapped = int(mat[0, 2])
    assert mapped > max_numeric


def test_all_missing(tmp_path):
    p = tmp_path / "missing.tsv"
    p.write_text("ST\tl1\tl2\nX\t0\t-NA\nY\t\t0\n")
    mat, names = prepare_mat(str(p))
    # All allele columns should be 0 (missing)
    assert (mat[:, 1:] == 0).all()


def test_smoke_first1000(tmp_path):
    # Use real small dataset: cgmlst_salmonella.first1000.csv.gz
    src = Path(__file__).resolve().parents[1] / "cgmlst_data" / "cgmlst_salmonella.first1000.csv.gz"
    if not src.exists():
        pytest.skip("first1000 file not present in cgmlst_data")

    dst = tmp_path / "first1000.tsv"
    # Convert commas to tabs (prepare_mat expects tab-delimited)
    with gzip.open(src, "rt", errors="replace") as fin, open(dst, "wt") as fout:
        for line in fin:
            fout.write(line.replace(",", "\t"))

    mat, names = prepare_mat(str(dst))
    # basic sanity checks
    assert mat.shape[0] > 0
    assert mat.shape[1] > 1
