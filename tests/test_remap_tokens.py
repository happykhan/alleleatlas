import importlib.util
import os
import numbers

# load module by path
spec = importlib.util.spec_from_file_location('convert_to_parquet', os.path.join(os.path.dirname(__file__), '..', 'approx', 'convert_to_parquet.py'))
ctp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctp)


def write_text(path: str, content: str):
    with open(path, 'w') as f:
        f.write(content)


def test_remap_md5_and_numeric(tmp_path):
    # MD5-like ST and one md5 allele, one numeric allele
    content = "ST,L1,L2\n" + \
              "a3f5c6d8e9f12345a6b7c8d9e0f12345,100,abcd\n" + \
              "ffffffffffffffffffffffffffffffff,200,123\n"
    p = tmp_path / "r.csv"
    write_text(str(p), content)
    # use loader (remapping is applied by default)
    df, _ = ctp.load_profiles(str(p))
    # STs should be integer-like (accept numpy integer scalars)
    assert isinstance(df.loc[0, 'ST'], numbers.Integral)
    assert isinstance(df.loc[1, 'ST'], numbers.Integral)
    # L1 should be numeric for the hex-like token and L2 numeric for the second row
    assert isinstance(df.loc[0, 'L1'], numbers.Integral)
    assert isinstance(df.loc[1, 'L2'], numbers.Integral)
