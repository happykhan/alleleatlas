import pandas as pd

from alleleatlas import convert_to_parquet as ctp


def write_text(path: str, content: str):
    with open(path, 'w') as f:
        f.write(content)


def test_tsv_input(tmp_path):
    content = """#Name\tST\tLocus1\tLocus2
sampleA\t120\t1\t2
sampleB\t121\t3\t4
"""
    p = tmp_path / "in.tsv"
    write_text(str(p), content)
    out = tmp_path / "out.parquet"
    ctp.main([str(p), str(out)])
    df = pd.read_parquet(str(out))
    assert list(df.columns)[0] == 'ST'
    assert df.shape == (2, 3)


def test_csv_gz_input_and_duplicates(tmp_path):
    import gzip
    content = "Name,ST,LocusA,LocusB\nS1,10,1,2\nS1,10,1,3\n"
    p = tmp_path / "in.csv.gz"
    with gzip.open(str(p), 'wt') as f:
        f.write(content)
    out = tmp_path / "out2.parquet"
    ctp.main([str(p), str(out)])
    df = pd.read_parquet(str(out))
    # duplicates should be preserved
    assert df.shape[0] == 2
    assert 'ST' in df.columns


def test_metadata_like_input(tmp_path):
    content = "ST Locus1 Locus2\n200 5 6\n"
    p = tmp_path / "meta.txt"
    write_text(str(p), content)
    out = tmp_path / "out3.parquet"
    ctp.main([str(p), str(out)])
    df = pd.read_parquet(str(out))
    assert int(df.iloc[0]['ST']) == 200


def test_name_map_output(tmp_path):
    # combined Name and ST: mapping should be produced when --name-map is given
    content = "Name,ST,L1\nA,50,1\nB,51,2\n"
    p = tmp_path / "nm.csv"
    write_text(str(p), content)
    out = tmp_path / "out_nm.parquet"
    map_out = tmp_path / "name_map.csv"
    ctp.main([str(p), str(out), '--name-map', str(map_out)])
    # mapping file exists and has correct pairs
    mp = pd.read_csv(str(map_out))
    assert set(mp['Name']) == {'A', 'B'}
    assert set(mp['ST'].astype(int)) == {50, 51}
