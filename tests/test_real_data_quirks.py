import pandas as pd

from approx import convert_to_parquet as ctp


def write_text(path: str, content: str):
    with open(path, 'w') as f:
        f.write(content)


def test_space_separated_csv_extension(tmp_path):
    # file has .csv extension but is space-separated (a real quirk seen in samples)
    content = "ST Locus1 Locus2\n300 5 6\n301 7 8\n"
    p = tmp_path / "space.csv"
    write_text(str(p), content)
    out = tmp_path / "out_space.parquet"
    # run conversion
    ctp.main([str(p), str(out)])
    df = pd.read_parquet(str(out))
    # ST should be integer and columns preserved
    assert list(df.columns)[:3] == ['ST', 'Locus1', 'Locus2']
    assert int(df.iloc[0]['ST']) == 300


def test_preserve_many_duplicate_sts(tmp_path):
    # create a CSV where the same ST repeats many times (like Salmonella.7geneMLST sample)
    lines = ["Name,ST,L1"]
    for i in range(20):
        lines.append(f"S{i},42,1")
    # add a different ST
    lines.append("Sx,43,2")
    content = "\n".join(lines) + "\n"
    p = tmp_path / "dupes.csv"
    write_text(str(p), content)
    out = tmp_path / "out_dupes.parquet"
    ctp.main([str(p), str(out)])
    df = pd.read_parquet(str(out))
    # all rows should be preserved
    assert df.shape[0] == 21
    # duplicated ST count should equal 20 for ST==42
    dup_count = df['ST'].duplicated(keep=False).sum()
    assert int(dup_count) >= 20
