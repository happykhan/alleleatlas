import pandas as pd
import importlib.util
import os

# Load the convert_to_parquet module directly from the approx directory to avoid
# import issues when pytest adjusts sys.path
spec = importlib.util.spec_from_file_location('convert_to_parquet', os.path.join(os.path.dirname(__file__), '..', 'approx', 'convert_to_parquet.py'))
ctp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctp)


def write_text(path: str, content: str):
    with open(path, 'w') as f:
        f.write(content)


def test_drop_organismid_and_split_cgmlst(tmp_path):
    # create a CSV with organismId and a cgmlst column
    content = "ST,organismId,cgmlst\n10,orgA,1_2_3\n11,orgB,4_5_6\n"
    p = tmp_path / "o.csv"
    write_text(str(p), content)
    out = tmp_path / "out.parquet"
    ctp.main([str(p), str(out)])
    df = pd.read_parquet(str(out))
    # organismId should be dropped
    assert 'organismId' not in df.columns
    # cgmlst expanded to Locus1..3
    assert all(c in df.columns for c in ['Locus1', 'Locus2', 'Locus3'])
    assert int(df.iloc[0]['Locus2']) == 2
