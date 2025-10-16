#!/usr/bin/env python3
"""
Split cgmlst_supp_species.csv.gz into per-organism gzipped CSVs.
Writes files to the same directory as the input, named: cgmlst_<organismid>.csv.gz

Usage:
    python scripts/split_cgmlst_supp_species.py 
    python scripts/split_cgmlst_supp_species.py --input cgmlst_data/cgmlst_supp_species.csv.gz --outdir cgmlst_data --overwrite
"""
import argparse
import gzip
import re
from pathlib import Path
import sys
import csv


def sanitize_filename(s: str) -> str:
    # Keep alphanumerics, dot, dash, underscore; replace others with underscore
    s = str(s)
    s = s.strip()
    s = re.sub(r"[\s/\\]+", "_", s)
    s = re.sub(r"[^\w.\-]", "_", s)
    s = re.sub(r"_+", "_", s)
    return s[:200]


def find_organism_col(cols):
    lowered = {c.lower(): c for c in cols}
    for candidate in ("organismid", "organism_id", "organism", "organismid1"):
        if candidate in lowered:
            return lowered[candidate]
    for c in cols:
        if 'organism' in c.lower():
            return c
    return None


def split_file_stream(input_path: Path, outdir: Path, overwrite: bool = False):
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"Reading (stream) {input_path} ...")
    outdir.mkdir(parents=True, exist_ok=True)

    open_writers = {}
    written = []
    total = 0

    # allow very large fields
    try:
        csv.field_size_limit(sys.maxsize)
    except Exception:
        pass

    with gzip.open(input_path, 'rt', newline='') as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("Could not determine CSV header/columns from input")
        print(f"Columns detected: {len(reader.fieldnames)}: {reader.fieldnames[:10]}")
        org_col = find_organism_col(reader.fieldnames)
        if org_col is None:
            raise ValueError(f"Could not find an organism column in input. Available columns: {reader.fieldnames[:20]}")
        print(f"Using organism column: '{org_col}'")

        # stream rows and write to per-organism gzip CSVs lazily
        for row in reader:
            total += 1
            orgid = row.get(org_col, '')
            if orgid is None:
                orgid = ''
            orgkey = str(orgid).strip()
            if orgkey == '':
                fname = outdir / "cgmlst_unknown_organism.csv.gz"
                safe = 'unknown'
            else:
                safe = sanitize_filename(orgkey)
                fname = outdir / f"cgmlst_{safe}.csv.gz"

            if fname in open_writers:
                writer, gzfh = open_writers[fname]
                # writer may be None when we decided to skip because file exists and overwrite=False
                if writer is None:
                    # skipping rows for this organism because file existed and overwrite=False
                    continue
                writer.writerow(row)
            else:
                # if file exists and not overwrite, skip writing any rows for this group
                if fname.exists() and not overwrite:
                    open_writers[fname] = (None, None)  # marker for skipped
                    written.append((orgkey, fname, 'skipped'))
                else:
                    gzfh = gzip.open(fname, 'wt', newline='')
                    writer = csv.DictWriter(gzfh, fieldnames=reader.fieldnames)
                    writer.writeheader()
                    writer.writerow(row)
                    open_writers[fname] = (writer, gzfh)

    # close open file handles
    for fname, (writer, gzfh) in open_writers.items():
        if gzfh:
            gzfh.close()
            # count rows for summary is expensive; indicate written
            written.append((fname.name, fname, 'written'))

    print(f"Processed {total:,} rows")
    return written


def main(argv=None):
    p = argparse.ArgumentParser(description="Split cgmlst_supp_species.csv.gz into per-organism files")
    p.add_argument('--input', '-i', default='cgmlst_data/cgmlst_supp_species.csv.gz')
    p.add_argument('--outdir', '-o', default='cgmlst_data')
    p.add_argument('--overwrite', action='store_true', help='Overwrite existing output files')
    args = p.parse_args(argv)

    input_path = Path(args.input)
    outdir = Path(args.outdir)

    try:
        written = split_file_stream(input_path, outdir, overwrite=args.overwrite)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(2)

    # summary
    print('\nSummary:')
    for orgid, fname, status in written:
        print(f"{status.upper():8} {orgid!r} -> {fname}")


if __name__ == '__main__':
    main()
