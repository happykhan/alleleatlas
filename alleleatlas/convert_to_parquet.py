"""approx.convert_to_parquet
=================================

Utilities to normalize allele profile files into a standardized compressed
Parquet table with the column order:

        ST, Locus1, Locus2, ...

This module supports a variety of real-world quirks found in allele profile
files (CSV, TSV, space-delimited "metadata" files, and compressed files
ending in ``.gz`` or ``.xz``). It also normalizes missing values and
deterministically remaps any non-numeric allele tokens (for example MD5-like
strings) to stable integers.

Features
- Detect input format by peeking at the first non-empty line(s).
- Remove leading ``#`` from header lines when present.
- Preserve duplicate ST rows and report a warning with the duplicate count.
- Split a single ``cgmlst`` underscore-delimited column into ``Locus1..N``.
- Normalize missing tokens (blank, ``_``, ``nan``, ``none``, ``-1``) to ``0``.
- Deterministic token -> int mapping using hex parsing or SHA256 fallback.

CLI and developer API
- The command-line user interface uses :mod:`typer` and prints colored
    messages using :mod:`rich`.
- A lightweight programmatic entrypoint ``main(argv)`` is preserved so tests
    and other code can call this module directly without invoking the Typer
    CLI.

Example
-------
Run on the command line::

        python -m approx.convert_to_parquet input.csv out.parquet

Or call programmatically in tests or scripts::

        from approx import convert_to_parquet as ctp
        ctp.main(["input.csv", "out.parquet"])  # kept for backward compatibility

"""
from __future__ import annotations

import os
import gzip
import lzma
import argparse
import tempfile
from typing import List, Optional, Tuple

import pandas as pd
import typer
from rich.console import Console


# Rich console for styled output (used across the module)
console = Console()

# Typer app used when running as a CLI (kept separate from programmatic main)
app = typer.Typer(help="Convert allele profile files into standardized Parquet")


def open_maybe_gz(path: str, encoding: str = 'utf-8'):
    """Open a file that may be plain, gzipped or lzma-compressed.

    Returns a text-mode file object. This helper centralizes compression
    handling so callers don't need to know about gzip/lzma.
    """
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding=encoding)
    if path.endswith('.xz') or path.endswith('.lzma'):
        return lzma.open(path, 'rt', encoding=encoding)
    return open(path, 'rt', encoding=encoding)


def read_first_nonempty_lines(path: str, n: int = 10) -> List[str]:
    """Return up to ``n`` first non-empty lines from ``path``.

    The function reads compressed or plain files using :func:`open_maybe_gz`.
    Useful for format detection without loading the whole file into memory.
    """
    lines: List[str] = []
    with open_maybe_gz(path) as fh:
        for _ in range(n):
            line = fh.readline()
            if not line:
                break
            lines.append(line.rstrip('\n'))
    return lines


def detect_format(path: str) -> str:
    """Detect a simple, best-effort format tag for the input file.

    Returns one of: ``'csv'``, ``'tsv'``, ``'metadata'`` (space-delimited),
    or ``'unknown'``. This function peeks at the first non-empty line and
    uses simple heuristics (presence of commas, tabs, or a leading ``ST``/``Name``
    token) to make a decision.
    """
    lines = read_first_nonempty_lines(path, 20)
    # remove comment markers and whitespace for detection
    stripped = [line.lstrip('#').strip() for line in lines if line.strip()]
    if not stripped:
        return 'unknown'
    header = stripped[0]
    # tab-delimited header -> tsv
    if '\t' in header:
        return 'tsv'
    # comma -> csv
    if ',' in header:
        return 'csv'
    # space-separated files often start with ST or Name
    if header.lower().startswith('st') or header.lower().startswith('name'):
        return 'metadata'
    # fallback
    return 'unknown'



def load_profiles(path: str) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """Load and normalize allele profiles from ``path``.

    Returns a tuple ``(df, name_map_df)`` where ``df`` is a DataFrame with all
    columns converted to integers (``int64``) and ``name_map_df`` is an
    optional DataFrame mapping names to STs when both ``Name`` and ``ST`` are
    present in the input.
    """
    fmt = detect_format(path)
    # report format to the rich console
    console.print(f"[blue]Detected format:[/blue] {fmt}")
    # If the header line itself begins with a '#', pandas will skip it when using
    # comment='#'. Preprocess such files into a temporary cleaned file where the
    # first header-like line has its leading '#' removed, then read via pandas.
    cleaned_path = None
    try:
        with open_maybe_gz(path) as fh:
            # peek at first few lines to see whether the header line starts with '#'
            lines = []
            for _ in range(5):
                ln = fh.readline()
                if not ln:
                    break
                lines.append(ln)

        first_nonempty = None
        for ln in lines:
            if ln.strip():
                first_nonempty = ln
                break

        if first_nonempty and first_nonempty.lstrip().startswith('#'):
            # create a temporary file with only the header cleaned (strip leading '#')
            tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt')
            cleaned_path = tmp.name
            done_header = False
            with open_maybe_gz(path) as fh_in:
                for ln in fh_in:
                    if not done_header and ln.lstrip().startswith('#'):
                        # heuristics: treat this as header if after removing '#' it looks like a header
                        after = ln.lstrip('#')
                        if '\t' in ln or ',' in ln or after.strip().lower().startswith('st') or 'st' in after.lower().split():
                            tmp.write(after)
                            done_header = True
                            continue
                    tmp.write(ln)
            tmp.close()
            read_path = cleaned_path
        else:
            read_path = path

    # read with appropriate pandas reader
        if fmt == 'tsv':
            df = pd.read_csv(read_path, sep='\t', comment='#', dtype=str)
        elif fmt == 'csv':
            # try reading csv; pandas can handle gz compressed directly
            df = pd.read_csv(read_path, comment='#', dtype=str)
        elif fmt == 'metadata':
            # try whitespace delim (some files use spaces) but if it fails fallback to csv
            try:
                df = pd.read_csv(read_path, sep=r"\s+", comment='#', dtype=str)
            except Exception:
                df = pd.read_csv(read_path, comment='#', dtype=str)
        else:
            # fallback: try pandas auto-detect
            try:
                df = pd.read_csv(read_path, comment='#', dtype=str)
            except Exception:
                df = pd.read_table(read_path, comment='#', dtype=str)
    finally:
        # remove temporary cleaned file if created
        if cleaned_path is not None:
            try:
                os.unlink(cleaned_path)
            except Exception:
                pass

    # Normalize column names: strip '#' and surrounding whitespace
    df.columns = [c.strip().lstrip('#').strip() for c in df.columns]

    # Ensure there's an ST column (case-insensitive match)
    col_lower = {c.lower(): c for c in df.columns}
    has_st = 'st' in col_lower
    has_name = 'name' in col_lower

    name_map_df = None

    if has_st:
        st_col = col_lower['st']
        # If both Name and ST exist, create mapping and drop the Name column
        if has_name:
            name_col = col_lower['name']
            # create mapping dataframe (preserve all rows)
            name_map_df = df[[name_col, st_col]].rename(columns={name_col: 'Name', st_col: 'ST'})
            # ensure mapping columns are strings
            name_map_df = name_map_df.astype(str)
            df = df.drop(columns=[name_col])
    elif has_name:
        st_col = col_lower['name']
        # rename chosen ST column to 'ST'
        if st_col != 'ST':
            df = df.rename(columns={st_col: 'ST'})
    else:
        # if ST isn't present, try first column
        st_col = df.columns[0]
        if st_col != 'ST':
            df = df.rename(columns={st_col: 'ST'})

    # ensure we have 'ST' in columns now
    if 'ST' not in df.columns:
        df = df.rename(columns={st_col: 'ST'})

    # Convert all data to string (alleles may be numbers or strings)
    df = df.astype(str)

    # Normalize missing data tokens to '0' for consistency.
    # Treat blank, underscore, 'nan', 'none', and '-1' as missing -> '0'.
    df = df.replace({
        r'^\s*$': '0',
        r'^_$': '0',
        r'(?i)^nan$': '0',
        r'(?i)^none$': '0',
        r'^-1$': '0'
    }, regex=True)

    # Drop organismId column if present (case-insensitive)
    lower_cols = {c.lower(): c for c in df.columns}
    if 'organismid' in lower_cols:
        df = df.drop(columns=[lower_cols['organismid']])

    # If a single 'cgmlst' column exists where all loci are encoded as
    # underscore-separated values, split it into multiple locus columns
    # named Locus1, Locus2, ... and insert them after 'ST'.
    cg_col = None
    for c in df.columns:
        if c.lower() == 'cgmlst':
            cg_col = c
            break
    if cg_col is not None:
        # split into separate columns
        split_df = df[cg_col].fillna('').astype(str).str.split('_', expand=True)
        # name the new columns Locus1..N
        new_names = [f'Locus{i+1}' for i in range(split_df.shape[1])]
        split_df.columns = new_names
        # drop original cgmlst column and concat the split columns
        df = df.drop(columns=[cg_col])
        # ensure split columns are strings
        split_df = split_df.astype(str)
        # insert split_df columns after ST
        cols = list(df.columns)
        if 'ST' in cols:
            st_index = cols.index('ST')
            before = cols[:st_index+1]
            after = cols[st_index+1:]
            df = pd.concat([df[before].reset_index(drop=True), split_df.reset_index(drop=True), df[after].reset_index(drop=True)], axis=1)
        else:
            df = pd.concat([split_df.reset_index(drop=True), df.reset_index(drop=True)], axis=1)

    # Reorder columns: ST first, then the remaining loci in original order
    cols = list(df.columns)
    others = [c for c in cols if c != 'ST']
    ordered = ['ST'] + others
    # If ST doesn't exist (shouldn't happen now), keep original order
    if 'ST' in df.columns:
        df = df[ordered]

    # Report duplicate STs count (but keep them). Use rich to make the
    # warning more visible in CLI runs.
    dup_count = df['ST'].duplicated(keep=False).sum()
    if dup_count:
        console.print(f"[yellow]Warning:[/yellow] {dup_count} duplicate ST rows detected (kept as-is)")

    # Remap/convert all allele tokens (and ST) to integer IDs in a vectorized
    # manner. For performance we map unique tokens per-column and then apply
    # the mapping using pandas' map (much faster than applymap on large frames).
    # Remap/convert all allele tokens (and ST) to integer IDs in a vectorized
    # manner. For performance we map unique tokens per-column and then apply
    # the mapping using pandas' map (much faster than applymap on large frames).
    import re
    import hashlib

    MOD = 2147483647
    hex_re = re.compile(r'^[0-9a-f]{8,}$', re.IGNORECASE)

    def _token_to_int(tok) -> int:
        """Deterministically map a single token to an integer.

        Rules:
        - Missing-like tokens -> 0
        - Pure decimal strings -> int(value)
        - Hex-like long tokens (e.g. MD5) -> parse leading 16 hex chars -> int % MOD
        - Fallback -> sha256(token) -> take first 8 bytes -> int % MOD
        """
        if tok is None:
            return 0
        s = str(tok).strip()
        if s == '' or s == '_' or s.lower() in ('nan', 'none') or s == '-1':
            return 0
        # decimal
        if s.isdigit():
            try:
                return int(s)
            except Exception:
                pass
        # hex-like tokens: parse a chunk to reduce collisions and size
        if hex_re.match(s):
            try:
                return int(s[:16], 16) % MOD
            except Exception:
                pass
        # deterministic hash fallback
        try:
            h = hashlib.sha256(s.encode('utf-8')).digest()
            return int.from_bytes(h[:8], 'big') % MOD
        except Exception:
            return 0

    # Build mapping per-column using unique values for speed (vectorized map)
    for col in df.columns:
        uniques = pd.Series(df[col].unique()).astype(str)
        mapping = {u: _token_to_int(u) for u in uniques}
        df[col] = df[col].map(mapping)

    # Ensure integer dtype (int64) for all columns
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0).astype('int64')

    return df, name_map_df


def write_parquet(df: pd.DataFrame, out: str) -> None:
    """Write DataFrame to a Parquet file using pyarrow and snappy compression.

    Uses :func:`rich.console.Console.print` for the final message so CLI
    output is styled.
    """
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    df.to_parquet(out, index=False, engine='pyarrow', compression='snappy')
    console.print(f"[green]Wrote parquet:[/green] {out} ({len(df)} rows, {len(df.columns)} columns)")


def convert_file(
    input_path: str,
    output_path: str,
    name_map: Optional[str] = None,
    dry_run: bool = False,
    verbose: bool = False,
) -> int:
    """High-level conversion routine used by both the CLI and tests.

    This function loads profiles, optionally prints verbose information to
    the console, writes the parquet file (unless ``dry_run``) and writes the
    optional name->ST mapping.
    """
    df, name_map_df = load_profiles(input_path)
    if verbose:
        # use a small rich table to summarize columns
        console.print(f"[blue]Inferred columns:[/blue] {list(df.columns)}")
        console.print(f"[blue]Rows:[/blue] {len(df)}")
    if dry_run:
        console.print("[italic]Dry-run:[/italic] skipping write")
        return 0
    write_parquet(df, output_path)
    if name_map and name_map_df is not None:
        out_map = name_map
        if out_map.lower().endswith('.parquet'):
            name_map_df.to_parquet(out_map, index=False, engine='pyarrow', compression='snappy')
        else:
            name_map_df.to_csv(out_map, index=False)
        console.print(f"[green]Wrote name->ST mapping:[/green] {out_map} ({len(name_map_df)} rows)")
    return 0


# typing imports are at top

def main(argv: Optional[List[str]] = None) -> int:
    """Backward-compatible programmatic entrypoint.

    When called with an ``argv`` list (as the tests in this repository do),
    this function will behave like the previous argparse-based main and
    return an exit code. When invoked with ``argv`` is ``None`` from the
    command-line, the Typer CLI will be used instead (see bottom of file).
    """
    parser = argparse.ArgumentParser(description='Convert allele profile to standardized parquet')
    parser.add_argument('input', help='Input profile file (csv, tsv, .gz, metadata)')
    parser.add_argument('output', help='Output parquet file (.parquet)')
    parser.add_argument('--name-map', help='Optional output path to write Name->ST mapping (parquet or csv)')
    parser.add_argument('--dry-run', action='store_true', help='Do not write output; print inferred columns and row counts')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose logging')
    args = parser.parse_args(argv)
    return convert_file(args.input, args.output, name_map=args.name_map, dry_run=args.dry_run, verbose=args.verbose)


@app.command()
def cli(
    input_path: str = typer.Argument(..., help='Input profile file (csv, tsv, .gz, metadata)'),
    output: str = typer.Argument(..., help='Output parquet file (.parquet)'),
    name_map: Optional[str] = typer.Option(None, help='Optional output path to write Name->ST mapping (parquet or csv)'),
    dry_run: bool = typer.Option(False, help='Do not write output; print inferred columns and row counts'),
    verbose: bool = typer.Option(False, '--verbose', '-v', help='Verbose logging'),
) -> None:
    """Command-line interface implemented with Typer.

    This function delegates to :func:`convert_file` for the real work and uses
    :mod:`rich` to print colored and user-friendly console messages.
    """
    raise SystemExit(convert_file(input_path, output, name_map=name_map, dry_run=dry_run, verbose=verbose))


if __name__ == '__main__':
    # Use Typer when running as a script so CLI help and completion are
    # available; programmatic callers should still invoke ``main([...])``.
    app()
