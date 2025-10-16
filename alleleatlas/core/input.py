"""Profile loading and format normalization.

Supports multiple cgMLST profile formats:
- Standard TSV/CSV with ID column and loci columns
- Gzipped/compressed files (.gz, .xz)
- PathogenWatch format (checksum + cgmlst columns)
- EnteroBase format (ST column)
- Raw matrix format (no header, whitespace-separated)
"""

import re
import pandas as pd
from rich.console import Console

console = Console()


def read_file_preview(cgmlst_profiles, nlines=10):
    """
    Reads the first few lines of a tabular file (CSV/TSV, optionally compressed) to detect columns or preview data.

    Parameters:
        cgmlst_profiles (str): Path to the input file. Supports .csv, .tsv, .tsv.gz, .tsv.xz, and .gz/.xz compressed files.
        nlines (int, optional): Number of lines to read from the start of the file. If set to -1, reads the entire file. Defaults to 10.

    Returns:
        pandas.DataFrame: DataFrame containing the previewed data from the file.
    """
    # Try common separators in order: tab, comma, whitespace (regex). Return the
    # first successful parse that yields more than one column or a header with
    # recognizable fields.
    read_base = {'dtype': str}
    if nlines != -1:
        read_base['nrows'] = nlines
    if cgmlst_profiles.endswith('.gz'):
        read_base['compression'] = 'gzip'
    elif cgmlst_profiles.endswith('.xz'):
        read_base['compression'] = 'xz'

    # 1) try tab
    try:
        df = pd.read_csv(cgmlst_profiles, **{**read_base, 'sep': '\t'})
        if df.shape[1] > 1 or 'ST' in df.columns:
            return df
    except Exception:
        pass

    # 2) try comma
    try:
        df = pd.read_csv(cgmlst_profiles, **{**read_base, 'sep': ','})
        if df.shape[1] > 1 or 'ST' in df.columns:
            return df
    except Exception:
        pass

    # 3) try whitespace regex (use python engine)
    try:
        df = pd.read_csv(cgmlst_profiles, **{**read_base, 'sep': r'\s+', 'engine': 'python'})
        return df
    except Exception:
        # last resort: let pandas try default inference
        return pd.read_csv(cgmlst_profiles, **read_base)


def detect_and_normalize_profile(cgmlst_profiles, remap_alleles=False):
    """
    Read the start of the file to detect input format and if it's PathogenWatch-style
    (contains `checksum` and a `cgmlst` column) normalize it into a table with columns:
      ST, L1, L2, L3, ..., <original other columns...>

    Parameters:
        cgmlst_profiles (str): Path to the input profile file
        remap_alleles (bool): If True, convert allele tokens to integers (for PathogenWatch format)

    Returns:
        tuple: (df, st_counts_dict) where:
            - df: Normalized DataFrame with ST and loci columns
            - st_counts_dict: Dictionary mapping ST values to their counts
    """
    st_counts_dict = None
    df = read_file_preview(cgmlst_profiles, nlines=10)
    if df.shape[0] == 0:
        raise ValueError('Empty input file')
    input_type = 'unknown'
    ST_col_name = 'ST'
    # remove legacy Name column if present
    if 'Name' in df.columns:
        input_type = 'regular'
    if 'checksum' in df.columns and ('cgmlst' in df.columns or 'cgMLST' in df.columns):
        input_type = 'pathogenwatch'
        ST_col_name = 'checksum'
    if 'ST' in df.columns and 'Name' not in df.columns:
        input_type = 'enterobase'

    if input_type == 'unknown':
        # Try one more time: some profile files are raw whitespace-separated
        # matrices with no header (rows = samples, columns = loci). Detect
        # that shape by reading a few lines with a whitespace separator.
        try:
            read_kwargs = {'dtype': str, 'sep': r'\s+', 'header': None, 'nrows': 10}
            if cgmlst_profiles.endswith('.gz'):
                read_kwargs['compression'] = 'gzip'
            elif cgmlst_profiles.endswith('.xz'):
                read_kwargs['compression'] = 'xz'
            test_df = pd.read_csv(cgmlst_profiles, **read_kwargs)
            # If we see multiple columns, assume it's a raw matrix file.
            if test_df.shape[1] > 1:
                input_type = 'matrix'
            else:
                raise ValueError('Unable to determine input type; please provide a valid cgMLST profile file.')
        except Exception:
            raise ValueError('Unable to determine input type; please provide a valid cgMLST profile file.')
        
    # At this point we know the input type and can re-read the full file
    if input_type == 'matrix':
        # Read full whitespace-separated matrix without header
        read_kwargs = {'dtype': str, 'sep': r'\s+', 'header': None}
        if cgmlst_profiles.endswith('.gz'):
            read_kwargs['compression'] = 'gzip'
        elif cgmlst_profiles.endswith('.xz'):
            read_kwargs['compression'] = 'xz'
        df_mat = pd.read_csv(cgmlst_profiles, **read_kwargs)
        # Create ST column and loci columns L1..Ln
        ncols = df_mat.shape[1]
        loci_cols = [f'L{i+1}' for i in range(ncols)]
        df_mat.columns = loci_cols
        df_mat.insert(0, 'ST', range(1, len(df_mat) + 1))
        # st_counts: each ST is unique here
        st_counts_dict = {str(i): 1 for i in df_mat['ST'].tolist()}
        df = df_mat
    else:
        df = read_file_preview(cgmlst_profiles, nlines=-1)
    st_counts_dict = df[ST_col_name].value_counts().to_dict()

    if 'Name' in df.columns:
        df = df.drop(columns=['Name'])

    # If ST column exists and has duplicates, keep first and record counts
    if df[ST_col_name].duplicated().any():
        df = df.drop_duplicates(subset=[ST_col_name])
        console.print('Warning: Duplicate STs found in header sample. Kept first occurrence; returning header-level ST counts in dictionary.')

    # PathogenWatch style: has 'checksum' and a cgmlst column (cgmlst / cgMLST)
    if 'checksum' in df.columns and ('cgmlst' in df.columns or 'cgMLST' in df.columns):
        cgcol = 'cgmlst' if 'cgmlst' in df.columns else 'cgMLST'

        # create a sequential ST column starting at 1 (user requested sequential ST numbers)
        df.insert(0, 'ST', range(1, len(df) + 1))

        # split cgmlst tokens into loci columns L1..Ln using '_' (also accept commas)
        toks = df[cgcol].fillna('').astype(str).map(lambda s: s.replace(',', '_').split('_'))

        # Build loci DataFrame in one step to avoid repeated column inserts (fragmentation)
        toks_list = toks.tolist()
        if len(toks_list) == 0:
            loci_df = pd.DataFrame()
        else:
            loci_df = pd.DataFrame(toks_list).fillna('').astype(str)

        # rename loci columns to L1..Ln
        max_len = loci_df.shape[1]
        loci_df.columns = [f'L{i+1}' for i in range(max_len)]

        # place ST and loci columns first
        df = pd.concat([df[['ST']].reset_index(drop=True), loci_df.reset_index(drop=True)], axis=1)
        loci_cols = list(loci_df.columns)

        # There are values that are MD5-like hex strings; convert those to stable integers
        hex_re = re.compile(r'^[0-9a-f]{8,}$', re.IGNORECASE)

        def _token_to_int(tok):
            """Convert a locus token to an integer when appropriate.

            - Numeric strings -> int
            - Long hex-like strings (MD5) -> reduced positive int
            - Otherwise -> stable hashed int
            """
            if tok is None:
                return ''
            s = str(tok).strip()
            if s == '' or s == '_' or s.lower() == 'nan':
                return ''
            # pure decimal
            if s.isdigit():
                try:
                    return int(s)
                except Exception:
                    pass
            # hex-like (md5 or long hex): convert leading portion to int
            if hex_re.match(s):
                try:
                    # take up to 16 hex chars (64 bits) then reduce to 31-bit positive
                    return int(s[:16], 16) % 2147483647
                except Exception:
                    pass
            # fallback: stable python hash reduced to positive int
            try:
                return abs(hash(s)) % 2147483647
            except Exception:
                return s

        # Apply token->int conversion across loci DataFrame in a vectorized manner
        if loci_cols and remap_alleles:
            # apply token->int conversion column-wise using Series.map to avoid deprecated applymap
            converted_loci = df[loci_cols].apply(lambda s: s.map(_token_to_int))
            # replace loci columns with converted integers (single assignment)
            df[loci_cols] = converted_loci

    return df, st_counts_dict


__all__ = ['read_file_preview', 'detect_and_normalize_profile']
