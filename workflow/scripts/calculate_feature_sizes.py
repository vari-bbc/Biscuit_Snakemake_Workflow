import pandas as pd

DEFAULT_2 = {'cov': [0], 'ol': [0]}
DEFAULT_3 = {'start': [0], 'end': [0], 'weight': [0]}

NAMES_2 = ['cov', 'ol']
NAMES_3 = ['start', 'end', 'weight']

# TODO: This can probably be adjusted to read files in chunks to reduce memory overhead
def calculate_feature_sizes(samp, infiles, paramfiles, outfile):
    # Load the asset files
    try:
        df_map_p = pd.read_csv(paramfiles['bismap'], sep='\t', header=None, usecols=[1, 2, 3], skiprows=1, names=NAMES_3)
    except pandas.errors.EmptyDataError:
        df_map_p = pd.DataFrame(DEFAULT_3)
    try:
        df_cpg_p = pd.read_csv(paramfiles['cpg'], sep='\t', header=None, usecols=[1, 2, 3], skiprows=1, names=NAMES_3)
    except pandas.errors.EmptyDataError:
        df_cpg_p = pd.DataFrame(DEFAULT_3)
    try:
        df_cgi_p = pd.read_csv(paramfiles['cgi'], sep='\t', header=None, usecols=[1, 2, 3], skiprows=1, names=NAMES_3)
    except pandas.errors.EmptyDataError:
        df_cgi_p = pd.DataFrame(DEFAULT_3)
    try:
        df_exn_p = pd.read_csv(paramfiles['exon'], sep='\t', header=None, usecols=[1, 2, 3], skiprows=1, names=NAMES_3)
    except pandas.errors.EmptyDataError:
        df_exn_p = pd.DataFrame(DEFAULT_3)
    try:
        df_gen_p = pd.read_csv(paramfiles['genic'], sep='\t', header=None, usecols=[1, 2, 3], skiprows=1, names=NAMES_3)
    except pandas.errors.EmptyDataError:
        df_gen_p = pd.DataFrame(DEFAULT_3)
    try:
        df_int_p = pd.read_csv(paramfiles['intergenic'], sep='\t', header=None, usecols=[1, 2, 3], skiprows=1, names=NAMES_3)
    except pandas.errors.EmptyDataError:
        df_int_p = pd.DataFrame(DEFAULT_3)
    try:
        df_msk_p = pd.read_csv(paramfiles['rmsk'], sep='\t', header=None, usecols=[1, 2, 3], skiprows=1, names=NAMES_3)
    except pandas.errors.EmptyDataError:
        df_msk_p = pd.DataFrame(DEFAULT_3)

    # Load the Snakemake generated files
    try:
        df_map_i = pd.read_csv(infiles['mapped'], sep='\t', header=None, usecols=[3, 8], names=NAMES_2)
    except pandas.errors.EmptyDataError:
        df_map_i = pd.DataFrame(DEFAULT_2)
    try:
        df_cpg_i = pd.read_csv(infiles['cpg'], sep='\t', header=None, usecols=[3, 8], names=NAMES_2)
    except pandas.errors.EmptyDataError:
        df_cpg_i = pd.DataFrame(DEFAULT_3)
    try:
        df_cgi_i = pd.read_csv(infiles['cgi'], sep='\t', header=None, usecols=[3, 8], names=NAMES_2)
    except pandas.errors.EmptyDataError:
        df_cgi_i = pd.DataFrame(DEFAULT_3)
    try:
        df_exn_i = pd.read_csv(infiles['exon'], sep='\t', header=None, usecols=[3, 8], names=NAMES_2)
    except pandas.errors.EmptyDataError:
        df_exn_i = pd.DataFrame(DEFAULT_3)
    try:
        df_gen_i = pd.read_csv(infiles['genic'], sep='\t', header=None, usecols=[3, 8], names=NAMES_2)
    except pandas.errors.EmptyDataError:
        df_gen_i = pd.DataFrame(DEFAULT_3)
    try:
        df_int_i = pd.read_csv(infiles['intergenic'], sep='\t', header=None, usecols=[3, 8], names=NAMES_2)
    except pandas.errors.EmptyDataError:
        df_int_i = pd.DataFrame(DEFAULT_3)
    try:
        df_msk_i = pd.read_csv(infiles['rmsk'], sep='\t', header=None, usecols=[3, 8], names=NAMES_2)
    except pandas.errors.EmptyDataError:
        df_msk_i = pd.DataFrame(DEFAULT_3)

    # Calculate sizes for finding the sum (weight * (end - start))
    df_map_p['size'] = df_map_p['weight'] * (df_map_p['end'] - df_map_p['start'])
    df_cpg_p['size'] = df_cpg_p['weight'] * (df_cpg_p['end'] - df_cpg_p['start'])
    df_cgi_p['size'] = df_cgi_p['weight'] * (df_cgi_p['end'] - df_cgi_p['start'])
    df_exn_p['size'] = df_exn_p['weight'] * (df_exn_p['end'] - df_exn_p['start'])
    df_gen_p['size'] = df_gen_p['weight'] * (df_gen_p['end'] - df_gen_p['start'])
    df_int_p['size'] = df_int_p['weight'] * (df_int_p['end'] - df_int_p['start'])
    df_msk_p['size'] = df_msk_p['weight'] * (df_msk_p['end'] - df_msk_p['start'])

    # Calculate sizes for finding the sum (coverage * number of overlapping bases)
    df_map_i['size'] = df_map_i['cov'] * df_map_i['ol']
    df_cpg_i['size'] = df_cpg_i['cov'] * df_cpg_i['ol']
    df_cgi_i['size'] = df_cgi_i['cov'] * df_cgi_i['ol']
    df_exn_i['size'] = df_exn_i['cov'] * df_exn_i['ol']
    df_gen_i['size'] = df_gen_i['cov'] * df_gen_i['ol']
    df_int_i['size'] = df_int_i['cov'] * df_int_i['ol']
    df_msk_i['size'] = df_msk_i['cov'] * df_msk_i['ol']

    # Calculate the necessary sums
    sum_map_p = df_map_p['size'].sum(); sum_map_i = df_map_i['size'].sum()
    sum_cpg_p = df_cpg_p['size'].sum(); sum_cpg_i = df_cpg_i['size'].sum()
    sum_cgi_p = df_cgi_p['size'].sum(); sum_cgi_i = df_cgi_i['size'].sum()
    sum_exn_p = df_exn_p['size'].sum(); sum_exn_i = df_exn_i['size'].sum()
    sum_gen_p = df_gen_p['size'].sum(); sum_gen_i = df_gen_i['size'].sum()
    sum_int_p = df_int_p['size'].sum(); sum_int_i = df_int_i['size'].sum()
    sum_msk_p = df_msk_p['size'].sum(); sum_msk_i = df_msk_i['size'].sum()

    with open(outfile, 'w') as f:
        f.write(
            '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                samp,
                sum_map_p, sum_cpg_p, sum_cgi_p, sum_msk_p, sum_exn_p, sum_gen_p, sum_int_p,
                sum_map_i, sum_cpg_i, sum_cgi_i, sum_msk_i, sum_exn_i, sum_gen_i, sum_int_i,
            )
        )

calculate_feature_sizes(
    snakemake.wildcards['sample'],
    snakemake.input,
    snakemake.params,
    snakemake.output['data']
)
