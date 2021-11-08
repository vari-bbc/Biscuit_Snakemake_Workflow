from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os

# Column names for Snakemake generated files
cols = [
    'sample',
    'refln', 'cpgln', 'cgiln', 'mskln', 'exnln', 'genln', 'intln',
    'mapln', 'fcpgs', 'fcgis', 'frmsk', 'fexon', 'fgene', 'fintr'
]

# Make it easier to parse the files based on their column names
ORDER = {
    'cpgs': ('cpgln', 'fcpgs'),
    'cgis': ('cgiln', 'fcgis'),
    'rmsk': ('mskln', 'frmsk'),
    'exon': ('exnln', 'fexon'),
    'gene': ('genln', 'fgene'),
    'intr': ('intln', 'fintr')
}

# Which columns to keep for generating plot
keepers = ['sample'] + list(ORDER.keys())

# Plot titles for each column
TITLE = {
    'cpgs': 'Observed / Expected Coverage for All CpGs',
    'cgis': 'Observed / Expected Coverage for CpG Islands',
    'rmsk': 'Observed / Expected Coverage for Repeat-Masked Space',
    'exon': 'Observed / Expected Coverage for Exons',
    'gene': 'Observed / Expected Coverage for Genic Space',
    'intr': 'Observed / Expected Coverage for Intergenic Space'
}

def create_plot(samples, outfiles):
    # Load the data and calculate the obs/exp ratios
    l_df = []
    for s in samples:
        df = pd.read_csv(s, sep='\t', header=None, names=cols)
        for k, v in ORDER.items():
            num = df[v[0]] / df['refln']
            den = df[v[1]] / df['mapln']

            df[k] = np.log2(num / den)
        l_df.append(df)

    # Combine dataframes for plotting
    df = pd.concat(l_df)
    df = df[keepers]
    df = df.sort_values(by = 'sample')

    y     = [i for i, s in enumerate(list(df['sample']))]
    names = [s for i, s in enumerate(list(df['sample']))]

    # Create a plot for each of the features
    for k in ORDER.keys():
        fig, ax = plt.subplots()
        plt.tight_layout()

        plt.plot(list(df[k]), y, 'kx', markersize=8)
        ax.axvline(0.0, alpha=0.6, color='grey', linestyle='--')

        plt.title(TITLE[k], fontsize=24)
        plt.xlabel('log2(obs / exp)', fontsize=20)
        plt.ylabel('', fontsize=20)

        plt.xlim(-3.5, 5.5)
        plt.ylim(-1, len(names))

        plt.xticks(
            [i for i in np.arange(-3, 6, 1)],
            [str(i) for i in np.arange(-3, 6, 1)],
            fontsize=18
        )
        plt.yticks(y, names, fontsize=18)

        plt.savefig(outfiles[k], bbox_inches='tight')
        plt.close('all')

create_plot(list(snakemake.input['files']), snakemake.output)
