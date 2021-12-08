from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os

# Create list of chromosomes (include enough to cover human and mouse)
CHROMS = [f'chr{i}' for i in range(1, 23)]
CHROMS.extend(['chrX', 'chrY', 'chrM'])

def find_x_ticks_labels(chrs):
    """Create list of chromosomes that are included in BED file

    Inputs -
        chrs - chr column from BED file loaded as a DataFrame
    Returns -
        dict, keys are the chromosomes and values are the first index location of the chromosome
    """
    out = {}
    for c in CHROMS:
        try:
            out[c] = list(chrs).index(c)
        except ValueError: # chromosome isn't in list, so ignore
            continue

    return out

def create_plot(files, params, outfile):
    # Setup the order that chromosomes need to be sorted into
    chr_sort_order = pd.api.types.CategoricalDtype(CHROMS, ordered=True)

    # Average bismap mappability scores for each bin
    # Treat bins with no score as having a weight of 0 since there's no data there
    weights = pd.read_csv(
        params['map_scores'],
        sep='\t',
        header=None,
        names=['chr', 'start', 'end', 'raw_weight'],
        na_values='.'
    )
    weights['weight'] = weights['raw_weight'].fillna(0)
    weights['ideal'] = weights['weight'] / (weights['end'] - weights['start'])

    # Setup figure for coverage across genome
    fig, ax = plt.subplots(figsize=(9,5))
    plt.tight_layout()

    plt.title('Weighted Coverage in 10kb Windows', fontsize=24)
    ax.set_xlabel('')
    ax.set_ylabel('Weighted Coverage / Window Width', fontsize=18)

    # Process data and fill coverage across genome figure
    x_values = []
    x_tk_lab = []
    plot_data = {
        'samp': [],
        'frac': [],
        'mean': [],
        'stdv': [],
        'lavg': [],
        'lstd': [],
    }

    for s in files:
        samp = os.path.basename(s).replace('.10kb_binned_coverage.bed.gz', '')

        # Read data
        df = pd.read_csv(s, sep='\t', header=None, names=['chr', 'start', 'end', 'covg'])

        # Merge average bismap mappability scores as a weight column
        df = df.merge(weights, on=['chr', 'start', 'end'])

        # Sort chromosomes to put chrM at the end
        df['chr'] = df['chr'].astype(chr_sort_order)
        df = df.sort_values(['chr', 'start'], ignore_index=True)

        # Calculate the weighted coverage [(avg mappability) * (coverage) / (window width)]
        df['weighted_covg'] = df['weight'] * df['covg'] / (df['end'] - df['start'])

        # Bins with both a raw weight and at least one base covered
        to_plot = df[(df['raw_weight'].notna()) & (df['covg'] > 0)]

        # Determine number of bins that had a bismap weight and at least one base covered
        n_weights = len(df[df['raw_weight'].notna()]['chr'])
        n_wgt_cov = len(to_plot[to_plot['weighted_covg'] > 0.01]['chr'])

        # Non-zero bins
        non_zero = [i if i > 0 else None for i in list(to_plot['weighted_covg'])]

        # Setup data that will be saved to output file and plotted in other figures
        plot_data['samp'].append(samp)
        plot_data['frac'].append(n_wgt_cov / n_weights)
        plot_data['mean'].append(np.average([i for i in non_zero if i]))
        plot_data['stdv'].append(np.std([i for i in non_zero if i]))
        plot_data['lavg'].append(np.average([np.log10(i) for i in non_zero if i]))
        plot_data['lstd'].append(np.std([np.log10(i) for i in non_zero if i]))

        # Fill x values and ticks/labels if that hasn't been done yet
        if len(x_values) == 0 or len(x_tk_lab) == 0:
            x_values = list(df.index)
            x_tk_lab = find_x_ticks_labels(df['chr'])

        # Add weighted coverage to coverage across genome figure
        ax.plot(x_values, df['weighted_covg'], 'k-')

    # Finish up coverage across genome figure
    ax.set_xlim(-500, len(x_values)+500)
    ax.set_xticks(list(x_tk_lab.values()))
    ax.set_xticklabels(x_tk_lab.keys(), rotation=90)

    plt.savefig(outfile['covg'], bbox_inches='tight')
    plt.close('all')

    # Create other figures
    df = pd.DataFrame(plot_data)

    fig, ax = plt.subplots()
    plt.tight_layout()

    sns.scatterplot(ax=ax, x='frac', y='lstd', data=df)

    line_loc = np.std([np.log10(i) for i in weights['ideal'] if i > 0])
    ax.axhline(line_loc, alpha=0.6, color='grey', linestyle='--')
    ax.text(x = 0.01, y = line_loc-0.01, s = 'Ideal Std. Dev.', ha='left', va='top', fontsize=14)

    ax.set_xlim(0, 1.1*max(plot_data['frac'])); ax.set_ylim(0, 1.1*max(plot_data['lstd']))

    plt.title('Whole Genome Coverage Uniformity', fontsize=24)
    plt.xlabel('(# Bins with Weighted Coverage > 0.01) /\n(# Bins with Nonzero Weights)', fontsize=18, va='top')
    plt.ylabel('std[ log10(weighted cov.) ]', fontsize=18)

    plt.savefig(outfile['frac'], bbox_inches='tight')
    plt.close('all')

    # Write data to output file
    df.to_csv(outfile['data'], index=False, sep='\t')

create_plot(list(snakemake.input['files']), snakemake.params, snakemake.output)
