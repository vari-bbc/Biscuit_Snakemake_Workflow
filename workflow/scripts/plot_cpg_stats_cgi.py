from matplotlib import pyplot as plt
import pandas as pd
import os

def create_plot(samples, outfiles):
    # Load data
    l_df = []
    for s in samples:
        df = pd.read_csv(s, sep='\t', header=None, names=['group', 'count']).set_index(['group'])

        d = {
            'sample': [os.path.basename(s).replace('.cpg_stats_cgi_table.txt', '')],
            '1': [100.0 * df['count']['one_cpg'] / df['count']['n_cpg_islands']],
            '3': [100.0 * df['count']['three_cpgs'] / df['count']['n_cpg_islands']],
            '5': [100.0 * df['count']['five_cpgs'] / df['count']['n_cpg_islands']],
            '7': [100.0 * df['count']['seven_cpgs'] / df['count']['n_cpg_islands']],
            'n_cgis': [df['count']['n_cpg_islands']]
        }

        l_df.append(pd.DataFrame(d))

    df = pd.concat(l_df)
    df = df.sort_values(by = 'sample')

    y = [i for i, s in enumerate(list(df['sample']))]
    n = [s for i, s in enumerate(list(df['sample']))]

    # Create plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex='col', sharey='row')
    plt.tight_layout()

    ax1.barh(y, df['1'], height=0.7, color='black')
    ax2.barh(y, df['3'], height=0.7, color='black')
    ax3.barh(y, df['5'], height=0.7, color='black')
    ax4.barh(y, df['7'], height=0.7, color='black')

    ax1.set_xlim(0, 110); ax1.set_ylim(-1, len(n)+1)
    ax2.set_xlim(0, 110); ax3.set_ylim(-1, len(n)+1)
    ax1.set_xticks(range(0, 110, 10)); ax1.set_xticklabels([str(i) for i in range(0, 110, 10)], fontsize=18)
    ax2.set_xticks(range(0, 110, 10)); ax2.set_xticklabels([str(i) for i in range(0, 110, 10)], fontsize=18)

    ax1.set_yticks(y); ax1.set_yticklabels(n, fontsize=18)
    ax3.set_yticks(y); ax3.set_yticklabels(n, fontsize=18)

    ax1.text(50, len(n), '1 Read', va='center', ha='center', size=16)
    ax2.text(50, len(n), '3 Reads', va='center', ha='center', size=16)
    ax3.text(50, len(n), '5 Reads', va='center', ha='center', size=16)
    ax4.text(50, len(n), '7 Reads', va='center', ha='center', size=16)

    plt.suptitle('# Reads Covering CpG Island', x=0.5, y=1.05, fontsize=24)
    ax1.set_title('', fontsize=16)
    ax2.set_title('', fontsize=16)
    ax3.set_title('', fontsize=16)
    ax4.set_title('', fontsize=16)

    ax1.set_xlabel('', fontsize=20)
    ax2.set_xlabel('', fontsize=20)
    ax3.set_xlabel('Percent Covered', fontsize=20)
    ax4.set_xlabel('Percent Covered', fontsize=20)

    ax1.set_ylabel('')
    ax2.set_ylabel('')
    ax3.set_ylabel('')
    ax4.set_ylabel('')

    plt.savefig(outfiles['cgi'], bbox_inches='tight')
    plt.close('all')

create_plot(snakemake.input['files'], snakemake.output)
