from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os

TITLE = {
    'AllCpGs': 'All CpGs in Genome',
    'CGICpGs': 'CpGs in CpG Islands',
    'ExonicCpGs': 'CpGs in Exons',
    'GenicCpGs': 'CpGs in Genes',
    'RepeatCpGs': 'CpGs in Repeat-masked Space'
}

def create_plot(samples, outfiles):
    # Load the data and rearrange by feature type
    d_df = {}
    for s in samples:
        df = pd.read_csv(s, sep='\t').set_index(['Feature'], drop=False)
        all_cnt = df['CpG_Count']['AllCpGs']
        all_q40 = df['Q40_Reads']['AllCpGs']
        all_all = df['All_Reads']['AllCpGs']

        for row in df.itertuples():
            data = {'sample': [], 'l2_q40': [], 'l2_all': [], 'q40_percent': [], 'all_percent': []}
            data['sample'].append(os.path.basename(s).replace('.cpg_stats_feature_table.tsv', ''))
            data['l2_q40'].append(np.log2( (row.Q40_Reads / all_q40) / (row.CpG_Count / all_cnt) ))
            data['l2_all'].append(np.log2( (row.All_Reads / all_all) / (row.CpG_Count / all_cnt) ))
            data['q40_percent'].append(100.0 * row.Q40_Reads / row.CpG_Count)
            data['all_percent'].append(100.0 * row.All_Reads / row.CpG_Count)

            if row.Feature not in d_df.keys():
                d_df[row.Feature] = pd.DataFrame(data)
            else:
                d_df[row.Feature] = d_df[row.Feature].append(pd.DataFrame(data), sort=True)

    for key, df in d_df.items():
        y = [i for i, s in enumerate(list(df['sample']))]
        n = [s for i, s in enumerate(list(df['sample']))]

        fig, ((ax1,ax2), (ax3,ax4)) = plt.subplots(nrows=2, ncols=2, sharex='col', sharey='row')
        plt.tight_layout()

        # Plots
        ax1.barh(y, df['all_percent'], height=0.7, color='black')
        ax2.plot(df['l2_all'], y, 'kx', markersize=8)
        ax2.axvline(0.0, alpha=0.6, color='grey', linestyle='--')

        ax3.barh(y, df['q40_percent'], height=0.7, color='black')
        ax4.plot(df['l2_q40'], y, 'kx', markersize=8)
        ax4.axvline(0.0, alpha=0.6, color='grey', linestyle='--')

        # Axis ticks and labels
        ax1.set_xlim(0, 110)
        ax2.set_xlim(-3.5, 5.5)

        ax1.set_xticks(range(0, 110, 10)); ax1.set_xticklabels([str(i) for i in range(0, 110, 10)], fontsize=18)
        ax2.set_xticks(range(-3, 6))     ; ax2.set_xticklabels([str(i) for i in range(-3, 6)], fontsize=18)

        ax1.set_ylim(-1, len(n))
        ax3.set_ylim(-1, len(n))

        ax1.set_yticks(y); ax1.set_yticklabels(n, fontsize=18)
        ax3.set_yticks(y); ax3.set_yticklabels(n, fontsize=18)

        # Titles and such
        plt.suptitle(TITLE[key], x=0.5, y=1.05, fontsize=24)
        ax1.set_title('', fontsize=20)
        ax2.set_title('', fontsize=20)
        ax3.set_title('', fontsize=20)
        ax4.set_title('', fontsize=20)

        ax1.set_xlabel('', fontsize=20)
        ax2.set_xlabel('', fontsize=20)
        ax3.set_xlabel('Percent Covered', fontsize=20)
        ax4.set_xlabel('log2(obs / exp)', fontsize=20)

        ax1.set_ylabel('All Reads', fontsize=20)
        ax2.set_ylabel('')
        ax3.set_ylabel('Q40 Reads', fontsize=20)
        ax4.set_ylabel('')

        plt.savefig(outfiles[key], bbox_inches='tight')
        plt.close('all')

create_plot(snakemake.input['files'], snakemake.output)
