from matplotlib import pyplot as plt
import os

def create_plot(data, outfile):
    # Load data
    all = {}
    for s in data['all']:
        with open(s, 'r') as f:
            fd = f.readlines()[2:]

        dd = {}
        for l in fd:
            fields = l.strip().split()
            dd[int(float(fields[0]))] = int(float(fields[1]))

        total = sum(dd.values())
        covrd = total - dd[0]

        samp = os.path.basename(s).replace('_covdist_all_base_table.txt', '')
        all[samp] = 100 * covrd / total

    q40 = {}
    for s in data['q40']:
        with open(s, 'r') as f:
            fd = f.readlines()[2:]

        dd = {}
        for l in fd:
            fields = l.strip().split()
            dd[int(float(fields[0]))] = int(float(fields[1]))

        total = sum(dd.values())
        covrd = total - dd[0]

        samp = os.path.basename(s).replace('_covdist_q40_base_table.txt', '')
        q40[samp] = 100 * covrd / total

    y = [i for i, s in enumerate(list(all.keys()))]
    n = [s for i, s in enumerate(list(all.keys()))]

    # Create plot
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey='row')
    plt.tight_layout()

    ax1.barh(y, all.values(), height=0.7, color='black')
    ax2.barh(y, q40.values(), height=0.7, color='black')

    ax1.set_xlim(0, 110);
    ax2.set_xlim(0, 110);
    ax1.set_xticks(range(0, 110, 20)); ax1.set_xticklabels([str(i) for i in range(0, 110, 20)], fontsize=18)
    ax2.set_xticks(range(0, 110, 20)); ax2.set_xticklabels([str(i) for i in range(0, 110, 20)], fontsize=18)

    ax1.set_ylim(-1, len(n))
    ax1.set_yticks(y); ax1.set_yticklabels(n, fontsize=18)

    plt.suptitle('Percent of Genome Covered', x=0.5, y=1.10, fontsize=24)
    ax1.set_title('All Reads', fontsize=16)
    ax2.set_title('Q40 Reads', fontsize=16)

    ax1.set_xlabel('Percent Covered', fontsize=20)
    ax2.set_xlabel('Percent Covered', fontsize=20)

    plt.savefig(outfile['out'], bbox_inches='tight')
    plt.close('all')

create_plot(snakemake.input, snakemake.output)
