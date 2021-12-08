###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in obs_exp.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
# SAMPLES          - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

rule covg_uniformity:
    input:
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
    output:
        uni = f'{output_directory}/analysis/covg_uniformity/{{sample}}.10kb_binned_coverage.bed.gz',
    params:
        outdir = f'{output_directory}/analysis/covg_uniformity',
        bins = f'{config["covg_uniformity"]["assets"]}/genome_window_10kb.bed.gz',
    log:
        f'{output_directory}/logs/covg_uniformity/{{sample}}.uniformity.log',
    benchmark:
        f'{output_directory}/benchmarks/covg_uniformity/{{sample}}.uniformity.txt',
    threads: 1
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['bedtools'],
    shell:
        """
        set +o pipefail
        mkdir -p {params.outdir}

        bedtools genomecov -bga -split -ibam {input.bam} | bedtools map -a {params.bins} -b - -c 4 | gzip > {output.uni}
        """

def get_covg_samples(wildcards):
    files = []
    for s in list(SAMPLES['sample']):
        files.append(f'{output_directory}/analysis/covg_uniformity/{s}.10kb_binned_coverage.bed.gz')
    files.sort()
    return files

rule uniformity_plot:
    input:
        files = get_covg_samples,
    output:
        covg = f'{output_directory}/analysis/covg_uniformity/binned_genome_coverage.pdf',
        frac = f'{output_directory}/analysis/covg_uniformity/coverage_uniformity.pdf',
        data = f'{output_directory}/analysis/covg_uniformity/sample_data.tsv',
    params:
        map_scores = f'{config["covg_uniformity"]["assets"]}/k100.bismap.10kb_avg.bed.gz',
    log:
        f'{output_directory}/logs/covg_uniformity/uniformity_plot.log',
    benchmark:
        f'{output_directory}/benchmarks/covg_uniformity/uniformity_plot.txt',
    threads: 1
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/python_packages.yaml'
    script:
        '../scripts/plot_covg_uniformity.py'
