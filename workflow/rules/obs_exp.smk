###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in obs_exp.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
# SAMPLES          - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

rule obs_exp_coverage_genomecov:
    input:
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
    output:
        cov = temp(f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.tmp.bed.gz'),
    params:
        outdir = f'{output_directory}/analysis/obs_exp',
    log:
        f'{output_directory}/logs/obs_exp/{{sample}}.genomecov.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/{{sample}}.genomecov.txt',
    threads: 1
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['bedtools'],
        config['envmodules']['samtools'],
    shell:
        """
        set +o pipefail
        mkdir -p {params.outdir}

        # Find genomic coverage
        samtools view -hb -F 0x4 -q 40 {input.bam} | \
        bedtools genomecov -bg -ibam stdin | \
        gzip -c > {output.cov}
        """

rule obs_exp_coverage_genomecov_cpg:
    input:
        cov = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.tmp.bed.gz',
    output:
        cpg = temp(f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.cpg.tmp.bed.gz'),
    params:
        cpg = f'{config["obs_exp"]["assets"]}/cpg_bismap.bed.gz',
    log:
        f'{output_directory}/logs/obs_exp/{{sample}}.cpg_genomecov.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/{{sample}}.cpg_genomecov.txt',
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

        # Find intersections
        bedtools intersect -a {params.cpg} -b {input.cov} -wo | gzip -c > {output.cpg}
        """

rule obs_exp_coverage_genomecov_cgi:
    input:
        cov = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.tmp.bed.gz',
    output:
        cgi = temp(f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.cgi.tmp.bed.gz'),
    params:
        cgi = f'{config["obs_exp"]["assets"]}/cgi_bismap.bed.gz',
    log:
        f'{output_directory}/logs/obs_exp/{{sample}}.cgi_genomecov.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/{{sample}}.cgi_genomecov.txt',
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

        # Find intersections
        bedtools intersect -a {params.cgi} -b {input.cov} -wo | gzip -c > {output.cgi}
        """

rule obs_exp_coverage_genomecov_intergenic:
    input:
        cov = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.tmp.bed.gz',
    output:
        intergenic = temp(f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.intergenic.tmp.bed.gz'),
    params:
        intergenic = f'{config["obs_exp"]["assets"]}/intergenic_bismap.bed.gz',
    log:
        f'{output_directory}/logs/obs_exp/{{sample}}.intergenic_genomecov.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/{{sample}}.intergenic_genomecov.txt',
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

        # Find intersections
        bedtools intersect -a {params.intergenic} -b {input.cov} -wo | gzip -c > {output.intergenic}
        """

rule obs_exp_coverage_genomecov_exon:
    input:
        cov = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.tmp.bed.gz',
    output:
        exon = temp(f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.exon.tmp.bed.gz'),
    params:
        exon = f'{config["obs_exp"]["assets"]}/exon_bismap.bed.gz',
    log:
        f'{output_directory}/logs/obs_exp/{{sample}}.exon_genomecov.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/{{sample}}.exon_genomecov.txt',
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

        # Find intersections
        bedtools intersect -a {params.exon} -b {input.cov} -wo | gzip -c > {output.exon}
        """

rule obs_exp_coverage_genomecov_genic:
    input:
        cov = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.tmp.bed.gz',
    output:
        genic = temp(f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.genic.tmp.bed.gz'),
    params:
        genic = f'{config["obs_exp"]["assets"]}/genic_bismap.bed.gz',
    log:
        f'{output_directory}/logs/obs_exp/{{sample}}.genic_genomecov.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/{{sample}}.genic_genomecov.txt',
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

        # Find intersections
        bedtools intersect -a {params.genic} -b {input.cov} -wo | gzip -c > {output.genic}
        """

rule obs_exp_coverage_genomecov_rmsk:
    input:
        cov = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.tmp.bed.gz',
    output:
        rmsk = temp(f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.rmsk.tmp.bed.gz'),
    params:
        rmsk = f'{config["obs_exp"]["assets"]}/rmsk_bismap.bed.gz',
    log:
        f'{output_directory}/logs/obs_exp/{{sample}}.rmsk_genomecov.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/{{sample}}.rmsk_genomecov.txt',
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

        # Find intersections
        bedtools intersect -a {params.rmsk} -b {input.cov} -wo | gzip -c > {output.rmsk}
        """

rule obs_exp_coverage_genomecov_mapped:
    input:
        cov = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.tmp.bed.gz',
    output:
        mapped = temp(f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.mapped.tmp.bed.gz'),
    params:
        mapped = f'{config["obs_exp"]["assets"]}/k100.bismap.bedgraph.gz',
    log:
        f'{output_directory}/logs/obs_exp/{{sample}}.mapped_genomecov.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/{{sample}}.mapped_genomecov.txt',
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

        # Find intersections
        bedtools intersect -a {params.mapped} -b {input.cov} -wo | gzip -c > {output.mapped}
        """

rule obs_exp_coverage_values:
    input:
        cpg = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.cpg.tmp.bed.gz',
        cgi = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.cgi.tmp.bed.gz',
        exon = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.exon.tmp.bed.gz',
        genic = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.genic.tmp.bed.gz',
        intergenic = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.intergenic.tmp.bed.gz',
        rmsk = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.rmsk.tmp.bed.gz',
        mapped = f'{output_directory}/analysis/obs_exp/{{sample}}.genomecov.mapped.tmp.bed.gz',
    output:
        data = f'{output_directory}/analysis/obs_exp/{{sample}}.feature_sizes.txt',
    params:
        cpg = f'{config["obs_exp"]["assets"]}/cpg_bismap.bed.gz',
        cgi = f'{config["obs_exp"]["assets"]}/cgi_bismap.bed.gz',
        exon = f'{config["obs_exp"]["assets"]}/exon_bismap.bed.gz',
        genic = f'{config["obs_exp"]["assets"]}/genic_bismap.bed.gz',
        intergenic = f'{config["obs_exp"]["assets"]}/intergenic_bismap.bed.gz',
        rmsk = f'{config["obs_exp"]["assets"]}/rmsk_bismap.bed.gz',
        bismap = f'{config["obs_exp"]["assets"]}/k100.bismap.bedgraph.gz',
    log:
        f'{output_directory}/logs/obs_exp/{{sample}}.feature_sizes.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/{{sample}}.feature_sizes.txt',
    threads: 1
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/python_packages.yaml'
    script:
        '../scripts/calculate_feature_sizes.py'

def get_all_samples(wildcards):
    files = []
    for s in list(SAMPLES['sample']):
        files.append(f'{output_directory}/analysis/obs_exp/{s}.feature_sizes.txt')
    files.sort()
    return files

rule obs_exp_coverage_plot:
    input:
        files = get_all_samples,
    output:
        cpgs = f'{output_directory}/analysis/obs_exp/obs_exp_cpgs.pdf',
        cgis = f'{output_directory}/analysis/obs_exp/obs_exp_cgis.pdf',
        exon = f'{output_directory}/analysis/obs_exp/obs_exp_exon.pdf',
        gene = f'{output_directory}/analysis/obs_exp/obs_exp_gene.pdf',
        intr = f'{output_directory}/analysis/obs_exp/obs_exp_intr.pdf',
        rmsk = f'{output_directory}/analysis/obs_exp/obs_exp_rmsk.pdf',
    params:
        outdir = f'{output_directory}/analysis/obs_exp',
    log:
        f'{output_directory}/logs/obs_exp/plot_obs_exp.log',
    benchmark:
        f'{output_directory}/benchmarks/obs_exp/plot_obs_exp.txt',
    threads: 1
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/python_packages.yaml'
    script:
        '../scripts/plot_obs_exp.py'
