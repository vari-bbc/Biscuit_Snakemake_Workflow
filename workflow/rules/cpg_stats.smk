###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in obs_exp.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

rule cpg_stats_genomecov:
    input:
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
    output:
        unf = f'{output_directory}/analysis/cpg_stats/{{sample}}.genomecov.cpg.q00.bed.gz',
        fil = f'{output_directory}/analysis/cpg_stats/{{sample}}.genomecov.cpg.q40.bed.gz',
    params:
        outdir = f'{output_directory}/analysis/cpg_stats',
        cpg = f'{config["cpg_stats"]["assets"]}/cpg.bed.gz',
    log:
        f'{output_directory}/logs/cpg_stats/{{sample}}.genomecov.cpg.log',
    benchmark:
        f'{output_directory}/benchmarks/cpg_stats/{{sample}}.genomecov.cpg.log',
    threads: 8
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    envmodules:
        config['envmodules']['bedtools'],
        config['envmodules']['samtools'],
    shell:
        """
        set +o pipefail
        mkdir -p {params.outdir}

        # Unfiltered coverage
        bedtools genomecov -bga -split -ibam {input.bam} | \
        LC_ALL=C sort --parallel={threads} -k1,1 -k2,2n -T {params.outdir} | \
        bedtools intersect -a {params.cpg} -b - -wo -sorted | \
        bedtools groupby -g 1-3 -c 7 -o min | \
        gzip -c > {output.unf}

        # MAPQ >= 40 filtered coverage
        samtools view -q 40 -b {input.bam} | \
        bedtools genomecov -bga -split -ibam stdin | \
        LC_ALL=C sort --parallel={threads} -k1,1 -k2,2n -T {params.outdir} | \
        bedtools intersect -a {params.cpg} -b - -wo -sorted | \
        bedtools groupby -g 1-3 -c 7 -o min | \
        gzip -c > {output.fil}
        """

rule cpg_stats_feature_table:
    input:
        unf = f'{output_directory}/analysis/cpg_stats/{{sample}}.genomecov.cpg.q00.bed.gz',
        fil = f'{output_directory}/analysis/cpg_stats/{{sample}}.genomecov.cpg.q40.bed.gz',
    output:
        f'{output_directory}/analysis/cpg_stats/{{sample}}.cpg_stats_feature_table.tsv',
    params:
        cgi = f'{config["cpg_stats"]["assets"]}/cpg_islands.bed.gz',
        exn = f'{config["cpg_stats"]["assets"]}/exon.bed.gz',
        msk = f'{config["cpg_stats"]["assets"]}/rmsk.bed.gz',
        gen = f'{config["cpg_stats"]["assets"]}/genic.bed.gz',
    log:
        f'{output_directory}/logs/cpg_stats/{{sample}}.feature_table.log',
    benchmark:
        f'{output_directory}/benchmarks/cpg_stats/{{sample}}.feature_table.log',
    threads: 8
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    envmodules:
        config['envmodules']['bedtools'],
    shell:
        """
        echo -e "Feature\tCpG_Count\tQ40_Reads\tAll_Reads" > {output}

        # All CpGs
        zcat {input.unf} | wc -l | awk '{{ printf("AllCpGs\\t%s", $1) }}' >> {output}
        zcat {input.fil} | awk '$4>0{{ a += 1 }} END{{ printf("\\t%d", a) }}' >> {output}
        zcat {input.unf} | awk '$4>0{{ a += 1 }} END{{ printf("\\t%d\\n", a) }}' >> {output}

        # CpG Island CpGs
        bedtools intersect -a {input.unf} -b <(bedtools merge -i {params.cgi}) -sorted | wc -l | \
        awk '{{ printf("CGICpGs\\t%s", $1) }}' >> {output}
        bedtools intersect -a {input.fil} -b <(bedtools merge -i {params.cgi}) -sorted | \
        awk '$4>0{{ a += 1 }} END{{ printf("\\t%d", a) }}' >> {output}
        bedtools intersect -a {input.unf} -b <(bedtools merge -i {params.cgi}) -sorted | \
        awk '$4>0{{ a += 1 }} END{{ printf("\\t%d\\n", a) }}' >> {output}

        # Exonic CpGs
        bedtools intersect -a {input.unf} -b <(bedtools merge -i {params.exn}) -sorted | wc -l | \
        awk '{{ printf("ExonicCpGs\\t%s", $1) }}' >> {output}
        bedtools intersect -a {input.fil} -b <(bedtools merge -i {params.exn}) -sorted | \
        awk '$4>0{{ a += 1 }} END{{ printf("\\t%d", a) }}' >> {output}
        bedtools intersect -a {input.unf} -b <(bedtools merge -i {params.exn}) -sorted | \
        awk '$4>0{{ a += 1 }} END{{ printf("\\t%d\\n", a) }}' >> {output}

        # Genic CpGs
        bedtools intersect -a {input.unf} -b <(bedtools merge -i {params.gen}) -sorted | wc -l | \
        awk '{{ printf("GenicCpGs\\t%s", $1) }}' >> {output}
        bedtools intersect -a {input.fil} -b <(bedtools merge -i {params.gen}) -sorted | \
        awk '$4>0{{ a += 1 }} END{{ printf("\\t%d", a) }}' >> {output}
        bedtools intersect -a {input.unf} -b <(bedtools merge -i {params.gen}) -sorted | \
        awk '$4>0{{ a += 1 }} END{{ printf("\\t%d\\n", a) }}' >> {output}

        # Repeat-masked CpGs
        bedtools intersect -a {input.unf} -b <(bedtools merge -i {params.msk}) -sorted | wc -l | \
        awk '{{ printf("RepeatCpGs\\t%s", $1) }}' >> {output}
        bedtools intersect -a {input.fil} -b <(bedtools merge -i {params.msk}) -sorted | \
        awk '$4>0{{ a += 1 }} END{{ printf("\\t%d", a) }}' >> {output}
        bedtools intersect -a {input.unf} -b <(bedtools merge -i {params.msk}) -sorted | \
        awk '$4>0{{ a += 1 }} END{{ printf("\\t%d\\n", a) }}' >> {output}
        """

rule cpg_stats_cgi_table:
    input:
        fil = f'{output_directory}/analysis/cpg_stats/{{sample}}.genomecov.cpg.q40.bed.gz',
    output:
        f'{output_directory}/analysis/cpg_stats/{{sample}}.cpg_stats_cgi_table.txt',
    params:
        cgi = f'{config["cpg_stats"]["assets"]}/cpg_islands.bed.gz',
        tmp = f'{output_directory}/analysis/cpg_stats/{{sample}}.cpg_dist.tmp.txt',
    log:
        f'{output_directory}/logs/cpg_stats/{{sample}}.cgi_table.log',
    benchmark:
        f'{output_directory}/benchmarks/cpg_stats/{{sample}}.cgi_table.log',
    threads: 8
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    envmodules:
        config['envmodules']['bedtools'],
    shell:
        """
        # Total number of CpG islands
        zcat {params.cgi} | wc -l | awk '{{ printf("n_cpg_islands\\t%s\\n", $1) }}' > {output}

        # Number of CpG islands with at least one read with MAPQ>=40 covering 1+, 3+, 5+, and 7+ CpGs in that CpG island
        bedtools intersect -a {input.fil} -b <(bedtools merge -i {params.cgi}) -sorted -wo | \
        awk '$4>0 {{ print $5":"$6"-"$7 }}' | uniq -c | \
        awk '{{ if ($1 >= 1) {{ a += 1; }} if ($1 >= 3) {{ b += 1; }} if ($1 >= 5) {{ c += 1; }} if ($1 >= 7) {{ d += 1; }} }} END {{ printf("%d %d %d %d", a, b, c, d) }}' | \
        awk '{{ printf("one_cpg\\t%s\\nthree_cpgs\\t%s\\nfive_cpgs\\t%s\\nseven_cpgs\\t%s\\n", $1, $2, $3, $4) }}' >> {output}
        """

def get_all_feature_samples(wildcards):
    files = []
    for s in list(SAMPLES['sample']):
        files.append(f'{output_directory}/analysis/cpg_stats/{s}.cpg_stats_feature_table.tsv')
    files.sort()
    return files

rule cpg_stats_feature_plot:
    input:
        files = get_all_feature_samples,
    output:
        AllCpGs = f'{output_directory}/analysis/cpg_stats/feature_all_cpgs.pdf',
        CGICpGs = f'{output_directory}/analysis/cpg_stats/feature_cpg_islands.pdf',
        ExonicCpGs = f'{output_directory}/analysis/cpg_stats/feature_exon_cpgs.pdf',
        GenicCpGs = f'{output_directory}/analysis/cpg_stats/feature_gene_cpgs.pdf',
        RepeatCpGs = f'{output_directory}/analysis/cpg_stats/feature_repeat_cpgs.pdf',
    params:
        outdir = f'{output_directory}/analysis/cpg_stats',
    log:
        f'{output_directory}/logs/cpg_stats/plot_cpg_stats_feature.log',
    benchmark:
        f'{output_directory}/benchmarks/cpg_stats/plot_cpg_stats_feature.log',
    threads: 1
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    script:
        '../scripts/plot_cpg_stats_features.py'

def get_all_cgi_samples(wildcards):
    files = []
    for s in list(SAMPLES['sample']):
        files.append(f'{output_directory}/analysis/cpg_stats/{s}.cpg_stats_cgi_table.txt')
    files.sort()
    return files

rule cpg_stats_cgi_plot:
    input:
        files = get_all_cgi_samples,
    output:
        cgi = f'{output_directory}/analysis/cpg_stats/cgi_stats.pdf',
    params:
        outdir = f'{output_directory}/analysis/cpg_stats',
    log:
        f'{output_directory}/logs/cpg_stats/plot_cpg_stats_cgi.log',
    benchmark:
        f'{output_directory}/benchmarks/cpg_stats/plot_cpg_stats_cgi.log',
    threads: 1
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    script:
        '../scripts/plot_cpg_stats_cgi.py'
