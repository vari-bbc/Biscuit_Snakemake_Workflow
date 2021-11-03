###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in qc.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
# SAMPLES          - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

rule biscuit_qc:
    input:
        ref = get_biscuit_reference,
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
        vcf = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz',
    output:
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_base_botgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_base_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_base_topgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_cpg_botgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_cpg_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_all_cpg_topgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_base_botgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_base_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_base_topgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_cpg_botgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_cpg_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_covdist_q40_cpg_topgc_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_cv_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_mapq_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_strand_table.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_dup_report.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_totalBaseConversionRate.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_totalReadConversionRate.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_CpHRetentionByReadPos.txt',
        f'{output_directory}/analysis/BISCUITqc/{{sample}}_CpGRetentionByReadPos.txt',
    params:
        assets = config['ref']['assets'],
        output_dir = f'{output_directory}/analysis/BISCUITqc',
    log:
        f'{output_directory}/logs/biscuit_qc/{{sample}}_QC.log'
    benchmark:
        f'{output_directory}/benchmarks/biscuit_qc/{{sample}}.txt'
    threads: 8
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['bedtools'],
        config['envmodules']['htslib'],
        config['envmodules']['samtools'],
        config['envmodules']['parallel'],
    shell:
        """
        set +o pipefail;
        QC.sh \
            -o {params.output_dir} \
            --vcf {input.vcf} \
            {params.assets} \
            {input.ref} \
            {wildcards.sample} \
            {input.bam} \
            2> {log}
        """

rule preseq:
    input:
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
    output:
        tsv = f'{output_directory}/analysis/preseq/{{sample}}.ccurve.txt',
    params:
        dir = f'{output_directory}/analysis/preseq',
        out = f'{output_directory}/analysis/preseq/{{sample}}.ccurve.txt',
    log:
        f'{output_directory}/logs/preseq/{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/preseq/{{sample}}.txt',
    threads: 1,
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    envmodules:
        config['envmodules']['preseq'],
    shell:
        """
        mkdir -p {params.dir}
        preseq c_curve -o {params.out} -P {input.bam} 2> {log}
        """

def get_multiQC_params(wildcards):
    raw = config['fastqs']
    out = output_directory
    indirs = f'{raw} {out}/analysis/BISCUITqc {out}/analysis/raw_fastqc {out}/analysis/align'
    if config['run_fastq_screen']:
        indirs += f' {out}/analysis/fastq_screen' # space needed at beginning to separate directories
    if config['trim_galore']['trim_before_BISCUIT']:
        indirs += f' {out}/analysis/trim_reads' # space needed at beginning to separate directories
    if config['preseq']:
        indirs += f' {out}/analysis/preseq' # space needed at beginning to separate directories
    return indirs
        
rule multiQC:
    input:
        # raw fastqc
        expand(f'{output_directory}/analysis/raw_fastqc/{{samples.sample}}-1-R{{read}}_fastqc.zip', read=[1,2], samples=SAMPLES.itertuples()),
        # flagstat
        expand(f'{output_directory}/analysis/align/{{samples.sample}}.sorted.markdup.bam.flagstat', samples=SAMPLES.itertuples()),
        # biscuit_qc
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_base_botgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_base_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_base_topgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_cpg_botgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_cpg_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_all_cpg_topgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_base_botgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_base_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_base_topgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_cpg_botgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_cpg_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_covdist_q40_cpg_topgc_table.txt', samples=SAMPLES.itertuples()),
        expand(f'{output_directory}/analysis/BISCUITqc/{{samples.sample}}_cv_table.txt', samples=SAMPLES.itertuples()),
        # fastq_screen
        expand(f'{output_directory}/analysis/fastq_screen/{{samples.sample}}-1-R{{read}}_screen.txt', read=[1,2], samples=SAMPLES.itertuples()) if config['run_fastq_screen'] else [],
        # trim_galore
        expand(f'{output_directory}/analysis/trim_reads/{{samples.sample}}-R{{read}}_val_{{read}}_merged.fq.gz', read=[1,2], samples=SAMPLES.itertuples()) if config['trim_galore']['trim_before_BISCUIT'] else [],
        # preseq
        expand(f'{output_directory}/analysis/preseq/{{samples.sample}}.ccurve.txt', samples=SAMPLES.itertuples()) if config['preseq'] else [],
    output:
        directory(f'{output_directory}/analysis/multiqc/multiqc_report_data',),
        f'{output_directory}/analysis/multiqc/multiqc_report.html',
    params:
        mqc_dirs = get_multiQC_params,
        output_dir = f'{output_directory}/analysis/multiqc',
    log:
        f'{output_directory}/logs/multiqc.log'
    benchmark:
        f'{output_directory}/benchmarks/multiQC.txt'
    threads: 1
    resources:
        mem_gb=8,
        walltime = config['walltime']['medium']
    envmodules:
        config['envmodules']['multiqc'],
    shell:
        """
        multiqc -f -o {params.output_dir} -n multiqc_report.html {params.mqc_dirs} 2> {log}
        """

