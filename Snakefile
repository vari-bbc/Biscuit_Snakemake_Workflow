import pandas as pd
import numpy as np
import os
import re
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.20.1")

samples = pd.read_table("bin/samples.tsv", dtype=str).set_index(["sample"], drop=False)
configfile: "bin/config.yaml"

rule all:
    input:
        # rename
        # expand("raw_data/{samples.sample}-R{read}.fastq.gz", read=[1,2], samples=samples.itertuples()),
        # trim_galore
        # expand("analysis/trim_reads/{samples.sample}-R{read}.fastq.gz_trimming_report.txt", read=[1,2], samples=samples.itertuples()),
        # expand("analysis/trim_reads/{samples.sample}-R{read}_val_{read}.fq.gz", read=[1,2], samples=samples.itertuples()),
        # biscuit
        # expand("analysis/align/{samples.sample}.sorted.markdup.bam", samples=samples.itertuples()),
        # expand("analysis/align/{samples.sample}.sorted.markdup.bam.bai", samples=samples.itertuples()),
        # expand("analysis/align/{samples.sample}.sorted.markdup.bam.flagstat", samples=samples.itertuples()),
        # biscuit_pileup
        # expand("analysis/pileup/{samples.sample}.vcf.gz", samples=samples.itertuples()),
        # expand("analysis/pileup/{samples.sample}.bed.gz", samples=samples.itertuples()),
        # expand("analysis/pileup/{samples.sample}.bed.gz.tbi", samples=samples.itertuples()),
        # mergecg
        expand("analysis/pileup/{samples.sample}_mergecg.bed.gz", samples=samples.itertuples()),
        expand("analysis/pileup/{samples.sample}_mergecg.bed.gz.tbi", samples=samples.itertuples()),
        # mergecg_combined
        # "analysis/pileup/combined_mergecg.bed.gz",
        # "analysis/pileup/combined_mergecg.bed.gz.tbi",
        # biscuit_qc
        # expand("analysis/BISCUITqc/{samples.sample}_mapq_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_strand_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_isize_score_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_covdist_cpg_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_covdist_cpg_q40_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_covdist_cpg_q40_botgc_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_covdist_cpg_q40_topgc_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_covdist_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_covdist_q40_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_covdist_q40_botgc_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_covdist_q40_topgc_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_all_cv_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_cpg_cv_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_cpg_dist_table.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_dup_report.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_CpHRetentionByReadPos.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_CpGRetentionByReadPos.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_CpGRetentionDist.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_freqOfTotalRetentionPerRead.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_totalBaseConversionRate.txt", samples=samples.itertuples()),
        # expand("analysis/BISCUITqc/{samples.sample}_totalReadConversionRate.txt", samples=samples.itertuples()),
        # multiQC
        "analysis/multiqc/multiqc_report.html",
        # biscuiteer

rule rename_fastq:
    output:
        "raw_data/{sample}-R1.fastq.gz",
        "raw_data/{sample}-R2.fastq.gz",
    log:
        "logs/rename/rename_{sample}.log"
    threads: 1
    resources:
        mem_gb=8
    envmodules:
        config["envmodules"]["R"],
        config["envmodules"]["snakemake"],
    script:
        "bin/rename.R"

rule trim_galore:
    input:
        "raw_data/{sample}-R1.fastq.gz",
        "raw_data/{sample}-R2.fastq.gz"
    output:
        "analysis/trim_reads/{sample}-R1_val_1.fq.gz",
        "analysis/trim_reads/{sample}-R1_val_1_fastqc.html",
        "analysis/trim_reads/{sample}-R1_val_1_fastqc.zip",
        "analysis/trim_reads/{sample}-R1.fastq.gz_trimming_report.txt",
        "analysis/trim_reads/{sample}-R2_val_2.fq.gz",
        "analysis/trim_reads/{sample}-R2_val_2_fastqc.html",
        "analysis/trim_reads/{sample}-R2_val_2_fastqc.zip",
        "analysis/trim_reads/{sample}-R2.fastq.gz_trimming_report.txt"
    log:
        stdout="logs/trim_reads/{sample}.o",
        stderr="logs/trim_reads/{sample}.e"
    benchmark:
        "benchmarks/trim_reads/{sample}.txt"
    params:
        outdir = "analysis/trim_reads/",
        quality = config["trim_galore"]["q"],
    envmodules:
        config["envmodules"]["trim_galore"],
        config["envmodules"]["fastqc"],
        config["envmodules"]["pigz"],
    threads: 4
    resources:
        mem_gb=80
    shell:
        """
        trim_galore \
        --paired \
        {input} \
        --output_dir {params.outdir} \
        --cores {threads} \
        -q {params.quality} \
        --fastqc \
        2> {log.stderr} 1> {log.stdout}
        """

rule biscuit_align:
    input:
        R1 = "analysis/trim_reads/{sample}-R1_val_1.fq.gz",
        R2 = "analysis/trim_reads/{sample}-R2_val_2.fq.gz",
    output:
        bam = "analysis/align/{sample}.sorted.markdup.bam",
        bai = "analysis/align/{sample}.sorted.markdup.bam.bai",
        disc = "analysis/align/{sample}.disc.sorted.bam",
        split = "analysis/align/{sample}.split.sorted.bam",
        unmapped = "analysis/align/{sample}.unmapped.fastq.gz",
        flagstat = "analysis/align/{sample}.sorted.markdup.bam.flagstat",
    params:
        # don't include the .fa/.fasta suffix for the reference biscuit idx.
        ref = config["ref"]["index"],
        LB = config["sam_header"]["LB"],
        ID = "{sample}",
        PL = config["sam_header"]["PL"],
        PU = config["sam_header"]["PU"],
        SM = "{sample}",
        lib_type = config["biscuit"]["lib_type"],
        disc = "analysis/align/{sample}.disc.sam",
        split = "analysis/align/{sample}.split.sam",
        unmapped = "analysis/align/{sample}.unmapped.fastq",
    log:
        biscuit = "logs/biscuit/biscuit_align.{sample}.log",
        samblaster = "logs/biscuit/samblaster.{sample}.log",
        samtools_view = "logs/biscuit/samtools_view.{sample}.log",
        samtools_sort = "logs/biscuit/samtools_sort.{sample}.log",
        samtools_index = "logs/biscuit/samtools_index.{sample}.log",
        samtools_flagstat = "logs/biscuit/samtools_flagstat.{sample}.log",
        sort_disc =  "logs/biscuit/sort_disc.{sample}.log",
        index_disc = "logs/biscuit/index_disc.{sample}.log",    
        sort_split = "logs/biscuit/sort_split.{sample}.log",
        index_split = "logs/biscuit/index_split.{sample}.log",
        bgzip_unmapped = "logs/biscuit/bgzip_unmapped.{sample}.log",       
    threads: 32
    resources:
        mem_gb=500
    benchmark:
        "benchmarks/biscuit_align/{sample}.txt"
    envmodules:
        config["envmodules"]["biscuit"],
        config["envmodules"]["samtools"],
        config["envmodules"]["samblaster"],
        config["envmodules"]["bedtools"],
        config["envmodules"]["htslib"],
    shell:
        """
        biscuit align -M -t {threads} -b {params.lib_type} \
        -R '@RG\tLB:{params.LB}\tID:{params.ID}\tPL:{params.PL}\tPU:{params.PU}\tSM:{params.SM}' \
        {params.ref} {input.R1} {input.R2} 2> {log.biscuit} | \
        samblaster -M -r --addMateTags -d {params.disc} -s {params.split} -u {params.unmapped} 2> {log.samblaster} | \
        samtools view -hbu -F 4 -q 20 2> {log.samtools_view} |
        samtools sort -@ {threads} -m 5G -o {output.bam} -O BAM - 2> {log.samtools_sort}
        samtools index -@ {threads} {output.bam} 2> {log.samtools_index}
        samtools flagstat {output.bam} 1> {output.flagstat} 2> {log.samtools_flagstat}
        samtools sort -o {output.disc} -O BAM {params.disc} 2> {log.sort_disc}
        samtools index -@ {threads} {output.disc} 2> {log.index_disc}
        samtools sort -o {output.split} -O BAM {params.split} 2> {log.sort_split}
        samtools index -@ {threads} {output.split} 2> {log.index_split}
        bgzip -@ {threads} {params.unmapped} 2> {log.bgzip_unmapped}
        rm {params.disc}
        rm {params.split}
        """

rule biscuit_pileup:
    input:
        bam="analysis/align/{sample}.sorted.markdup.bam",
    params:
        ref=config["ref"]["fasta"],
        vcf="analysis/pileup/{sample}.vcf",
        bed="analysis/pileup/{sample}.bed",
    output:
        vcf_gz="analysis/pileup/{sample}.vcf.gz",
        vcf_tabix="analysis/pileup/{sample}.vcf.gz.tbi",
        meth="analysis/pileup/{sample}.vcf_meth_average.tsv",
        bed_gz="analysis/pileup/{sample}.bed.gz",
        bed_tbi="analysis/pileup/{sample}.bed.gz.tbi",
    log:
        pileup = "logs/biscuit_pileup/pileup.{sample}.log",
        vcf_gz = "logs/biscuit_pileup/vcf_gz.{sample}.log",
        vcf_tbi = "logs/biscuit_pileup/vcf_tbi.{sample}.log",
        vcf2bed = "logs/biscuit_pileup/vcf2bed.{sample}.log",
        bed_gz="logs/biscuit_pileup/bed_gz.{sample}.log",
        bed_tbi="logs/biscuit_pileup/bed_tabix.{sample}.log",
    threads: 8
    resources:
        mem_gb=100
    benchmark:
        "benchmarks/biscuit_pileup/{sample}.txt"
    wildcard_constraints:
        sample=".*[^(_mergecg)]"
    envmodules:
        config["envmodules"]["biscuit"],
        config["envmodules"]["htslib"],
        config["envmodules"]["bedtools"],
    shell:
        """
        biscuit pileup -q {threads} -o {params.vcf} {params.ref} {input.bam} 2> {log.pileup}
        bgzip {params.vcf} 2> {log.vcf_gz}
        tabix -p vcf {output.vcf_gz} 2> {log.vcf_tbi}

        biscuit vcf2bed -t cg {output.vcf_gz} 1> {params.bed} 2> {log.vcf2bed}
        bgzip {params.bed} 2> {log.bed_gz}
        tabix -p bed {output.bed_gz} 2> {log.bed_tbi}
        """

rule biscuit_mergecg:
    input:
        bed="analysis/pileup/{sample}.bed.gz",
    params:
        ref=config["ref"]["fasta"],
        mergecg="analysis/pileup/{sample}_mergecg.bed",
    output:
        mergecg_gz="analysis/pileup/{sample}_mergecg.bed.gz",
        mergecg_tbi="analysis/pileup/{sample}_mergecg.bed.gz.tbi",
    log:
        mergecg="logs/biscuit_pileup/mergecg.{sample}.log",
        mergecg_gz="logs/biscuit_pileup/mergecg_gz.{sample}.log",
        mergecg_tbi="logs/biscuit_pileup/mergecg_tabix.{sample}.log",
    threads: 8
    resources:
        mem_gb=100
    benchmark:
        "benchmarks/biscuit_mergecg/{sample}.txt"
    wildcard_constraints:
        sample=".*[^(_mergecg)]"
    envmodules:
        config["envmodules"]["biscuit"],
        config["envmodules"]["htslib"],
    shell:
        """
        biscuit mergecg {params.ref} {input.bed} 1> {params.mergecg} 2> {log.mergecg}
        bgzip {params.mergecg} 2> {log.mergecg_gz}
        tabix -p bed {output.mergecg_gz} 2> {log.mergecg_tbi}
        """

rule biscuit_qc:
    input:
        vcf="analysis/pileup/{sample}.vcf.gz",
        bam="analysis/align/{sample}.sorted.markdup.bam"
    params:
        ref=config["ref"]["index"],
        assets=config["ref"]["assets"],
        sample="{sample}",
    output:
        "analysis/BISCUITqc/{sample}_covdist_all_base_botgc_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_all_base_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_all_base_topgc_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_all_cpg_botgc_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_all_cpg_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_all_cpg_topgc_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_q40_base_botgc_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_q40_base_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_q40_base_topgc_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_q40_cpg_botgc_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_q40_cpg_table.txt",
        "analysis/BISCUITqc/{sample}_covdist_q40_cpg_topgc_table.txt",
        "analysis/BISCUITqc/{sample}_cv_table.txt",
        "analysis/BISCUITqc/{sample}_mapq_table.txt",
        "analysis/BISCUITqc/{sample}_strand_table.txt",
        "analysis/BISCUITqc/{sample}_dup_report.txt",
        "analysis/BISCUITqc/{sample}_totalBaseConversionRate.txt",
        "analysis/BISCUITqc/{sample}_totalReadConversionRate.txt",
        "analysis/BISCUITqc/{sample}_CpHRetentionByReadPos.txt",
        "analysis/BISCUITqc/{sample}_CpGRetentionByReadPos.txt",
    threads: 8
    resources:
        mem_gb=100
    benchmark:
        "benchmarks/biscuit_qc/{sample}.txt"
    log:
        "logs/biscuit_qc/{sample}_QC.log"
    envmodules:
        config["envmodules"]["biscuit"],
        config["envmodules"]["samtools"],
        config["envmodules"]["htslib"],
        config["envmodules"]["bedtools"],
        config["envmodules"]["parallel"],
    shell:
        """
        set +o pipefail;
        QC.sh \
        -o analysis/BISCUITqc \
        --vcf {input.vcf} \
        {params.assets} \
        {params.ref} \
        {params.sample} \
        {input.bam} \
        2> {log}
        """

rule multiQC:
    input:
        # trim_galore
        expand("analysis/trim_reads/{samples.sample}-R{read}.fastq.gz_trimming_report.txt", read=[1,2], samples=samples.itertuples()),
        expand("analysis/trim_reads/{samples.sample}-R{read}_val_{read}.fq.gz", read=[1,2], samples=samples.itertuples()),
        # biscuit_qc
        expand("analysis/BISCUITqc/{samples.sample}_covdist_all_base_botgc_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_all_base_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_all_base_topgc_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_all_cpg_botgc_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_all_cpg_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_all_cpg_topgc_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_q40_base_botgc_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_q40_base_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_q40_base_topgc_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_q40_cpg_botgc_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_q40_cpg_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_covdist_q40_cpg_topgc_table.txt", samples=samples.itertuples()),
        expand("analysis/BISCUITqc/{samples.sample}_cv_table.txt", samples=samples.itertuples()),
    params:
        "raw_data/",
        "analysis/trim_reads/",
        "analysis/BISCUITqc/",
    output:
        directory("analysis/multiqc/multiqc_report_data",),
        "analysis/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log"
    threads: 1
    resources:
        mem_gb=32
    benchmark:
        "benchmarks/multiQC.txt"
    envmodules:
        config["envmodules"]["multiqc"],
        config["envmodules"]["fastqc"],
    shell:
        """
        multiqc -f {params} \
        -o analysis/multiqc \
        -n multiqc_report.html 2> {log}
        """
