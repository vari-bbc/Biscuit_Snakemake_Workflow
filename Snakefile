import pandas as pd
import numpy as np
import os
import re
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.20.1")

samples = pd.read_table("bin/samples.tsv", dtype=str).set_index(["sample"], drop=False)
configfile: "bin/config.yaml"
biscuitIndexFORMATS = ["bis.ann", "bis.amb","par.bwt","dau.bwt","bis.pac","par.sa","dau.sa","fai"]
wildcard_constraints:
   seqfile_index = '\d+'

rule all:
    input:
        # Sliding window / binned Averages
        # ~ "analysis/binnedAverages/" if config["SlidingWindow"]["run"] else [],
        # Download CpG Islands
        # ~ "analysis/cpgIslands/" if config["SubsetToFeatures"]["run"] else [],
        # rename_fastq
        expand("raw_data/{samples.sample}-1-R1.fastq.gz", samples=samples.itertuples()),
        # ~ # biscuit
        expand("analysis/align/{samples.sample}.sorted.markdup.bam", samples=samples.itertuples()),
        # ~ # mergecg
        expand("analysis/pileup/{samples.sample}_mergecg.bed.gz", samples=samples.itertuples()),
        expand("analysis/pileup/{samples.sample}_mergecg.bed.gz.tbi", samples=samples.itertuples()),
        # multiQC
        "analysis/multiqc/multiqc_report.html",
        # rule for checking that data is not trimmed
        
        # rules for quality control vectors
        expand("analysis/qc_vectors/lambda/{samples.sample}.bed", samples=samples.itertuples()) if config["control_vectors"] else [],
        expand("analysis/qc_vectors/puc19/{samples.sample}.bed", samples=samples.itertuples()) if config["control_vectors"] else [],
        "analysis/qc_vectors/control_vector_boxplot.pdf" if config["control_vectors"] else [], # for qc_vectors plot
        
        # rules for fastq_screen
        expand("analysis/fastq_screen/{samples.sample}-1-R{read}_screen.html", read=[1,2], samples=samples.itertuples()) if config["run_fastq_screen"] else [], # for fastQC screen
        
        # output of rule build_ref_with_methylation_controls
        expand("snakemake_built_reference_with_methylation_controls/merged.fa.gz.{ext}", ext=biscuitIndexFORMATS) if config["build_ref_with_methylation_controls"] else [],
        
        # epiread
        expand("analysis/snps/{samples.sample}.snp.bed.gz", samples=samples.itertuples()) if config["epiread"] else [],
        expand("analysis/epiread/{samples.sample}.epibed.gz", samples=samples.itertuples()) if config["epiread"] else [],
        # snps
        expand("analysis/snps/{samples.sample}.snp.bed.gz", samples=samples.itertuples()) if config["generate_snps"] else [],
               
rule get_R1_R2_files:
    output:
        "raw_data/{sample}-1-R1.fastq.gz", # only require 1 file / sample
        "raw_data/{sample}-1-R2.fastq.gz",
    log:
        "logs/rename/rename_{sample}.log"
    threads: 1
    resources:
        mem_gb=8,
        walltime = config["walltime"]["short"]
    envmodules:
        config["envmodules"]["R"],
        config["envmodules"]["snakemake"],
    script:
        "bin/rename.R"

if config["build_ref_with_methylation_controls"]:
    assert bool(re.match(".*\.gz$", config["ref"]["fasta"])), "Reference fasta in the config.yaml needs to be gzipped!"  
    newRef = "snakemake_built_reference_with_methylation_controls/merged.fa.gz"
    newIndex = newRef
    rule build_ref_with_methylation_controls:
        input:
            config["ref"]["fasta"],
        output:
            ref = newRef,
            newrefdir = directory('snakemake_built_reference_with_methylation_controls/'),
            indexes = expand("snakemake_built_reference_with_methylation_controls/merged.fa.gz.{ext}", ext=biscuitIndexFORMATS)
        log:
            "logs/control_vector_genome_build_index.log"
        threads: 
            2,
        resources:
            mem_gb=32
        benchmark:
            "benchmarks/build_ref_with_methylation_controls.txt"    
        envmodules:
            config["envmodules"]["R"],
            config["envmodules"]["samtools"],
            config["envmodules"]["biscuit"],
        shell:
            """
            mkdir -p {output.newrefdir}
            # should use gzip, NOT bgzip!
            cat bin/puc19.fa.gz bin/lambda.fa.gz {input} > {output.ref}
            biscuit index {output.ref}
            samtools faidx {output.ref}
            """
else:
    newRef = config["ref"]["fasta"],
    newIndex = config["ref"]["index"], 


def get_trim_reads_index(wildcards): # this is for getting files when there is more than 1 R1 and R2
    FILE_INDEX, = glob_wildcards("raw_data/" + wildcards.sample + "-{id}-R1.fastq.gz")
    return list(FILE_INDEX)
    
def get_trim_reads_input(wildcards):
    FILE_INDEX, = glob_wildcards("raw_data/" + wildcards.sample + "-{id}-R1.fastq.gz") # R1 and R2 must be the same
    files = list(expand("raw_data/" + wildcards.sample + "-{seqfile_index}-R{read}.fastq.gz", seqfile_index = FILE_INDEX, read = [1,2]))
    # ~ return FILE_INDEX
    return files

rule trim_reads:
    input:
        get_trim_reads_input
    output:
        "analysis/trim_reads/{sample}-R1_val_1_merged.fq.gz",
        # ~ "analysis/trim_reads/{sample}-1-R1_val_1_fastqc.html",
        # ~ "analysis/trim_reads/{sample}-1-R1_val_1_fastqc.zip",
        # ~ "analysis/trim_reads/{sample}-1-R1.fastq.gz_trimming_report.txt",
        "analysis/trim_reads/{sample}-R2_val_2_merged.fq.gz",
        # ~ "analysis/trim_reads/{sample}-1-R2_val_2_fastqc.html",
        # ~ "analysis/trim_reads/{sample}-1-R2_val_2_fastqc.zip",
        # ~ "analysis/trim_reads/{sample}-1-R2.fastq.gz_trimming_report.txt"
    log:
        stdout="logs/trim_reads/{sample}.o",
        stderr="logs/trim_reads/{sample}.e"
    benchmark:
        "benchmarks/trim_reads/{sample}.txt"
    params:
        sample = "{sample}",
        outdir = "analysis/trim_reads/",
        quality = config["trim_galore"]["quality"],
        hard_trim_R2 = config["trim_galore"]["hard_trim_R2"],
    envmodules:
        config["envmodules"]["trim_galore"],
        config["envmodules"]["fastqc"],
        config["envmodules"]["pigz"],
    threads: config["hpcParameters"]["trimThreads"]
    resources:
        mem_gb=config["hpcParameters"]["smallMemoryGb"],
        walltime = config["walltime"]["medium"]
    shell:
        """
        echo {input}
        if [ {params.hard_trim_R2} -ge 1 ]; then
            trim_galore \
            --paired \
            {input} \
            --output_dir {params.outdir} \
            --clip_R2 {params.hard_trim_R2} \
            --cores {threads} \
            -q {params.quality} \
            --fastqc \
            2> {log.stderr} 1> {log.stdout}
        else
            trim_galore \
            --paired \
            {input} \
            --output_dir {params.outdir} \
            --cores {threads} \
            -q {params.quality} \
            --fastqc \
            2> {log.stderr} 1> {log.stdout}
        fi
        
        cat {params.outdir}{wildcards.sample}-*-R1_val_1.fq.gz > {params.outdir}{wildcards.sample}-R1_val_1_merged.fq.gz # make the merged R1
        cat {params.outdir}{wildcards.sample}-*-R2_val_2.fq.gz > {params.outdir}{wildcards.sample}-R2_val_2_merged.fq.gz # make the merged R2
        rm {params.outdir}{wildcards.sample}-*-R1_val_1.fq.gz
        rm {params.outdir}{wildcards.sample}-*-R2_val_2.fq.gz
        """

if config["run_fastq_screen"]:
    rule fastq_screen:
        input:
            read1 = "raw_data/{sample}-1-R1.fastq.gz", # ***only the first if there are multiple files***
            read2 = "raw_data/{sample}-1-R2.fastq.gz",
        output:
            "analysis/fastq_screen/{sample}-1-R1_screen.html",
            "analysis/fastq_screen/{sample}-1-R2_screen.html",
            "analysis/fastq_screen/{sample}-1-R1_screen.txt",
            "analysis/fastq_screen/{sample}-1-R2_screen.txt",
        log:
            fastq_screen = "logs/fastq_screen/{sample}.log",
        params:
            config = config["fastq_screen_conf"],
        benchmark:
            "benchmarks/fastq_screen/{sample}.bmk"
        threads: 8
        resources:
            nodes = 1,
            mem_gb = 64,
        envmodules: 
            config["envmodules"]["fastq_screen"],
            config["envmodules"]["bismark"],
        shell:
            """
            fastq_screen --bisulfite --conf {params.config} --outdir analysis/fastq_screen/ {input} 2> {log.fastq_screen}
            """

def get_biscuit_align_reference(wildcards):
    if config["build_ref_with_methylation_controls"]:
        input = expand("snakemake_built_reference_with_methylation_controls/merged.fa.gz.{ext}", ext=biscuitIndexFORMATS)
        return input
    else:
        input = config["ref"]["fasta"] # else just require default reference
        return input

def get_rename_fastq_output_R1(wildcards):
    if config["trim_galore"]["trim_before_BISCUIT"]:
        files = "analysis/trim_reads/" + wildcards.sample + "-R1_val_1_merged.fq.gz"
        return files   
    else:
        FILE_INDEX, = glob_wildcards("raw_data/" + wildcards.sample + "-{id}-R1.fastq.gz")
        files = list(expand("raw_data/" + wildcards.sample + "-{seqfile_index}-R1.fastq.gz", seqfile_index = FILE_INDEX))
        files.sort()
        return files
        
def get_rename_fastq_output_R2(wildcards):
    if config["trim_galore"]["trim_before_BISCUIT"]:
        files = "analysis/trim_reads/" + wildcards.sample + "-R2_val_2_merged.fq.gz"
        return files   
    else:
        FILE_INDEX, = glob_wildcards("raw_data/" + wildcards.sample + "-{id}-R2.fastq.gz")
        files = list(expand("raw_data/" + wildcards.sample + "-{seqfile_index}-R2.fastq.gz", seqfile_index = FILE_INDEX))
        files.sort()
        return files

if config["SlidingWindow"]["run"]:  
    rule slidingWindow:
        input:
            # pileup
            # ~ expand("analysis/pileup/{samples.sample}.bed.gz", samples=samples.itertuples()),
            # ~ expand("analysis/pileup/{samples.sample}.bed.gz.tbi", samples=samples.itertuples()),
        params:
            covFilter = config["SlidingWindow"]["covFilter"],
        output:
            dir = directory("analysis/binnedAverages/",)
        log:
            "logs/slidingWindow.log"
        threads: 1
        resources:
            mem_gb = config["hpcParameters"]["intermediateMemoryGb"],
            walltime = config["walltime"]["medium"]
        benchmark:
            "benchmarks/slidingWindow.txt"
        envmodules:
            config["envmodules"]["multiqc"],
            config["envmodules"]["fastqc"],
        shell:
            """
            echo "hewo world"
            mkdir -p {output.dir}
            window=(1000 10000 100000 1000000 10000000)
            for i in "${{window[@]}}"
            do
                ./bin/slidingWindow.pl -indir analysis/pileup/ -windowSize ${{i}} -covFilter {params.covFilter} -out {output.dir}/${{i}}.tsv
            done
            """   

if config["SubsetToFeatures"]["run"]:  
    rule cpgIslands:
        input:
            reference = get_biscuit_align_reference,
            # ~ expand("analysis/pileup/{samples.sample}.bed.gz", samples=samples.itertuples()),
            # ~ expand("analysis/pileup/{samples.sample}.bed.gz.tbi", samples=samples.itertuples()),
        params:
        output:
            dir = directory("analysis/cpgIslands/",)
        log:
            "logs/cpgi.log"
        threads: 1
        resources:
            mem_gb = config["hpcParameters"]["intermediateMemoryGb"],
            walltime = config["walltime"]["medium"]
        benchmark:
            "benchmarks/cpgi.txt"
        envmodules:
            config["envmodules"]["bedtools"],
        shell:
            """
            mkdir -p {output.dir}
            cd {output.dir}
            
            # Reference location
            REFLOC={input.reference}.fai

            # CpG Islands from UCSC on only the canonical chromosomes
            wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz | \
            gunzip -c | \
            awk 'BEGIN{{ OFS="\\t"; }}{{ print $2, $3, $4; }}' | \
            awk '{{ if ($1 ~ /^chr[1234567890XYM]{{1,2}}$/) {{ print }} }}' | \
            bedtools sort -i - \
            > cpg_islands.bed

            # Create middles for finding locations
            bedtools slop \
                -i cpg_islands.bed \
                -g ${{REFLOC}} \
                -b 2000 | \
            bedtools merge -i - \
            > shores.tmp.bed

            bedtools slop \
                -i cpg_islands.bed \
                -g ${{REFLOC}} \
                -b 4000 | \
            bedtools merge -i - \
            > shelves.tmp.bed

            # CpG Open Seas (intervening locations)
            sort -k1,1 -k2,2n ${{REFLOC}} | \
            bedtools complement -L -i shelves.tmp.bed -g - \
            > cpg_open_seas.bed

            # CpG Shelves (shores +/- 2kb)
            bedtools subtract \
                -a shelves.tmp.bed \
                -b shores.tmp.bed \
            > cpg_shelves.bed

            # CpG Shores (island +/- 2kb)
            bedtools subtract \
                -a shores.tmp.bed \
                -b cpg_islands.bed \
            > cpg_shores.bed

            # Clean up temporary files
            rm -f shelves.tmp.bed shores.tmp.bed
            """       

rule biscuit_align:
    input:
        reference = get_biscuit_align_reference,
        R1 = get_rename_fastq_output_R1,
        R2 = get_rename_fastq_output_R2
    output:
        bam = "analysis/align/{sample}.sorted.markdup.bam",
        bai = "analysis/align/{sample}.sorted.markdup.bam.bai",
        disc = "analysis/align/{sample}.disc.sorted.bam" if config["biscuit"]["biscuit_blaster_version"] == "v2" else ["analysis/align/{sample}.no.disc"],
        disc_bai = "analysis/align/{sample}.disc.sorted.bam.bai" if config["biscuit"]["biscuit_blaster_version"] == "v2" else ["analysis/align/{sample}.no.disc.bai"],
        split = "analysis/align/{sample}.split.sorted.bam" if config["biscuit"]["biscuit_blaster_version"] == "v2" else ["analysis/align/{sample}.no.split"],
        split_bai = "analysis/align/{sample}.split.sorted.bam.bai" if config["biscuit"]["biscuit_blaster_version"] == "v2" else ["analysis/align/{sample}.no.split.bai"],
        unmapped = "analysis/align/{sample}.unmapped.fastq.gz" if config["biscuit"]["biscuit_blaster_version"] == "v2" else ["analysis/align/{sample}.no.unmapped"],
        flagstat = "analysis/align/{sample}.sorted.markdup.bam.flagstat",
    params:
        # don't include the .fa/.fasta suffix for the reference biscuit idx.
        ref = newRef,
        LB = config["sam_header"]["LB"],
        ID = "{sample}",
        PL = config["sam_header"]["PL"],
        PU = config["sam_header"]["PU"],
        SM = "{sample}",
        lib_type = config["biscuit"]["lib_type"],
        disc = "analysis/align/{sample}.disc.sam",
        split = "analysis/align/{sample}.split.sam",
        unmapped = "analysis/align/{sample}.unmapped.fastq",
        biscuit_version = config["biscuit"]["biscuit_blaster_version"]
    log:
        biscuit = "logs/biscuit/biscuit_align.{sample}.log",
        biscuit_blaster_version = "logs/biscuit/blaster_version.{sample}.log",
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
    threads: config["hpcParameters"]["maxThreads"]
    resources:
        mem_gb = config["hpcParameters"]["maxMemoryGb"],
        walltime = config["walltime"]["long"]
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
        if [ {params.biscuit_version} == "v2" ]; then
            echo "biscuit blaster v2" 2> {log.biscuit_blaster_version}
            biscuit align -@ {threads} -b {params.lib_type} \
            -R '@RG\tLB:{params.LB}\tID:{params.ID}\tPL:{params.PL}\tPU:{params.PU}\tSM:{params.SM}' \
            {params.ref} <(zcat {input.R1}) <(zcat {input.R2}) 2> {log.biscuit} | \
            samblaster -r --addMateTags -d {params.disc} -s {params.split} -u {params.unmapped} 2> {log.samblaster} | \
            samtools view -hbu -F 4 -q 30 2> {log.samtools_view} |
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
        elif [ {params.biscuit_version} == "v1" ]; then
            echo "biscuit blaster v1" 2> {log.biscuit_blaster_version}
            biscuit align -@ {threads} -b {params.lib_type} \
            -R '@RG\tLB:{params.LB}\tID:{params.ID}\tPL:{params.PL}\tPU:{params.PU}\tSM:{params.SM}' \
            {params.ref} <(zcat {input.R1}) <(zcat {input.R2}) 2> {log.biscuit} | \
            samblaster -r --addMateTags 2> {log.samblaster} | \
            samtools view -hbu -F 4 -q 30 2> {log.samtools_view} |
            samtools sort -@ {threads} -m 5G -o {output.bam} -O BAM - 2> {log.samtools_sort}
            samtools index -@ {threads} {output.bam} 2> {log.samtools_index}
            samtools flagstat {output.bam} 1> {output.flagstat} 2> {log.samtools_flagstat}
            touch {output.disc}
            touch {output.disc_bai}
            touch {output.split}
            touch {output.split_bai}
            touch {output.unmapped}
        else
            echo "biscuit: biscuit_blaster_version must be v1 or v2 in bin/config.yaml" 2> {log.biscuit_blaster_version}
        fi
        """
        
rule biscuit_pileup:
    input:
        bam="analysis/align/{sample}.sorted.markdup.bam",
    params:
        ref=newRef,
        vcf="analysis/pileup/{sample}.vcf",
        bed="analysis/pileup/{sample}.bed",
        nome=config["is_nome"]
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
    threads: config["hpcParameters"]["pileupThreads"]
    resources:
        mem_gb = config["hpcParameters"]["intermediateMemoryGb"],
        walltime = config["walltime"]["medium"]
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
        if [ {params.nome} == "TRUE" ]; then
            biscuit pileup -N -@ {threads} -o {params.vcf} {params.ref} {input.bam} 2> {log.pileup}
        else
            biscuit pileup -@ {threads} -o {params.vcf} {params.ref} {input.bam} 2> {log.pileup}
        fi
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
        ref=newRef,
        mergecg="analysis/pileup/{sample}_mergecg.bed",
        nome=config["is_nome"]
    output:
        mergecg_gz="analysis/pileup/{sample}_mergecg.bed.gz",
        mergecg_tbi="analysis/pileup/{sample}_mergecg.bed.gz.tbi",
    log:
        mergecg="logs/biscuit_pileup/mergecg.{sample}.log",
        mergecg_gz="logs/biscuit_pileup/mergecg_gz.{sample}.log",
        mergecg_tbi="logs/biscuit_pileup/mergecg_tabix.{sample}.log",
    threads: 8
    resources:
        mem_gb = config["hpcParameters"]["intermediateMemoryGb"],
        walltime = config["walltime"]["medium"]
    benchmark:
        "benchmarks/biscuit_mergecg/{sample}.txt"
    wildcard_constraints:
        sample=".*[^(_mergecg)]"
    envmodules:
        config["envmodules"]["biscuit"],
        config["envmodules"]["htslib"],
    shell:
        """
        if [ {params.nome} == "TRUE" ]; then
            biscuit mergecg -N {params.ref} {input.bed} 1> {params.mergecg} 2> {log.mergecg}
        else
            biscuit mergecg {params.ref} {input.bed} 1> {params.mergecg} 2> {log.mergecg}
        fi
        bgzip {params.mergecg} 2> {log.mergecg_gz}
        tabix -p bed {output.mergecg_gz} 2> {log.mergecg_tbi}
        """

rule biscuit_qc:
    input:
        vcf="analysis/pileup/{sample}.vcf.gz",
        bam="analysis/align/{sample}.sorted.markdup.bam"
    params:
        ref=newRef,
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
        mem_gb = config["hpcParameters"]["intermediateMemoryGb"],
        walltime = config["walltime"]["medium"]
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

def get_multiQC_params(wildcards):
    if config["run_fastq_screen"] and config["trim_galore"]["trim_before_BISCUIT"]:
        input = "raw_data/ analysis/trim_reads/ analysis/BISCUITqc/ analysis/fastq_screen"
        return input
    elif config["run_fastq_screen"]:
        input = "raw_data/ analysis/BISCUITqc/ analysis/fastq_screen"
        return input
    elif config["trim_galore"]["trim_before_BISCUIT"]:
        input = "raw_data/ analysis/BISCUITqc/ analysis/trim_reads/"
        return input
    else:
        input = "raw_data/ analysis/BISCUITqc/"
        return input
        
rule multiQC:
    input:
        # fastq_screen
        expand("analysis/fastq_screen/{samples.sample}-1-R{read}_screen.html", read=[1,2], samples=samples.itertuples()) if config["run_fastq_screen"] else [],
        expand("analysis/fastq_screen/{samples.sample}-1-R{read}_screen.txt", read=[1,2], samples=samples.itertuples()) if config["run_fastq_screen"] else [],
        # trim_galore
        # ~ expand("analysis/trim_reads/{samples.sample}-1-R{read}.fastq.gz_trimming_report.txt", read=[1,2], samples=samples.itertuples()) if config["trim_galore"]["trim_before_BISCUIT"] else [],
        expand("analysis/trim_reads/{samples.sample}-R{read}_val_{read}_merged.fq.gz", read=[1,2], samples=samples.itertuples()) if config["trim_galore"]["trim_before_BISCUIT"] else [],
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
        get_multiQC_params
    output:
        directory("analysis/multiqc/multiqc_report_data",),
        "analysis/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log"
    threads: 1
    resources:
        mem_gb=8,
        walltime = config["walltime"]["medium"]
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

if config["control_vectors"]:
    rule methylation_controls_qc:
        input:
           bed="analysis/pileup/{sample}_mergecg.bed.gz",
        envmodules:
           config["envmodules"]["samtools"],
           config["envmodules"]["htslib"],
        output:
           lambda_bed = "analysis/qc_vectors/lambda/{sample}.bed",
           puc19_bed = "analysis/qc_vectors/puc19/{sample}.bed",
        resources:
            mem_gb=32,
            walltime = config["walltime"]["medium"]
        benchmark:
            "benchmarks/methylation_controls_qc/{sample}.txt"
        log:
           lambda_log = "logs/qc_vectors/lambda.{sample}_QC.log",
           puc19_log = "logs/qc_vectors/puc19.{sample}_QC.log"
        shell:
           """
           mkdir -p "analysis/qc_vectors/lambda/"
           mkdir -p "analysis/qc_vectors/puc19/"
           
           # >J02459.1 Escherichia phage Lambda, complete genome - UNMETHYLATED CONTROL
           zcat {input.bed} | {{ grep '^J02459.1' || true; }} > {output.lambda_bed} 2> {log.lambda_log}
           
           # >M77789.2 Cloning vector pUC19, complete sequence - METHYLATED CONTROL
           zcat {input.bed} | {{ grep '^M77789.2' || true; }}  > {output.puc19_bed} 2> {log.puc19_log}
           """
       
    rule methylation_controls_figure:
        input:
           expand("analysis/qc_vectors/lambda/{samples.sample}.bed", samples=samples.itertuples()),
           expand("analysis/qc_vectors/puc19/{samples.sample}.bed", samples=samples.itertuples()),
        envmodules:
           config["envmodules"]["R"],
        resources:
            mem_gb=32,
            walltime = config["walltime"]["short"]
        output:
           control_vector_pdf = "analysis/qc_vectors/control_vector_boxplot.pdf"
        log:
           control_vector_pdf = "logs/qc_vectors/control_vector_pdf.log",
        script:
            "bin/control_vector.R"       
       
if config["generate_snps"] or config["epiread"]:
    rule biscuit_snps:
        input: 
           vcf_gz = "analysis/pileup/{sample}.vcf.gz",
           bam="analysis/align/{sample}.sorted.markdup.bam",
        output:
           snp_bed_gz = "analysis/snps/{sample}.snp.bed.gz",
           snp_bed_gz_tbi = "analysis/snps/{sample}.snp.bed.gz.tbi"
        envmodules:
           config["envmodules"]["biscuit"],
        benchmark:
            "benchmarks/biscuit_snps/{sample}.txt"        
        params:
           reference = newRef,
           snp_bed = "analysis/snps/{sample}.snp.bed"
        resources:
           mem_gb = config["hpcParameters"]["intermediateMemoryGb"],
           walltime = config["walltime"]["medium"]
        log:
           epiread = "logs/snps/snps.{sample}.log",
        shell:
           """
           biscuit vcf2bed -t snp {input.vcf_gz} > {params.snp_bed}
           bgzip {params.snp_bed}
           tabix -p bed {output.snp_bed_gz}
           """

if config["epiread"]:
    rule biscuit_epiread:
        input: 
            bam = "analysis/align/{sample}.sorted.markdup.bam",
            snps = "analysis/snps/{sample}.snp.bed.gz",
            snps_tbi = "analysis/snps/{sample}.snp.bed.gz.tbi",
        envmodules:
           config["envmodules"]["biscuit"],
           config["envmodules"]["htslib"],
        params:
            reference = newRef,
            epibed = "analysis/epiread/{sample}.epibed",
            nome=config["is_nome"]
        benchmark:
            "benchmarks/biscuit_epiread/{sample}.txt"    
        resources:
            mem_gb = config["hpcParameters"]["intermediateMemoryGb"],
            walltime = config["walltime"]["medium"]
        output:
            epibed_gz = "analysis/epiread/{sample}.epibed.gz",
            epibed_gz_tbi = "analysis/epiread/{sample}.epibed.gz.tbi",
        log:
           epiread = "logs/epiread/epiread.{sample}.log",
        shell:
           """
           if [[ "$(zcat {input.snps} | head -n 1 | wc -l)" == "1" ]]; then
               if [ {params.nome} -eq 1 ]; then
                   biscuit epiread -N -B <(zcat {input.snps}) {params.reference} {input.bam} | sort -k1,1 -k2,2n > {params.epibed}
               else
                   biscuit epiread -B <(zcat {input.snps}) {params.reference} {input.bam} | sort -k1,1 -k2,2n > {params.epibed}
               fi
           else
               if [ {params.nome} -eq 1 ]; then
                   biscuit epiread -N {params.reference} {input.bam} | sort -k1,1 -k2,2n > {params.epibed}
               else
                   biscuit epiread {params.reference} {input.bam} | sort -k1,1 -k2,2n > {params.epibed}
               fi
           fi
           bgzip {params.epibed}
           tabix -p bed {output.epibed_gz}
           """
