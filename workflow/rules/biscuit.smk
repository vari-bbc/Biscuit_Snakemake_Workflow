###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in biscuit.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
# BISCUIT_INDEX_FORMATS - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

if config['build_ref_with_methylation_controls']:
    rule build_ref_with_methylation_controls:
        input:
            config['ref']['fasta'],
        output:
            ref = 'merged_reference/merged.fa.gz',
            refdir = directory('merged_reference/'),
            idxfil = expand('merged_reference/merged.fa.gz.{ext}', ext=BISCUIT_INDEX_FORMATS),
        log:
            f'{output_directory}/logs/build_merged_reference_index.log',
        benchmark:
            f'{output_directory}/benchmarks/build_merged_reference_index.txt',
        threads: 2
        resources:
            mem_gb=32,
            walltime = config['walltime']['medium'],
        conda:
            '../envs/biscuit.yaml'
        envmodules:
            config['envmodules']['biscuit'],
            config['envmodules']['samtools'],
        shell:
            """
            mkdir -p {output.refdir}

            if (file {input} | grep -q "extra field"); then
                cat <(bgzip -d {input}) <(zcat bin/puc19.fa.gz) <(zcat bin/lambda.fa.gz) | bgzip > {output.ref}
            elif (file {input} | grep -q "gzip compressed data, was"); then
                cat <(zcat {input}) <(zcat bin/puc19.fa.gz) <(zcat bin/lambda.fa.gz) | bgzip > {output.ref}
            else
                cat {input} <(zcat bin/puc19.fa.gz) <(zcat bin/lambda.fa.gz) | bgzip > {output.ref}
            fi

            biscuit index {output.ref}
            samtools faidx {output.ref}
            """

def get_biscuit_reference(wildcards):
    if config['build_ref_with_methylation_controls']: # Currently not set up to be generated
        return 'merged_reference/merged.fa.gz'
    else:
        return config['ref']['fasta']

def get_rename_fastq_output_R1(wildcards):
    cp_output = checkpoints.rename_fastq_files.get().output.symlink_dir

    if config['trim_galore']['trim_before_BISCUIT']:
        return output_directory + '/analysis/trim_reads/' + wildcards.sample + '-R1_val_1_merged.fq.gz'
    else:
        IDX, = glob_wildcards(cp_output + '/' + wildcards.sample + '-{id}-R1.fastq.gz')
        files = list(expand(cp_output + '/' + wildcards.sample + '-{idx}-R1.fastq.gz', idx = IDX))
        files.sort()
        return files
        
def get_rename_fastq_output_R2(wildcards):
    cp_output = checkpoints.rename_fastq_files.get().output.symlink_dir

    if config['trim_galore']['trim_before_BISCUIT']:
        return output_directory + '/analysis/trim_reads/' + wildcards.sample + '-R2_val_2_merged.fq.gz'
    else:
        IDX, = glob_wildcards(cp_output + '/' + wildcards.sample + '-{id}-R2.fastq.gz')
        files = list(expand(cp_output + '/' + wildcards.sample + '-{idx}-R2.fastq.gz', idx = IDX))
        files.sort()
        return files

rule biscuit_blaster:
    input:
        reference = get_biscuit_reference,
        R1 = get_rename_fastq_output_R1,
        R2 = get_rename_fastq_output_R2,
    output:
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
        bai = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam.bai',
        disc = f'{output_directory}/analysis/align/{{sample}}.disc.sorted.bam' if config['biscuit']['biscuit_blaster_version'] == 'v2' else [f'{output_directory}/analysis/align/{{sample}}.no.disc'],
        disc_bai = f'{output_directory}/analysis/align/{{sample}}.disc.sorted.bam.bai' if config['biscuit']['biscuit_blaster_version'] == 'v2' else [f'{output_directory}/analysis/align/{{sample}}.no.disc.bai'],
        split = f'{output_directory}/analysis/align/{{sample}}.split.sorted.bam' if config['biscuit']['biscuit_blaster_version'] == 'v2' else [f'{output_directory}/analysis/align/{{sample}}.no.split'],
        split_bai = f'{output_directory}/analysis/align/{{sample}}.split.sorted.bam.bai' if config['biscuit']['biscuit_blaster_version'] == 'v2' else [f'{output_directory}/analysis/align/{{sample}}.no.split.bai'],
        unmapped = f'{output_directory}/analysis/align/{{sample}}.unmapped.fastq.gz' if config['biscuit']['biscuit_blaster_version'] == 'v2' else [f'{output_directory}/analysis/align/{{sample}}.no.unmapped'],
        flagstat = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam.flagstat',
    params:
        # don't include the .fa/.fasta suffix for the reference biscuit idx.
        LB = config['sam_header']['LB'],
        ID = '{sample}',
        PL = config['sam_header']['PL'],
        PU = config['sam_header']['PU'],
        SM = '{sample}',
        lib_type = config['biscuit']['lib_type'],
        disc = f'{output_directory}/analysis/align/{{sample}}.disc.sam',
        split = f'{output_directory}/analysis/align/{{sample}}.split.sam',
        unmapped = f'{output_directory}/analysis/align/{{sample}}.unmapped.fastq',
        biscuit_version = config['biscuit']['biscuit_blaster_version'],
    log:
        biscuit = f'{output_directory}/logs/biscuit/biscuit_blaster.{{sample}}.log',
        biscuit_blaster_version = f'{output_directory}/logs/biscuit/blaster_version.{{sample}}.log',
        samblaster = f'{output_directory}/logs/biscuit/samblaster.{{sample}}.log',
        samtools_view = f'{output_directory}/logs/biscuit/samtools_view.{{sample}}.log',
        samtools_sort = f'{output_directory}/logs/biscuit/samtools_sort.{{sample}}.log',
        samtools_index = f'{output_directory}/logs/biscuit/samtools_index.{{sample}}.log',
        samtools_flagstat = f'{output_directory}/logs/biscuit/samtools_flagstat.{{sample}}.log',
        sort_disc =  f'{output_directory}/logs/biscuit/sort_disc.{{sample}}.log',
        index_disc = f'{output_directory}/logs/biscuit/index_disc.{{sample}}.log',    
        sort_split = f'{output_directory}/logs/biscuit/sort_split.{{sample}}.log',
        index_split = f'{output_directory}/logs/biscuit/index_split.{{sample}}.log',
        bgzip_unmapped = f'{output_directory}/logs/biscuit/bgzip_unmapped.{{sample}}.log',       
    benchmark:
        f'{output_directory}/benchmarks/biscuit_blaster/{{sample}}.txt'
    threads: config['hpcParameters']['maxThreads']
    resources:
        mem_gb = config['hpcParameters']['maxMemoryGb'],
        walltime = config['walltime']['long'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['samtools'],
        config['envmodules']['samblaster'],
        config['envmodules']['htslib'],
    shell:
        """
        if [ {params.biscuit_version} == "v2" ]; then
            echo "biscuit blaster v2" 2> {log.biscuit_blaster_version}

            # biscuitBlaster pipeline
            # The out.bam##idx##out.bam.bai is the required format for using --write-index in samtools sort if you want
            # to get a bai index, rather than a csi index
            # NOTE: THIS ONLY WORKS WITH SAMTOOLS 1.12+
            biscuit align -@ {threads} -b {params.lib_type} \
                -R '@RG\tLB:{params.LB}\tID:{params.ID}\tPL:{params.PL}\tPU:{params.PU}\tSM:{params.SM}' \
                {input.reference} <(zcat {input.R1}) <(zcat {input.R2}) 2> {log.biscuit} | \
            samblaster --addMateTags -d {params.disc} -s {params.split} -u {params.unmapped} 2> {log.samblaster} | \
            samtools sort --write-index -@ {threads} -m 5G -o {output.bam}##idx##{output.bai} -O BAM - 2> {log.samtools_sort}

            # Get some initial stats
            samtools flagstat {output.bam} 1> {output.flagstat} 2> {log.samtools_flagstat}

            # Sort, compress, and index discordant read file
            samtools sort --write-index -@ {threads} -o {output.disc}##idx##{output.disc_bai} -O BAM {params.disc} 2> {log.sort_disc}
            #samtools index -@ {threads} {output.disc} 2> {log.index_disc}

            # Sort, compress, and index split read file
            samtools sort --write-index -@ {threads} -o {output.split}##idx##{output.split_bai} -O BAM {params.split} 2> {log.sort_split}
            #samtools index -@ {threads} {output.split} 2> {log.index_split}

            # Compress unmapped/clipped FASTQ
            bgzip -@ {threads} {params.unmapped} 2> {log.bgzip_unmapped}

            # Clean up
            rm {params.disc}
            rm {params.split}
        elif [ {params.biscuit_version} == "v1" ]; then
            echo "biscuit blaster v1" 2> {log.biscuit_blaster_version}

            # biscuitBlaster pipeline
            # The out.bam##idx##out.bam.bai is the required format for using --write-index in samtools sort if you want
            # to get a bai index, rather than a csi index
            # NOTE: THIS ONLY WORKS WITH SAMTOOLS 1.12+
            biscuit align -@ {threads} -b {params.lib_type} \
                -R '@RG\tLB:{params.LB}\tID:{params.ID}\tPL:{params.PL}\tPU:{params.PU}\tSM:{params.SM}' \
                {input.reference} <(zcat {input.R1}) <(zcat {input.R2}) 2> {log.biscuit} | \
            samblaster --addMateTags 2> {log.samblaster} | \
            samtools sort --write-index -@ {threads} -m 5G -o {output.bam}##idx##{output.bai} -O BAM - 2> {log.samtools_sort}

            # Get some initial stats
            samtools flagstat {output.bam} 1> {output.flagstat} 2> {log.samtools_flagstat}

            # Add files to ensure smooth output transition
            touch {output.disc}
            touch {output.disc_bai}
            touch {output.split}
            touch {output.split_bai}
            touch {output.unmapped}
        else
            echo "biscuit: biscuit_blaster_version must be v1 or v2 in config/config.yaml" 2> {log.biscuit_blaster_version}
        fi
        """

rule biscuit_pileup:
    input:
        ref = get_biscuit_reference,
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
    output:
        vcf_gz = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz',
        vcf_tabix = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz.tbi',
        meth = f'{output_directory}/analysis/pileup/{{sample}}.vcf_meth_average.tsv',
        bed_gz = f'{output_directory}/analysis/pileup/{{sample}}.bed.gz',
        bed_tbi = f'{output_directory}/analysis/pileup/{{sample}}.bed.gz.tbi',
    params:
        vcf = f'{output_directory}/analysis/pileup/{{sample}}.vcf',
        bed = f'{output_directory}/analysis/pileup/{{sample}}.bed',
        nome = config['is_nome'],
    log:
        pileup = f'{output_directory}/logs/biscuit_pileup/{{sample}}.pileup.log',
        vcf_gz = f'{output_directory}/logs/biscuit_pileup/{{sample}}.vcf_gz.log',
        vcf_tbi = f'{output_directory}/logs/biscuit_pileup/{{sample}}.vcf_tbi.log',
        vcf2bed = f'{output_directory}/logs/biscuit_pileup/{{sample}}.vcf2bed.log',
        bed_gz = f'{output_directory}/logs/biscuit_pileup/{{sample}}.bed_gz.log',
        bed_tbi = f'{output_directory}/logs/biscuit_pileup/{{sample}}.bed_tabix.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_pileup/{{sample}}.txt',
    threads: config['hpcParameters']['pileupThreads']
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    wildcard_constraints:
        sample = '.*[^(_mergecg)]',
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        if [ {params.nome} == "True" ]; then
            biscuit pileup -N -@ {threads} -o {params.vcf} {input.ref} {input.bam} 2> {log.pileup}
        else
            biscuit pileup -@ {threads} -o {params.vcf} {input.ref} {input.bam} 2> {log.pileup}
        fi

        bgzip {params.vcf} 2> {log.vcf_gz}
        tabix -p vcf {output.vcf_gz} 2> {log.vcf_tbi}

        biscuit vcf2bed -t cg {output.vcf_gz} 1> {params.bed} 2> {log.vcf2bed}
        bgzip {params.bed} 2> {log.bed_gz}
        tabix -p bed {output.bed_gz} 2> {log.bed_tbi}
        """

rule biscuit_mergecg:
    input:
        ref = get_biscuit_reference,
        bed = f'{output_directory}/analysis/pileup/{{sample}}.bed.gz',
    output:
        mergecg_gz = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz',
        mergecg_tbi = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz.tbi',
    params:
        mergecg = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed',
        nome = config['is_nome']
    log:
        mergecg = f'{output_directory}/logs/biscuit_pileup/mergecg.{{sample}}.log',
        mergecg_gz = f'{output_directory}/logs/biscuit_pileup/mergecg_gz.{{sample}}.log',
        mergecg_tbi = f'{output_directory}/logs/biscuit_pileup/mergecg_tabix.{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_mergecg/{{sample}}.txt',
    threads: 8
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    wildcard_constraints:
        sample = '.*[^(_mergecg)]'
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        if [ {params.nome} == "True" ]; then
            biscuit mergecg -N {input.ref} {input.bed} 1> {params.mergecg} 2> {log.mergecg}
        else
            biscuit mergecg {input.ref} {input.bed} 1> {params.mergecg} 2> {log.mergecg}
        fi

        bgzip {params.mergecg} 2> {log.mergecg_gz}
        tabix -p bed {output.mergecg_gz} 2> {log.mergecg_tbi}
        """

rule biscuit_snps:
    input: 
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
        vcf_gz = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz',
    output:
        snp_bed_gz = f'{output_directory}/analysis/snps/{{sample}}.snp.bed.gz',
        snp_bed_gz_tbi = f'{output_directory}/analysis/snps/{{sample}}.snp.bed.gz.tbi',
    params:
        snp_bed = f'{output_directory}/analysis/snps/{{sample}}.snp.bed',
    log:
        f'{output_directory}/logs/snps/snps.{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_snps/{{sample}}.txt',
    threads: 1
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        biscuit vcf2bed -t snp {input.vcf_gz} > {params.snp_bed} 2> {log}
        bgzip {params.snp_bed}
        tabix -p bed {output.snp_bed_gz}
        """

rule biscuit_epiread:
    input: 
        ref = get_biscuit_reference,
        bam = f'{output_directory}/analysis/align/{{sample}}.sorted.markdup.bam',
        snps = f'{output_directory}/analysis/snps/{{sample}}.snp.bed.gz',
        snps_tbi = f'{output_directory}/analysis/snps/{{sample}}.snp.bed.gz.tbi',
    output:
        epibed_gz = f'{output_directory}/analysis/epiread/{{sample}}.epibed.gz',
        epibed_gz_tbi = f'{output_directory}/analysis/epiread/{{sample}}.epibed.gz.tbi',
    params:
        epibed = f'{output_directory}/analysis/epiread/{{sample}}.epibed',
        nome = config['is_nome'],
    log:
        f'{output_directory}/logs/epiread/epiread.{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/biscuit_epiread/{{sample}}.txt'
    threads: config['hpcParameters']['pileupThreads']
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuit'],
        config['envmodules']['htslib'],
    shell:
        """
        if [[ "$(zcat {input.snps} | head -n 1 | wc -l)" == "1" ]]; then
            if [ {params.nome} == "True" ]; then
                biscuit epiread -N -@ {threads} -B <(zcat {input.snps}) {input.ref} {input.bam} | sort -k1,1 -k2,2n > {params.epibed} 2> {log}
            else
                biscuit epiread -@ {threads} -B <(zcat {input.snps}) {input.ref} {input.bam} | sort -k1,1 -k2,2n > {params.epibed} 2> {log}
            fi
        else
            if [ {params.nome} == "True" ]; then
                biscuit epiread -N {input.ref} {input.bam} | sort -k1,1 -k2,2n > {params.epibed} 2> {log}
            else
                biscuit epiread {input.ref} {input.bam} | sort -k1,1 -k2,2n > {params.epibed} 2> {log}
            fi
        fi

        bgzip {params.epibed}
        tabix -p bed {output.epibed_gz}
        """
