###-----------------------------------------------------------------------------------------------------------------###
# Notes on where variables are defined, if not defined in preprocessing.smk
#
# output_directory - workflow/Snakefile
# config           - workflow/Snakefile
#
###-----------------------------------------------------------------------------------------------------------------###

checkpoint rename_fastq_files:
    output:
        symlink_dir = directory(f'{output_directory}/analysis/renamed_fastqs'),
    params:
        samplesheet = config['samples'],
        fastq_dir = config['fastqs'],
    log:
        f'{output_directory}/logs/rename/rename.log',
    threads: 1,
    resources:
        mem_gb = config['hpcParameters']['smallMemoryGb'],
        walltime = config['walltime']['short'],
    envmodules:
        config['envmodules']['snakemake'],
        config['envmodules']['R'],
    script:
        '../scripts/rename.R'

def get_renamed_fastq_files(wildcards):
    cp_output = checkpoints.rename_fastq_files.get().output.symlink_dir

    # R1 and R2 will have the same id values
    IDX, = glob_wildcards(cp_output + '/' + wildcards.sample + '-{id}-R1.fastq.gz')
    files = list(
        expand(cp_output + '/' + wildcards.sample + '-{idx}-R{read}.fastq.gz', idx = IDX, read = [1, 2])
    )
    return files

rule raw_fastqc:
    input:
        get_renamed_fastq_files
    output:
        f'{output_directory}/analysis/raw_fastqc/{{sample}}-1-R1_fastqc.html',
        f'{output_directory}/analysis/raw_fastqc/{{sample}}-1-R1_fastqc.zip',
        f'{output_directory}/analysis/raw_fastqc/{{sample}}-1-R2_fastqc.html',
        f'{output_directory}/analysis/raw_fastqc/{{sample}}-1-R2_fastqc.zip',
    params:
        dir = f'{output_directory}/analysis/raw_fastqc',
    log:
        stdout = f'{output_directory}/logs/raw_fastqc/{{sample}}.o',
        stderr = f'{output_directory}/logs/raw_fastqc/{{sample}}.e',
    benchmark:
        f'{output_directory}/benchmarks/raw_fastqc/{{sample}}.txt',
    envmodules:
        config['envmodules']['fastqc'],
    threads: config['hpcParameters']['trimThreads'],
    resources:
        mem_gb = config['hpcParameters']['smallMemoryGb'],
        walltime = config['walltime']['medium'],
    shell:
        """
        mkdir -p {params.dir}
        fastqc --outdir {params.dir} --threads {threads} {input} 2> {log.stderr} 1> {log.stdout}
        """

rule trim_reads:
    input:
        get_renamed_fastq_files
    output:
        f'{output_directory}/analysis/trim_reads/{{sample}}-R1_val_1_merged.fq.gz',
        f'{output_directory}/analysis/trim_reads/{{sample}}-R2_val_2_merged.fq.gz',
    params:
        outdir = f'{output_directory}/analysis/trim_reads',
        quality = config['trim_galore']['quality'],
        hard_trim_R2 = config['trim_galore']['hard_trim_R2'],
    log:
        stdout = f'{output_directory}/logs/trim_reads/{{sample}}.o',
        stderr = f'{output_directory}/logs/trim_reads/{{sample}}.e',
    benchmark:
        f'{output_directory}/benchmarks/trim_reads/{{sample}}.txt',
    envmodules:
        config['envmodules']['trim_galore'],
        config['envmodules']['fastqc'],
    threads: config['hpcParameters']['trimThreads']
    resources:
        mem_gb = config['hpcParameters']['smallMemoryGb'],
        walltime = config['walltime']['medium'],
    shell:
        """
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

        # Create merged R1 and R2 FASTQs, clean up files that were merged
        cat {params.outdir}/{wildcards.sample}-*-R1_val_1.fq.gz > {params.outdir}/{wildcards.sample}-R1_val_1_merged.fq.gz
        cat {params.outdir}/{wildcards.sample}-*-R2_val_2.fq.gz > {params.outdir}/{wildcards.sample}-R2_val_2_merged.fq.gz
        rm {params.outdir}/{wildcards.sample}-*-R1_val_1.fq.gz
        rm {params.outdir}/{wildcards.sample}-*-R2_val_2.fq.gz
        """

rule fastq_screen:
    input:
        get_renamed_fastq_files
    output:
        f'{output_directory}/analysis/fastq_screen/{{sample}}-1-R1_screen.html',
        f'{output_directory}/analysis/fastq_screen/{{sample}}-1-R2_screen.html',
        f'{output_directory}/analysis/fastq_screen/{{sample}}-1-R1_screen.txt',
        f'{output_directory}/analysis/fastq_screen/{{sample}}-1-R2_screen.txt',
    params:
        conf = config['fastq_screen_conf'],
        output_dir = f'{output_directory}/analysis/fastq_screen/',
    log:
        fastq_screen = f'{output_directory}/logs/fastq_screen/{{sample}}.log',
    benchmark:
        f'{output_directory}/benchmarks/fastq_screen/{{sample}}.txt'
    threads: 8
    resources:
        nodes = 1,
        mem_gb = 64,
        walltime = config['walltime']['medium']
    envmodules: 
        config['envmodules']['fastq_screen'],
        config['envmodules']['bismark'],
    shell:
        """
        fastq_screen --bisulfite --conf {params.conf} --outdir {params.output_dir} {input} 2> {log.fastq_screen}
        """
