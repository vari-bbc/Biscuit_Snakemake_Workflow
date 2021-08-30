This Snakemake workflow takes paired-end whole-genome bisulfite sequencing (WGBS) data and processes it using BISulfite-seq CUI Toolkit (BISCUIT).

BISCUIT was written to perform alignment, DNA methylation and mutation calling, and allele specific methylation from bisulfite sequencing data (https://huishenlab.github.io/biscuit/).

Download BISCUIT here: https://github.com/huishenlab/biscuit/releases/latest.

# Components of the workflow
	1. Trim adapters
	2. Alignment, duplicate tagging, indexing, flagstat 
	3. Methylation information extraction (BED Format)
	4. Merge C annd G beta values in CpG dinucleotide context
	4. MultiQC with BICUIT QC modules specifically for methyaltion data
	5. [optional] QC methylated and unmethylated controls
	6. [optional] Build and index reference


# Running the workflow

+ Clone the repo `git clone https://github.com/vari-bbc/WGBS_Biscuit_Snakemake`


+ Place *gzipped* FASTQ files into `raw_data/`


+ Replace the example `bin/samples.tsv` with your own sample sheet containing:
	+ A row for each sample
	+ The following three columns
		A. `sample_1`
		B. `fq1` (name of R1 file for `sample_1` in `raw_data/`)
		C. `fq2` (name of R2 file for `sample_1` in `raw_data/`)
		D. Any other columns included are ignored
		
		
+ Modify the config.yaml to specify the appropriate 
	+ Reference genome 
	+ Biscuit index
	+ Biscuit QC assets (https://github.com/huishenlab/biscuit/releases/tag/v0.3.16.20200420)
	+ Environmental modules (If you do not have these modules, the executables for the required programs are in the path. It will give a warning but run)
	+ Additional parameters


+ Submit the workflow to an HPC using something similar to bin/run_snake.sh (e.g. qsub -q [queue_name] bin/run_snake.sh)

# After the workflow

+ The output files in analysis/pileup/ may be imported into a `BSseq` object using `bicuiteer::readBiscuit()`.
+ analysis/multiqc/multiqc_report.html contains the methylation-specific BISCUIT QC modules (https://huishenlab.github.io/biscuit/docs/alignment/QC.html)

# Diagram of workflow

![workflow diagram](bin/DAG.png)

# Helpful snakemake commands for debugging a workflow

snakemake -npR # test run

`snakemake --unlock --cores 1` # unlock after a manually aborted run

`snakemake --dag | dot -Tpng > my_dag.png` # create a workflow diagram for your run

`snakemake --use-envmodules --cores 1` # if running on the command line, need use-envmodules option

For more information on Snakemake: https://snakemake.readthedocs.io/en/stable/

