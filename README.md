# WGBS workflow

This workflow takes raw, paired-end whole-genome bisultite sequencing (WGBS) data and extracts methylation information and does bisulfite-specific QC.
The workflow centers around Biscuit, a tool for mapping and analyzing WGBS data: https://huishenlab.github.io/biscuit/.

## Running the workflow
+ Clone the repo to HPC
		git clone https://github.com/vari-bbc/WGBS_workflow
+ Load your gzipped reads into `raw_data/`
+ Setup a sample sheet in `src/samples.tsv` containing:
	+ a row for each sample
	+ three columns (feel free to include any other columns as needed, these are ignored):
		1. `sample`
		2. `fq1` (name of R1 file for sample in `raw_data/`)
		3. `fq2` (name of R2 file for sample in `raw_data/`)
+ Modify the config.yaml to specify the appropriate reference, biscuit index, biscuit QC assets, env modules & other parameters as needed
+ Submit the workflow to an HPC using something similar to src/run_snake.sh (will need to be changed)
		e.g. qsub -q [queue_name] src/run_snake.sh

## After the workflow
+ the following output files may be imported into a `BSseq` object using `bicuiteer::readBiscuit()`.
		analysis/pileup/combined.vcf.gz
		analysis/pileup/combined.bed.gz
+ proceed to differential methylation analysis with e.g. `dmrseq`

## Diagram of workflow
![workflow diagram](src/DAG.png)

## Helpful snakemake commands for debugging a workflow

This workflow currently depends on modules. If you do not have these modules, it will still run with a warning, but you need to have the executables for the dependencies in the path (all modules can be configered in config.yaml)

snakemake -np # test run

`snakemake --unlock --cores 1` # unlock after a manually aborted run

`snakemake --dag | dot -Tpng > my_dag.png` # create a workflow diagram for your run

`snakemake --use-envmodules --cores 1` # if running on the command line, need use-envmodules option

To test the workflow download (in the workflow dir) the hct116 WGBS data (you will need fastq-dump & seqtk): 

`fastq-dump -A SRR1536575 --split-files --gzip`. 

Then split it into 5 samples:
for i in {1..5}; do
   seqtk sample -s${i} SRR1536575_1.fastq.gz 1000 > raw_data/sample_${i}_R1.fq
   seqtk sample -s${i} SRR1536575_2.fastq.gz 1000 > raw_data/sample_${i}_R2.fq
done
