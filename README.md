This Snakemake workflow takes paired-end whole-genome bisulfite sequencing (WGBS) data and processes it using
BISulfite-seq CUI Toolkit (BISCUIT).

BISCUIT was written to perform alignment, DNA methylation and mutation calling, and allele specific methylation from
bisulfite sequencing data (https://huishenlab.github.io/biscuit/).

Download BISCUIT here: https://github.com/huishenlab/biscuit/releases/latest.

# Components of the workflow
  0. [default off] Modify and index genome reference to including methylation controls (currently not implemented)
  1. [default off] Trim adapters and/or hard clip R2
  2. [default off] Run Fastq Screen in bisulfite mode
  3. Run FastQC on raw FASTQ files
  4. Alignment, duplicate tagging, indexing, flagstat of input data (biscuitBlaster v1 and v2)
  5. Methylation information extraction (BED Format)
  6. Merge C and G beta values in CpG dinucleotide context
  7. [default off] SNP and Epiread extraction
  8. [default off] Run Preseq on aligned BAM
  9. MultiQC with BICUIT QC modules specifically for methyaltion data
  10. [default off] QC methylated and unmethylated controls (currently not implemented)

Many options can be easily specified in the `config.yaml`! Otherwise, the commands in the Snakefile can also be modified
to meet different needs.

# Dependencies
  + BISCUIT
  + R with packages tidyverse, patchwork, and viridis (only required for plotting methylation controls)
  + SAMTools 1.12+
  + Snakemake 6.0+
  + samblaster
  + htslib
  + bedtools
  + pigz
  + parallel
  + bismark (only required if running fastq_screen)
  + fastq_screen (only required if running fastq_screen)
  + fastQC
  + multiQC
  + parallel
  + preseq

# Running the workflow
+ Clone the repo `git clone https://github.com/vari-bbc/WGBS_Biscuit_Snakemake`.

+ Place *gzipped* FASTQ files into `raw_data/`. Alternatively, you can specify the location of your *gzipped* FASTQ
files in `config/config.yaml`.

+ Replace the example `config/samples.tsv` with your own sample sheet containing:
  + A row for each sample
  + The following three columns for each row:
    + A. `sample`
    + B. `fq1` (name of R1 file for `sample` in your raw data directory)
    + C. `fq2` (name of R2 file for `sample` in your raw data directory)
    + D. Any other columns included are ignored
Note, you can either edit `config/samples.tsv` in place or specify the path to your sample sheet in
`config/config.yaml`. If you create your own sample sheet, make sure to include the header line as is seen in the
example file.

+ Modify the config.yaml to specify the appropriate
  + Reference genome
  + Biscuit index
  + Biscuit QC assets (https://github.com/huishenlab/biscuit/releases/latest)
  + Environmental modules (If modules are not available, snakemake gives a warning but will run successfully
  *as long as the required executables are in the path*)
  + Toggle optional workflow components
  + Set other run parameters in `config/config.yaml`
  + Turn on optional rules in `config/config.yaml` (change from False to True)

+ Then submit the rest workflow to an HPC using something similar to bin/run_snakemake_workflow.sh (e.g.,
`qsub -q [queue_name] bin/run_snakemake_workflow.sh`)

# After the workflow

+ The output files in `analysis/pileup/` may be imported into a `BSseq` object using `bicuiteer::readBiscuit()`.
+ `analysis/multiqc/multiqc_report.html` contains the methylation-specific BISCUIT QC modules
(https://huishenlab.github.io/biscuit/docs/alignment/QC.html)

# Test dataset

This workflow comes with a working example dataset. To test the smakemake workflow on your system, place the 10
FASTQ files in `bin/working_example_dataset` into `raw_data/` and use the default `config/samples.tsv` sample sheet.
These example files can be mapped to the human genome.

# Example workflow - 1 sample
![workflow diagram](bin/DAGs/dag.png)


# Helpful snakemake commands for debugging a workflow

Do a test run: `snakemake -npr`

Unlock directory after a manually aborted run: `snakemake --unlock --cores 1`

Create a workflow diagram for your run: `snakemake --dag | dot -Tpng > my_dag.png`

Run pipeline from the command line: `snakemake --use-envmodules --cores 1`

For more information on Snakemake: https://snakemake.readthedocs.io/en/stable/

