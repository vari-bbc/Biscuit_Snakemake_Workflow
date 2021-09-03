#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

library(tidyverse)
units <- read_tsv("bin/samples.tsv")

for(samp in units$sample){
  print(paste("Renaming PE reads for:", samp))
  file_1 <- paste0("raw_data/",(units %>% dplyr::filter(sample==samp))$fq1)
  if (file.exists(file_1)){
    system(
      paste0("ln -sr ",toString(paste0("raw_data/",file_1))," raw_data/", samp, "-R1.fastq.gz")
    )
    print("R1 fq2 renamed")
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq1 file for", samp,'(',file_1,')', "listed in bin/samples.tsv not present in raw_reads/."))
  }
  # merge R2 reads
  file_2 <- paste0("raw_data/",(units %>% dplyr::filter(sample==samp))$fq2)
  if (file.exists(file_1)){
    system(
      paste0("ln -sr ",toString(paste0("raw_data/",file_1))," raw_data/", samp, "-R2.fastq.gz")
    )
    print("R2 fq2 renamed")
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq2 file for", samp,'(',file_2,')', "listed in bin/samples.tsv not present in raw_reads/."))
  }
}