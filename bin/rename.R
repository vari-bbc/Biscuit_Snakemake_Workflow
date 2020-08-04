log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

.libPaths("/primary/vari/software/BBC/R/build-3.6.0/lib64/R/library")
library(tidyverse)
units <- read_tsv("bin/samples.tsv")
samp = snakemake@output[[1]] %>% gsub("raw_data/","",.) %>% gsub("(-.{2}).fastq.gz","",.)

print(paste("merging PE reads for:", samp))
  if (file.exists(toString(paste0("raw_data/",(units %>% dplyr::filter(sample==samp))$fq1)))){
    system(paste0("ln -sr ",
                gsub(toString(paste0("raw_data/",(units %>% dplyr::filter(sample==samp))$fq1)),pattern=",", replacement = ""),
                " raw_data/", samp, "-R1.fastq.gz"))
    print("R1 units merged")
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq1 file for", samp, "listed in bin/samples.tsv not present in raw_reads/."))
  }
  # merge R2 reads
  if (file.exists(toString(paste0("raw_data/",(units %>% dplyr::filter(sample==samp))$fq2)))){
    system(paste0("ln -sr ",
                gsub(toString(paste0("raw_data/",(units %>% dplyr::filter(sample==samp))$fq2)),pattern=",", replacement = ""),
                " raw_data/", samp, "-R2.fastq.gz"))
    print("R2 units merged")
    save.image(file=paste0("logs/rename/rename-",samp,".RData"))
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq2 file for", samp, "listed in bin/samples.tsv not present in raw_reads/."))
  }
