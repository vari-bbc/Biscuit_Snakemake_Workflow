#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

# library(tidyverse)
# setwd('/Volumes/researchtemp/hpctmp/ian.beddows/canary_WGBS_snakemake/Biscuit_Snakemake_Workflow')
units <- read.delim("bin/samples.tsv",sep="\t")

for(samp in units$sample){
  print(paste("Renaming PE reads for:", samp))
  file_1 <- paste0("raw_data/",unlist(strsplit(subset(units,sample==samp)$fq1,split=',')))
  
  print(cat("Found",length(file_1),"R1 files for",samp,"\n"))
  for(i in 1:length(file_1)){
    if (file.exists(file_1[i])){
      system(
      # print(
        paste0("ln -sr ",file_1[i]," raw_data/", samp,"-", i,"-R1.fastq.gz")
      )
      print(cat("R1",file_1[i]," renamed (symlinked)"))
    } else {
      stop(paste("Error in mergeLanesAndRename.R: fq1 file for", samp,'(',file_1[i],')', "listed in bin/samples.tsv not present in raw_reads/."))
    }
  }
  
  # merge R2 reads
  file_2 <- paste0("raw_data/",unlist(strsplit(subset(units,sample==samp)$fq2,split=',')))
  print(cat("Found",length(file_2),"R2 files for",samp,"\n"))
  for(i in 1:length(file_2)){
    if (file.exists(file_2[i])){
      system(
      # print(
        paste0("ln -sr ",file_2[i]," raw_data/", samp,"-", i,"-R2.fastq.gz")
      )
      print(cat("R2",file_2[i]," renamed (symlinked)"))
    } else {
      stop(paste("Error in mergeLanesAndRename.R: fq2 file for", samp,'(',file_2[i],')', "listed in bin/samples.tsv not present in raw_reads/."))
    }
  }
}
