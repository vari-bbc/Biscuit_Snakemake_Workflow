rename_fastqs <- function(samplesheet, fastq_dir, samp) {
    sheet <- read.delim(samplesheet, sep="\t")
    
    print(paste("Renaming PE reads for:", samp))
    
    # Rename R1 reads
    file_1 <- paste0(fastq_dir, "/", unlist(strsplit(subset(sheet, sample == samp)$fq1, split = ",")))
    print(paste("Found", length(file_1), "R1 files for", samp))
    
    for (i in 1:length(file_1)) {
        if (file.exists(file_1[i])) {
            new_name <- paste0(fastq_dir, "/", samp, "-", i, "-R1.fastq.gz")
            system(paste0("ln -sr ", file_1[i], " ", new_name))
            print(paste(file_1[i], "successfully symlinked to", new_name))
        } else {
            stop(paste("Error:", file_1[i], "not found in", fastq_dir))
        }
    }
    
    # Rename R2 reads
    file_2 <- paste0(fastq_dir, "/", unlist(strsplit(subset(sheet, sample == samp)$fq2, split = ",")))
    print(paste("Found", length(file_2), "R2 files for", samp))
    
    for (i in 1:length(file_2)) {
        if (file.exists(file_2[i])) {
            new_name <- paste0(fastq_dir, "/", samp, "-", i, "-R2.fastq.gz")
            system(paste0("ln -sr ", file_2[i], " ", new_name))
            print(paste(file_2[i], "successfully symlinked to", new_name))
        } else {
            stop(paste("Error:", file_2[i], "not found in", fastq_dir))
        }
    }
}

rename_fastqs(snakemake@params[["samplesheet"]], snakemake@params[["fastq_dir"]], snakemake@params[["sample"]])
