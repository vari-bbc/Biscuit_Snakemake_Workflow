

require(tidyverse)
require(ggplot2)
require(patchwork)
require(viridis)

import_emseq <- function(dir){
  myReturnDF <- NULL
  files <- list.files(dir)
  for(f in files){
    if(file.info(f)$size>0){
      name <- gsub('\\.bed','',f)
      print(cat('File',f,name,'\n'))
      myBed <- read.delim(paste0(dir,f),sep="\t",header=FALSE)
      colnames(myBed) <- c('chr','start','end','beta','depth','info')
      myBed$sample <- rep(name,nrow(myBed))
      myReturnDF <- rbind(
        myReturnDF,
        myBed
      )
    }
  }
  return(myReturnDF)
}

# set the base directory
# baseDir <- '/Volumes/projects_secondary/shen/projects/SHEH_20210623_WGBS/dev2/WGBS_Biscuit_Snakemake/'
baseDir <- './'
posControlDir <- 'analysis/puc19/'
negControlDir <- 'analysis/lambda/'
# NEGATIVE CONTROL
dir <- paste0(baseDir,negControlDir)
# dir <- '../analysis/qc_vectors/lambda/'
lambda <- import_emseq(dir)
dim(lambda)
nSamples <- length(unique(lambda$sample))

if(exists("x")){
  x <- lambda
  
  topleft <- ggplot(x,aes(x=sample,y=depth)) +
    geom_boxplot(color='#357BA2FF') +
    theme_bw() +
    # ylim(c(0,max(x$depth,na.rm=TRUE))) +
    theme(
      axis.text.x = element_text(angle=60,hjust=1,size=12), 
      axis.text.y = element_text(size=12), 
      axis.title.y = element_text(size=25), 
      plot.title = element_text(size=25,hjust=0.5),
      plot.subtitle = element_text(size=15,hjust=0.5)
          
    ) +
    scale_color_manual() +
    ggtitle('Unmethylated',subtitle = paste('N =',nSamples,'Samples')) +
    ylab('Coverage') +
    xlab('')
  
  bottomleft <- ggplot(x,aes(x=sample,y=beta)) +
    geom_boxplot(color='#357BA2FF') +
    theme_bw() +
    ylim(c(0,1)) +
    theme(
      axis.text.x = element_text(angle=60,hjust=1,size=12), 
      axis.text.y = element_text(size=12), 
      axis.title.y = element_text(size=25), 
      plot.title = element_text(size=25,hjust=0.5),
      plot.subtitle = element_text(size=15,hjust=0.5)
          
    ) +
    scale_color_manual() +
    ylab('Beta') +
    xlab('')
} 

# POSITIVE CONTROL
dir2 <- paste0(baseDir,posControlDir)
# dir2 <- '../analysis/qc_vectors/puc19/'
puc19 <- import_emseq(dir2)
dim(puc19)
nSamples2 <- length(unique(puc19$sample))
x <- puc19

if(exists("x")){
  topright <- ggplot(x,aes(x=sample,y=depth)) +
    geom_boxplot(color='#357BA2FF') +
    theme_bw() +
    # ylim(c(0,max(x$depth,na.rm=TRUE))) +
    theme(
      axis.text.x = element_text(angle=60,hjust=1,size=12), 
      axis.text.y = element_text(size=12), 
      axis.title.y = element_text(size=27), 
      plot.title = element_text(size=27,hjust=0.5),
      plot.subtitle = element_text(size=15,hjust=0.5)
      
    ) +
    scale_color_manual() +
    ggtitle('Methylated',subtitle = paste('N =',nSamples,'Samples')) +
    xlab('') +
    ylab('')
  
  bottomright <- ggplot(x,aes(x=sample,y=beta)) +
    geom_boxplot(color='#357BA2FF') +
    theme_bw() +
    ylim(c(0,1)) +
    theme(
      axis.text.x = element_text(angle=60,hjust=1,size=12), 
      axis.text.y = element_text(size=12), 
      axis.title.y = element_text(size=25), 
      plot.title = element_text(size=25,hjust=0.5),
      plot.subtitle = element_text(size=15,hjust=0.5)
      
    ) +
    scale_color_manual() +
    ylab('') +
    xlab('')
}

if(exists("topright") & exists("topleft")){
  layout <- '
  AB
  CD
  '
  
  pdf(file=paste0(baseDir,'analysis/qc_vectors/control_vector_boxplot.pdf'),width = 7,height = 10)
  
  topleft + topright + bottomleft + bottomright + 
    patchwork::plot_annotation(tag_levels = 'A',title = 'Control Vectors') + 
    patchwork::plot_layout(design = layout)  & theme(plot.tag = element_text(face = 'bold'))  & theme(plot.title = element_text(size=27,hjust=0.5))
  
  dev.off()
}else{
  pdf(file=paste0(baseDir,'analysis/qc_vectors/control_vector_boxplot.pdf'),width = 7,height = 10)
  plot.new()
  text(x = 0.5,y=0.5,"NO PLOT CREATED.\nLIKELY SOLUTION: NO DEPTH FOR CONTROL VECTOR(S)\n IN SUPPLIED MERGECG FILES")  # first 2 numbers are xy-coordinates within [0, 1]
  dev.off()
  
}


