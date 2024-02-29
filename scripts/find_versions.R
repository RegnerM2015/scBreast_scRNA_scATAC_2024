###################
# scATAC-seq
###################
setwd("/datastore/nextgenout5/share/labs/francolab/scBreast_scRNA_scATAC_2023")
library(utils)
dirs <- read.table("scripts/copy_cellranger_outputs.sh")
dirs <- dirs$V2
paths <- gsub("(outs).*","\\1",dirs)
paths <- paths[1:42]

store.df <- data.frame(file=paths,
                       samples_processed_together=rep("fill",length(paths)),
                       version=rep("fill",length(paths)),
                       cellranger_info=rep("fill",length(paths)))                

for ( i in paths){
  setwd(i)
  setwd("..")
  text <- readLines("_cmdline")
  store.df$cellranger_info[which(store.df$file == i)] <- text
  
  # get other info
  if(length(grep(pattern = 'T47D|MCF7|SUM149PT|HCC1143', x = i)) >0){
    n <- 9
  }else{
    n <- 7
  }

  i.new <- vapply(strsplit(i, '/'), function(x)
    paste(x[seq.int(n)], collapse='/'), character(1L))

  setwd(i.new)
  dirs <- list.dirs(recursive = F)
  
  sample.idx <- grep("ATAC",dirs)
  sample.dirs <- dirs[sample.idx]
  
  sample.idx <- grep("date",sample.dirs)
  
  if (length(sample.idx) >0){
    sample.dirs <- sample.dirs[sample.idx]
  }else{
    sample.dirs <- sample.dirs
  }

  
  if(length(sample.dirs) ==2){
    text <- paste0(sample.dirs[1],sample.dirs[2])
    store.df$samples_processed_together[which(store.df$file == i)] <- text
  }else{
    print("sample sequenced or processed by itself")
    text <- sample.dirs
    store.df$samples_processed_together[which(store.df$file == i)] <- text
  }
 
  idx <- grep("./H",dirs)
  
  dirs <- dirs[idx]
  idx <- grep("-",dirs)
  
  if (length(idx) >0){
    dirs <- dirs[-idx]
  }else{
    idx <- grep("ATAC",dirs)
    if(length(idx) >0 ){
      dirs <- dirs[-idx]
    }else{
      dirs <- dirs
    }
    
    
    idx <- grep("date",dirs)
    
    if (length(idx) >0){
      dirs <- dirs[idx]
    }else{
      dirs <- dirs
    }
  }
  
  setwd(dirs)
  setwd("./MAKE_FASTQS_CS")
  setwd("./MAKE_FASTQS")
  setwd("./BCL2FASTQ_WITH_SAMPLESHEET")
  setwd("./fork0")
  
  text <- readLines("_outs", n=3,skip=2)
  text <- text[3]
  
  store.df$version[which(store.df$file == i)] <- text
  
  print(paste0("Completed:",i))
  setwd("/datastore/nextgenout5/share/labs/francolab/scBreast_scRNA_scATAC_2023")

}
setwd("/datastore/nextgenout5/share/labs/francolab/scBreast_scRNA_scATAC_2023")
write.table(store.df,"FilePathsAndVersions-scATAC.csv",col.names = T,row.names = F,sep = ",")

###################
# scRNA-seq
###################
library(utils)
dirs <- read.table("scripts/copy_cellranger_outputs.sh")
dirs <- dirs$V2
paths <- gsub("(outs).*","\\1",dirs)
paths <- paths[43:63]

store.df <- data.frame(file=paths,
                       samples_processed_together=rep("fill",length(paths)),
                       version=rep("fill",length(paths)),
                       cellranger_info=rep("fill",length(paths)))                

for ( i in paths){
  setwd(i)
  setwd("..")
  text <- readLines("_cmdline")
  store.df$cellranger_info[which(store.df$file == i)] <- text
  
  
  # get other info
  if(length(grep(pattern = 'T47D|MCF7|SUM149PT|HCC1143', x = i)) >0){
    n <- 9
  }else{
    n <- 7
  }
  
  i.new <- vapply(strsplit(i, '/'), function(x)
    paste(x[seq.int(n)], collapse='/'), character(1L))
  
  setwd(i.new)
  dirs <- list.dirs(recursive = F)

  if(length(grep(pattern = 'SUM149PT|HCC1143', x = dirs)) > 0){
    sample.idx <- grep("Mar2020",dirs)
    sample.dirs <- dirs[sample.idx]
  }else if(length(grep(pattern = 'MCF7|T47D', x = dirs)) > 0){
    sample.idx <- grep("July2020",dirs)
    sample.dirs <- dirs[sample.idx]
  }else{
    sample.idx <- grep("RNA",dirs)
    sample.dirs <- dirs[sample.idx]
  }
  
  sample.idx <- grep("date",sample.dirs)
  
  if (length(sample.idx) >0){
    sample.dirs <- sample.dirs[sample.idx]
  }else{
    sample.dirs <- sample.dirs
  }
  
  
  if(length(sample.dirs) ==2){
    text <- paste0(sample.dirs[1],sample.dirs[2])
    store.df$samples_processed_together[which(store.df$file == i)] <- text
  }else if (length(sample.dirs) ==1){
    print("sample sequenced or processed by itself")
    text <- sample.dirs
    store.df$samples_processed_together[which(store.df$file == i)] <- text
  }else if(length(grep(pattern = '_', x = sample.dirs)) > 0){
    idx <- grep(pattern = '_', x = sample.dirs)
    if(length(idx) >0 ){
      sample.dirs <- sample.dirs[-idx]
      text <- paste0(sample.dirs[1],sample.dirs[2])
      store.df$samples_processed_together[which(store.df$file == i)] <- text
    }else{
      sample.dirs <- sample.dirs
      text <- paste0(sample.dirs[1],sample.dirs[2])
      store.df$samples_processed_together[which(store.df$file == i)] <- text
    }
  }

  idx <- grep("./H",dirs)
  
  dirs <- dirs[idx]
  idx <- grep("-",dirs)
  
  if (length(idx) >0){
    dirs <- dirs[-idx]
  }else{
    idx <- grep("RNA",dirs)
    if(length(idx) >0 ){
      dirs <- dirs[-idx]
    }else{
      dirs <- dirs
    }
    
    idx <- grep("Mar2020",dirs)
    if(length(idx) >0 ){
      dirs <- dirs[-idx]
    }else{
      dirs <- dirs
    }
    
    idx <- grep("July2020",dirs)
    if(length(idx) >0 ){
      dirs <- dirs[-idx]
    }else{
      dirs <- dirs
    }
    
    idx <- grep("GFP_Luc",dirs)
    if(length(idx) >0 ){
      dirs <- dirs[-idx]
    }else{
      dirs <- dirs
    }
    
    idx <- grep("GFPLUC",dirs)
    if(length(idx) >0 ){
      dirs <- dirs[-idx]
    }else{
      dirs <- dirs
    }
    
    idx <- grep("date",dirs)
    
    if (length(idx) >0){
      dirs <- dirs[idx]
    }else{
      dirs <- dirs
    }
  }
  
  setwd(dirs)
  setwd("./MAKE_FASTQS_CS")
  setwd("./MAKE_FASTQS")
  setwd("./BCL2FASTQ_WITH_SAMPLESHEET")
  setwd("./fork0")
  
  text <- readLines("_outs", n=3,skip=2)
  text <- text[3]
  
  store.df$version[which(store.df$file == i)] <- text
  
  print(paste0("Completed:",i))
  setwd("/datastore/nextgenout5/share/labs/francolab/scBreast_scRNA_scATAC_2023")
  
}
setwd("/datastore/nextgenout5/share/labs/francolab/scBreast_scRNA_scATAC_2023")
write.table(store.df,"FilePathsAndVersions-scRNA.csv",col.names = T,row.names = F,sep = ",")
