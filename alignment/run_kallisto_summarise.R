#!/usr/bin/env Rscript
require(optparse)

option_list = list(
  make_option(c("--inputpath"), type="character", help="Full path to kallisto count tables per sample", metavar="character"),
  make_option(c("--protocol"), type="character", help="Lib prep protocol", metavar="character"),
  make_option(c("--annotation"), type="character", help="Annotation type", metavar="character"),
  make_option(c("--outputpath"), type="character", help="Full path to summarised output count tables per protocol/mapper/annotation setting", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

# SETTINGS ----------------------------------------------------------------

library(methods)
library(tidyverse)
library(data.table)

# INPUT -------------------------------------------------------------------

# opt=list(annotation='vega', 
# protocol='HGMM1kv3', 
# inputpath='counts/kallisto/')

UMIProtocols <- c("CELseq2",  "Dropseq",  "MARSseq",  "SCRBseq", "HGMM1kv3", "PBMC1kv3")
NoUMIProtocols <- c("Smartseq2", "SmartseqC1")

if(isTRUE(grepl(pattern = paste(UMIProtocols,collapse="|"), opt$protocol))) {
  #read in  count files
  files <- list.files(pattern = "tcc.rds",opt$inputpath, full.names = T, recursive = T)
  files <- files[grepl(pattern=paste0("/", opt$protocol, "/"), files)]
  files <- files[grepl(pattern=paste0("/", opt$annotation, "/"), files)]
  
  tmpumi <- readRDS(file=files[grepl(pattern="/umi/", files)])
  geneumis <- tmpumi[!grepl(tmpumi$gene, pattern=","),]
  geneumis <- geneumis %>% 
    dplyr::select(-transcript, -ec) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise_all(list(~sum(.)))
  geneid <- geneumis$gene
  geneumis <- geneumis[,grepl(colnames(geneumis), pattern= "[.]")]
  geneumis <- as.matrix(geneumis)
  rownames(geneumis) <- geneid
  
  transcriptumis <-  tmpumi[!grepl(tmpumi$transcript, pattern=","),]
  transcriptid <- transcriptumis$transcript
  transcriptumis <- transcriptumis[,grepl(colnames(transcriptumis), pattern= "[.]")]
  transcriptumis <- as.matrix(transcriptumis)
  rownames(transcriptumis) <- transcriptid
  
  tmpread <- readRDS(file=files[grepl(pattern="/no_umi/", files)])
  genereads <- tmpread[!grepl(tmpread$gene, pattern=","),]
  genereads <- genereads %>% 
    dplyr::select(-transcript, -ec) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise_all(list(~sum(.)))
  geneid <- genereads$gene
  genereads <- genereads[,grepl(colnames(genereads), pattern= "[.]")]
  genereads <- as.matrix(genereads)
  rownames(genereads) <- geneid
  
  transcriptreads <-  tmpread[!grepl(tmpread$transcript, pattern=","),]
  transcriptid <- transcriptreads$transcript
  transcriptreads <- transcriptreads[,grepl(colnames(transcriptreads), pattern= "[.]")]
  transcriptreads <- as.matrix(transcriptreads)
  rownames(transcriptreads) <- transcriptid
  
  readstats <- tmpread %>%
    dplyr::select(-c(ec, gene)) %>%
    tidyr::gather(key = "sample", value="value", -transcript) %>%
    dplyr::mutate(status=ifelse(grepl(',', transcript), 'RepeatedlyMapped', "UniquelyMapped")) %>%
    dplyr::group_by(sample, status) %>%
    dplyr::summarise(Total=sum(value)) %>%
    dplyr::ungroup()  %>%
    tidyr::spread(status, Total) %>%
    dplyr::mutate(Unmapped=NA) %>%
    data.frame() %>%
    tibble::column_to_rownames(var="sample")
  
  umistats <- tmpumi %>%
    dplyr::select(-c(ec, gene)) %>%
    tidyr::gather(key = "sample", value="value", -transcript) %>%
    dplyr::mutate(status=ifelse(grepl(',', transcript), 'RepeatedlyMapped', "UniquelyMapped")) %>%
    dplyr::group_by(sample, status) %>%
    dplyr::summarise(Total=sum(value)) %>%
    dplyr::ungroup()  %>%
    tidyr::spread(status, Total) %>%
    dplyr::mutate(Unmapped=NA) %>%
    data.frame() %>%
    tibble::column_to_rownames(var="sample")
  
  # combine the output matrices into list
  out.L <- list(GeneUmis=geneumis,
                TranscriptUmis=transcriptumis,
                GeneReads=genereads,
                TranscriptReads=transcriptreads,
                MapStatsReads=readstats,
                MapStatsUMIs=umistats)
  
}

if(isTRUE(grepl(pattern = paste(NoUMIProtocols,collapse="|"), opt$protocol))) {
  #read in  count files
  files <- list.files(pattern = "tcc.rds",opt$inputpath, full.names = T, recursive = T)
  files <- files[grepl(pattern=paste0("/", opt$protocol, "/"), files)]
  files <- files[grepl(pattern=paste0("/", opt$annotation, "/"), files)]
  
  # gene reads
  tmpread <- readRDS(file=files[grepl(pattern="/no_umi/", files)])
  genereads <- tmpread[!grepl(tmpread$gene, pattern=","),]
  genereads <- genereads %>% 
    dplyr::select(-transcript, -ec) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise_all(funs(sum))
  geneid <- genereads$gene
  genereads <- genereads[,grepl(colnames(genereads), pattern= "SmartSeq")]
  genereads <- as.matrix(genereads)
  rownames(genereads) <- geneid
  
  # transcript reads
  transcriptreads <-  tmpread[!grepl(tmpread$transcript, pattern=","),]
  transcriptid <- transcriptreads$transcript
  transcriptreads <- transcriptreads[,grepl(colnames(transcriptreads), pattern= "SmartSeq")]
  transcriptreads <- as.matrix(transcriptreads)
  rownames(transcriptreads) <- transcriptid
  
  # read stats
  readstats <- tmpread %>%
    dplyr::select(-c(ec, gene)) %>%
    tidyr::gather(key = "sample", value="value", -transcript) %>%
    dplyr::mutate(status=ifelse(grepl(',', transcript), 'RepeatedlyMapped', "UniquelyMapped")) %>%
    dplyr::group_by(sample, status) %>%
    dplyr::summarise(Total=sum(value)) %>%
    dplyr::ungroup()  %>%
    tidyr::spread(status, Total) %>%
    dplyr::mutate(Unmapped=NA) %>%
    data.frame() %>%
    tibble::column_to_rownames(var="sample")
  
  # combine the output matrices into list
  out.L <- list(GeneReads=genereads,
                TranscriptReads=transcriptreads,
                ReadStats=readstats)
  
}

# OUTPUT ------------------------------------------------------------------

saveRDS(out.L,file = paste0(opt$outputpath,"/", "tbls.rds"))

# FIN ---------------------------------------------------------------------

gc()

q()