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
  
   ## GENE UMIS
  
  tmpumi <- readRDS(file=files[grepl(pattern="/umi/", files)])
  UniqGene <- sapply(X = 1:length(tmpumi$gene), function(i) {
    tmp <- unlist(ifelse(grepl(pattern = ",", tmpumi$gene[i]),
           strsplit(tmpumi$gene[i], split = ","),
           tmpumi$gene[i]))
    tmp2 <- ifelse(length(unique(tmp)) == 1, unique(tmp), "multi")
    tmp2
  })
  
  geneumis <- tmpumi %>% 
    dplyr::mutate(gene = UniqGene) %>% 
    dplyr::filter(!gene == "multi") %>% 
    dplyr::select(-transcript, -ec) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise_all(list(~sum(.)))
  geneid <- geneumis$gene
  geneumis <- geneumis[,grepl(colnames(geneumis), pattern= "[.]")]
  geneumis <- as.matrix(geneumis)
  rownames(geneumis) <- geneid
  
  umistats.gene <- tmpumi %>%
    dplyr::select(-c(ec, transcript)) %>%
    dplyr::mutate(gene = UniqGene) %>% 
    tidyr::gather(key = "sample", value="value", -gene) %>%
    dplyr::mutate(status=ifelse(grepl('multi', gene), 'RepeatedlyMapped', "UniquelyMapped")) %>%
    dplyr::group_by(sample, status) %>%
    dplyr::summarise(Total=sum(value)) %>%
    dplyr::ungroup()  %>%
    tidyr::spread(status, Total) %>%
    dplyr::mutate(Unmapped=NA) %>%
    data.frame() %>%
    tibble::column_to_rownames(var="sample")
  
  ## TRANSCRIPT UMIS
  
  transcriptumis <-  tmpumi
  UniqTranscript <- sapply(X = 1:length(tmpumi$transcript), function(i) {
    tmp <- unlist(ifelse(grepl(pattern = ",", tmpumi$transcript[i]),
                         strsplit(tmpumi$transcript[i], split = ","),
                         tmpumi$transcript[i]))
    tmp2 <- ifelse(length(unique(tmp)) == 1, unique(tmp), "multi")
    tmp2
  })
  
  transcriptumis <- transcriptumis %>% 
    dplyr::mutate(transcript = UniqTranscript) %>% 
    dplyr::filter(!transcript == "multi") %>% 
    dplyr::select(-gene, -ec) %>% 
    dplyr::group_by(transcript) %>% 
    dplyr::summarise_all(list(~sum(.)))
  transcriptid <- transcriptumis$transcript
  transcriptumis <- transcriptumis[,grepl(colnames(transcriptumis), pattern= "[.]")]
  transcriptumis <- as.matrix(transcriptumis)
  rownames(transcriptumis) <- transcriptid

  umistats.transcript <- tmpumi %>%
    dplyr::mutate(transcript = UniqTranscript) %>% 
    dplyr::select(-c(ec, gene)) %>%
    tidyr::gather(key = "sample", value="value", -transcript) %>%
    dplyr::mutate(status=ifelse(grepl('multi', transcript), 'RepeatedlyMapped', "UniquelyMapped")) %>%
    dplyr::group_by(sample, status) %>%
    dplyr::summarise(Total=sum(value)) %>%
    dplyr::ungroup()  %>%
    tidyr::spread(status, Total) %>%
    dplyr::mutate(Unmapped=NA) %>%
    data.frame() %>%
    tibble::column_to_rownames(var="sample")
  
  ## GENE READS
  
  tmpread <- readRDS(file=files[grepl(pattern="/no_umi/", files)])
  UniqGene <- sapply(X = 1:length(tmpread$gene), function(i) {
    tmp <- unlist(ifelse(grepl(pattern = ",", tmpread$gene[i]),
                         strsplit(tmpread$gene[i], split = ","),
                         tmpread$gene[i]))
    tmp2 <- ifelse(length(unique(tmp)) == 1, unique(tmp), "multi")
    tmp2
  })
  
  genereads <- tmpread %>% 
    dplyr::mutate(gene = UniqGene) %>% 
    dplyr::filter(!gene == "multi") %>% 
    dplyr::select(-transcript, -ec) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise_all(list(~sum(.)))
  geneid <- genereads$gene
  genereads <- genereads[,grepl(colnames(genereads), pattern= "[.]")]
  genereads <- as.matrix(genereads)
  rownames(genereads) <- geneid
  
  readstats.gene <- tmpread %>%
    dplyr::mutate(gene = UniqGene) %>% 
    dplyr::select(-c(ec, transcript)) %>%
    tidyr::gather(key = "sample", value="value", -gene) %>%
    dplyr::mutate(status=ifelse(grepl('multi', gene), 'RepeatedlyMapped', "UniquelyMapped")) %>%
    dplyr::group_by(sample, status) %>%
    dplyr::summarise(Total=sum(value)) %>%
    dplyr::ungroup()  %>%
    tidyr::spread(status, Total) %>%
    dplyr::mutate(Unmapped=NA) %>%
    data.frame() %>%
    tibble::column_to_rownames(var="sample")
  
  ## TRANSCRIPT READS
  
  UniqTranscript <- sapply(X = 1:length(tmpread$transcript), function(i) {
    tmp <- unlist(ifelse(grepl(pattern = ",", tmpread$transcript[i]),
                         strsplit(tmpread$transcript[i], split = ","),
                         tmpread$transcript[i]))
    tmp2 <- ifelse(length(unique(tmp)) == 1, unique(tmp), "multi")
    tmp2
  })
  transcriptreads <-  tmpread %>% 
    dplyr::mutate(transcript = UniqTranscript) %>% 
    dplyr::filter(!transcript == "multi") %>% 
    dplyr::select(-gene, -ec) %>% 
    dplyr::group_by(transcript) %>% 
    dplyr::summarise_all(list(~sum(.)))
  transcriptid <- transcriptreads$transcript
  transcriptreads <- transcriptreads[,grepl(colnames(transcriptreads), pattern= "[.]")]
  transcriptreads <- as.matrix(transcriptreads)
  rownames(transcriptreads) <- transcriptid
  
  readstats.transcript <- tmpread %>%
    dplyr::mutate(transcript = UniqTranscript) %>% 
    dplyr::select(-c(ec, gene)) %>%
    tidyr::gather(key = "sample", value="value", -transcript) %>%
    dplyr::mutate(status=ifelse(grepl('multi', transcript), 'RepeatedlyMapped', "UniquelyMapped")) %>%
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
                MapStatsReads=readstats.gene,
                MapStatsUMIs=umistats.gene)
  
}

if(isTRUE(grepl(pattern = paste(NoUMIProtocols,collapse="|"), opt$protocol))) {
  #read in  count files
  files <- list.files(pattern = "tcc.rds",opt$inputpath, full.names = T, recursive = T)
  files <- files[grepl(pattern=paste0("/", opt$protocol, "/"), files)]
  files <- files[grepl(pattern=paste0("/", opt$annotation, "/"), files)]
  
  ## GENE READS
  
  tmpread <- readRDS(file=files[grepl(pattern="/no_umi/", files)])
  UniqGene <- sapply(X = 1:length(tmpread$gene), function(i) {
    tmp <- unlist(ifelse(grepl(pattern = ",", tmpread$gene[i]),
                         strsplit(tmpread$gene[i], split = ","),
                         tmpread$gene[i]))
    tmp2 <- ifelse(length(unique(tmp)) == 1, unique(tmp), "multi")
    tmp2
  })
  
  genereads <- tmpread %>% 
    dplyr::mutate(gene = UniqGene) %>% 
    dplyr::filter(!gene == "multi") %>% 
    dplyr::select(-transcript, -ec) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise_all(list(~sum(.)))
  geneid <- genereads$gene
  genereads <- genereads[,grepl(colnames(genereads), pattern= "[.]")]
  genereads <- as.matrix(genereads)
  rownames(genereads) <- geneid
  
  readstats.gene <- tmpread %>%
    dplyr::mutate(gene = UniqGene) %>% 
    dplyr::select(-c(ec, transcript)) %>%
    tidyr::gather(key = "sample", value="value", -gene) %>%
    dplyr::mutate(status=ifelse(grepl('multi', gene), 'RepeatedlyMapped', "UniquelyMapped")) %>%
    dplyr::group_by(sample, status) %>%
    dplyr::summarise(Total=sum(value)) %>%
    dplyr::ungroup()  %>%
    tidyr::spread(status, Total) %>%
    dplyr::mutate(Unmapped=NA) %>%
    data.frame() %>%
    tibble::column_to_rownames(var="sample")
  
  ## TRANSCRIPT READS
  
  UniqTranscript <- sapply(X = 1:length(tmpread$transcript), function(i) {
    tmp <- unlist(ifelse(grepl(pattern = ",", tmpread$transcript[i]),
                         strsplit(tmpread$transcript[i], split = ","),
                         tmpread$transcript[i]))
    tmp2 <- ifelse(length(unique(tmp)) == 1, unique(tmp), "multi")
    tmp2
  })
  transcriptreads <-  tmpread %>% 
    dplyr::mutate(transcript = UniqTranscript) %>% 
    dplyr::filter(!transcript == "multi") %>% 
    dplyr::select(-gene, -ec) %>% 
    dplyr::group_by(transcript) %>% 
    dplyr::summarise_all(list(~sum(.)))
  transcriptid <- transcriptreads$transcript
  transcriptreads <- transcriptreads[,grepl(colnames(transcriptreads), pattern= "[.]")]
  transcriptreads <- as.matrix(transcriptreads)
  rownames(transcriptreads) <- transcriptid
  
  readstats.transcript <- tmpread %>%
    dplyr::mutate(transcript = UniqTranscript) %>% 
    dplyr::select(-c(ec, gene)) %>%
    tidyr::gather(key = "sample", value="value", -transcript) %>%
    dplyr::mutate(status=ifelse(grepl('multi', transcript), 'RepeatedlyMapped', "UniquelyMapped")) %>%
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
                ReadStats=readstats.gene)
  
}

# OUTPUT ------------------------------------------------------------------

saveRDS(out.L,file = paste0(opt$outputpath,"/", "tbls.rds"))

# FIN ---------------------------------------------------------------------

gc()

q()
