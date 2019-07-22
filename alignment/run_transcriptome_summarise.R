# OPTPARSE ----------------------------------------------------------------

#!/usr/bin/env Rscript
require(optparse)

option_list = list(
  make_option(c("--inputpath"), type="character", help="Full path to bwa count tables for each sample", metavar="character"),
  make_option(c("--outputpath"), type="character", help="Full path to summarised output count tables per protocol/mapper/annotation setting", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

# SETTINGS ----------------------------------------------------------------

library(methods)
library(tidyverse)
library(data.table)

UMIProtocols <- c("CELseq2",  "Dropseq",  "MARSseq",  "SCRBseq")
NoUMIProtocols <- c("Smartseq2", "SmartseqC1")

# INPUT -------------------------------------------------------------------s

# opt=list(annotation='gencode',
# protocol='CELseq2',
# inputpath='counts/bwa/CELseq2/gencode')

if(isTRUE(grepl(pattern = paste(UMIProtocols,collapse="|"), opt$inputpath))) {
  
  #read in  count files
  files <- list.files(pattern = "counts.rds",opt$inputpath, full.names = T)
  filenames <- list.files(pattern = "counts.rds",opt$inputpath, full.names = F)
  filenames <- substr(filenames, start=1, stop = nchar(filenames)-11)
  names(files) <- filenames
  
  geneumis.L <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$GeneUmis
    names(tmp2) <- tmp$gene
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  geneumis <- do.call('cbind', geneumis.L)
  geneumis <- geneumis[!rownames(geneumis) %in% "*", ]
  
  transcriptumis.L  <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$TranscriptUmis
    names(tmp2) <- tmp$transcript
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  transcriptumis <- do.call('cbind', transcriptumis.L)
  transcriptumis <- transcriptumis[!rownames(transcriptumis) %in% "*", ]
  
  genecounts.L <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$GeneReads
    names(tmp2) <- tmp$gene
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  genecounts <- do.call('cbind', genecounts.L)
  genecounts <- genecounts[!rownames(genecounts) %in% "*", ]
  
  transcriptcounts.L  <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$TranscriptReads
    names(tmp2) <- tmp$transcript
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  transcriptcounts <- do.call('cbind', transcriptcounts.L)
  transcriptcounts <- transcriptcounts[!rownames(transcriptcounts) %in% "*", ]
  
  # read in filt counts tables
  files <- list.files(pattern = "counts.filt.rds",opt$inputpath, full.names = T)
  filenames <- list.files(pattern = "counts.filt.rds",opt$inputpath, full.names = F)
  filenames <- substr(filenames, start=1, stop = nchar(filenames)-16)
  names(files) <- filenames
  
  geneumis.L <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$GeneUmis
    names(tmp2) <- tmp$gene
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  filtgeneumis <- do.call('cbind', geneumis.L)
  filtgeneumis <- filtgeneumis[!rownames(filtgeneumis) %in% "*", ]
  
  transcriptumis.L  <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$TranscriptUmis
    names(tmp2) <- tmp$transcript
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  filttranscriptumis <- do.call('cbind', transcriptumis.L)
  filttranscriptumis <- filttranscriptumis[!rownames(filttranscriptumis) %in% "*", ]
  
  genecounts.L <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$GeneReads
    names(tmp2) <- tmp$gene
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  filtgenecounts <- do.call('cbind', genecounts.L)
  filtgenecounts <- filtgenecounts[!rownames(filtgenecounts) %in% "*", ]
  
  transcriptcounts.L  <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$TranscriptReads
    names(tmp2) <- tmp$transcript
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  filttranscriptcounts <- do.call('cbind', transcriptcounts.L)
  filttranscriptcounts <- filttranscriptcounts[!rownames(filttranscriptcounts) %in% "*", ]
  
  #read in  tbl files
  files <- list.files(pattern = ".tbl.rds", opt$inputpath, full.names = T)
  filenames <- list.files(pattern = ".tbl.rds",opt$inputpath, full.names = F)
  filenames <- substr(filenames, start=1, stop = nchar(filenames)-8)
  names(files) <- filenames
  
  readstats.L <- sapply(files, function(f) {
    print(names(f))
    tmp <- readRDS(file=f)
    unmapped <- sum(grepl(tmp$MapHitGene, pattern="Unmapped"))
    uniqmapped <- sum(grepl(tmp$MapHitGene, pattern="UniquelyMapped"))
    repmapped <- sum(grepl(tmp$MapHitGene, pattern="RepeatedlyMapped"))
    aligned <- uniqmapped + repmapped
    assigned <- sum(tmp$MapHitGene=="UniquelyMapped")
    tmp2 <- c(unmapped, uniqmapped, repmapped, aligned, assigned)
    names(tmp2) <- c("Unmapped", "UniquelyMapped", "RepeatedlyMapped", "Aligned", "Aligned+Assigned")
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  readstats <- do.call('rbind', readstats.L)
  readstats <- readstats[rowSums(readstats)>0,]
  
  umistats.L <- sapply(files, function(f) {
    print(names(f))
    tmp <- readRDS(file=f)
    tbl <- tmp %>% 
      dplyr::group_by(xm) %>% 
      dplyr::count(MapHitGene) %>% 
      tidyr::spread(MapHitGene, n, fill = 0) 
    unmapped <- sum(tbl$Unmapped)
    uniqmapped <- sum(tbl$UniquelyMapped)
    repmapped <- sum(tbl$RepeatedlyMapped)
    aligned <- uniqmapped + repmapped
    assigned <- uniqmapped
    tmp2 <- c(unmapped, uniqmapped, repmapped, aligned, assigned)
    names(tmp2) <- c("Unmapped", "UniquelyMapped", "RepeatedlyMapped", "Aligned", "Aligned+Assigned")
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  readstats <- do.call('rbind', readstats.L)
  readstats <- readstats[rowSums(readstats)>0,]
  
  # combine the output matrices into list
  out <- list(GeneUmis=geneumis,
              FiltGeneUmis=filtgeneumis,
              TranscriptUmis=transcriptumis,
              FiltTranscriptUmis=filttranscriptumis,
              GeneReads=genecounts,
              FiltGeneReads=filtgenecounts,
              TranscriptReads=transcriptcounts,
              FiltTranscriptReads=filttranscriptcounts,
              MapStatsReads=readstats,
              MapStatsUMIs=umistats)
}

if(isTRUE(grepl(pattern = paste(NoUMIProtocols,collapse="|"), opt$inputpath))) {
  
  #read in  count files
  files <- list.files(pattern = ".counts.rds",opt$inputpath, full.names = T)
  filenames <- list.files(pattern = ".counts.rds",opt$inputpath, full.names = F)
  filenames <- substr(filenames, start=1, stop = nchar(filenames)-11)
  names(files) <- filenames
  
  genecounts.L <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$GeneReads
    names(tmp2) <- tmp$gene
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  genecounts <- do.call('cbind', genecounts.L)
  genecounts <- genecounts[!rownames(genecounts) %in% "*", ]
  
  transcriptcounts.L  <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$TranscriptReads
    names(tmp2) <- tmp$transcript
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  transcriptcounts <- do.call('cbind', transcriptcounts.L)
  transcriptcounts <- transcriptcounts[!rownames(transcriptcounts) %in% "*", ]
  
  
  # read in filt counts tables
  files <- list.files(pattern = "counts.filt.rds",opt$inputpath, full.names = T)
  filenames <- list.files(pattern = "counts.filt.rds",opt$inputpath, full.names = F)
  filenames <- substr(filenames, start=1, stop = nchar(filenames)-16)
  names(files) <- filenames
  
  genecounts.L <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$GeneReads
    names(tmp2) <- tmp$gene
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  filtgenecounts <- do.call('cbind', genecounts.L)
  filtgenecounts <- filtgenecounts[!rownames(filtgenecounts) %in% "*", ]
  
  transcriptcounts.L  <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    tmp2 <- tmp$TranscriptReads
    names(tmp2) <- tmp$transcript
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  filttranscriptcounts <- do.call('cbind', transcriptcounts.L)
  filttranscriptcounts <- filttranscriptcounts[!rownames(filttranscriptcounts) %in% "*", ]
  
  #read in  tbl files
  files <- list.files(pattern = ".tbl.rds",opt$inputpath, full.names = T)
  filenames <- list.files(pattern = ".tbl.rds",opt$inputpath, full.names = F)
  filenames <- substr(filenames, start=1, stop = nchar(filenames)-8)
  names(files) <- filenames
  
  readstats.L <- sapply(files, function(f) {
    tmp <- readRDS(file=f)
    unmapped <- sum(grepl(tmp$MapHitGene, pattern="Unmapped"))
    uniqmapped <- sum(grepl(tmp$MapHitGene, pattern="UniquelyMapped"))
    repmapped <- sum(grepl(tmp$MapHitGene, pattern="RepeatedlyMapped"))
    aligned <- uniqmapped + repmapped
    assigned <- sum(tmp$MapHitGene=="UniquelyMapped")
    tmp2 <- c(unmapped, uniqmapped, repmapped, aligned, assigned)
    names(tmp2) <- c("Unmapped", "UniquelyMapped", "RepeatedlyMapped", "Aligned", "Aligned+Assigned")
    tmp2
  }, USE.NAMES = TRUE, simplify = FALSE)
  readstats <- do.call('rbind', readstats.L)
  readstats <- readstats[rowSums(readstats)>0,]
  
  # combine the output matrices into list
  out <- list(GeneReads=genecounts,
              FiltGeneReads=filtgenecounts,
              TranscriptReads=transcriptcounts,
              FiltTranscriptReads=filttranscriptcounts,
              MapStatsReads=readstats)
  
}


# OUTPUT ------------------------------------------------------------------

saveRDS(out,file = paste0(opt$outputpath,"/", "tbls.rds"))

# FIN ---------------------------------------------------------------------

gc()

q()
