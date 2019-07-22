# OPTPARSE ----------------------------------------------------------------

#!/usr/bin/env Rscript
require(optparse)

option_list = list(
  make_option(c("--samplename"), type="character", help="Sample name", metavar="character"),
  make_option(c("--bamfile"), type="character", help="Transcriptome mapped reads in read-sorted bam file with full path", metavar="character"),
  make_option(c("--GeneTranscript"), type="character", help="gene-transcript conversion table with full path", metavar="character"),
  make_option(c("--outputdir"), type="character", help="Path to output read/umi count rds file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

# SETTINGS ----------------------------------------------------------------

library(methods)
library(tidyverse)
library(data.table)
library(Rsamtools)

# INPUT -------------------------------------------------------------------

# opt = list(bamfile="mapping/bwa/Smartseq2/gencode/SmartSeq2A.C2.bam",
# samplename="SmartSeq2A.C2",
# GeneTranscript="annotation/mm10/gencode/transcript_gene_conv_spike_gencode.txt",
# outputdir="counts/bwa/Smartseq2/gencode")

trans_table=read.table(opt$GeneTranscript, 
                       col.names=c('gene', 'transcript'), stringsAsFactors = F, header=F)
trans_table$GeneTranscript=paste(trans_table$gene, trans_table$transcript, sep="_")
trans_table=rbind(c("*","*", "*_*"),trans_table)

orig.bam = Rsamtools::scanBam(opt$bamfile, param = Rsamtools::ScanBamParam(what=c("qname")))
orig.readid = orig.bam[[1]]$qname
sbamfile = substr(opt$bamfile, start = 1, stop = nchar(opt$bamfile)-4)
sbam = Rsamtools::sortBam(file = opt$bamfile, destination = paste0(sbamfile, ".sort"))
mbam = Rsamtools::scanBam(sbam, param = Rsamtools::ScanBamParam(what=c("qname", "rname", "mapq")))
fbam = Rsamtools::scanBam(sbam, param = Rsamtools::ScanBamParam(what="flag", tag=c("XT", "XA")))

fn <- paste0(sbamfile, ".sort.bam")
if (file.exists(fn)) {file.remove(fn)}

# quantification
transcript = as.character(mbam[[1]][["rname"]])
transcript = ifelse(is.na(transcript), "*", transcript)
#gather all data in one tibble
reads = tibble::tibble(id = mbam[[1]][["qname"]],
                       transcript = transcript,
                       gene = trans_table$gene[match(x = transcript,table = trans_table$transcript)],
                       AS = ifelse(!is.na(fbam[[1]][["tag"]]$XT), paste0("XT:A:", fbam[[1]][["tag"]]$XT), NA),
                       SAS = mbam[[1]][["mapq"]],
                       AH = ifelse(!is.na(fbam[[1]][["tag"]]$XA), paste0("XA:Z:", fbam[[1]][["tag"]]$XA), NA))

tbl = reads 
tbl$AS = ifelse(grepl(pattern="XT", tbl$AS), tbl$AS, NA)
tbl$AH = ifelse(grepl(pattern="XA:Z:", tbl$AH), tbl$AH, NA)
tbl$AH = gsub(pattern = "XA:Z:", replacement = "", x = tbl$AH)

head(tbl)

tmp <- strsplit(tbl$AH[!is.na(tbl$AH)], ";")
tmp2 <- sapply(tmp, function(i) {sapply(strsplit(i, ","), "[", 1)})
tmp3 <- lapply(1:length(tmp2), function(i) {sapply(1:length(tmp2[[i]]), function(j){
  tid <- tmp2[[i]][j]
  trans_table[trans_table$transcript == tid, "gene"]
})})

tbl.red <- tbl %>%
  dplyr::filter(!is.na(AH))
tmp4 <- sapply(1:length(tmp3), function(i) {
  tmp.res <- c(tbl.red$gene[i], tmp3[[i]])
  tmp.res2 <- unique(tmp.res)
  res <- ifelse(length(tmp.res2)==1, "UniquelyMapped", "RepeatedlyMapped")
  res
})

tmpres <- data.frame(id = tbl$id[!is.na(tbl$AH)], HitGene = tmp4, stringsAsFactors = F)
tbl <- tbl %>%
  dplyr::left_join(tmpres, by="id") %>%
  dplyr::mutate(bwaMapHitGene = recode(AS, 
                                       "XT:A:U" = "UniquelyMapped", 
                                       "XT:A:R" = "RepeatedlyMapped",
                                       .missing = "Unmapped")) %>%
  dplyr::mutate(CompareHit = ifelse(c(bwaMapHitGene == "RepeatedlyMapped" & 
                                        HitGene == "UniquelyMapped"), 
                                    "UniquelyMapped", NA)) %>% 
  dplyr::mutate(MapHitGene=ifelse(is.na(CompareHit), bwaMapHitGene, CompareHit)) 

# tables for counting
tbl4counts <- tbl
tbl4counts.filt <- tbl4counts %>%
  dplyr::filter(MapHitGene %in% c("UniquelyMapped", "Unmapped"))   


# count table
cnts= tbl4counts %>% 
  tidyr::unite("GeneTranscript", c('gene', 'transcript'), remove=FALSE) %>%
  dplyr::group_by(GeneTranscript, transcript, gene) %>%
  dplyr::summarise(TranscriptReads=length(transcript)) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(GeneReads=sum(TranscriptReads))

counts_full = dplyr::left_join(trans_table, cnts, by="GeneTranscript") %>% 
  dplyr::select(-contains("y")) %>%
  tidyr::replace_na(list(TranscriptReads = 0, GeneReads=0)) %>%
  dplyr::rename("gene"="gene.x", "transcript"="transcript.x")

cnts.filt = tbl4counts.filt %>% 
  tidyr::unite("GeneTranscript", c('gene', 'transcript'), remove=FALSE) %>%
  dplyr::group_by(GeneTranscript, transcript, gene) %>% 
  dplyr::summarise(TranscriptReads=length(transcript)) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(GeneReads=sum(TranscriptReads))

counts_full.filt = dplyr::left_join(trans_table, cnts.filt, by="GeneTranscript") %>% 
  dplyr::select(-dplyr::contains("y")) %>%
  tidyr::replace_na(list(TranscriptReads = 0, GeneReads=0)) %>%
  dplyr::rename("gene"="gene.x", "transcript"="transcript.x")

# OUTPUT ------------------------------------------------------------------

saveRDS(tbl,file = paste0(opt$outputdir,"/", opt$samplename,".tbl.rds"))
saveRDS(counts_full, file = paste0(opt$outputdir,"/", opt$samplename,".counts.rds"))
saveRDS(counts_full.filt, file = paste0(opt$outputdir,"/", opt$samplename,".counts.filt.rds"))

# FIN ---------------------------------------------------------------------

gc()

q()

 
