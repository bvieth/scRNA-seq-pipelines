# OPTPARSE ----------------------------------------------------------------

#!/usr/bin/env Rscript
require(optparse)

option_list = list(
  make_option(c("--inputfilepath"), type="character", help="Kallisto output files location with full path (matrix.tsv, matrix.cells, matrix.ec)", metavar="character"),
  make_option(c("--TranscriptIds"), type="character", help="Transcript IDs as given to kallisto with full path", metavar="character"),
  make_option(c("--TranscriptGene"), type="character", help="Transcript-Gene annotation file with full path", metavar="character"),
  make_option(c("--outputdir"), type="character", help="Path to output read/umi count rds file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

# SETTINGS ----------------------------------------------------------------

library(methods)
library(tidyverse)

# INPUT -------------------------------------------------------------------

# The matrix.tsv file contains the TCC counts in (row, col, value) format. 
# Rows are equivalence-class ids and columns are cell ids
# (a cell id corresponds to the order the cell appears in the matrix.cells file). 
# For example, "3 4 10" means that 10 reads from cell 4 are compatible with the equivalence class 3. 
# The set of transcripts that are in each eq. class are given in the matrix.ec file. 
# Each line in matrix.ec corresponds to an eq.class id and lists the ids of the transcripts it contains. 
# For example, "3\tab 4 5" means that eq class 3 contains tx4 and tx5. The transcript ids correspond 
# to their order in the specified (reference) fasta file. ***Note that everything is zero indexed, i.e., 0,1,2,...***.

# opt = list(inputfilepath="mapping/kallisto/Dropseq/vega/",
#            TranscriptIds="annotation/mm10/vega/transcript_spike_ids_vega.txt",
#            TranscriptGene="annotation/vega/transcript_gene_conv_spike_vega.txt")

# TCC in long format (row=equivalence class id, col=sample id, value=number of reads)
tsv=read.table(paste0(opt$inputfilepath, "matrix.tsv"), 
               stringsAsFactors = F)
# equivalence classes
ec=read.table(paste0(opt$inputfilepath, "matrix.ec"), 
              stringsAsFactors = F)
# cell ids (sample names)
cellnames=read.table(paste0(opt$inputfilepath, "matrix.cells"), 
                     stringsAsFactors = F)
#aquire the transcriptome ids
transcriptids=read.table(opt$TranscriptIds, stringsAsFactors = F, header = F)
transcriptgeneconv=read.table(opt$TranscriptGene, stringsAsFactors = F, header = F)
transcriptgeneconv = transcriptgeneconv[match(transcriptids$V1, transcriptgeneconv$V2),]
gids = transcriptgeneconv$V1
tids=transcriptids$V1

# transform  TCC to long table
tcc=tidyr::spread(data = tsv,key = V2,value = V3, fill = 0)

#name the abundance table according to the cellnames
colnames(tcc)[2:ncol(tcc)]=as.character(cellnames$V1)

# reduce ec matrix to used ecs of tcc matrix!!!
ec.used=ec[ec$V1 %in% tcc$V1,]

##replace the indexes of the ec file with the transcript ids (ec file is 0 based so use +1)
ecsplit=strsplit(x = as.character(ec.used$V2), split = ",")

#empty vector for result collection of loop
enst_ids=NULL
ensg_ids=NULL

#loop through all the comma seperated transcript indexes of the ec and 
for(i in 1:length(ecsplit)){
  temp_vector_t=NULL
  temp_vector_g=NULL
  for(j in 1:length(ecsplit[[i]])){
    temp_id_t = tids[as.numeric(ecsplit[[i]][j])+1]
    temp_vector_t=c(temp_vector_t,temp_id_t)
    temp_id_g = gids[as.numeric(ecsplit[[i]][j])+1]
    temp_vector_g=c(temp_vector_g,temp_id_g)
  }
  comb_vector_t=paste(temp_vector_t,collapse = ",")
  comb_vector_g=paste(temp_vector_g,collapse = ',')
  enst_ids=c(enst_ids,comb_vector_t)
  ensg_ids=c(ensg_ids, comb_vector_g)
}

#replace the ec indices with the ensmust ids
tcc$transcript=enst_ids
tcc$gene=ensg_ids
colnames(tcc)[1]="ec"

tcc[1:5, 1:10]

head(tcc$gene)



# OUTPUT ------------------------------------------------------------------

saveRDS(tcc,file = paste0(opt$outputdir,"/","tcc.rds"))

# FIN ---------------------------------------------------------------------

gc()

q()
