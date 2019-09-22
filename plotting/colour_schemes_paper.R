
# COLOUR SCHEMES ----------------------------------------------------------

# http://colorschemedesigner.com/csd-3.5/#

# slurm.input.simDe.full <- readRDS(file = "/data/share/htp/powsimR/scripts/simulateDE_params_full.rds")

# PROTOCOLS ---------------------------------------------------------------

# table(slurm.input.simDe.full$Protocol)

# per method
protocols.cols <- c("CELseq2" = "#FF6200", 
                    "MARSseq" = "#BF6370",
                    "SCRBseq" = "#FFA600", 
                    "mcSCRBseq" = "#BF8D30",
                    "Smartseq2" = "#34D800",
                    "SmartseqC1" = "#00B258",
                    "HGMM1kv3" = "#72217D", 
                    "PBMC1kv3" = "#A62A58",
                    "Dropseq" = "#A63F00")


# UMI vs non-UMI
umi.cols <- c("UMI" = "#A66C00", 
              "noUMI" = "#052A6E")


# MAPPER ------------------------------------------------------------------

 # table(slurm.input.simDe.full$Mapper)

mapper.cols <- c("bwa" = "#4671D5", 
                 "bowtie2" = "#6C8CD5",
                 "kallisto" = "#3914AF", 
                 "zumi" = "#009999", 
                 "star" = "#009999", 
                 "cellranger" = "009999")

# ANNOTATION --------------------------------------------------------------

# table(slurm.input.simDe.full$Annotation)

annotation.cols <- c("ensembl" = "#589900", 
                     "gencode" = "#589900", 
                     "refseq" ="#95A300", 
                     "vega" = "#A68F00")


# PREFILTERING ------------------------------------------------------------

# table(slurm.input.simDe.full$Prefilter)

prefilter.cols <- c("none" = "#bcbcbc",
                    "FreqFilter" = "#7908AA")


# IMPUTATION --------------------------------------------------------------

# table(slurm.input.simDe.full$Impute)

imputation.cols <- c("none" = "#8c8c8c",
                     "DrImpute" = "#CE36D3", 
                     "scImpute" = "#CF60D3", 
                     "scone" = "#F13D76", 
                     "SAVER" = "#F16D97")


# PREPROCESSING -----------------------------------------------------------

# a combination of prefiltering and imputation

preprocess.cols <- c("none" = "#8c8c8c", 
                     "Filtering" = "#fb9a99",
                     "DrImpute" = "#1f78b4", 
                     "scImpute" = "#b2df8a", 
                     "scone" = "#33a02c", 
                     "Seurat" = "#e31a1c",
                     "SAVER" = "#a6cee3")


# NORMALISATION -----------------------------------------------------------

# table(slurm.input.simDe.full$Normalisation)

# RColorBrewer::brewer.pal(n = 8, name = 'Dark2')
# 
# RColorBrewer::display.brewer.pal(n = 10, name = 'Dark2')

normalisation.cols <- c("SF" = "#666666", 
                        "MR" = "#1B9E77", 
                        "TMM" = "#66A61E",
                        "PosCounts" = "#D95F02", 
                        "Census" = "#E6AB02", 
                        "Linnorm" = "#A6761D",
                        "scran" = "#E7298A", 
                        "scran w/ groups" = "#F385BD",
                        "scran w/ cluster" = "#960D53",
                        "SCnorm" = "#7570B3",
                        "SCnorm w/ cluster" = "#2A2474")

spike.cols <- c("w/o Spike-Ins" = "#C46B31",  
                "with Spike-Ins" = "#2BAB92")


# DE-TOOL -----------------------------------------------------------------

# table(slurm.input.simDe.full$DEmethod)

detool.cols <- c("limma-trend" = "#A6611A", 
                 "edgeR-zingeR" = "#018571",
                 "MAST" = "#80CDC1", 
                 "T-Test" = "#DFC27D")


# SIMULATION SETUP --------------------------------------------------------

# percentage DE
percde.cols <- c("0.00" = "ivory3",
                 "0.05" = "indianred1",
                 "0.2" = "indianred3",
                 "0.6" = "indianred4")

# DE pattern
depattern.cols <- c("Symmetric" = "violet",
                   "Asymmetric" = "violetred",
                   "Completely Asymmetric" = "violetred4")

# sample size
samplesize.cols <- c("96 vs 96" = "olivedrab3",
                     "50 vs 200" = "olivedrab1",
                     "384 vs 384" = "olivedrab4")

