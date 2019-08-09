# SETTINGS ----------------------------------------------------------------

# library(powsimRDev)
library(tidyverse)
library(gridExtra)
library(grid)
library(cowplot)
library(data.table)

source(file  = "plotting/colour_schemes_paper.R")
source(file  = "plotting/renaming_paper.R")

# INPUT -------------------------------------------------------------------

mapres.dat <- readRDS(file = "MapStats_main_2019-05-27.rds")
scrbseq.locations.dat <- readRDS(file = "SCRBseq_locations_main.rds")

detectgenes.dat <- readRDS(file = "MargsDat.rds") %>% 
  dplyr::select(-Mean, -Dispersion)

fits.dat <- readRDS(file = "FitsDat.rds")
margs.dat <- readRDS(file = "MargsDat.rds") 
margs.dat <-   data.table::melt(setDT(margs.dat), 
                                id.vars = c("Name", "Protocol", "Mapper", "Annotation", "GeneID"),
                                measure.vars = c("Mean", "Dispersion", "Dropout"),
                                variable.factor = FALSE) %>% 
  data.frame()

perf.dat <- readRDS(file = "AlignerPerformance.rds")

# MAIN FIGURE -------------------------------------------------------------

##### Mapping statistics ##########

protocol_names <- c(CELseq2 = "CEL-seq2", SCRBseq = "SCRB-seq", Smartseq2 = "Smart-seq2", 
                    HGMM1kv3 = "10X - HGMM", PBMC1kv3 = "10X Genomics", `10XGenomics` = "10X - HGMM", 
                    Dropseq = "Drop-seq")
protocol.cols <- c(CELseq2 = "#FF6200", MARSseq = "#BF6370", SCRBseq = "#FFA600", 
                   mcSCRBseq = "#BF8D30", Smartseq2 = "#34D800", SmartseqC1 = "#00B258", 
                   HGMM1kv3 = "#72217D", PBMC1kv3 = "#A62A58", Dropseq = "#A63F00", 
                   `10X Genomics` = "#A62A58")
# read mapping statistics per protocol, aligner and annotation
mapres.plot <- ggplot(data = mapres.dat %>% 
                        dplyr::filter(! Protocol == "HGMM1kv3"), 
                      aes(x = Protocol, y= Proportion, 
                          fill=Protocol, group=Protocol, 
                          alpha=Feature)) +
  geom_bar(stat="identity", width=.8, position = "dodge", color="black") + 
  scale_y_continuous(labels = scales::percent, 
                     breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1)) +
  scale_x_discrete(labels = protocol_names) +
  scale_fill_manual(values = protocols.cols, labels = protocol_names) +
  scale_alpha_manual(values=c(0.1,0.6, 1)) +
  labs(y="Percentage of Total Reads", x=NULL, fill=NULL) + 
  facet_grid(Annotation ~ Mapper, 
             scales="fixed", 
             labeller = as_labeller(relabels)) + 
  theme_light() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.text = ggplot2::element_text(size=8, color='black'),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        legend.key.size = grid::unit(0.5, "lines"),
        axis.text.y = ggplot2::element_text(size=8,color='black'),
        axis.text.x = ggplot2::element_text(size=8, color='black'),
        axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
        strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
        strip.background = ggplot2::element_rect(fill="white")) +
  coord_flip()

# number of genes UMI aligns to
scrbseq.locations.plot <- ggplot2::ggplot(data = scrbseq.locations.dat %>%
                                            dplyr::filter(Annotation == "gencode"), 
                                          aes(x = genome_mapper_locations, 
                                              y = transcriptome_mapper_locations, 
                                              fill=cut(..count.., c(0, 5, 10, 50, 100, Inf))) ) +
  ggplot2::stat_bin2d(bins = 100) + 
  ggplot2::geom_abline(slope = 1) +
  ggplot2::scale_fill_hue(l = 70, h = c(0, 90)) +
  ggplot2::scale_x_continuous(limits = c(0,125)) + 
  scale_y_continuous(limits = c(0,125)) +
  ggplot2::geom_abline(slope=1, intercept=0, linetype = "dashed") +
  ggplot2::labs(y="BWA",
                x="STAR") +
  ggplot2::theme_light() +
  ggplot2::guides(fill = guide_legend(nrow = 2)) + 
  ggplot2::theme(legend.text = ggplot2::element_text(size=6, color='black'),
                 legend.position = "top",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(0.5, "lines"),
                 axis.text.y = ggplot2::element_text(size=8, color='black'),
                 axis.text.x = ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) 

# put them together into one figure
map.combine.plot <- cowplot::plot_grid(mapres.plot,
                                       scrbseq.locations.plot,
                                       rel_widths = c(0.7, 0.3),
                                       nrow = 1,
                                       ncol = 2,
                                       labels = c("A", "B"),
                                       label_size = 10)

#### Detected genes ######

detectgenes.sumdat <- detectgenes.dat %>% 
  dplyr::filter(Annotation %in% c("gencode", "refseq") &
                  Mapper %in% c("star", "bwa", "kallisto") &
                  Protocol %in% c("SCRBseq", "Smartseq2", "Dropseq", "HGMM1kv3", "PBMC1kv3", "CELseq2")) %>% 
  dplyr::filter(Dropout < 0.9) %>% 
  dplyr::group_by(Protocol, Mapper, Annotation) %>% 
  dplyr::tally() %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(! Protocol == "PBMC1kv3") %>% 
  dplyr::mutate(Protocol = ifelse(Protocol == "HGMM1kv3", "10X Genomics", Protocol))

detectgenes.sumdat$Protocol <- factor(detectgenes.sumdat$Protocol, 
                                      levels = c("10X Genomics", "CELseq2", "Dropseq", "SCRBseq", "Smartseq2"
                                      ))

relabels["10X Genomics"] <- "10X Genomics"

genedetect.plot <- ggplot2::ggplot(data = detectgenes.sumdat,
                                   ggplot2::aes(x = Protocol, 
                                                y = n, 
                                                fill = Mapper, 
                                                alpha = 0.5)) + 
  ggplot2::geom_col(position="dodge2", color = "black") + 
  ggplot2::theme_light() +
  ggplot2::scale_fill_manual(values = mapper.cols) +
  ggplot2::scale_x_discrete(labels = relabels) +
  ggplot2::facet_wrap(~Annotation, 
                      scales = "free", labeller = as_labeller(relabels)) +
  ggplot2::labs(x=NULL, y="Detected Genes") +
  ggplot2::guides(alpha=FALSE) + 
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "bottom",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(0.5, "lines"),
                 axis.text=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white"),
                 panel.spacing.y = grid::unit(1.5, "lines")) +
  ggplot2::coord_flip()

##### Parameters ##########

protocols <- c("Smartseq2", "SCRBseq")
annotation <- "gencode"
mapper <- c("bwa", "star", "kallisto")

param_names <- c(`Mean` = paste0("Log[2]", "~~", "Mean"), 
                 `Dispersion` = paste0("Log[2]", "~~", "Dispersion"), 
                 `Dropout` = paste0("Dropout"))
Label_names <- data.frame(param_names) %>% tibble::rownames_to_column(var = 'variable')

margs.dat.red <- margs.dat %>%
  dplyr::filter(Protocol %in% protocols, Annotation == annotation, Mapper %in% mapper) %>% 
  dplyr::mutate(Protocol = recode(Protocol, `Smartseq2` = "Smart-seq2", `SCRBseq` = "SCRB-seq")) %>% 
  dplyr::left_join(Label_names, by = 'variable') %>% 
  dplyr::rename(Parameter = param_names) %>% 
  dplyr::mutate(Parameter = factor(Parameter, 
                                   levels = c("Log[2]~~Mean", "Log[2]~~Dispersion", "Dropout")),
                Protocol = factor(Protocol, 
                                  levels = c("Smart-seq2", "SCRB-seq")))

fits.dat.red <- fits.dat %>%
  dplyr::filter(Protocol %in% protocols, Annotation == annotation, Mapper %in% mapper) %>% 
  dplyr::mutate(Protocol = recode(Protocol, `Smartseq2` = "Smart-seq2", `SCRBseq` = "SCRB-seq")) %>% 
  dplyr::mutate(Protocol = factor(Protocol, 
                                  levels = c("Smart-seq2", "SCRB-seq"))) 

marg.plot <- ggplot2::ggplot(data = margs.dat.red, 
                             ggplot2::aes(x = Mapper, y = value, fill = Mapper)) + 
  ggplot2::geom_violin(position=position_dodge(width=0.9), width=0.8, alpha = 0.5) +
  ggplot2::stat_summary(fun.y = median,
                        fun.ymin = median, 
                        fun.ymax = median,
                        color = "black",
                        position = position_dodge(width=0.9),
                        width = 0.5,
                        geom = "crossbar") +
  ggplot2::theme_light() +
  ggplot2::scale_x_discrete(labels = mapper_names) +
  ggplot2::scale_fill_manual(values = mapper.cols, labels = mapper_names) +
  ggplot2::scale_color_manual(values = mapper.cols, labels = mapper_names) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "top",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(0.5, "lines"),
                 axis.text=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white"),
                 panel.spacing.y = grid::unit(1.5, "lines")) +
  ggplot2::facet_grid(Protocol ~ Parameter, 
                      scales = "free_x", 
                      labeller = label_parsed) +
  ggplot2::labs(x=NULL, y=NULL) +
  ggplot2::coord_flip()

fit.plot <- ggplot2::ggplot(data=fits.dat.red, 
                            ggplot2::aes(x=Mean, 
                                         y=Dispersion, 
                                         fill = Mapper, 
                                         color = Mapper)) +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual(values = mapper.cols, labels = mapper_names) +
  ggplot2::scale_color_manual(values = mapper.cols, labels = mapper_names) +
  ggplot2::geom_ribbon(data=fits.dat.red, 
                       aes(ymin = DispersionLower, ymax = DispersionUpper), 
                       alpha = 0.3) +
  ggplot2::geom_line(data=fits.dat.red,
                     ggplot2::aes(x=Mean, y=Dispersion), 
                     linetype=2, size=1.5) +
  ggplot2::labs(y=expression(bold(paste(Log[2], " Dispersion", sep=""))),
                x=expression(bold(paste(Log[2], " (Mean)")))) +
  ggplot2::theme(legend.position = "top", 
                 legend.title = element_blank(),
                 axis.text=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 axis.line=element_line()) + 
  ggplot2::facet_wrap(~ Protocol, scales = "free", ncol = 1)


###### Performance #######

perf.dat.red <- perf.dat %>% 
  dplyr::filter(Protocol %in% protocols &
                  Mapper %in% c("bwa", "kallisto", "star") &
                  Annotation == annotation &
                  perc.DE == "0.05" &
                  LFC == "Symmetric" &
                  DEtool == "limma-trend" &
                  DEFilter == 'raw' &
                  Normalisation == 'scran' &
                  Preprocessing == 'none' &
                  Spike == "w/o Spike-Ins" &
                  SampleSize == "384 vs 384")

# stats
perf.sts <- perf.dat.red %>%
  dplyr::group_by(Protocol, Mapper) %>%
  dplyr::summarise_at(dplyr::vars(TPR), 
                      dplyr::funs(Min = min(., na.rm=T), 
                                  Max = max(., na.rm=T), 
                                  Mean = mean(., na.rm=T),
                                  Median = median(., na.rm=T), 
                                  N = sum(!is.na(.)),
                                  LowHinge = boxplot.stats(.)[[1]][2],  
                                  UpHinge = boxplot.stats(.)[[1]][4]))
low.whisk <- min(perf.sts$LowHinge, na.rm=T) 
up.whisk <- max(perf.sts$UpHinge, na.rm =T)

# performance plot
perf.dat.red$Mapper <- factor(perf.dat.red$Mapper,
                              levels = c("bwa" , "kallisto" , "star"))
perf.plot <- ggplot2::ggplot(data = perf.dat.red, 
                             aes(x = Mapper, y = TPR)) +
  ggplot2::geom_boxplot(data=perf.dat.red, 
                        aes(x=Mapper, y=TPR, fill = Protocol), 
                        alpha = 0.75, outlier.shape=NA, width = 0.75) + 
  ggplot2:: labs(y="Power (TPR)", x=NULL, fill=NULL) + 
  ggplot2::scale_fill_manual(values = protocols.cols, labels = protocol_names) +
  ggplot2::scale_y_continuous(expand =c(0.01,0.05)) +
  ggplot2::scale_x_discrete(labels = mapper_names) +
  ggplot2::theme_light() +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "top",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(0.5, "lines"),
                 axis.text=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) +
  ggplot2::coord_flip()


#### Combine plots and save ####

# put them together into one figure
fit.perf.combine.plot <- cowplot::plot_grid(fit.plot  + theme(legend.position = "none"),
                                            perf.plot,
                                            nrow = 1,
                                            ncol = 2,
                                            rel_widths = c(1, 0.6),
                                            labels = c("E", "F"),
                                            label_size = 12)

combine.plot <- cowplot::plot_grid(map.combine.plot,
                                   genedetect.plot,
                                   marg.plot + theme(legend.position = "none"),
                                   fit.perf.combine.plot,
                                   ncol = 1,
                                   rel_heights = c(0.8, 0.75, 0.9, 1.1),
                                   labels = c("", "C",  "D", ""),
                                   label_size = 12)

ggplot2::ggsave(filename=paste0("Quant_Results_", annotation, ".pdf"), 
                plot = combine.plot,
                width=210,
                height=297,
                units="mm")

ggplot2::ggsave(filename=paste0("Quant_Results_", annotation, ".png"), 
                plot = combine.plot,
                width=210,
                height=297,
                units="mm")


# SUPPLEMENTARY INPUT -----------------------------------------------------

mapres.si.dat <- readRDS(file = "MapStats_SI_2019-05-27.rds") %>% 
  dplyr::filter( Feature %in% c("Aligned", "Aligned+Assigned", "UMI") &
                   Protocol %in% c("CELseq2", "Dropseq", "SCRBseq", "Smartseq2", "PBMC1kv3", "HGMM1kv3") &
                   Mapper %in% c("star", 'bwa', 'kallisto'))

mapres.si.dat$Protocol <- factor(mapres.si.dat$Protocol , 
                                 levels = c("PBMC1kv3", "HGMM1kv3", "CELseq2", "Dropseq", "SCRBseq", "Smartseq2"))
mapres.si.dat$Mapper <- factor(mapres.si.dat$Mapper, levels = c("star", 'bwa', 'kallisto'))
mapres.si.dat$Annotation <- factor(mapres.si.dat$Annotation, levels = c("gencode", 'vega', 'refseq'))

detectgenes.si.dat <- readRDS(file = "MargsDat.rds") %>%
  dplyr::select(-Mean, -Dispersion) %>% 
  dplyr::filter(Protocol %in% c("SCRBseq", "Smartseq2", "Dropseq", "HGMM1kv3", "PBMC1kv3", "CELseq2") & Mapper %in% c("star", 'bwa', 'kallisto')) %>% 
  dplyr::filter(Dropout < 0.9)
detectgenes.si.dat$Protocol <- factor(detectgenes.si.dat$Protocol , 
                                      levels = c("PBMC1kv3", "HGMM1kv3", "CELseq2", "Dropseq", "SCRBseq", "Smartseq2"))
detectgenes.si.dat$Mapper <- factor(detectgenes.si.dat$Mapper, levels = c("star", 'bwa', 'kallisto'))
detectgenes.si.dat$Annotation <- factor(detectgenes.si.dat$Annotation, levels = c("gencode", 'vega', 'refseq'))

margs.si.dat <- readRDS(file = "MargsDat.rds")
margs.si.dat <-   data.table::melt(setDT(margs.si.dat), 
                                   id.vars = c("Name", "Protocol", "Mapper", "Annotation", "GeneID"),
                                   measure.vars = c("Mean", "Dispersion", "Dropout"),
                                   variable.factor = FALSE) %>% 
  data.frame()

fits.si.dat <- readRDS(file = "FitsDat.rds")

perf.si.dat <- readRDS(file = "AlignerPerformance.rds")

# SUPPLEMENT FIGURES ------------------------------------------------------

# (1) Mapping stats for all mappers, protocols and annotations
mapres.si.plot <- ggplot(data = mapres.si.dat, 
                         aes(x = Protocol, 
                             y = Proportion, fill = Protocol, group = Protocol, alpha = Feature)) +
  geom_bar(stat="identity", width=.8, position = "dodge", color="black") + 
  scale_y_continuous(labels = scales::percent, breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1)) +
  scale_x_discrete(labels = protocol_names) +
  labs(y="Percentage of Total Reads", x=NULL, fill=NULL) + 
  facet_grid(Annotation ~ Mapper, scales="fixed", labeller = as_labeller(relabels)) + 
  scale_fill_manual(values = protocols.cols) +
  scale_alpha_manual(values=c(0.1,0.6, 1)) +
  theme_light() +
  theme(legend.text = ggplot2::element_text(size=8, color='black'),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        legend.key.size = grid::unit(1, "lines"),
        axis.text.y = ggplot2::element_text(size=8,color='black'),
        axis.text.x = ggplot2::element_text(size=8, color='black'),
        axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
        strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
        strip.background = ggplot2::element_rect(fill="white")) +
  coord_flip()

ggplot2::ggsave(filename=paste0("Quant_Results_Mapping_SI", ".pdf"), 
                plot = mapres.si.plot,
                width=210,
                height=148,
                units="mm")

ggplot2::ggsave(filename=paste0("Quant_Results_Mapping_SI", ".png"), 
                plot = mapres.si.plot,
                width=210,
                height=148,
                units="mm")

# (2) Number of detected genes per annotation ; stratify over protocol and aligner
detectgenes.sum.dat <- detectgenes.si.dat %>% 
  dplyr::group_by(Protocol, Mapper, Annotation) %>% 
  dplyr::tally()
maxgenes <- max(detectgenes.sum.dat$n)

protocol_names <- c(`CELseq2` = "CEL-seq2", 
                    `SCRBseq` = "SCRB-seq", 
                    `Smartseq2` = "Smart-seq2", 
                    `HGMM1kv3` = "TenX-HGMM", 
                    `Dropseq` = "Drop-seq")

genedetect.si.plot <- ggplot2::ggplot(data = detectgenes.sum.dat,
                                      ggplot2::aes(x = Protocol, y = n, 
                                                   fill = Annotation, alpha = 0.75)) + 
  ggplot2::geom_col(position="dodge2", color = "black") + 
  ggplot2::theme_light() +
  ggplot2::scale_x_discrete(labels = protocol_names) +
  ggplot2::scale_y_continuous(breaks=c(0,5000, 10000, 15000, 20000), limits=c(0,maxgenes)) +
  ggplot2::scale_fill_manual(values = annotation.cols, labels = annotation_names) +
  ggplot2::guides(alpha = FALSE) + 
  ggplot2::labs(x=NULL, y="Detected Genes") +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "top",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(0.5, "lines"),
                 axis.text=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white"),
                 panel.spacing.y = grid::unit(1.5, "lines")) +
  ggplot2::facet_wrap(~Mapper, 
                      ncol = 1,
                      scales = "free", labeller = as_labeller(relabels)) +
  ggplot2::coord_flip()

ggplot2::ggsave(filename=paste0("Quant_Results_DetectGenes_SI", ".pdf"), 
                plot = genedetect.si.plot,
                width=210,
                height=210,
                units="mm")

ggplot2::ggsave(filename=paste0("Quant_Results_DetectGenes_SI", ".png"), 
                plot = genedetect.si.plot,
                width=210,
                height=210,
                units="mm")

# (3) Marginal distribution of parameters per protocol colored by annotation grouped by mapper
protocols <- c("Smartseq2", "SCRBseq", "CELseq2", "Dropseq", "HGMM1kv3")
annotations <- c("gencode", "vega", "refseq")
mappers <- c("star", "bwa", "kallisto")

param_names <- c(`Mean` = paste0("Log[2]", "~~", "Mean"), 
                 `Dispersion` = paste0("Log[2]", "~~", "Dispersion"), 
                 `Dropout` = paste0("Dropout"))
paramlabel_names <- data.frame(param_names) %>% 
  tibble::rownames_to_column(var = 'variable')
protocols_names <- c(`CELseq2` = "CEL-seq2", 
                     `SCRBseq` = "SCRB-seq", 
                     `Smartseq2` = "Smart-seq2", 
                     `HGMM1kv3` = "TenX-HGMM", 
                     `Dropseq` = "Drop-seq")
mapper_names <- c(`star` = "STAR", `bwa` = "BWA", `kallisto` = "kallisto")

protocollabel_names <- data.frame(protocols_names) %>% 
  tibble::rownames_to_column(var = 'Protocol')

margs.dat.red <- margs.si.dat %>%
  dplyr::filter(Mapper %in% mappers &
                  Annotation %in% annotations &
                  Protocol %in% protocols) %>% 
  dplyr::left_join(paramlabel_names, by = 'variable') %>% 
  dplyr::rename(Parameter = param_names) %>% 
  dplyr::left_join(protocollabel_names, by = 'Protocol') %>% 
  dplyr::rename(Protocols = protocols_names) %>%
  dplyr::mutate(Parameter = factor(Parameter, 
                                   levels = c("Log[2]~~Mean", 
                                              "Log[2]~~Dispersion", 
                                              "Dropout")),
                Protocols = factor(Protocols, 
                                   levels = c("Smart-seq2", 
                                              "SCRB-seq", 
                                              "Drop-seq", 
                                              "CEL-seq2", 
                                              "TenX-HGMM")),
                Mapper = factor(Mapper,
                                levels = c("kallisto", "bwa", "star")),
                Annotation = factor(Annotation,
                                    levels = c("gencode", 'vega', 'refseq')))

marg.si.plot <- ggplot2::ggplot(data = margs.dat.red, 
                                ggplot2::aes(x = Mapper, y = value, fill = Annotation)) + 
  ggplot2::geom_violin(position=position_dodge(width=0.9),
                       width=0.8, 
                       alpha = 0.75) +
  ggplot2::stat_summary(fun.y = median,
                        fun.ymin = median, 
                        fun.ymax = median,
                        color = "black",
                        position = position_dodge(width=0.9),
                        width = 0.5,
                        geom = "crossbar") +
  ggplot2::theme_light() +
  ggplot2::scale_x_discrete(labels = mapper_names) +
  ggplot2::scale_fill_manual(values = annotation.cols, labels = annotation_names) +
  ggplot2::scale_color_manual(values = annotation.cols, labels = annotation_names) +
  ggplot2::guides(color = FALSE) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "top",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white"),
                 panel.spacing.y = grid::unit(1.5, "lines")) +
  ggplot2::facet_wrap(Protocols ~ Parameter, 
                      scales = "free", 
                      ncol = 3,
                      labeller = label_parsed) +
  ggplot2::labs(x=NULL, y=NULL) +
  ggplot2::coord_flip()

ggplot2::ggsave(filename=paste0("Quant_Results_MargParams_SI", ".pdf"), 
                plot = marg.si.plot,
                width=210,
                height=297,
                units="mm")

ggplot2::ggsave(filename=paste0("Quant_Results_MargParams_SI", ".png"), 
                plot = marg.si.plot,
                width=210,
                height=297,
                units="mm")

# (4) mean-variance relationship Annotation ~ Protocol colored by mapper
protocols <- c("Smartseq2", "SCRBseq", "CELseq2", "Dropseq", "HGMM1kv3")
annotations <- c("gencode", "vega", "refseq")
mappers <- c("star", "bwa", "kallisto")

fits.dat.red <- fits.si.dat %>%
  dplyr::filter(Mapper %in% mappers &
                  Annotation %in% annotations &
                  Protocol %in% protocols) %>% 
  dplyr::mutate(Protocol = recode(Protocol, 
                                  `Smartseq2` = "Smart-seq2", 
                                  `SCRBseq` = "SCRB-seq",
                                  `Dropseq` = "Drop-seq",
                                  `CELseq2` = "CEL-seq2",
                                  `HGMM1kv3` = "10X-HGMM")) %>% 
  dplyr::mutate(Protocol = factor(Protocol, 
                                  levels = c("Smart-seq2", 
                                             "SCRB-seq",
                                             "Drop-seq",
                                             "CEL-seq2",
                                             "10X-HGMM"))) %>% 
  dplyr::mutate(Annotation = recode(Annotation, 
                                    `gencode` = "GENCODE",
                                    `vega` = "Vega",
                                    `refseq` = "RefSeq")) %>% 
  dplyr::mutate(Annotation = factor(Annotation, 
                                    levels = c("GENCODE",
                                               "Vega",
                                               "RefSeq")))

fit.si.plot <- ggplot2::ggplot(data=fits.dat.red, 
                               ggplot2::aes(x=Mean, 
                                            y=Dispersion, 
                                            fill = Mapper, 
                                            color = Mapper)) +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual(values = mapper.cols, labels = mapper_names) +
  ggplot2::scale_color_manual(values = mapper.cols, labels = mapper_names) +
  ggplot2::geom_ribbon(data=fits.dat.red, 
                       aes(ymin = DispersionLower, 
                           ymax = DispersionUpper), 
                       alpha = 0.3) +
  ggplot2::geom_line(data=fits.dat.red,
                     ggplot2::aes(x=Mean, 
                                  y=Dispersion), 
                     linetype=2, size=1.5) +
  ggplot2::labs(y=expression(bold(paste(Log[2], " Dispersion", sep=""))),
                x=expression(bold(paste(Log[2], " (Mean)")))) +
  ggplot2::theme(legend.position = "top", 
                 legend.title = element_blank(),
                 axis.text=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 axis.line=element_line()) + 
  ggplot2::facet_wrap(Protocol ~ Annotation, scales = "free", ncol = 3)


ggplot2::ggsave(filename=paste0("Quant_Results_Fit_SI", ".pdf"), 
                plot = fit.si.plot,
                width=210,
                height=297,
                units="mm")

ggplot2::ggsave(filename=paste0("Quant_Results_Fit_SI", ".png"), 
                plot = fit.si.plot,
                width=210,
                height=297,
                units="mm")

# (5) TPR per protocol stratified over aligner and annotation

perf.si.sts <- perf.si.dat %>%
  dplyr::filter(Protocol %in% c("10XGenomics", "CELseq2", "Dropseq", "SCRBseq", "Smartseq2") &
                  Mapper %in% c("bwa", "kallisto", "star") &
                  Annotation == c("gencode", "vega", "refseq") &
                  perc.DE == "0.05" &
                  LFC == "Symmetric" &
                  DEtool == "limma-trend" &
                  DEFilter == 'raw' &
                  Normalisation == 'scran' &
                  Preprocessing == 'none' &
                  Spike == "w/o Spike-Ins" &
                  SampleSize == "384 vs 384")

perf.si.plot <- ggplot2::ggplot(data = perf.si.sts,
                                aes(x = Protocol, 
                                    y = TPR, 
                                    fill = Annotation,
                                    alpha = 0.75)) +
  ggplot2::geom_boxplot(outlier.shape=NA, width = 1,
                        position=ggplot2::position_dodge2(1.25, preserve = "single", 
                                                          padding = 0.25)) +
  ggplot2::theme_light() +
  ggplot2::guides(alpha = FALSE) +
  ggplot2::scale_fill_manual(values = annotation.cols, labels = annotation_names) +
  scale_y_continuous(labels = scales::percent, 
                     breaks=c(0.5, 0.6, 0.7, 0.8, 0.9,1), 
                     limits=c(0.5,1)) +
  scale_x_discrete(labels = protocol_names) +
  ggplot2::labs(x = NULL) + 
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "top", 
                 legend.title = element_blank(),
                 axis.text=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title = ggplot2::element_text(size=8, face="bold", color='black', hjust = 0.5),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) +
  ggplot2::facet_wrap(~ Mapper, ncol = 1, scales = "free", 
                      labeller = as_labeller(relabels)) + 
  ggplot2::coord_flip()

ggplot2::ggsave(filename=paste0("Quant_Results_TPR_SI", ".pdf"), 
                plot = perf.si.plot,
                width=140,
                height=160,
                units="mm")

ggplot2::ggsave(filename=paste0("Quant_Results_TPR_SI", ".png"), 
                plot = perf.si.plot,
                width=210,
                height=297,
                units="mm")
