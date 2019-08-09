# SETTINGS ----------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(xtable)
# library(powsimRDev)

# relabelling annotations for plotting
source(file  = "renaming_paper.R")

# plotting colors
source(file  = "colour_schemes_paper.R")

# INPUT -------------------------------------------------------------------

evaljobs.dat <- readRDS(file = "evalfulldat.rds")

umi.info <- data.frame(Protocol = c("10XGenomics", "CELseq2", "Dropseq", "SCRBseq", "Smartseq2"),
                       Type = c("UMI", "UMI", "UMI", "UMI", "Full-Length"), 
                       stringsAsFactors = F)

# add UMI vs nonUMI method column
evaljobs.dat <- evaljobs.dat %>%
  dplyr::left_join(umi.info, by="Protocol") 

tpr <- readRDS(file = "TPR.rds")
fdr <- readRDS(file = "FDR.rds")
tf.pauc <- readRDS(file = "TPRvsFDR_pAUC_lib.rds")


# MAIN FIGURE -------------------------------------------------------------

# main figure with pAUC for x axis with DE-tools, 
# for 20 % asymmetric
# stratified over protocols ( Smartseq2, 10X, SCRB) and normalisation methods, 
# with alpha = spike-ins

main.plot.detool <- function(mapper, 
                             annotation,
                             normalisations,
                             detools,
                             protocols,
                             samplesize,
                             percde,
                             lfc) {
  evaljobs.dat2 <- evaljobs.dat
  evaljobs.dat2$Normalisation <- forcats::fct_rev(evaljobs.dat2$Normalisation)
  jobs2eval <- evaljobs.dat2 %>% dplyr::filter(Protocol %in% protocols &
                                                 Mapper == mapper &
                                                 Annotation ==  annotation &
                                                 Normalisation %in% normalisations &
                                                 DEtool %in% detools &
                                                 DEFilter == "raw" &
                                                 Preprocessing == "none" &
                                                 LFC == lfc &
                                                 perc.DE == percde) %>% 
    droplevels()
  
  # AUC
  auc.dat <- jobs2eval %>%
    dplyr::select(Name, Normalisation, DEtool, Spike, Protocol, LFC, perc.DE) %>%
    dplyr::left_join( tf.pauc , by = "Name") %>% 
    dplyr::filter(SampleSize == samplesize) %>%
    dplyr::rename(pAUC=`TPRvsFDR_AUC`)
  
  sts.auc <- auc.dat %>%
    dplyr::group_by(Normalisation, DEtool, Spike, Protocol, LFC, perc.DE) %>%
    dplyr::summarise_at(dplyr::vars(pAUC), 
                        dplyr::funs(Min = min(., na.rm=T), 
                                    Max = max(., na.rm=T), 
                                    Mean = mean(., na.rm=T), 
                                    MeanBoot = Hmisc::smean.cl.boot(.)[1],
                                    SE = sd(., na.rm=T)/n(),
                                    SD = sd(., na.rm=T),
                                    Median = median(., na.rm=T), 
                                    N = sum(!is.na(.)),
                                    LowHinge = boxplot.stats(.)[[1]][2],  
                                    UpHinge = boxplot.stats(.)[[1]][4],
                                    LowerBoot = Hmisc::smean.cl.boot(.)[2],
                                    UpperBoot = Hmisc::smean.cl.boot(.)[3]))
  
  # plot annotations
  percwritten <- c(`0.05` = "5%", `0.2` = "20%", `0.6` = "60 %")
  titleCap <- paste0("Performance of DE methods to detect ",  
                     percwritten[percde], " ", lfc)
  relabels["scran w/ cluster"] <- "scran \n with cluster"
  relabels["scran w/ groups"] <- "scran \n with groups"
  relabels["SCnorm"] <- "SCnorm \n with groups"
  relabels["SCnorm w/ cluster"] <- "SCnorm \n with cluster"
  relabels["SF"] <- "Simulated \nSize Factors"
  
  
  pdat <- sts.auc %>%  dplyr::filter(LFC == lfc & perc.DE == percde)
  
  # AUC plot
  auc.plot <- ggplot2::ggplot(data = pdat, 
                              aes(x = DEtool, 
                                  y = MeanBoot, 
                                  color = DEtool, 
                                  alpha = Spike)) +
    ggplot2::geom_pointrange(aes(x = DEtool, 
                                 ymin = LowerBoot, 
                                 ymax = UpperBoot, 
                                 color = DEtool, alpha = Spike),
                             size = 0.5, position = position_dodge(width = 1)) +
    ggplot2::geom_rect(data = pdat %>%  dplyr::filter(Normalisation == "SF"),
                       xmin = -Inf, xmax = Inf,
                       ymin = -Inf,ymax = Inf,
                       alpha = 0.05, fill ="grey") +
    ggplot2::scale_color_manual(values = detool.cols) +
    scale_alpha_manual(values=c(1,0.4)) +
    ggplot2::theme_light() + 
    ggplot2::facet_grid(Normalisation ~ Protocol, scales = "fixed",
                        labeller = as_labeller(relabels)) +
    ggplot2::labs(x = NULL,
                  y = "pAUC",
                  title = titleCap) +
    ggplot2::scale_y_continuous(limits=c(0, 1), breaks = c(0,0.25,0.5,0.75,1)) +
    ggplot2::scale_x_discrete(labels = detool_names) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                   legend.position = "none",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(1.5, "lines"),
                   axis.text.x=ggplot2::element_text(size=8, color='black'),
                   axis.text.y=ggplot2::element_text(size=8, color='black'),
                   axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white")) + 
    ggplot2::coord_flip()
  
  ggplot2::ggsave(filename=paste0("Performance_DE-tool_Main_", 
                                  mapper, "_", 
                                  annotation, "_", 
                                  gsub(" ", "", samplesize, fixed = TRUE), "_",
                                  percde, "_",
                                  lfc, 
                                  ".pdf"), 
                  plot = auc.plot,
                  width=180,
                  height=160,
                  units="mm")
  ggplot2::ggsave(filename=paste0("Performance_DE-tool_Main_", 
                                  mapper, "_", 
                                  annotation, "_", 
                                  gsub(" ", "", samplesize, fixed = TRUE), "_",
                                  percde, "_",
                                  lfc,
                                  ".png"), 
                  plot = auc.plot,
                  width=180,
                  height=160,
                  units="mm")
}

# grid of parameter combinations
mapper <- c('bwa', 'star', 'kallisto')
annotation <- c('gencode', 'refseq')
normalisations <- c("SF", "MR", "SCnorm",
                    "scran", "scran w/ cluster")
detools <- c('limma-trend', 
             "MAST", 
             'T-Test', 
             "edgeR-zingeR")
protocols <- c("10XGenomics", 
               "SCRBseq",
               "Smartseq2")
lfc <- c("Symmetric", "Asymmetric", "Completely Asymmetric")
percde <- c("0.05", "0.2", "0.6")
samplesize <- c("384 vs 384", "96 vs 96", "50 vs 200")
par.sim <- expand.grid(mapper=mapper, 
                       annotation=annotation,
                       percde = percde,
                       lfc=lfc,
                       samplesize = samplesize,
                       stringsAsFactors = F)

par.sim <- par.sim %>% 
  dplyr::filter(mapper == "star" &
                  annotation == "gencode" &
                  lfc == "Asymmetric" &
                  percde == "0.2" &
                  samplesize == "384 vs 384")

for(i in 1:nrow(par.sim)) {
  print(paste0(i , " out of ", nrow(par.sim)))
  tryCatch(main.plot.detool(mapper = par.sim[i, "mapper"], 
                            annotation = par.sim[i, "annotation"],
                            normalisations = normalisations,
                            detools = detools,
                            protocols = protocols,
                            lfc = par.sim[i,"lfc"],
                            percde = par.sim[i,"percde"],
                            samplesize = par.sim[i, "samplesize"]), 
           error = function() next)
}

# SUPPLEMENTS -------------------------------------------------------------

si.plot.detool <- function(mapper, 
                           annotation,
                           normalisations,
                           detools,
                           protocols,
                           samplesize,
                           percde,
                           lfc) {
  evaljobs.dat2 <- evaljobs.dat
  evaljobs.dat2$Normalisation <- forcats::fct_rev(evaljobs.dat2$Normalisation)
  jobs2eval <- evaljobs.dat2 %>% dplyr::filter(Protocol %in% protocols &
                                                 Mapper == mapper &
                                                 Annotation ==  annotation &
                                                 Normalisation %in% normalisations &
                                                 DEtool %in% detools &
                                                 DEFilter == "raw" &
                                                 Preprocessing == "none" &
                                                 LFC == lfc &
                                                 perc.DE == percde) %>% 
    droplevels()
  
  # shorten facet_wrap strip labels
  relabels["scran w/ cluster"] <- "scran \n with cluster"
  relabels["scran w/ groups"] <- "scran \n with groups"
  relabels["SCnorm"] <- "SCnorm \n with groups"
  relabels["SCnorm w/ cluster"] <- "SCnorm \n with cluster"
  relabels["SF"] <- "Simulated \nSize Factors"
  # title annotation
  percwritten <- c(`0.05` = "5%", `0.2` = "20%", `0.6` = "60 %")
  
  # AUC
  auc.dat <- jobs2eval %>%
    dplyr::select(Name, Normalisation, DEtool, Spike, Protocol, LFC, perc.DE) %>%
    dplyr::left_join( tf.pauc , by = "Name") %>% 
    dplyr::filter(SampleSize == samplesize) %>%
    dplyr::rename(pAUC=`TPRvsFDR_AUC`)
  
  sts.auc <- auc.dat %>%
    dplyr::group_by(Normalisation, DEtool, Spike, Protocol, LFC, perc.DE) %>%
    dplyr::summarise_at(dplyr::vars(pAUC), 
                        dplyr::funs(Min = min(., na.rm=T), 
                                    Max = max(., na.rm=T), 
                                    Mean = mean(., na.rm=T), 
                                    MeanBoot = Hmisc::smean.cl.boot(.)[1],
                                    SE = sd(., na.rm=T)/n(),
                                    SD = sd(., na.rm=T),
                                    Median = median(., na.rm=T), 
                                    N = sum(!is.na(.)),
                                    LowHinge = boxplot.stats(.)[[1]][2],  
                                    UpHinge = boxplot.stats(.)[[1]][4],
                                    LowerBoot = Hmisc::smean.cl.boot(.)[2],
                                    UpperBoot = Hmisc::smean.cl.boot(.)[3]))
  
  # plot annotation for AUC
  titleCap <- paste0("Performance of DE methods to detect ",  
                     percwritten[percde], " ", lfc)
  
  
  auc.pdat <- sts.auc %>%  dplyr::filter(LFC == lfc & perc.DE == percde)
  
  # AUC plot
  auc.plot <- ggplot2::ggplot(data = auc.pdat, 
                              aes(x = DEtool, 
                                  y = MeanBoot, 
                                  color = DEtool, 
                                  alpha = Spike)) +
    ggplot2::geom_pointrange(aes(x = DEtool, 
                                 ymin = LowerBoot, 
                                 ymax = UpperBoot, 
                                 color = DEtool, alpha = Spike),
                             size = 0.5, position = position_dodge(width = 1)) +
    ggplot2::geom_rect(data = auc.pdat %>%  dplyr::filter(Normalisation == "SF"),
                       xmin = -Inf, xmax = Inf,
                       ymin = -Inf,ymax = Inf,
                       alpha = 0.05, fill ="grey") +
    ggplot2::scale_color_manual(values = detool.cols) +
    scale_alpha_manual(values=c(1,0.4)) +
    ggplot2::theme_light() + 
    ggplot2::facet_grid(Normalisation ~ Protocol, scales = "fixed",
                        labeller = as_labeller(relabels)) +
    ggplot2::labs(x = NULL,
                  y = "pAUC",
                  title = titleCap) +
    ggplot2::scale_y_continuous(limits=c(0, 1), breaks = c(0,0.25,0.5,0.75,1)) +
    ggplot2::scale_x_discrete(labels = detool_names) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                   legend.position = "none",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(1.5, "lines"),
                   axis.text.x=ggplot2::element_text(size=8, color='black'),
                   axis.text.y=ggplot2::element_text(size=8, color='black'),
                   axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white")) + 
    ggplot2::coord_flip()
  
  ggplot2::ggsave(filename=paste0("Performance_DE-tool_SI_pAUC_", 
                                  mapper, "_", 
                                  annotation, "_", 
                                  gsub(" ", "", samplesize, fixed = TRUE), "_",
                                  percde, "_",
                                  lfc, 
                                  ".pdf"), 
                  plot = auc.plot,
                  width=210,
                  height=297,
                  units="mm")
  ggplot2::ggsave(filename=paste0("Performance_DE-tool_SI_pAUC_", 
                                  mapper, "_", 
                                  annotation, "_", 
                                  gsub(" ", "", samplesize, fixed = TRUE), "_",
                                  percde, "_",
                                  lfc,
                                  ".png"), 
                  plot = auc.plot,
                  width=210,
                  height=297,
                  units="mm")
  
  
  # FDR
  fdr.l <-  fdr %>% 
    tidyr::gather("SampleSize", "FDR", -Name) %>% 
    dplyr::filter(SampleSize == samplesize)
  
  fdr.dat <- jobs2eval %>%
    dplyr::select(Name, Normalisation, DEtool, Spike, Protocol, LFC, perc.DE) %>%
    dplyr::left_join( fdr.l , by = "Name") 
  
  sts.fdr <- fdr.dat %>%
    dplyr::group_by(Normalisation, DEtool, Spike, Protocol, LFC, perc.DE) %>%
    dplyr::summarise_at(dplyr::vars(FDR), 
                        dplyr::funs(Min = min(., na.rm=T), 
                                    Max = max(., na.rm=T), 
                                    Mean = mean(., na.rm=T), 
                                    MeanBoot = Hmisc::smean.cl.boot(.)[1],
                                    SE = sd(., na.rm=T)/n(),
                                    SD = sd(., na.rm=T),
                                    Median = median(., na.rm=T), 
                                    N = sum(!is.na(.)),
                                    LowHinge = boxplot.stats(.)[[1]][2],  
                                    UpHinge = boxplot.stats(.)[[1]][4],
                                    LowerBoot = Hmisc::smean.cl.boot(.)[2],
                                    UpperBoot = Hmisc::smean.cl.boot(.)[3]))
  
  # plot annotation for AUC
  titleCap <- paste0("FDR control of DE methods in the presence of ",  
                     percwritten[percde], " ", lfc)
  
  fdr.pdat <- sts.fdr %>%  dplyr::filter(LFC == lfc & perc.DE == percde)
  
  # FDR plot
  fdr.plot <- ggplot2::ggplot(data = fdr.pdat, 
                              aes(x = DEtool, 
                                  y = MeanBoot, 
                                  color = DEtool, 
                                  alpha = Spike)) +
    ggplot2::geom_pointrange(aes(x = DEtool, 
                                 ymin = LowerBoot, 
                                 ymax = UpperBoot, 
                                 color = DEtool, alpha = Spike),
                             size = 0.5, position = position_dodge(width = 1)) +
    ggplot2::geom_rect(data = fdr.pdat %>%  dplyr::filter(Normalisation == "SF"),
                       xmin = -Inf, xmax = Inf,
                       ymin = -Inf,ymax = Inf,
                       alpha = 0.05, fill ="grey") +
    ggplot2::scale_color_manual(values = detool.cols) +
    scale_alpha_manual(values=c(1,0.4)) +
    ggplot2::theme_light() + 
    ggplot2::facet_grid(Normalisation ~ Protocol, scales = "fixed",
                        labeller = as_labeller(relabels)) +
    ggplot2::labs(x = NULL,
                  y = "FDR",
                  title = titleCap) +
    ggplot2::geom_hline(yintercept = c(0.1), colour="darkgrey", linetype="dashed") + 
    ggplot2::scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4), 
                                labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_x_discrete(labels = detool_names) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                   legend.position = "none",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(1.5, "lines"),
                   axis.text.x=ggplot2::element_text(size=8, color='black'),
                   axis.text.y=ggplot2::element_text(size=8, color='black'),
                   axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white")) + 
    ggplot2::coord_flip()
  
  ggplot2::ggsave(filename=paste0("Performance_DE-tool_SI_FDR_", 
                                  mapper, "_", 
                                  annotation, "_", 
                                  gsub(" ", "", samplesize, fixed = TRUE), "_",
                                  percde, "_",
                                  lfc, 
                                  ".pdf"), 
                  plot = fdr.plot,
                  width=210,
                  height=297,
                  units="mm")
  ggplot2::ggsave(filename=paste0("Performance_DE-tool_SI_FDR_", 
                                  mapper, "_", 
                                  annotation, "_", 
                                  gsub(" ", "", samplesize, fixed = TRUE), "_",
                                  percde, "_",
                                  lfc,
                                  ".png"), 
                  plot = fdr.plot,
                  width=210,
                  height=297,
                  units="mm")
  
  
  # TPR
  tpr.l <-  tpr %>% 
    tidyr::gather("SampleSize", "TPR", -Name) %>% 
    dplyr::filter(SampleSize == samplesize)
  
  tpr.dat <- jobs2eval %>%
    dplyr::select(Name, Normalisation, DEtool, Spike, Protocol, LFC, perc.DE) %>%
    dplyr::left_join( tpr.l , by = "Name") 
  
  sts.tpr <- tpr.dat %>%
    dplyr::group_by(Normalisation, DEtool, Spike, Protocol, LFC, perc.DE) %>%
    dplyr::summarise_at(dplyr::vars(TPR), 
                        dplyr::funs(Min = min(., na.rm=T), 
                                    Max = max(., na.rm=T), 
                                    Mean = mean(., na.rm=T), 
                                    MeanBoot = Hmisc::smean.cl.boot(.)[1],
                                    SE = sd(., na.rm=T)/n(),
                                    SD = sd(., na.rm=T),
                                    Median = median(., na.rm=T), 
                                    N = sum(!is.na(.)),
                                    LowHinge = boxplot.stats(.)[[1]][2],  
                                    UpHinge = boxplot.stats(.)[[1]][4],
                                    LowerBoot = Hmisc::smean.cl.boot(.)[2],
                                    UpperBoot = Hmisc::smean.cl.boot(.)[3]))
  
  # plot annotation for AUC
  titleCap <- paste0("Power of DE methods to detect ",  
                     percwritten[percde], " ", lfc)
  
  tpr.pdat <- sts.tpr %>%  dplyr::filter(LFC == lfc & perc.DE == percde)
  
  # TPR plot
  tpr.plot <- ggplot2::ggplot(data = tpr.pdat, 
                              aes(x = DEtool, 
                                  y = MeanBoot, 
                                  color = DEtool, 
                                  alpha = Spike)) +
    ggplot2::geom_pointrange(aes(x = DEtool, 
                                 ymin = LowerBoot, 
                                 ymax = UpperBoot, 
                                 color = DEtool, alpha = Spike),
                             size = 0.5, position = position_dodge(width = 1)) +
    ggplot2::geom_rect(data = tpr.pdat %>%  dplyr::filter(Normalisation == "SF"),
                       xmin = -Inf, xmax = Inf,
                       ymin = -Inf,ymax = Inf,
                       alpha = 0.05, fill ="grey") +
    ggplot2::scale_color_manual(values = detool.cols) +
    scale_alpha_manual(values=c(1,0.4)) +
    ggplot2::theme_light() + 
    ggplot2::facet_grid(Normalisation ~ Protocol, scales = "fixed",
                        labeller = as_labeller(relabels)) +
    ggplot2::labs(x = NULL,
                  y = "TPR",
                  title = titleCap) +
    ggplot2::scale_y_continuous(limits = c(0.5,1),
                                breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1), 
                                labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_x_discrete(labels = detool_names) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                   legend.position = "none",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(1.5, "lines"),
                   axis.text.x=ggplot2::element_text(size=8, color='black'),
                   axis.text.y=ggplot2::element_text(size=8, color='black'),
                   axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white")) + 
    ggplot2::coord_flip()
  
  ggplot2::ggsave(filename=paste0("Performance_DE-tool_SI_TPR_", 
                                  mapper, "_", 
                                  annotation, "_", 
                                  gsub(" ", "", samplesize, fixed = TRUE), "_",
                                  percde, "_",
                                  lfc, 
                                  ".pdf"), 
                  plot = tpr.plot,
                  width=210,
                  height=297,
                  units="mm")
  ggplot2::ggsave(filename=paste0("Performance_DE-tool_SI_TPR_", 
                                  mapper, "_", 
                                  annotation, "_", 
                                  gsub(" ", "", samplesize, fixed = TRUE), "_",
                                  percde, "_",
                                  lfc,
                                  ".png"), 
                  plot = tpr.plot,
                  width=210,
                  height=297,
                  units="mm")
  
}

# grid of parameter combinations
mapper <- c('bwa', 'star', 'kallisto')
annotation <- c('gencode', 'refseq')
normalisations <- c("SF", "MR", "PosCounts", "TMM", "Linnorm", 
                    "SCnorm", "SCnorm w/ cluster", 
                    "scran", "scran w/ groups", "scran w/ cluster")
detools <- c('limma-trend', 
             "MAST", 
             'T-Test', 
             "edgeR-zingeR")
protocols <- c("SCRBseq", "Smartseq2", "10XGenomics")
lfc <- c("Symmetric", "Asymmetric", "Completely Asymmetric")
percde <- c("0.05", "0.2", "0.6")
samplesize <- c("384 vs 384", "96 vs 96", "50 vs 200")
par.sim <- expand.grid(mapper=mapper, 
                       annotation=annotation,
                       percde = percde,
                       lfc=lfc,
                       samplesize = samplesize,
                       stringsAsFactors = F)

par.sim <- par.sim %>% 
  dplyr::filter(mapper == "star" &
                  annotation == "gencode" &
                  lfc %in% c("Asymmetric", "Symmetric") &
                  percde %in% c("0.05", "0.2") &
                  samplesize == "384 vs 384")

for(i in 1:nrow(par.sim)) {
  print(paste0(i , " out of ", nrow(par.sim)))
  tryCatch(si.plot.detool(mapper = par.sim[i, "mapper"], 
                          annotation = par.sim[i, "annotation"],
                          normalisations = normalisations,
                          detools = detools,
                          protocols = protocols,
                          lfc = par.sim[i,"lfc"],
                          percde = par.sim[i,"percde"],
                          samplesize = par.sim[i, "samplesize"]), 
           error = function() next)
}
