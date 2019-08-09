# SETTINGS ----------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(xtable)
# library(powsimRDev)

# plotting colors
source(file  = "colour_schemes_paper.R")
# renaaming variables
source(file  = "renaming_paper.R")

# INPUT FOR MAIN FIGURE ---------------------------------------------------

# (A) RMSE for norm methods for UMI data for 20 % asymmetric, 384 vs 384, w/o spike-ins
rmse.norm.dat <- readRDS(file = "NormPreprocess_rmse_norm.rds")

# (B) pAUC for norm methods for UMI data stratified full pnael DE setups, 384 vs 384, w/o spike-ins
pauc.norm.dat <- readRDS(file = "NormPreprocess_pauc_norm.rds") %>% dplyr::filter(Spike == "w/o Spike-Ins")

# (C) pAUC for norm methods for UMI data stratified 20% asymmetric, 384 vs 384, with spike-ins
pauc.normspike.dat <- readRDS(file = "NormPreprocess_pauc_normspike.rds") 

# (D) pAUC for preprocessing methods over MR and scran with groups for 20% asymmetric, 384 vs 384 for 10X genomics, SCRBseq, Smartseq2
pauc.preprocess.dat <-readRDS(file = "NormPreprocess_pauc_preprocess.rds")

# MAIN FIGURE -------------------------------------------------------------

# (A)  RMSE for norm methods for UMI data for 20 % asymmetric, 384 vs 384, w/o spike-ins
rmse.norm.sts <- rmse.norm.dat %>%
  dplyr::group_by(Normalisation) %>%
  dplyr::summarise_at(dplyr::vars(RMSE), 
                      dplyr::funs(Min = min(., na.rm=T), 
                                  Max = max(., na.rm=T), 
                                  Mean = mean(., na.rm=T),
                                  Median = median(., na.rm=T), 
                                  N = sum(!is.na(.)),
                                  LowHinge = boxplot.stats(.)[[1]][2],  
                                  UpHinge = boxplot.stats(.)[[1]][4]))
low.whisk <- min(rmse.norm.sts$LowHinge, na.rm=T) 
up.whisk <- max(rmse.norm.sts$UpHinge, na.rm =T)

sf.plot <- ggplot2::ggplot() +
  ggplot2::geom_rect(data = rmse.norm.dat,
                     xmin = -Inf, xmax = Inf,
                     ymin = -Inf,ymax = Inf,
                     size = 2, fill = NA, color = "yellow") +
  ggplot2::geom_boxplot(data = rmse.norm.dat,
                        aes(x = Normalisation, y = RMSE, color = Normalisation),
                        outlier.shape=NA, width = 0.75,
                        position=ggplot2::position_dodge(1)) +
  scale_color_manual(values = normalisation.cols) +
  ggplot2::scale_y_continuous(limits=c(0, up.whisk*1.05), 
                              breaks = c(0, 0.1, 0.2),
                              labels = scales::number_format(accuracy = 0.1)) +
  ggplot2::theme_light() + 
  ggplot2::scale_x_discrete(labels = norm_names) +
  ggplot2::labs(x = NULL,
                y = "RMSE", 
                title = "Deviance between estimated and \nsimulated library size factors (RMSE)") +
  ggplot2::geom_hline(yintercept = c(0), colour="darkgrey", linetype="dashed") + 
  ggplot2::guides(fill = guide_legend(nrow = 2)) + 
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

# (B)  pAUC for norm methods for UMI data stratified full panel DE setups, 384 vs 384, w/o spike-ins
sts.norm.nospike.pauc <- pauc.norm.dat %>%
  dplyr::group_by(perc.DE, LFC, Normalisation) %>%
  dplyr::summarise_at(dplyr::vars(pAUC), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
  dplyr::ungroup()

sts.norm.nospike.pauc<-sts.norm.nospike.pauc %>% 
  dplyr::filter(Normalisation=="SF") %>% 
  dplyr::select(-Normalisation) %>%
  dplyr::right_join( dplyr::filter(sts.norm.nospike.pauc, Normalisation !="SF"),
                     by=c("perc.DE", "LFC"),suffix=c(".sf","")) %>%  
  droplevels() %>%
  dplyr::mutate(x=as.numeric(Normalisation)) %>%
  data.frame()

nospike.pauc.plot <-  ggplot2::ggplot()+
  ggplot2::geom_rect(data = subset(sts.norm.nospike.pauc,LFC == "Asymmetric" & perc.DE == "0.2") ,
                     xmin = -Inf, xmax = Inf,
                     ymin = -Inf,ymax = Inf,
                     size = 2, fill = NA, color = "yellow") +
  ggplot2::geom_rect(data = subset(sts.norm.nospike.pauc,LFC == "Completely Asymmetric" & perc.DE == "0.6") ,
                     xmin = -Inf, xmax = Inf,
                     ymin = -Inf,ymax = Inf,
                     size = 2, fill = NA, color = "tomato") +
  ggplot2::geom_ribbon(data = sts.norm.nospike.pauc, 
                       aes(x=x, ymin=Lower.sf, ymax=Upper.sf),
                       col="darkgrey",fill="darkgrey")+
  ggplot2::geom_pointrange(data = sts.norm.nospike.pauc,
                           aes(x = x, 
                               y = Mean,
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Normalisation),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = normalisation.cols[-1]) +
  ggplot2::scale_x_continuous(labels = norm_names[levels(sts.norm.nospike.pauc$Normalisation)], 
                              breaks = unique(sts.norm.nospike.pauc$x),
                              minor_breaks = unique(sts.norm.nospike.pauc$x)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y= "pAUC",
                title = "Trade-off between power and false discoveries (pAUC)") +
  ggplot2::facet_grid(perc.DE ~ LFC, labeller = as_labeller(relabels)) +
  ggplot2::scale_y_continuous(limits=c(0, 1), 
                              breaks = c(0,0.25,0.5,0.75,1),
                              minor_breaks = c(0,0.25,0.5,0.75,1)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

# (C) pAUC for norm methods for UMI data stratified 20% asymmetric, 384 vs 384, w/o spike-ins
sts.norm.spike.pauc <- pauc.normspike.dat %>%
  dplyr::group_by( perc.DE, LFC, Normalisation) %>%
  dplyr::summarise_at(dplyr::vars(pAUC), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% ungroup()

sts.norm.spike.pauc<-sts.norm.spike.pauc %>% 
  dplyr::filter(Normalisation=="SF") %>% 
  dplyr::select(-Normalisation) %>%
  dplyr::right_join( dplyr::filter(sts.norm.spike.pauc, Normalisation !="SF"),
                     by=c("perc.DE", "LFC"),suffix=c(".sf","")) %>%  
  droplevels() %>%
  dplyr::mutate(x=as.numeric(Normalisation)) %>%
  data.frame()


spike.pauc.plot <- ggplot2::ggplot() +
  ggplot2::geom_rect(data = sts.norm.spike.pauc,
                     xmin = -Inf, xmax = Inf,
                     ymin = -Inf,ymax = Inf,
                     size = 2, fill = NA, color = "tomato") +
  ggplot2::geom_ribbon(data = sts.norm.spike.pauc, 
                       aes(x=x, ymin=Lower.sf, ymax=Upper.sf),
                       alpha=0.5)+
  ggplot2::geom_pointrange(data = sts.norm.spike.pauc,
                           aes(x = x, y = Mean,
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Normalisation),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = normalisation.cols) +
  ggplot2::scale_x_continuous(labels = norm_names[levels(sts.norm.spike.pauc$Normalisation)],
                              breaks = unique(sts.norm.spike.pauc$x),
                              minor_breaks = unique(sts.norm.spike.pauc$x)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y = "pAUC",
                title = "Using spike-ins (pAUC)") +
  ggplot2::scale_y_continuous(limits=c(0, 1), 
                              breaks = c(0,0.25,0.5,0.75,1),
                              minor_breaks = c(0,0.25,0.5,0.75,1)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

# (D) pAUC for preprocessing methods over protocols with selection of normalisation methods
sts.preprocess.pauc <- pauc.preprocess.dat %>%
  dplyr::group_by(Protocol, Normalisation, Preprocessing, Thinning) %>%
  dplyr::summarise_at(dplyr::vars(pAUC), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>%
  dplyr::ungroup() %>% 
  data.frame()

relabels["scran w/ groups"]<- "scran \n with groups"
relabels["scran w/ cluster"]<- "scran \n with cluster"
preprocess.pauc.plot <- ggplot2::ggplot(data = sts.preprocess.pauc, 
                                        aes(x = Preprocessing, 
                                            y = Mean, color=Preprocessing, 
                                            shape = Thinning)) +
  ggplot2::geom_pointrange(aes(x = Preprocessing, 
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Preprocessing,
                               shape = Thinning),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = preprocess.cols, guide = FALSE) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y = "pAUC",
                title = "Performance of normalisation methods using filtered or imputed counts (pAUC)") +
  ggplot2::facet_grid(Normalisation ~ Protocol, labeller = as_labeller(relabels)) +
  ggplot2::scale_y_continuous(limits=c(0, 1), breaks = c(0,0.25,0.5,0.75,1)) +
  ggplot2::scale_alpha_manual(values=c(1,0.25)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "bottom",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()


# combine the plots into one figure
upperpanel.plot <- cowplot::plot_grid(sf.plot + theme(legend.position = "none"),
                                      spike.pauc.plot  + theme(legend.position = "none"),
                                      labels = c("A", "C"), 
                                      label_size = 10,
                                      ncol = 2)

lowerpanel.plot<-cowplot::plot_grid(nospike.pauc.plot,
                                    preprocess.pauc.plot,
                                    align="hv",axis="lr",
                                    nrow = 2,
                                    ncol = 1,
                                    labels = c("B", "D"),
                                    label_size = 10,
                                    rel_heights = c(1.3,0.85))

norm.preprocess.combine.plot <- cowplot::plot_grid(upperpanel.plot,
                                                   lowerpanel.plot,
                                                   nrow = 2,
                                                   label_size = 10,
                                                   rel_heights = c(0.5,2.15))

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess", ".pdf"), 
                plot = norm.preprocess.combine.plot,
                width=210,
                height=297,
                units="mm")

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess", ".png"), 
                plot = norm.preprocess.combine.plot,
                width=210,
                height=297,
                units="mm")

# INPUT FOR SUPPLEMENT FIGURES --------------------------------------------



# RMSE for all settings (3*3 Panel), 384 vs 384 cells, alpha for spike-ins
# (A) Smartseq2
# (B) UMI methods
rmse.norm.si.dat <- readRDS(file = "NormPreprocess_rmse_norm_SI.rds")

# full panel of TPR, FDR and pAUC (3*3 Panel), alpha for spike-ins
# (A) Smartseq2
# (B) UMI methods
fdr.norm.si.dat <- readRDS(file = "NormPreprocess_fdr_norm_SI.rds")
tpr.norm.si.dat <- readRDS(file = "NormPreprocess_tpr_norm_SI.rds")
auc.norm.si.dat <- readRDS(file = "NormPreprocess_pauc_norm_SI.rds")

# fuller panel of pAUC matching main figure
# matching TPR and FDR of these pAUCs.
fdr.preprocess.si.dat <- readRDS(file = "NormPreprocess_fdr_preprocess_SI.rds")
tpr.preprocess.si.dat <- readRDS(file = "NormPreprocess_tpr_preprocess_SI.rds")
pauc.preprocess.si.dat <- readRDS(file = "NormPreprocess_pauc_preprocess_SI.rds")

# SUPPLEMENT FIGURES ------------------------------------------------------

##### RMSE #######

# RMSE for all settings (3*3 Panel), 384 vs 384 cells, alpha for spike-ins
# (A) Smartseq2
# (B) UMI methods
rmse.norm.sts <- rmse.norm.si.dat %>%
  dplyr::group_by(Normalisation) %>%
  dplyr::summarise_at(dplyr::vars(RMSE), 
                      dplyr::funs(Min = min(., na.rm=T), 
                                  Max = max(., na.rm=T), 
                                  Mean = mean(., na.rm=T),
                                  Median = median(., na.rm=T), 
                                  N = sum(!is.na(.)),
                                  LowHinge = boxplot.stats(.)[[1]][2],  
                                  UpHinge = boxplot.stats(.)[[1]][4]))
low.whisk <- min(rmse.norm.sts$LowHinge, na.rm=T) 
up.whisk <- max(rmse.norm.sts$UpHinge, na.rm =T)

sf.smartseq2.plot <- ggplot2::ggplot(data =  rmse.norm.si.dat %>% 
                                       dplyr::filter(Type == "Full-Length"),
                                     aes(x = Normalisation, 
                                         y = RMSE, 
                                         color = Normalisation)) +
  ggplot2::geom_boxplot(data =  rmse.norm.si.dat %>% 
                          dplyr::filter(Type == "Full-Length"),
                        aes(x = Normalisation, 
                            y = RMSE, 
                            color = Normalisation),
                        outlier.shape=NA, width = 0.75,
                        position=ggplot2::position_dodge(1)) +
  scale_color_manual(values = normalisation.cols) +
  ggplot2::scale_y_continuous(limits = c(0, 0.28),
                              breaks = c(0,0.1,0.2),
                              labels = scales::number_format(accuracy = 0.1)) +
  ggplot2::theme_light() + 
  ggplot2::scale_x_discrete(labels = norm_names) +
  ggplot2::facet_grid(perc.DE ~ LFC, 
                      labeller = as_labeller(relabels)) +
  ggplot2::labs(x = NULL,
                y = "RMSE", 
                title = "Deviance between estimated and simulated library size factors (RMSE) for Smart-seq2 data") +
  ggplot2::geom_hline(yintercept = c(0), colour="darkgrey", linetype="dashed") + 
  ggplot2::guides(fill = guide_legend(nrow = 2)) + 
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

sf.umi.plot <- ggplot2::ggplot(data =  rmse.norm.si.dat %>% 
                                 dplyr::filter(Type == "UMI") %>% 
                                 droplevels(),
                               aes(x = Normalisation, 
                                   y = RMSE, 
                                   color = Normalisation)) +
  ggplot2::geom_boxplot(data =  rmse.norm.si.dat %>% 
                          dplyr::filter(Type == "UMI") %>% 
                          droplevels(),
                        aes(x = Normalisation, 
                            y = RMSE, 
                            color = Normalisation),
                        outlier.shape=NA, width = 0.75,
                        position=ggplot2::position_dodge(1)) +
  scale_color_manual(values = normalisation.cols) +
  ggplot2::scale_y_continuous(limits = c(0, 0.28),
                              breaks = c(0,0.1,0.2),
                              labels = scales::number_format(accuracy = 0.1)) +
  ggplot2::theme_light() + 
  ggplot2::scale_x_discrete(labels = norm_names) +
  ggplot2::facet_grid(perc.DE ~ LFC, 
                      labeller = as_labeller(relabels)) +
  ggplot2::labs(x = NULL,
                y = "RMSE", 
                title = "Deviance between estimated and simulated library size factors (RMSE) for UMI data") +
  ggplot2::geom_hline(yintercept = c(0), colour="darkgrey", linetype="dashed") + 
  ggplot2::guides(fill = guide_legend(nrow = 2)) + 
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

rmse.si.combine.plot<-cowplot::plot_grid(sf.smartseq2.plot,
                                         sf.umi.plot,
                                         align="hv",axis="lr",
                                         nrow = 2,
                                         ncol = 1,
                                         labels = c("A", "B"),
                                         label_size = 10)

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess_RMSE_SI", ".pdf"), 
                plot = rmse.si.combine.plot,
                width=210,
                height=297,
                units="mm")

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess_RMSE_SI", ".png"), 
                plot = rmse.si.combine.plot,
                width=210,
                height=297,
                units="mm")

##### TPR #######
# full panel of TPR, FDR and pAUC (3*3 Panel), alpha for spike-ins
# (A) Smartseq2
# (B) UMI methods
sts.norm.tpr <- tpr.norm.si.dat %>%
  dplyr::group_by(perc.DE, LFC, Normalisation, Type, Spike) %>%
  dplyr::summarise_at(dplyr::vars(TPR), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!Normalisation=="SF")

sts.sf.tpr <- tpr.norm.si.dat %>%
  dplyr::filter(Normalisation=="SF") %>% 
  dplyr::select(-Normalisation) %>%
  dplyr::group_by(perc.DE, LFC, Type) %>%
  dplyr::summarise_at(dplyr::vars(TPR), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
  dplyr::ungroup()

sts.all.tpr <- dplyr::left_join(sts.norm.tpr, sts.sf.tpr,
                                by=c("perc.DE", "LFC", "Type"),
                                suffix=c("", ".sf")) %>%  
  droplevels() %>%
  dplyr::mutate(x=as.numeric(Normalisation)) %>%
  data.frame()


si.tpr.smartseq2.plot <-  ggplot2::ggplot()+
  ggplot2::geom_ribbon(data = sts.all.tpr %>% 
                         dplyr::filter(Type == "Full-Length") %>% 
                         droplevels(), 
                       aes(x=x, ymin=Lower.sf, ymax=Upper.sf),
                       col="darkgrey",fill="darkgrey")+
  ggplot2::geom_pointrange(data = sts.all.tpr %>% 
                             dplyr::filter(Type == "Full-Length") %>% 
                             droplevels(),
                           aes(x = x, 
                               y = Mean,
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Normalisation,
                               alpha=Spike),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = normalisation.cols[-1]) +
  ggplot2::scale_alpha_manual(values=c(1,0.25)) +
  ggplot2::scale_y_continuous(limits = c(0.5,1),
                              breaks = c(0.5, 0.75, 1), 
                              labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_continuous(labels = norm_names[levels(sts.all.tpr$Normalisation)], 
                              breaks = unique(sts.all.tpr$x),
                              minor_breaks = unique(sts.all.tpr$x)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y= "TPR",
                title = "Power to detect true differential expression (TPR)  for Smart-seq2 data") +
  ggplot2::facet_grid(perc.DE ~ LFC, labeller = as_labeller(relabels)) +
  # ggplot2::scale_y_continuous(limits=c(0, 1), 
  #                             breaks = c(0,0.25,0.5,0.75,1),
  #                             minor_breaks = c(0,0.25,0.5,0.75,1)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

si.tpr.umi.plot <-  ggplot2::ggplot()+
  ggplot2::geom_ribbon(data = sts.all.tpr %>% 
                         dplyr::filter(Type == "UMI"), 
                       aes(x=x, ymin=Lower.sf, ymax=Upper.sf),
                       col="darkgrey",fill="darkgrey")+
  ggplot2::geom_pointrange(data = sts.all.tpr %>% 
                             dplyr::filter(Type == "UMI"),
                           aes(x = x, 
                               y = Mean,
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Normalisation,
                               alpha=Spike),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = normalisation.cols[-1]) +
  ggplot2::scale_alpha_manual(values=c(1,0.25)) +
  ggplot2::scale_y_continuous(limits = c(0.5,1),
                              breaks = c(0.5, 0.75, 1), 
                              labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_continuous(labels = norm_names[levels(sts.all.tpr$Normalisation)], 
                              breaks = unique(sts.all.tpr$x),
                              minor_breaks = unique(sts.all.tpr$x)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y= "TPR",
                title = "Power to detect true differential expression (TPR)  for UMI data") +
  ggplot2::facet_grid(perc.DE ~ LFC, labeller = as_labeller(relabels)) +
  # ggplot2::scale_y_continuous(limits=c(0, 1), 
  #                             breaks = c(0,0.25,0.5,0.75,1),
  #                             minor_breaks = c(0,0.25,0.5,0.75,1)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

tpr.si.combine.plot<-cowplot::plot_grid(si.tpr.smartseq2.plot,
                                        si.tpr.umi.plot,
                                        align="hv",axis="lr",
                                        nrow = 2,
                                        ncol = 1,
                                        labels = c("A", "B"),
                                        label_size = 10)

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess_TPR_SI", ".pdf"), 
                plot = tpr.si.combine.plot,
                width=210,
                height=297,
                units="mm")

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess_TPR_SI", ".png"), 
                plot = tpr.si.combine.plot,
                width=210,
                height=297,
                units="mm")


#### FDR #####

sts.norm.fdr <- fdr.norm.si.dat %>%
  dplyr::group_by(perc.DE, LFC, Normalisation, Type, Spike) %>%
  dplyr::summarise_at(dplyr::vars(FDR), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!Normalisation=="SF")

sts.sf.fdr <- fdr.norm.si.dat %>%
  dplyr::filter(Normalisation=="SF") %>% 
  dplyr::select(-Normalisation) %>%
  dplyr::group_by(perc.DE, LFC, Type) %>%
  dplyr::summarise_at(dplyr::vars(FDR), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
  dplyr::ungroup()

sts.all.fdr <- dplyr::left_join(sts.norm.fdr, sts.sf.fdr,
                                by=c("perc.DE", "LFC", "Type"),
                                suffix=c("", ".sf")) %>%  
  droplevels() %>%
  dplyr::mutate(x=as.numeric(Normalisation)) %>%
  data.frame()


si.fdr.smartseq2.plot <-  ggplot2::ggplot()+
  ggplot2::geom_ribbon(data = sts.all.fdr %>% 
                         dplyr::filter(Type == "Full-Length") %>% 
                         droplevels(), 
                       aes(x=x, ymin=Lower.sf, ymax=Upper.sf),
                       col="darkgrey",fill="darkgrey")+
  ggplot2::geom_pointrange(data = sts.all.fdr %>% 
                             dplyr::filter(Type == "Full-Length") %>% 
                             droplevels(),
                           aes(x = x, 
                               y = Mean,
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Normalisation,
                               alpha=Spike),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = normalisation.cols[-1]) +
  ggplot2::scale_alpha_manual(values=c(1,0.25)) +
  ggplot2::geom_hline(yintercept = c(0.1), colour="darkgrey", linetype="dashed") + 
  ggplot2::scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4), 
                              labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_continuous(labels = norm_names[levels(sts.all.fdr$Normalisation)], 
                              breaks = unique(sts.all.fdr$x),
                              minor_breaks = unique(sts.all.fdr$x)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y= "FDR",
                title = "Control over false discoveries (FDR)  for Smart-seq2 data") +
  ggplot2::facet_grid(perc.DE ~ LFC, labeller = as_labeller(relabels), scales = "free") +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

si.fdr.umi.plot <-  ggplot2::ggplot()+
  ggplot2::geom_ribbon(data = sts.all.fdr %>% 
                         dplyr::filter(Type == "UMI") %>% 
                         droplevels(), 
                       aes(x=x, ymin=Lower.sf, ymax=Upper.sf),
                       col="darkgrey",fill="darkgrey")+
  ggplot2::geom_pointrange(data = sts.all.fdr %>% 
                             dplyr::filter(Type == "UMI"),
                           aes(x = x, 
                               y = Mean,
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Normalisation,
                               alpha=Spike),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = normalisation.cols[-1]) +
  ggplot2::scale_alpha_manual(values=c(1,0.25)) +
  ggplot2::geom_hline(yintercept = c(0.1), colour="darkgrey", linetype="dashed") + 
  ggplot2::scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4), 
                              labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_continuous(labels = norm_names[levels(sts.all.fdr$Normalisation)], 
                              breaks = unique(sts.all.fdr$x),
                              minor_breaks = unique(sts.all.fdr$x)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y= "FDR",
                title = "Control over false discoveries (FDR)  for UMI data") +
  ggplot2::facet_grid(perc.DE ~ LFC, labeller = as_labeller(relabels), scales = "free") +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()


fdr.si.combine.plot<-cowplot::plot_grid(si.fdr.smartseq2.plot,
                                        si.fdr.umi.plot,
                                        align="hv",axis="lr",
                                        nrow = 2,
                                        ncol = 1,
                                        labels = c("A", "B"),
                                        label_size = 10)

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess_FDR_SI", ".pdf"), 
                plot = fdr.si.combine.plot,
                width=210,
                height=297,
                units="mm")

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess_FDR_SI", ".png"), 
                plot = fdr.si.combine.plot,
                width=210,
                height=297,
                units="mm")

##### pAUC ######

sts.norm.pauc <- auc.norm.si.dat %>%
  dplyr::group_by(perc.DE, LFC, Normalisation, Type, Spike) %>%
  dplyr::summarise_at(dplyr::vars(pAUC), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!Normalisation=="SF")

sts.sf.pauc <- auc.norm.si.dat %>%
  dplyr::filter(Normalisation=="SF") %>% 
  dplyr::select(-Normalisation) %>%
  dplyr::group_by(perc.DE, LFC, Type) %>%
  dplyr::summarise_at(dplyr::vars(pAUC), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
  dplyr::ungroup()

sts.all.pauc<- dplyr::left_join(sts.norm.pauc, sts.sf.pauc,
                                by=c("perc.DE", "LFC", "Type"),
                                suffix=c("", ".sf")) %>%  
  droplevels() %>%
  dplyr::mutate(x=as.numeric(Normalisation)) %>%
  data.frame()


si.pauc.smartseq2.plot <-  ggplot2::ggplot()+
  ggplot2::geom_ribbon(data = sts.all.pauc %>% 
                         dplyr::filter(Type == "Full-Length") %>% 
                         droplevels(), 
                       aes(x=x, ymin=Lower.sf, ymax=Upper.sf),
                       col="darkgrey",fill="darkgrey")+
  ggplot2::geom_pointrange(data = sts.all.pauc %>% 
                             dplyr::filter(Type == "Full-Length") %>% 
                             droplevels(),
                           aes(x = x, 
                               y = Mean,
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Normalisation,
                               alpha=Spike),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = normalisation.cols[-1]) +
  ggplot2::scale_alpha_manual(values=c(1,0.25)) +
  ggplot2::scale_x_continuous(labels = norm_names[levels(sts.all.pauc$Normalisation)], 
                              breaks = unique(sts.all.pauc$x),
                              minor_breaks = unique(sts.all.pauc$x)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y= "pAUC",
                title = "Trade-off between power and false discoveries (pAUC)  for Smart-seq2 data") +
  ggplot2::facet_grid(perc.DE ~ LFC, labeller = as_labeller(relabels)) +
  ggplot2::scale_y_continuous(limits=c(0, 1), 
                              breaks = c(0,0.25,0.5,0.75,1),
                              minor_breaks = c(0,0.25,0.5,0.75,1)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

si.pauc.umi.plot <-  ggplot2::ggplot()+
  ggplot2::geom_ribbon(data = sts.all.pauc %>% 
                         dplyr::filter(Type == "UMI"), 
                       aes(x=x, ymin=Lower.sf, ymax=Upper.sf),
                       col="darkgrey",fill="darkgrey")+
  ggplot2::geom_pointrange(data = sts.all.pauc %>% 
                             dplyr::filter(Type == "UMI"),
                           aes(x = x, 
                               y = Mean,
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Normalisation,
                               alpha=Spike),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = normalisation.cols[-1]) +
  ggplot2::scale_alpha_manual(values=c(1,0.25)) +
  ggplot2::scale_x_continuous(labels = norm_names[levels(sts.all.pauc$Normalisation)], 
                              breaks = unique(sts.all.pauc$x),
                              minor_breaks = unique(sts.all.pauc$x)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y= "pAUC",
                title = "Trade-off between power and false discoveries (pAUC)  for UMI data") +
  ggplot2::facet_grid(perc.DE ~ LFC, labeller = as_labeller(relabels)) +
  ggplot2::scale_y_continuous(limits=c(0, 1), 
                              breaks = c(0,0.25,0.5,0.75,1),
                              minor_breaks = c(0,0.25,0.5,0.75,1)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()


pauc.si.combine.plot<-cowplot::plot_grid(si.pauc.smartseq2.plot,
                                         si.pauc.umi.plot,
                                         align="hv",axis="lr",
                                         nrow = 2,
                                         ncol = 1,
                                         labels = c("A", "B"),
                                         label_size = 10)

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess_pAUC_SI", ".pdf"), 
                plot = pauc.si.combine.plot,
                width=210,
                height=297,
                units="mm")

ggplot2::ggsave(filename=paste0("Performance_Norm+Preprocess_pAUC_SI", ".png"), 
                plot = pauc.si.combine.plot,
                width=210,
                height=297,
                units="mm")


##### PREPROCESS  ######

sts.preprocess.pauc <- pauc.preprocess.si.dat %>%
  dplyr::group_by(Protocol, Normalisation, Preprocessing) %>%
  dplyr::summarise_at(dplyr::vars(pAUC), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>%
  dplyr::ungroup() %>% 
  data.frame()

preprocess.pauc.si.plot <- ggplot2::ggplot(data = sts.preprocess.pauc, 
                                           aes(x = Preprocessing,
                                               y = Mean, 
                                               color = Preprocessing)) +
  ggplot2::geom_pointrange(aes(x = Preprocessing, 
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Preprocessing),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = preprocess.cols) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y = "pAUC",
                title = "Performance of normalisation methods using filtered or imputed counts (pAUC)") +
  ggplot2::facet_grid(Normalisation ~ Protocol, labeller = as_labeller(relabels)) +
  ggplot2::scale_y_continuous(limits=c(0, 1), breaks = c(0,0.25,0.5,0.75,1)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()

ggplot2::ggsave(filename=paste0("Performance_Preprocess_pAUC_SI", ".pdf"), 
                plot = preprocess.pauc.si.plot,
                width=180,
                height=210,
                units="mm")

ggplot2::ggsave(filename=paste0("Performance_Preprocess_pAUC_SI", ".png"), 
                plot = preprocess.pauc.si.plot,
                width=180,
                height=210,
                units="mm")


sts.preprocess.tpr <- tpr.preprocess.si.dat %>%
  dplyr::group_by(Protocol, Normalisation, Preprocessing) %>%
  dplyr::summarise_at(dplyr::vars(TPR), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!Normalisation=="SF")


preprocess.tpr.si.plot <-  ggplot2::ggplot(data = sts.preprocess.tpr, 
                                           aes(x = Preprocessing,
                                               y = Mean, 
                                               color = Preprocessing)) +
  ggplot2::geom_pointrange(aes(x = Preprocessing, 
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Preprocessing),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = preprocess.cols) +
  ggplot2::scale_y_continuous(limits=c(0.5, 1), 
                              breaks = c(0.5,0.6,0.7,0.8,0.9,1), 
                              labels = scales::percent_format(accuracy = 1)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y = "TPR",
                title = "Power of normalisation methods using filtered or imputed counts") +
  ggplot2::facet_grid(Normalisation ~ Protocol, labeller = as_labeller(relabels)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()


ggplot2::ggsave(filename=paste0("Performance_Preprocess_TPR_SI", ".pdf"), 
                plot = preprocess.tpr.si.plot,
                width=180,
                height=210,
                units="mm")

ggplot2::ggsave(filename=paste0("Performance_Preprocess_TPR_SI", ".png"), 
                plot = preprocess.tpr.si.plot,
                width=180,
                height=210,
                units="mm")


sts.preprocess.fdr <- fdr.preprocess.si.dat %>%
  dplyr::group_by(Protocol, Normalisation, Preprocessing) %>%
  dplyr::summarise_at(dplyr::vars(FDR), 
                      dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                  Lower = Hmisc::smean.cl.boot(.)[2],
                                  Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!Normalisation=="SF")


preprocess.fdr.si.plot <-  ggplot2::ggplot(data = sts.preprocess.fdr, 
                                           aes(x = Preprocessing,
                                               y = Mean, 
                                               color = Preprocessing)) +
  ggplot2::geom_pointrange(aes(x = Preprocessing, 
                               ymin = Lower, 
                               ymax = Upper, 
                               color=Preprocessing),
                           size = 0.25, position = position_dodge(width = 1)) +
  ggplot2::scale_color_manual(values = preprocess.cols) +
  ggplot2::geom_hline(yintercept = c(0.1), colour="darkgrey", linetype="dashed") + 
  ggplot2::scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4), 
                              labels = scales::percent_format(accuracy = 1)) +
  ggplot2::theme_light() + 
  ggplot2::labs(x = NULL,
                y= "FDR",
                title = "FDR control of normalisation methods using filtered or imputed counts") +
  ggplot2::facet_grid(Normalisation ~ Protocol, labeller = as_labeller(relabels)) +
  ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                 legend.position = "none",
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = grid::unit(1.5, "lines"),
                 axis.text.x=ggplot2::element_text(size=8, color='black'),
                 axis.text.y=ggplot2::element_text(size=8, color='black'),
                 axis.title=ggplot2::element_text(size=8, face="bold", color='black'),
                 plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                 strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                 strip.background = ggplot2::element_rect(fill="white")) + 
  ggplot2::coord_flip()


ggplot2::ggsave(filename=paste0("Performance_Preprocess_FDR_SI", ".pdf"), 
                plot = preprocess.fdr.si.plot,
                width=180,
                height=210,
                units="mm")

ggplot2::ggsave(filename=paste0("Performance_Preprocess_FDR_SI", ".png"), 
                plot = preprocess.fdr.si.plot,
                width=180,
                height=210,
                units="mm")

