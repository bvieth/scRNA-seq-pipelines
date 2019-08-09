
# SETTINGS ----------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(xtable)
library(magick)
# library(powsimRDev)

# plotting colors
source(file  = "colour_schemes_paper.R")
# renaming variables for plotting
source(file  = "renaming_paper.R")

generateLabels <- function(type) {
  if(type == "NormType") {
    ModelComparison <- c("~ Protocol + Mapper + Annotation + Preprocessing + NormSpike + DEtool vs ~ 1", 
                         "~ Mapper + Annotation + Preprocessing + NormSpike + DEtool vs ~ 1", 
                         "~ Protocol + Annotation + Preprocessing + NormSpike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Preprocessing + NormSpike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + NormSpike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Preprocessing + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Preprocessing + NormSpike vs ~ 1", 
                         "Protocol", 
                         "Mapper", 
                         "Annotation", 
                         "Preprocessing", 
                         "NormSpike", 
                         "DEtool")
    
    ModelComparison <- factor(ModelComparison, levels = ModelComparison)
    
    ModelName <- c("Full Model", 
                   "Leave out Protocol", 
                   "Leave out Aligner", 
                   "Leave out Annotation", 
                   "Leave out Preprocessing", 
                   "Leave out NormSpike",
                   "Leave out DE-Tool",
                   "Protocol", 
                   "Aligner", 
                   "Annotation", 
                   "Preprocessing",
                   "Spike",
                   "DE-Tool")
    ModelName <- factor(ModelName, levels = ModelName)                     
    
    ModelFitType <- c(rep("Leave-one-out models", 7),
                      rep("Leave-one-out cross-validation", 6))
    ModelFitType <- factor(ModelFitType, levels = unique(ModelFitType))
    
    ModelAnnot <- data.frame(ModelComparison, ModelFitType, ModelName, stringsAsFactors = F)
  }
  if(type == "Norm+Norm:Spike") {
    ModelComparison <- c("~ Protocol + Mapper + Annotation + Preprocessing + Normalisation + Normalisation:Spike + DEtool vs ~ 1", 
                         "~ Mapper + Annotation + Preprocessing + Normalisation + Normalisation:Spike + DEtool vs ~ 1", 
                         "~ Protocol + Annotation + Preprocessing + Normalisation + Normalisation:Spike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Preprocessing + Normalisation + Normalisation:Spike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Normalisation + Normalisation:Spike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Preprocessing + Normalisation:Spike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Preprocessing + Normalisation + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Preprocessing + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Preprocessing + Normalisation + Normalisation:Spike vs ~ 1", 
                         "Protocol", 
                         "Mapper", 
                         "Annotation", 
                         "Preprocessing", 
                         "Normalisation", 
                         "Normalisation:Spike", 
                         "Normalisation+Normalisation:Spike", 
                         "DEtool")
    
    ModelComparison <- factor(ModelComparison, levels = ModelComparison)
    
    ModelName <- c("Full Model", 
                   "Leave out Protocol", 
                   "Leave out Aligner", 
                   "Leave out Annotation", 
                   "Leave out Preprocessing", 
                   "Leave out Normalisation",
                   "Leave out Normalisation:Spike", 
                   "Leave out Normalisation+Normalisation:Spike",
                   "Leave out DE-Tool",
                   "Protocol", 
                   "Aligner", 
                   "Annotation", 
                   "Preprocessing", 
                   "Normalisation", 
                   "Spike",
                   "Normalisation+Normalisation:Spike",
                   "DE-Tool")
    ModelName <- factor(ModelName, levels = ModelName)                     
    
    ModelFitType <- c(rep("Leave-one-out models", 9),
                      rep("Leave-one-out cross-validation", 8))
    ModelFitType <- factor(ModelFitType, levels = unique(ModelFitType))
    
    ModelAnnot <- data.frame(ModelComparison, ModelFitType, ModelName, stringsAsFactors = F)
  }
  if(type == "Norm+Spike") {
    ModelComparison <- c("~ Protocol + Mapper + Annotation + Preprocessing + Normalisation + Spike + DEtool vs ~ 1", 
                         "~ Mapper + Annotation + Preprocessing + Normalisation + Spike + DEtool vs ~ 1", 
                         "~ Protocol + Annotation + Preprocessing + Normalisation + Spike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Preprocessing + Normalisation + Spike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Normalisation + Spike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Preprocessing + Spike + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Preprocessing + Normalisation + DEtool vs ~ 1", 
                         "~ Protocol + Mapper + Annotation + Preprocessing + Normalisation + Spike vs ~ 1", 
                         "Protocol", 
                         "Mapper", 
                         "Annotation",
                         "Preprocessing", 
                         "Normalisation", 
                         "Spike", 
                         "DEtool")
    
    ModelComparison <- factor(ModelComparison, levels = ModelComparison)
    
    ModelName <- c("Full Model", 
                   "Leave out Protocol", 
                   "Leave out Aligner", 
                   "Leave out Annotation", 
                   "Leave out Preprocessing", 
                   "Leave out Normalisation",
                   "Leave out Spike", 
                   "Leave out DE-Tool",
                   "Protocol", 
                   "Aligner", 
                   "Annotation", 
                   "Preprocessing", 
                   "Normalisation", 
                   "Spike",
                   "DE-Tool")
    ModelName <- factor(ModelName, levels = ModelName)                     
    
    ModelFitType <- c(rep("Leave-one-out models", 8),
                      rep("Leave-one-out cross-validation", 7))
    ModelFitType <- factor(ModelFitType, levels = unique(ModelFitType))
    
    ModelAnnot <- data.frame(ModelComparison, ModelFitType, ModelName, stringsAsFactors = F)
  }
  return(ModelAnnot)
}

# BETA REGRESSION FOR NORM + NORM:SPIKE -----------------------------------

modelfit = "Norm+Norm:Spike"
average = TRUE
tenx = TRUE
measure = "MCC_lib"
samplesize = "384 vs 384"

aucdattype <- ifelse(isTRUE(average), "-avgauc.rds$", "-rawauc.rds$")
reflabel <- ifelse(isTRUE(tenx), "-10XGenomics", "-Smartseq2")
filext <- paste0(reflabel, aucdattype)

ModelAnnot <- generateLabels(type = modelfit)

# get the appropiate performance measure: MCC_lib
mcc <- readRDS(file = "MCC_lib.rds")
# setup annotation info
evaljobs.dat <- readRDS(file = "evalfulldat.rds")
umi.info <- data.frame(Protocol = c("10XGenomics", "CELseq2", "Dropseq", "SCRBseq", "Smartseq2"),
                       Type = c("UMI", "UMI", "UMI", "UMI", "Full-Length"), 
                       stringsAsFactors = F)
# add UMI vs nonUMI method column
evaljobs.dat <- evaljobs.dat %>%
  dplyr::left_join(umi.info, by="Protocol") 

# performance data must be in the format Name, SampleSize, Value
perfdat <- mcc %>% 
  dplyr::rename(Value = MCC)
perfdat$SampleSize <- factor(perfdat$SampleSize, 
                             levels = c("96 vs 96", "50 vs 200", "384 vs 384"))

favorite <- c("SCRBseq", "star", "gencode", 
              "SAVER", "scran w/ groups", "w/o Spike-Ins", "limma-trend")

naive <- c("SCRBseq", "bwa", "gencode", 
           "none", "MR", "w/o Spike-Ins", "T-Test")

perflabel <- "MCC"

plot.pseudoR.recommend <- function(measure, 
                                   samplesize, 
                                   tenx, 
                                   average, 
                                   modelfit,
                                   favorite, 
                                   naive, 
                                   perfdat, 
                                   perflabel) {
  
  # read in beta regression result files
  fitres.files <- list.files(path = paste0("beta_regression/", 
                                           measure, "/", modelfit, "/"), 
                             pattern=filext, full.names = T, include.dirs = T)
  fitres.names <-  substring(list.files(paste0("beta_regression/", 
                                               measure, "/", modelfit, "/"), 
                                        pattern=filext, full=FALSE), 
                             first = 1, 
                             last = nchar(list.files(paste0("beta_regression/", 
                                                            measure, "/", modelfit, "/"),  
                                                     pattern=filext, full=FALSE))-nchar(filext)+1)
  names(fitres.files) <- fitres.names
  
  # read in performance measure files
  res.fit.pseudoR.L <- lapply(1:length(fitres.files), function(f) {
    tmp <- readRDS(file=fitres.files[f])$Test
    tmp2 <- sapply(names(tmp), function(i) {
      tmp[[i]]$PseudoR
    }, USE.NAMES = TRUE, simplify = FALSE)
    tmp2
  })
  names(res.fit.pseudoR.L) <- fitres.names
  res.fit.pseudoR.L2 <- res.fit.pseudoR.L[!unlist(lapply(res.fit.pseudoR.L, is.null))]
  res.fit.pseudoR <- data.table::rbindlist(res.fit.pseudoR.L2, use.names = T, fill=T, idcol = 'Name')
  
  # process LRT results
  res.fit.lrt.L <- lapply(1:length(fitres.files), function(f) {
    tmp <- readRDS(file=fitres.files[f])$Test
    tmp2 <- sapply(names(tmp), function(i) {
      tmp[[i]]$LRT
    }, USE.NAMES = TRUE, simplify = FALSE)
    data.table::rbindlist(tmp2, use.names = T, fill=T, idcol = 'ModelComparison')
  })
  names(res.fit.lrt.L) <- fitres.names
  res.fit.lrt.L2 <- res.fit.lrt.L[!unlist(lapply(res.fit.lrt.L, is.null))]
  res.fit.lrt <- data.table::rbindlist(res.fit.lrt.L2, use.names = T, fill=T, idcol = 'Name') %>% 
    dplyr::select(Name, ModelComparison, P.Value) %>% 
    data.table::dcast(Name ~ ModelComparison, value.var = "P.Value") %>% 
    dplyr::mutate(Normalisation = `Normalisation+Normalisation:Spike`) %>%
    data.table::melt(id.vars = c("Name"),
                     variable.name = "ModelComparison", value.name = "P.Value",
                     variable.factor = FALSE)
  
  # process pseudoR and LRT values
  res.fit.processed <- res.fit.pseudoR %>% 
    dplyr::mutate(Normalisation = `Normalisation+Normalisation:Spike` - `Normalisation:Spike`) %>% 
    tidyr::gather(key = ModelComparison, value = Value, -Name)  %>%
    dplyr::rename(PseudoR=Value) %>%
    tidyr::separate(col=Name, into=c("SampleSize", "DEsetup"), sep="-", remove=F) %>%
    tidyr::separate(col=DEsetup, into=c("percDE", "LFC"), sep="_", remove=F) %>% 
    dplyr::left_join(res.fit.lrt, by = c("Name"="Name", "ModelComparison"="ModelComparison")) %>% 
    dplyr::left_join(ModelAnnot, by = "ModelComparison") %>%
    dplyr::mutate(ModelComparison = factor(ModelComparison, 
                                           levels = levels(ModelAnnot$ModelComparison)),
                  ModelFitType = factor(ModelFitType, 
                                        levels = levels(ModelAnnot$ModelFitType))) %>% 
    dplyr::group_by(SampleSize, percDE, LFC, ModelFitType) %>% 
    dplyr::mutate(adj.P.Value = p.adjust(P.Value, method = "bonferroni")) %>% 
    dplyr::mutate(Significant = ifelse(adj.P.Value <= 0.01, "p < 0.01", "p > 0.01")) %>% 
    dplyr::mutate(PseudoR = ifelse(PseudoR < 0, 0, PseudoR)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(ModelFitType = ifelse(ModelName == "Full Model", 
                                        "Leave-one-out cross-validation", 
                                        as.character(ModelFitType))) %>% 
    dplyr::filter(ModelFitType == "Leave-one-out cross-validation") %>% 
    dplyr::filter(ModelName != "Normalisation+Normalisation:Spike")
  
  res.fit.processed$LFC <-  dplyr::recode(res.fit.processed$LFC, lfc.gamma.narrow.sym = "Symmetric", 
                                          lfc.gamma.narrow.asym = "Asymmetric", 
                                          lfc.gamma.narrow.casym = "Completely Asymmetric")
  res.fit.processed$SampleSize <- factor(res.fit.processed$SampleSize, 
                                         levels=c("96 vs 96", "384 vs 384", "50 vs 200"))
  res.fit.processed$LFC <- factor(res.fit.processed$LFC, 
                                  levels=c("Symmetric", "Asymmetric", "Completely Asymmetric"))
  res.fit.processed$percDE <- factor(res.fit.processed$percDE, 
                                     levels=c("0.05", "0.2", "0.6"))
  
  ## plot data
  res.fit.pdat <- res.fit.processed %>% 
    dplyr::filter( SampleSize == samplesize)

  a <- c("0.05", "Symmetric")
  b <- c("0.2", "Symmetric")
  d <- c("0.2", "Asymmetric")
  e <- c("0.6", "Asymmetric")
  
  examplesetups <- list("5% Symmetric DE" = a,
                        "20% Symmetric DE" = b,
                        "20% Asymmetric DE" = d,
                        "60% Asymmetric DE" = e)
  
  alldat.plot.l <- setNames(vector("list", length(examplesetups)), names(examplesetups))
  for(setup in names(examplesetups)) {
    res.fit.dat.plot <- res.fit.pdat %>% 
      dplyr::filter(percDE == examplesetups[[setup]][1] &
                      LFC == examplesetups[[setup]][2] &
                      ModelName %in% c("Full Model", 
                                       "Protocol", 
                                       "Aligner", 
                                       "Annotation", 
                                       "Preprocessing", 
                                       "Normalisation",
                                       "Spike",
                                       "DE-Tool")) %>% 
      droplevels() %>% 
      dplyr::mutate(PseudoR = ifelse(PseudoR<0, 0, PseudoR))
    error.model <- data.frame(res.fit.dat.plot[1,], stringsAsFactors = F)
    error.model$PseudoR <- as.numeric(abs(res.fit.dat.plot[grepl("Full Model", res.fit.dat.plot$ModelName), "PseudoR"] - 
                                            sum(res.fit.dat.plot[!grepl("Full Model", res.fit.dat.plot$ModelName), "PseudoR"])))
    error.model$ModelName <- "Interaction"
    error.model$ModelFitType <- "Leave-one-out cross-validation"
    error.model$Significant <- "p > 0.01"
    unexplained.model <- data.frame(res.fit.dat.plot[grepl("Full Model", res.fit.dat.plot$ModelName),], stringsAsFactors = F)
    unexplained.model$PseudoR <- 1 - as.numeric(unexplained.model$PseudoR) - as.numeric(error.model$PseudoR)
    unexplained.model$ModelName <- "Unexplained"
    unexplained.model$ModelFitType <- "Leave-one-out cross-validation"
    unexplained.model$Significant <- "p > 0.01"
    
    dat.plot <- res.fit.dat.plot %>% 
      dplyr::filter(!ModelName == "Full Model") %>% 
      dplyr::bind_rows(error.model) %>% 
      dplyr::bind_rows(unexplained.model) %>%
      dplyr::mutate(RelContribute = PseudoR / sum(PseudoR))
    dat.plot$ModelName <- factor(dat.plot$ModelName, levels = c("Protocol", 
                                                                "Aligner", 
                                                                "Annotation", 
                                                                "Preprocessing", 
                                                                "Normalisation",
                                                                "Spike", 
                                                                "DE-Tool",
                                                                "Interaction",
                                                                "Unexplained"))
    dat.plot$ModelName <- forcats::fct_rev(dat.plot$ModelName)
    alldat.plot.l[[setup]] <- dat.plot
    
  }
  
  alldat.plot <- data.table::rbindlist(alldat.plot.l, use.names = T, fill=T, idcol = 'DEPattern')
  
  
  cols <- c("Protocol" = "#42c4c4", 
            "Aligner" = "khaki3", 
            "Annotation" = "khaki1", 
            "Preprocessing" = "#ab6e39", 
            "Normalisation" = "#F25151",
            "Spike" = "#F97D7D", 
            "DE-Tool" = "#4daf4a",
            "Interaction" = "gray50",
            "Unexplained" = "gray80")
  
  alldat.plot$DEPattern <- factor(alldat.plot$DEPattern, levels = c("5% Symmetric DE",
                                                                    "20% Symmetric DE",
                                                                    "20% Asymmetric DE",
                                                                    "60% Asymmetric DE"))
  
  res.pseudoR.plot <- ggplot2::ggplot(data = alldat.plot,
                                      ggplot2::aes(x = ModelFitType, 
                                                   y = RelContribute, 
                                                   fill = ModelName)) +
    ggplot2::geom_col(color = "black") + 
    ggplot2::scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1), 
                                limits = c(0, 1.001),
                                labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::guides(fill = guide_legend(nrow = 2, reverse=T, order = 1)) + 
    ggplot2::labs(x = NULL, y = perflabel,
                  title = bquote(bold('  Relative Contribution to R '^2))) +
    ggplot2::theme_light() +
    ggplot2::facet_wrap(~DEPattern, ncol = 1, nrow = 4) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                   legend.position = "bottom",
                   legend.box = "vertical",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(0.5, "lines"),
                   axis.text.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size=8, color='black'),
                   axis.title.x = ggplot2::element_blank(),
                   plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white")) +
    ggplot2::coord_flip()
  
  # panel (A): beta regression pseudo R-squared
  final.pseudoR <- res.pseudoR.plot
  
  # panel (B): performance measure of naive and favorite setup across sample sizes and DE settings
  jobs2eval.favorite <- evaljobs.dat %>% dplyr::filter(Protocol == favorite[1] &
                                                         Mapper == favorite[2] &
                                                         Annotation ==  favorite[3] &
                                                         Preprocessing == favorite[4] &
                                                         Normalisation  == favorite[5] &
                                                         Spike == favorite[6] &
                                                         DEtool == favorite[7] &
                                                         DEFilter == "raw") %>% 
    dplyr::mutate(Pipeline ="Good")
  jobs2eval.naive <- evaljobs.dat %>% dplyr::filter(Protocol == naive[1] &
                                                      Mapper == naive[2] &
                                                      Annotation ==  naive[3] &
                                                      Preprocessing == naive[4] &
                                                      Normalisation  == naive[5] &
                                                      Spike == naive[6] &
                                                      DEtool == naive[7] &
                                                      DEFilter == "raw") %>% 
    dplyr::mutate(Pipeline ="Naive")
  jobs2eval <- dplyr::bind_rows(jobs2eval.favorite, jobs2eval.naive) %>% 
    droplevels()
  
  pdat <- jobs2eval %>% 
    dplyr::left_join(perfdat, by = "Name") %>% 
    tidyr::replace_na(list(Value = 0)) %>% 
    dplyr::filter(! perc.DE == "0.00")
  
  
  sts.dat <- pdat %>%
    dplyr::group_by(perc.DE, LFC, SampleSize, Pipeline,
                    Protocol, Mapper, Annotation, Preprocessing, 
                    Normalisation, Spike, DEtool) %>%
    dplyr::summarise_at(dplyr::vars(Value), 
                        dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                    Lower = Hmisc::smean.cl.boot(.)[2],
                                    Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
    ungroup()
  
  pipeline.cols <- c("Naive" = "purple", 
                     "Good" = "orange")
  
  pdat <- sts.dat %>% 
    tidyr::unite("DEPattern", c("perc.DE", "LFC"), sep="_", remove = F)  %>% 
    dplyr::filter(DEPattern %in% c("0.05_Asymmetric", "0.05_Symmetric", "0.2_Asymmetric", "0.6_Asymmetric"))
  
  res.recommend.plot <- ggplot2::ggplot(data = pdat, 
                                        aes(x = SampleSize, 
                                            y = Mean, 
                                            color = Pipeline)) +
    ggplot2::geom_line(aes(x = SampleSize, 
                           y = Mean, 
                           color = Pipeline, 
                           group = Pipeline),
                       position = position_dodge(width = 0.25)) + 
    ggplot2::geom_pointrange(aes(x = SampleSize, 
                                 ymin = Lower, 
                                 ymax = Upper, 
                                 color = Pipeline),
                             size = 0.5, 
                             position = position_dodge(width = 0.25)) +
    ggplot2::scale_color_manual(values = pipeline.cols) +
    ggplot2::theme_light() + 
    ggplot2::scale_y_continuous(limits=c(0.25, 1), 
                                breaks = c(0.25,0.5,0.75,1)) +
    ggplot2::labs(x = NULL,
                  y = perflabel, title = "Performance of representative pipelines") +
    ggplot2::facet_wrap(~DEPattern, nrow = 2, ncol=2,  labeller = as_labeller(relabels)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                   legend.position = "bottom",
                   legend.box = "vertical",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(0.5, "lines"),
                   axis.text.x=ggplot2::element_text(size=8, color='black', angle = 45, hjust = 1),
                   axis.text.y=ggplot2::element_text(size=8, color='black'),
                   axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white"))
  
  
  final.recommend <- res.recommend.plot
  # combine A and B
  upperpanel <- cowplot::plot_grid(final.pseudoR, 
                                   final.recommend,  
                                   ncol = 2, 
                                   align = "hv", axis = "lr",
                                   labels = c("A", "B"),
                                   label_size = 10,
                                   rel_widths = c(1, 0.75))
  
  
  pbmc.pipeline.plot <- readRDS(file = "pipeline/Pipelines-BluePrint_MainTypes.rds")
  
  # panel (F): schematic recommendation 
  recommend.plot <- ggdraw() + 
    cowplot::draw_image(magick::image_read("pipeline.001.tiff", density = 1200))
  
  final.plot <- cowplot::plot_grid(upperpanel,
                                   pbmc.pipeline.plot,
                                   recommend.plot,
                                   labels = c("", "", "F"),
                                   label_size = 10,
                                   nrow = 3,
                                   ncol = 1,
                                   rel_heights = c(0.5,0.5,0.5))
  
  ggplot2::ggsave(filename = paste0("Pipeline_Comparison_",
                                    modelfit, "_",
                                    measure, "_",
                                    gsub(" ", "", samplesize, fixed = TRUE),
                                    ".pdf"),
                  plot = final.plot,
                  width=210,
                  height=300,
                  units="mm")
}

plot.pseudoR.recommend(measure = measure, 
                       samplesize = samplesize, 
                       tenx = tenx, 
                       average = average, 
                       modelfit = modelfit,
                       favorite = favorite, 
                       naive = naive, 
                       perfdat = perfdat, 
                       perflabel = perflabel)


# SUPPLEMENT --------------------------------------------------------------

modelfit = "Norm+Norm:Spike"
average = TRUE
tenx = TRUE
measure = "MCC_lib"
samplesize = "384 vs 384"

aucdattype <- ifelse(isTRUE(average), "-avgauc.rds$", "-rawauc.rds$")
reflabel <- ifelse(isTRUE(tenx), "-10XGenomics", "-Smartseq2")
filext <- paste0(reflabel, aucdattype)

ModelAnnot <- generateLabels(type = modelfit)

# get the appropiate performance measure: MCC_lib
mcc <- readRDS(file = "MCC_lib.rds")
# setup annotation info
evaljobs.dat <- readRDS(file = "evalfulldat.rds")
umi.info <- data.frame(Protocol = c("10XGenomics", "CELseq2", "Dropseq", "SCRBseq", "Smartseq2"),
                       Type = c("UMI", "UMI", "UMI", "UMI", "Full-Length"), 
                       stringsAsFactors = F)
# add UMI vs nonUMI method column
evaljobs.dat <- evaljobs.dat %>%
  dplyr::left_join(umi.info, by="Protocol") 

# performance data must be in the format Name, SampleSize, Value
perfdat <- mcc %>% 
  dplyr::rename(Value = MCC)
perfdat$SampleSize <- factor(perfdat$SampleSize, 
                             levels = c("96 vs 96", "50 vs 200", "384 vs 384"))

favorite <- c("SCRBseq", "star", "gencode", 
              "none", "scran", "scran w/ groups", "limma-trend")

naive <- c("SCRBseq", "bwa", "gencode", 
           "none", "MR", "T-Test")

perflabel <- "MCC"

plot.pseudoR.recommend.SI <- function(measure, 
                                      samplesize, 
                                      tenx, 
                                      average, 
                                      modelfit,
                                      favorite, 
                                      naive, 
                                      perfdat, 
                                      perflabel) {
  
  # panel of performance measure of naive and favorite setup across sample sizes and DE settings
  jobs2eval.favorite <- evaljobs.dat %>% dplyr::filter(Protocol == favorite[1] &
                                                         Mapper == favorite[2] &
                                                         Annotation ==  favorite[3] &
                                                         Preprocessing == favorite[4] &
                                                         Normalisation  %in% favorite[5:6] &
                                                         DEtool == favorite[7] &
                                                         DEFilter == "raw") %>% 
    dplyr::mutate(Pipeline ="Good") %>% 
    dplyr::mutate(Normalisation = "scran")
  jobs2eval.naive <- evaljobs.dat %>% dplyr::filter(Protocol == naive[1] &
                                                      Mapper == naive[2] &
                                                      Annotation ==  naive[3] &
                                                      Preprocessing == naive[4] &
                                                      Normalisation  == naive[5] &
                                                      DEtool == naive[6] &
                                                      DEFilter == "raw") %>% 
    dplyr::mutate(Pipeline ="Naive")
  jobs2eval <- dplyr::bind_rows(jobs2eval.favorite, jobs2eval.naive) %>% 
    droplevels()
  
  pdat <- jobs2eval %>% 
    dplyr::left_join(perfdat, by = "Name") %>% 
    tidyr::replace_na(list(Value = 0)) %>% 
    dplyr::filter(! perc.DE == "0.00")
  
  
  sts.dat <- pdat %>%
    dplyr::group_by(perc.DE, LFC, SampleSize, Pipeline,
                    Protocol, Mapper, Annotation, Preprocessing, 
                    Normalisation, Spike, DEtool) %>%
    dplyr::summarise_at(dplyr::vars(Value), 
                        dplyr::funs(Mean = Hmisc::smean.cl.boot(.)[1],
                                    Lower = Hmisc::smean.cl.boot(.)[2],
                                    Upper = Hmisc::smean.cl.boot(.)[3])) %>% 
    ungroup()
  
  pipeline.cols <- c("Naive" = "purple", 
                     "Good" = "orange")
  
  plot.dat <- sts.dat %>% 
    tidyr::unite("DEPattern", c("perc.DE", "LFC"), sep="_", remove = F)  %>% 
    droplevels()
  
  res.recommend.plot <- ggplot2::ggplot(data = plot.dat, 
                                        aes(x = SampleSize, 
                                            y = Mean, 
                                            color = Pipeline,
                                            group = interaction(Pipeline, Spike),
                                            alpha = Spike)) +
    ggplot2::geom_line(aes(x = SampleSize, 
                           y = Mean, 
                           color = Pipeline, 
                           group = interaction(Pipeline, Spike),
                           alpha = Spike),
                       position = position_dodge(width = 0.25)) + 
    ggplot2::geom_pointrange(aes(x = SampleSize, 
                                 ymin = Lower, 
                                 ymax = Upper, 
                                 group = interaction(Pipeline, Spike),
                                 color = Pipeline,
                                 alpha = Spike),
                             size = 0.5, 
                             position = position_dodge(width = 0.25)) +
    ggplot2::scale_color_manual(values = pipeline.cols) +
    ggplot2::scale_alpha_manual(values=c(1,0.25)) +
    ggplot2::theme_light() + 
    ggplot2::scale_y_continuous(limits=c(0, 1), 
                                breaks = c(0,0.25,0.5,0.75,1)) +
    ggplot2::labs(x = NULL,
                  y = perflabel, title = "Performance of representative pipelines") +
    ggplot2::facet_grid(perc.DE~LFC,  labeller = as_labeller(relabels)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size=8, color='black'),
                   legend.position = "bottom",
                   legend.box = "horizontal",
                   legend.title = ggplot2::element_blank(),
                   legend.key.size = grid::unit(0.5, "lines"),
                   axis.text.x=ggplot2::element_text(size=8, color='black', angle = 45, hjust = 1),
                   axis.text.y=ggplot2::element_text(size=8, color='black'),
                   axis.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   plot.title=ggplot2::element_text(size=10, face="bold", color='black'),
                   strip.text = ggplot2::element_text(size=8, face="bold", color='black'),
                   strip.background = ggplot2::element_rect(fill="white"))
  
  
  ggplot2::ggsave(filename = paste0("Pipeline_Comparison_SI_MCC_",
                                    modelfit, "_",
                                    measure, "_",
                                    gsub(" ", "", samplesize, fixed = TRUE),
                                    ".pdf"),
                  plot = res.recommend.plot,
                  width=210,
                  height=210,
                  units="mm")
}

plot.pseudoR.recommend.SI(measure = measure, 
                          samplesize = samplesize, 
                          tenx = tenx, 
                          average = average, 
                          modelfit = modelfit,
                          favorite = favorite, 
                          naive = naive, 
                          perfdat = perfdat, 
                          perflabel = perflabel)

