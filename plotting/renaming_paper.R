pattern_names <- c(`Symmetric` = "Symmetric", 
                   `Asymmetric` = "Asymmetric", 
                   `Completely Asymmetric` = "Completely \nAsymmetric")
de_names <- c(`0.05` = "5% DE", `0.2` = "20% DE", `0.6` = "60% DE")
sim_names <- c("0.05_lfc.gamma.narrow.asym" = "5% Asymmetric", 
               "0.05_lfc.gamma.narrow.casym" = "5% Completely \nAsymmetric", 
               "0.05_lfc.gamma.narrow.sym" = "5% Symmetric", 
               "0.2_lfc.gamma.narrow.asym" = "20% Asymmetric", 
               "0.2_lfc.gamma.narrow.casym"  = "20% Completely \nAsymmetric", 
               "0.2_lfc.gamma.narrow.sym" = "20% Symmetric", 
               "0.6_lfc.gamma.narrow.asym" = "60% Asymmetric", 
               "0.6_lfc.gamma.narrow.casym" = "60% Completely \nAsymmetric", 
               "0.6_lfc.gamma.narrow.sym" = "60% Symmetric")

de_pattern_names <- c("0.05_Asymmetric" = "5% Asymmetric DE", 
               "0.05_Completely Asymmetric" = "5% Completely \nAsymmetric DE", 
               "0.05_Symmetric" = "5% Symmetric DE", 
               "0.2_Asymmetric" = "20% Asymmetric DE", 
               "0.2_Completely Asymmetric"  = "20% Completely \nAsymmetric DE", 
               "0.2_Symmetric" = "20% Symmetric DE", 
               "0.6_Asymmetric" = "60% Asymmetric DE", 
               "0.6_Completely Asymmetric" = "60% Completely \nAsymmetric DE", 
               "0.6_Symmetric" = "60% Symmetric DE")

samplesize_names <- c(`50 vs 200` = "50 vs 200",
                      `96 vs 96` = "96 vs 96",
                      `384 vs 384` = "384 vs 384")
protocol_names <- c(`CELseq2` = "CEL-seq2", 
                    `SCRBseq` = "SCRB-seq", 
                    `Smartseq2` = "Smart-seq2", 
                    `HGMM1kv3` = "10X - HGMM",
                    `PBMC1kv3` = "10X - PBMC",
                    `10XGenomics` = "10X - HGMM", 
                    `Dropseq` = "Drop-seq")
type_names <- c(`UMI` = "UMI", `Full-Length` = "Full-Length")
mapper_names <- c(`bwa` = "BWA",
                  `kallisto` = "kallisto",
                  `bowtie2` = "Bowtie 2",
                  `star` = "STAR")
annotation_names <- c(`gencode` = "GENCODE",
                      `vega` = "VEGA",
                      `refseq` = "RefSeq")
preprocess_names <- c(`none` = "none", 
                      `Filtering` = "Filtering", 
                      `DrImpute` = "DrImpute", 
                      `scImpute` = "scImpute", 
                      `scone` = "scone", 
                      `Seurat` = "Seurat",
                      `SAVER` = "SAVER")
norm_names <- c(`SF`="Simulated Size Factors", 
                `MR`="MR", 
                `PosCounts`="PosCounts", 
                `TMM`="TMM", 
                `Census`="Census", 
                `Linnorm`="Linnorm", 
                `SCnorm`="SCnorm with groups", 
                `SCnorm w/ cluster`="SCnorm with cluster", 
                `scran`="scran", 
                `scran w/ groups`="scran with groups", 
                `scran w/ cluster`="scran with cluster")
spike_names <- c(`w/o Spike-Ins`="w/o Spike-Ins", 
                 `with Spike-Ins` = "with Spike-Ins")
detool_names <- c(`limma-trend` = "limma-trend", 
                  `T-Test` = "T-Test", 
                  `edgeR-zingeR` = "edgeR-zingeR", 
                  `MAST` = "MAST")

relabels <- c(pattern_names, 
              de_names, 
              de_pattern_names,
              sim_names,
              samplesize_names,
              protocol_names, 
              type_names, 
              mapper_names,
              annotation_names,
              preprocess_names,
              norm_names, 
              spike_names,
              detool_names)