rm(list = ls())

# Load packages
library(rstudioapi)
library(devtools)
library("flowCore")
library("flowWorkspace")
library(cytofCore)
library(FlowSOM)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(flowViz)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(uwot)
library(CATALYST)
library(diffcyt)
library(SummarizedExperiment)
library(stringr)
library(ggcyto)
library(SingleCellExperiment)
library(scran)
library(scater)
library(readxl)
library(flowStats)
library(FlowSOMworkshop)
library(tidyverse)
library(data.table)
library(ggpubr)
library(flowAI)
library(PeacoQC)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# set workingDir
workingDir <- "211021_WorkingDirectory"
workingDirPath <- paste(PrimaryDirectory, workingDir, sep = "/")
setwd(workingDirPath)

sce <- readRDS("SCE_BMvTvLN_concat_DR.rds")
CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta6", 
                by = "cluster_id", scale = "first", bars = TRUE, perc = TRUE)


# annotate clusters

annotation_table <- as.data.frame(cbind(c(1:6), c(1:6)))

colnames(annotation_table) <- c("meta6", "Clusters")
annotation_table$Clusters <- factor(annotation_table$Clusters)
sce <- mergeClusters(sce, k = "meta6", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)
sce$cluster_annotation <- cluster_ids(sce, "cluster_annotation")
# filtered_sce <- filterSCE(sce, cluster_id %in% c(paste0("C", c(1:12))), k = "cluster_annotation")

# store original_sce
# old_sce <- sce
# sce <- filtered_sce

FDR_cutoff <- 0.05
ei <- sce@metadata$experiment_info
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")

# DA using edgeR
design <- createDesignMatrix(ei,
                             cols_design = c("condition"))
contrastBMvTUM <- createContrast(c(0, 1, 0))
contrastBMvLN <-createContrast(c(0, 0, 1))

nrow(contrastBMvLN) == ncol(design)
nrow(contrastBMvTUM) == ncol(design)

# contrast <- createContrast(c(0,1))

# nrow(contrast) == ncol(design)

out_DA1 <- diffcyt(sce,
                  experiment_info = ei, design = design, contrast = contrastBMvTUM,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                  transform = FALSE, normalize = FALSE, min_samples = 1)

da1 <- rowData(out_DA1$res)
plotDiffHeatmap(sce, da1, top_n = 6, all = TRUE, fdr = FDR_cutoff)

out_DA2 <- diffcyt(sce,
                   experiment_info = ei, design = design, contrast = contrastBMvLN,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da2 <- rowData(out_DA2$res)
plotDiffHeatmap(sce, da2, top_n = 6, all = TRUE, fdr = FDR_cutoff)

# create a different design matrix for contrast of other groups, reference level switch
tissue <- factor(c("BM", "BM", "BM","LN", "LN", "LN", "TUM", "TUM","TUM"))
model.matrix(~tissue)
tissue <- relevel(tissue, ref = "LN")
design2<-model.matrix(~tissue)

# create another contrast
contrastLNvTUM <- createContrast(c(0, 0, 1))


out_DA3 <- diffcyt(sce,
                   experiment_info = ei, design = design2, contrast = contrastLNvTUM,
                   analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE, min_samples = 1)

da3 <- rowData(out_DA3$res)
plotDiffHeatmap(sce, da3, top_n = 6, all = TRUE, fdr = FDR_cutoff)

# design2 <-createDesignMatrix(ei, cols_design = c("condition", "patient_id"))
# design2



save(list = ls(), file = "workspaceDA2.rds")










