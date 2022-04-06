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

# SCE_AllSamples is all samples
# SCE_AllSamples2 is all samples with new markers as state
# SCE_BMvT is BM and TUM only, with new markers as state

sce <- readRDS("SCE_BMvT_concat.rds")
CATALYST::pbMDS(sce, color_by = "sample_id", features = type_markers(sce), fun = "median")

# Run FlowSOM and ConsensusClusterPlus
seed <- 123456
set.seed(seed)
sce <- cluster(sce, features = "type", xdim = 10, ydim = 10, maxK = 20, 
               verbose = TRUE, seed = seed)
delta_area(sce)
# Run dimensionality reduction
n_cells <- 5000
n_events <- min(n_cells(sce))
ifelse(n_cells > n_events, n_cells <- n_events, n_cells <- n_cells)
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
# sce <- runDR(sce, dr = "TSNE", cells = n_cells, features = "type", theta = 0.5, max_iter = 1000, 
#              distMethod = "euclidean",
#              PCA = TRUE, eta = eta, exaggeration_factor = 12.0)
sce <- runDR(sce, dr =  "UMAP", cells = NULL, features = "type")
# sce <- runDR(sce, dr = "DiffusionMap", cells = n_cells, features = "type", assay = "exprs")
sce <- runDR(sce, dr = "PCA", cells = n_cells, features = "type", assay = "exprs")

saveRDS(sce, file = "SCE_BMvTvLN_concat_DR.rds")

display.brewer.all(colorblindFriendly = TRUE)
delta_area(sce)
cdx2 <- type_markers(sce)
cdx <- state_markers(sce)
plotMultiHeatmap(sce, k = "meta6",
                 hm1 = cdx2, hm2 = "abundances", 
                 bars = TRUE, perc = TRUE, row_anno = FALSE, scale = "first")

plotMultiHeatmap(sce, k = "meta8",
                 hm1 = cdx2, hm2 = "abundances", 
                 bars = TRUE, perc = TRUE, row_anno = FALSE, scale = "first")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta6", by = "cluster_id",  fun = "mean", scale = "first")


plotDR(sce, dr = "UMAP", color_by = "meta6", facet_by = "condition")
plotDR(sce, dr = "PCA", color_by = "condition")

# plot expression of type markers by cluster
plotClusterHeatmap(sce,
                   hm2 = NULL, k = "meta6", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE)
plotClusterExprs(sce, k = "meta6", features = "state")


CATALYST::plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "sample_id")
CATALYST::plotDR(sce, dr = "UMAP", color_by = c("TCF1", "PD1", "TIM3", "EOMES", "Tbet", "Ly108", "CD101", "CX3CR1"), facet_by = "condition")
plotDR(sce, dr = "UMAP", color_by = "TCF1", facet_by = "condition") +  
  geom_density2d(binwidth = 0.006, colour = "grey")
plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "condition") +  
  geom_density2d(binwidth = 0.016, colour = "grey") + geom_point(size = 1)
