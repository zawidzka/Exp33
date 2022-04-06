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

load(file = "workspaceDA2.rds")

CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")

plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", 
                by = "cluster_id", scale = "first", bars = TRUE, perc = TRUE)

# set outputPath
plotsPath <- paste(getwd(), "plots2", sep = "/")
dir.create(plotsPath)

display.brewer.all(colorblindFriendly = FALSE)

# BM vs TUM
plotDiffHeatmap(sce, da1, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test1 <- as_tibble(da1)
p.adj.signif_1 <- c("ns", "*", "*", "ns", "ns", "ns")
group1 <- (rep("BM", nrow(stat_test1)))
group2 <- (rep("TUM", nrow(stat_test1)))
y.position <- c(90, 60, 90, 50, 60, 30)
stat.test1 <- cbind(stat_test1, group1, group2, p.adj.signif_1, y.position)
stat.test1 <- as_tibble(stat.test1)
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp1 <- bxp + stat_pvalue_manual(stat.test1, label = "p.adj.signif_1", tip.length = 0.01, size = 2.5)
bxp1

# BM v LN

plotDiffHeatmap(sce, da2, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test2 <- as_tibble(da2)
p.adj.signif_2 <- c("*", "*", "ns", "ns", "*", "*")
group1 <- (rep("BM", nrow(stat_test2)))
group2 <- (rep("LN", nrow(stat_test2)))
y.position <- c(90, 60, 90, 50, 60, 30)
stat.test2 <- cbind(stat_test2, group1, group2, p.adj.signif_2, y.position)
stat.test2 <- as_tibble(stat.test2)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp2 <- bxp1 + stat_pvalue_manual(stat.test2, label = "p.adj.signif_2", tip.length = 0.01, size = 2.5, step.increase = 0.05)
bxp2
# 
# # LN v TUM
plotDiffHeatmap(sce, da3, top_n = 6, all = TRUE, fdr = FDR_cutoff)
stat_test3 <- as_tibble(da3)
p.adj.signif_3 <- c("ns", "*", "*", "ns", "*", "*")
group1 <- (rep("LN", nrow(stat_test3)))
group2 <- (rep("TUM", nrow(stat_test3)))
y.position <- c(90, 60, 90, 50, 60, 30)
stat.test3 <- cbind(stat_test3, group1, group2, p.adj.signif_3, y.position)
stat.test3 <- as_tibble(stat.test3)
# bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
bxp3 <- bxp2 + stat_pvalue_manual(stat.test3, label = "p.adj.signif_3", tip.length = 0.01,
                                    size = 2.5, step.increase = 0.1)
bxp3

# view by sample to see outliers
# bxp2 <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "sample_id")
# bxp2 <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2.5)
# bxp2

ggsave("abundances_stat.svg", plot = last_plot(), dpi = 300)
ggsave("abundances_stat.pdf", plot = last_plot(), dpi = 300)



svg("ExprHeatmap.svg", width = 10, height = 10)
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "median",
                scale = "last", bars = FALSE, perc = FALSE)
dev.off()

pdf("ExprHeatmap.pdf", width = 10, height = 10)
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "median",
                scale = "first", bars = TRUE, perc = TRUE)
dev.off()


plotDR(sce, dr = "UMAP", color_by = "condition") 
# + geom_density2d(binwidth = 0.016, colour = "black")  # + geom_point(aes(size = 0.5))
ggsave("UMAP_wContours.svg", plot = last_plot())
ggsave("UMAP.pdf", plot = last_plot())
ggsave("UMAP2.pdf", height = 4, width = 4, plot = last_plot())


CATALYST::plotDR(sce, dr = "UMAP", color_by = c("TCF1", "PD1", "TIM3", "EOMES", "Tbet","CD44", "CD62l", "CD127", "Ly108", "CD45.1", "CD45.2"), facet_by = "condition")

plotDR(sce, dr = "UMAP", color_by = c("CTV", "CD107a", "IL2", "IFNg", "TNFa","OVATet", "Ly108", "TIM3", "PD1", "CD101", "CX3CR1"), 
                              facet_by = "condition") # +  geom_density2d(binwidth = 0.016, colour = "grey")

ggsave("UMAP_wContours_bymarker.svg", height = 4, width = 20, plot = last_plot())
ggsave("UMAP_bymarker.pdf", height = 4, width = 20, plot = last_plot())



CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "median")
ggsave("pbMDS.svg", plot = last_plot())
ggsave("pbMDS.pdf", plot = last_plot())


svg("Multiheatmap.svg", width = 12, height = 10)
plotMultiHeatmap(sce, k = "cluster_annotation", hm1 = type_markers(sce), 
                 hm2 = "abundances", bars = FALSE, perc = FALSE, 
                 row_anno = FALSE, scale = "last")
dev.off()

pdf("Multiheatmap.pdf", width = 12, height = 10)
plotMultiHeatmap(sce, k = "cluster_annotation", hm1 = type_markers(sce), 
                 hm2 = "abundances", bars = FALSE, perc = FALSE, 
                 row_anno = FALSE, scale = "first")
dev.off()

pdf("Multiheatmap2.pdf", width = 8, height = 10)
plotMultiHeatmap(sce, k = "cluster_annotation", hm1 = type_markers(sce), 
                 hm2 = "state", bars = FALSE, perc = FALSE, 
                 row_anno = FALSE, scale = "first")
dev.off()




plotClusterHeatmap(sce,
                   hm2 = NULL, k = "meta6", m = NULL, cluster_anno = TRUE, draw_freqs = TRUE)
plotClusterExprs(sce, k = "meta6", features = "state")
ggsave("clusterexp_CTV.pdf", plot = last_plot())

plotClusterExprs(sce, k = "meta8", features = "type")
ggsave("clusterexp_type.pdf", plot = last_plot())
