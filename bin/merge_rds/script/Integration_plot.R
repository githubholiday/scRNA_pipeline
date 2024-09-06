library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'outdir',	'o',	1,	"character"
	),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)

prefix <- opt$prefix
outdir <- opt$outdir

library(Seurat)
library(ggplot2)
library(cowplot)
library(SCP)

# data为seurat对象
setwd(paste(outdir,'1_Com_QC',sep='/'))
# 合并样本后
data <- readRDS("Combine_qc_before.rds")
qc_before <- FeatureStatPlot(data, stat.by = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "stim",palette = "Set2") & NoLegend()&theme(axis.title.y = element_blank())
pdf("Pre_Gene_UMI_mito_percent.pdf",h=5,w=8)
print(qc_before)
dev.off()
# 过滤后
data <- readRDS("Combine.rds")
qc_after <- FeatureStatPlot(data, stat.by = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "stim",palette = "Set2") & NoLegend() &theme(axis.title.y = element_blank())
pdf("Gene_UMI_mito_percent.pdf",h=5,w=8)
print(qc_after)
dev.off()

# UMAP图：cluster和样本分布
setwd(paste(outdir,'2_Com_clusters',sep='/'))
data <- readRDS("Combine_cluster.rds")
pdf(paste(prefix,"umap_cluster_samples.pdf",sep='_'),w=24,h=8)
p2 <- CellDimPlot(srt = data, group.by = c("seurat_clusters"), label = TRUE,label.size = 4,reduction = "UMAP",label_insitu=TRUE,palette = "Set2")
p1 <- CellDimPlot(srt = data, group.by = "stim", label = FALSE,reduction = "UMAP",label_insitu=TRUE,palette = "Set2")
p3<-plot_grid(p1, p2)
print(p3)
p4 <- CellDimPlot(srt = data, group.by = c("seurat_clusters"), split.by="stim",label = TRUE,label.size = 4,reduction = "UMAP",label_insitu=TRUE,palette = "Set2")
print(p4)
dev.off()

pdf(paste(prefix,"tsne_cluster_samples.pdf",sep='_'),w=24,h=8)
p1 <- CellDimPlot(srt = data, group.by = "stim", label = FALSE,reduction = "tSNE",label_insitu=TRUE,palette = "Set2")
p2 <- CellDimPlot(srt = data, group.by = c("seurat_clusters"), label = TRUE,label.size = 4,reduction = "tSNE",label_insitu=TRUE,palette = "Set2")
p3<-plot_grid(p1, p2)
print(p3)
p4 <- CellDimPlot(srt = data, group.by = c("seurat_clusters"), split.by="stim",label = TRUE,label.size = 4,reduction = "tSNE",label_insitu=TRUE,palette = "Set2")
print(p4)
dev.off()

# 比例分布柱状图
pdf(paste(prefix,"sample.clusters_stats.pdf",sep='_'),w=6,h=12)
p <- CellStatPlot(data,stat.by = "seurat_clusters",group.by="stim",plot_type = "bar",palette = "Set2")
print(p)
dev.off()

