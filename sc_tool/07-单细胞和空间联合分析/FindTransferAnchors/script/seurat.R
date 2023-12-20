#/annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/R
#整合映射分析
library('getopt')
para<- matrix(c(
	'prefix',	'p',	1,	"character",
	'RNArds',	'r',	1,	"character",
	'outdir',	'o',	1,	"character",
	'configini',	'c',	1,	"character",
	'spacerds',	's',	1,	"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	========================================================================================================================================
	Usage example:
	Rscript this.r -p prefix -r RNArds -o outdir -c config.ini -s spacerds
	Options:
	--help  h	NULL		get this help
	--prefix	p	character	the prefix for outputfiles [forced]
	--RNArds	r	character	RNArdsfile[forced]
	--outdir	o	character	The	resurt of out dir for analysis [forced]
	--configini	c	character	config.ini  [forced]
	--spacerds	s	character	spacerdsfile  [forced]
	\n")
	q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$RNArds) )	{ cat("Please input the RNArds data file ...\n\n") ; print_usage(para)}
if ( is.null(opt$prefix) )	{ cat("Please input the prefix for analysis ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para)}
if ( is.null(opt$configini) )	{ cat("Please give the configini file for analysis ...\n\n") ; print_usage(para)}
if ( is.null(opt$spacerds) )	{ cat("Please input the spacerds data file ...\n\n") ; print_usage(para)}

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(configr)
library(scibet)
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis)) 
suppressMessages(library(ggsci))
library(readr)
###创建目录
mkdirs <- function(outdir,fp) {
	if(!file.exists(file.path(outdir,fp))) {
		dir.create(file.path(outdir,fp))
	}else{
		print(paste(fp,"Dir already exists!",sep="     "))
		}
}

# 读取config.ini 文件，获取sample list 和Assay
prefix<-opt$prefix
outdir<- opt$outdir
mkdirs(outdir,prefix)
setwd(paste(outdir,prefix,sep='/'))
ini<-opt$configini
ini.list <- read.config(file = ini)


get_confusion_plot <- function(actual, predicted, align = 1){
  require(reshape2)
  require(ggplot2)

  if (class(actual) == "factor"){ actual_levels <- levels(actual) }
  else { actual_levels <- unique(actual) }
  
  if (class(predicted) == "factor"){ predicted_levels <- levels(predicted) }
  else { predicted_levels <- unique(predicted) }
  
  get_confusion_prop_tbl <- function(actual, predicted, align = 1){
    if (!(align %in% c(1,2))) {stop("align must be 1 or 2")}
    conf <- table(actual, predicted)
    if (align == 2) {
      conf <- t(conf)
      return (t(conf / rowSums(conf)))
    }
    return (conf / rowSums(conf))
  }
  
  conf_prop <- get_confusion_prop_tbl(actual, predicted, align)
  conf_prop_melt <- reshape2::melt(conf_prop, variable.name= predicted, id = actual)

  conf_prop_melt[,"actual"] = factor(conf_prop_melt$actual, levels=actual_levels)
  conf_prop_melt[,"predicted"] = factor(conf_prop_melt$predicted, levels=predicted_levels)
  
  p0<-theme(legend.text = element_text(face="plain", color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle=90),
        axis.text.y = element_text(color="black",size=10),
        axis.title.x = element_text(face="plain", color="black",size=15),
        axis.title.y = element_text(face="plain", color="black",size=15))
  plt <- ggplot2::ggplot(data = conf_prop_melt, aes(x=predicted,y=actual)) +
    geom_tile(aes(fill=value)) +theme_bw()+p0+scale_fill_gradientn(colours=c('#46085B','#453480','#20918C','#5AC763','#DFE318'), guide=guide_colorbar(reverse=F))
  
  #print(plt)
  return(plt)
}
#space
tmp<-readRDS(opt$spacerds)
#RNA
tmp1<-readRDS(opt$RNArds)

anchors_ims<-as.numeric(ini.list$Para$anchors_ims)
celltype<-ini.list$Para$RNAcelltype
seurat_clusters<-ini.list$Para$spacecluster
common_gene_num <-as.numeric(ini.list$Para$common_gene)
#检查celltype是否存在rds中
if(is.null(tmp1@meta.data[,celltype])){
	print(paste("celltype变量在",opt$RNArds,"文件中不存在，流程退出",sep="_"))
	q();
}
#检查seurat_clusters是否存在spacerds中
if(is.null(tmp@meta.data[,seurat_clusters])){
	print(paste("seurat_clusters变量在",opt$spacerds,"文件中不存在，流程退出",sep="_"))
	q();
}
#检查两者的基因名字是否有高于1000的重复
gene_space<-row.names(tmp@assays$Spatial@counts)
gene_RNA<-row.names(tmp1@assays$RNA@counts)
common_gene<-intersect(gene_space,gene_RNA)
if(length(common_gene)>common_gene_num){
	print("RNA和空间数据集中的gene重复超过1000个，继续分析")
}else{
	print(paste("RNA和空间数据集中的gene重复数为",length(common_gene),"流程退出",sep=""))
	q();
}

#输出每种细胞类型的细胞数量
cell_num<-as.data.frame(table(tmp1@meta.data[,celltype]))
names(cell_num)<-c("celltype","count")
cell_num<-cell_num[order(cell_num$count),]
write.table(cell_num,paste(prefix,"celltype_count.xls",sep="_"),sep="\t",quote=F,row.names=F)

cells <- unique(tmp1@meta.data[,celltype])
Idents(tmp1)<- tmp1@meta.data[,celltype]

#输出细胞类型的umap图
if(is.null(tmp1@reductions$umap)){
	print("RNA数据集中不存在umap slot，需要运行以下步骤")
	sce_scale <- ScaleData(tmp1, verbose = FALSE)
	sce_scale <- RunPCA(sce_scale, npcs = as.numeric(ini.list$Para$integration_pca_dims), verbose = FALSE)
	sce_umap <- RunUMAP(sce_scale, reduction = "pca",dims = 1:as.numeric(ini.list$Para$RunUMAP_dim))
}else{
	print("RNA数据集中存在umap slot,不再重新分析")
	sce_umap<-tmp1
}
DimPlot(sce_umap, reduction = "umap", group.by = celltype,pt.size = 1)
ggsave(paste(prefix,"celltype_umap.pdf",sep="_"),width = 12, height = 10)
#检查RNA数据集中是否有SCT的数据槽，如果没有就添加
if(is.null(tmp1@assays$SCT)){
	print("RNA数据集中不存在SCT，需要运行以下步骤")
	SCT_tmp1 = SCTransform(tmp1 , assay = "RNA" , return.only.var.genes = FALSE, verbose = FALSE)
	tmp1 = SCT_tmp1
}else{
	print("RNA数据集中存在SCT，可以直接进行后续分析")
}
#输出空间数据的umap图

graph <- paste(prefix,"cluster_umap.pdf",sep="_")
pdf(graph,w=12,h=8)
p1 <- DimPlot(tmp, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(tmp, label = TRUE, label.size = 4)
plot_grid(p1, p2,align = "h")
dev.off()

anchors <- FindTransferAnchors(reference = tmp1, query = tmp, normalization.method = ini.list$Para$normalization.method,dims = 1:anchors_ims)        
predictions <- TransferData(anchorset = anchors, refdata = tmp1@meta.data[,celltype], dims = 1:anchors_ims)
      
pancreas.query0 <- AddMetaData(tmp, metadata = predictions)

saveRDS(pancreas.query0, file = paste(prefix,"spot_cell.rds",sep='_'))
spot<-row.names(predictions)
predictions_all<-cbind(spot, predictions)
write.table(predictions_all,paste(prefix,"predictions.xls",sep="_"),sep="\t",quote=F,row.names=F)

write.table(t(predictions_all)[,1],paste(prefix,"predictions_example.xls",sep="_"),sep="\t",quote=F,row.names=T,col.names=F)

pancreas.query<-pancreas.query0
#统计spot每种细胞的数量
cell_num_space<-as.data.frame(table(pancreas.query$predicted.id))
names(cell_num_space)<-c("celltype","count")
cell_num_space<-cell_num_space[order(cell_num_space$count),]
cell_num_space1<-cell_num_space[cell_num_space$count>0,]
write.table(cell_num_space1,paste(prefix,"Spot_cell_count.xls",sep="_"),sep="\t",quote=F,row.names=F)

#将所有细胞展示在切片上
Idents(pancreas.query) <- pancreas.query$predicted.id
SpatialDimPlot(pancreas.query) + theme(legend.position = "right")
ggsave(file=paste(prefix,"SpatialDimPlot.pdf",sep="_"),width=10, height=8)

#将每种细胞单独展示
SpatialDimPlot(pancreas.query,cells.highlight = CellsByIdentities(object = pancreas.query, idents = unique(Idents(pancreas.query))),facet.highlight = TRUE)
ggsave(file=paste(prefix,"SpatialDimPlot_Pit1_eachcell.pdf",sep="_"),width=10, height=8)

#统计每个cluster中不同细胞的数量
cell_num_space<-as.data.frame(table(pancreas.query$predicted.id,pancreas.query@meta.data[,seurat_clusters]))
names(cell_num_space)<-c("celltype","cluster","count")
cell_num_space<-cell_num_space[order(cell_num_space$celltype),]
cell_num_space1<-cell_num_space[cell_num_space$count>0,]
write.table(cell_num_space1,paste(prefix,"celltype_count_space_cluster.xls",sep="_"),sep="\t",quote=F,row.names=F)
library(scibet)
p<-get_confusion_plot(pancreas.query@meta.data[,seurat_clusters], pancreas.query$predicted.id)+geom_text(aes(label = round(value,2)))+coord_flip()+ggtitle('')
pdf(paste(prefix,'heatmap_cluster.pdf',sep="_"),w=10,h=4)
print(p)
dev.off()

#write.table(pancreas.query0@meta.data,file=paste(prefix,'heatmap_spot_metadata.xls',sep="_"), quote=F,sep='\t')

my36colors <-c('#53A85F','#D6E7A3','#57C3F3','#476D87','#E95C59','#E59CC4','#AB3282','#23452F','#BD956A','#8C549C','#585658','#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3', '#E4C755','#F7F398', '#AA9A59','#E63863','#E39A35','#C1E6F3', '#6778AE','#91D0BE','#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6','#625D9E','#68A180','#3A6963','#968175')
prefix<-opt$prefix
pdf(paste0(prefix,'_predicted.id_umap1.pdf'),w=7,h=5)
p2 <- DimPlot(pancreas.query, reduction = "umap", label = F,group.by = "predicted.id",label.size = 4,pt.size = 0,raster=FALSE,cols=my36colors)+ggtitle("")
print(p2)
dev.off()

#https://www.jianshu.com/p/7b75b109cfb5
#https://www.jianshu.com/p/5335f7246c9a
#https://www.jianshu.com/p/3a535af7cf8b
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
#https://satijalab.org/seurat/articles/spatial_vignette.html
