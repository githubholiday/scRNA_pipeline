#export LD_LIBRARY_PATH=/opt/glibc-2.14/lib:D_LIBRARY_PATH:IBRARY_PATH && /annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/R
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

#输入space的rds
heart_ST = readRDS(opt$spacerds)
#输入scRNA的rds
sce = readRDS(opt$RNArds)

celltype<-ini.list$Para$RNAcelltype
seurat_clusters<-ini.list$Para$spacecluster
common_gene_num <-as.numeric(ini.list$Para$common_gene)

#检查celltype是否存在rds中
if(is.null(sce@meta.data[,celltype])){
	print(paste("celltype变量在",opt$RNArds,"文件中不存在，流程退出",sep="_"))
	q();
}
#检查seurat_clusters是否存在spacerds中
if(is.null(heart_ST@meta.data[,seurat_clusters])){
	print(paste("seurat_clusters变量在",opt$spacerds,"文件中不存在，流程退出",sep="_"))
	q();
}
#检查两者的基因名字是否有高于1000的重复
gene_space<-row.names(heart_ST@assays$Spatial@counts)
gene_RNA<-row.names(sce@assays$RNA@counts)
common_gene<-intersect(gene_space,gene_RNA)
if(length(common_gene)>common_gene_num){
	print("RNA和空间数据集中的gene重复超过1000个，继续分析")
}else{
	print(paste("RNA和空间数据集中的gene重复数为",length(common_gene),"流程退出",sep=""))
	q();
}

#输出每种细胞类型的细胞数量
cell_num<-as.data.frame(table(sce@meta.data[,celltype]))
names(cell_num)<-c("celltype","count")
cell_num<-cell_num[order(cell_num$count),]
write.table(cell_num,paste(prefix,"celltype_count.xls",sep="_"),sep="\t",quote=F,row.names=F)

cells <- unique(sce@meta.data[,celltype])
Idents(sce)<- sce@meta.data[,celltype]

#输出细胞类型的umap图
if(is.null(sce@reductions$umap)){
	print("RNA数据集中不存在umap slot，需要运行以下步骤")
	sce_scale <- ScaleData(sce, verbose = FALSE)
	sce_scale <- RunPCA(sce_scale, npcs = as.numeric(ini.list$Para$integration_pca_dims), verbose = FALSE)
	sce_umap <- RunUMAP(sce_scale, reduction = "pca",dims = 1:as.numeric(ini.list$Para$RunUMAP_dim))
}else{
	print("RNA数据集中存在umap slot,不再重新分析")
	sce_umap<-sce
}
DimPlot(sce_umap, reduction = "umap", group.by = celltype,pt.size = 1)
ggsave(paste(prefix,"celltype_umap.pdf",sep="_"),width = 12, height = 10)

#输出空间数据的umap图

graph <- paste(prefix,"cluster_umap.pdf",sep="_")
pdf(graph,w=12,h=8)
p1 <- DimPlot(heart_ST, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(heart_ST, label = TRUE, label.size = 4)
plot_grid(p1, p2,align = "h")
dev.off()


markers = FindAllMarkers(sce, only.pos = TRUE, min.pct = as.numeric(ini.list$Para$marker_gene_min.pct), logfc.threshold = as.numeric(ini.list$Para$marker_gene_logfc.threshold))
#markers = read.table("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/AddModuleScore/ASCP_1/ASCP_1_all_markers.xls",sep="\t",header=TRUE)
markers <- markers %>% dplyr::select(gene,everything()) %>% subset(p_val_adj < as.numeric(ini.list$Para$marker_gene_padj))
#write.table(markers,paste(prefix,"all_markers.xls",sep="_"),sep="\t",quote=F,row.names=F)
cells <- unique(markers$cluster)
print(cells)
#print(as.character(cells)[0:3])
for (cell in as.character(cells)){
	print("##########################################################################################")
	print(cell)
	marker_10 = markers %>% group_by(cluster) %>% top_n(n = as.numeric(ini.list$Para$marker_gene_num), wt = avg_log2FC) %>% filter(cluster==cell)
	gene_list=list(gene=marker_10$gene)
	print(gene_list)
	print("-AddModuleScore")
	heart_ST <-AddModuleScore(heart_ST,features = gene_list,name = cell)
    #print(str(heart_ST))
	#if (is.null(heart_ST@meta.data[paste(cell,"1",sep="")]){
		#heart_ST@meta.data[cell]<- 0
	#}else{
		heart_ST@meta.data[cell]<-heart_ST@meta.data[paste(cell,"1",sep="")]
	#}
	str(heart_ST@meta.data[[paste(cell,"1",sep="")]]) #这里会将细胞类型后面+1
}
# 打分结果保存在heart_ST@meta.data[["NK1"]] #这里会将细胞类型后面+1
print("#################################################################################################")
print("运行到这里了1")
write.table(heart_ST@meta.data,paste(prefix,"celltype_AddModuleScore.xls",sep="_"),sep="\t",quote=F,row.names=T)
meta<-heart_ST@meta.data
spot<-row.names(meta)
meta_all<-cbind(spot, meta)
meta_all_csv<-t(meta_all)[,1]
needcol<-c("spot","orig.ident","seurat_clusters",as.character(cells))
print(needcol)
meta_all_csv<- as.data.frame(meta_all_csv)
print(meta_all_csv[,1])
meta_all_csv_se<-meta_all_csv %>% filter(rownames(.) %in% needcol)
write.table(meta_all_csv_se,paste(prefix,"celltype_AddModuleScore_example.xls",sep="_"),sep="\t",quote=F,row.names=T,col.names=F)
#输出score值的Vlnplot图
VlnPlot(heart_ST,features = as.character(cells),split.by = seurat_clusters, ncol = 3)
ggsave(paste(prefix,"celltype_Vlnplot.pdf",sep="_"),width = 12, height = 10)

#输出score值的featureplot图
FeaturePlot(object = heart_ST, features = as.character(cells),ncol = 3)
ggsave(paste(prefix,"celltype_FeaturePlot.pdf",sep="_"),width = 12, height = 10)

#通过数值获取每个spot唯一的细胞类型，传入meta.data数据里面
metadata <-heart_ST@meta.data
cell_score<-metadata[,as.character(cells)]
cell_score$spot<-row.names(cell_score)
cell_score <- melt(cell_score,id.vars = "spot")
region_term = cell_score %>% group_by(spot) %>% slice(which.max(value))
heart_ST@meta.data$orig.spot<-row.names(heart_ST@meta.data)
labers=region_term[match(as.character(heart_ST@meta.data$orig.spot),as.character(region_term$spot)),2]
labers<-c(labers)
heart_ST$term <-labers
#统计空间中的细胞类型数量
cell_num_space<-as.data.frame(table(heart_ST@meta.data$term))
names(cell_num_space)<-c("celltype","count")
cell_num_space<-cell_num_space[order(cell_num_space$count),]
cell_num_space1<-cell_num_space[cell_num_space$count>0,]
write.table(cell_num_space1,paste(prefix,"celltype_count_space.xls",sep="_"),sep="\t",quote=F,row.names=F)
#按照空间的cluster统计空间中的细胞类型数量
cell_num_space<-as.data.frame(table(heart_ST@meta.data$term,heart_ST@meta.data[,seurat_clusters]))
names(cell_num_space)<-c("celltype","cluster","count")
cell_num_space<-cell_num_space[order(cell_num_space$celltype),]
cell_num_space1<-cell_num_space[cell_num_space$count>0,]
write.table(cell_num_space1,paste(prefix,"celltype_count_space_cluster.xls",sep="_"),sep="\t",quote=F,row.names=F)
#绘制百分比堆积图
ggplot(data=cell_num_space1,aes(x=cluster,y=count,fill=celltype)) +geom_bar(stat="identity",position="fill")+
  labs(x="cluster",y="percentage")+
  theme_bw()+
  theme(axis.title.y=element_text(size=16), axis.title.x=element_text(size=16))+
  theme(legend.text=element_text(size=10))+
  theme(axis.text.x = element_text(size = 12, color = "black"))+
  theme(axis.text.y = element_text(size = 12, color = "black"))
ggsave(paste(prefix,"celltype_count_space_cluster.pdf",sep="_"),width = 12, height = 10)
#绘制空间数据中细胞注释后的umap图
DimPlot(heart_ST, reduction = "umap", label = F,group.by = "term",pt.size = 1.5,raster=FALSE)+ggtitle("")
ggsave(paste(prefix,"celltype_space_umap.pdf",sep="_"),width = 12, height = 10)

#将所有细胞展示在切片上
Idents(heart_ST) <- heart_ST$term
SpatialDimPlot(heart_ST) + theme(legend.position = "right")
ggsave(file=paste(prefix,"SpatialDimPlot.pdf",sep="_"),width=10, height=8)

#将每种细胞单独展示
SpatialDimPlot(heart_ST,cells.highlight = CellsByIdentities(object = heart_ST, idents = unique(Idents(heart_ST))),facet.highlight = TRUE)
ggsave(file=paste(prefix,"SpatialDimPlot_Pit1_eachcell.pdf",sep="_"),width=10, height=8)

saveRDS(heart_ST, file = paste(prefix,"spot_cell.rds",sep='_'))



# 获取空间位置信息
#position_ST_ninth = GetTissueCoordinates(heart_ST)
#metadata <-heart_ST@meta.data
#metadata$x=position_ST_ninth$imagecol
#metadata$y=position_ST_ninth$imagerow

#ggplot(metadata, aes(x = x, y = y)) + 
#geom_point(data = metadata, aes(color= Epithelial1))+
#theme_bw()+
#xlab(NULL)+
#ylab(NULL)+
#scale_color_continuous(high="red",low="white")
#ggsave("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/limeng/16_Space_RNA/seurat/AddModuleScore/Shell/Epithelial_score.pdf",width = 12, height = 10)

