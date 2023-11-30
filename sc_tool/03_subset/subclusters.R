#前言：本分析主要用于亚群细分，rpca,harmony,cca,pca四种方法可供选择。
##nfeatures<-2000;npc<-30;dims<-20 #默认一般参数。
#===========================================================
#!/annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript
#时间：20230728
#版本：v1.0.0
#用途：亚群细分， 包含以下目录:
#1_QC  2_cluster  3_marker 
#===========================================================
library('getopt')
para<- matrix(c(
        'help',         'h',    0,      "logical",
        'rds',          'r',    1,      "character",
        'config',       'c',    2,      "character",
        'prefix',       'p',    2,      "character",
        'outdir',       'o',    1,      "character"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
        cat(getopt(para,usage=TRUE))
        cat("
        ==================================================================================================================
		选取特定的细胞群进行亚群细分， 包含以下结果#1_QC  2_cluster  3_marker ...
		亚群细分的方法有：rpca,harmony,cca,pca
		推荐rpca或者harmony; cca整合分析主要用于大群，pca是指不进行整合分群，也就是不去批次效应亚群细分（有些整合效果不好也可以尝试，一般不推荐）
		#nfeatures<-2000;npc<-30;dims<-20 #默认一般参数。
        ==================================================================================================================
        Usage example:
        Rscript this.r -r rds -c config.ini -p prefix -o outdir
        Options:
        --help          h       NULL            get this help
        --rds           r       character       rds file for analysis by seurat [forced]
        --config        c       character       ini file for analysis by seurat [forced]
        --prefix        p       character       prefix of output [forced]
        --outdir        o       character       The resurt of out dir for analysis [forced]
        \n")
        q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )       { print_usage(para) }
if ( is.null(opt$rds) )         { cat("Please input the rds data ...\n\n") ; print_usage(para)}
if ( is.null(opt$config) )      { cat("Please give the congfig file for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$outdir) )      { cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )      { opt$prefix <- 'test'}


# 读取配置文件
library(configr)

ini.list <- read.config(file = opt$config)
rdsfile <-opt$rds
outdir<-opt$outdir
prefix<-opt$prefix

#参数
ident					 <-ini.list$Subclusters$ident
sub_clusters			 <-ini.list$Subclusters$sub_clusters
sub_samples			 	 <-ini.list$Subclusters$sub_samples
reclusters_method		 <-ini.list$Subclusters$reclusters_method
nfeatures     			 <-as.numeric(ini.list$Subclusters$findvariablefeatures)
npc     			 	 <-as.numeric(ini.list$Subclusters$reduction_npc_num)
dims     			 	 <-as.numeric(ini.list$Subclusters$reduction_dims_num)
res     				 <-as.numeric(ini.list$Subclusters$reduction_resolution)

####------------------------------------------------------------------------
suppressPackageStartupMessages({
#library(FlexDotPlot) #dotplot_cluster
library(Seurat)
library(dplyr)
library(tidyverse)
library(viridis)
library(ggalluvial)
library(ggsci)
library(gridExtra)
library(ggplot2)
library(ComplexHeatmap)
require(Matrix)
require(magrittr)
library(scales)
library(configr)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(harmony)
#library(SeuratData) #转h5
library(ggpubr) #ggboxplot
})



mkdirs <- function(outdir,fp) {
        if(!file.exists(file.path(outdir,fp))) {
#               mkdirs(dirname(fp))
                dir.create(file.path(outdir,fp))
        }else{
                        print(paste(fp,"Dir already exists!",sep="     "))
                        unlink(file.path(outdir,fp), recursive=TRUE)
                        dir.create(file.path(outdir,fp))
                }
}

join_c <- function(ts,Connector=''){
        paste(ts,collapse=Connector)
}

Parse_abspath_c <- function(input_abspath){
  tmp <- c( dirname(input_abspath), basename(input_abspath))
  return(tmp)
}

#读取rds文件
print("###........................................")
print("开始对rds进行subset......")
print(Sys.time())
print("")
rds0<-readRDS(rdsfile)

# rds0<-readRDS('/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/B_fenxi-007/blood_celltype/qudouble/celltype2/subset/CellIdent2/CellIdent/SLEvsHCB_immune_combined_celltype.rds')
# outdir<-'/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/TET_PUBLIC/yaojiaying/sctools/scRNA_0009/v1.0.0/scRNA_0009/bin/script/'
# prefix<-'test'


cluster_cols <-c('#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87','#E95C59','#E59CC4','#AB3282','#23452F','#BD956A','#8C549C','#585658','#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3', '#E4C755','#F7F398', '#AA9A59','#E63863','#E39A35','#C1E6F3', '#6778AE','#91D0BE','#B53E2B','#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6','#625D9E','#68A180','#3A6963','#968175')

#对clusters进行subset
Idents(rds0)<-ident
if (sub_clusters!='all'){
	subclusters<-unique(unlist(strsplit(sub_clusters,split = ",",fixed=T)))
	print(paste0("所有细胞群：",paste(unique(Idents(rds0)),collapse=",")  ))
	print(paste0("subset的细胞群：",paste(subclusters,collapse=",")  ))
	if (!all(subclusters %in% unique(Idents(rds0))  )){print("抱歉，您选取了不存在的群，请仔细检查config配置文件......");quit()}
	rds1<-subset(rds0,idents=subclusters)
	}else{
	rds1<-rds0
	}
	
	
#对样本进行subset
if (sub_samples!='all'){
	subsamples<-unique(unlist(strsplit(sub_samples,split = ",",fixed=T)))
	print(paste0("所有细胞群：",paste(unique(rds1$orig.ident),collapse=",")   ))
	print(paste0("subset的细胞群：",paste(subsamples,collapse=",")   ))
	if (!all(subsamples %in% unique(rds1$orig.ident)  )){print("抱歉，您选取了不存在的样本，请仔细检查config配置文件......");quit()}
	rds<-subset(rds1,orig.ident==subsamples)
	}else{
	rds<-rds1
	}


#===========================================================
print("###........................................")
print("开始亚群细分......")
print(Sys.time())
print("")
#开始进行亚群细分
#nfeatures<-2000;npc<-30;dims<-20 #默认一般参数
if ('seurat_clusters' %in% colnames(rds@meta.data)){
    rds@meta.data$seurat_clusters_raw <- rds@meta.data$seurat_clusters
}
DefaultAssay(rds) <- "RNA"
#https://satijalab.org/seurat/articles/integration_rpca.html
if (reclusters_method == 'rpca'){
	ifnb.list <- SplitObject(rds, split.by = "orig.ident")
	# normalize and identify variable features for each dataset independently
	ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
		x <- NormalizeData(x)
		x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures)
	})
	features <- SelectIntegrationFeatures(object.list = ifnb.list)
	ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
		x <- ScaleData(x, features = features, verbose = FALSE)
		x <- RunPCA(x, features = features, verbose = FALSE)
	})
	immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
	# this command creates an 'integrated' data assay
	immune.combined <- IntegrateData(anchorset = immune.anchors)
	# specify that we will perform downstream analysis on the corrected data note that the
	# original unmodified data still resides in the 'RNA' assay
	DefaultAssay(immune.combined) <- "integrated"
	# Run the standard workflow for visualization and clustering
	immune.combined <- ScaleData(immune.combined, verbose = FALSE)
	immune.combined <- RunPCA(immune.combined, npcs = npc, verbose = FALSE)
	# p <- ElbowPlot(immune.combined)
	# ggsave(paste0(prefix,'_PCA_ElbowPlot.pdf'), w=6,h=5)
	immune.combined <- RunTSNE(immune.combined, reduction = "pca",dims = 1:dims)
	immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:dims)
	tmp <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:dims)
	
}else if (reclusters_method == 'harmony'){
	#nfeatures<-2000;npc<-30;dims<-20 #默认一般参数
	tmp <- NormalizeData(rds) %>% FindVariableFeatures(nfeatures =nfeatures) %>% ScaleData() %>% RunPCA(npcs = as.numeric(npc),verbose = FALSE)
	tmp <- RunHarmony(tmp, group.by.vars = "orig.ident")
	tmp <- RunTSNE(tmp, reduction = "harmony", dims = 1:dims)
	tmp <- RunUMAP(tmp, reduction = "harmony", dims = 1:dims)
	tmp <- FindNeighbors(tmp, reduction = "harmony", dims = 1:dims)#harmony
#CCA方法，亚群细分时一般不建议
}else if (reclusters_method == 'cca'){
	ifnb.list <- SplitObject(rds, split.by = "orig.ident")
	# normalize and identify variable features for each dataset independently
	ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
		x <- NormalizeData(x)
		x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures)
	})
	features <- SelectIntegrationFeatures(object.list = ifnb.list)
	# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
		# x <- ScaleData(x, features = features, verbose = FALSE)
		# x <- RunPCA(x, features = features, verbose = FALSE)
	# })
	immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
	# this command creates an 'integrated' data assay
	immune.combined <- IntegrateData(anchorset = immune.anchors)
	# specify that we will perform downstream analysis on the corrected data note that the
	# original unmodified data still resides in the 'RNA' assay
	DefaultAssay(immune.combined) <- "integrated"
	# Run the standard workflow for visualization and clustering
	immune.combined <- ScaleData(immune.combined, verbose = FALSE)
	immune.combined <- RunPCA(immune.combined, npcs = npc, verbose = FALSE)
	# p <- ElbowPlot(immune.combined)
	# ggsave(paste0(prefix,'_PCA_ElbowPlot.pdf'), w=6,h=5)
	immune.combined <- RunTSNE(immune.combined, reduction = "pca",dims = 1:dims)
	immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:dims)
	tmp <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:dims)

}else if (reclusters_method == 'pca'){
	tmp <- NormalizeData(rds) %>% FindVariableFeatures(nfeatures =nfeatures) %>% ScaleData() %>% RunPCA(npcs = as.numeric(npc),verbose = FALSE)
	#tmp <- RunHarmony(tmp, group.by.vars = "stim")
	tmp <- RunTSNE(tmp, reduction = "pca", dims = 1:dims)
	tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims)
	tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims)#harmony
}

	
#根据指定的分辨率分群
tmp <- FindClusters(tmp,resolution =res)
current.cluster.ids <- levels(Idents(tmp))
new.cluster.ids <- as.numeric(current.cluster.ids)+1
Idents(tmp) <- plyr::mapvalues(x = Idents(tmp), from = current.cluster.ids, to = new.cluster.ids)
tmp@meta.data$seurat_clusters<-Idents(tmp)[rownames(tmp@meta.data)]
table(tmp@meta.data$seurat_clusters)
	
# p <- clustree(immune.combined, prefix = "integrated_snn_res.")
# ggsave(paste0(prefix,'_clustree.pdf'), w=8,h=8)

setwd(outdir)	
saveRDS(tmp, paste0(prefix,'.rds'))
DefaultAssay(tmp) <- "RNA"
##===========================================================
mkdirs(outdir,'1_QC')
outdir_pre<-paste0(outdir,'/1_QC/')
setwd(outdir_pre)
#1_QC #绘制质控图
#percent.mt','percent.HB','nCount_RNA','nFeature_RNA',"S.Score","G2M.Score"
#qc_list<-c('percent.mt','percent.HB','nCount_RNA','nFeature_RNA',"S.Score","G2M.Score")
qc_list<-unique(unlist(strsplit(ini.list$Subclusters$qc_list,split = ",",fixed=T)))
for (qc in qc_list){
	#FeaturePlot
	p<-FeaturePlot(tmp,features=qc,cols = c("lightgrey", "red"),label=F,pt.size=0,label.size=3)
	pdf(paste0(qc,'_FeaturePlot.pdf'),w=5,h=4)
	print(p)
	dev.off()
	#VlnPlot
	p<-VlnPlot(tmp, features = qc,group.by='seurat_clusters',pt.size = -1,slot="data",assay="RNA",stack = F,flip =TRUE)+theme(legend.position='none',axis.text.y = element_text(color="black",size=8,angle=0),axis.text.x = element_text(color="black",size=10,angle=30),axis.title.x=element_text(size = 10),axis.title.y=element_text(size = 8))+xlab('')
	pdf(paste0(qc,'_clusters_VlnPlot.pdf'),w=2+0.2*length(unique(tmp$seurat_clusters)),h=2)
	print(p)
	dev.off()
}
##===========================================================
mkdirs(outdir,'2_clusters')
outdir_pre<-paste0(outdir,'/2_clusters/')
setwd(outdir_pre)	
#1_clusters #绘制分群结果
p1 <- DimPlot(tmp, reduction = "umap", label = T,group.by = "seurat_clusters",label.size = 4,pt.size = 0)+ggtitle('')
p2 <- DimPlot(tmp, reduction = "umap", label = F,group.by = "stim",label.size = 4,pt.size = 0)+ggtitle('')
p<-plot_grid(p1,p2,ncol = 2)
pdf(paste0(prefix,'_clusters_sample_umap.pdf'),w=12,h=5)
print(p)
dev.off()

# n1<-length(unique(tmp$stim))
# p3 <- DimPlot(tmp, reduction = "umap", label = T,group.by = "seurat_clusters",split.by='stim',label.size = 4,pt.size = 0,ncol=2)
# pdf(paste0(prefix,'_clusters_splitby-sample_umap.pdf'),w=12,h=5)
# print(p3)
# dev.off()

# p1 <- DimPlot(tmp, reduction = "tsne", label = T,group.by = "seurat_clusters",label.size = 4,pt.size = 0)+ggtitle('')
# p2 <- DimPlot(tmp, reduction = "tsne", label = T,group.by = "stim",label.size = 4,pt.size = 0)+ggtitle('')
# p<-plot_grid(p1,p2,ncol = 2)
# pdf(paste0(prefix,'clusters_sample_tsne.pdf'),w=12,h=5)
# print(p)
# dev.off()
	
#细胞占比分析
tmp$celltype<-tmp$seurat_clusters;tmp$Sample<-tmp$stim
a<-table(tmp$Sample,tmp$celltype)
write.csv(a,file=paste0(prefix,'_sample_cellnum.csv'),quote=F)

ratio<-round(table(tmp$Sample,tmp$celltype)/as.vector(table(tmp$Sample)),4)
write.csv(ratio,file=paste0(prefix,'_sample_cellratio.csv'),quote=F)

a<-data.frame(CellType=tmp@meta.data[,'celltype'],Samples=tmp@meta.data[,'Sample'])
a$CellType<-factor(a$CellType,levels=sort(unique(a$CellType)))
p0<-theme(legend.position='right',panel.grid=element_blank(), legend.background = element_rect(colour = NA),
		legend.title = element_blank(),legend.text = element_text(face="plain", color="black",size=8),
		axis.text.x = element_text(color="black",size=10,angle=30, hjust = 1),
		axis.text.y = element_text(color="black",size=10),
		axis.title.x = element_text(face="plain", color="black",size=12),
		axis.title.y = element_text(face="plain", color="black",size=12))
pdf(paste(prefix,"clusters.stats.pdf",sep='_'),w=6,h=5)
p<-ggplot(a, aes(Samples)) + geom_bar(aes(fill=CellType), position='fill',width=0.6)+labs(x=" ", y = "Cell type distribution",fill= "CellType")+theme_bw()+p0
print(p)
dev.off()

pdf(paste(prefix,"samples.stats.pdf",sep='_'),w=8,h=5)
p<-ggplot(a, aes(CellType)) + geom_bar(aes(fill=Samples), position='fill',width=0.6)+labs(x=" ", y = "Cell type distribution",fill= "Samples")+theme_bw()+p0
print(p)
dev.off()


##===========================================================
print("###........................................")
print("开始markers基因分析......")
print(Sys.time())
print("")
mkdirs(outdir,'3_markers')
outdir_pre<-paste0(outdir,'/3_markers/')
setwd(outdir_pre)	
DefaultAssay(tmp) <- "RNA"
#marker基因
Idents(tmp)<-'seurat_clusters'
cluster.averages <- AverageExpression(object = tmp, assays ='RNA',return.seurat = F)
write.csv(cluster.averages$RNA,file=paste(prefix,'clusters_averageExpression.csv',sep='_'),quote=F,row.names=T)
average.genes<-cluster.averages$RNA
all.markers<-FindAllMarkers(tmp, only.pos = T, min.pct =0.1, logfc.threshold = 0.25,test.use='wilcox')
markers<-subset(all.markers,p_val_adj<0.05)
da<-arrange(markers,cluster,desc(avg_log2FC))
write.csv(da,file=paste0(prefix,'_all.markers_FDR0.05.csv'),quote=F,row.names=T)
#da1<-da[-grep("^mt-",da$gene),]
topn<-10
top_gene <- da %>% group_by(cluster) %>% top_n(topn,avg_log2FC)
filename<-paste0(prefix,'_top10_markers_genes_pheatmap.pdf')
genes<-unique(as.vector(top_gene$gene))
data<-average.genes
da0<-log2(data[genes,]+1)
n1<-length(unique(tmp$seurat_clusters))
pheatmap::pheatmap(da0,filename=filename,width=4+0.3*n1,height=3+0.8*n1,cluster_rows=F,cluster_cols=F,display_numbers = F,number_format = "%.0f",fontsize = 8,fontsize_col = 10,border_color=NA,angle_col = "45",scale='row')

#top5
da %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
#options (repr.plot.width=28, repr.plot.height=5)
tmp <- ScaleData(tmp, features = genes, verbose = FALSE)
pdf(paste0(prefix,'_top5_markers_genes_DotPlot.pdf'), w=10+n1,h=6+0.2*n1)
DotPlot(tmp,group.by='seurat_clusters', features=unique(top5$gene),cols = c("lightgrey", "red"))+theme(axis.text.x = element_text(size=10,angle =60, hjust = 1),legend.text = element_text(face="plain", color="black",size=8),legend.title =element_text(face="plain", color="black",size=10))+xlab('')+ylab('')
dev.off()

##===========================================================
print("###........................................")
print("完成所有分析......")
print(Sys.time())
