#名称：Integration_rds.R
#作者：mengli
#邮箱：mengli@genome.cn
#时间：20231130
#版本：v0.0.3
#用途：将项目中所有的样本对过滤后的rds文件进行合并、聚类、分群。
#===========================================================
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'indir',	'i',	1,	"character",
	'config',	'c',	1,	"character",
	'outdir',	'o',	1,	"character"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	========================================================================================================================================
	indir数目目录:
	可以是两种形式：
	一种是我们的标准分析的目录：CellRanger_Count，比如/****/PM-JL190104-05/std/wangxiao/Analysis-test/Analysis/CellRanger_Count/
	另外一种是数目目录下全部为输入数据，需要注意文件名格式，必须paste(sample_names,'_all_UMI.csv',sep='')
	========================================================================================================================================
	prefix:the output prefix of files,such as pictures and excel.
	========================================================================================================================================
	config 配置文件，配置文件中包含样品分组，样品差异分析安排，以及其他的参数
	========================================================================================================================================
	outdir:outdir  of outputs,we will setwd(opt$outdir)
	Usage example:
	Rscript this.r -i1 BM-aLP_all_UMI.csv -i2 FL-alp_all_UMI.csv -o outdir -s1 BM-alp -s2 FL-alp -p alp
	Options:
	--help		h	NULL		get this help
	--indir	i	character	indir for expression file[forced]
	--config	c	character	config.ini file for group and other Para[forced]
	--outdir	o	character	The	resurt of out dir for analysis [forced]
	--prefix	p	character	the prefix for outputfiles [forced]
	\n")
	q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$indir) )	{ cat("Please input the data file1 ...\n\n") ; print_usage(para)}
if ( is.null(opt$config) )	{ cat("Please input the data file2 ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )	{ cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para) }
#if ( is.null(opt$species) )	{ cat("Please give the species ...\n\n") ; print_usage(para) }
##这个分析用最新的seurat包进行分析
require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(harmony)
mkdirs <- function(outdir,fp) {
	if(!file.exists(file.path(outdir,fp))) {
#		mkdirs(dirname(fp))
		dir.create(file.path(outdir,fp))
	}else{
			print(paste(fp,"Dir already exists!",sep="     "))
			#unlink(file.path(outdir,fp), recursive=TRUE)
			#dir.create(file.path(outdir,fp))
		}
}
get_table <- function(df1){
	row <- dim(df1)[1]
	col <- dim(df1)[2]
	df2 <- df1
	for (c in 1:col){
  		all <- sum(as.numeric(df1[,c]))
  		for (r in 1:row){
    			s <- round(as.numeric(df1[r,c])/all*100,2)
    			df2[r,c] <- paste(df1[r,c],'(',s,'%)',sep='')
 		}
	}
	return(df2)
}

##1_QC  2_clusters  3_marker  4_conserved_markers  5_diff_gene_condition
sample2SeuratObject_list <- function(indir,oudir,samplenames){
  print("start read")
  SeuratObject_list <- mapply(read_onesample2SeuratObject, indir, oudir, samplenames)
  print (length(SeuratObject_list))
  print("end read")
  names(SeuratObject_list) <- samplenames
  return(SeuratObject_list)
}
read_onesample2SeuratObject <- function(indir,outdir,sample_names,min.cells=as.numeric(ini.list$Para$object_list_min.cells),mitoName=ini.list$Para$object_list_mitoname,HB=ini.list$Para$object_list_hb,mt.percent=as.numeric(ini.list$Para$object_list_mt.percent),min_nFeature_RNA=as.numeric(ini.list$Para$object_list_min_nfeature_rna),max_nFeature_RNA=as.numeric(ini.list$Para$object_list_max_nfeature_rna),normalization.method=ini.list$Para$object_list_normalization.method,scale.factor=as.numeric(ini.list$Para$object_list_scale.factor),nfeatures_FindVariableFeatures=as.numeric(ini.list$Para$object_list_nfeatures_findvariablefeatures),FindVariableFeatures_method=ini.list$Para$object_list_findvariablefeatures_method){
  print(paste('Start read rawdata, the sample is :',sample_names,sep=' '))
  sample_indir<- file.path(indir, sample_names, paste(sample_names, '_filter_cell.rds', sep=''))
  print(sample_indir)
  SeuratObject <- readRDS(sample_indir)
  SeuratObject$stim <- sample_names
  #####
  print(paste('Finished read rawdata, the sample is :',sample_names,sep=' '))
  return(SeuratObject)
}
qc_pca_plot<-function(immune.combined,outdir=getwd(),pref='10x',w_h=c(12,8)){
    # plot pca
    pbmc <- JackStraw(immune.combined, num.replicate = 100)
    pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
    pdf(paste(outdir,paste(pref,"pca.qc.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
    p1<-DimPlot(immune.combined, reduction = "pca")
    # plot pca
#    pdf(paste(outdir,paste(pref,"pca_heatmap.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
    p2<-DimHeatmap(immune.combined, dims = 1:15, cells = 500, balanced = TRUE)
#    dev.off()
    # plot pca
#    pdf(paste(outdir,paste(pref,"pca_ElbowPlot.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
    p3<-VizDimLoadings(immune.combined, dims = 1:2, reduction = "pca")
    p4<-JackStrawPlot(pbmc, dims = 1:15)
    p5<-ElbowPlot(immune.combined)
	print(p1)
	print(p2)
	print(p3)
	print(p4)
	print(p5)
    dev.off()
}

reduction<-function(immune.combined,dims_num=20,resolution=0.8,outdir=getwd(),pref='10x',w_h=c(24,8)){
	x_umap <-c()
	y_umap <-c()
	if ( length(unique(immune.combined@meta.data$orig.ident)) > 1 ){
		immune.combined <- NormalizeData(immune.combined) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
		immune.combined <- RunHarmony(immune.combined, group.by.vars = "stim")
		immune.combined_umap <- RunUMAP(immune.combined, reduction = "harmony",dims = 1:as.numeric(dims_num))
		immune.combined_umap <- FindNeighbors(object = immune.combined_umap, reduction = "harmony", dims = 1:as.numeric(dims_num))
	}else{
		immune.combined_umap <- RunUMAP(immune.combined, reduction = "pca",dims = 1:as.numeric(dims_num))
		immune.combined_umap <- FindNeighbors(object = immune.combined_umap, reduction = "pca", dims = 1:as.numeric(dims_num))
	}
	immune.combined_umap <- FindClusters(immune.combined_umap, resolution = as.numeric(resolution))
	levels(Idents(immune.combined_umap))
	#改聚类号，从1开始
	current.cluster.ids <- levels(Idents(immune.combined_umap))
	new.cluster.ids <- as.numeric(current.cluster.ids)+1
	Idents(immune.combined_umap) <- plyr::mapvalues(x = Idents(immune.combined_umap), from = current.cluster.ids, to = new.cluster.ids)
	immune.combined_umap$seurat_clusters <- Idents(immune.combined_umap)
	#
	p_tmp<<-theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
        legend.title = element_blank(),legend.text =  element_text(color="black",size=30),
        axis.text.x = element_text(color="black",size=30),
        axis.text.y = element_text(color="black",size=30),
        axis.title.x = element_text(face="plain", color="black",size=30),
        axis.title.y = element_text(face="plain", color="black",size=30))
	pdf(paste(outdir,paste(pref,"umap_cluster_samples.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
	# Visualization
	p1 <- DimPlot(immune.combined_umap, reduction = "umap", group.by = "stim",pt.size = 1)+p_tmp+ggtitle("")
	p2 <- DimPlot(immune.combined_umap, reduction = "umap", label = TRUE,label.size = 5,pt.size = 1,repel = T)+p_tmp
	p3<-plot_grid(p1, p2)
	print(p3)
	p4 <- DimPlot(immune.combined_umap, reduction = "umap", label = TRUE,label.size = 5,pt.size = 1,split.by = "stim",repel = T)+p_tmp
	print(p4)
	dev.off()
	x_umap <-as.data.frame(immune.combined_umap@reductions$umap@cell.embeddings)
	x_umap$orig.ident<-rownames(x_umap)
	res <-immune.combined_umap@meta.data
	y_umap <-data.frame(orig.ident=rownames(res),sample=res$stim,cluster=res$seurat_clusters)
	c<-merge(x_umap,y_umap,by='orig.ident')
	print("保存细胞的umap坐标结果：cell.umap.csv")
	write.csv(c,paste(outdir,"cell.umap.csv",sep='/'),quote=F,row.names=F)
	
	#各样本cluster统计
	a<-data.frame(Clusters=immune.combined_umap@meta.data[,'seurat_clusters'],Samples=immune.combined_umap@meta.data[,'orig.ident'])
	a$Clusters<-factor(a$Clusters,levels=sort(unique(a$Clusters)))
	if (length(unique(a$Samples)) > 1) {
	pdf(paste(outdir,paste(pref,"sample.clusters_stats.pdf",sep='_'),sep='/'),w=12,h=8)
	p<-ggplot(a, aes(Clusters)) + geom_bar(aes(fill=Samples), position='fill',width=0.6)+labs(x=" ", y = "",fill= "Samples")+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.background = element_rect(colour = NA),plot.title=element_text(size = 25),axis.text.x=element_text(size=20,angle=60,hjust=1),axis.text.y=element_text(size=20),axis.title.x=element_text(size = 25),axis.title.y=element_text(size = 25),legend.text =element_text(size = 25),legend.title =element_text(size = 25))
	print(p)
	dev.off()
	}
	a1 <- get_table(as.data.frame.array(table(a)))
	a2 <- data.frame(Clusters=rownames(a1), a1)
	write.table(a2,paste(pref,'sample.clusters.xls',sep='_'),quote=F,sep="\t", row.names=F)


	#####tsne分析
	x_tsne <-c()
	y_tsne <-c()
	if ( length(unique(immune.combined@meta.data$orig.ident)) > 1 ){
		immune.combined_tsne <- RunTSNE(immune.combined_umap, reduction = "harmony",dims = 1:as.numeric(dims_num), check_duplicates = FALSE)
		immune.combined_tsne  <- FindNeighbors(object = immune.combined_tsne , reduction = "harmony", dims = 1:as.numeric(dims_num))
	}else{
		immune.combined_tsne <- RunTSNE(immune.combined_umap, reduction = "pca",dims = 1:as.numeric(dims_num))
		immune.combined_tsne  <- FindNeighbors(object = immune.combined_tsne , reduction = "pca", dims = 1:as.numeric(dims_num))
	}
	immune.combined_tsne  <- FindClusters(immune.combined_tsne , resolution = as.numeric(resolution))
	#改聚类号，从1开始
	current.cluster.ids <- levels(Idents(immune.combined_tsne))
	new.cluster.ids <- as.numeric(current.cluster.ids)+1
	Idents(immune.combined_tsne) <- plyr::mapvalues(x = Idents(immune.combined_tsne), from = current.cluster.ids, to = new.cluster.ids)
	immune.combined_tsne$seurat_clusters <- Idents(immune.combined_tsne)
	pdf(paste(outdir,paste(pref,"tsne_cluster_samples.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
	# Visualization
	p1 <- DimPlot(immune.combined_tsne, reduction = "tsne", group.by = "stim",pt.size = 1)+p_tmp
	p2 <- DimPlot(immune.combined_tsne, reduction = "tsne", label = TRUE,label.size = 5,pt.size = 1, repel = T)+p_tmp
	p3<-plot_grid(p1, p2)
	print(p3)
	p4 <- DimPlot(immune.combined_tsne, reduction = "tsne", label = TRUE,label.size = 5,pt.size = 1,split.by = "stim",repel = T)+p_tmp
	print(p4)
	dev.off()
	x_tsne <-as.data.frame(immune.combined_tsne@reductions$tsne@cell.embeddings)
	x_tsne$orig.ident<-rownames(x_tsne)
	res <-immune.combined_tsne@meta.data
	y_tsne <-data.frame(orig.ident=rownames(res),sample=res$stim,cluster=res$seurat_clusters)
	c<-merge(x_tsne,y_tsne,by='orig.ident')
	print("保存细胞的tsne坐标结果：cell.tsne.csv")
	write.csv(c,"cell.tsne.csv",quote=F,row.names=F)
	
	return (immune.combined_tsne)
}

prefix<-opt$prefix
outdir<-opt$outdir
indir<-opt$indir
ini<-opt$config
ini.list <- read.config(file = ini)
mkdirs(outdir,'1_Com_QC')
setwd(paste(outdir,'1_Com_QC',sep='/'))
sample_name<-unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T))

if (length(sample_name) < 2){
    print("一个样本，不做合并分析")
    q()
}
object_list<-sample2SeuratObject_list(indir,paste(outdir,'1_Com_QC',sep='/'), sample_name)
immune.combined=merge(object_list[[1]],object_list[2:length(object_list)])
saveRDS(immune.combined, file = paste(prefix,'qc_before.rds',sep='_'))


immune.combined <- NormalizeData(object = immune.combined, normalization.method = ini.list$Para$object_list_normalization.method, scale.factor = as.numeric(ini.list$Para$object_list_scale.factor), verbose = FALSE)
immune.combined <- FindVariableFeatures(object = immune.combined, selection.method = ini.list$Para$object_list_findvariablefeatures_method, nfeatures = as.numeric(ini.list$Para$object_list_nfeatures_findvariablefeatures))
saveRDS(immune.combined, file = paste(prefix,'nor.rds',sep='_'))ll
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = as.numeric(ini.list$Para$integration_pca_dims), verbose = FALSE)
qc_pca_plot(immune.combined,outdir=paste(outdir,'1_Com_QC',sep='/'),pref=prefix,w_h=c(as.numeric(unlist(strsplit(ini.list$Para$qc_pca_plot_w_h,split = ",",fixed=T))[1]),as.numeric(unlist(strsplit(ini.list$Para$qc_pca_plot_w_h,split = ",",fixed=T))[2])))

#聚类
mkdirs(outdir,'2_Com_clusters')
setwd(paste(outdir,'2_Com_clusters',sep='/'))
immune.combined$stim <- immune.combined$orig.ident
immune.combined<-reduction(immune.combined,dims_num=as.numeric(ini.list$Para$reduction_dims_num),resolution=as.numeric(ini.list$Para$reduction_resolution),outdir=paste(outdir,'2_Com_clusters',sep='/'),pref=prefix)
saveRDS(immune.combined, file = paste(prefix,'cluster.rds',sep='_'))


