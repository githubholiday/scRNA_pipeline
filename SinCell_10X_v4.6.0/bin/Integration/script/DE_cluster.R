
#名称：DE_cluster.R
#作者：mengli
#邮箱：mengli@genome.cn
#时间：20231130
#版本：v0.0.3
#用途：利用seurat对每个cluster对不同的比较组来源的细胞进行统计


library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'inrds',	'i',	1,	"character",
	'config',	'c',	1,	"character",
	'outdir',	'o',	1,	"character"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	==========================================================
	inrds: rds file
	prefix:the output prefix of files,such as pictures and excel.
	config 配置文件，配置文件中包含样品分组，样品差异分析安排，以及其他的参数
	==========================================================
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
if ( is.null(opt$inrds) )	{ cat("Please input the data file1 ...\n\n") ; print_usage(para)}
if ( is.null(opt$config) )	{ cat("Please input the data file2 ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )	{ cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para) }
require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(harmony)
mkdirs <- function(outdir, fp){
	if(!file.exists(file.path(outdir,fp))) {
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
reduction<-function(immune.combined,dims_num=20,resolution=0.8,outdir=getwd(),pref='10x',w_h=c(24,8)){
	x_umap <-c()
	y_umap <-c()
	p_tmp<<-theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
        legend.title = element_blank(),legend.text =  element_text(color="black",size=30),
        axis.text.x = element_text(color="black",size=30),
        axis.text.y = element_text(color="black",size=30),
        axis.title.x = element_text(face="plain", color="black",size=30),
        axis.title.y = element_text(face="plain", color="black",size=30))
	pdf(paste(outdir,paste(pref,"umap_cluster_samples.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
	# Visualization
	p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim",pt.size = 1)+p_tmp
	p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE,label.size = 5,pt.size = 1,repel = T)+p_tmp
	p3<-plot_grid(p1, p2)
	print(p3)
	p4 <- DimPlot(immune.combined, reduction = "umap", label = TRUE,label.size = 5,pt.size = 1,split.by = "stim",repel = T)+p_tmp
	print(p4)
	dev.off()
	pdf(paste(outdir,paste(pref,"umap_cluster_groups.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
	p5 <- DimPlot(immune.combined, reduction = "umap", label = TRUE,label.size = 5,pt.size = 1,split.by = "Group",repel = T)+p_tmp
	print(p5)
	dev.off()
	x_umap <-as.data.frame(immune.combined@reductions$umap@cell.embeddings)
	x_umap$orig.ident<-rownames(x_umap)
	res <-immune.combined@meta.data
	y_umap <-data.frame(orig.ident=rownames(res),sample=res$stim,cluster=res$seurat_clusters)
	c<-merge(x_umap,y_umap,by='orig.ident')
	print("保存细胞的umap坐标结果：cell.umap.csv")
	write.csv(c,paste(outdir,"cell.umap.csv",sep='/'),quote=F,row.names=F)
	
	#各样本cluster统计
	a<-data.frame(Clusters=immune.combined@meta.data[,'seurat_clusters'],Samples=immune.combined@meta.data[,'orig.ident'])
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

	#各比较组细胞类型统计
	b<-data.frame(Clusters=immune.combined@meta.data[,'seurat_clusters'],Group=immune.combined@meta.data[,'Group'])
	b$Clusters<-factor(b$Clusters,levels=sort(unique(b$Clusters)))
	if (length(unique(b$Group)) > 1) {
	pdf(paste(outdir,paste(pref,"group.clusters_stats.pdf",sep='_'),sep='/'),w=12,h=8)
	p<-ggplot(b, aes(Clusters)) + geom_bar(aes(fill=Group), position='fill',width=0.6)+labs(x=" ", y = "",fill= "Group")+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.background = element_rect(colour = NA),plot.title=element_text(size = 25),axis.text.x=element_text(size=20,angle=60,hjust=1),axis.text.y=element_text(size=20),axis.title.x=element_text(size = 25),axis.title.y=element_text(size = 25),legend.text =element_text(size = 25),legend.title =element_text(size = 25))
	print(p)
	dev.off()
	}
	b1 <- get_table(as.data.frame.array(table(b)))
	b2 <- data.frame(Clusters=rownames(b1), b1)
	write.table(b2,paste(pref,'group.clusters.xls',sep='_'),quote=F,sep="\t",row.names=F)

	#####tsne分析
	pdf(paste(outdir,paste(pref,"tsne_cluster_samples.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
	# Visualization
	p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "stim",pt.size = 1)+p_tmp
	p2 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE,label.size = 5,pt.size = 1, repel = T)+p_tmp
	p3<-plot_grid(p1, p2)
	print(p3)
	p4 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE,label.size = 5,pt.size = 1,split.by = "stim",repel = T)+p_tmp
	print(p4)
	dev.off()
	pdf(paste(outdir,paste(pref,"tsne_cluster_groups.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
	# Visualization
	p5 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE,label.size = 5,pt.size = 1,split.by = "Group",repel = T)+p_tmp
	print(p5)
	dev.off()
	return (immune.combined)

}



#主流程
#2_clusters  3_marker
prefix<-opt$prefix
outdir<-opt$outdir
indir<-opt$inrds
ini<-opt$config
ini.list <- read.config(file = ini)

#组合1	N	N2THY/N1THY/N3THY
#组合2	H	H2THY/H3THY/H1THY
#ini.list$sample$sample1  unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T))
sample_name<-unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T))
immune.combined1<-readRDS(indir)

#聚类
mkdirs(outdir,'2_clusters')
setwd(paste(outdir,'2_clusters',sep='/'))
sample1<-unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T)) #样本名称
sample2<-unlist(strsplit(ini.list$sample$sample2,split = "/",fixed=T)) #分组名称

#按照样本名称提取出对应的细胞-limeng

immune.combined= immune.combined1[,immune.combined1@meta.data$orig.ident %in% c(sample1)]

print(unique(immune.combined@meta.data$seurat_clusters))
immune.combined$stim <- immune.combined$orig.ident
immune.combined$Group <- immune.combined$orig.ident
for (n in 1:length(sample1)){
	immune.combined$Group<-gsub(paste('\\b',sample1[n],'\\b',sep=''),sample2[n],immune.combined$Group)
}
immune.combined<-reduction(immune.combined,dims_num=as.numeric(ini.list$Para$reduction_dims_num),resolution=as.numeric(ini.list$Para$reduction_resolution),outdir=paste(outdir,'2_clusters',sep='/'),pref=prefix)
setwd(paste(outdir,'2_clusters',sep='/'))
saveRDS(immune.combined, file = paste(prefix,'immune_combined.rds',sep='_'))
