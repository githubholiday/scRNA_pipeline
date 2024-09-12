#名称：cluster_multi.R
#作者：mengli
#邮箱：mengli@genome.cn
#时间：20231130
#版本：v0.0.3
#用途：对所有样本合并后的rds进行marker基因分析。
#===========================================================
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
#library(SeuratWrappers)

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


marker_gene <- function(immune.combined,min.pct = 0.1, logfc.threshold = 0.25,outdir=getwd(),pref='10x',test.use= "wilcox"){
	DefaultAssay(immune.combined) <- "RNA"
	cluster_num<-sort(unique(Idents(immune.combined)))
	print("Starting FindMarkers")
	all.markers <- FindAllMarkers(immune.combined, only.pos = FALSE, min.pct =min.pct, logfc.threshold = logfc.threshold,test.use=test.use)
	#all.markers <- RunPrestoAll(immune.combined,only.pos = FALSE, min.pct =min.pct, logfc.threshold = logfc.threshold,test.use=test.use)

	head(all.markers, n = 2)
	marker_result <- subset(all.markers, select=c(7,1,2,3,4,5,6))
	names(marker_result)[names(marker_result)=='gene'] <- 'gene_name'
	#write.csv(all.markers,paste(outdir,'3_marker',paste(pref,'all.markers.csv',sep='_'),sep='/'),quote=F)
	write.csv(marker_result,paste(outdir,'3_marker',paste(pref,'all.markers.csv',sep='_'),sep='/'),quote=F,row.names=F)
	print("Finished FindMarkers")
	return (all.markers)
}

marker_plot<-function(immune.combined,all.markers,outdir=getwd(),pref='10x',heatmap_n=10){
top1 <- all.markers %>% group_by(cluster) %>% top_n(1, avg_log2FC)
top1_gene<-unique(top1$gene)
pdf(paste(outdir,paste(pref,'top1markers_FeaturePlot.pdf',sep='_'),sep='/'),w=12,h=0.6*length(top1_gene)+3)
p<-FeaturePlot(immune.combined, features = as.vector(unique(top1_gene)), min.cutoff = "q9",cols = c("lightgrey", "red"))
print(p)
dev.off()
all.genes <- rownames(immune.combined)
immune.combined <- ScaleData(immune.combined, features = all.genes)
small<-immune.combined
#gene<-as.vector(unique(top10$gene))
exprs <- data.frame(FetchData(object = small, vars = as.vector(top1_gene)))
exprs$Barcod<-rownames(exprs)
ident<-data.frame()
#barcode与聚类信息提取
ident<-data.frame(Barcod=rownames(small@meta.data),orig.ident=small@meta.data$seurat_clusters,samples=small@meta.data$stim)
#通过merge函数，将表达量与聚类号对应起来
c<-merge(exprs,ident,by='Barcod')
#对其进行排序
#c$orig.ident<-factor(c$orig.ident,levels=new.ident)
a<-data.frame(Barcod=0,Clusters=0,samples=0,Exp=0,Gene=0)
for (g in colnames(c)[4:ncol(c)-2]){
	tmp<-c[,c('Barcod','orig.ident','samples',g)]
	tmp$Gene<-rep(g, times=nrow(tmp))
	colnames(tmp)<-c('Barcod','Clusters','samples','Exp','Gene')
	a<-rbind(a,tmp)
}
a<-a[2:nrow(a),]
noise <- rnorm(n = length(x = a[,c('Exp')])) / 100000
a[,c('Exp')] <- a[, c('Exp')] + noise
mkdirs(outdir,'marker_violin')
for (g in unique(a$Gene)){
	b<-a[a$Gene==g,]
	b$Clusters<-factor(b$Clusters,levels=as.character(sort(unique(b$Clusters))))
	pdf(paste(outdir,'marker_violin',paste(g,'clusters.pdf',sep='_'),sep='/'),w=0.4*length(unique(Idents(immune.combined))),h=0.2*length(top1_gene))
	p2 <- ggplot(b, aes(x = Clusters, y = Exp, fill = Clusters)) + geom_violin(scale = "width",adjust =1) +labs(x = '', y = 'Expression of Gene',title=g) +theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank(),legend.position="none",axis.text.x = element_text(color="black",size=0.7*length(unique(Idents(immune.combined))),angle = 60, hjust = 1),axis.text.y = element_text(color="black",size=0.7*length(unique(Idents(immune.combined)))),axis.title.x = element_text(face="plain", color="black",size=0.7*length(unique(Idents(immune.combined)))),axis.title.y = element_text(face="plain", color="black",size=0.7*length(unique(Idents(immune.combined)))))+theme(plot.title=element_text(size=1.1*length(unique(Idents(immune.combined))),color="black",hjust = 0.5))
	print(p2)
	dev.off()
}
	nrow=ceiling(length(unique(a$Gene))/4)
	a$Clusters<-factor(a$Clusters,levels=as.character(sort(as.numeric(unique(a$Clusters)))))
	pdf(paste(outdir,paste(pref,'all','clusters.pdf',sep='_'),sep='/'),w=2.2*length(unique(Idents(immune.combined))),h=1*length(top1_gene))
	p2 <- ggplot(a, aes(x = Clusters, y = Exp, fill = Clusters)) + geom_violin(scale = "width",adjust =1) + facet_wrap(~Gene, nrow =nrow,ncol =4,scales = 'free',strip.position='top') +labs(x = 'Clusters', y = 'Expression of Gene',title=' ') +theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank(),strip.text = element_text(size = 1.6*length(unique(Idents(immune.combined)))),strip.background = element_rect(fill = NA, colour = NA),strip.placement = "inside",legend.position="none",axis.text.x = element_text(color="black",size=1.2*length(unique(Idents(immune.combined))),angle = 60, hjust = 1),axis.text.y = element_text(color="black",size=1.2*length(unique(Idents(immune.combined)))),axis.title.x = element_text(face="plain", color="black",size=1.2*length(unique(Idents(immune.combined)))),axis.title.y = element_text(face="plain", color="black",size=1.2*length(unique(Idents(immune.combined)))))
	print(p2)
	dev.off()
	top10 <- all.markers %>% group_by(cluster) %>% top_n(as.numeric(heatmap_n), avg_log2FC)
	pdf(paste(outdir,paste(pref,'top10_marker_heatmap.pdf',sep='_'),sep='/'),w=1*length(unique(Idents(immune.combined))),h=0.15*length(unique(top10$gene)))
	p<-DoHeatmap(small, features = unique(top10$gene)) + NoLegend()
	print(p)
	dev.off()
	
	pdf(paste(outdir,paste(pref,'top_dotplot.pdf',sep='_'),sep='/'),w=0.5*length(unique(Idents(immune.combined)))+4,h=0.5*length(top1_gene)+4)
	p<-DotPlot(object = small, features=top1_gene,cols = c("blue", "red"))+RotatedAxis()+coord_flip()+
	labs(x="Marker",y="")+
	theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
	print(p)
	dev.off()
}


#主流程
#2_clusters  3_marker
prefix<-opt$prefix
outdir<-opt$outdir
indir<-opt$inrds
ini<-opt$config
ini.list <- read.config(file = ini)

sample_name<-unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T))
immune.combined<-readRDS(indir)

#Marker基因分析
DefaultAssay(immune.combined) <- "RNA"
mkdirs(outdir,'3_marker')
setwd(paste(outdir,'3_marker',sep='/'))
print(Sys.time())
all.markers<-marker_gene(immune.combined,min.pct =as.numeric(ini.list$Para$marker_gene_min.pct), logfc.threshold = as.numeric(ini.list$Para$marker_gene_logfc.threshold),outdir=outdir,test.use=ini.list$Para$marker_gene_test.use,pref=prefix)
immune.combined@misc$markers <- all.markers
saveRDS(immune.combined, file = paste(prefix,"marker.rds",sep='_'))

marker_plot(immune.combined,all.markers,outdir=paste(outdir,'3_marker',sep='/'),pref=prefix)
print(Sys.time())
mkdirs(outdir,'2_Com_clusters')
setwd(paste(outdir,'2_Com_clusters',sep='/'))
saveRDS(immune.combined, file = paste(prefix,"immune_combined_marker.rds",sep='_'))
