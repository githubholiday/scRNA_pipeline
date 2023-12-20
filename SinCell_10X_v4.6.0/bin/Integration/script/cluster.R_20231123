#!/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/Rscript
#名称：DE_10xGenomics.R
#作者：姚盟成
#邮箱：mengchengyao@genome.cn
#时间：201901011
#版本：v0.0.2
#用途：利用seurat进行10x 进行不同条件比较分析,需要输入配置文件，配置可以设置具有生物学重复的分析，指定分组。如果有
#两个以上样品，则需要分别取做差异分析。
###说明：
#程序开发环境/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/Rscript，需要指定R，指定包的路径,使用为seurat3.0版本以上
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
	p1 <- DimPlot(immune.combined_umap, reduction = "umap", group.by = "stim",pt.size = 1)+p_tmp
	p2 <- DimPlot(immune.combined_umap, reduction = "umap", label = TRUE,label.size = 5,pt.size = 1,repel = T)+p_tmp
	p3<-plot_grid(p1, p2)
	print(p3)
	p4 <- DimPlot(immune.combined_umap, reduction = "umap", label = TRUE,label.size = 5,pt.size = 1,split.by = "stim",repel = T)+p_tmp
	print(p4)
	dev.off()
	pdf(paste(outdir,paste(pref,"umap_cluster_groups.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
	p5 <- DimPlot(immune.combined_umap, reduction = "umap", label = TRUE,label.size = 5,pt.size = 1,split.by = "Group",repel = T)+p_tmp
	print(p5)
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

	#各比较组细胞类型统计
	b<-data.frame(Clusters=immune.combined_umap@meta.data[,'seurat_clusters'],Group=immune.combined_umap@meta.data[,'Group'])
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
	x_tsne <-c()
	y_tsne <-c()
	if ( length(unique(immune.combined@meta.data$orig.ident)) > 1 ){
		immune.combined_tsne <- RunTSNE(immune.combined, reduction = "harmony",dims = 1:as.numeric(dims_num), check_duplicates = FALSE)
		immune.combined_tsne  <- FindNeighbors(object = immune.combined_tsne , reduction = "harmony", dims = 1:as.numeric(dims_num))
	}else{
		immune.combined_tsne <- RunTSNE(immune.combined, reduction = "pca",dims = 1:as.numeric(dims_num))
		immune.combined_tsne  <- FindNeighbors(object = immune.combined_tsne , reduction = "pca", dims = 1:as.numeric(dims_num))
	}
	immune.combined_tsne  <- FindClusters(immune.combined_tsne , resolution = as.numeric(resolution))
	#改聚类号，从1开始
	current.cluster.ids <- levels(Idents(immune.combined_tsne))
	new.cluster.ids <- as.numeric(current.cluster.ids)+1
	Idents(immune.combined_tsne) <- plyr::mapvalues(x = Idents(immune.combined_tsne), from = current.cluster.ids, to = new.cluster.ids)
	
	pdf(paste(outdir,paste(pref,"tsne_cluster_samples.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
	# Visualization
	p1 <- DimPlot(immune.combined_tsne, reduction = "tsne", group.by = "stim",pt.size = 1)+p_tmp
	p2 <- DimPlot(immune.combined_tsne, reduction = "tsne", label = TRUE,label.size = 5,pt.size = 1, repel = T)+p_tmp
	p3<-plot_grid(p1, p2)
	print(p3)
	p4 <- DimPlot(immune.combined_tsne, reduction = "tsne", label = TRUE,label.size = 5,pt.size = 1,split.by = "stim",repel = T)+p_tmp
	print(p4)
	dev.off()
	pdf(paste(outdir,paste(pref,"tsne_cluster_groups.pdf",sep='_'),sep='/'),w=w_h[1],h=w_h[2])
	# Visualization
	p5 <- DimPlot(immune.combined_tsne, reduction = "tsne", label = TRUE,label.size = 5,pt.size = 1,split.by = "Group",repel = T)+p_tmp
	print(p5)
	dev.off()
	x_tsne <-as.data.frame(immune.combined_tsne@reductions$tsne@cell.embeddings)
	x_tsne$orig.ident<-rownames(x_tsne)
	res <-immune.combined_tsne@meta.data
	y_tsne <-data.frame(orig.ident=rownames(res),sample=res$stim,cluster=res$seurat_clusters)
	c<-merge(x_tsne,y_tsne,by='orig.ident')
	print("保存细胞的tsne坐标结果：cell.tsne.csv")
	write.csv(c,"cell.tsne.csv",quote=F,row.names=F)
	
	
	cluster_num<-length(unique(Idents(immune.combined_tsne)))
	color<-hue_pal()(cluster_num)
	for (i in 1:cluster_num){
		c_umap<-merge(x_umap,y_umap,by='orig.ident')
		c_umap$cluster<-as.character(as.numeric(c_umap$cluster))
		#c_umap$cluster <-as.character(c_umap$cluster)
		c_umap[c_umap$cluster!=i,]$cluster <- 'Not in'
		mycolo<-c(color[as.numeric(i+1)],'grey') 
		legend<-c(i,'Not in')
		graph <- paste(i,'grey.clusters.pdf',sep='_')
		pdf(graph, w=12, h=8)
		#tsne 高亮图
		c_umap %>%dplyr::group_by(cluster) %>% summarize(x = median(x = UMAP_1), y = median(x = UMAP_2)) -> centers
		p<-ggplot(c_umap,aes(x=c_umap$UMAP_1, y=c_umap$UMAP_2, colour=factor(c_umap$cluster,levels=legend)))+ geom_point(size=0.5)+scale_colour_manual(values=mycolo) +theme_bw(base_size = 20)+p_tmp+
		labs(x='UMAP1', y= 'UMAP2')+guides(colour = guide_legend(override.aes = list(size=8)))+geom_text(data=centers[centers$cluster!='Not in',],aes(x=x,y=y,label=cluster),colour='black',size = 10)
		lab <-paste("正在绘制",i,"高亮图。。。",sep=" ")
		print(lab)
		print(p)
#		dev.off()
		#tsne 高亮图
		c_tsne<-merge(x_tsne,y_tsne,by='orig.ident')
		c_tsne$cluster<-as.character(as.numeric(c_tsne$cluster))
		#c_tsne$cluster <-as.character(c_tsne$cluster)
		c_tsne[c_tsne$cluster!=i,]$cluster <- 'Not in'
		c_tsne %>%dplyr::group_by(cluster) %>% summarize(x = median(x = tSNE_1), y = median(x = tSNE_2)) -> centers
		p<-ggplot(c_tsne,aes(x=c_tsne$tSNE_1, y=c_tsne$tSNE_2, colour=factor(c_tsne$cluster,levels=legend)))+ geom_point(size=0.5)+scale_colour_manual(values=mycolo) +theme_bw(base_size = 20)+p_tmp+
		labs(x='TSNE1', y= 'TSNE2')+guides(colour = guide_legend(override.aes = list(size=8)))+geom_text(data=centers[centers$cluster!='Not in',],aes(x=x,y=y,label=cluster),colour='black',size = 10)
		print(lab)
		print(p)
		dev.off()
		print("已完成全部cluster的高亮图绘制。")
	}
	return (immune.combined_umap)
}

marker_gene <- function(immune.combined,min.pct = 0.1, logfc.threshold = 0.25,outdir=getwd(),pref='10x',test.use= "wilcox"){
	DefaultAssay(immune.combined) <- "RNA"
	cluster_num<-sort(unique(Idents(immune.combined)))
	print("Starting FindMarkers")
	all.markers <- FindAllMarkers(immune.combined, only.pos = FALSE, min.pct =min.pct, logfc.threshold = logfc.threshold,test.use=test.use)
	head(all.markers, n = 2)
	marker_result <- subset(all.markers, select=c(7,1,2,3,4,5,6))
	names(marker_result)[names(marker_result)=='gene'] <- 'gene_name'
	#write.csv(all.markers,paste(outdir,'3_marker',paste(pref,'all.markers.csv',sep='_'),sep='/'),quote=F)
	write.csv(marker_result,paste(outdir,'3_marker',paste(pref,'all.markers.csv',sep='_'),sep='/'),quote=F,row.names=F)
	print("Finished FindMarkers")
	return (all.markers)
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

#组合1	N	N2THY/N1THY/N3THY
#组合2	H	H2THY/H3THY/H1THY
#ini.list$sample$sample1  unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T))
sample_name<-unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T))
immune.combined<-readRDS(indir)
mkdirs(outdir,'1_QC')
setwd(paste(outdir,'1_QC',sep='/'))

immune.combined <- NormalizeData(object = immune.combined, normalization.method = ini.list$Para$object_list_normalization.method, scale.factor = as.numeric(ini.list$Para$object_list_scale.factor), verbose = FALSE)
immune.combined <- FindVariableFeatures(object = immune.combined, selection.method = ini.list$Para$object_list_findvariablefeatures_method, nfeatures = as.numeric(ini.list$Para$object_list_nfeatures_findvariablefeatures))
saveRDS(immune.combined, file = paste(prefix,'nor.rds',sep='_'))
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = as.numeric(ini.list$Para$integration_pca_dims), verbose = FALSE)
qc_pca_plot(immune.combined,outdir=paste(outdir,'1_QC',sep='/'),pref=prefix,w_h=c(as.numeric(unlist(strsplit(ini.list$Para$qc_pca_plot_w_h,split = ",",fixed=T))[1]),as.numeric(unlist(strsplit(ini.list$Para$qc_pca_plot_w_h,split = ",",fixed=T))[2])))

#聚类
mkdirs(outdir,'2_clusters')
setwd(paste(outdir,'2_clusters',sep='/'))
sample1<-unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T)) #样本名称
sample2<-unlist(strsplit(ini.list$sample$sample2,split = "/",fixed=T)) #分组名称
immune.combined$Group <- immune.combined$orig.ident
immune.combined$stim <- immune.combined$orig.ident
for (n in 1:length(sample1)){
	immune.combined$Group<-gsub(paste('\\b',sample1[n],'\\b',sep=''),sample2[n],immune.combined$Group)
}
immune.combined<-reduction(immune.combined,dims_num=as.numeric(ini.list$Para$reduction_dims_num),resolution=as.numeric(ini.list$Para$reduction_resolution),outdir=paste(outdir,'2_clusters',sep='/'),pref=prefix)

#Marker基因分析
DefaultAssay(immune.combined) <- "RNA"
mkdirs(outdir,'3_marker')
setwd(paste(outdir,'3_marker',sep='/'))
print(Sys.time())
saveRDS(immune.combined, file = paste(prefix,'cluster.rds',sep='_'))
all.markers<-marker_gene(immune.combined,min.pct =as.numeric(ini.list$Para$marker_gene_min.pct), logfc.threshold = as.numeric(ini.list$Para$marker_gene_logfc.threshold),outdir=outdir,test.use=ini.list$Para$marker_gene_test.use,pref=prefix)
immune.combined@misc$markers <- all.markers
marker_plot(immune.combined,all.markers,outdir=paste(outdir,'3_marker',sep='/'),pref=prefix)
print(Sys.time())

setwd(paste(outdir,'2_clusters',sep='/'))
saveRDS(immune.combined, file = paste(prefix,'immune_combined.rds',sep='_'))
