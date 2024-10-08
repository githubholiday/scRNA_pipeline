#名称：marker_celltype_multi.R
#作者：mengli
#邮箱：mengli@genome.cn
#时间：20231130
#版本：v0.0.3
#用途：对所有样本合并后的rds,基于scibet细胞注释后的结果，对不同分组的细胞进行统计
#===========================================================
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'inrds',	'i',	1,	"character",
	'config',	'c',	1,	"character",
	'outdir',	'o',	1,	"character",
	'marker',	'm',	'2',	"character",
	'celltype',	't',	'2',	"character"
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
	--marker	m	character	第4列为marker基因的表格
	--celltype	t	character	第2列为cluster，第3列为细胞类型的表格
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
mkdirs <- function(outdir, fp){
	if(!file.exists(file.path(outdir,fp))) {
		dir.create(file.path(outdir,fp))
	}else{
		print(paste(fp,"Dir already exists!",sep="     "))
		#unlink(file.path(outdir,fp), recursive=TRUE)
		#dir.create(file.path(outdir,fp))
		}
}

marker_plot<-function(immune.combined,marker_gene,outdir=getwd(),pref='10x',heatmap_n=10){
#all.genes <- rownames(immune.combined)
#immune.combined <- ScaleData(immune.combined, features = all.genes)
small<-immune.combined
Idents(small) <- small$Cell_type
print(levels(factor(Idents(small))))
all.genes <- rownames(small)
#small<-ScaleData(small, features = all.genes)
#gene<-as.vector(unique(marker_gene))
exprs <- data.frame(FetchData(object = small, vars = as.vector(marker_gene)))
exprs$Barcod<-rownames(exprs)
ident<-data.frame()
#barcode与聚类信息提取
ident<-data.frame(Barcod=rownames(small@meta.data),orig.ident=small@meta.data$Cell_type,samples=small@meta.data$stim)
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
head(a)
mkdirs(outdir,'marker_plot')
for (g in unique(a$Gene)){
	b<-a[a$Gene==g,]
	
	#marker基因的高量图
	pdf(paste(outdir,'marker_plot',paste(g,'FeaturePlot.pdf',sep='_'),sep='/'),w=12,h=8)
	p<-FeaturePlot(immune.combined, features = as.vector(g), min.cutoff = "q9",cols = c("lightgrey", "red"))
	print(p)
	dev.off()
	
}
	nrow=ceiling(length(unique(a$Gene))/1)

	pdf(paste(outdir,paste(pref,'all','clusters.pdf',sep='_'),sep='/'),w=4*length(unique(Idents(small))),h=2*length(marker_gene))
	p2<-VlnPlot(small, features = marker_gene,split.by = "Group", group.by = "Cell_type", combine = FALSE,pt.size = -1,slot="data",assay="RNA",stack = T,flip =T) + theme(text = element_text(size = 30))+theme(axis.text.y = element_text(color="black",size=26,angle=0),axis.text.x = element_text(color="black",size=26,angle=30))+xlab('')
	print(p2)
	dev.off()

	print('dotplot画图：')
	pdf(paste(outdir,paste(pref,'top_dotplot.pdf',sep='_'),sep='/'),w=0.5*length(unique(Idents(small)))+4,h=0.5*length(unique(marker_gene)))
	#p<-DotPlot(object =seurat_data1, features=unique(gene_exist),cols = c("lightgrey", "red"))+RotatedAxis()
	p<-DotPlot(object =small, features=unique(marker_gene),cols = c("blue", "red"))+RotatedAxis()+coord_flip()+
	labs(x="Marker",y="Celltype")+
	theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18)) +
	theme(axis.text.x = element_text(angle = 45))
	print(p)
	dev.off()
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

Rename_cluster<-function(seurat_data,label_names,outdir=getwd(),pref='10x'){
	labers=label_names[match(as.numeric(as.character(seurat_data@active.ident)),label_names[,1]),2]
	seurat_data$Cell_type <-labers
	p_tmp <- theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
		legend.title = element_blank(),legend.text =  element_text(color="black",size=12),
		axis.text.x = element_text(color="black",size=12),
		axis.text.y = element_text(color="black",size=12),
		axis.title.x = element_text(face="plain", color="black",size=12),
		axis.title.y = element_text(face="plain", color="black",size=12))
	print("细胞聚类注释后umap图：")
	pdf(paste(outdir,paste(pref,"umap_cluster_anno.pdf",sep='_'),sep='/'),w=24,h=8)
	p1<-DimPlot(seurat_data, reduction = "umap",pt.size = 1)+p_tmp
	p2<-DimPlot(seurat_data, reduction = "umap", group.by = "Cell_type", label=T, label.size=5, pt.size=1, repel = T) + p_tmp
	p3<-plot_grid(p1, p2)
	print(p3)
	dev.off()
	pdf(paste(outdir,paste(pref,"umap_group_anno.pdf",sep='_'),sep='/'),w=12,h=8)
	p4<-DimPlot(seurat_data, reduction = "umap", group.by = "Cell_type", label=F, pt.size=1, split.by = "Group", ncol=2) + p_tmp
	print(p4)
	dev.off()
	pdf(paste(outdir,paste(pref,"umap_sample_anno.pdf",sep='_'),sep='/'),w=12,h=8)
	p5<-DimPlot(seurat_data, reduction = "umap", group.by = "Cell_type", label=F, pt.size=1, split.by = "orig.ident", ncol=2) + p_tmp
	print(p5)
	dev.off()

	
	#各样本细胞类型统计
	a<-data.frame(CellType=seurat_data@meta.data[,'Cell_type'],Samples=seurat_data@meta.data[,'orig.ident'])
	a$CellType<-factor(a$CellType,levels=sort(unique(a$CellType)))
	if (length(unique(a$Samples)) > 1) {	
	pdf(paste(outdir,paste(pref,"cellType.sample_stats.pdf",sep='_'),sep='/'),w=length(unique(a$CellType))*1+4,h=max(nchar(as.character(unique(a$CellType))))*0.1+6)
	p<-ggplot(a, aes(CellType)) + geom_bar(aes(fill=Samples), position='fill',width=0.6)+labs(x=" ", y = "",fill= "Samples")+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.background = element_rect(colour = NA),plot.title=element_text(size = 25),axis.text.x=element_text(size=20,angle=60,hjust=1),axis.text.y=element_text(size=20),axis.title.x=element_text(size = 25),axis.title.y=element_text(size = 25),legend.text =element_text(size = 25),legend.title =element_text(size = 25))
	print(p)
	dev.off()
	}
	a1 <- get_table(as.data.frame.array(table(a)))
	a2 <- data.frame(celltype=rownames(a1), a1)
	write.table(a2,paste(pref,'sample.celltype.xls',sep='_'),quote=F,sep="\t", row.names=F)
	
	#各比较组细胞类型统计
	b<-data.frame(CellType=seurat_data@meta.data[,'Cell_type'],Group=seurat_data@meta.data[,'Group'])
	b$CellType<-factor(b$CellType,levels=sort(unique(b$CellType)))
	if (length(unique(b$Group)) > 1) {
	pdf(paste(outdir,paste(pref,"cellType.group_stats.pdf",sep='_'),sep='/'),w=length(unique(a$CellType))*1+4,h=max(nchar(as.character(unique(a$CellType))))*0.1+6)
	p<-ggplot(b, aes(CellType)) + geom_bar(aes(fill=Group), position='fill',width=0.6)+labs(x=" ", y = "",fill= "Group")+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.background = element_rect(colour = NA),plot.title=element_text(size = 25),axis.text.x=element_text(size=20,angle=60,hjust=1),axis.text.y=element_text(size=20),axis.title.x=element_text(size = 25),axis.title.y=element_text(size = 25),legend.text =element_text(size = 25),legend.title =element_text(size = 25))
	print(p)
	dev.off()
	}
	b1 <- get_table(as.data.frame.array(table(b)))
	b2 <- data.frame(celltype=rownames(b1), b1)
	write.table(b2,paste(pref,'group.celltype.xls',sep='_'),quote=F,sep="\t",row.names=F)
	print("细胞聚类注释后统计完成")	
	return(seurat_data)
}


#主流程
#2_clusters  3_marker  4_conserved_markers  5_diff_gene_condition
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
sample1<-unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T)) #样本名称
sample2<-unlist(strsplit(ini.list$sample$sample2,split = "/",fixed=T)) #分组名称
immune.combined$Group <- immune.combined$orig.ident
immune.combined$stim <- immune.combined$orig.ident
all_sample <- unique(immune.combined$orig.ident)
for (n in 1:length(sample1)){
	if (sample1[n] %in% all_sample){
		immune.combined$Group<-gsub(paste('\\b',sample1[n],'\\b',sep=''),sample2[n],immune.combined$Group)
	} else {
		stop(paste(sample1[n], 'not in the rds!!', sep=' '))
	}
}
#按照样本名称提取出对应的细胞-limeng

immune.combined= immune.combined[,immune.combined@meta.data$orig.ident %in% c(sample1)]

#Marker基因分析
DefaultAssay(immune.combined) <- "RNA"
print(Sys.time())


mkdirs(outdir,'CellIdent')
setwd(paste(outdir,'CellIdent',sep='/'))
if ( is.null(opt$celltype) ) {
	print('not celltype')
} else {
	a <- read.table(opt$celltype, header = T, sep = '\t')
	#a <- a[a[,1]==prefix,]
	label_names <- a[,c(2,3)]
	immune.combined <- Rename_cluster(immune.combined,label_names,outdir=getwd(),pref=prefix)
}

## 获取原始矩阵的assay
ori_assay <- ""
if("RNA" %in% Assays(object = immune.combined)){
	print("原始矩阵是 RNA Assay")
	ori_assay <- "RNA"
} else if("SCT" %in% Assays(object = immune.combined)) {
	print("原始矩阵是 SCT Assay")
	ori_assay <- "SCT"
} else{
	print("请确认原始矩阵的assay 类型")
}
saveRDS(immune.combined, file = paste(prefix,'immune_combined_celltype.rds',sep='_'))

if ( is.null(opt$marker) ) {
	print('not marker')
} else {
	mkdirs(outdir,'select_marker')
	setwd(paste(outdir,'select_marker',sep='/'))
	a <- read.table(opt$marker, header = T, sep = '\t')
	a <- a[a[,1]==prefix,]
	#Marker_gene <- unique(unlist(strsplit(a[,4], ',')))
	Marker_gene <- unique(unlist(strsplit(a[,4], '[,，]')))
	#select_gene <- marker_realname(Marker_gene)
	genes <- rownames(GetAssayData(object = immune.combined[[ori_assay]], slot = "counts"))
	gene_exist <- c()
	gene_notexist <- c()
	for (gene in unique(Marker_gene)){
		if(length(grep(paste("^",gene,"$",sep=""), genes, ignore.case = T)) > 0){ #区分大小写
			Mark_real = genes[grep(paste("^",gene,"$",sep=""), genes, ignore.case = T)[1]]
			gene_exist<-c(gene_exist,Mark_real)
		}else{
			gene_notexist<-c(gene_notexist,gene)
		}
	}
	print(gene_exist)
	print('以下基因无法识别：')
	print(gene_notexist)
	marker_plot(immune.combined,gene_exist,outdir=getwd(),pref=prefix)
}



