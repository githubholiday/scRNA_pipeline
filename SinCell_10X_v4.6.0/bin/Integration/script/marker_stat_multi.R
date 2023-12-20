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
	p1<-DimPlot(seurat_data, reduction = "umap",pt.size = 1,label=T,repel = T)+p_tmp
	p2<-DimPlot(seurat_data, reduction = "umap", group.by = "Cell_type", label=T, label.size=5, pt.size=1, repel = T) + p_tmp + ggtitle("")
	p3<-plot_grid(p1, p2)
	print(p3)
	dev.off()
	pdf(paste(outdir,paste(pref,"umap_sample_anno.pdf",sep='_'),sep='/'),w=12,h=8)
	p5<-DimPlot(seurat_data, reduction = "umap", group.by = "Cell_type", label=F, pt.size=1, split.by = "orig.ident", ncol=2) + p_tmp
	print(p5)
	dev.off()
	#print("细胞聚类注释后tsne图：")
	#seurat_tsne<-RunTSNE(seurat_data,reduction = "pca",dims = 1:20)
	#seurat_tsne <- seurat_data
	#pdf(paste(outdir,paste(pref,"tsne_cluster_anno.pdf",sep='_'),sep='/'),w=24,h=8)
	#p1 <- DimPlot(seurat_tsne, reduction = "tsne", pt.size = 1)+p_tmp
	#p2 <- DimPlot(seurat_tsne, reduction = "tsne", group.by = "Cell_type", label=T, label.size=5, pt.size=1) + p_tmp
	#p3<-plot_grid(p1, p2)
	#print(p3)
	#dev.off()
	
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
immune.combined$stim <- immune.combined$orig.ident
all_sample <- unique(immune.combined$orig.ident)
DefaultAssay(immune.combined) <- "RNA"
print(Sys.time())


mkdirs(outdir,'4_CellIdent')
setwd(paste(outdir,'4_CellIdent',sep='/'))
if ( is.null(opt$celltype) ) {
	print('not celltype')
} else {
	a <- read.table(opt$celltype, header = T, sep = '\t')
	#a <- a[a[,1]==prefix,]
	label_names <- a[,c(2,3)]
	immune.combined <- Rename_cluster(immune.combined,label_names,outdir=getwd(),pref=prefix)
}

saveRDS(immune.combined, file = paste(prefix,'immune_combined_celltype.rds',sep='_'))

