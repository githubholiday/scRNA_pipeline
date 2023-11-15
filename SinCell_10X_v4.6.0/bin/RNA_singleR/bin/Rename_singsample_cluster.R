#!/annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript
#名称：Rename_cluster.R
#作者：赵倩
#邮箱：mengchengyao@genome.cn
#时间：20210615
#版本：v0.0.3
#用途：利用seurat v4.0进行10x 进行单细胞类群注释上
###说明：
#程序开发环境/annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript，需要指定R，指定包的路径,使用为seurat4.0版本以上
#===========================================================
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'rds',	'r',	1,	"character",
	'anno',	'a',	1,	"character",
	'outdir',	'o',	1,	"character"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	========================================================================================================================================
	rds:the output of seurat analysis,rdsfile.
	========================================================================================================================================
	prefix:the output prefix of files,such as pictures and excel.
	========================================================================================================================================
	anno:anno for the cluster from file and scibet's result.
	========================================================================================================================================
	outdir:outdir  of outputs,we will setwd(opt$outdir)
	Usage example:
	Rscript this.r -r sample.rds -a annolist.xls -o outdir -p alp
	Options:
	--help		h	NULL		get this help
	--rds	r	character	indir for rds file[forced]
	--anno	a	character	anno file for cluster[forced]
	--outdir	o	character	The	resurt of out dir for analysis [forced]
	--prefix	p	character	the prefix for outputfiles [forced]
	\n")
	q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$rds) )	{ cat("Please input the rdsfile ...\n\n") ; print_usage(para)}
if ( is.null(opt$anno) )	{ cat("Please input the anno_list ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )	{ cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para) }
##这个分析用最新的seurat包进行分析
require(Seurat)
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(cowplot)

mkdirs <- function(outdir,fp) {
	if(!file.exists(file.path(outdir,fp))) {
#		mkdirs(dirname(fp))
		dir.create(file.path(outdir,fp))
	}else{
			print(paste(fp,"Dir already exists!",sep="     "))
			unlink(file.path(outdir,fp), recursive=TRUE)
			dir.create(file.path(outdir,fp))
		}
}

sub_cells <- readRDS(file =opt$rds)
label_names <- read.csv(opt$anno,header=T,sep='\t')
outdir <- opt$outdir
pref <- opt$prefix
seurat_data <-sub_cells
labers=label_names[match(as.numeric(as.character(seurat_data@active.ident)),label_names[,1]),2]
seurat_data$lables <-labers

p_tmp <- theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
		legend.title = element_blank(),legend.text =  element_text(color="black",size=18),
		axis.text.x = element_text(color="black",size=18),
		axis.text.y = element_text(color="black",size=18),
		axis.title.x = element_text(face="plain", color="black",size=18),
		axis.title.y = element_text(face="plain", color="black",size=18))
print("细胞聚类注释后umap图：")
pdf(paste(outdir,paste(pref,"umap_cluster_anno.pdf",sep='_'),sep='/'),w=24,h=8)
p1<-DimPlot(seurat_data, reduction = "umap",pt.size = 1)+p_tmp
p2<-DimPlot(seurat_data, reduction = "umap", group.by = "lables", label=T, label.size=5, pt.size=1, repel = T) + p_tmp
p3<-plot_grid(p1, p2)
print(p3)
dev.off()
print("细胞聚类注释后tsne图：")
seurat_tsne<-RunTSNE(seurat_data,reduction = "pca",dims = 1:20)
pdf(paste(outdir,paste(pref,"tsne_cluster_anno.pdf",sep='_'),sep='/'),w=24,h=8)
p1 <- DimPlot(seurat_tsne, reduction = "tsne", pt.size = 1)+p_tmp
p2 <- DimPlot(seurat_tsne, reduction = "tsne", group.by = "lables", label=T, label.size=5, pt.size=1, repel = T) + p_tmp
p3<-plot_grid(p1, p2)
print(p3)
dev.off()



