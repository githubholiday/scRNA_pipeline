#!/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/envs/cellassign/bin/Rscript
#名称：singleR1.0.6.r
#作者：姚盟成
#邮箱：mengchengyao@genome.cn
#时间：2020-04-22
#版本：v0.0.2
#用途：利用新版的singleR对细胞亚群进行自动鉴定，目前仅仅支持人和小鼠两个物种，7个数据集。
###说明：
#程序开发环境/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/envs/cellassign/bin/Rscript，需要指定R，指定包的路径,使用为seurat3.0版本以上
#主要支持多个样品的，列名必须有stim
#===========================================================
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'marker',	'm',	2,	"character",
	'rds',	'r',	1,	"character",
	'db',	'd',	1,	"character",
	'species',	's',	1,	"character",
	'config',	'c',	2,	"character",
	'outdir',	'o',	1,	"character"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	========================================================================================================================================
	marker:
	是否进行marker基因差异分析，默认不进行差异分析，这一步比较耗时
	========================================================================================================================================
	prefix:the output prefix of files,such as pictures and excel.
	========================================================================================================================================
	rds 主要是seurat分析的结果rds文件，这里只支持seurat3.0以上版本的结果
	========================================================================================================================================
	outdir:outdir  of outputs,we will setwd(opt$outdir)
	Usage example:
	Rscript this.r -r rds -s mouse or human  -o outdir -d ImmGenData -p alp [-m marker ]
	Options:
	--help		h	NULL		get this help
	--marker	m	character	marker analysis or not [option]
	--rds	r	character	result for analysis by seurat3.0 [forced]
	--db	e	character	ref_data for analysis, only 7 data sets, chose one of follows  [forced]
		for example:
		[1] ImmGenData                       MouseRNAseqData                 
		[3] HumanPrimaryCellAtlasData        BlueprintEncodeData             
		[5] DatabaseImmuneCellExpressionData NovershternHematopoieticData    
		[7] MonacoImmuneData
	--species	s	character	species for celltype indenty[mouse or human]  [forced]
	--outdir	o	character	The	resurt of out dir for analysis [forced]
	--prefix	p	character	the prefix for outputfiles [forced]
	\n")
	q(status=1)
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
#if ( is.null(opt$model) )	{ cat("Please input the model data  ...\n\n") ; print_usage(para)}
if ( is.null(opt$rds) )	{ cat("Please input the rds data ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para) }
if ( is.null(opt$prefix) )	{ cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para) }
#if ( is.null(opt$db) && is.null(opt$model) &&  is.null(opt$species) )	{ cat("Please give the train data file  ...\n\n") ; print_usage(para) }
if ( is.null(opt$marker) )	{opt$marker<-""}
if ( is.null(opt$db) )	{cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para)}
if ( is.null(opt$config) )	{opt$config<-"/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/B_TET-024/std/result/changjiyang/yaomengcheng/singleR/v1.0.6/database/config.ini"}
if ( is.null(opt$species) )	{cat("Please give the prefix for outputfiles ...\n\n") ; print_usage(para)}

suppressMessages(library(scRNAseq))
suppressMessages(library(SingleR))
suppressMessages(library(pheatmap))
suppressMessages(library(configr))
suppressMessages(library(cowplot))
suppressMessages(library(Seurat))
suppressMessages(library(scater))
suppressMessages(library(dplyr))

outdir<-opt$outdir
pre<-opt$prefix
rds<-opt$rds
marker<-opt$model
db<-opt$db
species<-opt$species

as_matrix <- function(mat){
 
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

#seurat对象处理以及细胞亚群预测
RG<-function(rds,ref_data,prefix,w_h=c(12,8)){
#	immune.combined <- readRDS(rds)
	querry <- as_matrix(rds@assays$RNA@data)
	print(paste("开始分析，时间为：",Sys.time(),sep=''))
	pred.hesc <- SingleR(test = querry, ref =ref_data, labels = ref_data$label.fine) ###这里选择的是用label.fine，而没有用label.main，主要是考虑结果的精细化
	#pred.hesc_label.fine <- SingleR(test = data, ref =ref_data, labels = ref_data$label.fine)
	print(paste("完成分析，时间为：",Sys.time(),sep=''))
	print(table(pred.hesc$labels))

	####画热图列注释信息
	annotation_col<-data.frame(clusters=immune.combined@meta.data$seurat_clusters,stim=immune.combined@meta.data$stim)
	rownames(annotation_col)<-rownames(immune.combined@meta.data)
	order_cells<-annotation_col %>% dplyr::mutate(.,barcod=rownames(.)) %>% dplyr::arrange(.,clusters,stim)
	pred.hesc$clusters<-annotation_col[as.vector(rownames(pred.hesc)),]$clusters
	pred.hesc$stim<-annotation_col[as.vector(rownames(pred.hesc)),]$stim
	pred.hesc<-pred.hesc[as.vector(order_cells$barcod),]
	pdf(paste(prefix,'stim_celltype.pdf',sep="_"),w=w_h[1],h=w_h[2]) 
	p<-plotScoreHeatmap(pred.hesc,annotation_col=annotation_col,cells.order=order(pred.hesc$clusters))
	print(p)
	dev.off()
	write.table(dplyr :: bind_cols(Barcode=rownames(pred.hesc),as.data.frame(pred.hesc)),paste(prefix,'singleR_celltype.xls',sep="_"),sep="\t",row.names=F) 
	saveRDS(pred.hesc,paste(prefix,'singleR.rds',sep="_"))
	print(paste("完成singleR细胞亚群自动鉴定分析，时间为：",Sys.time(),sep=''))
	return (pred.hesc)
}

MK<-function(rds,pred.hesc,prefix){
	to.remove <- pruneScores(pred.hesc)
	print("remove的细胞具体情况如下：")
	print(summary(to.remove))
	##Based on the deltas across cells
	pdf(paste(prefix,'ScoreDistribution.pdf',sep="_"),w=12,h=3+2.5*length(unique(pred.hesc$labels))) 
	p<-plotScoreDistribution(pred.hesc, show = "delta.med", ncol = 5, show.nmads = 3)
	print(p)
	dev.off()
	#By default, SingleR() will report pruned labels in the pruned.labels field where low-quality assignments are replaced with NA
	new.pruned <- pred.hesc$labels
	new.pruned[pruneScores(pred.hesc, nmads=5)] <- NA
	print("NA的细胞的情况如下：")
	print(table(new.pruned, useNA="always"))
	y<-do.call(SummarizedExperiment, c(list(assays=as_matrix(rds@assays$RNA@data))))
	y$labels <- pred.hesc$labels
	names(assays(y))<-'logcounts'
	all.markers <- metadata(pred.hesc)$de.genes
	y$labels <- pred.hesc$labels
	for (lab in unique(pred.hesc$labels)) {
		print(lab)
		pdf(paste(prefix,gsub(" ",'_',lab),'celltype_marker.pdf',sep="_"),w=12,h=8)
		p<-plotHeatmap(y, order_columns_by=list(I(pred.hesc$labels)),features=unique(unlist(all.markers[[lab]])),colour_columns_by='labels',show_colnames =FALSE,) 
		print(p)
		dev.off()
	}
}
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
print(paste("读取配置文件，时间为：",Sys.time(),sep=''))
#logcounts(sce) <- log2(normcounts(sce)+1)
ini<-'/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/B_TET-024/std/result/changjiyang/yaomengcheng/singleR/v1.0.6/database/config.ini'
ini<- read.config(file = opt$config)
#ini<- read.config(file = ini)
#all_db<-c('HumanPrimaryCellAtlasData','BlueprintEncodeData',)
if ( length(intersect(species,c('human','mouse')))==0 ){ print(paste("目前只支持人和小鼠两个物种，请注意输入的物种，当前输入的物种为：",species,sep="  "))}
if ( length(intersect(db,c(names(ini$mouse),names(ini$human))))==0 ){ 
	print(paste("目前只支持人和小鼠7个数据集，请注意输入的参数数据集名称，当前输入的数据集为：",db,sep="    "))
	print(paste("目前能支持的参考数据集为：",c(names(ini$mouse),names(ini$human)),sep="  "))
	q()
	}
print(paste("读取rds文件，时间为：",Sys.time(),sep=''))
rds <- readRDS(rds)
print(paste("读取ref_data文件，时间为：",Sys.time(),sep=''))
ref_data<-readRDS(paste(ini$dir$ref_data,ini[[species]][[db]],sep="/")) #file.path
print(paste(ini$dir$ref_data,ini[[species]][[db]],sep="/"))
###开始分析
mkdirs(outdir,'result/')
setwd(paste(outdir,'result/',sep='/'))
pred.hesc<-RG(rds,ref_data,prefix)
if (marker==""){
	print(paste("完成singleR细胞亚群自动鉴定分析,不进行marker基因分析，时间为：",Sys.time(),sep=''))
	q()
}else{
	print(paste("开始进行marker基因分析，时间为：",Sys.time(),sep=''))
	MK(rds,pred.hesc,prefix)
	print(paste("完成marker基因分析，时间为：",Sys.time(),sep=''))
}
