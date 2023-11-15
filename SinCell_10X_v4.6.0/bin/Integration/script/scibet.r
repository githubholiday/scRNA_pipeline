#!/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/Rscript
#名称：scibet.r
#作者：姚盟成
#邮箱：mengchengyao@genome.cn
#时间：2020-04-16
#版本：v0.0.1
#用途：利用scibet包对单细胞转录组数据进行细胞亚群鉴定，默认的数据模型，都是细胞较大的类型，如果细胞亚群分的特别的细的话，可能不适合。
###说明：
#程序开发环境/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/Rscript，需要指定R，指定包的路径,使用为seurat3.0版本以上
#===========================================================
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'model',	'm',	2,	"character",
	'rds',	'r',	1,	"character",
	'exp',	'e',	2,	"character",
	'species',	's',	2,	"character",
	'outdir',	'o',	1,	"character"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	========================================================================================================================================
	model文件:
	这个文件是从scibet包中官方下载，目前官方提供了100中数据集，可以根据自己感兴趣的数据集进行细胞亚群鉴定
	========================================================================================================================================
	prefix:the output prefix of files,such as pictures and excel.
	========================================================================================================================================
	rds 主要是seurat分析的结果rds文件，这里只支持seurat3.0以上版本的结果
	========================================================================================================================================
	outdir:outdir  of outputs,we will setwd(opt$outdir)
	Usage example:
	Rscript this.r -i1 BM-aLP_all_UMI.csv -i2 FL-alp_all_UMI.csv -o outdir -s1 BM-alp -s2 FL-alp -p alp
	Options:
	--help		h	NULL		get this help
	--model	m	character	model for celltype identity[option]
	--rds	r	character	result for analysis by seurat3.0 [forced]
	--exp	e	character	TPM file for single cell RNA,must be a csv file,and last colum should be cell label  [option]
		for example:
			gene1	gene2	gene3	.......	label
		cell-1	0	1	1.5	.....	T cell
		cell-2	1.3	1.98	1.5	.....	B cell
		......
		cell-n	2.0	3.28	4.8	....	B cell
	--species	s	character	species for celltype indenty[mouse or human]  [option]
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
if ( is.null(opt$exp) && is.null(opt$model) &&  is.null(opt$species) )	{ cat("Please give the train data file  ...\n\n") ; print_usage(para) }
if ( is.null(opt$model) )	{opt$model<-""}
if ( is.null(opt$exp) )	{opt$exp<-""}
if ( is.null(opt$species) )	{opt$species<-""}

outdir<-opt$outdir
pre<-opt$prefix
rds<-opt$rds
model<-opt$model
train<-opt$exp
species<-opt$species
##这个分析用最新的seurat包进行分析
suppressMessages(require(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
library(dplyr)
library(tibble)

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


####
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
###
#' Process scibet.core
#' @name pro.core
#' @usage pro.core(scibet.core)
#' @param scibet.core A SciBet core
#' @return A processed SciBet core
#' @export
pro.core <- function(scibet.core){
  cell.type <- unname(unlist(scibet.core[,1]))
  scibet.core <- as.data.frame(t(scibet.core[,-1]))
  colnames(scibet.core) <- cell.type
  return(as.matrix(scibet.core))
}

#' Heatmap of classification result.
#' @name Confusion_heatmap
#' @usage Confusion_heatmap(ori, prd)
#' @param ori A vector of the original labels for each cell in the test set.
#' @param prd A vector of the predicted labels for each cell in the test set..
#' @return A heatmap for the confusion matrix of the classification result.
#' @export
Confusion_heatmap <- function(ori, prd){
  tibble(
    ori = ori,
    prd = prd
  ) %>%
    dplyr::count(ori, prd) %>%
    tidyr::spread(key = prd, value = n) -> cross.validation.filt

  cross.validation.filt[is.na(cross.validation.filt)] = 0
  cross.validation.filt[,-1] <- round(cross.validation.filt[,-1]/rowSums(cross.validation.filt[,-1]),2)
  cross.validation.filt <- cross.validation.filt %>%
    tidyr::gather(key = 'prd', value = 'Prob', -ori)

  cross.validation.filt %>%
    ggplot(aes(ori,prd,fill = Prob)) +
    geom_tile() +
    theme(axis.title = element_text(size = 0)) +
    theme(axis.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 0)) +
    theme(legend.text = element_text(size = 10)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black",angle = 45, hjust = 1)) +
    scale_fill_viridis() -> p

  return(p)
}
#' Heatmap for the confusion matrix of the classification with the false postive control.
#' @name Confusion_heatmap_negctrl
#' @usage Confusion_heatmap_negctrl(res, cutoff = 0.4)
#' @param res Classification result.
#' @param cutoff The cutoff of confifence score C.
#' @return A heatmap for the confusion matrix of the classification result with the false postive control.
#' @export
Confusion_heatmap_negctrl <- function(res, cutoff = 0.4){
  res %>%
    dplyr::mutate(prd = ifelse(c_score < cutoff, 'unassigned', prd)) %>%
    dplyr::count(ori, prd) %>%
    tidyr::spread(key = prd, value = n) -> cla.res

  cla.res[is.na(cla.res)] = 0
  cla.res[,-1] <- round(cla.res[,-1]/rowSums(cla.res[,-1]),2)
  cla.res <- cla.res %>% tidyr::gather(key = 'prd', value = 'Prob', -ori)
  label <- cla.res$ori %>% unique()
  cla.res %>%
    ggplot(aes(prd, factor(ori, levels = c(label[-3],'Neg.cell')), fill = Prob)) +
    geom_tile(colour = 'white', lwd = 0.5) +
    theme(axis.title = element_text(size = 12)) +
    theme(axis.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 0)) +
    theme(legend.text = element_text(size = 12)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black", angle = 50, hjust = 1)) +
    scale_fill_material('blue')
}
#' Expression patterns of informative genes across cell types.
#' @name Marker_heatmap
#' @usage Marker_heatmap(expr, gene)
#' @param expr The expression dataframe. Rows should be cells, columns should be genes and last column should be "label".
#' @param gene A vector of informative genes.
#' @return A figure.
#' @export
Marker_heatmap_ymc <- function(expr, gene){
  expr <- expr[,c(gene,'label')]
  type_expr <- expr %>%
    tidyr::nest(-label) %>%
    dplyr::rename(expr = data) %>%
    dplyr::mutate(colmeans = purrr::map(
      .x = expr,
      .f = function(.x){colMeans(.x)}))

  type_expr$colmeans %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    t() %>%
    as.data.frame() %>%
    tibble::remove_rownames() -> type_mean_expr

  rownames(type_mean_expr) <- type_expr$label
  colnames(type_mean_expr) <- colnames(expr)[-ncol(expr)]

  sub_expr <- type_mean_expr
  sub_expr <- sub_expr %>%
    as_tibble() %>%
    dplyr::mutate_all(funs((. - mean(.))/sd(.))) %>%
    t()
  colnames(sub_expr) <- type_expr$label
  get_label <- function(num){
    v <- sub_expr[num,]
    colnames(sub_expr)[which(v == max(v))]
  }
  sub_expr <- sub_expr %>%
    tibble::as.tibble() %>%
    dplyr::mutate(group = purrr::map_chr(1:length(gene), get_label))
  sub_expr <- as.data.frame(sub_expr)
  rownames(sub_expr) <- gene
  sub_expr <- sub_expr %>%
    dplyr::mutate(gene = gene) %>%
    tidyr::gather(key = 'cell_type', value = 'zscore', -group, -gene) %>%
    dplyr::arrange(group, desc(zscore))
  if (sum(is.na(as.numeric(sub_expr$cell_type)))>0){
	  print("The celltype is not numeric,将安装字母顺序排序")
   }else{
	  print(unique(sub_expr$cell_type))
	  sub_expr$cell_type<-as.numeric(as.character(sub_expr$cell_type))
#	  sub_expr$cell_type<-factor(sub_expr$cell_type,levels=1:20)
	  print(sort(unique(sub_expr$cell_type), decreasing = T))
   }

  sub_expr %>%
    ggplot(aes(factor(gene, levels = unique(sub_expr$gene)),
               factor(cell_type, levels = sort(unique(sub_expr$cell_type), decreasing = T)))) +
    geom_point(aes(size = zscore, colour = zscore)) +
    theme(
      strip.text.x = element_blank(),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black", angle = -90, hjust = 0),
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      )
    ) +
    facet_grid(. ~ group, scales = "free", space = "free") +
    scale_colour_distiller(palette = "RdYlBu") +
    labs(
      x = '',
      y = ''
    ) -> p

  return(p)
}

celltype<-function(rdsfile,pre,model,train="",type="model"){
###获取处理的数据TPM值矩阵
	immune.combined <- readRDS(rdsfile)
	immune.combined@meta.data$seurat_clusters<-Idents(immune.combined)[rownames(immune.combined@meta.data)]
#	DefaultAssay(immune.combined) <- "RNA"
	query<-as.data.frame(t(as_matrix(immune.combined@assays$RNA@data)))
	print(query[1:10, 1:10])
	label<-immune.combined@meta.data[as.vector(rownames(query)),]$seurat_clusters
	label2<-immune.combined@meta.data[as.vector(rownames(query)),]$seurat_clusters
	query$label<-label
	etest_gene <- SelectGene(query, k = 100)
	print(etest_gene)
	pdf(paste(pre,"_marker_heatmap.pdf",sep=""),w=25,h=(1.5+0.3*length(unique(query$label)))) 
	p<-Marker_heatmap_ymc(query, etest_gene)
	print(p)
	dev.off()
	###
	if ( type =="exp"){
		train<-read.csv(train,sep=",",row.names=1)
		prd <- SciBet(train, query[,-ncol(query)])
		pdf(paste(pre,"_query_model_heatmap.pdf",sep=""),w=(8+5*length(unique(query$label))),h=(2+0.15*length(unique(prd))))
		p<-Confusion_heatmap(query$label, prd)
		print(p)
		dev.off()
		output<-tibble( Barcode=rownames(query),seurat_clusters= label2,prd_celltype = prd)#test_set$label, prd
		write.csv(output,paste(pre,"_celltype_predict",'.csv',sep=""),quote=F)
	}else {
		model <- readr::read_csv(model) 
		model <- pro.core(model)
		#query <- readr::read_rds("test.rds.gz")
		#query <- readr::read_rds("TEST.rds.gz")
		print(query[1:10, 1:10])
		ori_label <- query$label
		query <- query[,-ncol(query)]
		prd <- LoadModel(model)
		label <- prd(query)
		pdf(paste(pre,"_query_model_heatmap.pdf",sep=""),w=(8+5*length(unique(query$label))),h=(2+0.15*length(unique(label))))
		p<-Confusion_heatmap(ori_label,label)
		print(p)
		dev.off()
		output<-tibble( Barcode=rownames(query),seurat_clusters= label2,prd_celltype = label)#test_set$label, prd
		write.csv(output,paste(pre,"_celltype_predict",'.csv',sep=""),quote=F)
	}
}

#### 1 建输出文件夹
#mkdirs(outdir,'1_result')
#setwd(paste(outdir,'1_result',sep='/'))
setwd(outdir)

binpath=commandArgs(trailingOnly = F)
scriptdir=dirname(sub("--file=","",binpath[grep("--file",binpath)]))  ###获取脚本路径
#### 2 数据的准备
if (model=="" &&  species=="mouse" ){
	modelfile=paste(scriptdir,"/database/GSE109774_scibet_core.mouse20.csv",sep="/")
	#modelfile="/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/luhuiping/Pipline/SinCell_10X_Genomics/SinCell_10X/bin/CellIdent/database/GSE109774_scibet_core.mouse20.csv"
}else if ( model=="" &&  species=="human" ) {
	#modelfile="/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/luhuiping/Pipline/SinCell_10X_Genomics/SinCell_10X/bin/CellIdent/database/major_human_cell_types.csv"
	modelfile=paste(scriptdir,"/database/major_human_cell_types.csv",sep="/")
}else {
	print("靓仔或者美女，你输入的不是人或者小鼠物种，你选择了自己的道路，请注意！")
	modelfile=model
}
if ( train==""){
	type="model"
}else {
	type="exp"
}
print(paste("开始scibet分析，时间为：",Sys.time(),sep=''))
celltype(rds,pre,modelfile,train=train,type=type)
print(paste("完成scibet分析，时间为：",Sys.time(),sep=''))

