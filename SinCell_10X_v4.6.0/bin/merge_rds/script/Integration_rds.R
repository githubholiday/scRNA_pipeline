# 组合1	N	N2THY/N1THY/N3THY
# 组合2	H	H2THY/H3THY/H1THY
# SCH1THY	Homo sapiens	样品2	SCH
# 1)提供组间比较结果
# N、H：N为对照组，H为处理组
# N、H、SCH：N为对照组 处理组H、SCH

# 2）基因在各组的哪些细胞中是否有差异:NLRP3   AIM2   NLRP1A   NLRP1B   NLRC4   caspase-1   caspase-11    Gasdermin-D     IL18   IL-1β   ，想看看在哪些组织中有表达，在各组中是否有差异，尤其是想看看NLRP3
#!/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/Rscript
#名称：DE_10xGenomics.R
#作者：姚盟成
#邮箱：mengchengyao@genome.cn
#时间：201901011
#版本：v0.0.2
#用途：利用seurat v3.0进行10x 进行不同条件比较分析,需要输入配置文件，配置可以设置具有生物学重复的分析，指定分组。如果有
#两个以上样品，则需要分别取做差异分析。
###说明：
#程序开发环境/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/Rscript，需要指定R，指定包的路径,使用为seurat3.0版本以上
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
#library(harmony)
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
integration <- function(object_list,cca_dims=20,pca_dims=20,RunPCA_npcs=30){
  all_samples<-names(object_list)
  if(length(all_samples)>1){
  immune.anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:as.numeric(cca_dims))
  immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:as.numeric(pca_dims))
  DefaultAssay(immune.combined) <- "integrated"
  }else{
  print(paste('You just have only one sample: hahahahah',all_samples,sep=' '))
  immune.combined=object_list[[1]]
  }
  
  return (immune.combined)
}

#1_QC  2_clusters  3_marker  4_conserved_markers  5_diff_gene  
prefix<-opt$prefix
outdir<-opt$outdir
indir<-opt$indir
ini<-opt$config
ini.list <- read.config(file = ini)
mkdirs(outdir,'1_QC')
setwd(paste(outdir,'1_QC',sep='/'))
#组合1	N	N2THY/N1THY/N3THY
#组合2	H	H2THY/H3THY/H1THY
#ini.list$sample$sample1  unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T))
sample_name<-unlist(strsplit(ini.list$sample$sample1,split = "/",fixed=T))
object_list<-sample2SeuratObject_list(indir,paste(outdir,'1_QC',sep='/'), sample_name)
#immune.combined<-integration(object_list,cca_dims=as.numeric(ini.list$Para$integration_cca_dims),pca_dims=as.numeric(ini.list$Para$integration_pca_dims),RunPCA_npcs=as.numeric(ini.list$Para$integration_runpca_npcs))
immune.combined=merge(object_list[[1]],object_list[2:length(object_list)])
#immune.combined <- NormalizeData(immune.combined) %>% FindVariableFeatures()
#immune.combined <- RunHarmony(immune.combined, group.by.vars = "stim")
saveRDS(immune.combined, file = paste(prefix,'.rds',sep=''))


