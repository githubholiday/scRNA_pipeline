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

DE_gene<-function(immune.combined,outdir,logfc.threshold = 0.25, test.use = "wilcox",min.pct = 0.1){
	tmp <- immune.combined
	tmp$celltype.stim <- paste(Idents(tmp), tmp$Group, sep = "_")
	tmp$celltype <- Idents(tmp)
	Idents(tmp) <- "celltype.stim"
	cluster_num<-length(sort(unique(Idents(immune.combined))))
	for (i in 1:cluster_num){
			dir.create(paste(outdir,"/",i,sep=""))
			for (m in 1:length(ini.list$cmp)){
					name <- paste(i,unlist(strsplit(ini.list$cmp[[m]],split = "/",fixed=T))[1],"VS",unlist(strsplit(ini.list$cmp[[m]],split = "/",fixed=T))[2],sep="_")
					diff_g.dir <- paste(outdir,"/",i,"/",name,sep="")
					dir.create(diff_g.dir,showWarnings = TRUE, recursive = TRUE,mode = "0777")
					###差异分析中，处理组_VS_对照组，因此需要ini文件中，顺序正确，正确的顺序为：处理组/对照组
		cmp_list = unlist(strsplit(ini.list$cmp[[m]],split = "/",fixed=T))
		case_cmp = cmp_list[1]
		control_cmp = cmp_list[2]
		cmp_name = paste(case_cmp,"_VS_",control_cmp,sep="")
		cmp_prefix = paste(outdir,"/",cmp_name,"/",name,sep="")

					tryCatch({b.interferon.response <- FindMarkers(tmp, ident.1 = paste(i,case_cmp,sep="_"), ident.2 = paste(i,control_cmp,sep="_"),print.bar = FALSE,logfc.threshold =logfc.threshold ,test.use =test.use ,min.pct = min.pct)

		gene_tsv = paste(cmp_prefix,"_diff_gene.csv",sep="")

					write.csv(b.interferon.response,gene_tsv,quote=F)
					b.interferon.response$gene<-rownames(b.interferon.response)
					b.interferon.response<-b.interferon.response[order(-abs(as.numeric(b.interferon.response$avg_log2FC))),]
					gene<-c()
		//
					if(length(b.interferon.response$gene) >=10){
							gene<-b.interferon.response$gene[1:10]
							}else{
							gene<-b.interferon.response$gene
					}//
					tmp2<-subset(tmp,cells=rownames(tmp@meta.data[tmp@meta.data$stim==case_cmp|tmp@meta.data$stim==control_cmp,])) #by yaomengcheng 20191115
					pdf(paste(cmp_prefix,"_VlnPlot.pdf",sep=""),w=0.8*cluster_num,h=3*length(gene))
					plots <- VlnPlot(tmp2, features = as.vector(gene), split.by = "stim", group.by = "celltype",pt.size = 0, combine = FALSE)
					p<-CombinePlots(plots = plots, ncol = 1)
					print(p)
					dev.off()
					pdf(paste(cmp_prefix,"_FeaturePlot.pdf",sep=""),w=12,h=3*length(gene))
					p<-FeaturePlot(tmp, features = as.vector(gene), split.by = "stim", max.cutoff = 3,cols = c("grey", "red"))
					print(p)
					dev.off()
					a<-FetchData(tmp2,as.vector(gene))
					a$Barcod<-rownames(a)
					b<-tmp@meta.data[,c('orig.ident','stim','celltype')]
					b$Barcod<-rownames(tmp@meta.data)
					a_b<-merge(a,b,by='Barcod')
					signif_data<-data.frame(Expression=0,stim=0,celltype=0,orig.ident=0,Gene=0)
					for (g1 in gene){
							tmp_data<-a_b[,c(g1,'stim','celltype','orig.ident')]
							tmp_data$Gene<-rep(g1,nrow(tmp_data)) #MT1X stim celltype orig.ident Gene
							colnames(tmp_data)<-c('Expression','stim','celltype','orig.ident','Gene')
							signif_data<-rbind(signif_data,tmp_data)
					}
					signif_data<-signif_data[2:nrow(signif_data),]
					library(ggsignif)
					compaired <- list(cmp_list)
					pdf(paste(cmp_prefix,"_signif_Plot.pdf",sep=""),w=12,h=10)

					p<-ggplot(signif_data[signif_data$celltype==i,],aes(stim,Expression,fill=stim))+geom_boxplot(width=0.5)+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test,vjust=5)+facet_wrap(~Gene, nrow =5,ncol =4,scales = 'free',strip.position='bottom')+theme(panel.border=element_blank(),panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'),plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23),strip.text = element_text(size = 10),strip.background = element_rect(fill = NA, colour = NA),strip.placement = "inside")+labs(x='Gene', y= 'Expression')
					print(p)
					dev.off()
					print(paste(i,cmp_name,"diff_gene.csv",sep="_"))},error=function(e){print(paste(i,cmp_name,"diff_gene.csv","--聚类不在某个组中,出现错误",sep="_"))} )
					}
			}
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

DefaultAssay(immune.combined) <- "RNA"

#多样本差异分析
if (length(sample_name) < 2){
	print(paste("您只有一个样品，不能做差异分析，样品名为：", sample_name, sep=' '))
	q()
}
mkdirs(outdir,'5_diff_gene_condition')
setwd(paste(outdir,'5_diff_gene_condition',sep='/'))
DE_gene(immune.combined, paste(outdir,'5_diff_gene_condition',sep='/'), logfc.threshold = as.numeric(ini.list$Para$de_gene_logfc.threshold), test.use = ini.list$Para$de_gene_test.use, min.pct = as.numeric(ini.list$Para$de_gene_min.pct))
