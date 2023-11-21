library(getopt)
para<-matrix(c(
	'help',	    'h',	0,	"logical",
	'datadir',	'd',	1,	"character",
	'outdir',	'o',	1,	"character",
	'sample',	's',	1,	"character",
	'minumi',	'mu',	2,	"integer",
	'mincell',	'mc',	2,	"integer",
	'min',	    'min',	2,	"integer",
	'max',	    'max',	2,	"integer",
	'mtp',	    'mtp',	2,	"integer",
	'hbp',	    'hbp',	2,	"integer",
	'species',	'sp',	1,	"character",
	'mt',	    'mt',	2,	"character",
	'testmethod','tm',	2,	"character"
),byrow=TRUE,ncol=4)
print("1")
opt <- getopt(para,debug=TRUE)
print("4")
#=============
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	Usage example:
	Rscript securat_parameter.R -d hg19/ -p ANproject -o outdir/
	Options:
	--help		h	NULL		get this help
	--datadir	d	character	the cellRanger output datadir [forced]
	--sample	s	character	SampleName[forced]
	--outdir	o	character	output file dir [forced]
	--minumi	mu	character	Filter genes have minium umi count [default: 200]
	--mincell	mc	character	Filter genes expession based on minium cell numbers [default: 3]
	--min	min	character	Filter cells based on minium gene numbers [default: 200]
	--max	max	character	 Filter cells based on max gene numbers [default: 10000]
	--mtp	mtp	character	Filter cells based on mito percent [default:20]
	--hbp	hbp	character	Filter cells based on HB percent [default:100]
	--mtt	mtt	character	mitochondria Gene symbol [default : MT]
	--species	sp	character	species[human:9606 mouse:10090]
	--testmethod	tm	character	test method for findMakers [default: wilcox]
	\n")
	q(status=1)
}
print("3")
if ( is.null(opt$datadir) || is.null(opt$sample) || is.null(opt$outdir) || is.null(opt$species)){ print_usage(para) } 
if ( is.null(opt$mincell))	{ opt$mincell <- 3 }
if ( is.null(opt$minumi))	{ opt$minumi <- 200 }
if ( is.null(opt$min))	{ opt$min <- 200 }
if ( is.null(opt$max))	{ opt$max <- 10000 }
if ( is.null(opt$mtp))	{ opt$mtp <- 20 }
if ( is.null(opt$hbp))	{ opt$hbp <- 100}
if ( is.null(opt$mtt))	{ opt$mtt <- c("mt") }
if ( is.null(opt$testmethod))	{ opt$testmethod <- c("wilcox") }

#======================
library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)

print("2")
prefix <- paste(opt$outdir,opt$sample,sep='/')
## creat object
expression_matrix <- Read10X(data.dir = opt$datadir) ## 数据读入路径
#object <- CreateSeuratObject(counts = expression_matrix, project = opt$sample, min.cells = opt$mincell, min.features = opt$minumi) ### 设置最小细胞中表达的基因和最小基因表达数的细胞
object <- CreateSeuratObject(counts = expression_matrix, project = opt$sample, min.cells = opt$mincell)
print("Rawgene:")
RawGene <- nrow(object@meta.data)
print(RawGene)
#threshold <- unname(quantile(object@meta.data$nFeature_RNA,0.25)) #细胞中表达基因数目的下四分位数，probs参数可调
#print(threshold)
#threshold <-max(threshold,200)
#file <- paste(prefix,'filter_UMI.csv',sep='_')
#print(file)
species<<-opt$species
if (species == '9606'){ 
  HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # 人类血液常见红细胞基因
  HB_m <- match(HB.genes_total,rownames(object@assays$RNA))
  HB.genes <- rownames(object@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  object[["percent.HB"]]<-PercentageFeatureSet(object,features=HB.genes)
 }else if (species == '10090'){
   HB.genes_total <- c("Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y","Hba-x","Hba-a1","Hba-a2") # 小鼠血液常见红细胞基因
   HB_m <- match(HB.genes_total,rownames(object@assays$RNA))
   HB.genes <- rownames(object@assays$RNA)[HB_m]
   HB.genes <- HB.genes[!is.na(HB.genes)]
   object[["percent.HB"]]<-PercentageFeatureSet(object,features=HB.genes)
 }else{
   object@meta.data$percent.HB<-rep(0,length(rownames(object@meta.data)))
}
if (species == '9606'){
   mitoName<<-'MT'
}else{
  mitoName<<-opt$mtt
}
print(mitoName)
A<<-opt$min
B<<-opt$max
mt_percent<<-opt$mtp
HB_percent<<-opt$hbp
print(A)
print(B)
print(mt_percent)
print(HB_percent)
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = paste("^", mitoName,"-",sep=""))
print('start filter')
object1 <- subset(object,  subset = nFeature_RNA > A & nFeature_RNA < B )
print("filter:min<200,max>10000:")
filterGene <- nrow(object1@meta.data)
print(filterGene)
object <- subset(object1,  percent.mt < mt_percent & percent.HB < HB_percent)
print("filter:mt>20%,HB>5%:")
nGene <- nrow(object@meta.data)
print(nGene)

gene_stat <- c(RawGene,filterGene,nGene)
#gene_result <- data.matrix(gene_stat)
#print(str(gene_result))
names(gene_stat)=c("RawCellbarcodes","Cells_filterGene_200_10000","Cells_mt20_HB5")
#row<-c("RawCellbarcodes","filterGene_Cells(200<x<10000)","Cells(filter:mt>20%,HB>5%)")
#column <- c(opt$sample)
#dimnames(gene_result)=list(row,column)
gene_stats=t(gene_stat)
file <- paste(prefix,'filter_cells_stat.csv',sep='_')
print(file)
print(gene_stats)
write.csv(gene_stats,file,quote=F,row.names=FALSE)
object

## 改为fwrite写入
#raw_data<-as.matrix(object@assays$RNA@counts)
#raw_data  <-  data.frame(data.frame(Gene=rownames(raw_data)),raw_data)
#data.table::fwrite(data.table::data.table(raw_data),file = file,sep = ",",col.names = TRUE,row.names = FALSE,quote=F)
file <- paste(prefix,'filter_cell.csv',sep='_')
print(file)
write.csv(object@meta.data,file,quote=F)
#cell_for_heatmap = round(length(object@meta.data[,1])*0.05)
#print(cell_for_heatmap)
#gene_stat <- c(RawGene,filterGene,nGene)
#gene_result <- data.matrix(gene_stat)
#print(str(gene_result))
#row<-c("RawCells","filterGene_Cells","Cells")
#column <- c(opt$sample)
#dimnames(gene_result)=list(row,column)
#file <- paste(prefix,'filter_gene_stat.csv',sep='_')
#print(file)
#write.table(gene_result,file,sep='\t',quote=F)

##QC
# Expression Violin plot
gene_umi_mito <- paste(prefix,'Gene_UMI_mito_percent.pdf',sep='_')
#unlink('Rplots.pdf')
pdf(gene_umi_mito,width=16,height=8) ## 注意画图输出
#VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), nCol = 3)
plot1 <- VlnPlot(object = object, features = "nFeature_RNA", ncol = 1, combine = "False")
plot2 <- VlnPlot(object = object, features = "nCount_RNA", ncol = 1, combine = "False")
plot3 <- VlnPlot(object = object, features = "percent.mt", ncol = 1, combine = "False")
CombinePlots(plots=c(plot1,plot2,plot3), ncol=3)
print(gene_umi_mito)
dev.off()
# UMI correlation scatter plot
graph<-paste(prefix,'GenePlot.pdf',sep='_')
pdf(graph,width=12,height=8)
par(mfrow = c(1, 2)) 
plot1 <- FeatureScatter(object = object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
print(graph)
dev.off()
